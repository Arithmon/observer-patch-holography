"""Full-register source encoding, native correction and hostile interface checks."""
from fractions import Fraction as F
import math

import numpy as np
import pytest

from ..observer import native
from ..observer.source_discrimination import circuit
from . import encoded_memory as memory, native_routing as routing


def _ideal_code():
    matrix = np.zeros((512, 2), complex)
    matrix[[2, 73], 0] = 1/math.sqrt(2)
    matrix[[74, 1], 1] = 1/math.sqrt(2)
    return matrix


def _execute_integer_edges(edges):
    images = list(range(512))
    for target, pattern in edges:
        mask = 511 ^ (1 << target)
        images = [value ^ (1 << target) if value & mask == pattern else value for value in images]
    return tuple(images)


def test_generic_permutations_have_a_complete_phase_exact_source_schedule():
    assignments = ((0, 2), (2, 73), (1, 74), (3, 1))
    images = routing.complete_permutation(assignments)
    assert all(images[source] == target for source, target in assignments)
    assert set(images) == set(range(512))
    assert _execute_integer_edges(routing.permutation_edges(images)) == images
    inverse = _execute_integer_edges(routing.permutation_edges(images, inverse=True))
    assert all(inverse[images[index]] == index for index in range(512))
    expected = np.zeros((512, 512), complex)
    expected[images, np.arange(512)] = 1
    np.testing.assert_allclose(routing.permutation_unitary(images), expected, atol=3e-11, rtol=0)
    np.testing.assert_allclose(routing.permutation_unitary(images, inverse=True), expected.T,
                               atol=3e-11, rtol=0)


def test_source_encoder_preserves_unknown_logical_phase_and_reference_without_projection():
    expected = _ideal_code()
    np.testing.assert_allclose(memory.code_matrix(), expected, atol=3e-11, rtol=0)
    np.testing.assert_array_equal(memory.ideal_code_matrix(), expected)
    logical = np.array([[.3+.4j, -.2j], [.7j, -.5+.1j]])
    actual = memory.encoding(logical)
    np.testing.assert_allclose(actual, expected@logical, atol=3e-11, rtol=0)
    decoded = memory.decoding(actual)
    target = np.zeros((512, 2), complex)
    target[:2] = logical
    assert decoded.shape == (512, 2)
    np.testing.assert_allclose(decoded, target, atol=3e-11, rtol=0)
    for result in (memory.code_matrix(), actual, decoded):
        with pytest.raises(ValueError):
            result.setflags(write=True)
    # Decode acts unitarily on arbitrary syndrome input, not just the code.
    full = memory.decoding(np.eye(512))
    assert np.linalg.norm(full.conj().T@full-np.eye(512), 'fro') < 2e-10


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('seam', ((0, 1), (0, 5), (2, 3)))
def test_every_raw_native_outcome_has_its_complete_logical_recovery(law, seam):
    branches = memory.repair_branches(law, seam, memory.code_matrix())
    count = 4 if law == 'C' else 8
    assert [(b['coin'], b['record_code']) for b in branches] == [
        (coin, raw) for coin in (0, 1) for raw in range(count)]
    gram = np.zeros((2, 2), complex)
    weight = F(0)
    for branch in branches:
        probability = memory.branch_weight(law, seam, branch['coin'], branch['record_code'])
        assert branch['code_branch_weight'] == probability
        expected = math.sqrt(float(probability))*_ideal_code()
        np.testing.assert_allclose(branch['amplitude'], expected, atol=4e-11, rtol=0)
        gram += branch['amplitude'].conj().T@branch['amplitude']
        weight += probability
        if branch['abort']:
            assert not branch['recovery_applied']
        else:
            unitary = memory.recovery_unitary(law, seam, branch['coin'], branch['record_code'])
            assert np.linalg.norm(unitary.conj().T@unitary-np.eye(512), 'fro') < 2e-10
            if not probability:
                np.testing.assert_array_equal(unitary, np.eye(512))
                assert branch['recovery_resources']['native_pulses'] == 0
    assert weight == 1
    np.testing.assert_allclose(gram, np.eye(2), atol=4e-11, rtol=0)


def test_new_seam_route_executes_all_columns_and_fixes_every_padding_address():
    expected = np.eye(512, dtype=complex)
    for state in native.states():
        target = list(state)
        target[0], target[2] = target[2], target[0]
        target[1], target[3] = target[3], target[1]
        column = circuit.encode(state)
        expected[:, column] = 0
        expected[circuit.encode(tuple(target)), column] = 1
    for inverse in (False, True):
        actual = routing.routing_unitary(seam=(2, 3), inverse=inverse)
        np.testing.assert_allclose(actual, expected, atol=4e-11, rtol=0)
    assert len(routing.routing_pairs((2, 3))) == 20
    assert len(routing.routing_edges(seam=(2, 3))) == 156


@pytest.mark.parametrize('recover', (False, True))
def test_arbitrary_native_padding_and_reference_amplitudes_are_retained(recover):
    state = np.zeros((512, 3), complex)
    state[0, 0], state[511, 1], state[400, 2] = 2+.3j, -.7j, .6
    branches = memory.repair_branches('E', (2, 3), state, recover=recover)
    gram = sum((b['amplitude'].conj().T@b['amplitude'] for b in branches),
               np.zeros((3, 3), complex))
    np.testing.assert_allclose(gram, state.conj().T@state, atol=3e-10, rtol=0)
    assert all(b['amplitude'].shape == (512, 3) for b in branches)
    assert all(not b['recovery_applied'] for b in branches if b['abort'] or not recover)
    # These are unnormalized off-code data, so projection or branchwise
    # normalization would materially change this reference Gram matrix.
    assert np.trace(gram).real > 4


def test_encoder_and_record_recovery_costs_charge_source_calls_without_fresh_quantum_memory():
    encoding, decoding = memory.encoding_resources(), memory.decoding_resources()
    assert encoding['blank_code_qubits'] == 8 and encoding['supplied_logical_qubits'] == 1
    assert decoding['blank_code_qubits'] == 0 and decoding['output_rows_retained'] == 512
    for key in routing.COUNT_KEYS:
        assert encoding[key] == decoding[key]
    assert encoding['native_pulses'] > 11782
    for law in ('C', 'E'):
        for seam in routing.SEAMS:
            for coin in (0, 1):
                for raw in range(4 if law == 'C' else 8):
                    costs = memory.recovery_resources(law, seam, coin, raw)
                    assert costs['blank_code_qubits'] == costs['archive_record_rewrites'] == 0
                    assert costs['new_classical_records'] == 0
                    assert costs['proper_code_service_events'] == 2*costs['native_pulses']
                    assert bool(costs['native_pulses']) == bool(memory.branch_weight(law, seam, coin, raw))


@pytest.mark.parametrize('bad', [True, np.array([True, False]), [1., False],
                                [np.array([1.]), np.array([False])],
                                np.zeros(3), np.empty((2, 0)), [math.nan, 0.],
                                np.ma.array([1., 0.], mask=False)])
def test_bad_logical_inputs_are_rejected_before_source_execution(bad, monkeypatch):
    monkeypatch.setattr(memory, 'code_matrix', lambda: pytest.fail('invalid input reached source execution'))
    with pytest.raises(ValueError):
        memory.encoding(bad)


@pytest.mark.parametrize('law,seam,coin,raw', [
    ('Q', (0, 1), 0, 0), ('C', (1, 0), 0, 0), ('C', [0, 1], 0, 0),
    ('C', (False, 1), 0, 0), ('C', (0, 1), True, 0),
    ('C', (0, 1), 0, False), ('C', (0, 1), 0, 4), ('E', (0, 1), 0, -1)])
def test_malformed_record_labels_are_not_coerced(law, seam, coin, raw):
    for operation in (memory.branch_weight, memory.recovery_unitary, memory.recovery_resources):
        with pytest.raises(ValueError):
            operation(law, seam, coin, raw)


@pytest.mark.parametrize('assignments', [((0, 1), (0, 2)), ((0, 1), (2, 1)),
                                       ((False, 1),), ((0, 512),), [(0, 1)]])
def test_ambiguous_or_malformed_permutation_completions_are_rejected(assignments):
    with pytest.raises(ValueError):
        routing.complete_permutation(assignments)


@pytest.mark.parametrize('bad', [tuple(range(511)), tuple([0]*512),
                               (False, *range(1, 512)), list(range(512))])
def test_invalid_full_permutations_are_rejected_before_source_execution(bad):
    with pytest.raises(ValueError):
        routing.permutation_unitary(bad)


def test_abort_is_not_a_recovery_label_and_boolean_control_is_explicit():
    assert memory.branch_weight('C', (0, 1), 0, 3) == 0
    with pytest.raises(ValueError, match='halt'):
        memory.recovery_unitary('C', (0, 1), 0, 3)
    with pytest.raises(ValueError, match='Boolean'):
        memory.repair_branches('C', (0, 1), np.zeros(512), recover=1)

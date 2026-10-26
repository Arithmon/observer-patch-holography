"""Independent whole-output tests of the narrow source instrument adapter."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.observer.source_scalar import compiler as c
from . import history as h


X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def embed(local, targets, work_bits):
    """Independent index embedding, with all unmeasured bits unchanged."""
    size = 1 << work_bits
    result = np.zeros((size, size), complex)
    mask = sum(1 << bit for bit in targets)
    for column in range(size):
        incoming = sum(((column >> bit) & 1) << j for j, bit in enumerate(targets))
        for outgoing in range(len(local)):
            row = (column & ~mask) | sum(((outgoing >> j) & 1) << bit
                                         for j, bit in enumerate(targets))
            result[row, column] = local[outgoing, incoming]
    return result


def reference_columns(size, columns=3):
    values = np.arange(size * columns).reshape(size, columns)
    state = np.sin(values + .2) + 1j * np.cos(values * 1.3 + .7)
    return state / np.linalg.norm(state)


@pytest.mark.parametrize("vector", [np.array([1., 2., 0., -3.]) / math.sqrt(14),
                                    np.array([1., 2j, 0., -3. + 1j]) / math.sqrt(15)])
def test_known_complex_preparation_is_an_executed_word(vector):
    zero = np.zeros((len(vector), 1), complex)
    zero[0, 0] = 1
    word = h.preparation_word(vector)
    output, resources = h.execute_word(word, 2, zero)
    np.testing.assert_allclose(output[:, 0], vector, rtol=0, atol=2e-13)
    assert resources["native_pulses"] == len(word)
    assert not resources["physical_service_error_certified"]


def test_weak_read_preserves_all_reference_columns_and_clock_spectator():
    state = reference_columns(8)
    effects = np.array([0., .2, .7, 1.])
    basis_word = [c.pulse("Y", .37, [0]), c.pulse("Z", -.21, [1])]
    local_basis = np.kron(expm(-.21j * Z), expm(.37j * Y))
    result = h.diagonal_read(state, effects, targets=(2, 0), basis_word=basis_word)
    for bit, probability in enumerate((effects, 1 - effects)):
        local = local_basis.conj().T @ np.diag(np.sqrt(probability)) @ local_basis
        expected = embed(local, (2, 0), 3) @ state
        np.testing.assert_allclose(result.outputs[bit], expected, rtol=0, atol=3e-13)
    assert sum(np.vdot(v, v).real for v in result.outputs) == pytest.approx(1, abs=2e-13)
    assert result.resources["fresh_pointer_blanks"] == 1
    assert result.resources["committed_classical_bits"] == 1


def test_projector_restores_arbitrary_complex_vector_and_keeps_spectators():
    vector = np.array([1., 2j, -1j, -.7 + .3j])
    vector /= np.linalg.norm(vector)
    state = reference_columns(8, 2)
    result = h.projector_read(state, vector)
    projector = embed(np.outer(vector, vector.conj()), (0, 1), 3)
    np.testing.assert_allclose(result.outputs[0], projector @ state, rtol=0, atol=5e-13)
    np.testing.assert_allclose(result.outputs[1], (np.eye(8) - projector) @ state,
                               rtol=0, atol=5e-13)


def test_zero_probability_outcome_is_retained_and_no_branch_is_renormalized():
    vector = np.array([1., 1j, -2., 0.]) / math.sqrt(6)
    state = vector[:, None] / 7
    result = h.projector_read(state, vector)
    assert len(result.outputs) == 2
    np.testing.assert_allclose(result.outputs[0], state, rtol=0, atol=2e-13)
    assert np.linalg.norm(result.outputs[1]) < 2e-13
    assert sum(np.vdot(v, v).real for v in result.outputs) == pytest.approx(1 / 49, abs=2e-13)
    zero = h.diagonal_read(np.zeros((4, 2), complex), [.1, .2, .3, .4])
    assert len(zero.outputs) == 2 and all(np.count_nonzero(v) == 0 for v in zero.outputs)


@pytest.mark.parametrize("axis,matrix", [("X", X), ("Y", Y), ("Z", Z)])
def test_clock_pauli_reads_are_complete_qnd_instruments(axis, matrix):
    state = reference_columns(8, 4)
    result = h.pauli_read(state, 1, axis)
    for bit, sign in enumerate((1, -1)):
        projector = embed((np.eye(2) + sign * matrix) / 2, (1,), 3)
        np.testing.assert_allclose(result.outputs[bit], projector @ state, rtol=0, atol=3e-13)
    # No measurement branch substitutes or traces out a reference.
    assert all(value.shape == state.shape for value in result.outputs)


def test_sequential_complete_reads_preserve_old_record_marginals():
    state = reference_columns(8, 2)
    first = h.diagonal_read(state, [.2, .8], targets=(0,))
    records = []
    for old_bit, branch in enumerate(first.outputs):
        # A source operation on the same retained branch precedes another read.
        continued, _ = h.execute_word([c.pulse("Y", .31, [2])], 3, branch)
        later = h.pauli_read(continued, 2, "X")
        assert sum(np.vdot(v, v).real for v in later.outputs) == pytest.approx(
            np.vdot(branch, branch).real, abs=3e-13)
        records.extend(((old_bit, bit), value) for bit, value in enumerate(later.outputs))
    assert [record for record, _ in records] == [(0, 0), (0, 1), (1, 0), (1, 1)]
    assert sum(np.vdot(value, value).real for _, value in records) == pytest.approx(1, abs=4e-13)


def test_unitary_word_preserves_reference_reduction():
    state = reference_columns(8, 3)
    output, _ = h.execute_word([c.pulse("X", .24, [0]),
                                c.pulse("J", -.29, [1, 2])], 3, state)
    np.testing.assert_allclose(output.conj().T @ output, state.conj().T @ state,
                               rtol=0, atol=2e-13)


@pytest.mark.parametrize("effect", [[-.01, 1.], [0., 1.01], [0., math.nan],
                                     [True, False], [1., 0., 0.]])
def test_invalid_effects_fail_closed(effect):
    with pytest.raises(ValueError):
        h.diagonal_read(np.eye(4, dtype=complex)[:, :1], effect)


def test_invalid_vectors_targets_and_qubit_counts_fail_closed():
    state = reference_columns(4)
    with pytest.raises(ValueError):
        h.preparation_word([1., 1.])
    with pytest.raises(ValueError):
        h.diagonal_read(state, [0., 1.], targets=(True,))
    with pytest.raises(ValueError):
        h.diagonal_read(state, [0., 1.], basis_word=[c.pulse("X", .2, [1])])
    with pytest.raises(ValueError):
        h.pauli_read(state, 2, "X")
    with pytest.raises(ValueError):
        h.pauli_read(state, 1, "A")
    with pytest.raises(ValueError):
        h.execute_word([], 3, state)


def test_eight_bit_projector_keeps_all_scratch_outputs_on_a_reference():
    vector = np.zeros(256, complex)
    vector[5], vector[194] = 1 / math.sqrt(2), 1j / math.sqrt(2)
    # Eight measured work bits, two distinct clock spectators, and an
    # untouched purification reference. The first clock/reference pair
    # carries a Bell correlation; another component populates clock bit9.
    state = np.zeros((1024, 2), complex)
    state[5, 0], state[194 + 256, 1] = .5, .5j
    state[17 + 512, 0] = 1 / math.sqrt(2)
    result = h.projector_read(state, vector)
    projector = np.outer(vector, vector.conj())
    selected = np.vstack([projector @ block for block in np.split(state, 4)])
    np.testing.assert_allclose(result.outputs[0], selected, rtol=0, atol=3e-12)
    np.testing.assert_allclose(result.outputs[1], state - selected, rtol=0, atol=3e-12)
    assert len(result.failure_outputs) == len(result.failure_records) == 6
    assert result.failure_records == ((0, 1, 0), (1, 1, 0), (0, 0, 1),
                                      (1, 0, 1), (0, 1, 1), (1, 1, 1))
    assert result.scratch_return_error < 3e-12
    # The long executed word retains its floating trace defect. Use the
    # instrument's declared mass tolerance; the branch/amplitude and scratch
    # checks above remain stricter and no output is renormalized.
    assert sum(np.vdot(branch, branch).real for branch in
               (*result.outputs, *result.failure_outputs)) == pytest.approx(
                   np.vdot(state, state).real, rel=0, abs=h.TOLERANCE)
    assert result.resources["fresh_pointer_blanks"] == 3
    assert result.resources["committed_classical_bits"] == 3
    assert result.resources["maximum_work_qubits"] == 13
    assert result.resources["scratch_failure_outputs_retained"]


def test_overflowed_mass_masks_and_mixed_booleans_are_rejected():
    for state in (np.full((2, 1), 1e308 + 1e308j),
                  np.ma.array([[1.], [0.]], mask=[[False], [True]]),
                  [[True], [0.]]):
        with pytest.raises(ValueError):
            h.execute_word([], 1, state)
    with pytest.raises(ValueError):
        h.preparation_word([1e308 + 1e308j, 0.])
    with pytest.raises(ValueError):
        h.preparation_word([True, 0.])
    with pytest.raises(ValueError):
        h.diagonal_read(np.array([[1.], [0.]]), [True, .2])


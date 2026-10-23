"""Independent native Kraus and arbitrary-reference pointer regressions."""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import pointer as p


def _full_native_oracle(factor, law):
    """Execute full78 native branches after the raw word, without effect tables."""
    states = p.native.states()
    index = {state: i for i, state in enumerate(states)}
    reference, rank = factor.shape[0]//2, factor.shape[1]
    coupled = np.zeros((78, 2, reference, rank), complex)
    matrix = p.local_matrix()
    for z in range(4):
        for a in (0, 1):
            for c in (0, 1):
                coupled[index[p.compiled.circuit.decode(z)], a] += (
                    matrix[z+4*a, p.READY_CODE+4*c]
                    *factor.reshape(2, reference, rank)[c])
    groups, visible = defaultdict(list), {}
    for initial, state in enumerate(states):
        for coin in (0, 1):
            # Independent endpoint floor/ceiling rule, including collisions.
            total = state[0]+state[1]
            low, high = total//2, (total+1)//2
            out = list(state)
            out[0], out[1] = (low, high) if coin == 0 else (high, low)
            current = state[0]-out[0]
            key = coin, current if law == 'C' else state[0]-state[1]
            visible[key] = int(2*current+2*coin-1 < 0)
            groups[key].append((initial, index[tuple(out)]))
    result = [np.zeros((2*reference, 2*reference), complex) for _ in (0, 1)]
    for key, mapping in groups.items():
        transformed = np.zeros_like(coupled)
        for initial, output in mapping:
            transformed[output] += coupled[initial]/np.sqrt(2)
        tensor = transformed.reshape(78, 2*reference, rank)
        result[visible[key]] += np.einsum('air,ajr->ij', tensor, tensor.conj())
    return tuple(result)


def _trace_norm(matrix):
    return float(np.sum(abs(np.linalg.eigvalsh(matrix))))


def test_exact_native_instrument_and_raw_source_envelope():
    certificate = p.certificate()
    audit = certificate['exact_instrument']
    assert audit['logical_superoperator_entries_checked'] == 64
    assert audit['ready_native_configuration'] == [1, 0, 1]+[0]*9
    assert certificate['source_pulses'] == 6
    assert certificate['proper_code_service_events'] == 12
    assert certificate['native_pointer_attempts'] == 1
    assert certificate['half_diamond_upper'] == p.WORD_ERROR_GATE*(2+p.WORD_ERROR_GATE)/2
    for law in ('C', 'E'):
        assert p.native_effects(law) == ((1, 0, 2, 0), (1, 2, 0, 2))
    assert len(p.source_word()) == 6
    assert not p.local_matrix().flags.writeable
    # Wrong pointer ready orientation reverses the intended input label.
    wrong = p.target()[:, 1]
    assert sum(p.native_effects('C')[1][z]*wrong[z]**2 for z in range(4)) == 2


@pytest.mark.parametrize('law', ['C', 'E'])
def test_full78_kraus_execution_matches_reduced_map_with_qutrit_reference(law):
    rng = np.random.default_rng(729)
    factor = rng.normal(size=(6, 3))+1j*rng.normal(size=(6, 3))
    factor *= np.sqrt(.37)/np.linalg.norm(factor)
    density = factor@factor.conj().T
    actual = p.read_bit(density, law, bit=0, work_qubits=1, reference_dimension=3)
    oracle = _full_native_oracle(factor, law)
    for observed, expected in zip(actual, oracle, strict=True):
        np.testing.assert_allclose(observed, expected, atol=2e-16, rtol=0)
    assert sum(np.trace(value).real for value in actual) == pytest.approx(.37, abs=2e-15)


@pytest.mark.parametrize('bit', [0, 1, 2])
def test_all_untouched_work_and_reference_indices_survive(bit):
    rng = np.random.default_rng(100+bit)
    vector = rng.normal(size=24)+1j*rng.normal(size=24)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    masks = np.repeat([(i >> bit) & 1 for i in range(8)], 3)
    for law in ('C', 'E'):
        actual = p.read_bit(density, law, bit=bit, work_qubits=3, reference_dimension=3)
        expected = [density*(masks[:, None] == value)*(masks[None, :] == value)
                    for value in (0, 1)]
        distance = sum(_trace_norm(a-b) for a, b in zip(actual, expected, strict=True))/2
        assert distance < float(p.certificate()['half_diamond_upper'])


def test_full512_register_includes_padding_and_preserves_coherence():
    report = p.run()
    for witness in report['full_register_reference_executions']:
        assert witness['work_dimension'] == 512 and witness['reference_dimension'] == 2
        assert witness['complete_output_max_entry_error_diagnostic'] < 2e-14
        assert witness['mass_error_diagnostic'] < 2e-14
    assert report['scope']['proper_code_service_blanks_not_derived']


def test_correlated_auxiliary_reset_can_steer_but_complete_timeout_bound_holds():
    # The auxiliary label is classically correlated with two orthogonal WR
    # states. A source reset with label-dependent success steers the WR input.
    plus = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    minus = np.array([0, 1, 1, 0], complex)/np.sqrt(2)
    first, second = np.outer(plus, plus.conj())/2, np.outer(minus, minus.conj())/2
    original = first+second
    failure = second/4
    success = original-failure
    assert _trace_norm(success/np.trace(success)-original)/2 > .05
    masks = np.repeat([0, 1], 2)
    ideal = [original*(masks[:, None] == b)*(masks[None, :] == b) for b in (0, 1)]
    for law in ('C', 'E'):
        actual = p.read_bit(success, law, bit=0, work_qubits=1, reference_dimension=2)
        complete = (sum(_trace_norm(a-b) for a, b in zip(actual, ideal, strict=True))
                    +_trace_norm(failure))/2
        assert complete == pytest.approx(.125, abs=2e-15)
        assert complete < float(p.comparison_bound(F(1, 4)))


def test_zero_substate_and_tiny_source_mass_roundoff_are_not_renormalized():
    empty = np.zeros((2, 2), complex)
    for value in p.read_bit(empty, 'C', bit=0, work_qubits=1):
        assert not np.any(value)
    density = np.diag([1+2e-14, 0]).astype(complex)
    outputs = p.read_bit(density, 'E', bit=0, work_qubits=1)
    assert sum(np.trace(value).real for value in outputs) > 1+1e-14


@pytest.mark.parametrize('diagonal', [[1., -9e-11], [1.]+[-5e-14]*255])
def test_negative_input_mass_cannot_hide_in_per_eigenvalue_tolerance(diagonal):
    density = np.diag(diagonal)
    with pytest.raises(ValueError, match='aggregate roundoff'):
        p.read_bit(density, 'C', bit=0, work_qubits=(len(diagonal)-1).bit_length())


def test_dense_rank_one_reference_roundoff_passes_without_changing_input():
    rng = np.random.default_rng(1125)
    vector = rng.normal(size=1024)+1j*rng.normal(size=1024)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    before = density.copy()
    outputs = p.read_bit(density, 'C', reference_dimension=2)
    np.testing.assert_array_equal(density, before)
    assert abs(sum(np.trace(value) for value in outputs)-1) < 1e-13


@pytest.mark.parametrize('field,value', [
    ('mass_error_diagnostic', -1e-10),
    ('mass_error_diagnostic', 1e-10),
    ('complete_output_max_entry_error_diagnostic', -1e-10),
    ('complete_output_max_entry_error_diagnostic', 1e-10),
])
def test_portable_replay_cannot_accept_invalid_diagnostic_domains(field, value):
    actual = json.loads(p.OUTPUT.read_text(encoding='utf-8'))
    changed = deepcopy(actual)
    changed['full_register_reference_executions'][0][field] = value
    with pytest.raises(ValueError, match='execution domain'):
        p.compare_receipt(changed, actual)


@pytest.mark.parametrize('kind', ['negative', 'inflated', 'nan', 'inf', 'object',
                                  'integer', 'bool', 'masked', 'subclass', 'nonhermitian'])
def test_hostile_density_is_rejected(kind):
    density = np.diag([.5, .5]).astype(complex)
    if kind == 'negative':
        density = np.diag([2., -1.])
    elif kind == 'inflated':
        density *= 1.001
    elif kind in ('nan', 'inf'):
        density[0, 0] = float(kind)
    elif kind in ('object', 'integer', 'bool'):
        dtype = {'object': object, 'integer': int, 'bool': bool}[kind]
        density = np.eye(2, dtype=dtype)
    elif kind == 'masked':
        density = np.ma.array(density)
    elif kind == 'subclass':
        class Forged(np.ndarray):
            pass
        density = density.view(Forged)
    else:
        density[0, 1] = .1j
    with pytest.raises(ValueError):
        p.read_bit(density, 'C', bit=0, work_qubits=1)


@pytest.mark.parametrize('keywords', [{'law': 'X'}, {'bit': True}, {'work_qubits': True},
                                     {'reference_dimension': True}, {'reference_dimension': 4},
                                     {'bit': 1}, {'work_qubits': 0}])
def test_malformed_interface_addresses_are_rejected(keywords):
    options = dict(law='C', bit=0, work_qubits=1)
    options.update(keywords)
    with pytest.raises(ValueError):
        p.read_bit(np.eye(2)/2, **options)


def test_raw_source_mutation_cannot_pass_certificate(monkeypatch):
    changed = np.array(p.local_matrix())
    changed[0, 0] += 1e-7
    monkeypatch.setattr(p, 'local_matrix', lambda: changed)
    with pytest.raises(ArithmeticError, match='operator envelope'):
        p.certificate()


@pytest.mark.parametrize('timeout', [True, .1, F(-1, 2), F(3, 2)])
def test_comparison_budget_requires_exact_valid_timeout(timeout):
    with pytest.raises(ValueError):
        p.comparison_bound(timeout)

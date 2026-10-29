"""Complete source ancestry, correlated memory and record-only coefficient recovery."""
from unittest.mock import patch

import numpy as np
import pytest
import sympy as s

from . import native_family_record_memory as m


@pytest.fixture(scope='module')
def data():
    return m.architecture()


@pytest.fixture(scope='module')
def execution():
    return m.execute()


def test_native_labels_and_actual_event_union_are_complete(data):
    assert (data['p'], data['q'], data['coincidence']) == (s.Rational(1, 12), s.Rational(1, 6), s.Rational(1, 72))
    assert len(data['records']) == 128 and len(data['labels']) == 1536
    assert {r[4] for r in data['labels']} == set(range(12))
    assert tuple(r[:4] for r in data['records']) == m.RAW_LABELS
    assert any(r[-1] == 0 for r in data['records'])
    assert sum(r[-1] for r in data['records'] if r[2:4] == (0, 4)) == s.Rational(1, 12)
    assert sum(r[-1] for r in data['records'] if r[2:4] == (0, 5)) == s.Rational(1, 12)


def test_core_reset_and_exact_stationary_law(data):
    expected = s.zeros(12, 1)
    for port, weight in [(2, s.Rational(1206, 1291)), (3, s.Rational(335, 5164)),
                         (9, s.Rational(25, 397628)), (10, s.Rational(90, 99407))]:
        expected[port] = weight
    assert data['stationary'] == expected
    assert data['T']*expected == expected and sum(expected) == 1
    assert data['T'][0, :] == data['T'][1, :] == s.zeros(1, 12)
    assert all(sum(data['T'][:, j]) == 1 for j in range(12))
    support = [2, 3, 9, 10]
    assert data['T'].extract(support, support) == s.Matrix([
        [67, 71, 67, 67], [5, 0, 0, 0], [0, 0, 0, 5], [0, 1, 5, 0]])/72


def test_exact_kernel_invariants_and_uncalibrated_diagnostics(data):
    Y = data['Y']
    assert Y == Y.T and s.trace(Y) == 1 and Y.det() > 0
    for P in data['projectors']:
        assert objects_zero(P*P-P) and s.trace(P) == 1
    values = np.linalg.eigvalsh(np.asarray(Y, float))
    np.testing.assert_allclose(values/values[0], [1, 195.11703757, 3579.82624571], rtol=2e-10)
    diagnostics = m.kernel_diagnostics()
    assert not diagnostics['physical_pole_masses_inferred']
    assert s.simplify(s.sympify(diagnostics['determinant_exact'])-Y.det()) == 0


def objects_zero(matrix):
    return all(s.simplify(x) == 0 for x in matrix)


@pytest.mark.parametrize('change', ['missing_zero', 'relabel', 'flag', 'renormalize'])
def test_source_interface_rejects_dropped_or_changed_raw_inputs(data, change):
    rows = list(data['records'])
    if change == 'missing_zero':
        rows.pop(next(i for i, row in enumerate(rows) if row[-1] == 0))
    elif change == 'relabel':
        rows[0] = (0, 0, 0, 1, *rows[0][4:])
    elif change == 'flag':
        rows[0] = (*rows[0][:4], True, rows[0][-1])
    else:
        rows[0] = (*rows[0][:-1], rows[0][-1]+s.Rational(1, 100))
    with pytest.raises(ValueError):
        m.transition_from_records(rows)


def test_supplied_interaction_factors_preserve_exact_expressions_and_scope(data):
    p = s.sin(s.pi/5)**2/12
    q = s.sin(s.sqrt(2)*s.pi/5)**2/6
    other = m.architecture_from_factors(p, q)
    assert other['p'] == p and other['q'] == q
    assert not other['physical_fermion_mass_map_derived']
    assert 'caller-supplied' in other['event_factor_origin']
    T, pi = np.asarray(other['T'].evalf(30), float), np.asarray(other['stationary'].evalf(30), float).ravel()
    np.testing.assert_allclose(T@pi, pi, atol=1e-15, rtol=0)
    assert m.architecture_from_factors(data['p'], data['q'])['Y'] == data['Y']
    with pytest.raises(ValueError):
        m.architecture_from_factors(1, q)


@pytest.mark.parametrize('p,q', [
    (-s.Rational(1, 10**400), s.Rational(1, 2)),
    (-s.Rational(1, 10**20), s.Rational(1, 2)),
    (1, s.Rational(1, 2)),
    (1+s.Rational(1, 10**20), s.Rational(1, 2)),
    (s.Rational(1, 2), -s.Rational(1, 10**400)),
    (s.Rational(1, 2), 1+s.Rational(1, 10**20)),
    (s.Rational(1, 2), 1+s.Rational(1, 10**400)),
])
def test_exact_out_of_domain_factors_cannot_construct_signed_transitions(p, q):
    with pytest.raises(ValueError, match='declared domain'):
        m.architecture_from_factors(p, q)


def test_exact_near_boundary_values_keep_their_domain_without_float_coercion():
    tiny = s.Rational(1, 10**400)
    assert m._probability(tiny, strict_upper=True) == tiny
    assert m._probability(1-tiny, strict_upper=True) == 1-tiny
    assert m._probability(0, strict_upper=True) == 0
    assert m._probability(1) == 1
    assert m._probability(np.float64('0.125')) == s.Rational(1, 8)
    assert m._probability(s.sqrt(2)/2) == s.sqrt(2)/2


@pytest.mark.parametrize('p,q', [
    (0, 1),
    (s.Rational(1, 2), 1),
    (1-s.Rational(1, 10**20), s.Rational(1, 2)),
])
def test_admitted_endpoint_and_near_endpoint_kernels_are_exactly_stochastic(p, q):
    result = m.architecture_from_factors(p, q)
    assert result['p'] == p and result['q'] == q
    assert all(entry.is_nonnegative is True for entry in result['T'])
    assert all(sum(result['T'][:, column]) == 1 for column in range(12))
    assert result['T']*result['stationary'] == result['stationary']
    assert sum(result['stationary']) == 1


@pytest.mark.parametrize('value', [
    s.nan, s.oo, -s.oo, s.zoo, s.I, s.Symbol('p', real=True),
    s.Function('unknown_factor', real=True, finite=True)(0),
    np.float64('nan'), np.float64('inf'), True, None,
])
def test_nonfinite_nonreal_or_undetermined_factors_fail_closed(value):
    with pytest.raises(ValueError, match='finite real event factor'):
        m.architecture_from_factors(s.Rational(1, 2), value)


def test_raw_record_weights_reject_negative_underflow_even_when_normalized(data):
    rows = list(data['records'])
    zero = next(i for i, row in enumerate(rows) if row[-1] == 0)
    positive = next(i for i, row in enumerate(rows) if row[-1] > 0)
    tiny = s.Rational(1, 10**400)
    rows[zero] = (*rows[zero][:-1], -tiny)
    rows[positive] = (*rows[positive][:-1], rows[positive][-1]+tiny)
    assert sum(row[-1] for row in rows) == 1
    with pytest.raises(ValueError, match='declared domain'):
        m.transition_from_records(rows)


def test_all1536_effective_kraus_factors_preserve_every_memory_input(data):
    factors, labels = m.raw_instrument_factors(1/64)
    assert factors.shape == (1536, 6, 6) and labels.shape == (1536, 6)
    grams = np.zeros((12, 6, 6), complex)
    for K, label in zip(factors, labels, strict=True):
        grams[label[4]] += K.conj().T@K
    np.testing.assert_allclose(grams, np.broadcast_to(np.eye(6), grams.shape), atol=1e-14, rtol=0)
    zero = [i for i in range(1536) if data['records'][i//12][-1] == 0]
    assert np.max(abs(factors[zero])) == 0


def test_complete_finite_CQ_channel_is_CP_TP_and_correlated(data):
    # Apply the finite channel to all36 quantum matrix units with the full
    # stationary product memory, retaining all12 terminal memory outputs.
    N = 8
    joint = m.joint_channel(1/N)
    initial = np.vstack([float(pi)*np.eye(36) for pi in data['stationary']])
    terminal = (np.linalg.matrix_power(joint, N)@initial).reshape(12, 36, 36)
    marginal = terminal.sum(axis=0)
    choi = marginal.reshape(6, 6, 6, 6).transpose(0, 2, 1, 3).reshape(36, 36)
    np.testing.assert_allclose(choi, choi.conj().T, atol=1e-14, rtol=0)
    assert np.linalg.eigvalsh(choi).min() > -1e-13
    np.testing.assert_allclose(np.eye(6).reshape(36)@marginal, np.eye(6).reshape(36), atol=1e-13, rtol=0)
    from scipy.linalg import expm
    U = expm(-1j*np.kron([[0, 1], [1, 0]], np.asarray(data['Y'], float)))
    assert 6*np.linalg.norm(marginal-np.kron(U, U.conj())) < m.averaging_bound(N, initial_memory_stationary=True)
    assert all(sum(abs(x) for x in data['T'][:, i]-data['T'][:, j]) <= 2*data['p']
               for i in range(12) for j in range(12))


def test_memory_has_stationary_current_but_no_coherent_phase_claim(data):
    diagnostic = m.stationary_currents()
    assert not diagnostic['detailed_balanced']
    assert diagnostic['current'][3, 2] == s.Rational(335, 371808)
    assert data['T'][3, 2]*data['T'][10, 3]*data['T'][2, 10] > 0
    assert data['T'][10, 2]*data['T'][3, 10]*data['T'][2, 3] == 0
    assert not diagnostic['classical_current_is_coherent_U1_phase']


def test_fixed_basis_inverse_recovers_coefficient_from_complete_records(execution, data):
    report, records = execution
    law = records['probabilities']
    assert law.shape == (2, 9, 1536, 3)
    np.testing.assert_allclose(law.sum(axis=(2, 3)), np.ones((2, 9)), atol=1e-14, rtol=0)
    assert np.any(law[:, :, :, 2] > 0)
    with patch.object(m, 'architecture', side_effect=AssertionError('inverse read source')):
        inferred = m.decode_records(law, steps=records['steps'], operator_error=report['error_budget']['operator_error_upper'])
    assert np.linalg.norm(inferred['estimate']-np.asarray(data['Y'], float), 2) <= inferred['operator_error_upper']
    for value, (lo, hi) in zip(np.linalg.eigvalsh(np.asarray(data['Y'], float)), inferred['eigenvalue_intervals'], strict=True):
        assert lo <= value <= hi
    assert inferred['eigenvalue_intervals'][0][0] > 0
    assert not report['scope']['literal_native_pulse_history_compiled']
    assert not report['scope']['physical_pole_mass_map_derived']
    assert not report['scope']['measured_spectrum_claimed']


def test_forward_path_does_not_read_kernel_or_spectrum(execution, data):
    without_kernel = {key: value for key, value in data.items() if key != 'Y'}
    with patch.object(m, 'architecture', return_value=without_kernel), \
            patch.object(np.linalg, 'eigvalsh', side_effect=AssertionError('forward read spectrum')):
        other = m.forward_records()
    np.testing.assert_array_equal(other['probabilities'], execution[1]['probabilities'])


def test_inverse_rejects_lost_terminal_record_without_renormalization(execution):
    law = execution[1]['probabilities'].copy()
    # Choose a definitely positive branch instead of assuming a particular
    # raw label is supported: zero labels must also remain in the array.
    i = int(np.argmax(execution[1]['probabilities'][0, 0].sum(axis=1)))
    law[0, 0, i] = 0
    with pytest.raises(ValueError, match='complete'):
        m.decode_records(law, steps=m.READOUT_STEPS, operator_error=1e-5)


def test_two_clock_inverse_cancels_bias_before_symmetrization(execution, data):
    report, records = execution
    assert report['error_budget']['finite_rank_one_operator_bias_upper'] == 0
    h, k = m.READOUT_STEPS
    u, v = report['error_budget']['two_clock_inverse_weights']
    assert abs(u+v-1) < 1e-14
    assert abs(u*(np.cos(h)-1)+v*(np.cos(k)-1)) < 1e-14
    finite_pi = np.asarray(data['T'], float)@records['burnin_memories'][-1]
    actual = sum(finite_pi[j]*np.asarray(data['projectors'][j], float) for j in range(12))
    assert np.linalg.norm(report['inferred_kernel']['estimate']-actual, 2) < 1e-12


@pytest.mark.parametrize('change', ['zero_label', 'clock', 'memory', 'normalized_replacement'])
def test_rehashed_or_self_consistent_replacement_laws_fail_source_replay(execution, change):
    records = {key: value.copy() if isinstance(value, np.ndarray) else value
               for key, value in execution[1].items()}
    if change == 'zero_label':
        records['labels'][0, 0] += 1
    elif change == 'clock':
        records['steps'] = (1/9, 2/9)
    elif change == 'memory':
        records['burnin_memories'][0] = np.eye(12)[3]
    else:
        # Preserve normalization and allow a blind inverse to decode its own
        # result; the replacement still fails the independent source replay.
        law = records['probabilities']
        i = int(np.argmax(law[0, 0, :, 0]))
        law[0, 0, i, 0] -= 1e-6
        law[0, 0, i, 1] += 1e-6
        m.decode_records(law, steps=m.READOUT_STEPS, operator_error=1e-5)
    with pytest.raises(ValueError):
        m.validate_records(records)

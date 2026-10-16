"""Independent stopped-map, modular convention and continuing-record controls."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as c


@pytest.mark.parametrize('law', ['C', 'E'])
def test_stopped_channel_trace_and_protocol_reduction_on_every_matrix_unit(law):
    columns = c.event_columns(law)
    assert len(columns) == 66**2
    for i in range(66):
        for j in range(66):
            column = columns[66*i+j]
            assert all(value > 0 for _, _, value in column)
            assert sum(value for a, b, value in column if a == b) == int(i == j)
            # Both S(A) subset A and S†(A) subset A: no matrix coefficient
            # crosses either side of the Hilbert-Schmidt orthogonal split.
            assert all(((i, j) in c.ALGEBRA) == ((a, b) in c.ALGEBRA)
                       for a, b, _ in column)
    with pytest.raises(TypeError):
        columns[0][0] = (0, 0, F(0))


@pytest.mark.parametrize('law', ['C', 'E'])
def test_geometric_sum_matches_actual_native_kraus_matrices(law):
    source = c.marked(law)
    palm = c.events.embedded()[1]
    # Summing the exact Schur geometric series numerically independently of
    # event_columns, then contracting the actual native T1 Kraus matrices.
    w = np.asarray(source.no_current_counts)/60
    states = [np.diag([float(p) for p in palm]), c.dense(c.prepared_states(law)[1])]
    for rho in states:
        actual = sum(k@(rho/(1-w))@k.conj().T for k in source.instrument.kraus[1])
        sparse = {(i, j): c.z(F(float(rho[i, j].real)), F(float(rho[i, j].imag)))
                  for i in range(66) for j in range(66) if rho[i, j]}
        np.testing.assert_allclose(actual, c.dense(c.apply(law, sparse)), rtol=0, atol=5e-16)
    # No-current first, movement second: order is T1 T0, not T0 T1.
    initial = c.dense(c.prepared_states(law)[1])
    waited = sum(k@initial@k.conj().T for k in source.instrument.kraus[0])
    second = sum(k@waited@k.conj().T for k in source.instrument.kraus[1])
    np.testing.assert_allclose(second, c.dense(c.apply(law, c.prepared_states(law)[1], first_at=2)),
                               rtol=0, atol=5e-17)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_stationary_reference_is_event_palm_not_attempt_stationarity(law):
    degree, palm, _, transition = c.events.embedded()
    sigma = c.diagonal(palm)
    assert c.apply(law, sigma) == sigma
    attempt = c.apply(law, sigma, no_current=1)
    for key, value in c.apply(law, sigma, first_at=1).items():
        c.accumulate(attempt, key, value)
    assert sum(abs(attempt[i, i][0]-palm[i]) for i in range(66))/2 == F(1, 150)
    assert c.apply(law, c.diagonal([F(1, 66)]*66)) != sigma
    # Native one-attempt movement herald, not a supplied Palm preparation.
    herald = c.apply(law, c.diagonal([F(1, 66)]*66), first_at=1)
    assert c.trace(herald) == F(5, 33)
    assert herald == c.diagonal([F(5, 33)*p for p in palm])
    assert all(F(degree[i], 600)*transition[i][j] == F(degree[j], 600)*transition[j][i]
               for i in range(66) for j in range(66))


def test_modular_dual_state_sign_and_exact_protocol_domain():
    degree, palm, _, _ = c.events.embedded()
    full = np.diag([1j if d == 8 else 1 for d in degree])
    for key in c.ALGEBRA:
        unit = {key: c.z(1)}
        np.testing.assert_array_equal(c.dense(c.modular_action(unit)), full@c.dense(unit)@full.conj().T)
        assert c.modular_action(unit) == c.modular_action(unit, full=False)
    # Omitted degree-eight phases cancel only in the stated algebra.
    outside = {(0, 1): c.z(1)}
    assert (0, 1) not in c.ALGEBRA
    assert c.modular_action(outside) != c.modular_action(outside, full=False)
    unequal = next((i, j) for i, j in c.PAIRS if palm[i] < palm[j])
    unit = {unequal: c.z(1)}
    assert c.modular_action(unit)[unequal] == c.z(0, 1)
    assert c.modular_action(unit, -1)[unequal] == c.z(0, -1)


def test_full_matrix_kms_boundary_has_inverse_observable_convention():
    palm = c.events.embedded()[1]
    for i in range(66):
        for j in range(66):
            # A=Eij, B=Eji. alpha_(-i)(B)=sigma B sigma^-1.
            lhs = palm[i]*palm[j]/palm[i]
            rhs = palm[j]
            assert lhs == rhs
    i, j = next((i, j) for i, j in c.PAIRS if palm[i] != palm[j])
    assert palm[i]*palm[i]/palm[j] != palm[j]  # wrong KMS sign fails


@pytest.mark.parametrize('law', ['C', 'E'])
def test_actual_native_branches_and_apparatus_reduce_the_algebra(law):
    result = c.algebra_audit(law)
    assert result['source_map_and_adjoint_preserve_protocol_algebra']
    assert result['nonzero_full_channel_modular_frequency_mixing_entries'] == 3360
    assert result['Palm_reference_one_native_attempt_TV_change'] == '1/150'
    assert result['ideal_active_phase_target_equals_full_modular_state_action_on_protocol_algebra']
    # Direct dense apparatus cross-check, independent of exact read formula.
    phase, degree, read = c.native_targets()
    for state in c.prepared_states(law):
        rho = c.dense(state)
        for instrument, expected in ((phase, [c.modular_action(state, full=False)]),
                                    (degree, [c.degree_read(state, a) for a in (0, 1)]),
                                    (read, [c.y_read(state, a) for a in (0, 1)])):
            outputs = c.runtime._apply(instrument, rho, 1)
            for actual, exact in zip(outputs, expected, strict=True):
                np.testing.assert_allclose(actual, c.dense(exact), rtol=0, atol=5e-17)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_two_event_joint_law_is_continuing_and_differs_from_reset(law):
    degree, palm, _, transition = c.events.embedded()
    # Independently use the diagonal source jump chain; no sparse CP apply.
    joint = [[sum(palm[i]*transition[i][j] for i in range(66) for j in range(66)
                  if degree[i] == a and degree[j] == b) for b in (8, 10)] for a in (8, 10)]
    assert joint == [[F(1, 5), F(1, 5)], [F(1, 5), F(2, 5)]]
    assert joint[0][0]-sum(joint[0])*sum(row[0] for row in joint) == F(1, 25)
    result = c.predictions(law)
    assert result['two_event_degree_joint_probabilities'] != result['reset_control_joint_probabilities']
    for row in result['two_event_finite_cap_enclosures']:
        cap = row['attempt_cap_per_event']
        mass = F(row['conditional_unresolved_mass'])
        assert 0 <= mass <= min(F(1), 2*F(13, 15)**cap)
        assert F(row['unconditional_unresolved_mass']) == F(5, 33)*mass


@pytest.mark.parametrize('law,expected,defect', [('C', F(69, 140), F(1, 924)),
                                               ('E', F(109, 220), F(1, 1452))])
def test_observable_covariance_defect_keeps_original_herald_yield(law, expected, defect):
    result = c.predictions(law)
    assert F(result['phase_before_event_conditional_plus']) == expected
    assert F(result['phase_after_event_conditional_plus']) == F(1, 2)
    assert F(result['pooled_herald_joint_probability_defect']) == defect
    assert defect == F(5, 33)*(F(1, 2)-expected)
    # An individual coherence already exposes the noncommuting channel.
    unit = {(55, 57): c.z(1)}
    assert c.apply(law, c.modular_action(unit)) != c.modular_action(c.apply(law, unit))


@pytest.mark.parametrize('law', ['C', 'E'])
def test_two_event_runtime_retains_both_timeouts_and_failed_preparation(law):
    result = c.degree_experiment(law)
    rows = result['complete_records']
    assert len(rows) == 22
    assert sum(F(row['probability']) for row in rows) == 1
    assert rows[0]['probability'] == '28/33'
    assert rows[0]['native_attempts'] == 1
    assert rows[1]['degree_reads'] == 0
    assert any(row['record'][2] >= 0 and row['record'][4] == -1 for row in rows)
    assert result['full_state_residual_diagnostic'] < 1e-14
    assert not result['reset_between_events']
    for row in rows:
        flag, n, a, m, b = row['record']
        assert row['native_attempts'] == 1+(min(n, 2) if flag else 0)+(min(m, 2) if a >= 0 else 0)
        assert row['degree_reads'] == int(a >= 0)+int(b >= 0)


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('order', ['phase_before_event', 'phase_after_event'])
def test_covariance_runtime_retains_phase_order_resources_and_timeouts(law, order):
    result = c.covariance_experiment(law, order)
    rows = result['complete_records']
    assert len(rows) == 19
    assert sum(F(row['probability']) for row in rows) == 1
    assert rows[0]['probability'] == '28/33'
    assert rows[0]['phase_calls'] == rows[0]['balanced_mixer_calls'] == 0
    assert result['full_state_residual_diagnostic'] < 1e-14
    for row in rows[1:]:
        assert row['balanced_mixer_calls'] == 1
        if row['record'][2] == -1:
            assert row['phase_calls'] == int(order == 'phase_before_event')
            assert row['phase_reads'] == 0
            assert F(row['probability']) > 0


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_tail_preserves_arbitrary_coherence_and_uniform_trace_bound(law):
    # This state lies outside A. Tail and event maps still act on full M66;
    # source trace escape therefore also bounds arbitrary untouched references.
    state = {(0, 0): c.z(F(1, 2)), (1, 1): c.z(F(1, 2)),
             (0, 1): c.z(0, F(1, 2)), (1, 0): c.z(0, F(-1, 2))}
    for cap in (1, 3, 16):
        tail = c.apply(law, state, no_current=cap)
        stopped = c.apply(law, state, through=cap)
        assert c.trace(tail)+c.trace(stopped) == 1
        assert c.trace(tail) <= F(13, 15)**cap
        assert (0, 1) in tail
        np.testing.assert_allclose(c.dense(tail), c.dense(tail).conj().T, rtol=0, atol=0)
        assert np.linalg.eigvalsh(c.dense(tail)).min() >= -1e-16


@pytest.mark.parametrize('bad', [True, 1, 0, 'c', None, [], {}])
def test_candidate_validation_precedes_cache_dispatch(bad):
    for operation in (c.marked, c.event_columns):
        with pytest.raises(ValueError, match='candidate C or E'):
            operation(bad)


@pytest.mark.parametrize('bad', [True, 0, -1, 257, 1.5, '1'])
def test_attempt_caps_and_exclusive_controls_reject_aliases(bad):
    with pytest.raises(ValueError):
        c.apply('C', {}, through=bad)
    with pytest.raises(ValueError):
        c.degree_records(bad)
    with pytest.raises(ValueError):
        c.covariance_records(bad)


def test_exact_boundary_rejects_off_domain_or_ambiguous_inputs():
    for state in ({(True, 1): c.z(1)}, {(66, 1): c.z(1)}, {(1, 1): (1, 0)}, []):
        with pytest.raises(ValueError):
            c.sparse(state)
    for quarter in (True, 5, '1', 1.0):
        with pytest.raises(ValueError):
            c.modular_action({}, quarter)
    with pytest.raises(ValueError):
        c.modular_action({}, full=1)
    with pytest.raises(ValueError):
        c.apply('C', {}, through=2, no_current=2)
    with pytest.raises(ValueError):
        c.y_read({(0, 1): c.z(1)}, 0)
    with pytest.raises(ValueError):
        c.finite_preparation('C', coherent=1)
    with pytest.raises(ValueError):
        c.covariance_experiment('C', 'native_modular_evolution')


def test_receipt_exact_quantities_cannot_hide_behind_float_tolerance():
    correct = {'probability': '1/5', 'residual_diagnostic': 1e-17}
    c.compare_receipt(correct, {**correct, 'residual_diagnostic': 2e-17})
    with pytest.raises(ValueError):
        c.compare_receipt(correct, {**correct, 'probability': '1/4'})

"""Complete-source closure, exact contraction and multi-time continuation controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
import sympy as s

from . import repair_continuation as r, parameter_selection as selection, record_code as code


@pytest.mark.parametrize('count', [1, 4, 5, 6])
def test_complete_existing_gray_mixer_not_only_its_prepared_column(count):
    target = r._mixer_exact(count)
    assert (target.T*target-s.eye(8)).applyfunc(s.simplify) == s.zeros(8)
    change = r._mixer_exact(6)*target.T
    assert change[6:8, :count] == s.zeros(2, count)
    # The basis completion is material: simply keeping the first count
    # syndrome coordinates would be an incorrect explanation of closure.
    if count in (5, 6):
        assert target[count:, :count] != s.zeros(8-count, count)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_all_record_kraus_maps_are_tp_and_fix_arbitrary_code_operators(law):
    model = r._model(law)
    assert len(model['records']) == (160 if law == 'C' else 320)
    assert model['source_integer_input_checks'] == 20480
    assert sum(row['zero'] for row in model['records']) == (80 if law == 'C' else 200)
    for logical in range(6):
        gram = sum((row[logical].T@row[logical] for row in model['maps']), np.zeros((6, 6)))
        np.testing.assert_allclose(gram, np.eye(6), atol=2e-15, rtol=0)
    target = np.ones((6, 6))/6
    for left in range(6):
        for right in range(6):
            np.testing.assert_allclose((r._pair_transfer(law, left, right)@target.ravel()).reshape(6, 6),
                                       target, atol=3e-16, rtol=0)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_exact_block_contraction_matches_independent_full512_completion(law):
    certificate = r.contraction_certificate(law)
    assert certificate['outward_rational_arithmetic']
    assert max(certificate['directed_Frobenius_upper_diagnostics']) < float(F(certificate['strict_Frobenius_upper_exact']))
    V = code.ideal_code_matrix()
    Q = np.eye(512)-V@V.T.conj()
    full = Q@r.source_cross_transfer(law)@Q
    from ..observer.source_discrimination import circuit
    for logical, support in enumerate(code.SUPPORTS):
        indices = [circuit.encode(state) for state in support]
        np.testing.assert_allclose(full[np.ix_(indices, indices)],
                                   np.asarray(r._model(law)['cross'][logical], float), atol=8e-16, rtol=0)
        outside = full[:, indices].copy(); outside[indices] = 0
        assert np.linalg.norm(outside, 2) < 1e-15


@pytest.mark.parametrize('law', ['C', 'E'])
def test_dressed_poisson_identity_closes_the_full_off_code_first_derivative(law):
    model = r._model(law)
    geometry = selection.code_geometry()
    energies = selection.diagonal_cost(2)['values']
    u = np.ones(6)/math.sqrt(6)
    projector = np.eye(6)-np.outer(u, u)
    H = [np.diag([float(energies[i]) for i in row]) for row in geometry['supports']]
    h = [float(value) for value in r.logical_hamiltonian().diagonal()]
    Y = []
    for j in range(6):
        M = np.asarray(model['cross'][j], float)
        K = projector@H[j]@u
        Y.append(np.linalg.solve(np.eye(6)-M, M@K))
    for left in range(6):
        for right in range(6):
            D = -1j*(np.outer(Y[left], u)-np.outer(u, Y[right]))
            LJ = -1j*(H[left]@np.outer(u, u)-np.outer(u, u)@H[right])
            PhiLJ = (r._pair_transfer(law, left, right)@LJ.ravel()).reshape(6, 6)
            PhiD = (r._pair_transfer(law, left, right)@D.ravel()).reshape(6, 6)
            expected = PhiLJ+1j*(h[left]-h[right])*np.outer(u, u)
            np.testing.assert_allclose(D-PhiD, expected, atol=2e-16, rtol=0)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_finite_composition_with_reference_converges_and_matches_universal_bound(law):
    errors = []
    for steps in (256, 1024, 4096):
        time = F(256)
        state = r._choi(r.continuation_blocks(law, time, steps))
        assert abs(np.trace(state)-1) < 2e-11
        assert np.linalg.eigvalsh(state).min() > -2e-12
        error = r._trace_norm(state-r._ideal_choi(time))
        assert error < float(F(r.error_bound(law, time, steps)['trace_norm_upper_exact']))+2e-11
        errors.append(error)
    assert errors[2] < errors[1]/3 < errors[0]/9


def test_exact_macroscopic_bound_requires_per_interval_service_error_accounting():
    for law, gate in (('C', F(1, 512)), ('E', F(1, 1024))):
        ideal = r.error_bound(law, F(1024), 1 << 24)
        assert F(ideal['trace_norm_upper_exact']) < gate
        amplified = r.implementation_bound(law, F(1024), 1 << 24, F(1, 10**6))
        assert amplified['combined_trace_norm_upper_exact'] == '2'
        assert not amplified['literal_source_diagnostics_or_physical_errors_certify_this_allowance']


@pytest.mark.parametrize('law', ['C', 'E'])
def test_public_contraction_receipt_mutation_cannot_weaken_a_warm_cached_bound(law):
    before = r.error_bound(law, F(1024), 65536)
    public = r.contraction_certificate(law)
    public['sum_norm_powers_from_one_upper_exact'] = '0'
    public['directed_Frobenius_upper_diagnostics'][0] = 0.
    public['strict_Frobenius_upper_exact'] = '0'
    after = r.error_bound(law, F(1024), 65536)
    assert after == before
    fresh = r.contraction_certificate(law)
    assert fresh['sum_norm_powers_from_one_upper_exact'] == ('15' if law == 'C' else '5')
    assert fresh['directed_Frobenius_upper_diagnostics'][0] > 0


@pytest.mark.parametrize('law', ['C', 'E'])
def test_finite_rate_corrections_have_damping_and_approach_the_derived_gaps(law):
    result = r.finite_rate_modes(law)
    limits = (1/256, -1/1024, -1/3072)
    errors = []
    for row in result['controls']:
        error = 0.
        for channel, gap in zip(row['channels'], limits, strict=True):
            assert channel['attenuation_rate'] > 0
            error = max(error, abs(channel['wrapped_phase_rate']-gap))
        errors.append(error)
    assert errors[0] < 1e-5
    assert errors[0] < errors[1]/5 < errors[2]/25
    assert not result['wrapped_phase_rates_are_unwrapped_physical_masses']


@pytest.mark.parametrize('law', ['C', 'E'])
def test_two_time_complete_instruments_include_outside_syndrome_outcomes(law):
    result = r.two_time_control(law)
    probabilities = np.array(result['joint_probabilities'])
    assert abs(probabilities.sum()-1) < 3e-11
    assert probabilities.min() > -2e-12
    assert probabilities[2].sum()+probabilities[:, 2].sum() > 1e-8
    assert result['L1_difference'] < float(F(result['universal_L1_upper_exact']))
    assert result['L1_difference'] < 0.004


@pytest.mark.parametrize('law', ['C', 'E'])
def test_old_terminal_timeout_control_freezes_instead_of_claiming_continuation(law):
    first = r.absorbing_control(law, F(16), 16)
    later = r.absorbing_control(law, F(16), 256)
    assert F(first['survival_probability_exact']) == F(27, 32)**16
    assert F(later['survival_probability_exact']) == F(27, 32)**256
    assert later['Choi_trace_norm_to_unchanged_input_diagnostic'] < first['Choi_trace_norm_to_unchanged_input_diagnostic']/5
    assert later['Choi_trace_norm_to_unchanged_input_diagnostic'] < float(F(later['trace_norm_to_unchanged_input_upper_exact']))
    assert not later['timeout_restarted']


def test_clock_winding_has_same_endpoint_and_different_intermediate_record():
    result = r.clock_winding_control()
    endpoint, middle = result['intermediate_quadrature_records_distinguish_it']
    assert endpoint['X_mean'] == pytest.approx(1., abs=1e-14)
    assert middle['X_mean'] == pytest.approx(-1., abs=1e-14)
    assert not result['an_identity_endpoint_selects_a_unique_unwrapped_generator']


@pytest.fixture(scope='module')
def executed():
    return r.source_execution()


def test_literal_perturbed_corrected_source_retains_all_rows_and_matches_actual_completion(executed):
    report, arrays = executed
    assert max(report['retained_source_operator_errors'].values()) < 2e-10
    assert min(report['actual512_versus_alternative66_output_operator_difference'].values()) > 1e-7
    assert arrays['native_C_outputs'].shape == (8, 512, 42)
    assert arrays['native_E_outputs'].shape == (16, 512, 42)
    assert not report['full_N_interval_pulse_or_protection_history_executed']


def test_retained_source_arrays_must_meet_independent_targets(executed):
    _, arrays = executed
    changed = dict(arrays)
    changed['mixer_5'] = np.eye(8, dtype=complex)
    with pytest.raises(ValueError, match='admission'):
        r.verify_source_arrays(changed)
    changed = dict(arrays)
    changed['phase_outputs'] = arrays['phase_outputs'].copy()
    changed['phase_outputs'][511, 0] += 1e-4
    with pytest.raises(ValueError, match='admission'):
        r.verify_source_arrays(changed)


@pytest.mark.parametrize('time,steps', [(True, 16), (1., 16), (F(-1), 16), (F(1), True), (F(1), 0), (F(1), 1 << 25)])
def test_parameter_domain_refuses_ambiguous_clock_and_schedule(time, steps):
    with pytest.raises(ValueError):
        r.continuation_blocks('C', time, steps)


def test_exact_contract_types_and_admissions_cannot_drift():
    with pytest.raises(ValueError):
        r._compare({'numerical_source_admission': 1e-8}, {'numerical_source_admission': 1.01e-8})
    with pytest.raises(ValueError):
        r._compare({'derived': False}, {'derived': 0})
    with pytest.raises(ValueError):
        r.error_bound('bad', F(1), 16)
    with pytest.raises(ValueError):
        r.implementation_bound('C', F(1), 16, 1e-9)

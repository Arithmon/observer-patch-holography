"""Common native continuum generator, time parameters and wave-path controls."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from model.observer.scalar_response.experiment import compare_receipt
from . import continuum_dynamics as d


@pytest.mark.parametrize('level', (0, 1, 3, 8))
def test_native_integer_step_and_coherent_clock_errors_use_actual_normalizations(level):
    bounds = d.propagation_bounds(level, band=16, horizon=4)
    q, n = 630*16**level, 1260*16**level
    for a in (F(1, 10), F(1), F(7), F(16)):
        for parameter in (F(1, 13), F(7, 10), F(4)):
            with mp.workdps(60):
                a0, t = mp.mpf(a.numerator)/a.denominator, mp.mpf(parameter.numerator)/parameter.denominator
                attempts = int(mp.floor(t*n))
                heat_error = abs((1-a0/n)**attempts-mp.exp(-t*a0))
                steps = int(mp.floor(t*mp.sqrt(q)))
                theta = 2*mp.asin(mp.sqrt(a0/(4*q)))
                wave_error = abs(mp.cos(steps*theta)-mp.cos(t*mp.sqrt(a0)))
            assert heat_error <= float(bounds['integer_attempt_heat_error_on_band'])
            assert wave_error <= float(bounds['integer_coherent_wave_error_on_band'])
    assert bounds['spatial_subsequence_convergence_rate_supplied'] is False


@pytest.mark.parametrize('level', (0, 2, 5, 10))
def test_actual_positive_length_polynomial_freezes_at_fixed_cap(level):
    cap, band = 8, F(16)
    value = d.fixed_cap_collapse(level, cap, band)
    coefficients = d.wave_reference.coefficients(cap)
    with mp.workdps(80):
        q = 630*mp.mpf(16)**level
        x = 1-mp.mpf(16)/(2*q)
        weights = tuple(mp.mpf(c.numerator)/c.denominator for c in coefficients)
        theta = -mp.log(sum(c*x**i for i, c in enumerate(weights))/sum(weights))
        scaled = mp.sqrt(q)*theta
        ideal = mp.sqrt(q)*mp.acos(x)
    assert 0 < scaled <= float(value['scaled_modular_frequency_upper'])
    assert ideal >= 4-mp.mpf('1e-60')
    if level >= 5:
        assert scaled < .01


def test_growing_cap_budget_has_finite_scaled_error_and_distinct_clock_parameter():
    old = None
    for level in range(6):
        value = d.propagation_bounds(level)
        assert value['sufficient_native_return_cap'] == 64*256**level
        if old is not None:
            assert value['joint_cap_modular_cosine_error_upper'] <= old
        old = value['joint_cap_modular_cosine_error_upper']
    # The passive attempt and coherent step units diverge under refinement.
    assert F(1260*16**5, 1260*16**4) == 16
    assert F(4**5, 4**4) == 4


@pytest.mark.parametrize('level', (0, 1, 4, 12))
def test_weighted_native_heat_bound_has_correct_continuum_constant(level):
    for kappa in (F(1, 3), F(1), F(2)):
        value = d.finite_speed_bounds(level, kappa, separation=3, heat_parameter=2)
        with mp.workdps(80):
            k = mp.mpf(kappa.numerator)/kappa.denominator
            exact = 2*4**level*(mp.cosh(k/2**level)-1)
            upper = value['native_weighted_heat_log_norm_upper']
            assert exact <= mp.mpf(upper.numerator)/upper.denominator
            assert exact >= k*k
        assert value['continuum_weighted_heat_log_norm_upper'] == kappa*kappa
        assert value['continuum_separated_heat_log_norm_upper'] == -F(9, 4)
        assert value['continuum_zero_wave_matrix_element_when_time_squared_below'] == 18
        assert value['continuum_wave_speed_squared_upper'] == F(1, 2)
        assert value['constant_zero_mode_restored_before_support_localization']


def test_weighted_form_uses_actual_full_native_edges_and_rejects_internal_weights():
    value = d.native_weighted_form_control()
    assert value['exact_weighted_form_decomposition']
    assert F(value['maximum_port_penalty']) == F(value['maximum_cross_seam_cosh_minus_one'])
    assert F(value['internal_weight_mutation_cross_seam_penalty']) == 0
    assert F(value['internal_weight_mutation_constant_quadratic_form']) < 0
    assert value['finite_continuous_wave_has_instantaneous_small_tail']
    assert F(value['finite_cross_seam_wave_t_squared_coefficient']) == F(1, 2)


def test_davies_gaffney_optimization_and_wave_distance_normalization():
    k, r, tau = sp.symbols('k r tau', positive=True)
    exponent = tau*k*k/2-k*r
    assert sp.diff(exponent, k).subs(k, r/tau) == 0
    assert exponent.subs(k, r/tau) == -r*r/(2*tau)
    assert -(sp.sqrt(2)*r)**2/(4*tau) == -r*r/(2*tau)
    # A mean-zero projection spreads a localized input; support claims need
    # the full operator with the constant zero mode restored.
    localized = sp.Matrix([1, 0, 0, 0])
    projected = (sp.eye(4)-sp.ones(4)/4)*localized
    assert all(projected[i] != 0 for i in range(4))


def test_weighted_heat_geometric_series_requires_a_sufficiently_fine_level():
    with pytest.raises(ValueError):
        d.finite_speed_bounds(0, lipschitz=4)
    assert d.finite_speed_bounds(1, lipschitz=4)['native_weighted_heat_log_norm_upper'] == 12


def test_native_stirring_and_random_canonical_covariance_have_independent_oracles():
    value = d.native_stirring_control()
    assert value['canonical_configurations'] == 495
    assert value['native_endpoint_branches'] == 29700
    assert value['independent_two_attempt_stirring_branches'] == 3600
    assert value['native_unread_repair_equals_identity_or_swap']
    assert value['conditional_covariance_varies_before_averaging']
    assert value['conditional_canonical_covariance_verified_on_every_branch']
    assert value['stirring_and_full_native_two_time_covariance_exactly_equal']
    assert 0 < value['two_attempt_overlap_variance'] <= value['uniform_overlap_variance_upper']
    assert value['two_time_sum_variance'] == 2+2*value['two_attempt_mean_overlap']


def test_reusing_a_seam_keeps_stationarity_but_changes_multitime_density_law():
    value = d.native_stirring_control()
    assert value['repeated_seam_mutation_mean_overlap'] != value['two_attempt_mean_overlap']
    states, _, edges = d._native_trace_counts(4)
    samples = np.asarray(states, dtype=np.int64)
    # Every individual seam permutation preserves the complete uniform sector.
    for a, b in edges:
        order = list(range(12))
        order[a], order[b] = order[b], order[a]
        assert set(map(tuple, samples[:, order])) == set(states)
    assert value['single_time_unit_variance_fourth_moment'] != 3
    assert not value['finite_native_snapshot_is_exactly_Gaussian']


@pytest.mark.parametrize('level', (0, 2, 5, 10))
def test_actual_density_process_uses_simultaneous_diffusive_attempt_scaling(level):
    value = d.density_process_bounds(level)
    v, q = 240*4**level, 80*4**level
    assert value['ports'] == v and value['particles'] == q
    assert value['native_attempts_per_comparison_time'] == 2*(630*4**level)*4**level
    assert value['susceptibility'] == F(q*(v-q), v*(v-1))
    assert value['stirring_overlap_variance_upper'] == F(4, v)
    assert value['overlap_bound_uniform_in_native_attempt_count']
    assert value['conditional_Hajek_ratio_upper_under_variance_assumption'] == F(18, q)
    assert value['simultaneous_refinement_and_attempt_time_finite_dimensional_limit']
    assert not value['count_first_long_block_limit_required']
    assert not value['path_space_tightness_asserted']
    assert not value['clock_is_laboratory_proper_time']


def test_stirring_variance_telescopes_only_for_a_positive_contraction():
    x = sp.symbols('x', real=True)
    # The discarded remainder is nonnegative exactly on the positive band.
    assert sp.expand((1-x*x)-(1-x)-x*(1-x)) == 0
    for value in (F(0), F(1, 3), F(2, 3), F(1)):
        for n in (1, 3, 20):
            cost = sum(value**(2*k)*(1-value) for k in range(n))
            assert cost <= 1-value**(2*n)
    assert 1-F(-1, 2) > 1-F(-1, 2)**2


def test_reflected_ou_covariance_reconstructs_decay_generator_not_wave_frequency():
    # An exact positive transfer tests the reflected covariance factorization.
    transfer = sp.diag(sp.Rational(1, 2), sp.Rational(1, 3))
    fields = [sp.Matrix([1, 0]), sp.Matrix([1, 2]), sp.Matrix([-1, 1])]
    times = (1, 2, 4)
    gram = sp.Matrix(3, 3, lambda i, j:
        (fields[i].T*transfer**(times[i]+times[j])*fields[j])[0])
    columns = sp.Matrix.hstack(*(transfer**t*f for t, f in zip(times, fields)))
    assert gram == columns.T*columns
    for value in (sp.Matrix([1, 0, -1]), sp.Matrix([1, -2, 3])):
        assert (value.T*gram*value)[0] >= 0
    # Spatial filtering changes the covariance, not the eigenmode decay rate.
    eigenvalues = np.array([1., 4., 9.])
    tau = .5
    actual = np.exp(-tau*eigenvalues)/eigenvalues
    relativistic = np.exp(-tau*np.sqrt(eigenvalues))/eigenvalues
    assert actual[0] == relativistic[0]
    assert np.all(actual[1:] < relativistic[1:])


@pytest.mark.parametrize('bad', (True, 0, -1, 1., '1'))
def test_density_limit_requires_an_original_positive_test_count(bad):
    with pytest.raises(ValueError):
        d.density_process_bounds(0, tests=bad)


def test_native_green_wave_covariance_has_identity_velocity_not_quantum_vacuum():
    value = d.exact_wave_synthesis()
    assert value['constant_mode_exactly_removed']
    _, _, _, _, _, source = d.dynamics._base_blocks()
    eigen, vectors = np.linalg.eigh(np.asarray(source, dtype=float))
    positive = eigen[1:]
    covariance = 1/positive-1/2520
    assert np.all(covariance > 0)
    np.testing.assert_allclose(positive*covariance, 1-positive/2520, atol=1e-15)
    assert np.linalg.norm(1/positive-1/(2*np.sqrt(positive))) > 1
    for t, u in ((0., .3), (.5, 1.2), (2., 2.)):
        actual = covariance*(np.cos(t*np.sqrt(positive))*np.cos(u*np.sqrt(positive))+
                             np.sin(t*np.sqrt(positive))*np.sin(u*np.sqrt(positive)))
        expected = covariance*np.cos((t-u)*np.sqrt(positive))
        np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=2e-14)


def test_reusing_one_gaussian_block_breaks_stationary_wave_covariance():
    # Two disjoint CLT increments are essential for the declared Gaussian law.
    a, t = 3., .4
    one_block = (np.cos(t*np.sqrt(a))+np.sin(t*np.sqrt(a)))**2/a
    two_blocks = (np.cos(t*np.sqrt(a))**2+np.sin(t*np.sqrt(a))**2)/a
    assert abs(two_blocks-1/a) < 1e-15
    assert abs(one_block-1/a) > .1


def test_weak_compressed_static_data_alone_do_not_determine_the_resolvent():
    # In a growing basis the second coordinate can escape every fixed dictionary.
    # G_n=I+(|1><n|+|n><1|)/2 converges weakly to I, but not strongly.
    block = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    assert block[0, 0] == 1
    actual = (block*(sp.eye(2)+block).inv())[0, 0]
    assert actual == sp.Rational(7, 15)
    assert actual != sp.Rational(1, 2)
    # Source compactness excludes the escaping off-diagonal mass in this control.
    assert (block-sp.eye(2))[:, 0].dot((block-sp.eye(2))[:, 0]) == sp.Rational(1, 4)


def test_actual_full_native_coherent_oracle_keeps_both_Ramsey_outcomes():
    value = d.native_propagator_control()
    assert value['source_ports'] == 240
    assert value['full_native_seams'] == 630
    rows = value['readouts']
    for row in rows:
        assert abs(row['Ramsey_plus_diagnostic']+row['Ramsey_minus_diagnostic']-1) < 1e-15
        assert row['spectral_oracle_residual_diagnostic'] < 3e-12
    assert rows[0]['Ramsey_plus_diagnostic'] == 1
    assert rows[-1]['Ramsey_minus_diagnostic'] > .01


def test_wave_path_tail_controls_the_supremum_and_tends_to_zero():
    for cutoff in (0, 10, 100, 1000):
        value = d.wave_path_bounds(cutoff, F(3, 2))
        expected = 4*(d.spatial.GREEN_BOUND+F(9, 4))/(cutoff+1)**2
        assert value['expected_squared_supremum_tail_upper'] == expected
        assert value['expected_fourth_increment_upper_per_dt_fourth'] == 12
    assert d.wave_path_bounds(1000)['expected_squared_supremum_tail_upper'] < F(1, 100)


@pytest.mark.parametrize('particles', (2, 4))
def test_complete_native_branch_trace_matches_permutation_formula(particles):
    value = d.native_reference_state_control(particles)
    assert value['native_and_permutation_trace_moments_exactly_equal']
    assert value['configuration_dimension'] == (66 if particles == 2 else 495)
    assert value['native_endpoint_branches'] == 60*value['configuration_dimension']
    assert value['cap_one_commutator_expectation_minus_i'] == (F(1, 366) if particles == 2 else F(2, 435))
    assert value['ideal_reference_half_trace_distance_diagnostic'] > 0
    assert not value['finite_reference_exactly_tracial']
    # A different supplied nontracial state retains a nonzero commutator.
    assert value['different_supplied_two_configuration_pure_state_commutator'] == 1


def test_permutation_cycle_formula_has_an_independent_subset_oracle():
    for particles in (2, 4):
        states, _, edges = d._native_trace_counts(particles)
        moments = d._permutation_moments(12, particles, 30, 240)
        for cycle, key in (((1, 0, *range(2, 12)), 'swap_trace'),
                           ((1, 2, 0, *range(3, 12)), 'three_cycle_trace'),
                           ((1, 0, 3, 2, *range(4, 12)), 'disjoint_swap_trace')):
            fixed = sum(tuple(state[i] for i in cycle) == state for state in states)
            assert F(fixed, len(states)) == moments[key]
    assert d._native_trace_counts(2)[1] == tuple(tuple(map(int, row))
        for row in d.dynamics.native_density()['counts'].tolist())


@pytest.mark.parametrize('level', (0, 1, 2))
def test_actual_source_seams_satisfy_all_level_trace_census(level):
    value = d.source_trace_census(level)
    assert value['source_census_matches_all_level_formula']
    assert value['degree_five_ports']+value['degree_six_ports'] == 240*4**level


def test_trace_variance_bound_is_an_all_level_rational_inequality():
    v = sp.symbols('v', positive=True)
    exact = 2*(365*v*v-1631*v+1560)/(11907*(v-2)*(v-1)**2)
    difference = sp.factor(1/(16*v)-exact)
    expected = (227*v**3+4564*v*v+9615*v-23814)/(190512*v*(v-2)*(v-1)**2)
    assert sp.simplify(difference-expected) == 0
    # Every numerator term is positive for v>=3 after this constant grouping.
    assert 9615*3-23814 > 0
    for level in (0, 1, 4, 10, 100):
        value = d.reference_state_bounds(level)
        direct = d._permutation_moments(240*4**level, 80*4**level,
                                        630*4**level, 5400*4**level)
        assert value['transfer_trace_variance'] == direct['transfer_trace_variance']
        assert value['ideal_reference_half_trace_distance_to_normalized_trace_upper'] == F(3, 8*2**level)


def test_density_mode_weight_bound_charges_full_configuration_normalization():
    import math
    for level in (0, 1):
        value = d.reference_state_bounds(level)
        v, q = value['ports'], value['particles']
        exact_upper = F(5*(v-1), math.comb(v, q))
        exponent = value['density_mode_reference_probability_log2_upper']
        assert exact_upper <= F(1, 2**(-exponent))
        assert exponent < -60
    assert d.reference_state_bounds(4)['density_mode_reference_probability_log2_upper'] < -20000


def test_fixed_parameter_modular_bound_does_not_bound_rescaled_wave_time():
    fixed = d.reference_state_bounds(12)
    rescaled = d.reference_state_bounds(12, modular_parameter=26*4**12)
    assert fixed['ideal_modular_automorphism_normalized_trace_L2_error_upper'] < F(1, 10000)
    assert rescaled['ideal_modular_automorphism_normalized_trace_L2_error_upper'] == 2
    assert not fixed['modular_operator_norm_triviality_claimed']
    assert not fixed['refinement_rescaled_modular_time_triviality_claimed']


def test_local_fluctuation_commutator_bound_accounts_for_overlapping_terms():
    value = d.reference_state_bounds(10, terms_per_vertex=3, overlaps=5, norm=2)
    assert value['bounded_overlap_fluctuation_commutator_norm_upper'] == 120
    epsilon = value['growing_cap_reference_half_trace_distance_to_normalized_trace_upper']
    assert value['growing_cap_fluctuation_commutator_expectation_upper'] == 240*epsilon
    assert not value['every_noncommutative_or_nontracial_limit_excluded']


def test_cached_native_trace_rows_cannot_be_mutated():
    _, rows, _ = d._native_trace_counts(2)
    with pytest.raises(TypeError):
        rows[0][0] = 0


@pytest.mark.parametrize('bad', (True, -1, 1., '1'))
def test_ambiguous_levels_caps_and_cutoffs_rejected(bad):
    with pytest.raises(ValueError):
        d.propagation_bounds(bad)
    with pytest.raises(ValueError):
        d.fixed_cap_collapse(0, bad)
    with pytest.raises(ValueError):
        d.wave_path_bounds(bad)


def test_retained_continuum_receipt_replays_and_rejects_physical_promotion():
    result = d.run()
    path = Path(__file__).with_name('continuum_dynamics_results.json')
    compare_receipt(json.loads(path.read_text()), result)
    changed = copy.deepcopy(result)
    changed['scope']['observed_one_plus_three_spacetime_identified'] = True
    with pytest.raises(ValueError):
        compare_receipt(changed, result)
    changed = copy.deepcopy(result)
    changed['scope']['density_process_path_space_limit_proved'] = True
    with pytest.raises(ValueError):
        compare_receipt(changed, result)


def test_fresh_runtime_source_closure_is_pinned():
    command = '''
import sys
from pathlib import Path
from model.observer_time_fields import continuum_dynamics as d
d.run()
missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path is None:continue
    try:name=Path(path).resolve().relative_to(d.ROOT).as_posix()
    except ValueError:continue
    if name.endswith('.py') and name not in d.SOURCE_PATHS:missing.append(name)
assert not missing,missing
'''
    subprocess.run([sys.executable, '-c', command], cwd=d.ROOT, check=True)

"""Complete native moments, coupled clocks and continuous-path controls."""
from fractions import Fraction as F
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

from . import observation_process as o


@pytest.mark.parametrize('points', (1, 3, 5))
def test_exact_product_dictionary_diagonalizes_declared_radial_energy(points):
    spacing = sp.Rational(2, points+1)
    matrix = sp.diag(*([2]*points))
    for j in range(points-1):
        matrix[j, j+1] = matrix[j+1, j] = -1
    for k in range(1, points+1):
        row = o.product_mode_dictionary(6, k, points)
        vector = sp.Matrix([sp.sin(sp.pi*k*j/(points+1)) for j in range(1, points+1)])
        assert all(sp.simplify(value) == 0 for value in (
            matrix*vector/spacing**2-row['radial_eigenvalue']*vector))
        assert row['product_eigenvalue'] == row['radial_eigenvalue']+6
        assert row['native_eigenvector_membership_requires_a_separate_certificate']


@pytest.mark.parametrize('beta', (F(1, 16), F(1, 2), F(1), F(2), F(4)))
def test_record_filter_distinguishes_classical_and_quantum_product_covariances(beta):
    row = o.product_thermal_control(beta)
    for key in ('direct_product_matrix_diagonalization_error',
                'sampled_bridge_DST_covariance_error',
                'derived_record_filter_covariance_error',
                'ordered_covariance_commutator_identity_error',
                'ordered_covariance_KMS_identity_error'):
        assert row[key] < 1e-10
    assert row['sampled_bridge_vs_grid_inverse_discrepancy'] > .001
    assert row['classical_inverse_vs_quantum_thermal_covariance_discrepancy'] > .001
    assert row['smallest_mode_uncertainty_product'] >= .25
    assert row['omitted_angular_constant_disjoint_readout_response_derivative'] == -.05
    assert max(row['source_shifted_product_resolvent_errors']) < .125
    assert np.shape(row['all_thermal_radial_readouts_from_complete_96_event_tape']) == (3, 20)
    assert not row['classical_record_filter_establishes_quantum_CCR']
    assert not row['finite_tape_certifies_Gaussian_thermal_ensemble']


def test_retained_product_covariance_and_bound_share_inverse_temperature():
    assert o.product_thermal_control()['inverse_temperature'] == '1'
    assert o.product_thermal_bounds()['inverse_temperature'] == F(1)


@pytest.mark.parametrize('points', (3, 8, 16))
def test_radial_frequency_and_thermal_errors_obey_the_band_bound(points):
    beta, length, horizon = F(3, 2), F(5, 2), F(2)
    bounds = o.product_thermal_bounds(beta, length, points, 2, horizon)
    for k in (1, 2):
        row = o.product_mode_dictionary(0, k, points, length)
        rho, continuous = map(float, (row['radial_eigenvalue'], row['continuum_radial_eigenvalue']))
        assert 0 <= continuous-rho <= float(bounds['radial_squared_frequency_error_upper'])
        for lam in (0., .01, 1., 100.):
            frequencies = np.sqrt(np.array([rho, continuous])+lam)
            thermal = 1/(2*frequencies*np.tanh(float(beta)*frequencies/2))
            assert abs(np.diff(thermal)[0]) <= float(bounds['thermal_covariance_grid_error_upper'])
            for time in (0., .5, 2.):
                response = np.sin(time*frequencies)/frequencies
                ordered = thermal*np.cos(time*frequencies)-.5j*response
                assert abs(np.diff(response)[0]) <= float(bounds['sine_response_grid_error_upper'])
                assert abs(np.diff(ordered)[0]) <= float(bounds['ordered_thermal_covariance_grid_error_upper'])


def test_finite_count_covariance_budget_follows_the_actual_source_gap():
    graph = o.spatial.source.source.federation(0)
    laplacian = np.zeros((graph.ports, graph.ports))
    for a, b in graph.seams:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    eigenvalues = np.linalg.eigvalsh(laplacian)[1:]
    assert eigenvalues[0] >= 1/4097
    density = 1-eigenvalues/(2*len(graph.seams))
    for count in (2, 17, 100):
        bound = o.native_block_covariance_bound(0, count, 6, F(3))
        # Exact finite-density eigenmode expressions, independently summed.
        cross = sum(min(i+1, 2*count-1-i)*density**(i+1)
                    for i in range(2*count-1))/(4*len(graph.seams)*count)
        one = (1+2*sum((1-j/count)*density**j for j in range(1, count)))/(4*len(graph.seams))
        limiting = (1+density)/(1-density)/(4*len(graph.seams))
        assert max(cross) <= float(bound['adjacent_block_cross_covariance_upper'])
        assert max(abs(one-limiting)) <= float(bound['one_block_covariance_defect_upper'])
        assert bound['filtered_joint_covariance_defect_upper'] == 36*bound['adjacent_block_cross_covariance_upper']
        assert not bound['finite_count_Gaussian_total_variation_bound']


def test_product_thermal_control_detects_a_changed_radial_mode_map(monkeypatch):
    original = o.product_mode_dictionary
    def changed(*args, **kwargs):
        row = original(*args, **kwargs)
        row['radial_eigenvalue'] *= sp.Rational(11, 10)
        return row
    monkeypatch.setattr(o, 'product_mode_dictionary', changed)
    with pytest.raises(ArithmeticError, match='product thermal'):
        o.product_thermal_control()


@pytest.mark.parametrize('eigenvalue', (True, 1., sp.Float(6), '6', -1, sp.oo, sp.Symbol('lambda'), sp.I))
def test_product_mode_dictionary_requires_exact_certifiable_eigenvalues(eigenvalue):
    with pytest.raises(ValueError):
        o.product_mode_dictionary(eigenvalue)


def test_product_contract_inputs_fail_closed():
    for kwargs in ({'radial_mode': 0}, {'radial_mode': 4}, {'radial_points': True},
                   {'length': 2.}, {'angular_scale': F(0)}):
        with pytest.raises(ValueError): o.product_mode_dictionary(6, **kwargs)
    for kwargs in ({'radial_band': 4}, {'radial_points': 0}, {'horizon': F(0)},
                   {'inverse_temperature': 2.}):
        with pytest.raises(ValueError): o.product_thermal_bounds(**kwargs)
    with pytest.raises(ValueError): o.product_thermal_control(F(5))
    with pytest.raises(ValueError): o.native_block_covariance_bound(0, 0, 1)


@pytest.mark.parametrize('level', (0, 1))
def test_native_innovations_supply_full_rank_radial_markov_gluing(level):
    row = o.radial_record_gluing(level)
    assert row['joint_covariance_rank'] == row['expected_joint_covariance_rank'] == 60*4**level
    assert row['smallest_conditional_innovation_eigenvalue'] > 0
    assert row['native_filter_covariance_error'] < 1e-10
    assert row['joint_cross_radius_covariance_error'] < 1e-10
    assert row['joint_precision_identity_error'] < 1e-10
    assert row['nonneighbor_radial_precision_norm'] == 0
    assert row['conditional_past_future_covariance_norm'] < 1e-10
    assert row['inserted_layer_gluing_error'] < 1e-10
    assert row['reused_block_covariance_discrepancy'] > 0.01
    assert row['independent_CLT_blocks_of_one_history_required'] == 6
    assert row['zero_angular_mode_is_independent_encoded_readout']
    executed = row['finite_native_record_execution']
    assert len(executed['complete_signed_current_records_including_zero']) == 96
    assert 0 in executed['complete_signed_current_records_including_zero']
    assert np.shape(executed['all_radial_slice_readouts']) == (3, 20*4**level)
    assert executed['decoded_field_checked_after_every_native_event']
    assert not executed['finite_tape_is_a_Gaussian_ensemble_certificate']
    assert not row['native_conserved_total_load_fluctuates']
    assert not row['Gaussian_independence_holds_at_finite_native_count']
    assert not row['physical_radial_transport_or_Lorentzian_clock_identified']


@pytest.mark.parametrize('eigenvalue', (0., 5-np.sqrt(5), 6., 5+np.sqrt(5)))
def test_radial_covariance_inverse_matches_independent_dirichlet_boundary_energy(eigenvalue):
    times, length = np.array([0.2, 0.7, 1.4]), 2.0
    covariance = np.array([[o._radial_dirichlet_covariance(
        np.array([eigenvalue]), r, s, length)[0] for s in times] for r in times])
    gaps = np.diff(np.concatenate(([0.], times, [length])))
    root = np.sqrt(eigenvalue)
    diagonal = 1/gaps if root == 0 else root/np.tanh(root*gaps)
    crossing = 1/gaps if root == 0 else root/np.sinh(root*gaps)
    precision = np.diag(diagonal[:-1]+diagonal[1:])
    for i in range(len(times)-1):
        precision[i, i+1] = precision[i+1, i] = -crossing[i+1]
    np.testing.assert_allclose(precision@covariance, np.eye(3), atol=2e-14, rtol=0)
    # A deterministic harmonic extension of one boundary value has no innovation.
    boundary = np.exp(-root*times)
    harmonic_covariance = np.outer(boundary, boundary)
    assert np.linalg.matrix_rank(harmonic_covariance) == 1
    assert np.linalg.matrix_rank(covariance) == 3
    conditional_variance = covariance[1, 1]-covariance[1, 0]**2/covariance[0, 0]
    assert conditional_variance > 0


def test_finite_native_blocks_do_not_receive_exact_gaussian_independence():
    row = o.finite_radial_block_control()
    assert row['native_branch_binding'] == 3960
    assert F(row['normalized_adjacent_block_cross_covariance_trace']) > 0
    assert row['complete_66_state_and_density_oracles_agree_exactly']
    assert not row['adjacent_blocks_independent_at_finite_count']
    assert not row['bound_is_uniform_over_source_refinement']


@pytest.mark.parametrize('energy', (3, 4, 12, 48))
def test_product_spectral_count_bounds_follow_both_extremal_source_spectra(energy):
    row = o.radial_product_bounds(F(energy))
    # These exact spectral oracles saturate the existing source eigenvalue
    # inequalities; the counting proof needs no guessed eigenfunctions.
    high_spectrum_count = sum((2*(energy-k*k))//3
                              for k in range(1, int(np.sqrt(energy))+1))
    low_spectrum_count = sum(49920*(energy-k*k)+1
                             for k in range(1, int(np.sqrt(energy))+1))
    lower = row['certified_positive_angular_mode_count_lower']
    upper = row['certified_including_zero_angular_mode_count_upper']
    assert lower <= high_spectrum_count
    assert low_spectrum_count <= upper
    assert 288*lower**2 >= energy**3
    assert upper**2 <= 49921**2*energy**3
    assert row['spectral_count_power'] == '3/2'
    assert not row['interval_boundary_length_or_metric_selected_by_repair']


def test_radial_gluing_rejects_changed_innovation_and_invalid_contracts(monkeypatch):
    for kwargs in ({'level': True}, {'level': 2}, {'slices': 2},
                   {'length': 2.0}, {'length': F(0)}, {'length': F(5)}):
        with pytest.raises(ValueError):
            o.radial_record_gluing(**kwargs)
    original = o._radial_bridge_step

    def wrong_innovation(*args):
        transfer, variance = original(*args)
        return transfer, variance*0.9

    monkeypatch.setattr(o, '_radial_bridge_step', wrong_innovation)
    with pytest.raises(ArithmeticError, match='radial innovation'):
        o.radial_record_gluing()


def test_radial_readout_rejects_a_forged_native_current_record(monkeypatch):
    original = o.spatial.source.source.binary_step

    def wrong_current(*args):
        occupied, current = original(*args)
        return occupied, current+1

    monkeypatch.setattr(o.spatial.source.source, 'binary_step', wrong_current)
    with pytest.raises(ArithmeticError, match='complete field readback'):
        o.radial_record_gluing()


def test_each_actual_native_seam_separately_preserves_the_full_binary_sector():
    states, targets, marks, score, count = o.marked_source()
    assert count == 3960 and len(states) == 66
    for seam in range(30):
        assert np.array_equal(np.bincount(targets[:, 2*seam:2*seam+2].ravel(), minlength=66),
                              np.full(66, 2))
    assert marks.sum() == 10
    assert np.all(score.sum(axis=0) == 0) and np.all(score.sum(axis=1) == 0)
    assert not targets.flags.writeable and not marks.flags.writeable and not score.flags.writeable


def test_full_native_two_step_law_matches_every_exact_joint_covariance_entry():
    actual = o.direct_two_steps()
    expected = o.finite_covariances(2)
    assert actual['branches'] == 237600
    assert set(actual['terminal_counts']) == {3600}
    for key, matrix in zip(('attempt_second_moment', 'difference_second_moment',
                            'cross_second_moment'), expected, strict=True):
        assert actual[key] == matrix
    assert actual['difference_lag_one'] == sp.zeros(12)
    assert actual['difference_mean'] == sp.zeros(12, 1)
    assert expected[2] != sp.zeros(12)  # The two clocks are not independent.


def test_one_event_and_joint_clt_covariance_have_distinct_quantifiers():
    a, r, cross = o.finite_covariances(1)
    assert a == o.occupation.CHI*o.occupation.P0
    assert r == 5*a and cross == sp.zeros(12)
    cross_limit, delta = o.joint_clt()
    assert cross_limit != cross and delta == r
    a = o.occupation.field('attempt')['covariance']
    joint = np.array(sp.BlockMatrix([[a, cross_limit], [cross_limit.T, delta]]), dtype=float)
    eigenvalues = np.linalg.eigvalsh(joint)
    assert eigenvalues.min() > -2e-13
    assert (eigenvalues > 1e-10).sum() == 22
    assert a+cross_limit+cross_limit.T+delta == 6*o.occupation.field('incident')['covariance']


def test_all_native_joint_innovations_are_centered_and_have_the_predicted_covariance():
    certificate = o.martingale_certificate()
    assert certificate['conditional_states'] == 66
    assert certificate['native_innovations'] == 3960
    assert certificate['all_conditional_means_zero']
    assert certificate['joint_asymptotic_covariance_checked']


def test_stationary_state_marginals_do_not_replace_independent_seam_scheduling():
    sticky = o.direct_two_steps(sticky=True)
    assert sticky['branches'] == 7920
    assert set(sticky['terminal_counts']) == {120}
    assert sticky['difference_mean'] == sp.zeros(12, 1)
    assert sp.trace(sticky['difference_lag_one']) == sp.Rational(250, 33)
    _, delta = o.joint_clt()
    lag = sticky['difference_lag_one']
    assert sticky['difference_second_moment'] == 2*delta+lag+lag.T
    assert sticky['difference_second_moment'] != 2*delta


def test_nonstationary_native_preparation_invalidates_the_reverse_mean_identity():
    biased = o.direct_two_steps(stationary=False)
    assert biased['branches'] == 3600
    assert biased['difference_mean'] != sp.zeros(12, 1)
    assert sum(biased['difference_mean']) == 0  # Total load remains conserved.
    assert biased['difference_second_moment'] != o.finite_covariances(2)[1]


def test_uncentered_clock_sums_do_not_admit_the_centered_random_time_change():
    # Constant readout1 at a two-attempt prefix: N~Binomial(2,1/6).
    probabilities = (F(25, 36), F(10, 36), F(1, 36))
    noise = tuple(6*k-2 for k in range(3))
    assert sum(p*x for p, x in zip(probabilities, noise, strict=True)) == 0
    assert sum(p*x*x for p, x in zip(probabilities, noise, strict=True)) == 10
    # At a prescribed tick count the uncentered sum is deterministic instead.
    assert len({sum([1]*2) for _ in range(3)}) == 1


@pytest.mark.parametrize('level', (0, 1))
def test_path_clock_constants_match_the_actual_source_seams(level):
    raw, _, _, incident, _, _ = o.continuum.uniform.ancestry.problem(level)
    seams = len(raw.seam_a)
    assert seams == 630*4**level and len(incident) == 33
    assert o.clock_difference(level) == (1-F(len(incident), seams))/(4*len(incident)*4**level)


def test_all_level_hilbert_path_bounds_vanish_with_refinement_and_cutoff():
    for level in (0, 3, 10, 100):
        report = o.path_bounds(level, horizon=F(3, 2), count=3, degree=25)
        d = o.clock_difference(level)
        assert report['finite_prefix_clock_sup_squared_upper'] == 30*d  # ceil(9/2)=5
        assert report['joint_Gaussian_clock_sup_squared_upper'] == 12*d
        assert report['joint_Gaussian_clock_sup_squared_upper'] < F(1, 11*4**level)
        assert report['Gaussian_harmonic_sup_tail_squared_upper'] == F(16389, 26**2)
    assert o.path_bounds(100)['joint_Gaussian_clock_sup_squared_upper'] < F(1, 10**60)
    assert o.path_bounds(0, degree=1000)['Gaussian_harmonic_sup_tail_squared_upper'] < F(11, 1000)


@pytest.mark.parametrize('times', ((F(1),), (F(1, 3), F(2), F(7, 3))))
def test_gaussian_cylinder_action_matches_independent_joint_inverse_and_increments(times):
    result = o.cylinder_action(times)
    covariance = sp.kronecker_product(result['time_covariance'], result['spatial_covariance'])
    inverse = covariance.inv(method='DM')
    assert inverse == sp.kronecker_product(result['time_precision'], result['spatial_precision'])
    difference = sp.eye(len(times))
    for i in range(1, len(times)):
        difference[i, i-1] = -1
    assert difference*result['time_covariance']*difference.T == sp.diag(*result['increments'])
    assert result['dimension'] == 11*len(times)


def test_native_finite_count_observations_are_not_independent_brownian_increments():
    native = o.occupation.field('attempt')['kernel']
    adjacent = o.occupation.CHI*(native*o.occupation.P0)[0, 0]
    assert adjacent == sp.Rational(25, 198) > 0
    assert o.finite_covariances(2)[0] != 2*o.finite_covariances(1)[0]


@pytest.mark.parametrize('bad', (True, -1, 1., '2', F(2)))
def test_integer_hypotheses_fail_closed(bad):
    with pytest.raises(ValueError): o.clock_difference(bad)
    with pytest.raises(ValueError): o.path_bounds(0, count=bad)
    with pytest.raises(ValueError): o.finite_covariances(bad)


@pytest.mark.parametrize('bad', (True, 1., F(), F(-1)))
def test_exact_horizon_hypothesis_fails_closed(bad):
    with pytest.raises(ValueError): o.path_bounds(0, horizon=bad)


def test_noncanonical_fractions_and_time_order_are_rejected():
    bad = F(1, 2); bad._numerator = 2; bad._denominator = 4
    with pytest.raises(ValueError): o.path_bounds(0, horizon=bad)
    for times in ((), (F(1), F(1)), (F(2), F(1)), [F(1)]):
        with pytest.raises(ValueError): o.cylinder_action(times)
    with pytest.raises(ValueError): o.path_bounds(0, count=0)
    with pytest.raises(ValueError): o.direct_two_steps(sticky=1)


def test_uniform_bounds_alone_do_not_force_a_unique_whole_sequence_covariance():
    # A control for the inference, not a substitute native source family.
    sequence = [sp.Rational(5+n % 2) for n in range(8)]
    assert min(sequence) > 4 and max(sequence) < o.continuum.uniform.GREEN_BOUND
    assert abs(sequence[-1]-sequence[-2]) == 1


def test_all_runtime_python_sources_are_pinned_in_a_fresh_process():
    code = '''
from pathlib import Path
import sys
from model.physical_space_time import observation_process as o
o.run()
missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if not path: continue
    try: name=Path(path).resolve().relative_to(o.ROOT).as_posix()
    except ValueError: continue
    if not name.startswith('.') and name.endswith('.py') and name not in o.SOURCE_PATHS: missing.append(name)
assert not missing,sorted(set(missing))
'''
    subprocess.run([sys.executable, '-c', code], cwd=o.ROOT, check=True, timeout=120)

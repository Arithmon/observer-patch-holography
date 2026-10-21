from copy import deepcopy
from collections import Counter
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as sp

from ..observer_reconstruction import native_geometry as ng
from . import intrinsic_geometry as g


@pytest.fixture(scope='module')
def radial_records():
    seams = tuple(tuple(edge) for edge in g.recover_structure(ng.execute_calibration())['seams'])
    return seams, g.execute_radial_records(seams)


def test_native_radial_quotient_has_exact_source_spectrum_and_clock_obstruction():
    result = g.native_radial_record_control()
    radial = sp.Matrix(result['radial_generator'])
    mass = sp.diag(1, 5, 5, 1)
    symmetric = mass**sp.Rational(1, 2)*radial*mass**sp.Rational(-1, 2)
    assert symmetric == sp.Matrix([[5, -sp.sqrt(5), 0, 0], [-sp.sqrt(5), 3, -2, 0],
                                  [0, -2, 3, -sp.sqrt(5)], [0, 0, -sp.sqrt(5), 5]])
    assert symmetric.eigenvals() == {0: 1, 6: 1, 5-sp.sqrt(5): 1, 5+sp.sqrt(5): 1}
    stopped = sp.Matrix(result['first_incident_row_kernel'])
    assert stopped*sp.ones(4, 1) == sp.ones(4, 1)
    assert (stopped*radial-radial*stopped)[0, 1] == sp.Rational(13, 37)
    assert result['exchanged_shell_ports_control_rejected']
    assert not result['radial_shells_are_an_independent_product_coordinate']
    assert not result['product_radial_factor_or_physical_length_selected']
    assert not result['source_pulses_for_this_radial_input_reexecuted']
    # Independent spectral cosine evaluates the emitted arc-walk response.
    kernel = np.array(sp.Matrix(result['native_row_kernel']), dtype=float)
    values, vectors = np.linalg.eigh(kernel)
    phases = np.arccos(np.clip(values, -1, 1))
    for row in result['active_wave_complete_Ramsey_outputs']:
        amplitude = np.sum(vectors[0, :]**2*np.cos(row['incident_triggers']*phases))
        assert abs((1+amplitude)/2-float(F(row['Ramsey_plus_probability']))) < 2e-15
        assert F(row['Ramsey_plus_probability'])+F(row['Ramsey_minus_probability']) == 1


def test_native_radial_histories_match_every_two_attempt_mark_and_intermediate_shell(radial_records):
    seams, records = radial_records
    recovered = g.recover_radial_records(records)
    shells = recovered['shells']
    shell_of = {p: r for r, block in enumerate(shells) for p in block}
    kernels = [sp.Matrix(recovered[name]) for name in ('remote_row_kernel', 'incident_row_kernel')]
    for origin in range(12):
        before = tuple(int(p == origin) for p in range(12))
        observed = Counter()
        for edge in seams:
            for coin in (0, 1):
                middle, _ = g.native.native(before, edge, coin)
                first = int(0 in edge), shell_of[middle.index(1)]
                for next_edge in seams:
                    for next_coin in (0, 1):
                        final, _ = g.native.native(middle, next_edge, next_coin)
                        observed[*first, int(0 in next_edge), shell_of[final.index(1)]] += 1
        for m1 in (0, 1):
            for s1 in range(4):
                for m2 in (0, 1):
                    for s2 in range(4):
                        expected = kernels[m1][shell_of[origin], s1]*kernels[m2][s1, s2]
                        assert sp.Rational(observed[m1, s1, m2, s2], 3600) == expected


@pytest.mark.parametrize('observer', range(12))
def test_native_radial_generator_is_observer_relative_with_the_same_law(radial_records, observer):
    _, records = radial_records
    result = g.recover_radial_records(records, observer)
    assert result['shells'][0] == [observer]
    assert result['shell_sizes'] == [1, 5, 5, 1]
    assert sp.Matrix(result['radial_generator']) == sp.Matrix([
        [5, -5, 0, 0], [-1, 3, -2, 0], [0, -2, 3, -1], [0, 0, -5, 5]])


def test_radial_first_incident_resolvent_keeps_remote_history_and_timeout(radial_records):
    _, records = radial_records
    result = g.recover_radial_records(records)
    remote, incident = (sp.Matrix(result[name]) for name in ('remote_row_kernel', 'incident_row_kernel'))
    stopped = sp.Matrix(result['first_incident_row_kernel'])
    completed = sp.zeros(4)
    for n in range(1, 9):
        completed += remote**(n-1)*incident
        assert completed+remote**n*stopped == stopped
        assert remote**n*sp.ones(4, 1) == sp.ones(4, 1)*sp.Rational(5, 6)**n
    assert stopped != 6*incident
    assert remote*incident != incident*remote


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'no_zero_current', 'wrong_coin', 'hidden_endpoint'])
def test_radial_record_reconstruction_rejects_incomplete_or_false_records(radial_records, mutation):
    _, original = radial_records
    records = deepcopy(original)
    if mutation == 'missing': records.pop()
    elif mutation == 'duplicate': records.append(records[0])
    elif mutation == 'no_zero_current': records = [row for row in records if row['current']]
    elif mutation == 'wrong_coin': records[0]['coin'] = 1-records[0]['coin']
    else: records[0]['hidden_endpoint'] = 0
    with pytest.raises(ValueError):
        g.recover_radial_records(records)


def test_averaging_an_inequitable_partition_does_not_define_shell_histories(radial_records):
    _, records = radial_records
    result = g.recover_radial_records(records)
    blocks = result['shells']
    blocks[1][0], blocks[2][0] = blocks[2][0], blocks[1][0]
    with pytest.raises(ValueError, match='not strongly lumpable'):
        g._lumped_kernel(sp.Matrix(result['native_row_kernel']), blocks)


@pytest.fixture(scope='module')
def radial_extensions():
    return g.radial_nonidentifiability_control()


def test_radial_extensions_share_angular_readouts_but_have_distinguishing_local_probes(radial_extensions):
    result = radial_extensions
    assert result['native_ports'] == 240 and result['product_ports'] == 720
    assert [row['radial_spectrum'] for row in result['models']] == [
        ['0', '1', '3'], ['0', '4', '12'], ['0', '3', '3']]
    assert [row['cycle_rank'] for row in result['models']] == [0, 0, 1]
    assert [row['radial_contrast_wave_second_derivative_at_zero']
            for row in result['models']] == ['-5/2', '-10', '-3']
    for row in result['models']:
        radial = sp.Matrix(row['radial_laplacian'])
        incidence = sp.zeros(len(row['radial_edges']), 3)
        for index, (a, b) in enumerate(row['radial_edges']):
            incidence[index, a] = 1; incidence[index, b] = -1
        assert radial == row['radial_edge_weight']*incidence.T*incidence
        assert row['exact_product_intertwiner_residual_nonzeros'] == 0
        # A separately constructed cosine readout detects the radial contrast;
        # its Taylor remainder is bounded by ||R||^2 t^4/24.
        values, vectors = np.linalg.eigh(np.array(radial, float))
        values = np.maximum(values, 0)
        probe = np.array([1., -1., 0.])/np.sqrt(2)
        t = .01
        observed = np.sum((vectors.T@probe)**2*np.cos(t*np.sqrt(values)))
        derivative = float(F(row['radial_contrast_wave_second_derivative_at_zero']))
        assert abs(observed-1-t*t*derivative/2) <= max(values)**2*t**4/24+1e-15
        lo, hi = map(F, row['nonconstant_radial_cosine_interval'])
        assert float(lo) <= observed <= float(hi)
    assert not result['radial_average_is_a_boundary_Dirichlet_to_Neumann_readout']
    assert not result['additional_radial_contrast_is_available_in_current_native_record_menu']
    assert not result['radial_models_are_executed_native_repair_extensions']
    assert not result['existing_Dirichlet_product_or_general_OPH_is_refuted']


def test_radial_extensions_keep_complete_adaptive_clock_reference_histories(radial_extensions):
    result = radial_extensions['complete_adaptive_angular_histories']
    assert result['complete_native_tapes'] == 3600
    assert result['complete_tape_and_field_outcomes'] == 4900
    assert result['reference_dimension'] == 2
    assert result['zero_current_incidents_retained'] > 0
    assert [row['count'] for row in result['clock_pattern_counts']] == [2500, 500, 500, 100]
    rows = result['complete_angular_record_law']
    assert len(rows) == 7 and sum(row['probability'] for row in rows) == pytest.approx(1., abs=1e-12)
    assert rows[0]['field_bits'] == []
    assert rows[0]['probability'] == pytest.approx(25/36, abs=1e-12)
    assert result['second_duration_by_first_bit'] == ['1/4', '1/2']
    assert result['feedback_removal_record_tv'] > 1e-4
    for row in result['radial_model_comparisons']:
        assert row['complete_work_reference_record_distance_upper_numerical'] < 1e-11
        assert row['complete_output_probability'] == pytest.approx(1., abs=1e-12)
    assert not result['source_pulses_for_comparison_operations_executed']
    assert not result['full_product_Gibbs_or_nonconstant_radial_inputs_covered']
    intervals = sorted(tuple(map(F, row['nonconstant_radial_cosine_interval']))
                       for row in radial_extensions['models'])
    assert all(a[1] < b[0] for a, b in zip(intervals, intervals[1:]))


def test_radial_history_equivalence_needs_the_constant_radial_code():
    # The angular superposition of eigenvalues zero and lambda acquires a
    # different reduced coherence under a nonconstant radial preparation.
    # No full-input channel equality follows from the constant-code theorem.
    probe = np.array([1., -1., 0.])/np.sqrt(2)
    radial = np.array([[1., -1., 0.], [-1., 2., -1.], [0., -1., 1.]])
    values, axes = np.linalg.eigh(radial)
    values[0] = 0
    lam, time = 6., .5
    amplitude = np.sum((axes.T@probe)**2*np.exp(-1j*time*(np.sqrt(values+lam)-np.sqrt(values))))
    angular_only = np.exp(-1j*time*np.sqrt(lam))
    assert abs(amplitude-angular_only) > .1


def test_radial_control_rejects_loss_of_the_native_angular_constant(monkeypatch):
    graph, rows, denominator = g._native_density_rows(0)
    altered = [list(row) for row in rows]
    q, count = altered[0][0]
    altered[0][0] = q, count-1
    monkeypatch.setattr(g, '_native_density_rows', lambda level: (graph, altered, denominator))
    with pytest.raises(ArithmeticError, match='angular lifting'):
        g.radial_nonidentifiability_control()


def test_native_optical_front_retains_clock_coordinate_and_detector_scope():
    result = g.native_optical_front_control()
    distance = F(result['constant_coefficient_optical_distance_squared'])
    rate = F(result['wave_coordinate_change_t_equals_rate_times_tau'])
    assert F(result['rescaled_optical_distance_squared']) == distance/rate**2
    assert result['optical_dual_energy_saturated_exactly']
    assert result['metric_null_vector_identity_exact']
    assert result['spatial_chart_covariance_exact']
    dark = result['dark_native_mode']
    graph, _, denominator = g._native_density_rows(0)
    profile = dark['first_carrier_profile']+[0]*(graph.ports-12)
    # Reconstruct the generator directly from the seam incidence, independently
    # of the endpoint-record kernel used by the control.
    laplacian_profile = [0]*graph.ports
    for a, b in graph.seams:
        difference = profile[a]-profile[b]
        laplacian_profile[a] += difference
        laplacian_profile[b] -= difference
    assert laplacian_profile == [6*v for v in profile]
    assert sum(profile) == 0 and sum(v*v for v in profile) == 4
    assert F(dark['uniform_first_carrier_one_attempt_probability_at_detector']) == F(1, 12*denominator)
    assert not dark['finite_matrix_wave_has_a_positive_hard_arrival_delay']
    changed = result['null_set_representative_control']
    assert F(changed['naive_curve_length_squared_on_that_line']) == distance/4
    scope = result['continuum_front_scope']
    assert not scope['finite_coefficient_average_is_the_limiting_metric']
    assert not scope['one_fixed_detector_must_see_every_signal']
    assert not scope['smooth_null_geodesics_or_fastest_native_transport_identified']


@pytest.mark.parametrize('time', [F(1, 4), F(1, 2), F(5, 8), F(3, 4), F(1)])
def test_radial_front_against_independent_dirichlet_spectral_response(time):
    result = g.radial_front_readout(time)
    n = np.arange(1, 4097, dtype=float)
    omega = n*np.pi/4
    source = (np.cos(omega)-np.cos(1.5*omega))/(np.sqrt(2)*omega)
    detector = (np.cos(2*omega)-np.cos(2.5*omega))/(np.sqrt(2)*omega)
    cosine = np.sum(source*detector*np.cos(float(time)*omega))
    sine = np.sum(source*detector*np.sin(float(time)*omega)/omega)
    # Each interval coefficient is at most 4*sqrt(2)/(n*pi). Summing
    # the absolute omitted products gives these independent tail bounds.
    assert abs(cosine-float(F(result['cosine_pairing']))) <= 32/(np.pi**2*4096)
    assert abs(sine-float(F(result['retarded_sine_pairing']))) <= 64/(np.pi**3*4096**2)


def test_optical_duality_for_nonunit_time_weight_and_zero_displacement():
    a = sp.Matrix([[2, 1, 0], [1, 3, 0], [0, 0, 5]])
    v = sp.Matrix([2, -1, 3])
    w = sp.Rational(7, 3)
    result = g.constant_optical_separation(a.tolist(), v, w)
    dual = result['maximizing_covector']
    assert result['optical_tensor'] == w*a.inv()
    assert sp.simplify((dual.T*a*dual)[0]) == w
    assert sp.simplify((dual.T*v)[0]**2) == result['distance_squared']
    for density in (F(1, 10), F(7, 3)):
        metric = g.wave_density_metric(a.tolist(), density, w)['metric']
        assert sp.simplify(metric[1:, 1:]/(-metric[0, 0])) == result['optical_tensor']
    zero = g.constant_optical_separation(a.tolist(), [0, 0, 0], w)
    assert zero['distance_squared'] == 0
    assert zero['maximizing_covector'] == sp.zeros(3, 1)


@pytest.mark.parametrize('vector', [None, [1, 2], [1, 2, 3, 4], [1., 2, 3], [True, 0, 0]])
def test_optical_separation_rejects_inexact_or_malformed_displacements(vector):
    with pytest.raises(ValueError, match='displacement'):
        g.constant_optical_separation(sp.eye(3).tolist(), vector)


@pytest.mark.parametrize('time', [True, .5, F(-1), F(9, 8)])
def test_radial_front_rejects_parameters_outside_its_exact_unreflected_window(time):
    with pytest.raises(ValueError, match='wave parameter'):
        g.radial_front_readout(time)


def test_native_density_kernel_counts_every_endpoint_record_and_hold():
    graph, rows, denominator = g._native_density_rows(1)
    origin = 36
    before = tuple(int(p == origin) for p in range(graph.ports))
    counts = Counter()
    for seam in graph.seams:
        for coin in (0, 1):
            after, (_, current) = g.native.native(before, seam, coin)
            decoded = list(before)
            decoded[seam[0]] -= current
            decoded[seam[1]] += current
            assert tuple(decoded) == after
            counts[after.index(1)] += 1
    assert dict(rows[origin]) == counts
    assert sum(counts.values()) == denominator
    assert counts[origin] == denominator-len(graph.neighbours[origin])


def test_native_increment_readout_retains_chart_and_affine_probe_covariance(monkeypatch):
    before = g.native_chart_increment(1, 2)
    original = g._native_chart_dictionary
    change = sp.Matrix([[2, 1], [-1, 3]])
    def changed(level):
        probes, weights = original(level)
        return tuple((2*x+y+F(5), -x+3*y-F(7)) for x, y in probes), tuple(7*w for w in weights)
    monkeypatch.setattr(g, '_native_chart_dictionary', changed)
    after = g.native_chart_increment(1, 2)
    first = sp.Matrix(before['averaged_screen_tensor'])
    assert sp.Matrix(after['averaged_screen_tensor']) == change*first*change.T


def test_native_increment_control_rejects_cell_reset_and_reversed_limits():
    result = g.native_increment_control()
    assert result['two_attempt_endpoint_record_oracle_exact']
    assert result['zero_current_records_included']
    reset = result['false_cell_reset']
    assert reset['actual_probability'] == '0'
    assert F(reset['reset_probability']) > 0
    assert not result['false_Schur_kernel']['finite_cell_mean_generator_is_a_Markov_generator']
    bare = sp.Matrix(result['bare_chart_tensor'])
    assert (bare-sp.Rational(15, 19)*bare).det() > 0
    values = result['fixed_density_parameter_readouts']
    assert values[0]['density_parameter'] == values[1]['density_parameter']
    assert values[0]['native_attempts'] == 1 and values[1]['native_attempts'] == 16
    assert values[0]['averaged_screen_tensor'] != values[1]['averaged_screen_tensor']
    a = sp.Matrix(result['finite_averaged_product_tensor'])
    metric = sp.Matrix(result['finite_averaged_comparison_metric'])
    assert sp.simplify(sp.sqrt(-metric.det())*metric.inv()) == sp.diag(-1, a)/10
    assert not result['limiting_tensor_pointwise_error_from_finite_run_bounded']
    assert not result['native_time_is_identified_with_wave_or_lab_time']


@pytest.mark.parametrize('level,attempts', [(True, 1), (0, 1), (3, 1),
                                          (1, False), (1, 0), (1, 33), (1, 1.)])
def test_native_increment_rejects_undeclared_execution_parameters(level, attempts):
    with pytest.raises(ValueError):
        g.native_chart_increment(level, attempts)


def test_wave_metric_reproduces_action_and_anisotropic_characteristics():
    a = sp.Matrix([[1, 0, 0], [0, 2, 1], [0, 1, 3]])
    result = g.wave_density_metric(a.tolist(), F(3, 2), F(4, 3))
    metric = result['metric']
    assert metric.det() < 0
    assert metric[0, 0] < 0
    assert all(metric[1:k, 1:k].det() > 0 for k in (2, 3, 4))
    density = sp.sqrt(-metric.det())*metric.inv()
    assert sp.simplify(density-result['action_density']) == sp.zeros(4)
    # Independent dispersion test in a non-eigenvector spatial direction.
    xi = sp.Matrix([1, 1, 1])
    omega = sp.sqrt((xi.T*a*xi)[0]/sp.Rational(4, 3))
    covector = sp.Matrix([-omega, *xi])
    assert sp.simplify((covector.T*metric.inv()*covector)[0]) == 0
    assert (covector.T*sp.diag(-1, 1, 1, 1)*covector)[0] != 0


def test_wave_metric_is_covariant_under_a_spatial_chart_change():
    a = sp.Matrix([[2, 1, 0], [1, 3, 1], [0, 1, 2]])
    change = sp.Matrix([[2, 1, 0], [0, 1, 1], [0, 0, 1]])
    before = g.wave_density_metric(a.tolist(), F(2, 3))
    # x=change*y: gradients transform contragrediently; volume is a density.
    after = g.wave_density_metric(
        (change.inv()*a*change.inv().T).tolist(),
        sp.Rational(2, 3)*abs(change.det()))
    jacobian = sp.diag(1, change)
    assert sp.simplify(after['metric']-jacobian.T*before['metric']*jacobian) == sp.zeros(4)


def test_wave_metric_retains_action_normalization_and_clock_coordinates():
    control = g.product_wave_metric_control()
    assert control['doubling_measure_preserves_cone_but_doubles_metric']
    assert control['clock_reparameterization_is_tensor_covariant']
    assert not control['smooth_coefficients_or_curvature_proved']
    assert not control['observed_spacetime_or_Einstein_dynamics_proved']
    # Both tensors have rank three. Rank alone does not identify their cones.
    first = g.wave_density_metric(sp.eye(3).tolist())['inverse_metric']
    second = g.wave_density_metric(sp.diag(1, 2, 3).tolist())['inverse_metric']
    light_covector = sp.Matrix([-1, 0, 0, 1])
    assert (light_covector.T*first*light_covector)[0] == 0
    assert (light_covector.T*second*light_covector)[0] != 0


@pytest.mark.parametrize('tensor,measure,weight', [
    ([[1, 0], [0, 1]], 1, 1),
    ([[1, 0, 0], [0, 1, 0], [0, 0, 0]], 1, 1),
    ([[1, 0, 0], [0, -1, 0], [0, 0, 1]], 1, 1),
    ([[1, 1, 0], [0, 1, 0], [0, 0, 1]], 1, 1),
    ([[1., 0, 0], [0, 1, 0], [0, 0, 1]], 1, 1),
    (sp.eye(3).tolist(), 0, 1),
    (sp.eye(3).tolist(), 1, -1),
    (sp.eye(3).tolist(), True, 1),
])
def test_wave_metric_rejects_degenerate_or_inadmissible_action(tensor, measure, weight):
    with pytest.raises(ValueError):
        g.wave_density_metric(tensor, measure, weight)


@pytest.fixture(scope='module')
def calibrations():
    records = ng.execute_calibration()
    seams = tuple(sorted({tuple(row['seam']) for row in records}))
    return records, seams, g.execute_biased_calibration(seams)


def test_structural_geometry_survives_preparation_dependent_schedule_and_deterministic_coins(calibrations):
    uniform, _, biased = calibrations
    before = ng.recover(uniform)
    after = g.recover_structure(biased)
    assert len(biased) == 540
    assert {row['coin'] for row in biased if row['preparation_port'] == 0} == {0}
    assert after['seam_exposures_exact'][0] != after['seam_exposures_exact'][1]
    assert after['laplacian_exact'] == before['laplacian_exact']
    assert after['gram_qsqrt5'] == before['gram_qsqrt5']
    assert after['rank_three_structural_gram_certified']


def test_missing_exposure_does_not_inherit_the_icosahedral_certificate(calibrations):
    _, seams, _ = calibrations
    result = g.recover_structure(ng.execute_calibration({seams[0]: F(1)}))
    assert not result['rank_three_structural_gram_certified']


@pytest.mark.parametrize('coin_value', [None, 0, 1])
def test_one_native_shot_per_setting_recovers_gram_without_probabilities(calibrations, coin_value):
    records, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams, coin_value)
    assert len(shots) == 360
    assert all('probability' not in row for row in shots)
    result = g.recover_from_shots(shots)
    assert result['gram_qsqrt5'] == ng.recover(records)['gram_qsqrt5']
    assert not result['probability_law_oracle_used']
    assert result['conditional_current_variance_zero_under_native_model']
    assert not result['laboratory_noise_or_preparation_error_bounded']


def test_unequal_finite_repetition_counts_leave_structural_recovery_unchanged(calibrations):
    _, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams)
    repeated = [row for index, row in enumerate(shots) for _ in range(1+index % 5)]
    assert g.recover_from_shots(shots)['gram_qsqrt5'] == g.recover_from_shots(repeated)['gram_qsqrt5']


def test_finite_shots_reject_a_wrong_native_current_or_probability_oracle(calibrations):
    _, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams)
    changed = deepcopy(shots)
    changed[0]['current'] = 0
    with pytest.raises(ValueError, match='even-load'):
        g.recover_from_shots(changed)
    with pytest.raises(ValueError, match='not probability laws'):
        g.recover_from_shots(calibrations[0])
    with pytest.raises(ValueError, match='inconsistent currents'):
        g.recover_from_shots(shots+changed[:1])


def test_finite_shot_missing_coverage_cannot_recover_the_full_metric(calibrations):
    _, seams, _ = calibrations
    shots = g.execute_single_shot_calibration(seams)
    with pytest.raises(ValueError, match='same positive seam support'):
        g.recover_from_shots(shots[1:])
    incomplete = [row for row in shots if tuple(row['seam']) != seams[0]]
    assert not g.recover_from_shots(incomplete)['rank_three_structural_gram_certified']


def test_relabeling_covariance_without_source_coordinates(calibrations):
    _, _, records = calibrations
    permutation = (7, 2, 10, 0, 5, 11, 8, 4, 1, 9, 3, 6)
    relabelled = deepcopy(records)
    for row in relabelled:
        row['preparation_port'] = permutation[row['preparation_port']]
        a, b = [permutation[x] for x in row['seam']]
        row['seam'] = sorted((a, b))
        if a > b:
            row['current'] *= -1
            row['coin'] = 1-row['coin']
    first, second = g.recover_structure(records), g.recover_structure(relabelled)
    assert second['rank_three_structural_gram_certified']
    for i in range(12):
        for j in range(12):
            assert first['gram_qsqrt5'][i][j] == second['gram_qsqrt5'][permutation[i]][permutation[j]]


@pytest.mark.parametrize('mutate', [
    lambda rows: rows.append(deepcopy(rows[0])),
    lambda rows: rows[0].update(probability='0'),
    lambda rows: rows[0].update(current=0),
    lambda rows: rows[0].update(coin=True),
    lambda rows: rows[0].update(preparation_port=False),
    lambda rows: rows[0].update(seam=[1, 0]),
    lambda rows: rows[0].update(hidden_loads=[]),
    lambda rows: rows.pop(),
])
def test_corrupted_or_incomplete_calibration_rejected(calibrations, mutate):
    rows = deepcopy(calibrations[0])
    mutate(rows)
    with pytest.raises(ValueError):
        g.recover_structure(rows)


def test_preparation_specific_missing_support_rejected(calibrations):
    rows = deepcopy(calibrations[0])
    seam = rows[0]['seam']
    rows = [row for row in rows if row['preparation_port'] != 0 or row['seam'] != seam]
    for row in rows:
        if row['preparation_port'] == 0:
            row['probability'] = '1/58'
    with pytest.raises(ValueError, match='same positive seam support'):
        g.recover_structure(rows)


def test_doubling_any_one_seam_splits_slowest_native_response_from_three_to_two(calibrations):
    records, seams, _ = calibrations
    laplacian = [[F(x) for x in row] for row in g.recover_structure(records)['laplacian_exact']]
    # Check all edges, without relying on an assumed symmetry action.
    for seam in seams:
        result = g.scheduling_control(laplacian, seam)
        assert result['smallest_positive_weighted_laplacian_multiplicity_exact'] == 2
        assert result['normalized_long_attempt_response_gram_rank'] == 2
        assert result['structural_gram_rank_after_exposure_removal'] == 3


def test_full_response_inverse_metric_has_exact_rank_eleven_and_matches_poisson_response(calibrations):
    records, _, _ = calibrations
    laplacian = [[F(x) for x in row] for row in g.recover_structure(records)['laplacian_exact']]
    result = g.inverse_response_metric(laplacian)
    assert result['centered_gram_rank_exact'] == 11
    assert result['squared_distance_shell_counts'] == {'11/30': 30, '7/15': 30, '1/2': 6}
    gram = [[F(x) for x in row] for row in result['gram_exact']]
    resistance = [[F(x) for x in row] for row in result['squared_euclidean_distance_exact']]
    for a in range(12):
        for b in range(a+1, 12):
            potential = [row[a]-row[b] for row in gram]
            forcing = [sum(laplacian[i][j]*potential[j] for j in range(12)) for i in range(12)]
            assert forcing == [F(i == a)-F(i == b) for i in range(12)]
            assert potential[a]-potential[b] == resistance[a][b]


def test_exact_quadratic_field_rank_handles_a_nonrational_pivot():
    z, r = (F(0), F(0)), (F(0), F(1))
    assert g.qrank([[r, z], [z, r]]) == 2
    assert g.qrank([[r, r], [r, r]]) == 1


def test_complete_single_carrier_receipt(calibrations):
    result = g.run()
    assert result['all_pass']
    assert result['same_structural_gram_as_uniform_protocol_exact']
    assert result['finite_shot_recovery']['native_branch_calls'] == 360
    assert not result['scope']['complete_ensemble_probability_law_needed_for_even_preparation_structure']
    assert not result['scope']['federation_needed_for_finite_carrier_result']
    assert not result['scope']['physical_length_time_or_lab_calibration_established']

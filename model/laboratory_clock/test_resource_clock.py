from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.optimize import linprog

from . import resource_clock as r


def point(name, n, p, seconds, error=F(0)):
    return r.ResourceEndpoint(name, n, p, seconds-error, seconds+error)


@pytest.fixture(scope='module')
def calibration():
    return (point('zero', 0, 0, F(0)), point('native', 100, 0, F(1, 10), F(1, 10**6)),
            point('mixed', 100, 100, F(51, 500), F(1, 10**6)))


@pytest.fixture(scope='module')
def receipt():
    return r.run()


def test_two_operation_durations_are_identified_without_equating_them(calibration):
    fit = r.calibration(calibration)
    assert fit['native_seconds_bounds'] == [F(99999, 10**8), F(100001, 10**8)]
    assert fit['pulse_seconds_bounds'] == [F(999, 50000000), F(1001, 50000000)]
    assert fit['both_durations_bounded_and_bounded_away_from_zero']
    assert not fit['physical_timing_hypothesis_authenticated']


def test_fixed_resource_ratio_leaves_a_continuum_of_calibrations():
    points = (point('zero', 0, 0, F(0)), point('mixed', 100, 200, F(13, 125)))
    fit = r.calibration(points)
    assert fit['strict_positive_durations_feasible']
    assert not fit['both_durations_bounded_and_bounded_away_from_zero']
    assert r.duration_range(points, F(100), F(0)) == (F(0), F(13, 125))
    # Both positive duration pairs reproduce all calibration observations.
    for tau, kappa in ((F(1, 1000), F(1, 50000)), (F(1, 1250), F(3, 25000))):
        assert 100*tau+200*kappa == F(13, 125)


def test_missing_pulse_exposure_is_unbounded_not_a_zero_rate_estimate():
    points = (point('zero', 0, 0, F(0)), point('native', 100, 0, F(1, 10)))
    fit = r.calibration(points)
    assert fit['pulse_seconds_bounds'] == [F(0), None]
    assert fit['strict_positive_durations_feasible']
    assert r.duration_range(points, F(0), F(51)) == (F(0), None)


def test_zero_duration_boundary_is_not_positive_feasibility():
    points = (point('zero', 0, 0, F(0)), point('native', 1, 0, F(0)),
              point('pulse', 0, 1, F(1)))
    fit = r.calibration(points)
    assert fit['closed_outer_feasible']
    assert not fit['strict_positive_durations_feasible']


def test_joint_validation_can_fail_despite_every_individual_overlap():
    calibration = (point('zero', 0, 0, F(0)), point('mixed', 1, 1, F(2)))
    held = (point('v_native', 1, 0, F(1, 2)), point('v_pulse', 0, 1, F(1, 2)))
    result = r.validate(calibration, held)
    assert result['individually_compatible'] == [True, True]
    assert not result['jointly_compatible_with_some_original_calibration']
    assert not result['held_out_used_in_calibration']


def test_held_out_shift_rejected_without_refitting(receipt):
    assert receipt['held_out']['jointly_compatible_with_some_original_calibration']
    assert not receipt['held_out_changed']['jointly_compatible_with_some_original_calibration']
    assert receipt['held_out']['calibration'] == receipt['held_out_changed']['calibration']


def test_resource_correlation_is_retained_in_expected_duration(receipt):
    joint = list(map(F, receipt['expected_mean_duration_seconds_bounds']))
    box = list(map(F, receipt['independent_marginal_box_mean_duration_seconds_bounds']))
    assert box[0] < joint[0] < joint[1] < box[1]
    assert joint == [F(25689909, 6250000000), F(25700091, 6250000000)]


@pytest.mark.parametrize('gamma,cap', [(0, 1), (1, 1), (0, 2), (1, 2), (0, 3), (1, 3)])
def test_counted_history_forgets_to_every_native_output_matrix(gamma, cap):
    law = r.resource_law(gamma, cap)  # Checks the full matrix identity internally.
    q = 1-F(3, 5)**cap
    assert sum(law.values(), F(0)) == 1
    assert sum(v for (n, p, s, a, b, bit), v in law.items() if s.startswith('timeout')) == 1-q*q
    assert sum(n*v for (n, p, s, a, b, bit), v in law.items()) == q*(1+q)/F(2, 5)+q
    assert sum(p*v for (n, p, s, a, b, bit), v in law.items()) == 27+24*q


def test_resource_marginal_and_joint_setting_response(receipt):
    assert receipt['resource_marginal'] == [
        {'attempts': 2, 'pulses': 27, 'probability': '9/25'},
        {'attempts': 3, 'pulses': 51, 'probability': '4/25'},
        {'attempts': 4, 'pulses': 51, 'probability': '42/125'},
        {'attempts': 5, 'pulses': 51, 'probability': '18/125'}]
    assert receipt['joint_resource_clock_detector_setting_response_tv'] == '1848/15625'


def test_deadline_retains_timeouts_and_coherent_response(calibration):
    outputs = [r.deadline(r.resource_law(g), calibration, F(9, 2000),
                         timestamp_error=F(1, 10**6), detector=1) for g in (0, 1)]
    assert outputs[0]['event_probability_bounds'] != outputs[1]['event_probability_bounds']
    assert outputs[0]['event_probability_bounds'] == [F(2451, 10000)]*2
    assert outputs[1]['event_probability_bounds'] == [F(543, 2000)]*2
    assert all(x['bounds_sharp_for_ordered_cost_support'] and x['timeouts_retained'] for x in outputs)
    total = r.deadline(r.resource_law(), calibration, F(9, 2000), timestamp_error=F(1, 10**6))
    assert total['event_probability_bounds'] == [F(13, 25), F(13, 25)]
    assert r.deadline(r.resource_law(), iter(calibration), F(9, 2000),
                      timestamp_error=F(1, 10**6)) == total
    edge = r.deadline(r.resource_law(), calibration, F(201, 50000))
    assert edge['event_probability_bounds'] == [F(9, 25), F(13, 25)]


def test_mean_duration_cannot_replace_the_stopped_history(receipt):
    control = receipt['mean_replacement_negative_control']
    assert control['exact_cadence_deadline']['event_probability_bounds'] == ['9/25', '9/25']
    assert control['probability_from_deterministic_mean_duration'] == '0'


def test_raw_mass_defects_are_rejected_without_probability_normalization(calibration):
    law = r.resource_law()
    changed = {key: value*(1+F(1, 10**12)) for key, value in law.items()}
    with pytest.raises(ValueError, match='normalized'):
        r.deadline(changed, calibration, F(1))


def test_exact_projection_matches_independent_three_variable_linear_program():
    rng = np.random.default_rng(408125)
    for case in range(60):
        tau, kappa, offset = [F(int(x), 100) for x in rng.integers(1, 100, 3)]
        exposures = [(0, 0), (10, 0), (0, 10)]
        exposures += [tuple(map(int, pair)) for pair in rng.integers(0, 20, (3, 2))]
        points = tuple(point(f'p{i}', n, p, offset+n*tau+p*kappa,
                             F(int(rng.integers(1, 20)), 100))
                       for i, (n, p) in enumerate(exposures))
        n, p = map(int, rng.integers(0, 30, 2))
        rows, bounds = [], []
        for endpoint in points:
            row = [1, endpoint.attempts, endpoint.pulses]
            rows.extend([row, [-v for v in row]])
            bounds.extend([float(endpoint.seconds_hi), -float(endpoint.seconds_lo)])
        actual = r.duration_range(points, F(n), F(p))
        for sign, expected in ((1, actual[0]), (-1, actual[1])):
            lp = linprog([0, sign*n, sign*p], A_ub=rows, b_ub=bounds,
                         bounds=[(None, None), (0, None), (0, None)], method='highs')
            assert lp.success
            assert sign*lp.fun == pytest.approx(float(expected), abs=1e-10)


def test_counter_cache_is_copy_safe_and_does_not_accept_boolean_aliases():
    history = r.counted_history(1, 2)
    key = next(key for key, state in history.items() if state)
    entry = next(iter(history[key]))
    original = F(history[key][entry].numerator, history[key][entry].denominator)
    history[key][entry]._numerator = 999999
    assert r.counted_history(1, 2)[key][entry] == original
    history[key].clear()
    assert r.counted_history(1, 2)[key]
    with pytest.raises(ValueError):
        r.counted_history(True, 2)
    with pytest.raises(ValueError):
        r.counted_history(1, True)


def test_bad_bounds_and_duplicate_calibration_validation_ids_rejected(calibration):
    with pytest.raises(ValueError):
        replace(calibration[0], pulses=True)
    with pytest.raises(ValueError):
        replace(calibration[0], seconds_lo=F(1))
    with pytest.raises(ValueError):
        r.validate(calibration, (calibration[0],))


@pytest.mark.parametrize('numerator,denominator', [(True, 2), (2, 4), (1, 0), (1, -2)])
def test_noncanonical_mutated_rationals_are_rejected(numerator, denominator, calibration):
    bad = F(1, 2)
    bad._numerator, bad._denominator = numerator, denominator
    with pytest.raises(ValueError):
        r.ResourceEndpoint('bad', 1, 1, bad, F(2))
    with pytest.raises(ValueError):
        r.duration_range(calibration, bad, F(1))
    with pytest.raises(ValueError):
        r.deadline(r.resource_law(), calibration, bad)
    corrupt = point('mutated', 1, 1, F(1))
    corrupt.seconds_lo._numerator, corrupt.seconds_lo._denominator = numerator, denominator
    with pytest.raises(ValueError):
        r.calibration((corrupt,))


def test_endpoint_detaches_caller_owned_fraction_slots():
    value = F(1, 2)
    endpoint = r.ResourceEndpoint('detached', 1, 1, value, value)
    value._numerator = 77
    assert endpoint.seconds_lo == endpoint.seconds_hi == F(1, 2)


@pytest.mark.parametrize('cap', [1, 2, 3])
def test_first_hidden_wait_factorizes_from_every_complete_incident_output(cap):
    control = r.first_wait_obstruction(cap)
    assert control['supported_matrix_units_checked'] == 9
    for reference in (False, True):
        initial = r.oracle.preparation(reference=reference)
        first, timeout = r._interval(initial, cap)
        by_mark = {}
        for count, mark, state in first:
            normalized = {key: value/(1-r.target.P)**(count-1) for key, value in state.items() if value}
            if mark in by_mark:
                assert normalized == by_mark[mark]
            by_mark[mark] = normalized
        assert timeout == {key: value*(1-r.target.P)**cap for key, value in initial.items()}
    if cap == 2:
        assert control['first_wait_given_success'] == [F(5, 8), F(3, 8)]
        assert control['total_attempt_classification_error_lower'] == F(6, 25)
        assert control['total_attempt_squared_error_lower'] == F(3, 20)


@pytest.mark.parametrize('gamma', [0, 1])
def test_visible_exposure_bayes_risks_and_no_hidden_count_service(gamma):
    result = r.visible_exposure(gamma)
    assert not result['hidden_counter_readout_added']
    assert result['ambiguous_record_mass'] == F(16, 25)
    assert result['minimum_record_only_classification_error'] == F(38, 125)
    assert result['minimum_record_only_squared_error'] > F(3, 20)
    # Recover both the unconditional attempt distribution and its Bayes risk
    # directly from posterior records, including the known first-timeout cost.
    recovered = {}
    direct_risk = F(0)
    for row in result['posterior_rows']:
        posterior = dict(row['attempt_posterior'])
        assert sum(posterior.values(), F(0)) == 1
        for count, mass in posterior.items():
            recovered[count] = recovered.get(count, F(0))+row['mass']*mass
        losses = {guess: 1-posterior.get(guess, F(0)) for guess in range(2, 6)}
        direct_risk += row['mass']*min(losses.values())
    assert recovered == {2: F(9, 25), 3: F(4, 25), 4: F(42, 125), 5: F(18, 125)}
    assert direct_risk == result['minimum_record_only_classification_error']


def test_equal_quantum_output_blocks_exact_attempt_inference():
    witness = r.indistinguishable_exposure()
    assert witness['native_attempt_masses'] == [(4, F(9, 1000)), (5, F(27, 5000))]
    assert witness['same_terminal_occupied_ports'] == [1, 5]
    assert not witness['terminal_quantum_measurement_can_identify_the_count']
    # Arbitrary effects have equal response in the two normalized output
    # states, so changing the terminal measurement cannot reveal the count.
    record = witness['visible_record']
    states = [state for (n, word, status, bit), state in r.counted_history().items()
              if (word, status, bit) == record and r.oracle.trace(state)]
    normalized = [{key: value/r.oracle.trace(state) for key, value in state.items()} for state in states]
    assert normalized[0] == normalized[1]


def test_deadline_forecast_marginalizes_latent_cost_without_granting_read_access(calibration):
    record = r.indistinguishable_exposure()['visible_record']
    prediction = r.visible_deadline(record, calibration, F(11, 2000), timestamp_error=F(1, 10**6))
    assert prediction['record_probability'] == F(9, 625)
    assert prediction['joint_record_and_deadline_probability_bounds'] == [F(9, 1000)]*2
    assert prediction['deadline_given_record_probability_bounds'] == [F(5, 8)]*2
    assert prediction['latent_costs_summed_without_counter_readout']
    assert prediction['elapsed_clock_observation_is_additional_hypothesis']
    assert prediction['timeouts_retained']
    with pytest.raises(ValueError, match='positive-probability'):
        r.visible_deadline(('not-a-record',), calibration, F(1))

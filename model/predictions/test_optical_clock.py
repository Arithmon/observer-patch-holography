"""Exact controls for extraction, concentration bounds and conditional scope."""
from datetime import datetime, timedelta
from fractions import Fraction as F
from math import comb

import pytest

from . import optical_clock as o


def fixture(seconds=900, omit=()):
    start = datetime(2021, 1, 1)
    rows = ['Date Time(UTC) Ratio_36ar_yb']
    for index in range(seconds):
        if index not in omit:
            stamp = start+timedelta(seconds=index)
            rows.append(stamp.strftime('%Y-%m-%d\t%H:%M:%S')+'\t1.050000000000000000000000001')
    return ('\n'.join(rows)+'\n').encode()


def test_complete_calendar_blocks_do_not_bridge_missing_rows():
    full = o.summarize(fixture(), '36ar_yb_ratio.txt')
    missing = o.summarize(fixture(omit=(450,)), '36ar_yb_ratio.txt')
    assert full['complete_blocks'] == 3 and full['adjacent_complete_block_pairs'] == 2
    assert missing['complete_blocks'] == 2 and missing['possible_calendar_blocks'] == 3
    assert missing['adjacent_complete_block_pairs'] == 0
    assert missing['nonoverlapping_Allan_variance_diagnostic'] is None
    assert full['minimum_block_ratio'] == str(F('1.050000000000000000000000001'))
    assert full['nonoverlapping_Allan_variance_diagnostic'] == '0'


@pytest.mark.parametrize('change', ['duplicate', 'nan', 'reorder', 'fractional_second', 'wrong_header'])
def test_reject_malformed_original_observations(change):
    lines = fixture().decode().splitlines()
    if change == 'duplicate':
        lines.insert(3, lines[2])
    elif change == 'nan':
        lines[2] = lines[2].rsplit('\t', 1)[0]+'\tNaN'
    elif change == 'reorder':
        lines[2], lines[3] = lines[3], lines[2]
    elif change == 'fractional_second':
        lines[2] = lines[2].replace('00:00:01', '00:00:01.5')
    else:
        lines[0] = 'other apparatus'
    with pytest.raises(ValueError):
        o.summarize(('\n'.join(lines)+'\n').encode(), '36ar_yb_ratio.txt')


@pytest.mark.parametrize('delta', [0, 1])
def test_geometric_parity_and_native_renewal_factorization(delta):
    parameters = o.event_geometry.renewal_parameters(delta)
    p = parameters['geometric_success_probability']; q = 1-p
    parity = q/(1+q)
    for value in range(12):
        # Independent parity Bernoulli and a doubled Geom0(1-q²).
        split = (parity if value % 2 else 1-parity)*(1-q*q)*q**(2*(value//2))
        assert split == p*q**value
    terminal = parameters['terminal_B_probability']
    for value in range(1, 12):
        probability = (1-terminal)*p*q**value+terminal*p*q**(value-1)
        assert probability == F(25, (10-delta)**2)*q**(value-1)
    assert (1-terminal)*p == F(5-delta, 10-delta)
    for count in (16, 32, 64):
        maximum = max(F(comb(count, k))*parity**k*(1-parity)**(count-k)
                      for k in range(count+1))
        assert maximum <= o.max_point_mass_bound(count)


def test_selection_certificate_is_an_upper_bound_on_anchored_union():
    total, retained, mass = 100, 80, F(1, 10000)
    result = o.selection_bound(total, retained, mass)
    actual_union_bound = total*comb(total-1, retained-1)*mass**(retained-1)
    binary = total*F(1, 2)**(result['binary_power_per_retained_draw']*(retained-1))
    decimal = F(10)**result['conservative_probability_upper_power_of_ten']
    assert actual_union_bound <= binary <= decimal < F(1, 100)
    assert o.selection_bound(100, 2, F(1))['conservative_probability_upper_power_of_ten'] == 0
    selected = o.selection_bound(total, retained, mass, choice_budget=10**6)
    assert selected['conservative_probability_upper_power_of_ten'] == (
        result['conservative_probability_upper_power_of_ten']+6)


def test_calendar_search_is_charged_and_does_not_include_unlimited_stopping():
    row = o.strict_json.load(o.SUMMARY)['datasets'][0]
    result = o.evaluate(row)
    assert result['calendar_and_count_selection_union_budget'] == o.CALENDAR_BLOCKS**3
    with pytest.raises(ValueError, match='calendar'):
        o.evaluate({**row, 'first_UTC_second': o.CALENDAR_START-1})
    with pytest.raises(ValueError, match='calendar'):
        o.summarize(fixture().replace(b'2021-', b'2020-'), '36ar_yb_ratio.txt')


def test_public_comparison_charges_readout_error_and_retains_its_scope():
    row = o.strict_json.load(o.SUMMARY)['datasets'][0]
    result = o.evaluate(row)
    assert result['attachment_rejected_at_one_percent_under_stated_hypotheses']
    enlarged = o.evaluate(row, error=F(1, 1000))
    assert not enlarged['attachment_rejected_at_one_percent_under_stated_hypotheses']
    assert F(o.evaluate(row, gain_lower=F(1, 2))['single_window_interval_mass_upper']) > F(
        result['single_window_interval_mass_upper'])


@pytest.mark.parametrize('count', [True, 1, -1, 2.0])
def test_exact_count_contract(count):
    with pytest.raises(ValueError):
        o.max_point_mass_bound(count)


def test_selection_bound_multi_draw_binary_exponent_when_per_draw_base_exceeds_half():
    total, retained, mass = 6978, 210, F(885342, 168243169)
    result = o.selection_bound(total, retained, mass, choice_budget=o.CALENDAR_BLOCKS**3)
    assert result['binary_power_per_retained_draw'] == 0
    assert F(1, 2) < F(result['interval_occupancy_union_base_upper']) < F(1)
    actual_union_bound = total*o.CALENDAR_BLOCKS**3*comb(total-1, retained-1)*mass**(retained-1)
    decimal = F(10)**result['conservative_probability_upper_power_of_ten']
    assert actual_union_bound <= decimal <= F(1, 10**38)


def test_servo_and_refinement_analysis_oracles_and_decision_boundary():
    receipt = o.run()
    high = (o.REFERENCE_HZ+1)*o.BLOCK_SECONDS
    expected_floors = {
        '36ar_yb_ratio.txt': (2544097785170, 3180122231462, -123),
        '40ar_yb_ratio.txt': (2407037601662, 3008797002078, -38),
    }
    for row, servo in zip(receipt['data']['datasets'], receipt['servo_and_refinement_tests']):
        obs = F(row['nonoverlapping_Allan_variance_diagnostic'])
        adj_floor, nonadj_floor, dec_exp = expected_floors[row['file']]
        assert servo['servo_filtered_multi_replica_attachment_rejected_at_one_percent']
        assert servo['conservative_probability_upper_power_of_ten'] == dec_exp
        assert servo['adjacent_to_observed_Allan_variance_ratio_floor'] == adj_floor
        assert servo['nonadjacent_to_observed_Allan_variance_ratio_floor'] == nonadj_floor
        r_min = servo['moment_matching_independent_native_replica_count_lower']
        assert r_min == adj_floor+1
        assert F(8, 5*high*r_min) <= obs < F(8, 5*high*(r_min-1))
        nu_min = servo['moment_matching_unfiltered_A_tick_frequency_Hz_lower']
        assert F(8, 5*o.BLOCK_SECONDS*nu_min) <= obs < F(8, 5*o.BLOCK_SECONDS*(nu_min-1))
        assert servo['moment_target_requires_unfiltered_zero_error_population_variance_matching']
        assert not servo['sample_Allan_value_is_a_population_upper_confidence_bound']
        assert not servo['continuum_refinement_level_inferred_from_public_data']
        enlarged = o.servo_and_refinement_analysis(row, error=F(1, 1000))
        assert not enlarged['servo_filtered_multi_replica_attachment_rejected_at_one_percent']


def test_servo_and_refinement_rejects_invalid_inputs():
    row = o.strict_json.load(o.SUMMARY)['datasets'][0]
    malformed = F(1, 2); malformed._denominator = -2
    for bad_replicas in (0, -1, True, 8.0):
        with pytest.raises(ValueError):
            o.servo_and_refinement_analysis(row, replicas=bad_replicas)
    for bad_servo in (-1, o.BLOCK_SECONDS, True, 30.0):
        with pytest.raises(ValueError):
            o.servo_and_refinement_analysis(row, servo_seconds=bad_servo)
    for bad_frac in (malformed, True, -F(1), 0.5):
        with pytest.raises(ValueError):
            o.servo_and_refinement_analysis(row, error=bad_frac)
        with pytest.raises(ValueError):
            o.servo_and_refinement_analysis(row, gain_lower=bad_frac)
    for bad_row in ({**row, 'adjacent_complete_block_pairs': 0},
                    {**row, 'nonoverlapping_Allan_variance_diagnostic': None},
                    {**row, 'nonoverlapping_Allan_variance_diagnostic': '0'}):
        with pytest.raises(ValueError):
            o.servo_and_refinement_analysis(bad_row)


def test_finite_filter_conditioning_is_block_diagonal_only_for_interior_ticks():
    # Three six-tick blocks and a causal two-tick-memory unit-mass filter.
    kernel = (F(1, 4), F(1, 2), F(1, 4))
    def block_weights(tick):
        return tuple(sum((weight for lag, weight in enumerate(kernel)
                          if 6*block < tick+lag <= 6*(block+1)), F(0))
                     for block in range(3))
    for block in range(3):
        for tick in range(6*block+1, 6*block+5):
            assert block_weights(tick) == tuple(F(j == block) for j in range(3))
    assert block_weights(6) == (F(1, 4), F(3, 4), F(0))
    # A decaying response with any finite time constant still has an infinite tail.
    assert sum((F(1, 2)**(lag+1) for lag in range(3)), F(0)) == F(7, 8)
    hypothesis = o.servo_hypothesis()
    assert not hypothesis['published_servo_time_constant_establishes_compact_support']
    assert not hypothesis['eight_published_servo_outputs_establish_independent_native_replicas']
    assert hypothesis['A_ticks_equally_spaced_in_declared_calendar_coordinate']


def test_refinement_clock_normalization_counts_actual_source_seams():
    from ..relative_observer_geometry import temporal_records as temporal
    for level in (0, 1):
        graph = temporal.clocks.screen.source.federation(level)
        for degree in (5, 6):
            port = next(i for i, row in enumerate(graph.neighbours) if len(row) == degree)
            incident = sum(port in seam for seam in graph.seams)
            scale = o.refined_incident_clock_scale(level, degree)
            assert scale['source_seams'] == len(graph.seams)
            assert F(scale['incident_probability']) == F(incident, len(graph.seams))
            expected = F(1260*16**level*incident, len(graph.seams))
            assert F(scale['mean_incident_ticks_per_density_continuum_unit']) == expected == 2*degree*4**level
            assert expected != 210*16**level
    for level, degree in ((True, 5), (-1, 5), (0, True), (0, 30), (1.0, 5)):
        with pytest.raises(ValueError):
            o.refined_incident_clock_scale(level, degree)


def test_complete_incident_clock_statistics_fix_the_native_field_transfer():
    certificate = o.clock_field_consistency()
    assert certificate['native_binary_histories_checked'] == 3960
    assert certificate['port_mean_tick_probability'] == ['1/6']*12
    assert not certificate['second_field_rate_fitted']
    # Dropping the shared-event covariance predicts a false direct coupling
    # between nonadjacent ports and loses the native seam coupling as well.
    p = F(1, 6)
    assert p*p/2 != F(0)
    assert (F(1, 180)+p*p)/2 == F(1, 60)
    assert (-F(1, 36)+p*p)/2 == F(0)


def test_clock_field_sampling_bound_preserves_shared_diagonal_and_covariance():
    certificate = o.clock_field_sampling_bound(3900000, F(1, 100))
    assert certificate['failure_probability_upper'] == '1/20'
    assert certificate['density_transfer_operator_error_upper'] == '3/50'
    # An error delta in every second-moment entry changes every matrix entry
    # by delta/2 in magnitude, because its diagonal is reused for p.
    delta = F(1, 100)
    for i in range(12):
        for j in range(12):
            error = -delta*int(i == j)+delta/2
            assert abs(error) == delta/2
    for events, delta in ((True, F(1, 10)), (0, F(1, 10)), (1, F(0)), (1, 0.1)):
        with pytest.raises(ValueError):
            o.clock_field_sampling_bound(events, delta)

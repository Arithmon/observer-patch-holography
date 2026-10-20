"""Exact concentration, finite selection cost, and attachment boundary controls."""
from fractions import Fraction as F
from math import comb

import pytest

from . import regional_optical as r


@pytest.fixture(scope='module')
def receipt():
    return r.run()


@pytest.mark.parametrize('counts', [(5, 5, 1), (5, 30, 5), (30, 5, 5), (29, 30, 29), (29, 29, 28)])
def test_independent_geometric_parity_factor_and_point_mass_bound(counts):
    law = r.regional.parameters(*counts)
    p, terminal = law['geometric_success_probability'], law['terminal_B_probability']
    q = 1-p; parity = q/(1+q)
    for value in range(20):
        assert p*q**value == (parity if value % 2 else 1-parity)*(1-q*q)*q**(2*(value//2))
    choice = max((terminal, parity), key=lambda x: x*(1-x))
    assert r.bernoulli_variance(*counts) == choice*(1-choice)
    for count in (16, 32, 64):
        maximum = max(F(comb(count, k))*choice**k*(1-choice)**(count-k) for k in range(count+1))
        assert maximum <= r.point_mass_bound(*counts, count)


def test_all_distinct_incidence_classes_have_positive_independent_noise_factor():
    for a, b, common, _ in r.regional.pair_classes():
        assert r.bernoulli_variance(a, b, common) > 0
        bound = r.point_mass_bound(a, b, common, 10**12)
        assert 0 < bound < 1
        assert bound**2 * 10**12 >= F(22, 56)/r.bernoulli_variance(a, b, common)
    with pytest.raises(ValueError, match='identical'):
        r.point_mass_bound(30, 30, 30, 100)


def test_exact_selection_bound_keeps_probability_before_decimal_rounding():
    total, retained, mass, choices = 100, 80, F(1, 10000), 10**9
    result = r.selection_bound(total, retained, mass, choice_budget=choices)
    exact = total*choices*comb(total-1, retained-1)*mass**(retained-1)
    exponent = result['conservative_probability_upper_power_of_ten']
    assert F(10)**(exponent-1) < exact <= F(10)**exponent
    assert r.selection_bound(10, 2, F(1))['conservative_probability_upper_power_of_ten'] == 0
    assert r.decimal_ceiling(F(1, 10**70)) == -70
    assert r.decimal_ceiling(F(2, 10**70)) == -69


def test_public_bound_tests_every_distinct_regional_pair_but_not_aliases(receipt):
    assert receipt['family']['distinct_event_sets'] == 3433
    assert receipt['family']['clock_aliases_excluded_from_distinct_pair_test']
    for row in receipt['conditional_tests']:
        assert row['all_distinct_regional_direct_attachments_rejected_at_one_percent']
        assert row['ordered_distinct_event_set_pair_budget'] == 11782056
        assert row['finite_selection_union_budget'] == r.optical.CALENDAR_BLOCKS**3*11782056
        assert row['worst_incidence_counts'] == [5, 30, 5]
    assert receipt['hypothesis']['extra_temporal_filtering_or_independent_replica_averaging_included'] is False
    assert receipt['assessment']['atomic_clocks_or_OPH_excluded'] is False


def test_larger_error_or_smaller_gain_weakens_the_actual_count_window():
    row = r.optical.strict_json.load(r.optical.SUMMARY)['datasets'][1]
    baseline = r.evaluate(row)
    wider = r.evaluate(row, error=F(1, 10**6))
    assert not wider['all_distinct_regional_direct_attachments_rejected_at_one_percent']
    assert F(r.evaluate(row, gain_lower=F(1, 2))['single_window_interval_mass_upper_for_every_distinct_pair']) > F(
        baseline['single_window_interval_mass_upper_for_every_distinct_pair'])


def test_strict_contract_rejects_mutated_fraction_bool_and_undeclared_calendar():
    malformed = F(1, 2); malformed._denominator = -2
    for value in (malformed, True, 0.5, F(-1)):
        with pytest.raises(ValueError):
            r.selection_bound(10, 5, value)
    for count in (True, 1, 2.0):
        with pytest.raises(ValueError):
            r.point_mass_bound(5, 5, 1, count)
    row = r.optical.strict_json.load(r.optical.SUMMARY)['datasets'][0]
    with pytest.raises(ValueError):
        r.evaluate({**row, 'first_UTC_second': r.optical.CALENDAR_START-1})
    with pytest.raises(ValueError):
        r.evaluate({**row, 'complete_blocks': True})


def test_receipt_rejects_undeclared_summary_fields(monkeypatch):
    summary = r.optical.strict_json.load(r.optical.SUMMARY)
    monkeypatch.setattr(r.optical.strict_json, 'load', lambda _: {**summary, 'extra': True})
    with pytest.raises(ValueError, match='retained summaries'):
        r.run()


@pytest.mark.parametrize('regions', [
    ((0,), (1,)),
    ((0,), tuple(range(12))),
    (tuple(p for p in range(12) if p not in (0, 1)), tuple(range(12))),
    (tuple(p for p in range(12) if p not in (0, 1)), tuple(p for p in range(12) if p not in (0, 5))),
])
def test_binary_geometric_digits_and_exact_discrete_convolution_oracle(regions):
    first, second = regions
    a, b, common = r.regional.incidence_counts(first, second)
    law = r.regional.parameters(a, b, common)
    p = law['geometric_success_probability']; q = 1-p
    for z in (F(0), F(1, 3), F(2, 3)):
        # Euler binary product identity for Geom0(p) up to 2^3.
        binary_prefix = F(1)
        for j in (0, 1, 2):
            qj = q**(1 << j)
            rj = qj/(1+qj)
            binary_prefix *= (1-rj+rj*z**(1 << j))
        remainder = (1-q**8)/(1-(q*z)**8)
        assert binary_prefix*remainder == p/(1-q*z)
    for count in (8, 16):
        bins, tail = r.regional.count_law(first, second, count=count, maximum=64)
        for width in (0, 1, 4, 9):
            # Any window of width w in [0, 64] has mass at most max window + tail.
            max_window = max(sum(bins[k:k+width+1], F(0)) for k in range(len(bins)-width))
            bound = r.servo_interval_mass_bound(a, b, common, count, width)
            assert max_window <= bound


def test_regional_servo_and_allan_analysis_oracles_and_controls(receipt):
    summary = r.optical.strict_json.load(r.optical.SUMMARY)['datasets']
    high = (r.optical.REFERENCE_HZ+1)*r.optical.BLOCK_SECONDS
    expected = {
        '36ar_yb_ratio.txt': (53002037192, 109659387292, -123),
        '40ar_yb_ratio.txt': (50146616702, 103751620762, -25),
    }
    for row, servo in zip(summary, receipt['servo_and_allan_tests']):
        obs = F(row['nonoverlapping_Allan_variance_diagnostic'])
        r_min_exp, r_eq_exp, dec_exp = expected[row['file']]
        assert servo['all_distinct_regional_servo_attachments_rejected_at_one_percent']
        assert servo['conservative_probability_upper_power_of_ten'] == dec_exp
        assert servo['worst_incidence_counts'] == [29, 30, 29]
        r_min = servo['moment_matching_independent_regional_replica_count_lower']
        r_eq = servo['moment_matching_independent_equal_rate_regional_replica_count_lower']
        assert r_min == r_min_exp and r_eq == r_eq_exp
        assert F(1, 30*high*r_min) <= obs < F(1, 30*high*(r_min-1))
        assert F(2, 29*high*r_eq) <= obs < F(2, 29*high*(r_eq-1))
        assert servo['moment_target_requires_unfiltered_zero_error_population_variance_matching']
        assert not servo['sample_Allan_value_is_a_population_upper_confidence_bound']
        assert not servo['thirty_seam_regional_floor_transferred_to_refining_graph']
        wider = r.servo_and_allan_analysis(row, error=F(1, 10**6))
        assert not wider['all_distinct_regional_servo_attachments_rejected_at_one_percent']


def test_regional_servo_rejects_invalid_inputs():
    row = r.optical.strict_json.load(r.optical.SUMMARY)['datasets'][0]
    malformed = F(1, 2); malformed._denominator = -2
    for bad_count, bad_width in ((1, 0), (True, 0), (10, -1), (10, True), (10, 1.0)):
        with pytest.raises(ValueError):
            r.servo_interval_mass_bound(5, 5, 1, bad_count, bad_width)
    for bad_replicas in (0, -1, True, 8.0):
        with pytest.raises(ValueError):
            r.servo_and_allan_analysis(row, replicas=bad_replicas)
    for bad_servo in (-1, r.optical.BLOCK_SECONDS, True, 30.0):
        with pytest.raises(ValueError):
            r.servo_and_allan_analysis(row, servo_seconds=bad_servo)
    for bad_frac in (malformed, True, -F(1), 0.5):
        with pytest.raises(ValueError):
            r.servo_and_allan_analysis(row, error=bad_frac)
        with pytest.raises(ValueError):
            r.servo_and_allan_analysis(row, gain_lower=bad_frac)
    for bad_row in ({**row, 'adjacent_complete_block_pairs': 0},
                    {**row, 'nonoverlapping_Allan_variance_diagnostic': None},
                    {**row, 'nonoverlapping_Allan_variance_diagnostic': '0'}):
        with pytest.raises(ValueError):
            r.servo_and_allan_analysis(bad_row)

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

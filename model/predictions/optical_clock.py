"""A conditional public-data test of a direct native-event clock attachment.

A clock is defined internally by its counted events. The tested hypothesis
identifies two such counts with reported optical phase, up to one constant
gain and a bounded readout error. It is deliberately stronger than defining
a clock. Rejection concerns this attachment, not every OPH clock model.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import re

from ..laboratory_clock import event_geometry
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT/'model/predictions/data/optical_clock_summary.json'
SOURCES = {
    '36ar_yb_ratio.txt': '4852780f7d506f811c1bd5eaa52596949d89107d395f6a40c549db97f71e701f',
    '40ar_yb_ratio.txt': '21e297071cacf567fc7ce0f9041fabe21b8a1365c97d93b08d55e93351d5e32d',
}
BLOCK_SECONDS = 300
REFERENCE_HZ = 642121496772645
CALENDAR_START = int(datetime(2021, 1, 1, tzinfo=timezone.utc).timestamp())
CALENDAR_END = int(datetime(2022, 1, 1, tzinfo=timezone.utc).timestamp())
CALENDAR_BLOCKS = (CALENDAR_END-CALENDAR_START)//BLOCK_SECONDS
SOURCE_PATHS = tuple(sorted(set((*event_geometry.SOURCE_PATHS,
    'model/predictions/optical_clock.py', 'model/predictions/__init__.py',
    'model/predictions/data/optical_clock_summary.json',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py'))))


def summarize(content, name):
    """Reduce all published rows; retain only complete fixed UTC-grid blocks.

    No value-based filter is applied. Published removals remain gaps. The
    range test below allows arbitrary removal of whole candidate blocks.
    """
    if name not in SOURCES or type(content) is not bytes:
        raise ValueError('a named public ratio file as original bytes is required')
    lines = content.decode('ascii').splitlines()
    isotope = name[:2]
    if not lines or lines[0] != f'Date Time(UTC) Ratio_{isotope}ar_yb':
        raise ValueError('the published ratio header is required')
    bins, previous, count, first = defaultdict(list), None, 0, None
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != 3 or not re.fullmatch(r'[0-9]+\.[0-9]+', fields[2]):
            raise ValueError('each row needs a UTC second and a finite decimal ratio')
        stamp = fields[0]+'T'+fields[1]
        if not re.fullmatch(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}', stamp):
            raise ValueError('original whole-second UTC timestamps are required')
        second = int(datetime.fromisoformat(stamp).replace(tzinfo=timezone.utc).timestamp())
        ratio = F(fields[2])
        if not CALENDAR_START <= second < CALENDAR_END:
            raise ValueError('the tested acquisition calendar is the fixed year 2021')
        if previous is not None and second <= previous:
            raise ValueError('timestamps must be unique and strictly increasing')
        if not 1 <= ratio <= 2:
            raise ValueError('the declared optical-ratio range is one through two')
        if first is None:
            first = second
        bins[second//BLOCK_SECONDS].append((second, ratio))
        previous, count = second, count+1
    full = {i: sum((r for _, r in rows), F(0))/BLOCK_SECONDS
            for i, rows in bins.items() if len(rows) == BLOCK_SECONDS
            and rows[0][0] == BLOCK_SECONDS*i
            and rows[-1][0] == BLOCK_SECONDS*i+BLOCK_SECONDS-1}
    if len(full) < 2:
        raise ValueError('at least two complete UTC-grid blocks are required')
    adjacent = [(full[i+1]-value)**2 for i, value in full.items() if i+1 in full]
    mean = sum(full.values(), F(0))/len(full)
    return {'file': name, 'sha256': hashlib.sha256(content).hexdigest(),
            'published_rows': count, 'first_UTC_second': first, 'last_UTC_second': previous,
            'block_seconds': BLOCK_SECONDS, 'complete_blocks': len(full),
            'possible_calendar_blocks': max(bins)-min(bins)+1,
            'minimum_block_ratio': str(min(full.values())),
            'maximum_block_ratio': str(max(full.values())),
            'mean_block_ratio': str(mean),
            'adjacent_complete_block_pairs': len(adjacent),
            'nonoverlapping_Allan_variance_diagnostic': (
                str(sum(adjacent, F(0))/(2*len(adjacent)*mean**2)) if adjacent else None),
            'additional_value_based_filter': False}


def verify_raw(directory):
    rows = []
    for name, expected in SOURCES.items():
        content = (Path(directory)/name).read_bytes()
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError('the public dataset bytes differ from the verified source')
        rows.append(summarize(content, name))
    return {'schema': 'oph.optical_clock.public_summary.v1', 'datasets': rows}


def max_point_mass_bound(count):
    """Uniform bound for either native pair: sup_k P[B(count)=k] <= 2/sqrt(count).

    X=Bernoulli(delta/5)+Geom0(5/(10-delta)). The geometric parity bit is
    independent of its even remainder and has success r=q/(1+q), where
    q is 4/9 or 1/2. Thus B(count) is a convolution with Binomial(count,r).
    Fourier inversion, 1-x<=exp(-x), and sin(|t|/2)>=|t|/pi give its maximum
    mass <= sqrt(pi)/(2*sqrt(2*count*r*(1-r))) < 2/sqrt(count).
    Replacing sqrt(count) by its integer floor preserves an exact upper bound.
    """
    if type(count) is not int or count < 2:
        raise ValueError('an original integer count of at least two is required')
    for delta in (0, 1):
        parameters = event_geometry.renewal_parameters(delta)
        p = parameters['geometric_success_probability']
        r = (1-p)/(2-p)
        if (parameters['terminal_B_probability'] != F(delta, 5)
                or p != F(5, 10-delta) or parameters['mean'] != 1
                or parameters['variance'] != 2-F(2*delta, 5) or 8*r*(1-r) <= 1):
            raise ArithmeticError('the native renewal decomposition lost its uniform bound')
    return min(F(1), F(2, isqrt(count)))


def selection_bound(total, retained, interval_mass, *, choice_budget=1):
    """Any retained subset of IID draws, in any interval of the declared width.

    Anchor the interval at its smallest retained value, union over all total
    anchors and (retained-1)-subsets, then use binomial(n,k)<=(3n/k)^k.
    The constant gain/centre may be chosen after seeing the values. A further
    union budget can cover choosing the calendar window and retained count.
    """
    if (type(total) is not int or type(retained) is not int
            or not 2 <= retained <= total or type(interval_mass) is not F
            or not 0 < interval_mass <= 1 or type(choice_budget) is not int or choice_budget < 1):
        raise ValueError('valid integer window counts and exact positive mass are required')
    base = min(F(1), 3*(total-1)*interval_mass/(retained-1))
    power = 0
    while base*2**(power+1) <= 1:
        power += 1
    exponent = power*(retained-1)
    prefactor = total*choice_budget
    # 2**10>10**3; exact integer digit count bounds the prefactor.
    decimal_exponent = min(0, len(str(prefactor))-3*(exponent//10))
    return {'interval_occupancy_union_base_upper': str(base),
            'calendar_and_count_selection_union_budget': choice_budget,
            'binary_power_per_retained_draw': power,
            'probability_upper_form': f'min(1,{prefactor}*2^(-{exponent}))',
            'conservative_probability_upper_power_of_ten': decimal_exponent}


def evaluate(row, *, error=F(1, 10**12), gain_lower=F(1)):
    if type(error) is not F or error < 0 or type(gain_lower) is not F or gain_lower <= 0:
        raise ValueError('exact nonnegative error and positive gain lower bound are required')
    low = (REFERENCE_HZ-1)*BLOCK_SECONDS
    high = (REFERENCE_HZ+1)*BLOCK_SECONDS
    span = F(row['maximum_block_ratio'])-F(row['minimum_block_ratio'])
    if span < 0 or row['block_seconds'] != BLOCK_SECONDS:
        raise ValueError('ordered extrema for the fixed observation window are required')
    width = high*(span+2*error)/gain_lower
    integer_width = -(-width.numerator//width.denominator)
    mass = min(F(1), (integer_width+1)*max_point_mass_bound(low))
    if (not CALENDAR_START <= row['first_UTC_second'] <= row['last_UTC_second'] < CALENDAR_END
            or not 2 <= row['possible_calendar_blocks'] <= CALENDAR_BLOCKS):
        raise ValueError('the candidate calendar window must lie within the fixed acquisition year')
    # At most Y² calendar endpoint pairs and Y retained-count choices. This
    # charges selection of the observed span, start, stop and subset size.
    result = selection_bound(row['possible_calendar_blocks'], row['complete_blocks'], mass,
                             choice_budget=CALENDAR_BLOCKS**3)
    return {'file': row['file'], 'complete_blocks': row['complete_blocks'],
            'possible_calendar_blocks': row['possible_calendar_blocks'],
            'observed_ratio_range': str(span), 'count_interval_width_upper': integer_width,
            'single_window_interval_mass_upper': str(mass), **result,
            'attachment_rejected_at_one_percent_under_stated_hypotheses':
                result['conservative_probability_upper_power_of_ten'] <= -2}


def run():
    summary = strict_json.load(SUMMARY)
    if (set(summary) != {'schema', 'datasets'} or summary['schema'] != 'oph.optical_clock.public_summary.v1'
            or [row['file'] for row in summary['datasets']] != list(SOURCES)
            or any(row['sha256'] != SOURCES[row['file']] for row in summary['datasets'])):
        raise ValueError('the retained summaries must bind both named public datasets')
    return {'schema': 'oph.optical_clock.direct_attachment.v1',
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            'public_dataset': 'https://zenodo.org/records/6901524',
            'measurement_protocol': 'https://arxiv.org/abs/2205.13053',
            'reference_frequency': 'https://arxiv.org/abs/2010.06620',
            'data': summary,
            'hypothesis': {'clock_definition': 'each incident repair ticks, including zero current',
                'comparison': 'count B between successive A ticks, including coincident terminal ticks',
                'scheduler': 'independent uniform native seam selection',
                'optical_attachment': 'reported block ratio = constant_gain*B(N)/N + bounded_error',
                'reference_frequency_Hz_interval': [REFERENCE_HZ-1, REFERENCE_HZ+1],
                'A_ticks_per_block': 'one fixed integer N in the frequency interval times 300 seconds',
                'constant_gain_lower': '1', 'absolute_ratio_error_per_block_upper': '1/1000000000000',
                'arbitrary_missing_or_discarded_blocks_allowed': True,
                'window_endpoints_and_retained_count_selection_charged_within_calendar_year': 2021,
                'searches_outside_declared_calendar_year_covered': False,
                'extra_temporal_filtering_or_collective_clock_mechanism_included': False,
                'bounded_error_and_attachment_are_tested_hypotheses_not_established_physics': True},
            'conditional_tests': [evaluate(row) for row in summary['datasets']],
            'predicted_relative_variance_per_N_A_ticks': {'adjacent': '8/(5*N)', 'nonadjacent': '2/N'},
            'assessment': {'retrospective_conditional_attachment_test': True,
                'OPH_or_every_native_clock_excluded': False,
                'native_atomic_transition_or_interacting_QFT_derived': False,
                'new_distinction_from_SM_GR_established': False}, 'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', type=Path, help='independently re-extract both original public files')
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.raw_dir:
        compare_receipt(strict_json.load(SUMMARY), verify_raw(args.raw_dir))
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Conditional optical event-clock attachment comparison replayed.')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()

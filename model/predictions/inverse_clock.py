"""Visible current correlations predict the same transfer's coherent clock.

Four lag covariances identify two positive simple modes within the declared
two-mode domain. Exact interval propagation permits inconclusive outcomes.
The input tape retains zero-current records at every attempt; it is not an
isolated local-port tape. No physical dataset is identified with that input.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from . import modular_optical as inherited
from ..record_laws import burn_in

source, exact, ROOT = inherited.source, inherited.exact, inherited.ROOT
SOURCE_PATHS = tuple(sorted(set((
    *(p for p in inherited.SOURCE_PATHS if p != inherited.DATA),
    *burn_in.SOURCE_PATHS, 'model/predictions/inverse_clock.py'))))
RATE_UPPER = F(9, 10)


def _integer(value, low=0, high=10**40):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original integer in the stated range is required')
    return value


def _fraction(value, low=None, high=None):
    value = exact._fraction(value)
    if (low is not None and value < low) or (high is not None and value > high):
        raise ValueError('an exact rational in the stated range is required')
    return value


def interval(value, *, probability=False):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError('a pair of exact ordered endpoints is required')
    lo, hi = (_fraction(x, 0 if probability else None, 1 if probability else None) for x in value)
    if lo > hi:
        raise ValueError('interval endpoints are reversed')
    return lo, hi


def rounded(box, digits=18):
    lo, hi = interval(box)
    scale = 10**_integer(digits, 1, 30)
    lo, hi = lo*scale, hi*scale
    return F(lo.numerator//lo.denominator, scale), F(-((-hi).numerator//(-hi).denominator), scale)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def sub(a, b):
    return a[0]-b[1], a[1]-b[0]


def mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def div(a, b):
    if b[0] <= 0 <= b[1]:
        raise ValueError('division interval contains zero')
    return mul(a, (1/b[1], 1/b[0]))


def sqrt_box(a):
    if a[0] < 0:
        raise ValueError('a square-root interval is negative')
    return rounded((source.exact._sqrt_interval(a[0])[0], source.exact._sqrt_interval(a[1])[1]))


def minus_log(a):
    if not F(1, 2) <= a[0] <= a[1] < 1:
        raise ValueError('the certified log regime is [1/2,1)')
    return rounded((-source._log_interval(a[1])[1], -source._log_interval(a[0])[0]))


def identify(mean, joint):
    """Outer enclosure of every compatible positive two-exponential model.

    joint[n] encloses E[B_t B_(t+n+1)], n=0..3. The rank-two spectral form
    is an admitted domain, not a conclusion about arbitrary hidden models.
    Repeated uncertain quantities remain shared in the true model; ordinary
    interval arithmetic only widens their image and never assumes independence.
    """
    mean = interval(mean, probability=True)
    if type(joint) is not tuple or len(joint) != 4:
        raise ValueError('four consecutive lag-product probability intervals are required')
    joint = tuple(interval(x, probability=True) for x in joint)
    covariance = tuple(sub(x, mul(mean, mean)) for x in joint)
    c0, c1, c2, c3 = covariance
    determinant = sub(mul(c0, c2), mul(c1, c1))
    result = {'certified': False,
              'certificate_means': 'outer enclosure of every compatible model in the admitted domain',
              'existence_of_model_in_probability_box_certified': False,
              'covariance_intervals': covariance,
              'Hankel_determinant_interval': determinant}
    def fail(reason):
        return {**result, 'reason': reason}
    if c0[0] <= 0 or determinant[0] <= 0:
        return fail('positive two-mode Hankel determinant is not certified')
    trace = div(sub(mul(c0, c3), mul(c1, c2)), determinant)
    product = div(sub(mul(c1, c3), mul(c2, c2)), determinant)
    discriminant = sub(mul(trace, trace), mul((F(4), F(4)), product))
    if trace[0] <= 0 or product[0] <= 0 or discriminant[0] <= 0:
        return fail('positive simple real roots are not certified')
    separation = sqrt_box(discriminant)
    plus = rounded(div(add(trace, separation), (F(2), F(2))))
    minus = rounded(div(sub(trace, separation), (F(2), F(2))))
    if not F(1, 2) <= minus[0] <= minus[1] < plus[0] <= plus[1] < 1:
        return fail('ordered separated roots inside the certified logarithm regime are not certified')
    weight_plus = rounded(div(sub(c1, mul(minus, c0)), sub(plus, minus)))
    weight_minus = rounded(div(sub(mul(plus, c0), c1), sub(plus, minus)))
    if min(weight_plus[0], weight_minus[0]) <= 0:
        return fail('positive residues are not certified')
    energies = minus_log(plus), minus_log(minus)
    return {**result, 'certified': True,
        'recurrence_trace_interval': rounded(trace), 'recurrence_product_interval': rounded(product),
        'transfer_root_intervals': (plus, minus), 'positive_weight_intervals': (weight_plus, weight_minus),
        'energy_intervals': energies,
        'unit_free_frequency_ratio_interval': rounded(div(energies[1], energies[0])),
        'model_domain': 'two positive exponential modes; no completeness claim for arbitrary hidden processes'}


def _model_boxes(identification):
    """Validate interval inputs, not provenance or existence of a fitted model."""
    if type(identification) is not dict or identification.get('certified') is not True:
        raise ValueError('a certified two-mode enclosure is required')
    result = []
    for key in ('transfer_root_intervals', 'positive_weight_intervals'):
        values = identification.get(key)
        if type(values) is not tuple or len(values) != 2:
            raise ValueError('two exact root and residue intervals are required')
        result.append(tuple(interval(x) for x in values))
    roots, weights = result
    plus, minus = roots
    if not F(1, 2) <= minus[0] <= minus[1] < plus[0] <= plus[1] < 1:
        raise ValueError('ordered separated roots in [1/2,1) are required')
    if min(weight[0] for weight in weights) <= 0:
        raise ValueError('positive residue intervals are required')
    return roots, weights


def predict_joint(mean, identification, lag):
    """Predict an unused lag from the same inferred roots and residues."""
    mean = interval(mean, probability=True)
    _integer(lag, 0, 32)
    roots, weights = _model_boxes(identification)
    covariance = (F(0), F(0))
    for root, weight in zip(roots, weights, strict=True):
        covariance = add(covariance, mul(weight, (root[0]**lag, root[1]**lag)))
    return rounded(add(mul(mean, mean), covariance))


def predict_readout(identification, mode, parameter=4, quadrature='x', *, raw_error=F(1, 10**10)):
    """All four coherent code outcomes, with raw-array event-weight allowance.

    Parameter uses the same logarithm convention as the source compiler.
    A common time unit rescales both energies together, never independently.
    raw_error bounds complete half-diamond distance; an unnormalized output
    can change an event weight by twice this amount.
    """
    _integer(mode, 1, 2)
    _integer(parameter, 0, 8)
    roots, _ = _model_boxes(identification)
    if type(quadrature) is not str or quadrature not in ('x', 'y'):
        raise ValueError('a certified identification and x/y quadrature are required')
    raw_error = _fraction(raw_error, 0, 1)
    # Derive phases from the supplied root boxes again: callers cannot decouple
    # the phase and record scales by replacing a cached energy field.
    energy = minus_log(roots[mode-1])
    phase = source._phase_box(mul((F(parameter), F(parameter)), energy))
    lo, hi = div(add((F(1), F(1)), phase[int(quadrature == 'y')]), (F(2), F(2)))
    lo, hi = max(F(0), lo), min(F(1), hi)
    ideal = tuple((lo, hi) if k == 0 else (1-hi, 1-lo) if k == mode else (F(0), F(0)) for k in range(4))
    return {'ideal_complete_probability_intervals': tuple(rounded(x, 15) for x in ideal),
        'raw_complete_weight_intervals': tuple(rounded((max(F(0), a-2*raw_error), b+2*raw_error), 15)
                                                for a, b in ideal),
        'raw_half_diamond_error': raw_error, 'event_weight_error': 2*raw_error}


def summarize_tape(bits):
    """Read only retained bits, including zeros; all six estimates use M starts."""
    if (type(bits) is not tuple or len(bits) < 6
            or any(type(x) is not int or x not in (0, 1) for x in bits)):
        raise ValueError('a complete tuple of at least six original current bits is required')
    samples = len(bits)-5
    return {'windows': samples, 'mean_ones': sum(bits[:samples]),
            'joint_ones': tuple(sum(bits[t]*bits[t+lag+1] for t in range(samples)) for lag in range(5))}


def variance_coefficient(window_length):
    """Stationary variance(mean of M overlapping bounded windows) <= A_L/M.

    The marked strong lumping has reversible nonconstant L2 rate <9/10.
    Overlap covariances are <=1/4; separated-window covariances are bounded
    by (9/10)^gap/4. No independent-trial or independent-window assumption.
    """
    length = _integer(window_length, 1, 32)
    return (2*length-1+2/(1-RATE_UPPER))/4


def tape_confidence(summary, radius, *, warmup=221):
    """Joint coverage for mean and five products under the declared source law.

    An arbitrary initial binary state is allowed after unread warm-up. The
    complete six-tube failure bound adds the upstream total-variation error.
    The fifth product is available for validation and is not used by identify.
    """
    if type(summary) is not dict or set(summary) != {'windows', 'mean_ones', 'joint_ones'}:
        raise ValueError('exactly the complete tape-count fields are required')
    count = _integer(summary['windows'], 1)
    mean = _integer(summary['mean_ones'], 0, count)
    products = summary['joint_ones']
    if type(products) is not tuple or len(products) != 5:
        raise ValueError('five complete lag counts are required')
    products = tuple(_integer(x, 0, count) for x in products)
    radius = _fraction(radius, 0, 1)
    if radius == 0:
        raise ValueError('a positive confidence radius is required')
    _integer(warmup, 0, 10000)
    values = (F(mean, count), *(F(x, count) for x in products))
    boxes = tuple((max(F(0), x-radius), min(F(1), x+radius)) for x in values)
    coefficient = sum((variance_coefficient(n) for n in range(1, 7)), F(0))
    return {'mean_interval': boxes[0], 'joint_intervals': boxes[1:],
        'stationary_union_variance_coefficient': coefficient,
        'coverage_failure_upper': min(F(1), coefficient/(count*radius**2)+burn_in.binary_bound(warmup)),
        'overlapping_windows_treated_as_independent': False,
        'includes_unused_fifth_lag': True}


def source_moments():
    """Independent full-native marked path sums; no fitted spectral data."""
    kernel = source.transfer.native_counts()
    mean = F(sum(sum(row) for row in kernel.marked_counts[1]), 66*60)
    joint = tuple(source.transfer.exact_record_law(kernel, lag)[1, 1][0] for lag in range(6))
    return mean, joint


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: _json(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    return value


def run():
    mean, joint = source_moments()
    c = tuple(x-mean**2 for x in joint)
    determinant = c[0]*c[2]-c[1]**2
    trace = (c[0]*c[3]-c[1]*c[2])/determinant
    product = (c[1]*c[3]-c[2]**2)/determinant
    if (determinant, trace, product) != (F(1, 9702990000), F(5, 3), F(311, 450)):
        raise ArithmeticError('the record-derived recurrence changed')
    exact_fit = identify((mean, mean), tuple((x, x) for x in joint[:4]))
    if not exact_fit['certified'] or exact_fit['transfer_root_intervals'][0][1] >= RATE_UPPER:
        raise ArithmeticError('the native inverse spectral or rate certificate failed')
    radius = F(1, 10**12)
    mean_box = mean-radius, mean+radius
    tubes = tuple((x-radius, x+radius) for x in joint[:4])
    fitted = identify(mean_box, tubes)
    if not fitted['certified']:
        raise ArithmeticError('the declared source-probability tubes do not identify separated modes')
    if source._sequence_allowance(6) > F(1, 10**10):
        raise ArithmeticError('the complete source compiler exceeds the declared raw-output allowance')
    reads = []
    for mode in (1, 2):
        for quadrature in ('x', 'y'):
            prediction = predict_readout(fitted, mode, 4, quadrature)
            program = source.program(mode, (4,), quadrature=quadrature)
            actual = source.factor_execute(program)
            weights = tuple(branch.probability for branch in actual.branches)
            if len(weights) != 4 or any(not lo <= F(weight) <= hi for weight, (lo, hi) in
                    zip(weights, prediction['raw_complete_weight_intervals'], strict=True)):
                raise ArithmeticError('the shared source readout falls outside the record-derived forecast')
            reads.append({'mode': mode, 'quadrature': quadrature, 'parameter': 4,
                **prediction, 'raw_source_record_weights_diagnostic': weights,
                'all_padding_outcomes_retained': True})
    predicted = predict_joint(mean_box, fitted, 4)
    if not predicted[0] <= joint[4] <= predicted[1] or predicted[0] <= joint[4]+F(1, 1000) <= predicted[1]:
        raise ArithmeticError('the unused-lag positive and altered-data controls failed')
    one_mode_defect = c[2]-c[1]**2/c[0]
    required = sum((variance_coefficient(n) for n in range(1, 7)), F())/(radius**2*F(9, 1000))
    sufficient_windows = -((-required.numerator)//required.denominator)
    return _json({'schema': 'oph.inverse_record_modular_clock.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'visible_input': 'complete binary current tape including every zero-current attempt record',
        'exact_source_mean': mean, 'exact_joint_products': joint, 'exact_covariances': c,
        'Hankel_determinant': determinant, 'recurrence_trace': trace, 'recurrence_product': product,
        'exact_identification': exact_fit,
        'declared_probability_tubes': {'synthetic_source_centered_example': True, 'absolute_radius': radius,
            'mean': mean_box, 'training_joint_products': tubes},
        'uncertain_identification': fitted, 'coherent_readout_forecasts': reads,
        'unused_lag_validation': {'unread_attempts': 4, 'prediction': predicted,
            'actual_source_joint_product': joint[4], 'shift_by_1_over_1000_rejected': True},
        'single_exponential_control': {'matches_first_two_covariances': True,
            'third_covariance_defect': one_mode_defect,
            'claimed_native_alternative_instrument': False},
        'finite_tape_coverage': {'six_statistic_variance_sum': F(39),
            'sufficient_windows_for_declared_radius_and_failure_1_over_100': sufficient_windows,
            'discarded_binary_warmup_attempts': 221,
            'bound_is_sufficient_not_minimal_sample_requirement': True,
            'independent_windows_assumed': False},
        'scope': {'positive_rank_two_spectral_domain_admitted': True,
            'unknown_hidden_dimension_identified_by_four_moments': False,
            'native_C_E_instrument_identified_from_classical_bits': False,
            'coherent_code_preparation_and_readout_are_declared_services': True,
            'different_clock_or_field_frequency_scales_fitted': False,
            'matched_physical_record_dataset_available': False,
            'new_discriminator_against_SM_GR_established': False}, 'all_pass': True})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        source.compare_receipt(result, source.strict_json.load(args.check))
        print('Visible-record inverse clock and complete coherent forecasts replayed.')
    else:
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

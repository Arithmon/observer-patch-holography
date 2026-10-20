"""Visible current correlations predict the same transfer's coherent clock.

Four lag covariances identify two positive simple modes within the declared
two-mode domain. Exact interval propagation permits inconclusive outcomes.
The input tape retains zero-current records at every attempt; it is not an
isolated local-port tape. No physical dataset is identified with that input.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from . import modular_optical as inherited
from ..observer import coherent_continuation as continuation
from ..record_laws import burn_in

source, exact, ROOT = inherited.source, inherited.exact, inherited.ROOT
SOURCE_PATHS = tuple(sorted(set((
    *(p for p in inherited.SOURCE_PATHS if p != inherited.DATA),
    *burn_in.SOURCE_PATHS, *continuation.SOURCE_PATHS,
    'model/predictions/inverse_clock.py'))))
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


def edge_gibbs_bound(parameter, heralds):
    """Ideal lazy-filter Gibbs limit, in a preparation parameter, not seconds.

    lambda_max=5+sqrt(5)<8 and p=30*t/n<=1. For x=t*lambda/(2*n),
    0<=-2*n*log(1-x)-t*lambda<=t²*lambda²/[4*n*(1-x)].
    Thus finite weights are the Gibbs weights times factors in [exp(-e),1].
    Rejection coupling gives TV<=min(1,e); normalized log ratios give D<=e.
    These are analytic bounds, not numerical source-pulse error certificates.
    """
    t = _fraction(parameter, 0)
    n = _integer(heralds, 1)
    p = 30*t/n
    if p > 1:
        raise ValueError('the lazy event preparation requires 30*t/n <= 1')
    error = 16*t*t/(n*(1-4*t/n))
    return {'preparation_parameter': t, 'successful_heralds': n, 'lazy_weight': p,
            'spectral_log_error_upper': error,
            'trace_distance_to_edge_Gibbs_upper': min(F(1), error),
            'relative_entropy_to_edge_Gibbs_upper': error,
            'free_entropy_deficit_upper': error,
            'all_success_probability_lower': F(1, 12),
            'physical_seconds_or_modular_time': False}


def _lazy_filter_program(weight):
    """Source-only dilation of (1-p)I+pK12; its event bias is supplied."""
    import math
    import mpmath as mp
    import numpy as np
    from ..observer.source_scalar import compiler as c
    p = _fraction(weight, 0, 1)
    program = continuation.density_walk_program()
    gray = [j ^ (j >> 1) for j in range(64)]
    counts = [((p.numerator if j & 1 else 2*p.denominator-p.numerator)
               if j < 60 else 0) for j in gray]
    rotations = []
    with mp.workdps(80):
        for row in range(63, 0, -1):
            x, y = counts[row-1], counts[row]
            theta = 0. if x == y == 0 else float(mp.atan2(mp.sqrt(y), mp.sqrt(x)))
            rotations.append((gray[row-1], gray[row], theta))
            counts[row-1], counts[row] = x+y, 0
    prepare, target = [], np.eye(64)
    for a, b, theta in reversed(rotations):
        prepare.extend(c._gray_rotation(a, b, -theta, 6))
        x, y = target[a].copy(), target[b].copy()
        target[a] = math.cos(theta)*x-math.sin(theta)*y
        target[b] = math.sin(theta)*x+math.cos(theta)*y
    prepare = [c.pulse(g['kind'], round(g['angle']*2**44)/2**44, g['targets'])
               for g in prepare]
    program['routines'] = {name: program['routines'][name] for name in ('ccx', 'flip')}
    program['routines'].update(prepare=prepare, unprepare=c.inverse(prepare))
    program['walk'] = ([{'routine': 'prepare', 'targets': list(range(4, 10))}]
                       + program['shift']
                       + [{'routine': 'unprepare', 'targets': list(range(4, 10))}])
    return program, target


@lru_cache(maxsize=1)
def prepared_edge_gibbs():
    """Execute one full filter on all twelve logical input columns.

    Reuse that raw executed linear map for two heralds, retaining every
    first-failure output. The supplied input is I12/12; it is not the native
    66-state stationary density. Exact rational matrices independently check
    the success operator, sector ranks, coherence and all outcome weights.
    """
    import math
    import numpy as np
    import sympy as sp
    from ..observer.source_scalar import compiler as c
    program, ideal_prepare = _lazy_filter_program(F(1, 2))
    permutation = continuation.density_shift_permutation(program)
    edge = continuation.edge_generator_audit()
    laplacian = sp.zeros(12)
    for a, b in continuation.quantum_mixing.seams():
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    lazy = sp.eye(12)-laplacian/120
    if laplacian.eigenvals() != {sp.Integer(0): 1, 5-sp.sqrt(5): 3,
                                sp.Integer(6): 5, 5+sp.sqrt(5): 3}:
        raise ArithmeticError('native edge-generator spectrum changed')
    exact_success = [(lazy**(2*n)).trace()/12 for n in range(3)]
    exact_outcomes = (1-exact_success[1], exact_success[1]-exact_success[2], exact_success[2])
    rho = lazy**4/(12*exact_success[2])
    if rho[0, 1] <= 0 or (rho*rho).trace() <= sp.Rational(1, 12):
        raise ArithmeticError('the heralded density lost its nontracial coherence')
    matrices = {name: c.execute(tape, 6 if name.endswith('prepare') else
                                3 if name == 'ccx' else 1)
                for name, tape in program['routines'].items()}
    target_ccx = np.zeros((8, 8))
    for j in range(8):
        target_ccx[j ^ (4 if j & 3 == 3 else 0), j] = 1
    wanted = {'prepare': ideal_prepare, 'unprepare': ideal_prepare.T,
              'ccx': target_ccx, 'flip': np.array([[0., 1.], [1., 0.]])}
    guards = {'prepare': F(1, 10**8), 'unprepare': F(1, 10**8),
              'ccx': F(1, 10**10), 'flip': F(1, 10**10)}
    for name, matrix in matrices.items():
        if np.linalg.norm(matrix-wanted[name], 'fro') > float(guards[name]):
            raise ArithmeticError(f'raw source {name} exceeds its numerical guard')
    amplitudes = np.array([math.sqrt(float((F(1, 2) if j & 1 else F(3, 2))/60))
                           if j < 60 else 0. for j in range(64)])
    if np.linalg.norm(ideal_prepare[:, 0]-amplitudes) > 1e-13:
        raise ArithmeticError('biased event preparation has incorrect amplitudes')
    inputs = np.zeros((1 << 15, 12), complex)
    inputs[:12] = np.eye(12)
    actual = continuation.execute_density_walk(program, matrices, inputs)
    ideal = c._apply(inputs, ideal_prepare, list(range(9, 3, -1)), 15)
    shifted = np.empty_like(ideal); shifted[permutation] = ideal
    ideal = c._apply(shifted, ideal_prepare.T, list(range(9, 3, -1)), 15)
    wanted_filter = np.array(lazy, dtype=float)
    if np.linalg.norm(ideal[:12]-wanted_filter, 'fro') > 1e-12:
        raise ArithmeticError('source dilation does not compress to the native lazy transfer')
    uses = {name: sum(op['routine'] == name for op in program['walk']) for name in matrices}
    guard_sum = sum((uses[name]*guards[name] for name in uses), F())
    delta = guard_sum/(1-guard_sum)
    if np.linalg.norm(actual-ideal, 'fro') > math.sqrt(12)*float(delta):
        raise ArithmeticError('complete source filter output exceeds its product bound')
    # Every nonherald row is retained, including site/event padding and dirty
    # scratch. No normalization occurs before the complete outcome checks.
    first = actual/math.sqrt(12)
    second = actual@actual[:12]/math.sqrt(12)
    outcomes = (np.linalg.norm(first[12:])**2, np.linalg.norm(second[12:])**2,
                np.linalg.norm(second[:12])**2)
    two_delta = 2*delta+delta*delta
    weight_guard = 2*two_delta+two_delta*two_delta
    if (any(abs(value-float(expected)) > float(weight_guard)
            for value, expected in zip(outcomes, exact_outcomes, strict=True))
            or abs(sum(outcomes)-1) > 3*float(weight_guard)
            or np.linalg.norm(second-ideal@wanted_filter/math.sqrt(12), 'fro') > float(two_delta)):
        raise ArithmeticError('complete two-herald source outcomes exceed their guards')
    read_word = [c.pulse(g['kind'], round(g['angle']*2**44)/2**44, g['targets'])
                 for g in c._gray_rotation(0, 1, math.pi/4, 4)]
    read_matrix = c.execute(read_word, 4)
    read_target = np.eye(16)
    read_target[:2, :2] = np.array([[1., 1.], [-1., 1.]])/math.sqrt(2)
    read_guard = F(1, 10**9)
    if np.linalg.norm(read_matrix-read_target, 'fro') > float(read_guard):
        raise ArithmeticError('source coherence readout exceeds its numerical guard')
    successful = np.zeros((16, 12), complex)
    successful[:12] = second[:12]
    read_weights = np.sum(abs(read_matrix@successful)**2, axis=1)
    unnormalized = lazy**4/12
    wanted_reads = [unnormalized[j, j] for j in range(12)]+[sp.Integer(0)]*4
    wanted_reads[:2] = [unnormalized[0, 0]+unnormalized[0, 1],
                       unnormalized[0, 0]-unnormalized[0, 1]]
    if unnormalized[0, 0] != unnormalized[1, 1]:
        raise ArithmeticError('the exact transitive site symmetry changed')
    read_delta = (1+two_delta)*(1+read_guard)-1
    read_weight_guard = 2*read_delta+read_delta*read_delta
    if any(abs(value-float(expected)) > float(read_weight_guard)
           for value, expected in zip(read_weights, wanted_reads, strict=True)):
        raise ArithmeticError('retained source coherence-readout weights changed')
    word = continuation.density_walk_word(program)
    resources = c.resources(word)
    stable_resources = {key: resources[key] for key in (
        'native_pulses', 'single_code_pulses', 'pair_code_pulses', 'identity_phase_pulses',
        'proper_code_load_events', 'proper_code_export_events', 'proper_code_service_events',
        'maximum_work_qubits', 'service_failure_outcomes_retained',
        'valid_code_failure_probability_in_exact_arithmetic')}
    as_fraction = lambda x: F(int(sp.numer(x)), int(sp.denom(x)))
    return {'native_density_generator': 'L=6I-(1/2)sum_{g in C} U_g; K12=I-L/60',
        'sector_dimensions': [1, 3, 5, 3],
        'sector_eigenvalues': ['0', '5-sqrt(5)', '6', '5+sqrt(5)'],
        'representation_algebra_dimension': edge['group_image_algebra_dimension'],
        'prepared_register_algebra': 'B(C^12); includes operators between representation sectors',
        'restricted_group_image_modular_flow': 'identity: exp(-tL)/Z is central in M1+M3+M5+M3',
        'event_probability_even': F(1, 40), 'event_probability_odd': F(1, 120),
        'one_herald_operator': 'F=(1-p)I+pK12=I-L/120, p=1/2',
        'success_herald': 'valid site address 0..11, event zero and clean scratch',
        'initial_density': 'supplied I12/12 on a protected logical register',
        'successful_heralds': 2, 'executed_full_register_input_columns': 12,
        'retained_first_failure_and_success_probabilities': tuple(map(as_fraction, exact_outcomes)),
        'normalized_success_density': 'F^4/Tr(F^4)',
        'normalized_success_matrix_element_0_1': as_fraction(rho[0, 1]),
        'normalized_success_purity': as_fraction((rho*rho).trace()),
        'coherence_readout': {
            'operation': 'exp(i*pi*Y_01/4) on successful site register, then all 16 address outcomes',
            'ideal_unnormalized_address_probabilities': tuple(map(as_fraction, wanted_reads)),
            'normalized_zero_minus_one_contrast': as_fraction(2*rho[0, 1]),
            'dephased_input_contrast': F(0),
            'source_native_pulses': len(read_word),
            'source_word_sha256': hashlib.sha256(json.dumps(read_word, sort_keys=True,
                separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
            'complete_outcome_weight_error_upper': read_weight_guard,
            'sector_projector_detector_compiled': False},
        'event_or_scratch_failure_outputs_discarded_from_total_probability': False,
        'compiled_word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True,
            separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
        'one_filter_resources': stable_resources, 'maximum_filters_per_trial': 2,
        'routine_error_guards': guards, 'one_filter_operator_error_upper': delta,
        'two_filter_outcome_weight_error_upper': weight_guard,
        'raw_source_routines_executed_on_all_columns': True,
        'execution': 'raw routine matrices reused on all 32768 outputs; full filter reused for second herald',
        'expanded_word_executed_literally_on_full_register': False,
        'numerical_guards_are_interval_or_physical_error_certificates': False,
        'ideal_family': edge_gibbs_bound(F(1, 30), 2),
        'sector_law': 'p_R=d_R*(1-t*lambda_R/(2n))^(2n)/sum_S d_S*(1-t*lambda_S/(2n))^(2n)',
        'limit_sector_ratio': '(p_R/d_R)/(p_S/d_S)=exp[-t*(lambda_R-lambda_S)]',
        'limit_free_entropy_identity': 'S(rho)-t*Tr(rho L)=log Z-D(rho||exp(-tL)/Z)',
        'dephasing_control': 'diagonalizing the heralded density preserves site probabilities but removes rho[0,1]>0',
        'edge_and_wave_modular_generators_identified': False,
        'full_group_regular_or_compact_group_sewing_established': False,
        'native_passive_equilibrium_prepares_this_state': False,
        'growing_source_compiler_error_budget_certified': False,
        'physical_dataset_with_this_joint_preparation_and_readout': False}


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
        'actively_prepared_edge_reference': prepared_edge_gibbs(),
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

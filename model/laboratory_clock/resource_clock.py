"""Exact shared timing calibration for two distinct operation counters.

Native attempts and compiled source-pulse invocations need not take equal
physical durations. The affine timing law is an explicit external hypothesis;
its finite calibration never identifies a native event with an atomic clock.
"""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache

from ..quantum_observer_correspondence import clock_response as oracle
from ..quantum_observer_correspondence import construction as target
from ..compiled_clock_history import receipt as compiled
from ..repair_identification import core as exact_core


def _rational(value):
    """Reject corrupted rational slots and detach caller-owned values."""
    if type(value) is not F:
        raise ValueError('canonical rational input required')
    return exact_core._fraction(value)


@dataclass(frozen=True)
class ResourceEndpoint:
    name: str
    attempts: int
    pulses: int
    seconds_lo: F
    seconds_hi: F

    def __post_init__(self):
        lo, hi = _rational(self.seconds_lo), _rational(self.seconds_hi)
        if (type(self.name) is not str or not self.name
                or any(type(v) is not int or v < 0 for v in (self.attempts, self.pulses))
                or lo > hi):
            raise ValueError('named exact operation counters and rational timestamp bounds required')
        object.__setattr__(self, 'seconds_lo', lo)
        object.__setattr__(self, 'seconds_hi', hi)


def inequalities(points):
    """Eliminate the common origin; return a*tau+b*kappa <= c.

    Closed nonnegative resource durations are used for the outer polytope.
    Strict-positive feasibility and a positive lower bound are distinguished.
    """
    points = tuple(points)
    if (not points or any(type(p) is not ResourceEndpoint for p in points)
            or len({p.name for p in points}) != len(points) or len(points) > 32):
        raise ValueError('one through 32 distinct validated endpoints required')
    points = tuple(ResourceEndpoint(p.name, p.attempts, p.pulses, p.seconds_lo, p.seconds_hi)
                   for p in points)
    rows = {(F(-1), F(0), F(0)), (F(0), F(-1), F(0))}
    for left in points:
        for right in points:
            rows.add((F(right.attempts-left.attempts), F(right.pulses-left.pulses),
                      right.seconds_hi-left.seconds_lo))
    return tuple(sorted(rows))


def linear_range(rows, attempts=F(1), pulses=F(0)):
    """Exact infimum/supremum of a resource duration over the closed region.

    Substitute y=attempts*tau+pulses*kappa and eliminate the other duration
    with Fourier--Motzkin inequalities. No independent marginal box is used.
    None denotes an unbounded endpoint, never a fitted replacement value.
    """
    attempts, pulses = _rational(attempts), _rational(pulses)
    rows = tuple(tuple(_rational(value) for value in row) for row in rows)
    if any(len(row) != 3 for row in rows):
        raise ValueError('three exact rational inequality coefficients required')
    if not attempts and not pulses:
        feasible, _, _ = linear_range(rows)
        return feasible, F(0), F(0)
    if not attempts:
        rows = tuple((b, a, c) for a, b, c in rows)
        attempts, pulses = pulses, attempts
    positive, negative, reduced = [], [], []
    for a, b, c in rows:
        y, other = a/attempts, b-a*pulses/attempts
        if other > 0:
            positive.append((y/other, c/other))
        elif other < 0:
            negative.append((y/(-other), c/(-other)))
        else:
            reduced.append((y, c))
    reduced.extend((a+b, c+d) for a, c in positive for b, d in negative)
    lower = upper = None
    for coefficient, bound in reduced:
        if coefficient > 0:
            value = bound/coefficient
            upper = value if upper is None else min(upper, value)
        elif coefficient < 0:
            value = bound/coefficient
            lower = value if lower is None else max(lower, value)
        elif bound < 0:
            return False, None, None
    return (not (lower is not None and upper is not None and lower > upper), lower, upper)


def calibration(points):
    points = tuple(points)
    rows = inequalities(points)
    feasible, native_lo, native_hi = linear_range(rows)
    _, pulse_lo, pulse_hi = linear_range(rows, F(0), F(1))
    positive = feasible and (native_hi is None or native_hi > 0) and (pulse_hi is None or pulse_hi > 0)
    identified = (positive and native_lo is not None and native_lo > 0
                  and pulse_lo is not None and pulse_lo > 0
                  and native_hi is not None and pulse_hi is not None)
    return {'closed_outer_feasible': feasible, 'strict_positive_durations_feasible': positive,
            'both_durations_bounded_and_bounded_away_from_zero': identified,
            'native_seconds_bounds': [native_lo, native_hi] if feasible else None,
            'pulse_seconds_bounds': [pulse_lo, pulse_hi] if feasible else None,
            'calibration_ids': [point.name for point in points],
            'one_common_origin_and_duration_pair': True,
            'physical_timing_hypothesis_authenticated': False}


def validate(calibration_points, validation_points):
    calibration_points, validation_points = tuple(calibration_points), tuple(validation_points)
    if len(calibration_points) < 2 or not validation_points:
        raise ValueError('at least two calibration and one validation endpoints required')
    inequalities(calibration_points+validation_points)  # Also checks split disjointness.
    return {'calibration': calibration(calibration_points),
            'validation_ids': [p.name for p in validation_points],
            'held_out_used_in_calibration': False,
            'individually_compatible': [calibration(calibration_points+(p,))[
                'strict_positive_durations_feasible'] for p in validation_points],
            'jointly_compatible_with_some_original_calibration': calibration(
                calibration_points+validation_points)['strict_positive_durations_feasible'],
            'replacement_fitted_duration_pair_exported': False}


def _interval(state, cap):
    waiting, marked = dict(state), []
    for attempt in range(1, cap+1):
        step = oracle.native_step(waiting)
        marked.extend((attempt, mark, rho) for mark, rho in step.items() if mark is not None)
        waiting = step[None]
    return marked, waiting


@lru_cache(maxsize=6)
def _counted_history(gamma, cap):
    """Native endpoint replay with previously hidden scheduled counts retained.

    Values are full unnormalized source matrices. The old observer output is
    recovered exactly by forgetting the count, without assigning it a mean.
    """
    result = {}

    def finish(count, word, status, state):
        for bit in (0, 1):
            result[count, word, status, bit] = {
                key: value for key, value in state.items()
                if oracle.STATES[key[0]][1] == bit == oracle.STATES[key[1]][1]}

    first, timeout = _interval(oracle.preparation(), cap)
    finish(cap, (), 'timeout_first', timeout)
    for first_count, mark, state in first:
        seam = (0, 1) if mark[1][0] > 0 else (0, 5)
        for feedback, controlled in oracle.native_step(state, seam).items():
            word = (('scheduled', mark), ('feedback', feedback))
            second, timeout = _interval(oracle.apparatus(controlled, gamma), cap)
            finish(first_count+1+cap, word, 'timeout_second', timeout)
            for second_count, last, final in second:
                # Different first/second splits can have the same total count
                # and visible word; sum those CP alternatives below.
                key_count = first_count+1+second_count
                key_word = word+(('scheduled', last),)
                for bit in (0, 1):
                    key = key_count, key_word, 'complete', bit
                    output = result.setdefault(key, {})
                    for entry, value in final.items():
                        if oracle.STATES[entry[0]][1] == bit == oracle.STATES[entry[1]][1]:
                            output[entry] = output.get(entry, F(0))+value
    return result


def counted_history(gamma=0, cap=2):
    """Validated, copy-safe counted source output, including every failure."""
    if type(gamma) is not int or gamma not in (0, 1) or type(cap) is not int or not 1 <= cap <= 3:
        raise ValueError('binary setting and cap one through three required')
    return {key: {entry: _rational(value) for entry, value in state.items()}
            for key, state in _counted_history(gamma, cap).items()}


def resource_law(gamma=0, cap=2):
    costs = compiled.costs(cap=cap)
    result, collapsed = defaultdict(F), {}
    for (attempts, word, status, bit), state in counted_history(gamma, cap).items():
        output = collapsed.setdefault((word, status, bit), {})
        for entry, value in state.items():
            output[entry] = output.get(entry, F(0))+value
        pulses = costs['preparation_pulses']+costs['terminal_pointer_pulses']
        if status != 'timeout_first':
            pulses += costs['conditional_intervention_pulses']
        a, b = (sum(mark[0][i] for _, mark in word) for i in (0, 1))
        result[attempts, pulses, status, a, b, bit] += oracle.trace(state)
    collapsed = {key: {entry: value for entry, value in state.items() if value}
                 for key, state in collapsed.items()}
    if collapsed != target.execute(gamma, cap=cap):
        raise ArithmeticError('forgetting native costs changed a complete conditional state')
    if sum(result.values(), F(0)) != 1:
        raise ArithmeticError('ideal complete resource law must be normalized')
    return dict(sorted(result.items()))


def duration_range(points, attempts, pulses):
    attempts, pulses = _rational(attempts), _rational(pulses)
    if any(value < 0 for value in (attempts, pulses)):
        raise ValueError('nonnegative rational resource exposure required')
    feasible, lower, upper = linear_range(inequalities(points), attempts, pulses)
    if not feasible:
        raise ValueError('incompatible timing calibration')
    return lower, upper


def deadline(law, points, seconds, *, timestamp_error=F(0), detector=None):
    """Outer event probability from one common calibrated timing region.

    Complete ideal probabilities only: raw compiler weights are not silently
    normalized. Time is elapsed from the declared protocol start: the common
    origin cancels. Errors bound the elapsed reading, including both endpoints.
    """
    points = tuple(points)
    seconds, timestamp_error = _rational(seconds), _rational(timestamp_error)
    if type(law) is dict:
        law = {key: _rational(value) for key, value in law.items()}
    if (type(law) is not dict or seconds < 0
            or timestamp_error < 0 or (detector is not None and (type(detector) is not int or detector not in (0, 1)))
            or any(type(key) is not tuple or len(key) != 6
                   or any(type(key[i]) is not int or key[i] < 0 for i in (0, 1, 3, 4))
                   or type(key[2]) is not str or type(key[5]) is not int or key[5] not in (0, 1)
                   for key in law)
            or any(type(value) is not F or value < 0 for value in law.values())
            or sum(law.values(), F(0)) != 1):
        raise ValueError('normalized exact law, nonnegative seconds/errors and declared detector event required')
    if not calibration(points)['strict_positive_durations_feasible']:
        raise ValueError('a compatible positive timing calibration is required')
    lower = upper = F(0)
    costs = sorted({key[:2] for key, value in law.items() if value})
    ordered = all(n <= m and p <= q for (n, p), (m, q) in zip(costs, costs[1:]))
    ranges = {cost: duration_range(points, F(cost[0]), F(cost[1])) for cost in costs}
    for (n, p, status, a, b, bit), value in law.items():
        if not value or (detector is not None and bit != detector):
            continue
        lo, hi = ranges[n, p]
        if hi is not None and hi <= seconds-timestamp_error:
            lower += value
        if lo is None or lo <= seconds+timestamp_error:
            upper += value
    return {'deadline_seconds': seconds, 'timestamp_error_bound': timestamp_error,
            'detector': detector, 'event_probability_bounds': [lower, upper],
            'calibration_used_jointly': True,
            'bounds_sharp_for_ordered_cost_support': ordered and calibration(points)[
                'both_durations_bounded_and_bounded_away_from_zero'],
            'timeouts_retained': True, 'raw_compiler_weights_normalized': False}


def run():
    """Source-derived resource law and explicitly synthetic timing controls."""
    point = ResourceEndpoint
    calibration_points = (point('origin', 0, 0, F(0), F(0)),
        point('native', 100, 0, F(1, 10)-F(1, 10**6), F(1, 10)+F(1, 10**6)),
        point('mixed', 100, 100, F(51, 500)-F(1, 10**6), F(51, 500)+F(1, 10**6)))
    held_out = (point('validation', 200, 50, F(201, 1000)-F(1, 10**6), F(201, 1000)+F(1, 10**6)),)
    changed = (point('validation_changed', 200, 50, F(201, 1000)+F(9, 10**6),
                     F(201, 1000)+F(11, 10**6)),)
    laws = [resource_law(gamma) for gamma in (0, 1)]
    marginals = []
    for law in laws:
        marginal = defaultdict(F)
        for (n, p, status, a, b, bit), value in law.items():
            marginal[n, p] += value
        marginals.append(dict(marginal))
    if marginals[0] != marginals[1]:
        raise ArithmeticError('the coherent setting unexpectedly changed resource exposure')
    marginal = marginals[0]
    expected_n = sum((F(n)*value for (n, p), value in marginal.items()), F(0))
    expected_p = sum((F(p)*value for (n, p), value in marginal.items()), F(0))
    if (expected_n, expected_p) != (F(408, 125), F(1059, 25)):
        raise ArithmeticError('source history resource expectations changed')
    fit = calibration(calibration_points)
    actual_mean = duration_range(calibration_points, expected_n, expected_p)
    nlo, nhi = fit['native_seconds_bounds']
    plo, phi = fit['pulse_seconds_bounds']
    independent_box = (expected_n*nlo+expected_p*plo, expected_n*nhi+expected_p*phi)
    deficient = (point('d0', 0, 0, F(0), F(0)), point('d1', 100, 200, F(13, 125), F(13, 125)))
    q = F(16, 25)
    exact = (point('e0', 0, 0, F(0), F(0)), point('e1', 100, 0, F(1, 10), F(1, 10)),
             point('e2', 100, 100, F(51, 500), F(51, 500)))
    response = oracle.total_variation(laws[0], laws[1])
    return serial({'schema': 'oph.two_resource_clock.v1',
        'timing_hypothesis': 't=b+tau_native*N+tau_pulse*P; counts include every declared operation',
        'source_duration_selected': False, 'data_kind': 'synthetic_fixture',
        'calibration': fit, 'held_out': validate(calibration_points, held_out),
        'held_out_changed': validate(calibration_points, changed),
        'deficient_fixed_resource_ratio_calibration': calibration(deficient),
        'deficient_held_out_native_only_duration_bounds': duration_range(deficient, F(100), F(0)),
        'complete_resource_laws': [[{'attempts': n, 'pulses': p, 'status': status,
            'clock_ticks': [a, b], 'detector': bit, 'probability': probability}
            for (n, p, status, a, b, bit), probability in law.items()] for law in laws],
        'resource_marginal': [{'attempts': n, 'pulses': p, 'probability': mass}
                              for (n, p), mass in sorted(marginal.items())],
        'same_resource_marginal_both_settings': True,
        'complete_conditional_matrices_coarse_grain_to_native_oracle': True,
        'mean_native_attempts_including_feedback': expected_n,
        'mean_source_pulses': expected_p,
        'expected_mean_duration_seconds_bounds': actual_mean,
        'independent_marginal_box_mean_duration_seconds_bounds': independent_box,
        'joint_resource_clock_detector_setting_response_tv': response,
        'calibrated_deadline_predictions': [deadline(law, calibration_points, F(9, 2000),
             timestamp_error=F(1, 10**6), detector=bit)
             for law in laws for bit in (None, 1)],
        'mean_replacement_negative_control': {
            'exact_cadence_deadline': deadline(laws[0], exact, F(1, 250)),
            'probability_from_deterministic_mean_duration': F(int(
                expected_n*F(1, 1000)+expected_p*F(1, 50000) <= F(1, 250))),
            'first_timeout_probability': 1-q,
            'total_timeout_probability': 1-q*q},
        'scope': {'all_resource_cost_classes_and_timeouts_retained': True,
            'additional_attempt_counter_is_declared_readout': True,
            'raw_compiled_history_refinement_executed': False,
            'coarse_current_calibration_bound_transferred_to_resource_records': False,
            'physical_service_duration_or_fidelity_established': False,
            'compiled_raw_array_mass_used_as_probability': False,
            'primitive_source_pulses_are_native_incident_ticks': False,
            'arbitrary_modular_or_atomic_clock_identification': False,
            'unrecorded_service_latency_must_be_bounded_or_added_as_a_resource': True}})


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value

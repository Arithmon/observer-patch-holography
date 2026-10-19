"""One native history gives jointly coupled continuum observation paths.

The time parameter is rescaled cumulative observation count. The field is a
long-observation Gaussian limit, not instantaneous matter or a wave amplitude.
The proof and access contract are in OCCUPATION_FIELD.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as sp

from . import occupation_field as occupation
from . import continuum_field as continuum
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*occupation.SOURCE_PATHS, *continuum.SOURCE_PATHS,
    'model/physical_space_time/observation_process.py'))))


def _integer(value, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError('an original integer in the admitted range is required')
    return value


def _positive(value):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.numerator <= 0
            or value.denominator <= 0 or math.gcd(value.numerator, value.denominator) != 1):
        raise ValueError('a positive canonical Fraction is required')
    return value


def clock_difference(level):
    """Exact covariance coefficient of the normalized same-prefix difference."""
    scale = 4**_integer(level)
    m, seams = 33, 630*scale
    probability = F(m, seams)
    return (1-probability)/(4*m*scale)


def path_bounds(level, *, horizon=F(1), count=None, degree=0):
    """All-level path bounds, requiring no allocation of a refined graph.

    A count supplies the exact finite-prefix reverse-martingale bound. Without
    it, the joint fixed-graph functional CLT has already been taken.
    """
    level, degree = _integer(level), _integer(degree)
    horizon = _positive(horizon)
    difference = clock_difference(level)
    result = {'joint_Gaussian_clock_sup_squared_upper': 8*horizon*difference,
        'Gaussian_harmonic_sup_tail_squared_upper': 4*horizon*continuum.harmonic_tail(degree),
        'Gaussian_fourth_increment_coefficient_upper': 3*continuum.harmonic_tail(0)**2,
        'clock_difference_covariance_coefficient': difference}
    if count is not None:
        _integer(count, 1)
        endpoint = (count*horizon.numerator+horizon.denominator-1)//horizon.denominator
        result['finite_prefix_clock_sup_squared_upper'] = 18*F(endpoint, count)*difference
    return result


@lru_cache(maxsize=1)
def marked_source():
    """Every native endpoint branch, with stationary marking verified exactly."""
    states, remote, incident, z, checked = occupation.native_counts()
    seams, _ = occupation.passage.recovered_carrier()
    index = {state: i for i, state in enumerate(states)}
    targets = np.empty((66, 60), dtype=np.int64)
    marks = np.array([int(0 in edge) for edge in seams for _ in (0, 1)], dtype=np.int64)
    for i, state in enumerate(states):
        for e, edge in enumerate(seams):
            for coin in (0, 1):
                final, (_, current) = native.native(state, edge, coin)
                decoded = tuple(state[p]+current*(int(p == edge[1])-int(p == edge[0]))
                                for p in range(12))
                if decoded != final:
                    raise ArithmeticError('the complete current tape lost field readback')
                targets[i, 2*e+coin] = index[final]
        # The instrument on each seam is exactly half identity plus half swap.
        for e, (a, b) in enumerate(seams):
            swapped = list(state); swapped[a], swapped[b] = swapped[b], swapped[a]
            if sorted(targets[i, 2*e:2*e+2]) != sorted((i, index[tuple(swapped)])):
                raise ArithmeticError('a native seam no longer preserves the uniform sector')
    uniform = sp.ones(66, 1)
    signed = 5*incident-remote
    if signed*uniform != sp.zeros(66, 1) or sp.ones(1, 66)*signed != sp.zeros(1, 66):
        raise ArithmeticError('a stationary marked output has nonzero conditional mean')
    score = np.array([[6*x-1 for x in state] for state in states], dtype=np.int64)
    if sp.Matrix(score.T) != 6*z:
        raise ArithmeticError('the integer current-decoded score changed')
    for array in (targets, marks, score):
        array.flags.writeable = False
    return states, targets, marks, score, checked


def finite_covariances(count):
    """Exact covariances of attempt sum A and clock-difference sum R.

    R=sum_j (I_j/p-1) eta_j uses post-repair eta_j; p=1/6.
    These formulas are valid for every count; this exact executable is bounded.
    """
    _integer(count, 1, 16)
    source = occupation.field('attempt')
    p = source['kernel']
    _, _, incident, _, _ = occupation.native_counts()
    z = sp.Matrix(marked_source()[3].T)/6
    # Density-mark matrix B/p-P is reconstructed from the complete state map.
    local = source['laplacian']-occupation.field('incident')['remote_laplacian']
    b = sp.eye(12)/6-local/60
    r = 6*b-p
    if z*(incident/10-(occupation.native_counts()[1]+incident)/60) != r*z:
        raise ArithmeticError('marked density closure disagrees with full native configurations')
    cross, power = sp.zeros(12), sp.eye(12)
    for lag in range(1, count):
        cross += (count-lag)*power*r
        power = power*p
    return tuple(sp.ImmutableMatrix(x) for x in (
        count*occupation.covariance_at_count('attempt', count),
        count*5*occupation.CHI*occupation.P0,
        occupation.CHI*cross))


def direct_two_steps(*, sticky=False, stationary=True):
    """Independent full endpoint-path moments, with no density-transfer use.

    sticky=True reuses the first seam at the second attempt. Each marginal
    remains stationary, but the IID-clock theorem need not survive.
    """
    if type(sticky) is not bool or type(stationary) is not bool:
        raise ValueError('original Boolean countercontrol switches are required')
    states, targets, marks, score, _ = marked_source()
    weights = 6*marks-1
    aa = np.zeros((12, 12), dtype=np.int64)
    rr = np.zeros_like(aa); ar = np.zeros_like(aa); lag = np.zeros_like(aa)
    mean = np.zeros(12, dtype=np.int64)
    terminal = np.zeros(66, dtype=np.int64)
    origins = range(66) if stationary else (0,)
    branches = 0
    for origin in origins:
        first = targets[origin]
        if sticky:
            columns = np.array([[2*(i//2), 2*(i//2)+1] for i in range(60)])
            final = targets[first[:, None], columns]
            second_weight = weights[columns]
        else:
            final = targets[first]
            second_weight = np.broadcast_to(weights, final.shape)
        one = np.broadcast_to(score[first, None, :], (*final.shape, 12))
        two = score[final]
        r1 = weights[:, None, None]*one
        r2 = second_weight[:, :, None]*two
        attempted = (one+two).reshape(-1, 12)
        difference = (r1+r2).reshape(-1, 12)
        aa += attempted.T@attempted
        rr += difference.T@difference
        ar += attempted.T@difference
        lag += r2.reshape(-1, 12).T@r1.reshape(-1, 12)
        mean += difference.sum(axis=0)
        terminal += np.bincount(final.ravel(), minlength=66)
        branches += final.size
    denominator = 36*branches
    return {'branches': branches, 'attempt_second_moment': sp.ImmutableMatrix(aa)/denominator,
        'difference_second_moment': sp.ImmutableMatrix(rr)/denominator,
        'cross_second_moment': sp.ImmutableMatrix(ar)/denominator,
        'difference_lag_one': sp.ImmutableMatrix(lag)/denominator,
        'difference_mean': sp.ImmutableMatrix(mean)/(6*branches),
        'terminal_counts': tuple(map(int, terminal)), 'initial_state': states[0]}


def joint_clt():
    """Cross-covariance derives the incident covariance on the same history."""
    attempt = occupation.field('attempt')
    incident = occupation.field('incident')
    local = attempt['laplacian']-incident['remote_laplacian']
    r = attempt['laplacian']/60-local/10
    cross = 60*occupation.CHI*attempt['green']*r
    delta = 5*occupation.CHI*occupation.P0
    if attempt['covariance']+cross+cross.T+delta != 6*incident['covariance']:
        raise ArithmeticError('joint clock covariance does not recover the exact stopped law')
    return sp.ImmutableMatrix(cross), sp.ImmutableMatrix(delta)


def martingale_certificate():
    """Exact forward innovations for the JOINT reward functional CLT.

    This is independent of using reverse orthogonality to obtain the clock
    difference covariance. It checks each conditional mean on all 66 states.
    """
    _, targets, marks, score, _ = marked_source()
    attempt = occupation.field('attempt')
    local = attempt['laplacian']-occupation.field('incident')['remote_laplacian']
    r = attempt['laplacian']/60-local/10
    h = 60*attempt['green']
    qa, qr = h-occupation.P0, r*h
    denominator = int(sp.ilcm(*(x.q for matrix in (qa, qr) for x in matrix)))
    qa = np.array(qa*denominator, dtype=object)
    qr = np.array(qr*denominator, dtype=object)
    after = score[targets].astype(object)
    jump = after-score[:, None, :]
    a = denominator*after+jump@qa.T
    z = denominator*(6*marks-1)[None, :, None]*after+jump@qr.T
    if np.any(a.sum(axis=1)) or np.any(z.sum(axis=1)):
        raise ArithmeticError('a joint reward innovation is not conditionally centered')
    a, z = a.reshape(-1, 12), z.reshape(-1, 12)
    normalizer = 36*denominator**2*len(a)
    blocks = tuple(sp.ImmutableMatrix(x)/normalizer for x in (a.T@a, z.T@z, a.T@z))
    cross, delta = joint_clt()
    if blocks != (attempt['covariance'], delta, cross):
        raise ArithmeticError('native joint martingale covariance disagrees with the Green formulas')
    return {'native_innovations': len(a), 'conditional_states': len(score),
            'all_conditional_means_zero': True, 'joint_asymptotic_covariance_checked': True}


def cylinder_action(times=(F(1, 2), F(1), F(2))):
    """Gaussian increment action AFTER the native functional CLT."""
    if type(times) is not tuple or not times or len(times) > 8:
        raise ValueError('one through eight increasing exact observation times required')
    for time in times:
        _positive(time)
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError('observation times must strictly increase')
    size = len(times)
    time_cov = sp.Matrix(size, size, lambda i, j: min(times[i], times[j]))
    increments = tuple(b-a for a, b in zip((F(), *times[:-1]), times, strict=True))
    difference = sp.eye(size)
    for i in range(1, size):
        difference[i, i-1] = -1
    time_precision = difference.T*sp.diag(*(1/x for x in increments))*difference
    spatial = occupation.field('attempt')['covariance'][:11, :11]/(120*occupation.CHI)
    precision = spatial.inv()
    covariance = sp.kronecker_product(time_cov, spatial)
    inverse = sp.kronecker_product(time_precision, precision)
    if covariance*inverse != sp.eye(11*size):
        raise ArithmeticError('the native Gaussian path action is not its joint precision')
    return {'times': times, 'increments': increments,
        'time_covariance': time_cov, 'time_precision': time_precision,
        'spatial_covariance': spatial, 'spatial_precision': precision,
        'dimension': 11*size}


def _rows(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def run():
    _, _, _, _, branches = marked_source()
    actual = direct_two_steps()
    predicted = finite_covariances(2)
    for key, expected in zip(('attempt_second_moment', 'difference_second_moment',
                              'cross_second_moment'), predicted, strict=True):
        if actual[key] != expected:
            raise ArithmeticError('the complete native two-step oracle differs from the theorem')
    if actual['difference_lag_one'] != sp.zeros(12) or actual['difference_mean'] != sp.zeros(12, 1):
        raise ArithmeticError('the reverse martingale covariance check failed')
    sticky, biased = direct_two_steps(sticky=True), direct_two_steps(stationary=False)
    if (sticky['difference_lag_one'] == sp.zeros(12)
            or len(set(sticky['terminal_counts'])) != 1
            or biased['difference_mean'] == sp.zeros(12, 1)):
        raise ArithmeticError('a hypothesis countercontrol stopped distinguishing the laws')
    cross, delta = joint_clt()
    action = cylinder_action()
    return {'schema': 'oph.native_joint_observation_process.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'native_exact_controls': {'one_step_branches': branches,
            'two_step_complete_branches': actual['branches'],
            'post_repair_mark_mean_zero_given_output_checked': True,
            'same_history_clock_difference_lag_one_zero': True,
            'difference_covariance_per_attempt': _rows(delta),
            'attempt_difference_asymptotic_cross_covariance': _rows(cross),
            'joint_covariance_recovers_exact_incident_clock': True},
        'joint_functional_CLT_innovations': martingale_certificate(),
        'countercontrols': {
            'sticky_seam_preserves_uniform_one_time_state': True,
            'sticky_seam_two_step_branches': sticky['branches'],
            'sticky_seam_difference_lag_one_trace': str(sp.trace(sticky['difference_lag_one'])),
            'nonstationary_preparation_difference_mean': _rows(biased['difference_mean']),
            'finite_native_adjacent_occupation_covariance': str(occupation.CHI*(
                occupation.field('attempt')['kernel']*occupation.P0)[0, 0]),
            'uncentered_constant_readout_random_time_change_noise_variance_per_attempt': '5'},
        'Gaussian_cylinder_action': {'dimension': action['dimension'],
            'times': list(map(str, action['times'])),
            'time_covariance': _rows(action['time_covariance']),
            'time_precision': _rows(action['time_precision']),
            'joint_precision_equals_time_precision_tensor_native_spatial_precision': True,
            'applies_after_functional_CLT_not_at_finite_native_count': True},
        'all_level_path_bounds': {str(level): {k: str(v) for k, v in path_bounds(
            level, count=100, degree=20).items()} for level in (0, 1, 3, 10, 30)},
        'theorem': {'joint_limit_in_C_of_time_into_H_minus_2_of_S2': True,
            'time_covariance': 'min(t,u) times the native spatial subsequence covariance',
            'normalized_native_clocks_are_jointly_coupled_to_the_same_path': True,
            'finite_count_same_prefix_clock_error_bound': True,
            'incident_index_time_change_at_each_fixed_source_level': True,
            'limit_order': 'native observation count to infinity, then common spatial subsequence',
            'full_field_requires_initial_load_and_complete_current_read_contract': True,
            'instantaneous_density_process_identified_with_Brownian_path': False,
            'uniform_joint_count_refinement_CLT_error': False,
            'unique_spatial_limit_or_physical_radial_coordinate': False,
            'Lorentzian_propagator_interacting_QFT_or_modular_flow_identity': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', type=Path)
    modes.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Joint native clock coupling, continuum observation paths and Gaussian action verified.')
    elif args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()

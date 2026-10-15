"""Exact diagonal-action response, with the global-energy gauge retained.

The response is twice the initial probability derivative of a Ramsey probe.
It is an identification/design object, not finite-shot derivative data.  The
binomial basis exposes which occupation levels each candidate can affect.
"""
from __future__ import annotations

from fractions import Fraction
import math
from numbers import Integral, Real

import sympy as sp


def integer(value, name, minimum=0, maximum=64):
    if (isinstance(value, bool) or not isinstance(value, Integral)
            or not minimum <= value <= maximum):
        raise ValueError(f'{name} must be an integer in [{minimum}, {maximum}]')
    return int(value)


def finite(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number')
    return float(value)


def response_matrix(levels, order=5):
    """Exact gap responses for |0> to |n>, including unobservable k=0."""
    order = integer(order, 'polynomial order', 1, 16)
    levels = tuple(integer(n, 'occupation', 1, 64) for n in levels)
    if len(levels) != len(set(levels)):
        raise ValueError('distinct occupations required')
    return sp.Matrix([[math.comb(n, k)-math.comb(0, k)
                       for k in range(order+1)] for n in levels]) if levels else sp.zeros(0, order+1)


def rank_report(levels, order=5):
    levels = tuple(levels)
    matrix = response_matrix(levels, order)
    return {
        'occupations': list(levels),
        'columns': [f'binom(N,{k})' for k in range(order+1)],
        'response_exact': [[str(x) for x in row] for row in matrix.tolist()],
        'rank': matrix.rank(),
        'nullspace_exact': [[str(x) for x in vector] for vector in matrix.nullspace()],
        'global_energy_gauge_column': 0,
    }


def complete_design(initial=(1, 2, 3), candidates=range(4, 8), order=5):
    """Greedily select the lowest available occupation that increases rank.

    This minimizes maximum occupation for this nested binomial family.  It
    does not claim optimal finite-shot conditioning among all instruments.
    """
    initial, candidates = tuple(initial), tuple(candidates)
    current = list(initial)
    matrix = response_matrix(current, order)
    accepted = []
    for level in sorted(integer(n, 'candidate occupation', 1, 64) for n in candidates):
        if level in current:
            continue
        trial = response_matrix((*current, level), order)
        if trial.rank() > matrix.rank():
            current.append(level)
            accepted.append(level)
            matrix = trial
        if matrix.rank() == order:
            break
    return {'before': rank_report(initial, order), 'added_occupations': accepted,
            'after': rank_report(tuple(current), order),
            'remaining_non_gauge_directions': order-matrix.rank()}


def finite_difference_inverse(order=5):
    """Exact inverse from measured energy gaps to nonconstant coefficients."""
    matrix = response_matrix(range(1, order+1), order)[:, 1:]
    return [[str(x) for x in row] for row in matrix.inv().tolist()]


def bernoulli_design(time, energy_difference, trials=19, probability_error=0.):
    """Midpoint Ramsey discrimination and a conditional Hoeffding guarantee.

    Each trial is a fresh independently prepared Bernoulli experiment.  The
    error argument is an *assumed absolute probability bound*, not inferred
    from floating source residuals.  Each candidate's probability can move
    adversarially by that amount while the decision threshold stays 1/2.
    """
    time, difference = finite(time, 'time'), finite(energy_difference, 'energy difference')
    error = finite(probability_error, 'probability error')
    trials = integer(trials, 'trial count', 1, 1_000_000)
    if time < 0 or not 0 <= error <= 1:
        raise ValueError('nonnegative time and probability error in [0,1] required')
    gap = abs(math.sin(time*difference/2))
    margin = max(0., gap/2-error)
    return {
        'time_parameter': 'supplied model coordinate', 'model_time': time,
        'energy_difference': difference, 'ideal_binary_tv': gap,
        'ideal_probabilities': [(1-gap)/2, (1+gap)/2],
        'fresh_independent_trials': trials,
        'assumed_probability_error': error,
        'worst_case_one_hypothesis_error_upper': math.exp(-2*trials*margin**2),
        'decision_rule': 'choose the higher-probability candidate when the sample mean exceeds 1/2',
        'bound_scope': 'ideal or externally bounded Bernoulli probabilities; no numerical roundoff certification',
    }


def derivative_noise_radius(time, trials, failure_probability):
    """Sampling-only radius for 2*(sample_mean-1/2)/time.

    This does not bound finite-time derivative bias.  The initial probability
    1/2 is assumed exactly known; estimating it separately adds uncertainty.
    """
    time = finite(time, 'time')
    trials = integer(trials, 'trial count', 1, 10**18)
    failure_probability = finite(failure_probability, 'failure probability')
    if time <= 0 or not 0 < failure_probability < 1:
        raise ValueError('positive time and failure probability in (0,1) required')
    return 2*math.sqrt(math.log(2/failure_probability)/(2*trials))/time


def exact_blind_compression(cutoff=4, order=5):
    cutoff = integer(cutoff, 'cutoff', 1, 32)
    order = integer(order, 'polynomial order', 1, 32)
    return {f'binom(N,{k})': [str(Fraction(math.comb(n, k))) for n in range(cutoff)]
            for k in range(cutoff, order+1)}

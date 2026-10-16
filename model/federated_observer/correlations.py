"""Exact detector correlations without enumerating the many-body state space.

The parameter is a number of uniform global seam attempts.  The observable
is the net signed incident current, accessible without an occupation oracle.
No physical duration or relativistic identification is assigned to attempts.
"""
from __future__ import annotations

from fractions import Fraction as F

from . import source


SOURCE_PATHS = (
    'model/federated_observer/correlations.py',
    'Lean/Thermodynamics/NativeRepairFourCycle.lean',
    'Lean/Screen/RegionalContinuity.lean',
)


def _integer(value, lower, upper, name):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f'{name} must be an integer between {lower} and {upper}')
    return value


def _graph(neighbours):
    if (type(neighbours) is not tuple or len(neighbours) < 2
            or any(type(row) is not tuple for row in neighbours)):
        raise ValueError('an immutable simple undirected neighbour table is required')
    vertices = len(neighbours)
    for i, row in enumerate(neighbours):
        if (any(type(j) is not int or not 0 <= j < vertices or j == i for j in row)
                or len(set(row)) != len(row)):
            raise ValueError('invalid simple graph row')
    if any(i not in neighbours[j] for i, row in enumerate(neighbours) for j in row):
        raise ValueError('undirected seams are required')
    denominator = sum(map(len, neighbours))
    if denominator == 0:
        raise ValueError('at least one seam is required')
    return vertices, denominator


def one_point_counts(neighbours, start, steps):
    """Integer row of (2M I-L)^steps and denominator (2M)^steps.

    A local native fair-coin instrument has conditional mean
    E[n_a'|n]=(n_a+n_b)/2.  Summing over its M possible seams proves
    E[n'|n]=(I-L/(2M))n for every binary occupation state.  This closes
    a marginal, not the complete many-body law or quantum channel.
    """
    vertices, denominator = _graph(neighbours)
    _integer(start, 0, vertices-1, 'start')
    _integer(steps, 0, 64, 'attempt count')
    counts = [int(i == start) for i in range(vertices)]
    for _ in range(steps):
        following = [count*(denominator-len(neighbours[i]))
                     for i, count in enumerate(counts)]
        for i, count in enumerate(counts):
            if count:
                for j in neighbours[i]:
                    following[j] += count
        counts = following
    total = denominator**steps
    if min(counts) < 0 or sum(counts) != total:
        raise ArithmeticError('one-point propagation lost exact probability mass')
    return tuple(counts), total


def covariance(neighbours, particles, first, second, gap):
    """Cov(n_first(0),n_second(gap)) in uniform fixed-number preparation."""
    vertices, _ = _graph(neighbours)
    _integer(particles, 0, vertices, 'particle count')
    _integer(second, 0, vertices-1, 'second detector')
    row, total = one_point_counts(neighbours, first, gap)
    density = F(particles, vertices)
    return density*(1-density)/(vertices-1)*(vertices*F(row[second], total)-1)


def increment_covariance(neighbours, particles, first, second, interval_a, interval_b):
    """Cov of two net incident-current reads over declared attempt intervals.

    n_i(v)-n_i(u) is the sum of signed recorded currents incident on i.
    Stationarity and reversibility give the four-covariance identity below.
    An absolute initial occupation is unnecessary for this difference read.
    """
    for interval in (interval_a, interval_b):
        if (type(interval) is not tuple or len(interval) != 2
                or any(type(t) is not int or not 0 <= t <= 64 for t in interval)
                or interval[0] > interval[1]):
            raise ValueError('ordered integer attempt intervals in [0,64] are required')
    u, v = interval_a
    r, s = interval_b
    def c(gap):
        return covariance(neighbours, particles, first, second, abs(gap))
    return c(v-s)-c(v-r)-c(u-s)+c(u-r)


def run():
    rows = []
    for level in (0, 1, 2):
        graph = source.federation(level)
        a, b = graph.inter_seams[0]
        # Both preparations are disclosed; neither is inferred from an experiment.
        for particles in (2, graph.ports//2):
            for attempts in (1, 2, 4, 8):
                variance = increment_covariance(graph.neighbours, particles, a, a,
                                                (0, attempts), (0, attempts))
                shared = increment_covariance(graph.neighbours, particles, a, b,
                                              (0, attempts), (0, attempts))
                successive = increment_covariance(graph.neighbours, particles, a, a,
                                                  (0, attempts), (attempts, 2*attempts))
                if variance <= 0 or shared >= 0 or successive >= 0:
                    raise ArithmeticError('the declared current correlation controls changed')
                rows.append({'level': level, 'elementary_observers': graph.cells,
                    'ports': graph.ports, 'particles': particles,
                    'cross_seam_detectors': [a, b], 'uniform_seam_attempts': attempts,
                    'net_current_variance': str(variance),
                    'opposite_endpoint_net_current_covariance': str(shared),
                    'successive_window_net_current_covariance': str(successive),
                    'cross_endpoint_covariance_over_local_variance': str(shared/variance)})
                if attempts == 1 and shared/variance != F(-1, 6):
                    raise ArithmeticError('the degree-six cross-seam one-attempt ratio changed')
    return {
        'preparation': 'uniform over all binary configurations with the disclosed global particle count',
        'parameter': 'number of uniformly selected global seam attempts with a fair endpoint coin',
        'observable': 'signed incident-current sums over the stated attempt windows',
        'initial_occupation_access_needed_for_increment_reads': False,
        'one_point_kernel': 'P=I-L/(2M), M=number of seams',
        'static_covariance': 'rho*(1-rho)*V/(V-1)*(I-11^T/V)',
        'lag_covariance': 'rho*(1-rho)/(V-1)*(V*P^n_ij-1)',
        'same_window_increment_covariance': '2*(C_ij(0)-C_ij(n))',
        'successive_window_increment_covariance': 'C_ij(a)+C_ij(b)-C_ij(a+b)-C_ij(0)',
        'integer_arithmetic_predictions': rows,
        'exponential_many_body_matrix_constructed': False,
        'closes_all_higher_correlations': False,
        'matches_new_empirical_measurements': False,
        'physical_clock_or_lorentz_metric_inferred': False,
    }

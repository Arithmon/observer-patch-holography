"""Exact ideal-operation bounds; free spatial estimates are not interacting ones."""
from fractions import Fraction as F
from math import factorial

from .exact import EDGES, cutoffs, sqrt_upper


def _rational(value, name, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError(f'{name} must be an exact nonnegative rational, positive where required')
    return F(value)


def local_centered_matrix(d, coupling=1):
    """Bargmann entries of the true compressed onsite polynomial, Omega=2.

    The removed scalar 1+lambda/128 occurs equally in the exact evolution
    and every symmetric product, so its phase cancels in their comparison.
    """
    if type(d) is not int or d < 1 or type(coupling) is not int or coupling not in (0, 1):
        raise ValueError('positive integer dimension and coupling zero or one required')
    out = [[F(0) for _ in range(d)] for _ in range(d)]
    for n in range(d):
        out[n][n] = 2*n+F(coupling*6*n*(n+1), 384)
        for delta, factor in ((4, 1), (2, 4*n+6),
                              (-2, n*(n-1)*(4*n-2)),
                              (-4, n*(n-1)*(n-2)*(n-3))):
            if 0 <= n+delta < d:
                out[n+delta][n] = F(coupling*factor, 384)
    return out


def onsite_norm_upper(d, coupling=1):
    matrix = local_centered_matrix(d, coupling)
    # Change to the normalized Fock basis before applying the Hermitian
    # absolute-row-sum bound. Every square root is rounded outwards.
    rows = [sum((abs(matrix[i][j])*sqrt_upper(F(factorial(i), factorial(j)))
                 for j in range(d)), F(0)) for i in range(d)]
    return max(rows)


def generator_sum_bound(dimensions, coupling=1):
    dims = cutoffs(dimensions)
    if dims is None:
        raise ValueError('finite local dimensions required')
    # ||(a+a*)/2|| <= sqrt(d-1); the twelve bond terms commute, but
    # the conservative triangle bound does not need that fact.
    q = [sqrt_upper(F(d-1)) for d in dims]
    return (sum((onsite_norm_upper(d, coupling) for d in dims), F(0))
            + sum((q[i]*q[j] for i, j in EDGES), F(0)))


def product_error(dimensions, time, steps, coupling=1):
    """Vector/operator error of n symmetric onsite/bond steps, <=S^3 t^3/(3n^2)."""
    time = _rational(time, 'model time')
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer symmetric steps required')
    s = generator_sum_bound(dimensions, coupling)
    return min(F(2), s**3*time**3/(3*steps**2))


def free_window(spacing, momentum_radius, model_time, mass=1, speed=1, particles=1):
    """Common band-limited free momentum representation, away from aliasing.

    K*a <=1 is a convenient rational sufficient condition for the common
    Brillouin-zone support. This compares the infinite spatial free lattice
    with Klein--Gordon, not the eight-site experiment or interacting flow.
    """
    a = _rational(spacing, 'spacing', True)
    k = _rational(momentum_radius, 'momentum radius')
    t = _rational(model_time, 'model time')
    m, c = _rational(mass, 'mass', True), _rational(speed, 'speed', True)
    if a*k > 1 or type(particles) is not int or particles < 0:
        raise ValueError('require spacing*momentum_radius <=1 and nonnegative integer particle cap')
    frequency = c*c*a*a*k**4/(24*m)
    return {'spacing': a, 'momentum_radius': k, 'model_time': t,
            'maximum_particles': particles,
            'frequency_error_upper': frequency,
            'free_evolution_vector_error_upper': min(F(2), particles*t*frequency),
            'free_equal_time_covariance_multiplier_error_upper': frequency/(2*m*m),
            'interaction_included': False, 'eight_site_boundary_included': False}


def gap_statement():
    return {
        'target': 'Untruncated canonical oscillators in every finite spatial box',
        'lower_bound_exact': '1',
        'premises': 'Supplied m=1; nonnegative quartic coupling and edge weights; unit hbar',
        'method': 'Strong convexity and harmonic-oscillator comparison, followed by configuration-space exhaustion',
        'primary_reference': 'https://arxiv.org/abs/1006.1686',
        'theorem': 'Andrews--Clutterbuck, Theorem 1.3 and Section 5',
        'finite_occupation_gap_claimed': False,
        'infinite_spatial_volume_representation_constructed': False,
        'gap_measured_by_protocol': False,
    }


def serializable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    return value

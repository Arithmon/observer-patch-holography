"""Exact Fourier coefficients, compressed oscillator polynomials and bounds.

The new modes are cos(3x), sin(3x), not a complete momentum shell. All
occupation powers are formed before compression. Time is model time.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import factorial

from ..scalar_observations.bounds import sqrt_upper


LOW_DIMENSIONS = (2, 4, 4, 2, 2, 2, 2)
DIMENSIONS = LOW_DIMENSIONS + (4, 4)
MODES = ('zero', 'cos_x', 'sin_x', 'cos_y', 'sin_y', 'cos_z', 'sin_z',
         'cos_3x', 'sin_3x')
WAVES = ((0, 0, 0), (1, 0, 0), (1, 0, 0), (0, 1, 0), (0, 1, 0),
         (0, 0, 1), (0, 0, 1), (3, 0, 0), (3, 0, 0))


def integer(value, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError('an original integer in the admitted range is required')
    return value


def rational(value, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('an exact nonnegative rational is required')
    return F(value)


@lru_cache(maxsize=2)
def quartic_coefficients(extended=True):
    """Multinomial coefficients of V integral phi**4, in exact rationals.

    Nonconstant real modes have Fourier coefficients 1/sqrt(2) and
    +/-i/sqrt(2). Constant Fourier coefficients are selected explicitly.
    This is an independent integral calculation, not the model's formula.
    """
    if type(extended) is not bool:
        raise ValueError('boolean band selector required')
    count = 9 if extended else 7
    result = {}
    for indices in combinations_with_replacement(range(count), 4):
        factors = []
        nonconstant = 0
        for index in indices:
            if index == 0:
                factors.append([((0, 0, 0), 1)])
            else:
                nonconstant += 1
                wave = WAVES[index]
                sine = index % 2 == 0
                factors.append([(tuple(sign*x for x in wave), (-sign*1j if sine else 1))
                                for sign in (-1, 1)])
        total = 0j
        for choices in product(*factors):
            if all(sum(row[0][axis] for row in choices) == 0 for axis in range(3)):
                coefficient = 1
                for _, number in choices:
                    coefficient *= number
                total += coefficient
        if total == 0:
            continue
        if total.imag or nonconstant % 2:
            raise ArithmeticError('unexpected nonrational Fourier coefficient')
        powers = tuple(indices.count(i) for i in range(count))
        multiplicity = factorial(4)
        for power in powers:
            multiplicity //= factorial(power)
        result[powers] = F(int(total.real)*multiplicity, 2**(nonconstant//2))
    return result


def declared_quartic():
    """The explicit nine-mode action, independently checked by Fourier sums."""
    result = {}

    def term(coefficient, powers):
        key = tuple(powers.get(i, 0) for i in range(9))
        result[key] = result.get(key, F(0)) + F(coefficient)

    for i in range(9):
        term(1 if i == 0 else F(3, 2), {i: 4})
    for i in range(9):
        for j in range(i+1, 9):
            same_pair = (i, j) in ((1, 2), (3, 4), (5, 6), (7, 8))
            term(3 if same_pair else 6, {i: 2, j: 2})
    term(2, {1: 3, 7: 1})
    term(-6, {1: 1, 2: 2, 7: 1})
    term(6, {1: 2, 2: 1, 8: 1})
    term(-2, {2: 3, 8: 1})
    return result


@lru_cache(maxsize=32)
def polynomial(dimension, power, momentum=False):
    """Full Q**power (Q=(a+a*)/2) or P**2 in Bargmann coordinates."""
    d, power = integer(dimension), integer(power, 0)
    if type(momentum) is not bool or (momentum and power != 2):
        raise ValueError('only the squared momentum polynomial is admitted')
    matrix = [[F(0) for _ in range(d)] for _ in range(d)]
    for column in range(d):
        entries = {column: F(1)}
        for _ in range(power):
            following = {}
            for n, value in entries.items():
                following[n+1] = following.get(n+1, F(0)) + value*(1 if momentum else F(1, 2))
                if n:
                    following[n-1] = following.get(n-1, F(0)) + value*n*(-1 if momentum else F(1, 2))
            entries = following
        for row, value in entries.items():
            if row < d:
                matrix[row][column] = -value if momentum else value
    return tuple(tuple(row) for row in matrix)


def row_norm(matrix):
    return max(sum((abs(value)*sqrt_upper(F(factorial(i), factorial(j)))
                    for j, value in enumerate(row)), F(0))
               for i, row in enumerate(matrix))


def local_norm(d, power):
    return row_norm(polynomial(d, power))


def residual_bound(coupling=1):
    """R**2 upper bound for Q H J on the declared 9-mode occupation code.

    Only high labels (2,0),(0,2),(2,2),(1,0),(0,1) occur. Their ranges
    are orthogonal. The first two blocks include free high-blank mismatch;
    the odd blocks come exclusively from the resonant quartic terms.
    """
    g = rational(coupling)
    squares = sum((local_norm(d, 2) for d in LOW_DIMENSIONS), F(0))
    even_without_root_two = F(3, 4)+g/32+g*squares/16
    x, y = LOW_DIMENSIONS[1:3]
    odd_z = g*(local_norm(x, 3)+3*local_norm(x, 1)*local_norm(y, 2))/24
    odd_w = g*(3*local_norm(x, 2)*local_norm(y, 1)+local_norm(y, 3))/24
    squared = 4*even_without_root_two**2+(g/64)**2+odd_z**2+odd_w**2
    return {'squared_residual_norm_upper': squared,
            'residual_norm_upper': sqrt_upper(squared),
            'even_high_blocks_without_sqrt_two_bound': even_without_root_two,
            'odd_high_block_bounds': (odd_z, odd_w),
            'two_high_excitation_block_exact': g/64,
            'matched_mass_squared_increment': g/4,
            'matched_scalar_increment': F(7, 2)+g/32}


def history_bound(total_time, coupling=1):
    time = rational(total_time)
    result = residual_bound(coupling)
    squared = result['squared_residual_norm_upper']*time*time
    return {**result, 'total_model_time': time,
            'reduced_complete_output_trace_distance_upper': min(F(1), squared),
            'high_retained_complete_output_trace_distance_upper': min(F(1), sqrt_upper(squared)),
            'high_blank_survival_failure_probability_upper': min(F(1), squared),
            'postselection_used': False,
            'all_low_input_states_and_untouched_references': True,
            'arbitrary_identical_low_only_adaptive_instruments': True,
            'unrestricted_high_occupation_or_uv_limit_included': False}


def grid_history_bound(grid, total_time):
    """Same nine-mode code: finite differences versus exact gradients.

    Quartic Fourier quadrature is exact for grid>=13. For harmonic k,
    0<=k**2-kappa_grid(k)<=k**4*pi**2/(3*grid**2), using pi<22/7.
    """
    grid = integer(grid, 13)
    time = rational(total_time)
    weighted = sum((sum(k*k for k in wave)**2*local_norm(d, 2)
                    for wave, d in zip(WAVES, DIMENSIONS, strict=True)), F(0))
    generator = F(484, 294*grid*grid)*weighted
    return {'grid_per_axis': grid,
            'same_nine_mode_hamiltonian_norm_difference_upper': generator,
            'complete_history_trace_distance_upper': min(F(1), time*generator),
            'total_model_time': time,
            'quartic_quadrature_exact': True,
            'number_of_independent_modes': 9,
            'grid_increases_bandwidth': False}


def serializable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    return value

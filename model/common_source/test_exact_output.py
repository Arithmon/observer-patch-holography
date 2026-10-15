"""Independent algebraic controls for exact reported-output comparisons."""
from fractions import Fraction as F
import math

import mpmath
import numpy as np
import pytest

from . import exact_output as output, protocol, smooth_history as smooth


@pytest.mark.parametrize('value', [F(0), F(1, 9), F(2), F(1, 10**700),
                                   F(10**700), F(2**192+1, 2**192)])
def test_rational_square_root_brackets_original_values_without_float_conversion(value):
    lower, upper = output._sqrt_interval(value)
    assert type(lower) is type(upper) is F
    assert 0 <= lower <= upper
    assert lower*lower <= value <= upper*upper
    if value:
        assert lower > 0
        assert (upper-lower)/lower < F(1, 2**180)


def test_rank_one_distance_matches_closed_original_component_formula():
    program = protocol.Program()
    enclosure = smooth.enclose(program, dimension=2)
    vector = np.zeros(program.work_dimension, complex)
    vector[0], vector[1] = .25+.75j, .5-.25j
    answer = output.compare_vectors(program, {(): vector}, enclosure)
    # ||v||^2=15/16, ||p||^2=1, |<v,p>|^2=5/8, p=|0>.
    exact_squared = F(321, 1024)
    assert answer['polynomial_distance_lower']**2 <= exact_squared
    assert answer['polynomial_distance_upper']**2 >= exact_squared
    assert answer['source_mass_real'] == F(15, 16)
    assert answer['source_mass_imag'] == 0
    assert answer['distance_upper']-answer['distance_lower'] < F(1, 10**50)


def _mp(value, ctx):
    value = F(value)
    return ctx.mpf(value.numerator)/value.denominator


def _independent_embedding(program, enclosure, ctx):
    """Construct a finite matrix containing every nonzero polynomial entry."""
    d = program.dimension
    keys = [(n0, n1, g) for n1 in range(d) for n0 in range(d) for g in range(4)]
    polynomial = enclosure.branches[()]
    keys += sorted(set(polynomial.coefficients)-set(keys))
    vector = []
    for n0, n1, g in keys:
        real, imag = polynomial.coefficients.get((n0, n1, g), (0, 0))
        vector.append(ctx.mpc(real, imag)*ctx.sqrt(math.factorial(n0)*math.factorial(n1))
                      /polynomial.denominator/ctx.sqrt(_mp(enclosure.normalization, ctx)))
    return vector


@pytest.mark.parametrize('dimension', [2, 4])
def test_rank_one_formula_includes_factorial_weights_and_omitted_rectangle(dimension):
    program = protocol.Program(dimension=dimension)
    resolution = smooth.Resolution(cutoff=4, preparation=2, degree=4,
                                   steps_per_unit=2, bits=64)
    enclosure = smooth.enclose(program, resolution=resolution)
    source = np.zeros(program.work_dimension, complex)
    source[0], source[4], source[-1] = .75+.125j, -.25j, .5
    answer = output.compare_vectors(program, {(): source}, enclosure)
    ctx = mpmath.mp.clone(); ctx.dps = 120
    target = _independent_embedding(program, enclosure, ctx)
    original = [ctx.mpc(_mp(z.real, ctx), _mp(z.imag, ctx)) for z in source]
    mass = sum(abs(z)**2 for z in original)
    norm = sum(abs(z)**2 for z in target)
    inner = sum(ctx.conj(z)*target[i] for i, z in enumerate(original))
    expected = ctx.sqrt((mass+norm)**2-4*abs(inner)**2)/2
    assert _mp(answer['polynomial_distance_lower'], ctx) <= expected
    assert expected <= _mp(answer['polynomial_distance_upper'], ctx)
    assert answer['distance_upper'] == answer['polynomial_distance_upper']+enclosure.output_error


@pytest.mark.parametrize('dimension,infinite', [(4, False), (2, True)])
def test_matrix_frobenius_bound_matches_direct_embedding_without_hermitian_assumption(dimension, infinite):
    program = protocol.Program(dimension=dimension)
    resolution = smooth.Resolution(cutoff=4, preparation=2, degree=4,
                                   steps_per_unit=2, bits=64)
    enclosure = smooth.enclose(program, dimension=None if infinite else dimension,
                               resolution=resolution)
    rho = np.zeros((program.work_dimension,)*2, complex)
    for index in (0, 4, 8, 12):
        rho[index, index] = .25
    rho[0, 8], rho[8, 0], rho[4, 0] = .125+.25j, -.125j, .5
    assert not np.array_equal(rho, rho.conj().T)
    answer = output.compare_densities(program, {(): rho}, enclosure)
    ctx = mpmath.mp.clone(); ctx.dps = 120
    target = _independent_embedding(program, enclosure, ctx)
    squared = ctx.mpf(0)
    for i, a in enumerate(target):
        for j, b in enumerate(target):
            z = rho[i, j] if i < len(rho) and j < len(rho) else 0j
            value = ctx.mpc(_mp(z.real, ctx), _mp(z.imag, ctx))-a*ctx.conj(b)
            squared += abs(value)**2
    frobenius = ctx.sqrt(squared)/2
    rank = program.work_dimension+int(infinite)
    assert answer['per_record'][()]['rank_upper'] == rank
    assert _mp(answer['polynomial_distance_lower'], ctx) <= frobenius
    assert frobenius*ctx.sqrt(rank) <= _mp(answer['polynomial_distance_upper'], ctx)
    assert _mp(answer['polynomial_distance_upper'], ctx)-frobenius*ctx.sqrt(rank) < ctx.mpf('1e-50')
    assert answer['source_mass_real'] == 1
    assert answer['source_mass_imag'] == 0


def test_complex_trace_is_reported_instead_of_repairing_a_nonhermitian_matrix():
    program = protocol.Program()
    enclosure = smooth.enclose(program, dimension=2)
    rho = np.zeros((program.work_dimension,)*2, complex)
    rho[0, 0] = 1+.25j
    result = output.compare_densities(program, {(): rho}, enclosure)
    assert result['source_mass_real'] == 1
    assert result['source_mass_imag'] == F(1, 4)
    assert result['distance_lower'] == F(1, 8)
    assert result['distance_upper'] == F(1, 2)


def test_complete_target_error_is_added_once_without_probability_normalization():
    program = protocol.Program(operations=(protocol.Operation('evolve', F(1, 64)),))
    resolution = smooth.Resolution(cutoff=2, preparation=0, degree=2,
                                   steps_per_unit=2, bits=32)
    enclosure = smooth.enclose(program, dimension=2, resolution=resolution)
    assert enclosure.output_error > 0
    vector = np.zeros(program.work_dimension, complex)
    vector[0] = 2
    answer = output.compare_vectors(program, {(): vector}, enclosure)
    assert answer['source_mass_real'] == 4
    assert answer['distance_lower'] == max(0, answer['polynomial_distance_lower']-enclosure.output_error)
    assert answer['distance_upper'] == answer['polynomial_distance_upper']+enclosure.output_error
    assert answer['distance_lower'] > 1


@pytest.mark.parametrize('value', [True, False, 1, 0, 1., 0.])
@pytest.mark.parametrize('warm', [False, True])
def test_square_root_rejects_nonfractions_before_any_cache_hit(value, warm):
    # Clear the private computation cache, including the pre-fix decorated
    # public helper, so each cold/warm pair tests the same original input.
    cached = getattr(output, '_cached_sqrt_interval', output._sqrt_interval)
    cached.cache_clear()
    if warm:
        assert output._sqrt_interval(F(int(value))) == (F(int(value)), F(int(value)))
    with pytest.raises(ValueError, match='exact rational'):
        output._sqrt_interval(value)

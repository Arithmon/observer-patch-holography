"""Independent local proof controls for complete smooth-history certificates.

Dense full-history and public-boundary controls live in the companion suite.
These tests check the generator domain, factorial metric, omitted support,
pointer phase/normalization and the exact complete-output comparison.
"""
from fractions import Fraction as F
import math

import mpmath as mp
import pytest

from . import infinite_history as ih, smooth_history as sh


def _rational(vector):
    return {key: (F(a, vector.denominator), F(b, vector.denominator))
            for key, (a, b) in vector.coefficients.items()}


def _inner(left, right):
    real = imaginary = F(0)
    for (a, b, gauge), (lr, li) in left.items():
        rr, ri = right.get((a, b, gauge), (F(0), F(0)))
        weight = math.factorial(a)*math.factorial(b)
        real += weight*(lr*rr+li*ri)
        imaginary += weight*(lr*ri-li*rr)
    return real, imaginary


def _distance_squared(left, right):
    keys = left.keys() | right.keys()
    difference = {key: tuple(x-y for x, y in zip(
        left.get(key, (F(0), F(0))), right.get(key, (F(0), F(0))), strict=True))
        for key in keys}
    return _inner(difference, difference)[0]


@pytest.mark.parametrize('dimension', [None, 2, 4])
def test_coordinate_is_self_adjoint_in_the_factorial_metric(dimension):
    high = 3 if dimension is None else dimension-1
    left = {(0, 0, 0): (2, 3), (high, 1, 2): (-3, 1), (1, 1, 2): (5, -2)}
    right = {(1, 0, 0): (-1, 4), (high-1, 1, 2): (3, 2)}
    assert _inner(left, sh._coordinate(right, dimension)) == _inner(
        sh._coordinate(left, dimension), right)


@pytest.mark.parametrize('dimension', [2, 4])
def test_finite_generator_is_compressed_between_powers(dimension):
    top = {(dimension-1, 0, 0): (1, 0)}
    infinite_square = sh._coordinate(sh._coordinate(top))
    finite_square = sh._coordinate(sh._coordinate(top, dimension), dimension)
    # Pi B^2 Pi - (Pi B Pi)^2 = d |d-1><d-1| at the upper boundary.
    assert infinite_square[dimension-1, 0, 0][0]-finite_square[dimension-1, 0, 0][0] == dimension


def test_infinite_hamiltonian_powers_keep_excursions_outside_the_finite_rectangle():
    vacuum = {(0, 0, 0): (1, 0)}
    first = ih._hamiltonian_numerator(vacuum)
    assert first[2, 0, 0] == first[0, 2, 0] == (9, 0)
    infinite_square = ih._hamiltonian_numerator(first)
    finite_square = sh._compress(ih._hamiltonian_numerator(sh._compress(first, 2)), 2)
    # Two distinct excursions contribute 9 * (2*9) each to (128H)^2.
    assert infinite_square[0, 0, 0][0]-finite_square[0, 0, 0][0] == 324
    assert sh._compress(first, 2)[0, 0, 0] == (146, 0)


def test_taylor_remainder_is_the_full_next_vector_with_original_normalization():
    vector = ih._Polynomial({(0, 0, 0): (3, 4)}, 7)
    polynomial, remainder = sh._series(vector, F(2, 5), 1, sh._coordinate)
    assert _rational(polynomial) == {
        (0, 0, 0): (F(3, 7), F(4, 7)),
        (1, 0, 0): (F(8, 35), F(-6, 35))}
    assert _rational(remainder) == {
        (0, 0, 0): (F(-6, 175), F(-8, 175)),
        (2, 0, 0): (F(-6, 175), F(-8, 175))}
    assert sh._norm(remainder, F(25, 49)) == F(12, 625)


def test_record_mixing_obeys_the_complete_parallelogram_identity():
    left = ih._Polynomial({(0, 0, 0): (2, -3), (2, 1, 1): (1, 5)}, 7)
    right = ih._Polynomial({(0, 0, 0): (-4, 1), (1, 2, 3): (2, -1)}, 11)
    outputs = [sh._combine(left, right, sign) for sign in (-1, 1)]
    assert sum((sh._norm(v, F(1)) for v in outputs), F(0)) == (
        sh._norm(left, F(1))+sh._norm(right, F(1)))/2
    # Both arms affect both outputs; a normalized or selected branch breaks this.
    assert all(sh._norm(v, F(1)) != 1 for v in outputs)


def _mp_fraction(value):
    return mp.mpf(value.numerator)/value.denominator


def _vacuum_pointer_distance(outputs, *, dimension=None):
    """Full omitted-tail distance from an independent normalized analytic state."""
    alpha = 1/mp.sqrt(8)
    overlap = mp.mpc(0)
    norm = mp.mpf(0)
    for bit, vector in enumerate(outputs):
        sign = -1 if bit == 0 else 1
        for (n0, n1, gauge), (real, imag) in _rational(vector).items():
            actual = mp.mpc(_mp_fraction(real), _mp_fraction(imag))
            norm += math.factorial(n0)*math.factorial(n1)*abs(actual)**2
            if n1 or gauge:
                continue
            if dimension is None:
                common = mp.exp(-alpha**2/2)/math.factorial(n0)
                expected = common*((-1j*alpha)**n0+sign*1j*(1j*alpha)**n0)/2
            else:
                assert dimension == 2 and n0 in (0, 1)
                minus = mp.cos(alpha) if n0 == 0 else -1j*mp.sin(alpha)
                plus = mp.cos(alpha) if n0 == 0 else 1j*mp.sin(alpha)
                expected = (minus+sign*1j*plus)/2
            overlap += math.factorial(n0)*mp.conj(expected)*actual
    # Exact complete analytic output norm is one. This includes the entire
    # infinite tail without pretending a large finite truncation is exact.
    squared = 1+norm-2*mp.re(overlap)
    assert squared >= 0
    return mp.sqrt(squared)


@pytest.mark.parametrize('dimension', [None, 2])
def test_pointer_vacuum_error_covers_exact_displacement_or_finite_rotation(dimension):
    resolution = sh.Resolution(cutoff=4, preparation=0, degree=10,
                               steps_per_unit=2, bits=64)
    outputs, errors, maximum = sh._pointer(
        ih._Polynomial({(0, 0, 0): (1, 0)}), F(1), resolution, dimension)
    with mp.workdps(100):
        distance = _vacuum_pointer_distance(outputs, dimension=dimension)
        assert distance <= _mp_fraction(sum(errors.values(), F(0)))
        assert distance > 0  # The polynomial approximant is not the exact channel.
    assert errors['Taylor'] > 0 and errors['coordinate_coefficient'] > 0
    assert errors['coefficient_rounding'] > 0
    assert (errors['discarded_support'] > 0) == (dimension is None)
    assert maximum == (11 if dimension is None else 1)


def test_discarding_and_rounding_are_charged_without_renormalizing():
    vector = ih._Polynomial({(0, 0, 0): (2, 1), (1, 0, 0): (-1, 2),
                            (3, 0, 1): (1, -1)}, 3)
    normalization = F(7, 5)
    rounded, discarded, rounding = ih._cut_and_round(vector, 1, 32, normalization)
    kept = {k: v for k, v in _rational(vector).items() if sum(k[:2]) <= 1}
    removed = {k: v for k, v in _rational(vector).items() if sum(k[:2]) > 1}
    assert discarded**2 >= _inner(removed, removed)[0]/normalization > 0
    assert rounding**2 >= _distance_squared(kept, _rational(rounded))/normalization > 0
    assert sh._norm(rounded, normalization) != 1


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_native_kicks_preserve_complex_unnormalized_spectator_sectors(current):
    vector = ih._Polynomial({(0, 0, 0): (2, 3), (0, 0, 1): (-4, 1),
                            (2, 1, 2): (3, -2)}, 7)
    kicked = sh._kick(vector, current)
    assert sh._norm(kicked, F(1)) == sh._norm(vector, F(1))
    restored = sh._kick(kicked, -current)
    assert _rational(restored) == _rational(vector)


def test_rank_one_overlap_uses_both_normalizations_and_is_branch_phase_invariant():
    left = ih._Polynomial({(0, 0, 0): (1, 0), (1, 0, 0): (0, 1)}, 1)
    # Same normalized ket up to a global record phase, using another common
    # polynomial scale. Committed record phases cannot change this distance.
    phased = ih._Polynomial({(0, 0, 0): (0, 3), (1, 0, 0): (-3, 0)}, 1)
    opposite = ih._Polynomial({(0, 0, 0): (3, 0), (1, 0, 0): (0, -3)}, 1)
    normalizations = (F(2), F(18))
    for right, expected in ((phased, F(0)), (opposite, F(4))):
        p, q = sh._norm(left, normalizations[0]), sh._norm(right, normalizations[1])
        overlap = sh._inner_squared(left, right)/(normalizations[0]*normalizations[1])
        assert (p+q)**2-4*overlap == expected
    assert sh._norm(phased, F(18)) == sh._norm(opposite, F(18)) == 1

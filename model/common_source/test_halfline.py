"""Independent Hermite integrals and exact mass checks for infinite sign readout."""
from fractions import Fraction as Q
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_hermite

from . import halfline as h


def _bargmann_wavefunction(n, x):
    # Independent coordinate-space definition: sqrt(n!) times normalized
    # Hermite function. The factorial therefore cancels from this formula.
    return eval_hermite(n, x) * math.exp(-x*x/2) / (math.pi**.25 * 2**(n/2))


def test_machin_and_integer_root_enclosures_have_correct_orientation():
    low, high = h.pi_interval()
    known_low = Q("3.14159265358979323846264338327950288419716939937510")
    known_high = Q("3.14159265358979323846264338327950288419716939937511")
    assert known_low < low < high < known_high
    a, b = h.inverse_sqrt_two_pi_interval()
    assert 0 < a < b
    assert a*a * 2*high <= 1
    assert b*b * 2*low >= 1
    assert b-a <= Q(2, h.ROOT_SCALE)
    assert abs(float(a) - 1/math.sqrt(2*math.pi)) < 1e-16


@pytest.mark.parametrize("m,n", ((0, 0), (0, 1), (0, 3), (2, 1), (2, 3),
                                 (3, 2), (4, 3), (3, 3), (2, 4), (5, 6)))
def test_halfline_coefficients_match_independent_coordinate_integrals(m, n):
    constant, cross = h.halfline_coefficient(m, n)
    expected, quadrature_error = quad(
        lambda x: _bargmann_wavefunction(m, x) * _bargmann_wavefunction(n, x),
        0, np.inf, epsabs=1e-11, epsrel=1e-12)
    actual = float(constant) + float(cross)/math.sqrt(2*math.pi)
    # Quadrature is only an independent test, never used by the certificate.
    assert abs(actual-expected) < 10*quadrature_error + 1e-11


def test_every_coefficient_is_hermitian_and_complement_has_full_bargmann_norm():
    for m in range(9):
        for n in range(9):
            constant, cross = h.halfline_coefficient(m, n)
            assert (constant, cross) == h.halfline_coefficient(n, m)
            assert 2*constant == (math.factorial(m) if m == n else 0)
            if (m-n) % 2 == 0:
                assert cross == 0
    assert h.halfline_coefficient(0, 1) == (Q(0), Q(1))
    assert h.halfline_coefficient(0, 3) == (Q(0), Q(-1))
    assert h.halfline_coefficient(2, 3) == (Q(0), Q(3))


def test_d2_infinite_effect_is_not_the_finite_coordinate_sign_projector():
    # (|0>+|1>)/sqrt(2) is the positive eigenvector of finite Q_2,
    # hence its finite positive-readout probability is exactly one.
    vector = {(0, 0, 0): (Q(1), Q(0)), (1, 0, 0): (Q(1), Q(0))}
    assert h.positive_weight_coefficients(vector, Q(2)) == (Q(1, 2), Q(1))
    low, high = h.positive_weight_interval(vector, Q(2))
    assert Q(89, 100) < low < high < Q(9, 10)
    assert abs(float(1-high) - .1010577195985673) < 1e-15


def test_physical_postread_output_has_a_nonzero_tail_outside_d2():
    # On vacuum, the true plus branch has norm^2=1/2. Its two retained
    # Hermite coefficients are 1/2 and 1/sqrt(2*pi), so its omitted norm
    # is 1/4-1/(2*pi). The postread state must NOT be truncated to d=2.
    vacuum = {(0, 0, 0): (Q(1), Q(0))}
    assert h.positive_weight_interval(vacuum) == (Q(1, 2), Q(1, 2))
    _, high = h.inverse_sqrt_two_pi_interval()
    assert Q(1, 4)-high*high > Q(9, 100)


def test_complex_coherence_and_untouched_gauge_labels_have_exact_weights():
    vector = {(0, 0, 0): (Q(1), Q(2)), (1, 0, 0): (Q(3), Q(-4))}
    assert h.norm_squared(vector) == 30
    assert h.positive_weight_coefficients(vector) == (Q(15), Q(-10))
    # Changing only the gauge label removes this cross term, without
    # changing any occupation probability or total norm.
    vector[(1, 0, 1)] = vector.pop((1, 0, 0))
    assert h.positive_weight_coefficients(vector) == (Q(15), Q(0))


def test_factorial_weights_and_spectator_labels_are_not_implicitly_normalized():
    vector = {(2, 3, 0): (Q(2), Q(3))}
    assert h.norm_squared(vector) == 156
    assert h.positive_weight_interval(vector) == (Q(78), Q(78))
    assert h.positive_weight_interval(vector, Q(156)) == (Q(1, 2), Q(1, 2))
    vector[(3, 2, 0)] = (Q(1), Q(0))
    assert h.norm_squared(vector) == 168
    assert h.positive_weight_coefficients(vector) == (Q(84), Q(0))


@pytest.mark.parametrize("sign", (-1, 1))
def test_raw_complement_intervals_sum_to_mass_for_either_cross_term_sign(sign):
    vector = {(0, 0, 0): (Q(1), Q(2)), (1, 0, 0): (Q(sign*3), Q(-sign*4)),
              (2, 1, 3): (Q(2, 3), Q(1, 5)), (3, 1, 3): (Q(-1, 4), Q(2, 7))}
    normalization = Q(7, 3)
    positive = h.weight_interval(vector, 1, normalization)
    negative = h.weight_interval(vector, 0, normalization)
    mass = h.norm_squared(vector)/normalization
    assert positive[0]+negative[1] == mass
    assert positive[1]+negative[0] == mass
    assert positive[0] >= 0 and negative[0] >= 0
    assert positive[1] > 1  # Raw weights are not capped to a unit probability.


def test_arbitrary_complex_polynomial_against_independent_halfline_quadrature():
    vector = {(0, 0, 0): (Q(1, 3), Q(2, 5)), (1, 0, 0): (Q(1, 2), Q(-1, 7)),
              (2, 0, 0): (Q(1, 4), Q(1, 6)), (3, 0, 0): (Q(-1, 5), Q(2, 9)),
              (0, 2, 1): (Q(1, 6), Q(-1, 3)), (1, 2, 1): (Q(1, 8), Q(1, 9))}
    expected = 0.
    for n1, gauge in {(key[1], key[2]) for key in vector}:
        coefficients = [(n0, complex(float(z[0]), float(z[1])))
                        for (n0, other, g), z in vector.items() if (other, g) == (n1, gauge)]
        value, _ = quad(lambda x: abs(sum(z*_bargmann_wavefunction(n, x)
                                         for n, z in coefficients))**2,
                        0, np.inf, epsabs=1e-12, epsrel=1e-12)
        expected += math.factorial(n1)*value
    low, high = h.positive_weight_interval(vector)
    assert abs(float(low)-expected) < 2e-12
    assert abs(float(high)-expected) < 2e-12


def test_exact_finite_dynamics_has_nonzero_response_for_the_infinite_sign_effect():
    # This is a proved connection to the existing d=2 signal, with the
    # corrected physical observable. It makes no occupation-dynamics claim.
    from . import finite_signal

    configurations = []
    for enabled in (True, False):
        branches, _, _, eta = finite_signal.branch_vectors(enabled)
        before_sign = {}
        for (first, _), vector in branches:
            if first not in before_sign:
                before_sign[first] = vector
            else:
                before_sign[first] = tuple((a[0]+b[0], a[1]+b[1])
                                           for a, b in zip(before_sign[first], vector, strict=True))
        low = high = Q(0)
        for vector in before_sign.values():
            # B0+B1=I reconstructs the exact Taylor polynomial before the
            # finite pointer; only its readout is replaced in this comparison.
            polynomial = {((index >> 2) & 1, (index >> 3) & 1, index & 3): z
                          for index, z in enumerate(vector)}
            a, b = h.positive_weight_interval(polynomial)
            low, high = low+a, high+b
        # Each pre-sign branch has vector error <=eta. Both old records
        # remain present, so their effect-weight errors add conservatively.
        radius = 2*(2*eta+eta*eta)
        configurations.append((low-radius, high+radius))
    low = configurations[0][0]-configurations[1][1]
    high = configurations[0][1]-configurations[1][0]
    assert Q("-0.025262207") < low <= high < Q("-0.025262206")


def test_empty_vector_and_zero_coefficients_are_valid_zero_branches():
    for vector in ({}, {(0, 0, 0): (Q(0), Q(0))}):
        assert h.norm_squared(vector) == 0
        assert h.positive_weight_coefficients(vector) == (Q(0), Q(0))
        assert h.positive_weight_interval(vector) == (Q(0), Q(0))
        assert h.weight_interval(vector, 0) == (Q(0), Q(0))


@pytest.mark.parametrize("vector", (
    [], None, {(True, 0, 0): (Q(1), Q(0))}, {(-1, 0, 0): (Q(1), Q(0))},
    {(0, h.MAX_OCCUPATION+1, 0): (Q(1), Q(0))}, {(0, 0, True): (Q(1), Q(0))},
    {(0, 0, 4): (Q(1), Q(0))}, {(0, 0): (Q(1), Q(0))},
    {(0, 0, 0): (1, Q(0))}, {(0, 0, 0): (True, Q(0))},
    {(0, 0, 0): (float("nan"), Q(0))}, {(0, 0, 0): (Q(1), float("inf"))},
    {(0, 0, 0): [Q(1), Q(0)]}, {(0, 0, 0): (Q(1),)},
))
def test_invalid_vectors_rejected(vector):
    with pytest.raises(ValueError):
        h.norm_squared(vector)
    with pytest.raises(ValueError):
        h.positive_weight_interval(vector)


@pytest.mark.parametrize("normalization", (0, 1, True, 1., Q(0), Q(-1), float("nan")))
def test_normalization_requires_an_explicit_positive_rational(normalization):
    with pytest.raises(ValueError):
        h.positive_weight_interval({}, normalization)


@pytest.mark.parametrize("value", (True, -1, 1., h.MAX_OCCUPATION+1, float("nan")))
def test_invalid_matrix_indices_rejected_even_after_cache_is_warm(value):
    h.halfline_coefficient(1, 0)
    with pytest.raises(ValueError):
        h.halfline_coefficient(value, 0)
    with pytest.raises(ValueError):
        h.halfline_coefficient(0, value)


@pytest.mark.parametrize("outcome", (True, 2, -1, 0.))
def test_invalid_sign_records_rejected(outcome):
    with pytest.raises(ValueError):
        h.weight_interval({}, outcome)

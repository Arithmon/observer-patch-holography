"""Original-span controls for projection and cancellation-sensitive distance."""

from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
import pytest

from quantum_information.algebras import (
    ALGEBRA_ATOL, FiniteAlgebra, algebra_intersection, matrix_units,
)


def _sqrt_fraction(value):
    # An independent decimal evaluation is ample for these non-midpoint
    # controls; the producer's norm and rounding routines are not reused.
    with localcontext() as context:
        context.prec = 120
        return float((Decimal(value.numerator)/Decimal(value.denominator)).sqrt())


@pytest.mark.parametrize("phase", (1., 1j))
def test_tiny_distance_survives_a_large_original_scalar_component(phase):
    algebra = FiniteAlgebra([np.eye(3)])
    scalar, perturbation = 2.**600, 2.**-600
    x = phase*scalar*np.eye(3, dtype=complex)
    x[0, 1] = phase*perturbation
    # The original span is exactly CI. Rounded QR factors do not define it:
    # their small orthogonality error would dominate this true residual.
    np.testing.assert_array_equal(algebra.project(x), phase*scalar*np.eye(3))
    assert algebra.distance(x) == perturbation


@pytest.mark.parametrize("phase", (1., 1j))
def test_distance_is_not_the_norm_of_a_separately_rounded_projection(phase):
    algebra = FiniteAlgebra([np.eye(3)])
    upper = np.nextafter(1., np.inf)
    x = phase*np.diag([1., upper, 1.])
    delta = Fraction(float(upper))-1
    # The exact projection is (1+delta/3)I, which rounds to I. Its exact
    # residual is delta*(-1,2,-1)/3 and has squared norm 2*delta**2/3.
    expected = _sqrt_fraction(Fraction(2, 3)*delta**2)
    np.testing.assert_array_equal(algebra.project(x), phase*np.eye(3))
    assert algebra.distance(x) == expected
    assert expected < float(delta)


@pytest.mark.parametrize("phase", (1., 1j))
def test_subnormal_multicomponent_distance_rounds_after_the_norm(phase):
    tiny = np.nextafter(0., 1.)
    algebra = FiniteAlgebra([np.eye(4)])
    x = np.zeros((4, 4), dtype=complex)
    for i, j in ((0, 1), (0, 2), (1, 0), (1, 2)):
        x[i, j] = phase*tiny
    # Four orthogonal entries give sqrt(4*tiny**2)=2*tiny. Sequential
    # unscaled hypot instead rounds back to tiny after every addition.
    np.testing.assert_array_equal(algebra.project(x), np.zeros((4, 4)))
    assert algebra.distance(x) == 2*tiny


@pytest.mark.parametrize("phase", (1., 1j))
def test_distance_succeeds_when_only_intermediate_components_underflow(phase):
    tiny = np.nextafter(0., 1.)
    algebra = FiniteAlgebra([np.eye(2)])
    x = phase*np.diag([tiny, 0.])
    # Projection and residual coordinates have size tiny/2, but the final
    # distance tiny/sqrt(2) rounds to tiny. Refuse the unrepresentable matrix,
    # while retaining the representable scalar diagnostic.
    assert algebra.distance(x) == tiny
    with pytest.raises(ValueError, match="underflow|precision|represent|range"):
        algebra.project(x)


def test_intersection_projection_cache_uses_the_retained_original_span():
    identity = np.eye(2)
    n = np.array([[1., 3.-1j], [3.+1j, -1.]])
    algebra = FiniteAlgebra([identity, n])
    full = FiniteAlgebra(matrix_units(2))
    # N^2=11I, so CI+CN is an exact algebra. Converting its intersection
    # basis to a float chart introduces thirds; that chart is not the exact
    # retained intersection. Constructor-time cached projections must not
    # survive replacement by the retained exact span.
    np.testing.assert_array_equal(n@n, 11*identity)
    common = algebra_intersection([algebra, full])
    assert common.distance(n) == 0.
    np.testing.assert_array_equal(common.project(n), n)
    repeated = algebra_intersection([common, full, algebra])
    assert repeated.distance(n) == 0.
    np.testing.assert_array_equal(repeated.project(n), n)


@pytest.mark.parametrize("ulps, accepted", ((900719, True), (900720, False)))
def test_constructor_identity_tolerance_retains_both_sides(ulps, accepted):
    a = np.diag([1., 1.+ulps*2.**-52])
    entry = Fraction(float(a[1, 1]))
    # In span(diag(1,b)), minimization of (1-c)^2+(1-b*c)^2
    # gives identity-distance squared (b-1)^2/(1+b^2). These adjacent
    # supplied values straddle the documented constructor tolerance.
    squared_distance = (entry-1)**2/(1+entry**2)
    threshold = Fraction(float(ALGEBRA_ATOL*np.sqrt(2.)))
    assert (squared_distance <= threshold**2) is accepted
    if accepted:
        assert FiniteAlgebra([a]).dimension == 1
    else:
        with pytest.raises(ValueError, match="contain the identity"):
            FiniteAlgebra([a])

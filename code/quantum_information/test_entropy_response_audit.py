"""Maintainer-audit regressions: retain coherent reference data and slopes."""

import math
from fractions import Fraction

import numpy as np
import pytest

from quantum_information.entropy_response import (
    entropy_tangent, finite_entropy_balance, first_law_diagnostic,
    gibbs_entropy_response,
)

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])


@pytest.mark.parametrize("epsilon", [1e-10, 1e-20, 1e-100])
def test_coherent_reference_retains_its_analytic_tangent(epsilon):
    sigma = I/2+epsilon*X
    # Exact 2-by-2 divided difference: log(sigma)'s X coefficient.
    expected = -math.log1p(4*epsilon/(1-2*epsilon))
    assert entropy_tangent(sigma, X) == pytest.approx(expected, rel=3e-14, abs=0)


@pytest.mark.parametrize("epsilon", [1e-10, 1e-20, 1e-100])
def test_constant_generator_cannot_pass_for_a_coherent_reference(epsilon):
    sigma = I/2+epsilon*X
    expected = math.log1p(4*epsilon/(1-2*epsilon))/math.sqrt(2)
    report = first_law_diagnostic(sigma, np.zeros((2, 2)))
    assert report["all_tangent_defect"] == pytest.approx(expected, rel=3e-14, abs=0)


@pytest.mark.parametrize("epsilon", [1e-10, 1e-20, 1e-100])
def test_finite_entropy_difference_keeps_the_reference_modular_response(epsilon):
    report = finite_entropy_balance(I/2+2*epsilon*X, I/2+epsilon*X)
    # S(I/2+2eX)-S(I/2+eX) = -6e^2+O(e^4); relative correction <1e-18 here.
    assert report["entropy_change"] == pytest.approx(-6*epsilon**2, rel=3e-14, abs=0)


@pytest.mark.parametrize("scale", [1e-160, 1e-161])
def test_subnormal_variance_cannot_report_a_false_gibbs_slope(scale):
    with pytest.raises(ValueError, match="resolved|precision"):
        gibbs_entropy_response(scale*np.diag([1., -1.]), 1.3)


@pytest.mark.parametrize("probability", [.01, .4, .5])
@pytest.mark.parametrize("epsilon", [1e-20, 1e-100, 1e-300])
def test_unequal_reference_weights_keep_tiny_complex_coherence(probability, epsilon):
    y = np.array([[0., -1j], [1j, 0.]])
    sigma = np.diag([probability, 1-probability])+epsilon*y
    if probability == .5:
        coefficient = 2.
    else:
        difference = 2*probability-1
        coefficient = math.log1p(difference/(1-probability))/difference
    # The divided difference is exact at first order; O(epsilon^3) is below
    # the tolerance for these inputs. No production matrix log is reused.
    expected = -2*epsilon*coefficient
    assert entropy_tangent(sigma, y) == pytest.approx(expected, rel=3e-14, abs=0)


def test_reference_log_does_not_change_global_mpmath_precision():
    import mpmath
    original = mpmath.mp.dps
    entropy_tangent(I/2+1e-100*X, X)
    assert mpmath.mp.dps == original


@pytest.mark.parametrize("epsilon", [1e-20, 1e-100, 1e-150])
def test_three_level_reference_keeps_its_quadratic_diagonal_response(epsilon):
    # Adjacency of a three-site path has eigenvalues 0,+sqrt(2),-sqrt(2).
    # For D=diag(1,-2,1), Tr(log(sigma) D)=-log(1-2(e/a)^2)/2.
    a = 1/3
    adjacency = np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]])
    sigma = a*np.eye(3)+epsilon*adjacency
    expected = math.log1p(-2*(epsilon/a)**2)/2
    assert entropy_tangent(sigma, np.diag([1., -2., 1.])) == pytest.approx(expected, rel=3e-14, abs=0)


def test_longer_coherent_path_retains_a_representable_fourth_order_response():
    a, epsilon = .2, 1e-80
    adjacency = np.diag(np.ones(4), 1)+np.diag(np.ones(4), -1)
    direction = np.zeros((5, 5))
    direction[0, 4] = direction[4, 0] = 1
    # The unique length-four path gives log(sigma)[0,4]=-(epsilon/a)^4/4
    # to leading order. Higher orders are below the binary64 range here.
    expected = float(Fraction(epsilon)**4/(2*Fraction(a)**4))
    actual = entropy_tangent(a*np.eye(5)+epsilon*adjacency, direction)
    assert actual == pytest.approx(expected, abs=1e-323, rel=0)


def test_exact_polynomial_zero_is_not_confused_with_a_tiny_nonzero_coherence():
    q = np.kron(X, I)+np.kron(np.diag([1., -1.]), X)
    np.testing.assert_array_equal(q@q, 2*np.eye(4))
    sigma = np.eye(4)/4+q/16
    direction = np.zeros((4, 4))
    direction[0, 3] = direction[3, 0] = 1
    # q^2=2I implies log(sigma)=aI+bq. The selected entry is exactly zero,
    # even though sigma's coordinate graph is connected.
    assert entropy_tangent(sigma, direction) == 0
    # Adding an actual tiny entry invalidates that exact zero certificate.
    epsilon = 1e-100
    assert entropy_tangent(sigma+epsilon*direction, direction) < -epsilon


def test_unrepresentable_true_response_is_not_declared_an_exact_zero():
    adjacency = np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]])
    sigma = np.eye(3)/3+1e-200*adjacency
    # The diagonal response is of order 1e-400, outside binary64. Its exact
    # power moments are nonzero, so the polynomial-zero certificate fails.
    with pytest.raises(ValueError, match="resolved|precision"):
        entropy_tangent(sigma, np.diag([1., -2., 1.]))

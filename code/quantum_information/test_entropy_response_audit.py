"""Maintainer-audit regressions: retain coherent reference data and slopes."""

import math

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


def test_subnormal_variance_cannot_report_a_false_gibbs_slope():
    with pytest.raises(ValueError, match="resolved|precision"):
        gibbs_entropy_response(1e-160*np.diag([1., -1.]), 1.3)

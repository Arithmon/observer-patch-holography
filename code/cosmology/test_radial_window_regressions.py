"""Independent controls for the finite-window theorem's numerical evaluation."""
from __future__ import annotations

import math

import mpmath
import numpy as np
import pytest

from oph_radial_lift_330 import (
    derivative_mellin_norm,
    finite_window_stability_bound,
    normalized_radial_window,
    window_transfer,
)


def high_precision_norms(ell: int, theta: float):
    """Use the original gamma integrals, not the producer's recurrence."""
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    t = ctx.mpf(theta)

    def integral(exponent):
        return (ctx.sqrt(ctx.pi) / 4 * ctx.gamma(1 + exponent / 2)
                / ctx.gamma(ctx.mpf(3) / 2 + exponent / 2)
                * ctx.gamma(ell - exponent / 2)
                / ctx.gamma(ell + 2 + exponent / 2))

    i = integral(t)
    j = integral(t - 2) - (ell * (ell + 1) - t * (t + 1) / 2) * i
    return ctx, i, j


@pytest.mark.parametrize("theta", [1e-20, 1e-16, 1e-14, 0.03397850436258248, 1.9, 3.9])
def test_derivative_norm_retains_positive_near_scale_invariant_domain(theta):
    _, _, expected = high_precision_norms(2, theta)
    assert derivative_mellin_norm(2, theta) == pytest.approx(float(expected), rel=2e-13)


@pytest.mark.parametrize("theta", [1e-20, 1e-16, 1e-14, 0.03397850436258248])
def test_window_bound_preserves_nonzero_radius_displacement(theta):
    ctx, i, j = high_precision_norms(2, theta)
    r = 1.01
    t = ctx.mpf(theta)
    # Direct high-precision subtraction independently checks expm1 in the
    # producer; the near-zero float result was a false zero uncertainty.
    eta = 2 * ctx.sqrt(j) / t * abs(ctx.mpf(r) ** (t / 2) - 1)
    expected = 4 * ctx.pi * eta * (2 * ctx.sqrt(i) + eta)
    actual = finite_window_stability_bound(
        2, theta, A_zeta=1, Z_q=1, k_pivot=1, R_star=1,
        radii=[r], radial_weights=[1],
    )
    assert actual.absolute_cl_bound > 0
    assert actual.absolute_cl_bound == pytest.approx(float(expected), rel=2e-13)


def test_window_weight_normalization_and_transfer_are_scale_invariant():
    radii = [0.9, 1.1]
    _, normalized = normalized_radial_window(radii, [1e308, 1e308])
    np.testing.assert_array_equal(normalized, [0.5, 0.5])
    k = np.geomspace(0.01, 100, 100)
    expected = window_transfer(2, k, radii, [1, 1])
    actual = window_transfer(2, k, radii, [1e308, 1e308])
    np.testing.assert_array_equal(actual, expected)
    assert np.any(actual != 0)


def test_shifted_shell_gives_a_strictly_positive_analytic_spectrum_difference():
    # For a one-point window, the exact shell law scales as r**theta.
    # Evaluate that actual difference without subtracting two rounded C_l's.
    theta, radius = 1e-14, 1.01
    ctx, i, _ = high_precision_norms(2, theta)
    exact_difference = 4 * ctx.pi * i * ctx.expm1(ctx.mpf(theta) * ctx.log(radius))
    bound = finite_window_stability_bound(
        2, theta, A_zeta=1, Z_q=1, k_pivot=1, R_star=1,
        radii=[radius], radial_weights=[1],
    )
    assert 0 < exact_difference < bound.absolute_cl_bound
    assert math.isfinite(bound.absolute_cl_bound)

"""Independent Gibbs-response controls: no derivative-ratio self-checks."""

from fractions import Fraction

import mpmath
import numpy as np
import pytest

from maxent.information_projection import duhamel_covariance


@pytest.mark.parametrize("small", [1e-20, 1e-100, 1e-300])
def test_uniform_qutrit_retains_small_cross_response(small):
    a = np.diag([1., -1., small])
    b = np.diag([1., 1., -2.])
    result = duhamel_covariance([a, b], [0., 0.])
    # Exact uniform-state trace pairing, independent of thermal numerics.
    expected = float(-2*Fraction(small)/3)
    assert result[0, 1] == pytest.approx(expected, rel=2e-14, abs=0)
    assert result[1, 0] == result[0, 1]


@pytest.mark.parametrize("gap", [700., 710., 744., 745., 750., 1000.])
@pytest.mark.parametrize("scale", [1., 1e150])
@pytest.mark.parametrize("sign", [-1., 1.])
def test_rare_population_moments_are_evaluated_before_binary64_rounding(gap, scale, sign):
    beta = sign*gap/scale
    ctx = mpmath.mp.clone()
    ctx.dps = 500
    e, b = ctx.mpf(scale), ctx.mpf(beta)
    w = ctx.exp(-abs(b*e))
    expected = e*e*w/(1+w)**2
    rounded = float(expected)
    if rounded == 0 or abs(ctx.mpf(rounded)-expected) > 32*np.finfo(float).eps*expected:
        with pytest.raises(ValueError, match="precision|underflow"):
            duhamel_covariance([np.diag([0., scale])], [beta])
    else:
        result = duhamel_covariance([np.diag([0., scale])], [beta])
        assert result[0, 0] == pytest.approx(rounded, rel=2e-13, abs=0)


def test_exact_zero_cross_response_is_distinct_from_lost_output():
    result = duhamel_covariance([np.diag([1., -1., 0.]), np.diag([1., 1., -2.])], [0., 0.])
    assert result[0, 1] == result[1, 0] == 0
    with pytest.raises(ValueError, match="precision|underflow"):
        duhamel_covariance([1e-200*np.diag([1., -1.])], [0.])

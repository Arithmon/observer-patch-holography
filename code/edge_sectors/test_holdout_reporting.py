"""Public fit and reporting boundaries, separately from ground-state solvers."""
from fractions import Fraction
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge_sectors.heat_kernel_holdout_validation import fit_t, predict, report_zn


@pytest.mark.parametrize("weights", [(1., 1e-320, 1., 4.),
    (.5, np.nextafter(.5, 0.), 1., 3.), (1e-200, 1e-200, 2., 3.)])
def test_fitted_time_preserves_extreme_and_near_equal_weights(weights):
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    p0, p1, dim, eigenvalue = map(ctx.mpf, weights)
    expected = float(ctx.log(p0*dim/p1)/eigenvalue)
    assert fit_t(*weights) == pytest.approx(expected, rel=5e-13, abs=0)


def test_uniform_weights_have_zero_fitted_time_without_false_division():
    assert fit_t(.5, .5, 1., 4.) == 0
    # The full strictly positive-h model is almost uniform; its computed
    # probability contrast cannot support a normalized residual here.
    with pytest.raises(ValueError, match="log-gap is unresolved"):
        report_zn(5, [1e-20])


@pytest.mark.parametrize("bad", [True, np.bool_(False), -1., 0., float("inf"),
    float("nan"), 1+0j, np.ma.array(1., mask=True), 2**53+1,
    np.uint64(2**53+1), Fraction(1, 3), 10**400])
def test_fit_rejects_invalid_or_lossy_original_parameters(bad):
    with pytest.raises(ValueError):
        fit_t(1., bad, 1., 4.)


def test_prediction_avoids_premature_exponential_underflow():
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    expected = float(ctx.mpf(1e-100)*ctx.mpf(1e200)*ctx.exp(-800))
    assert predict(1e-100, 1e200, 1., 800.) == pytest.approx(expected, rel=5e-13, abs=0)


@pytest.mark.parametrize("time", [1e308, -1e308])
def test_unrepresentable_predictions_are_not_zero_or_infinite(time):
    with pytest.raises(ValueError, match="reporting range"):
        predict(.5, 1., 4., time)

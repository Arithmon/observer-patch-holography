"""Checks for calibration, signed residuals, and the retained data comparison."""
import numpy as np
import pytest

import measured_comparison as comparison


def test_calibration_is_applied_once_to_theory_only():
    # Deliberately large calibration distinguishes division from multiplication,
    # and one division from a repeated conversion.
    data = np.array([[30., 7., 2., 3.], [31., 12., 4., 5.]])
    original = data.copy()
    result = comparison.exact_residuals(data, np.array([30, 31]), np.array([20., 64.]), 2.)
    np.testing.assert_array_equal(data, original)
    np.testing.assert_array_equal(result[:, 1], [7., 12.])
    np.testing.assert_array_equal(result[:, 2], [5., 16.])
    np.testing.assert_array_equal(result[:, 3:5], original[:, 2:4])
    np.testing.assert_allclose(result[:, -1], [1., -.8])


def test_negative_te_uses_signed_difference_and_error_toward_model():
    data = np.array([[30., -7., 2., 5.], [31., -6., 2., 5.]])
    result = comparison.exact_residuals(data, np.array([30, 31]), np.array([-20., -40.]), 2.)
    # Model [-5,-10]; the first residual points upward to the model (error_plus),
    # the second downward (error_minus). Dividing by TE itself would be wrong.
    np.testing.assert_array_equal(result[:, 2], [-5., -10.])
    np.testing.assert_allclose(result[:, -1], [-.4, 2.])


@pytest.mark.parametrize("multipoles", [[30.5], [31.]])
def test_effective_bin_centers_or_missing_integer_multipoles_are_rejected(multipoles):
    data = np.column_stack((multipoles, np.ones((len(multipoles), 3))))
    with pytest.raises(ValueError, match="exact matching multipoles"):
        comparison.exact_residuals(data, np.array([30, 32]), np.array([4., 4.]), 1.)


@pytest.mark.parametrize("calibration", [0., -1., np.nan, np.inf])
def test_invalid_calibration_is_rejected(calibration):
    with pytest.raises(ValueError, match="Calibration"):
        comparison.exact_residuals(np.array([[30., 1., 1., 1.]]), np.array([30]), np.array([1.]), calibration)


def test_retained_comparison_matches_independently_computed_ranges_and_residual_rms():
    columns, tables, calibration, residuals, receipt = comparison.prepare()
    expected = {
        "TT": (2471, 2500, 1.0136332508133403),
        "TE": (1967, 1996, 1.0197569815778813),
        "EE": (1967, 1996, 1.019438700917485),
    }
    assert calibration == 1.000442
    for channel, (count, maximum, rms) in expected.items():
        rows = residuals[channel]
        assert rows.shape == (count, 6)
        np.testing.assert_array_equal(rows[:, 0], np.arange(30, maximum + 1))
        assert receipt["diagnostics"][channel]["rms_in_quoted_error_units"] == pytest.approx(rms, abs=1e-12)
        # Independent frozen coverage uses ell-2, not the producer's searchsort.
        raw = columns[f"mean_D_{channel}_uK2"][28:maximum-1]
        np.testing.assert_allclose(rows[:, 2], raw * .9991165857467871, rtol=1e-14)
        full = tables[channel, "full"]
        selected = full[(full[:, 0] >= 30) & (full[:, 0] <= maximum)]
        np.testing.assert_array_equal(rows[:, 1], selected[:, 1])
        assert not np.array_equal(rows[:, 0], tables[channel, "binned"][:, 0])

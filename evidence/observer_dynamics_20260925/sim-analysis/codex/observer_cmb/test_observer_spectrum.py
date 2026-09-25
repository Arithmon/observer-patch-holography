"""Independent statistical checks and false-green mutations for the readout."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("conditional_observer_spectrum", HERE / "observer_spectrum.py")
observer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(observer)


def test_d_to_c_normalization_and_temperature_variance():
    ell = np.array([2, 3])
    d = np.array([[6., -2., 3.], [12., 4., 8.]])
    c = observer.d_to_c(ell, d)
    np.testing.assert_allclose(c[:, 0], [2 * np.pi, 2 * np.pi])
    assert c[0, 1] == pytest.approx(-2 * np.pi / 3)
    variance = np.sum((2 * ell + 1) * c[:, 0]) / (4 * np.pi)
    assert variance == pytest.approx(6.)


@pytest.mark.parametrize("ell,power", [
    ([2, 4], [[4, 0, 9], [4, 0, 9]]),
    ([2, 2], [[4, 0, 9], [4, 0, 9]]),
    ([2.5], [[4, 0, 9]]),
    ([1], [[4, 0, 9]]),
    ([2], [[4, 6.1, 9]]),
    ([2], [[4, -6.1, 9]]),
    ([2], [[-4, 0, 9]]),
    ([2], [[0, 1e-30, 9]]),
    ([2], [[4, np.nan, 9]]),
])
def test_invalid_multipoles_and_covariance_rejected(ell, power):
    with pytest.raises(ValueError):
        observer.validate_spectra(ell, power)


def test_covariance_hand_fixture_preserves_negative_te():
    # Sigma=[[4,-3],[-3,9]], 2ell+1=5. Isserlis' fourth moments.
    expected = np.array([[32, -24, 18], [-24, 45, -54], [18, -54, 162]]) / 5
    actual = observer.gaussian_estimator_covariance([2], [[4, -3, 9]])[0]
    np.testing.assert_allclose(actual, expected)
    assert np.linalg.eigvalsh(actual).min() > 0


def test_covariance_against_independent_real_harmonic_monte_carlo():
    # Generate the actual five independent real harmonic pairs, without the
    # producer's Bartlett algorithm, then form the three quadratic estimators.
    rng = np.random.Generator(np.random.PCG64(1729))
    u, v = rng.standard_normal((2, 60000, 5))
    t = 2 * u
    e = -1.5 * u + np.sqrt(6.75) * v
    estimates = np.column_stack(((t*t).mean(axis=1), (t*e).mean(axis=1), (e*e).mean(axis=1)))
    np.testing.assert_allclose(estimates.mean(axis=0), [4, -3, 9], rtol=.012)
    expected = np.array([[32, -24, 18], [-24, 45, -54], [18, -54, 162]]) / 5
    np.testing.assert_allclose(np.cov(estimates, rowvar=False, ddof=1), expected, rtol=.035)


def test_bartlett_draw_joint_moments_and_psd():
    # At all ell, sqrt(2ell+1)*(estimator-mean) has the same covariance.
    # Combining different df checks the normalization without huge sky arrays.
    ell = np.arange(2, 20002)
    mean = np.tile([4., -3., 9.], (len(ell), 1))
    sample = observer.draw_wishart_spectra(ell, mean, seed=9182)
    deviations = (sample - mean) * np.sqrt(2 * ell + 1)[:, None]
    expected = np.array([[32, -24, 18], [-24, 45, -54], [18, -54, 162]])
    assert np.all(np.abs(deviations.mean(axis=0)) < .04 * np.sqrt(expected.diagonal()))
    np.testing.assert_allclose(np.cov(deviations, rowvar=False), expected, rtol=.065)
    assert (sample[:, 1]**2 <= sample[:, 0] * sample[:, 2]).all()
    np.testing.assert_array_equal(sample, observer.draw_wishart_spectra(ell, mean, seed=9182))


def test_singular_positive_semidefinite_input_is_supported():
    ell = [2, 3, 4]
    mean = np.array([[4., -6., 9.], [0., 0., 9.], [0., 0., 0.]])
    result = observer.draw_wishart_spectra(ell, mean)
    assert result[0, 1] == pytest.approx(-1.5 * result[0, 0])
    assert result[0, 2] == pytest.approx(2.25 * result[0, 0])
    assert result[1, 0] == result[1, 1] == 0
    assert result[1, 2] > 0
    assert not result[2].any()


def test_low_ell_marginal_interval_is_asymmetric_and_not_parameter_posterior():
    bounds = observer.auto_sampling_interval([2], [10.], .95)[0]
    # Independent published chi-square(df=5) quantiles: .025 and .975.
    np.testing.assert_allclose(bounds, [1.662423226973325, 25.665003988060054], rtol=1e-12)
    assert bounds[1] - 10 > 10 - bounds[0]
    assert bounds[0] > 0


def test_imported_blackbody_peak_is_frequency_radiance_peak():
    output = observer.blackbody_summary(2.7255)
    assert output["frequency_radiance_peak_GHz"] == pytest.approx(160.2301215, abs=.0000001)
    assert "Imported" in output["temperature_status"]


def test_pin_rejects_mutated_spectrum_even_with_matching_mutated_receipt(tmp_path):
    spectrum = tmp_path / "spectrum.txt"
    spectrum.write_bytes(observer.SPECTRUM.read_bytes().replace(b"1.016163750489", b"1.116163750489", 1))
    receipt = json.loads(observer.RECEIPT.read_text())
    receipt["spectrum_sha256"] = observer.sha256(spectrum)
    receipt_path = tmp_path / "receipt.json"
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="Pinned spectrum"):
        observer.load_inputs(spectrum, receipt_path, observer.INI)


def test_outputs_recompute_and_reject_mutations(tmp_path):
    output = observer.write_output(tmp_path)
    observer.verify_output(tmp_path)
    table = np.loadtxt(tmp_path / "observer_spectrum.txt")
    np.testing.assert_allclose(table, np.column_stack(list(output["columns"].values())), rtol=1e-13)
    assert len(output["columns"]["ell"]) == 2499
    assert output["summary"]["temperature_rms_uK"] == pytest.approx(112.46628, abs=.00001)
    assert len(output["summary"]["TT_peaks"]) == 3
    assert output["summary"]["monopole"]["temperature_K"] == 2.7255
    output["columns"]["mean_D_TT_uK2"][0] *= 1.01
    (tmp_path / "observer_spectrum.json").write_text(json.dumps(output))
    with pytest.raises(ValueError, match="numeric mismatch"):
        observer.verify_output(tmp_path)
    observer.write_output(tmp_path)
    path = tmp_path / "observer_spectrum.txt"
    path.write_text(path.read_text().replace("mean_D_TT_uK2", "mean_C_TT_uK2", 1))
    with pytest.raises(ValueError, match="header mismatch"):
        observer.verify_output(tmp_path)


def test_retained_outputs_pass_recomputation():
    observer.verify_output()

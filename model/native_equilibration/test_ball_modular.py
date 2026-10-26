"""Independent controls for the proper-ball scalar comparison."""
from fractions import Fraction as F
from copy import deepcopy
import hashlib
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import ball_modular as ball


@pytest.fixture(scope="module")
def result():
    return ball.run()


def test_proper_frame_is_same_native_read_action(result):
    read = result["proper_frame"]
    assert read["native_chart_read"]["ports"] == 960
    assert read["native_chart_read"]["native_seams"] == 2520
    assert read["native_chart_read"]["averaged_screen_tensor"] == [
        ["5039/181440", "-5039/362880"], ["-5039/362880", "5039/181440"]]
    matrix = lambda key: s.Matrix([[s.sympify(v) for v in row] for row in read[key]])
    E = matrix("proper_coordinate_matrix")
    assert s.simplify(E.T*s.diag(-1, 1, 1, 1)*E-matrix("constant_comparison_metric")) == s.zeros(4)
    assert not read["finite_average_identified_with_native_continuum_tensor"]
    assert read["proper_action_density_exact"] == "diag(-1,1,1,1)"


def test_odd_image_vacuum_equals_independent_sine_spectrum():
    X, P, sites, _ = ball.covariance(F(1, 2), F(4))
    otherX, otherP = ball.direct_covariance(F(1, 2), F(4))
    np.testing.assert_allclose(X, otherX, atol=2e-14, rtol=0)
    np.testing.assert_allclose(P, otherP, atol=2e-14, rtol=0)
    assert len(sites) == 27
    assert np.linalg.eigvalsh(X)[0] > 0
    assert np.linalg.eigvalsh(P)[0] > 0


def test_partial_trace_does_not_condition_the_exterior():
    X, P, sites, _ = ball.covariance(F(1, 2), F(4))
    distances = np.sum(abs(sites[:, None, :]-sites[None, :, :]), axis=2)
    V = 24*np.eye(len(sites))-4*(distances == 1)
    values, U = np.linalg.eigh(V)
    conditional_X, conditional_P = (U/(2*np.sqrt(values)))@U.T, (U*(np.sqrt(values)/2))@U.T
    assert np.linalg.norm(X-conditional_X) > .01
    assert np.linalg.norm(P-conditional_P) > .01
    with pytest.raises(ValueError, match="purity gaps unresolved"):
        ball.modular_momentum(conditional_X, conditional_P)
    assert ball.modular_momentum(X, P)["minimum_symplectic_gap_squared"] > 1e-6


def test_single_mode_modular_form_and_exact_coherent_relative_entropy():
    nu, omega = .7, 2.
    X, P = np.array([[nu/omega]]), np.array([[nu*omega]])
    epsilon = math.log((nu+.5)/(nu-.5))
    hp = ball.modular_momentum(X, P)["momentum_form"][0, 0]
    assert hp == pytest.approx(epsilon/omega, abs=2e-14)
    assert ball.charge(X, P, [1.]) == pytest.approx(epsilon/(2*omega), abs=2e-14)
    # Independent finite-Fock check of D(U sigma U*||sigma)=Delta<K>:
    # conjugation preserves the entire reference spectrum and entropy.
    cutoff, amplitude = 30, .2
    lowering = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    q = (lowering+lowering.T)/math.sqrt(2*omega)
    thermal = np.exp(-epsilon*np.arange(cutoff)); thermal /= sum(thermal)
    sigma = np.diag(thermal)
    U = expm(1j*amplitude*q)
    rho = U@sigma@U.conj().T
    log_sigma = np.diag(np.log(thermal))
    log_rho = U@log_sigma@U.conj().T
    relative_entropy = float(np.trace(rho@(log_rho-log_sigma)).real)
    assert relative_entropy/amplitude**2 == pytest.approx(ball.charge(X, P, [1.]), abs=2e-12)


def test_noise_entropy_secants_bracket_actual_modular_derivative(result):
    fine = next(row for row in result["grid_controls"] if row["proper_spacing"] == "2/7")
    for row in fine["profiles"]:
        tangent = row["modular_entropy_derivative"]
        controls = row["entropy_secant_control"]
        assert controls["step"] == 1e-6
        assert controls["negative_step_below_one_twentieth_numerical_positivity_radius"]
        # Entropy is concave along the affine covariance-noise family. The
        # slack covers floating entropy evaluation, not a physical error fit.
        assert controls["positive_noise_slope"] <= tangent+3e-7
        assert tangent <= controls["regional_negative_noise_slope"]+3e-7
        assert abs(controls["positive_noise_slope"]-tangent) < .001


def test_all_grid_controls_retained_with_precision_refusal(result):
    grids = result["grid_controls"]
    assert [row["proper_spacing"] for row in grids] == ["2/3", "1/2", "2/5", "1/3", "2/7", "1/4"]
    assert grids[-1]["status"] == "REFUSED"
    assert "purity gaps unresolved" in grids[-1]["reason"]
    for row in grids[:-1]:
        assert row["all_symplectic_modes_retained"]
        assert row["proper_lattice_spacing_isotropic"]
        assert row["minimum_regional_symplectic_gap_squared"] >= ball.MIN_PURITY_GAP
    assert abs(grids[-2]["profiles"][0]["relative_discrepancy"]) < .04
    assert abs(grids[-2]["profiles"][1]["relative_discrepancy"]) < .07
    # Pixel boundaries give both signs; the controls do not claim monotonic
    # or proved continuum convergence of these numerical discrepancies.
    assert grids[0]["profiles"][0]["relative_discrepancy"] > 0
    assert grids[1]["profiles"][0]["relative_discrepancy"] < 0


def test_infrared_box_effect_is_measured(result):
    controls = result["infrared_controls"]
    charges = [row["profiles"][0]["modular_entropy_derivative"] for row in controls]
    assert [row["proper_box_length"] for row in controls] == ["4", "8", "12", "16"]
    assert all(a > b for a, b in zip(charges, charges[1:]))
    assert charges[0]-charges[1] > charges[1]-charges[2] > charges[2]-charges[3] > 0
    assert not result["scope"]["interval_arithmetic_or_proven_continuum_error_rate"]


def test_nested_regions_use_one_identical_vacuum(result):
    rows = result["same_vacuum_nested_ball_relative_entropy_controls"]
    values = [row["constant_coherent_relative_entropy_per_amplitude_squared"] for row in rows]
    assert all(a < b for a, b in zip(values, values[1:]))
    # Independently recover the smaller covariance by actual partial trace.
    largeX, largeP, large_sites, _ = ball.covariance(F(2, 5), F(8), F(6, 5))
    smallX, smallP, small_sites, _ = ball.covariance(F(2, 5), F(8), F(4, 5))
    lookup = {tuple(site): i for i, site in enumerate(large_sites)}
    indices = [lookup[tuple(site)] for site in small_sites]
    np.testing.assert_array_equal(largeX[np.ix_(indices, indices)], smallX)
    np.testing.assert_array_equal(largeP[np.ix_(indices, indices)], smallP)


def test_conformal_kernel_and_area_coefficient_independently_derived(result):
    row = result["exact_ball_identities"]
    assert row["constant_momentum_noise_entropy_kernel"] == "4*pi**2*R**4/15"
    assert row["linear_momentum_noise_entropy_kernel"] == "4*pi**2*R**4/105"
    assert row["fixed_volume_area_derivative_per_k"] == "-4*pi*R**4/5"
    assert row["fixed_volume_area_derivative_per_G_tau_tau"] == "-4*pi*R**4/15"
    assert row["space_form_G_tau_tau"] == "3*k"
    assert not row["geometry_defined_from_stress_or_entropy"]
    assert not row["Newton_constant_or_renormalized_area_coefficient_derived"]
    # The2*pi weight and factor1/2 in independently measured kinetic energy
    # are fixed before the Gaussian modular charges are computed.
    constant = result["grid_controls"][-2]["profiles"][0]
    assert constant["conformal_continuum_stress_integral"] == pytest.approx(4*math.pi**2/15)


@pytest.mark.parametrize("bad", [True, 1., 1+0j, F(0), F(-1), F(1, 10**400), F(10**100)])
def test_original_exact_lengths_and_precision_domain(bad):
    with pytest.raises(ValueError):
        ball.covariance(bad, F(8))


@pytest.mark.parametrize("args", [(F(2, 5), F(6), F(1)), (F(1, 10), F(8), F(1)),
                                  (F(1), F(8), F(1)), (F(2, 5), F(2), F(1))])
def test_grid_domain_and_region_shape_refusals(args):
    with pytest.raises(ValueError):
        ball.ball_sites(*args)


@pytest.mark.parametrize("bad", [np.array([[.7+0j]]), np.ma.array([[.7]], mask=[[False]]),
                                  np.array([[float("nan")]]), [[-.7]], [[.7, .1]], [[1e300]]])
def test_covariance_domain_refusals(bad):
    with pytest.raises(ValueError):
        ball.modular_momentum(bad, [[.7]])


@pytest.mark.parametrize("bad", [[1+0j], [float("inf")], [[1.]], np.ma.array([1.])])
def test_direction_domain_refusals(bad):
    with pytest.raises(ValueError):
        ball.charge([[.7]], [[.7]], bad)


def test_nonfinite_requested_charge_refused():
    with pytest.raises(ValueError, match="charge overflowed"):
        ball.charge([[.7]], [[.7]], [1e300])
    assert ball.charge([[.7]], [[.7]], [0.]) == 0


def test_tiny_covariance_asymmetry_not_hidden_by_absolute_tolerance():
    with pytest.raises(ValueError, match="symmetric"):
        ball.modular_momentum([[1e-40, 1e-41], [0., 1e-40]], np.eye(2)*1e40)


def test_source_bindings_and_scope(result):
    for name, digest in result["source_bindings"].items():
        assert hashlib.sha256((ball.ROOT/name).read_bytes()).hexdigest() == digest
    assert result["scope"]["finite_Gaussian_partial_trace_and_first_law"]
    assert not result["scope"]["native_smooth_tangent_modular_identification"]
    assert not result["scope"]["finite_regulator_rotation_invariant"]
    assert not result["scope"]["renormalized_generalized_entropy_or_area_density"]
    assert not result["scope"]["joint_state_geometry_stationarity_or_Einstein_dynamics"]
    assert not result["scope"]["massive_full_matter_or_calibrated_particle_spectrum"]


def test_receipt_roundtrip_and_corruption_controls(result, tmp_path, monkeypatch):
    receipt = tmp_path/"ball.json"
    monkeypatch.setattr(ball, "run", lambda: result)
    monkeypatch.setattr("sys.argv", ["ball_modular", "--write", str(receipt)])
    ball.main()
    loaded = json.loads(receipt.read_text())
    assert isinstance(loaded, dict)
    ball.compare_ball_receipt(result, loaded)
    monkeypatch.setattr("sys.argv", ["ball_modular", "--check", str(receipt)])
    ball.main()
    corrupt = deepcopy(loaded)
    corrupt["grid_controls"][0]["regional_sites"] += 1
    with pytest.raises(ValueError):
        ball.compare_ball_receipt(result, corrupt)
    corrupt = deepcopy(loaded)
    corrupt["grid_controls"][0]["profiles"][0]["modular_entropy_derivative"] += .001
    with pytest.raises(ValueError):
        ball.compare_ball_receipt(result, corrupt)
    corrupt = deepcopy(loaded)
    corrupt["source_bindings"][next(iter(corrupt["source_bindings"]))] = "0"*64
    with pytest.raises(ValueError):
        ball.compare_ball_receipt(result, corrupt)
    corrupt = deepcopy(loaded)
    corrupt["entropy_secant_replay_absolute_tolerance"] = 0.
    with pytest.raises(ValueError):
        ball.compare_ball_receipt(result, corrupt)
    corrupt = deepcopy(loaded)
    corrupt["grid_controls"][4]["minimum_regional_symplectic_gap_squared"] = -1e-7
    with pytest.raises(ValueError, match="purity/precision"):
        ball.compare_ball_receipt(result, corrupt)
    corrupt = deepcopy(loaded)
    corrupt["grid_controls"][4]["modular_eigen_residual_norm"] = 1e-11
    with pytest.raises(ValueError, match="purity/precision"):
        ball.compare_ball_receipt(result, corrupt)

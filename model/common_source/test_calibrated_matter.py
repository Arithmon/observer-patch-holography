"""Independent spectrum, current and complex-convention calibration checks."""
from copy import deepcopy
import json

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import calibrated_matter as c
from .chiral_candidate import action, kinetic


def numeric(matrix):
    return np.asarray(matrix.evalf(20) if hasattr(matrix, "evalf") else matrix, complex)


@pytest.fixture(scope="module")
def report():
    return c.report()


def test_one_common_scale_is_calibrated_and_not_pole_or_held_out_validation(report):
    packet = report["calibration"]
    assert packet["common_scale"]["scheme"] == "MSbar Standard Model"
    assert 248 < float(packet["matched_SM"]["v_SM_GeV"]) < 249
    assert packet["calibration_contract"]["all_nine_Yukawas_CKM_and_v_are_empirical_inputs"]
    assert not packet["calibration_contract"]["same_input_mass_recovery_is_held_out_validation"]
    assert not report["scope"]["quark_running_parameters_are_free_asymptotic_poles"]
    assert not report["scope"]["physical_parameters_derived_from_native_repair"]
    b = c.bosonic_parameters()
    assert b["gY"]**2 == s.Rational(3, 5)*b["g1_GUT"]**2
    assert b["lambda_SM"] == b["lambda_paper"]/4
    assert b["neutral_higgs_norm_GeV"]**2 == b["v_SM_GeV"]**2/2


def test_dense_complete_45_singular_values_recover_all_fitted_channels_and_color(report):
    k = numeric(c.mass_kernel())
    found = np.sort(np.linalg.svd(k, compute_uv=False))
    expected = [0.]*3
    for name in c.CHANNELS:
        for mass in c.mass_coordinates(name):
            expected.extend([float(mass)]*(2 if name == "charged_lepton" else 6))
    np.testing.assert_allclose(found, np.sort(expected), rtol=4e-11, atol=4e-14)
    assert np.linalg.matrix_rank(k, tol=1e-10) == 42
    masses = [float(x) for x in c.mass_coordinates("charged_lepton")]
    assert masses[1]/masses[0] > 200 and masses[2]/masses[0] > 3000
    assert report["identities"]["source_current_Higgs_identities_exact"] == 72


def test_physical_charged_current_and_cp_orientation_follow_the_same_mass_basis():
    w, k = numeric(c.mass_basis_change()), numeric(c.mass_kernel())
    np.testing.assert_allclose(w.conj().T@w, np.eye(45), atol=4e-16, rtol=0)
    d = w.T@k@w
    for name in c.CHANNELS:
        for (i, j), sign in c._oriented_components(name)[1].todok().items():
            np.testing.assert_allclose(d[3*i:3*i+3, 3*j:3*j+3], np.diag(list(map(float, c.mass_coordinates(name)))), atol=4e-14, rtol=0)
    plus = numeric(c.weak_raising_current())
    transformed = w.conj().T@plus@w
    for up, down in ((2, 3), (5, 6), (7, 8)):
        np.testing.assert_allclose(transformed[3*up:3*up+3, 3*down:3*down+3], numeric(c.ckm_matrix()), atol=1e-16, rtol=0)
    v = numeric(c.ckm_matrix())
    j = np.imag(v[0, 0]*v[1, 1]*np.conj(v[0, 1]*v[1, 0]))
    assert 3.0e-5 < j < 3.3e-5
    assert np.linalg.norm(v-v.conj()) > .005
    # Conjugating only the family coefficient reverses the weak-current CP
    # convention; the paired transpose/adjoint distinction is observable.
    wrong = w.conj().T@k@w.conj()
    assert np.max(abs(wrong-d)) > .01


def test_actual_finite_weak_basis_transform_preserves_a_complex_higgs_vertex():
    name, component = "down", 0
    coefficients = c.higgs_components()[name]
    h = np.array([.17+.11j, -.23+.07j])
    k = sum((x*numeric(g) for x, g in zip(h, coefficients)), np.zeros((45, 45), complex))
    primitive = action._representation()[1][10]
    current = action.generators()[10]
    higgs_generator = primitive[3:, 3:]
    if action._representation()[4][name]["conjugate_Higgs"]:
        higgs_generator = -higgs_generator.T
    uf = expm(.19*numeric(current)); uh = expm(.19*numeric(higgs_generator))
    rotated = uh@h
    target = sum((x*numeric(g) for x, g in zip(rotated, coefficients)), np.zeros((45, 45), complex))
    np.testing.assert_allclose(uf.T@target@uf, k, atol=2e-18, rtol=1e-14)


def test_full_spin_CAR_electric_and_Weyl_mass_shell_share_the_calibrated_kernel():
    k = c.mass_kernel(); delta = c.pairing_kernel(); q = action.electric_charge()
    assert k == k.T and delta == -delta.T
    assert q.T*k+k*q == s.zeros(45)
    np.testing.assert_allclose(numeric(delta), numeric(kinetic.hamiltonian_pairing(k)), atol=1e-16, rtol=1e-15)
    momentum = (.11, -.07, .19)
    speed = 1.
    symbol = c.weyl_symbol(momentum)
    eigen = np.linalg.eigvalsh(symbol)
    masses = np.linalg.svd(numeric(k), compute_uv=False)
    expected_positive = np.repeat(np.sqrt(np.dot(momentum, momentum)+masses*masses), 2)
    np.testing.assert_allclose(eigen, np.sort(np.r_[-expected_positive, expected_positive]), atol=2e-12, rtol=3e-13)
    charge = kinetic.nambu_charge(q)
    np.testing.assert_allclose(symbol@charge, charge@symbol, atol=0, rtol=0)


def test_higgs_off_and_calibrated_scale_controls_are_not_new_mass_predictions():
    assert c.mass_kernel(0) == s.zeros(45)
    assert c.pairing_kernel(0) == s.zeros(90)
    assert c.mass_kernel(2*c.neutral_higgs_norm()) == 2*c.mass_kernel()
    for value in (False, .1, s.Float(.1), -1, s.Symbol("v")):
        with pytest.raises(ValueError):
            c.mass_kernel(value)
    for component in (False, -1, 2, 1.):
        with pytest.raises(ValueError):
            c.mass_kernel(component=component)


def test_declared_bosonic_tree_comparison_exposes_its_extra_m2_prescription():
    result = c.bosonic_tree_control()
    assert result["Goldstone_zero_directions_exact"] == 3
    assert result["photon_zero_direction_exact"] == 1
    assert result["unbroken_color_gluon_zero_directions_exact"] == 8
    b = c.bosonic_parameters()
    g, y, v = b["gW"], b["gY"], b["v_SM_GeV"]
    neutral = s.Matrix([[g*g, -g*y], [-g*y, y*y]])*v*v/4
    assert neutral*s.Matrix([y, g]) == s.zeros(2, 1)
    assert s.simplify(s.trace(neutral)-(g*g+y*y)*v*v/4) == 0
    assert not result["physical_pole_masses_or_widths_are_these_tree_coefficients"]
    assert result["calibrated_v_is_SMDR_loop_minimum_VEV"]
    assert result["tree_stationarity_is_an_additional_comparison_prescription"]
    assert not result["tree_comparison_m2_is_SMDR_loop_matched_m2"]
    assert not result["tree_Goldstone_zeros_are_claimed_for_full_SMDR_loop_action_Hessian"]
    assert s.sympify(result["tree_comparison_m2_GeV2_exact"]) == -b["lambda_SM"]*v*v


def test_optional_same_Higgs_Weinberg_extension_has_calibrated_PMNS_and_splittings():
    report = c.neutrino_report()
    packet = report["calibration"]
    assert packet["source"]["selection"].startswith("normal ordering")
    assert not packet["matching_contract"]["renormalization_group_transport_to_common_MZ_boundary_executed"]
    u = numeric(c.pmns_matrix()); m = np.array(list(map(float, c.neutrino_mass_coordinates())))
    np.testing.assert_allclose(u.conj().T@u, np.eye(3), atol=4e-16, rtol=0)
    np.testing.assert_allclose(u.T@numeric(c.neutrino_family_kernel())@u, np.diag(m), atol=1e-26, rtol=4e-15)
    assert m[0] == 0
    assert m[1]**2*1e18 == pytest.approx(.00007537, abs=2e-20)
    assert m[2]**2*1e18 == pytest.approx(.002511, abs=2e-18)
    # Majorana operator is already in the actual neutrino slot of M15; no
    # right neutrino is appended, and its Higgs law is quadratic.
    a = c.mass_kernel(include_neutrinos=True)-c.mass_kernel()
    assert a.shape == (45, 45)
    assert {i//3 for i, j in a.todok()} == {10}
    scaled = c.mass_kernel(2*c.neutral_higgs_norm(), include_neutrinos=True)-c.mass_kernel(2*c.neutral_higgs_norm())
    assert c._zero(scaled-4*a)
    w = numeric(c.mass_basis_change(include_neutrinos=True))
    current = w.conj().T@numeric(c.weak_raising_current())@w
    np.testing.assert_allclose(current[30:33, 33:36], u.conj().T, atol=1e-16, rtol=0)
    # Using U rather than U dagger here would reverse the all-left basis
    # convention and generally changes the neutral mass response.
    assert np.max(abs(u.T@numeric(c.neutrino_family_kernel()).conj()@u-np.diag(m))) > 1e-12


def test_public_metadata_is_detached_and_exact_types_are_checked(report):
    packet = c.parameters(); packet["matched_SM"]["yukawa"]["up"][0] = "9"
    assert 6e-6 < float(c.parameters()["matched_SM"]["yukawa"]["up"][0]) < 8e-6
    forged = deepcopy(report); forged["identities"]["charged_K45_rank_exact"] = True
    with pytest.raises(ValueError, match="type"):
        c._compare(forged, report)
    forged = deepcopy(report); forged["calibration"]["CKM"]["delta_radians"] = "-1.147"
    with pytest.raises(ValueError, match="changed"):
        c._compare(forged, report)
    assert c.source_bindings()["model/common_source/calibrated_matter_parameters.json"]


def test_retained_action_replays_exactly():
    path = c.ROOT/"model/common_source/calibrated_matter_results.json"
    retained = json.loads(path.read_text())
    c._compare(retained, json.loads(json.dumps(c.report(), allow_nan=False)))

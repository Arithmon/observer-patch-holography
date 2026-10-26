"""Independent algebra, source, finite-cut and refusal controls for the wedge."""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
import sys

import numpy as np
import pytest
from scipy.linalg import eigh_tridiagonal
from scipy.special import ellipk, eval_legendre
import sympy as s

from . import wedge_modular as w


@lru_cache(maxsize=1)
def prepared():
    return w.source_preparation()


def test_exact_native_parameters_and_transverse_representation():
    from ..common_source.native_staggered_field import clifford_matrices
    x, z = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    assert clifford_matrices() == (s.kronecker_product(z, s.eye(4)), s.kronecker_product(x, z, s.eye(2)), s.kronecker_product(x, x, z), s.kronecker_product(x, x, x))
    assert tuple(w.parameters(j)[0] for j in range(3)) == (F(1, 256), F(1, 1024), F(1, 3072))
    assert w.parameters()[1] == F(125, 512)
    assert w.transverse_algebra()["scalar_chains_per_reduced_transverse_momentum"] == 4
    assert w.transverse_algebra()["Dirac_tastes"] == 2


def test_lapse_stress_fixed_before_modular_prefactor():
    m, c, a = F(3, 8), F(5, 4), F(2, 3)
    h = w.chain_matrix(8, m, c, a)
    x = np.diag(float(a)*(np.arange(1, 9)-.5))
    b = w.boost_matrix(8, m, c, a)
    np.testing.assert_allclose(b, (x@h+h@x)/2, atol=1e-15, rtol=0)
    np.testing.assert_allclose(np.diag(b, 1), -float(c/2)*np.arange(1, 8), atol=1e-15, rtol=0)
    mu = m*a/c
    t = w.boost_matrix(8, mu, F(1), F(1))
    np.testing.assert_allclose(b, float(c)*t, atol=1e-15, rtol=0)


def test_onsite_gauge_preserves_cut_and_oriented_imaginary_hopping():
    h = w.chain_matrix(8, *w.parameters())
    gauge = (-1j)**np.arange(1, 9)
    out = gauge[:, None]*h*gauge.conj()[None, :]
    np.testing.assert_allclose(np.diag(out, 1), -1j*float(w.parameters()[1]/2), atol=0, rtol=0)
    assert np.count_nonzero(np.diag(gauge)[:4, 4:]) == 0


@pytest.mark.parametrize("mu2", [F(1, 562500), F(1, 62500), F(4, 15625), F(1, 16), F(1), F(16)])
def test_exact_elliptic_intervals_and_continuum_prefactor(mu2):
    data = w.elliptic_data(mu2, F(125, 512))
    kappa = 1/math.sqrt(1+float(mu2))
    gamma = 4*kappa*ellipk(float(mu2/(1+mu2)))/float(F(125, 512))
    delta = 2*math.pi*ellipk(float(mu2/(1+mu2)))/ellipk(float(1/(1+mu2)))
    assert math.isclose(float(sum(data["gamma_interval"])/2), gamma, rel_tol=4e-15)
    assert math.isclose(float(sum(data["entanglement_spacing_interval"])/2), delta, rel_tol=3e-11)
    low, high = data["continuum_prefactor_ratio_interval"]
    assert high >= 1-mu2/4
    assert low <= 1
    # Enclosure of nontrivial subleading coefficient, not only mu²/4 bound.
    assert abs(float((1-(low+high)/2)/mu2)-.25) < .03 if mu2 <= F(1, 16) else True


def test_band_bound_is_uniform_and_requires_declared_band():
    m, c = w.parameters()
    for a in (F(1, 4), F(1, 8), F(1, 16)):
        bound = w.band_bound(m, c, a, F(1))
        for py, pz in ((0., 0.), (.2, -.6), (1., -1.)):
            mu2 = (float(m*a/c))**2+math.sin(float(a)*py)**2+math.sin(float(a)*pz)**2
            ratio = w.elliptic_data(F(mu2))["continuum_prefactor_ratio_interval"]
            assert float(1-ratio[0]) <= float(bound)+1e-15
    with pytest.raises(ValueError):
        w.band_bound(m, c, F(1), F(1))


def test_legendre_generating_function_normalization_and_tail():
    mass, limit, degree = F(1, 4), F(2), 32
    r = float((limit-mass)/(limit+mass))
    values = np.linspace(float(mass*mass), float(limit*limit), 101)
    z = (float(limit*limit+mass*mass)-2*values)/float(limit*limit-mass*mass)
    polynomial = 2/float(limit+mass)*sum(r**j*eval_legendre(j, z) for j in range(degree+1))
    inverse_tail = r**(degree+1)/float(mass)
    assert np.max(abs(1/np.sqrt(values)-polynomial)) <= inverse_tail*(1+1e-12)
    # Missing factor2 in generating function would fail even at the gap.
    assert abs(polynomial[0]-1/float(mass)) < inverse_tail
    assert w.projector_bound(mass, mass, 0) == 0
    with pytest.raises(ValueError):
        w.projector_bound(F(2), F(1), 8)


def test_finite_walk_agreement_and_projector_bound():
    # Independent finite matrices with identical paths from one central site.
    m, c = F(1, 2), F(1)
    small, large = w.chain_matrix(16, m, c), w.chain_matrix(32, m, c)
    for degree in range(1, 8):
        assert math.isclose(np.linalg.matrix_power(small, degree)[8, 8], np.linalg.matrix_power(large, degree)[16, 16], abs_tol=1e-15)
    def projector(h):
        vals, vecs = np.linalg.eigh(h)
        return vecs[:, vals < 0]@vecs[:, vals < 0].T
    actual = abs(projector(small)[8, 8]-projector(large)[16, 16])
    bound = w.projector_bound(m, F(6, 5), 3)
    assert actual < float(bound)
    # Equal window entries alone do not determine the vacuum covariance.
    changed = small.copy()
    changed[7, 8] = changed[8, 7] = 0
    assert changed[8, 8] == small[8, 8]
    assert abs(projector(changed)[8, 8]-projector(small)[8, 8]) > .01


def test_analytic_sine_mass_frame_matches_independent_diagonalization():
    m, c = w.parameters()
    occupied = w.orbitals(8, m, c)
    positive = w.orbitals(8, m, c, positive=True)
    frame = np.column_stack((occupied, positive))
    h = w.chain_matrix(16, m, c)
    # Accumulated sine-phase/matrix arithmetic scales with the finite
    # dimension; these diagnostics do not replace the exact sine identity.
    arithmetic = 64*len(frame)*np.finfo(float).eps
    assert np.linalg.norm(frame.T@frame-np.eye(16), 2) < arithmetic
    vals, vecs = eigh_tridiagonal(np.diag(h), np.diag(h, 1))
    assert np.linalg.norm(occupied@occupied.T-vecs[:, vals < 0]@vecs[:, vals < 0].T, 2) < arithmetic
    assert np.linalg.eigvalsh(occupied.T@h@occupied).max() < 0


@pytest.mark.parametrize("sector", [0, 1, 2])
def test_infinite_kernel_independent_of_BW_and_certified_local_IR_control(sector):
    m, c = w.parameters(sector)
    infinite, quadrature = w.infinite_window(m, c)
    half = 32768
    rows = list(range(half+1, half+5))
    frame = w.orbitals(half, m, c, rows=rows)
    error = np.linalg.norm(frame@frame.T-infinite, 2)
    bound = w.chain_bound(half, 4, m, c)
    assert bound["maximum_walk_length"] <= bound["distance_to_omitted_boundary"]
    assert error < float(bound["local_covariance_operator_norm_upper"]+quadrature)+w.NUMERIC_GATE
    assert float(bound["local_covariance_operator_norm_upper"]) < 1e-10


def test_exact_staggered_level_offset_and_single_multiplicity():
    absolute = sorted(abs(2*j-F(1, 2)) for j in range(-10, 11))[:12]
    assert absolute == [j+F(1, 2) for j in range(12)]
    control = w.entropy_series(F(1, 16))
    finite = w.finite_entropy(32, F(1, 4), F(1))
    assert abs(finite-control["entropy_diagnostic"]) < 1e-8
    assert float(control["analytic_series_tail_upper"]) < 1e-30
    # Wrong half-level/multiplicity inherited from Ising would be material.
    assert abs(2*finite-control["entropy_diagnostic"]) > .4


def test_native_finite_entropy_retains_all_IR_and_nonmonotonic_controls():
    m, c = w.parameters()
    exact = w.entropy_series((m/c)**2)["entropy_diagnostic"]
    finite = [w.finite_entropy(n, m, c) for n in (8, 32, 128, 512)]
    assert np.max(np.diff(abs(np.array(finite)-exact))) < 0
    assert abs(finite[-1]-exact) < 1e-8
    m, c = w.parameters(1)
    assert w.finite_entropy(512, m, c) > w.entropy_series((m/c)**2)["entropy_diagnostic"]


def test_bare_area_and_taste_count_are_independently_geometric():
    result = w.planar_entropy(0, 8, F(1, 2))
    assert result["proper_area_exact"] == 16
    assert result["Dirac_tastes"] == 2
    assert result["scalar_chains_per_reduced_transverse_momentum"] == 4
    assert not result["renormalized_area_coefficient_or_Newton_constant"]
    assert result["bare_entropy_per_proper_area_diagnostic"] == result["entropy_diagnostic"]/16
    for count in (6, 10):
        with pytest.raises(ValueError):
            w.planar_entropy(0, count)


def test_independent3D_projector_confirms_four_chain_taste_factor():
    result = w.finite_planar_control()
    assert result["full_one_body_modes"] == 128
    assert result["scalar_chains"] == 16
    assert abs(result["direct_native_projector_entropy_diagnostic"]-result["four_chain_reduction_entropy_diagnostic"]) < 1e-12
    assert not result["large_Fock_source_preparation_executed"]


def test_source_literal_vacuum_and_exact_global_energy_certificate():
    state, target, resources = prepared()
    np.testing.assert_allclose(state, target, atol=1e-11, rtol=0)
    certificate = w.certify_source(state)
    assert certificate["global_raw_half_trace_norm_upper"] < F(1, 10**10)
    assert certificate["gap_lower"] > F(1, 100)
    assert resources["proper_code_service_events"] == 2*resources["native_pulses"]
    assert resources["native_pulses"] == 796
    h = w.chain_matrix(8, *w.parameters())
    value = np.vdot(state, w._fock(h)@state).real/np.vdot(state, state).real
    assert abs(float(certificate["Rayleigh_energy_exact"])-value) < 1e-14
    roots = np.linalg.eigvalsh(h)
    assert abs(sum(roots[:4])-float(sum(certificate["ground_energy_interval"])/2)) < 1e-14


def test_source_retained_arrays_are_not_normalized_or_phase_repaired():
    state = prepared()[0]
    for changed in (state*(1+1e-4), np.eye(1, 256, 0, dtype=complex)[0], np.roll(state, 1)):
        with pytest.raises(ValueError):
            w.certify_source(changed)
    # Energy proof is phase independent, but the literal source target gate
    # is not: receipt checker also compares the original complex amplitudes.
    assert w.certify_source(1j*state)["global_raw_half_trace_norm_upper"] < w.TRACE_GATE
    assert np.linalg.norm(1j*state-prepared()[1]) > 1
    with pytest.raises(ValueError):
        w.certify_source(np.ma.array(state, mask=False))
    bad = state.copy()
    bad[0] = np.nan
    with pytest.raises(ValueError):
        w.certify_source(bad)


def test_full_Fock_partial_trace_nonGaussian_mixture_first_law():
    result = w.source_first_law(prepared()[0])
    assert abs(result["full_Fock_first_law_derivative"]-result["quadratic_first_law_derivative"]) < 1e-9
    errors = abs(np.asarray(result["mixture_secants"])-result["full_Fock_first_law_derivative"])
    assert np.all(np.diff(errors) < 0)
    assert abs(result["full_Fock_first_law_derivative"]-result["infinite_wedge_weight_finite_chain_response"]) > .001
    assert not result["finite_interval_BW_equality_claimed"]


def test_local_phase_stress_change_is_fixed_by_only_first_bond():
    m, c = w.parameters()
    infinite, _ = w.infinite_window(m, c)
    u = np.diag([np.exp(.25j), 1, 1, 1])
    dc = u@infinite@u.conj().T-infinite
    direct = float(np.trace(w.boost_matrix(4, m, c)@dc).real)
    local = float(c)*(1-math.cos(.25))*infinite[0, 1]
    assert local > 0
    assert abs(direct-local) < 1e-16


def test_finite_modular_refuses_unresolved_purity_and_bad_domains():
    for value in (np.eye(4), np.zeros((4, 4)), np.diag([1e-16, .5]), np.diag([-.01, .5]), np.ma.array(np.eye(4)/2, mask=False), np.array([[.5, 1j], [0, .5]])):
        with pytest.raises(ValueError):
            w.finite_modular(value)
    m, c = w.parameters()
    small = w.orbitals(4, m, c)[4:]
    assert np.isfinite(w.finite_modular(small@small.T)).all()
    large = w.orbitals(32, m, c)[32:]
    with pytest.raises(ValueError):
        w.finite_modular(large@large.T)


@pytest.mark.parametrize("bad", [True, .1, 0, F(0), F(-1)])
def test_exact_positive_domain_refusals(bad):
    with pytest.raises(ValueError):
        w.elliptic_data(bad)


def test_numerical_profiles_refuse_overflow_underflow_before_evaluation():
    for mass in (F(10)**308, F(1, 10**400)):
        with pytest.raises(ValueError, match="numerical chain scales"):
            w.boost_matrix(8, mass, F(1))
        with pytest.raises(ValueError, match="numerical chain scales"):
            w.orbitals(4, mass, F(1))


def test_receipt_exact_gates_scope_and_numeric_domains():
    fresh = {"admission": w.SOURCE_GATE, "diagnostic": .5, "scope": False, "tail": "1/1000"}
    for key, value in (("admission", 0.), ("scope", True), ("tail", "0"), ("diagnostic", float("nan"))):
        old = deepcopy(fresh)
        old[key] = value
        with pytest.raises(ValueError):
            w.compare(old, fresh)
    w.compare({**fresh, "diagnostic": .5+1e-10}, fresh)


def test_public_exact_data_return_detached_not_mutable_cached_claims():
    mu2 = F(4, 15625)
    first = w.elliptic_data(mu2)
    first["relative_prefactor_error_upper"] = F(0)
    assert w.elliptic_data(mu2)["relative_prefactor_error_upper"] == mu2/4


def _receipt(tmp_path, monkeypatch):
    path = tmp_path/"wedge.json"
    monkeypatch.setattr(sys, "argv", ["wedge", "--output", str(path)])
    w.main()
    monkeypatch.setattr(sys, "argv", ["wedge", "--check", str(path)])
    return path


def test_cli_roundtrip_and_forged_exact_binary_certificate(tmp_path, monkeypatch):
    path = _receipt(tmp_path, monkeypatch)
    w.main()
    report = json.loads(path.read_text())
    report["source_preparation"]["own_binary_certificate"]["infidelity_upper"] = "0"
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="exact replay"):
        w.main()


@pytest.mark.parametrize("change, admitted", [(1e-10, True), (1e-7, False)])
def test_cli_own_binary_gates_allow_roundoff_but_reject_phase_drift(tmp_path, monkeypatch, change, admitted):
    path = _receipt(tmp_path, monkeypatch)
    artifact = path.with_suffix(".npz")
    with np.load(artifact, allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in saved.files}
    arrays["source_state"] = np.exp(1j*change)*arrays["source_state"]
    arrays["source_covariance"] = w._covariance(arrays["source_state"])
    report = json.loads(path.read_text())
    report["source_preparation"]["own_binary_certificate"] = w._json(w.certify_source(arrays["source_state"]))
    with artifact.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    report["artifact"]["sha256"] = hashlib.sha256(artifact.read_bytes()).hexdigest()
    path.write_text(json.dumps(report))
    if admitted:
        w.main()
    else:
        with pytest.raises(ValueError, match="retained binary replay|independent target"):
            w.main()

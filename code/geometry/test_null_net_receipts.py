#!/usr/bin/env python3
"""Tests for the null-net receipt instrumentation (#503, #524)."""

from __future__ import annotations

import sys
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from null_net_receipts import (  # noqa: E402
    hsm_compression_receipt,
    instrument_null_net,
    lie_closure_receipt,
    mixed_gns_cauchy_receipt,
    nti_receipt,
    separating_modulus_receipt,
    weak_additivity_receipt,
    correlation_matrix,
    bond_limit_error_bound,
    embed,
    modular_subspace_diagnostic,
    momentum_profile,
    _shape_fit,
)


def test_nti_and_weak_additivity_are_witnessed():
    for n in (16, 32):
        assert nti_receipt(n)["nontrivial"]
        assert weak_additivity_receipt(n)["covers_ring"]


def test_separating_modulus_positive_at_every_stage():
    r = separating_modulus_receipt(32)
    assert r["strictly_faithful"]
    # the gap is astronomically small but strictly positive: this is why
    # the instrumentation needs extended precision
    assert r["occupation_gap_log10"] < -10


def test_mixed_gns_cauchy():
    r = mixed_gns_cauchy_receipt((16, 32, 64))
    assert r["cauchy_decreasing"]


def test_hsm_compression_has_consistent_sign_asymmetry():
    r32 = hsm_compression_receipt(32)
    r64 = hsm_compression_receipt(64)
    assert r32["asymmetry_ratio"] > 3.0
    assert r64["asymmetry_ratio"] > 3.0
    assert r32["lower_packet_leakage_sign"] == r64["lower_packet_leakage_sign"]
    for r in (r32, r64):
        assert r["max_subspace_leakage_plus"] > 0.98
        assert r["max_subspace_leakage_plus"] == pytest.approx(
            r["max_subspace_leakage_minus"], abs=1e-12)
        assert not r["hsm_inclusion_certified"]


def test_hsm_real_packet_is_the_null_control():
    # a real (non-chiral) packet is exactly time-reversal symmetric: the
    # asymmetry must vanish. Chirality selects directional packet transport,
    # not half-sided inclusion of the entire subspace.
    import numpy as np
    from scipy.linalg import expm
    from modular_clock_instrumentation import arc_entanglement_hamiltonian
    from null_net_receipts import embed
    n, m_a, m_b = 32, 16, 8
    h_a = arc_entanglement_hamiltonian(n, m_a)
    h_b = embed(arc_entanglement_hamiltonian(n, m_b), m_a)
    xs = np.arange(m_a)
    psi = np.exp(-((xs - m_b / 2.0) ** 2) / (2 * (m_b / 6.0) ** 2)).astype(complex)
    psi[m_b:] = 0.0
    psi /= np.linalg.norm(psi)

    def leakage(t):
        out = expm(1j * h_a * t) @ expm(-1j * h_b * t) @ psi
        return float(np.sum(np.abs(out[m_b:]) ** 2))

    assert abs(leakage(0.12) - leakage(-0.12)) < 1e-12


def test_lie_closure_percent_level_and_resummation_control():
    r = lie_closure_receipt(32)
    assert r["relative_residual"] < 0.02
    # the unresummed range-2 truncation is several times worse: the
    # resummation is load-bearing
    assert (r["relative_residual_unresummed_control"]
            > 5.0 * r["relative_residual"])


def test_instrumented_null_net_verdicts():
    report = instrument_null_net(rings=(16, 32))
    w = report["receipts_witnessed"]
    assert w["nti"] and w["weak_additivity"]
    assert w["separating_faithfulness"] and w["single_bond_observable_limit"]
    assert not w["mixed_gns_cauchy"]
    assert not w["hsm_compression_one_particle"]
    assert not w["modular_lie_closure_percent_level"]
    assert report["verdicts"]["lie_closure_percent_level"]
    # the rate clause must be recorded as open in the report
    assert any("convergence rate" in p or "rate" in p
               for p in report["receipts_pending"])


def test_full_ring_covariance_is_the_independent_fourier_projector():
    n = 16
    momenta = (2 * np.arange(n) + 1) * np.pi / n
    occupied = momenta[np.cos(momenta) > 0]
    modes = np.exp(1j * np.outer(np.arange(n), occupied)) / np.sqrt(n)
    expected = modes @ modes.conj().T
    actual = correlation_matrix(n)
    np.testing.assert_allclose(actual, expected, atol=3e-15)
    np.testing.assert_allclose(actual @ actual, actual, atol=3e-15)


def test_packet_asymmetry_cannot_witness_half_sided_inclusion():
    from scipy.linalg import expm
    from modular_clock_instrumentation import arc_entanglement_hamiltonian
    h = arc_entanglement_hamiltonian(16, 8)
    leakage = [np.linalg.norm(expm(1j * t * h)[4:, :4], 2) ** 2
               for t in (0.12, -0.12)]
    assert min(leakage) > 0.65
    assert leakage[0] == pytest.approx(leakage[1], abs=1e-14)
    assert not instrument_null_net((16,))["receipts_witnessed"][
        "hsm_compression_one_particle"]


def test_single_stage_cannot_supply_mixed_gns_or_an_empty_closure_witness():
    report = instrument_null_net((16,))
    assert not report["verdicts"]["lie_closure_percent_level"]
    assert not report["receipts_witnessed"]["mixed_gns_cauchy"]


def test_lie_fit_rejects_a_cutoff_that_erases_the_signal():
    with pytest.raises(ValueError, match="cutoff|rmax"):
        lie_closure_receipt(32, rmax=0)


@pytest.mark.parametrize("n", [8, 12, 16, 24, 32])
def test_covariance_spectrum_seam_and_arc_faithfulness(n):
    # Separate discrete Fourier implementation, including unoccupied modes.
    k = (2 * np.arange(n) + 1) * np.pi / n
    basis = np.exp(1j * np.outer(np.arange(n), k)) / np.sqrt(n)
    occupied, empty = basis[:, np.cos(k) > 0], basis[:, np.cos(k) < 0]
    c = correlation_matrix(n)
    np.testing.assert_allclose(c, occupied @ occupied.conj().T, atol=6e-15)
    np.testing.assert_allclose(np.eye(n) - c, empty @ empty.conj().T, atol=6e-15)
    np.testing.assert_allclose(c @ c, c, atol=2e-15)
    assert np.trace(c) == n / 2
    assert c[0, -1] == pytest.approx(-c[0, 1], abs=1e-16)
    for m in (1, n // 4, n // 2):
        # Gram rank checks the Vandermonde argument without subtracting tiny
        # eigenvalues from one in binary64 (which fails by N=64).
        assert np.linalg.matrix_rank(occupied[:m]) == m
        assert np.linalg.matrix_rank(empty[:m]) == m


def test_precision_context_is_local():
    import mpmath as mp
    c = correlation_matrix(16)
    expected = separating_modulus_receipt(32)
    with mp.workdps(15):
        np.testing.assert_array_equal(correlation_matrix(16), c)
        assert separating_modulus_receipt(32) == expected
        assert mp.mp.dps == 15


def test_bond_limit_against_independent_quadrature_and_fourier_sum():
    from scipy.integrate import quad
    from math import exp, pi
    rings = (8, 12, 32, 64, 128)
    report = mixed_gns_cauchy_receipt(rings)
    f = lambda x: exp(-200 * (x - 0.25) ** 2)
    limit, quadrature_error = quad(f, 0, 1, epsabs=1e-13, epsrel=1e-13)
    limit *= 2 / pi
    assert quadrature_error < 1e-12
    assert report["limit_value"] == pytest.approx(limit, abs=1e-15)
    for n, value, bound in zip(rings, report["values"], report["limit_error_bounds_exact"]):
        momenta = (2 * np.arange(n) + 1) * pi / n
        bond = np.sum(np.exp(1j * momenta[np.cos(momenta) > 0])).real / n
        independent = 2 * bond * sum(f(j / n) for j in range(n - 1)) / n
        assert value == pytest.approx(independent, abs=1e-15)
        assert abs(value - limit) < Fraction(bound)
    for difference, bound in zip(report["cauchy_differences"], report["pairwise_cauchy_bounds_exact"]):
        assert difference < Fraction(bound)
    # The bound supplies a cofinal rate for arbitrary large stages, without
    # evaluating larger matrices or extrapolating from a finite fit.
    for n in (8, 32, 128, 10**6, 10**20):
        assert 0 < bond_limit_error_bound(n) < Fraction(3, n)
        assert bond_limit_error_bound(2 * n) < bond_limit_error_bound(n)
    assert not report["mixed_gns_certified"]


def test_exact_cyclic_counterexample_to_packet_inclusion():
    import sympy as sp
    # U e0=e1 stays inside B, but U e1=e2 escapes. Reverse flow sends
    # e0 outside. Thus a perfect packet asymmetry coexists with maximal
    # leakage of the entire subspace in BOTH directions.
    u = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    p = sp.diag(1, 1, 0)
    q = sp.eye(3) - p
    e0 = sp.eye(3)[:, 0]
    assert q * u * e0 == sp.zeros(3, 1)
    assert (q * u.T * e0).norm() == 1
    for flow in (u, u.T):
        block = q * flow * p
        assert max((block.T * block).eigenvals()) == 1
    # A Hermitian logarithm realizes this unitary as a finite continuous flow.
    eigenvalues, eigenvectors = np.linalg.eig(np.array(u).astype(complex))
    h = (eigenvectors * np.angle(eigenvalues)) @ eigenvectors.conj().T
    r = modular_subspace_diagnostic(h, np.zeros((2, 2)), 1, np.array([1, 0]))
    assert r["leakage_plus"] < 1e-28
    assert r["leakage_minus"] == pytest.approx(1, abs=2e-14)
    assert r["max_subspace_leakage_plus"] == pytest.approx(1, abs=2e-14)
    assert r["max_subspace_leakage_minus"] == pytest.approx(1, abs=2e-14)
    assert not r["hsm_inclusion_certified"]


@pytest.mark.parametrize("ma,mb", [(3, 1), (4, 3), (7, 2), (7, 5)])
def test_full_subspace_leakage_by_independent_complex_spectral_evolution(ma, mb):
    rng = np.random.default_rng(70 + ma + mb)
    a = rng.normal(size=(ma, ma)) + 1j * rng.normal(size=(ma, ma))
    b = rng.normal(size=(mb, mb)) + 1j * rng.normal(size=(mb, mb))
    ha, hb = a + a.conj().T, b + b.conj().T
    packet = rng.normal(size=mb) + 1j * rng.normal(size=mb)
    t = 0.31
    report = modular_subspace_diagnostic(ha, hb, t, packet)
    eig, basis = np.linalg.eigh(ha)
    for sign, label in ((1, "plus"), (-1, "minus")):
        va = (basis * np.exp(1j * sign * t * eig)) @ basis.conj().T
        # No h_B in this calculation: it only rotates the input subspace.
        block = va[mb:, :mb]
        maximum = np.linalg.eigvalsh(block.conj().T @ block)[-1]
        assert report[f"max_subspace_leakage_{label}"] == pytest.approx(maximum, abs=2e-14)
        assert report[f"leakage_{label}"] <= maximum + 2e-14
    assert report["max_subspace_sign_difference"] < 2e-14


def test_reducing_subspace_is_not_promoted_to_proper_inclusion():
    r = modular_subspace_diagnostic(np.diag([1, 2, 3]), np.eye(2), 0.1, np.ones(2))
    assert r["leakage_plus"] == r["leakage_minus"] == 0
    assert r["max_subspace_leakage_plus"] == 0
    assert r["asymmetry_ratio"] is None
    assert r["lower_packet_leakage_sign"] is None
    assert not r["proper_finite_subspace_inclusion_possible"]
    json.dumps(r, allow_nan=False)


def test_imaginary_commutator_obstruction_is_exact():
    import sympy as sp
    x = sp.Matrix([[0, 1], [1, 0]])
    z = sp.diag(1, -1)
    comm = sp.I * (x * z - z * x)
    a, b, c = sp.symbols("a b c", real=True)
    candidate = sp.Matrix([[a, b], [b, c]])
    squared = sp.trace((comm - candidate).adjoint() * (comm - candidate)).expand()
    assert squared == 8 + a * a + 2 * b * b + c * c
    assert sp.trace(comm.adjoint() * comm) == 8
    r = lie_closure_receipt(32)
    assert r["real_symmetric_span_relative_distance"] == pytest.approx(1, abs=1e-15)
    assert not r["operator_lie_closure_certified"]


def test_even_envelope_has_an_exact_nonzero_invisible_commutator():
    import sympy as sp
    invisible = sp.zeros(6)
    invisible[1, 3], invisible[3, 1] = 2, -2
    invisible[0, 4], invisible[4, 0] = 1, -1
    # Both entries have midpoint 2 and weighted sum 2 - 2*1 = 0.
    diagonal = sp.diag(*range(6))
    symmetric = sp.Matrix(6, 6, lambda i, j: invisible[i, j] / (i - j) if i != j else 0)
    assert symmetric == symmetric.T
    assert diagonal * symmetric - symmetric * diagonal == invisible
    assert invisible.norm() ** 2 == 10
    profile, tail = momentum_profile(np.array(invisible).astype(float))
    np.testing.assert_array_equal(profile, np.zeros(4))
    np.testing.assert_array_equal(tail, np.zeros(4))
    base = np.zeros((6, 6))
    base[1, 3], base[3, 1] = 1, -1
    contaminated = base + 1e6 * np.array(invisible).astype(float)
    np.testing.assert_array_equal(momentum_profile(base)[0], momentum_profile(contaminated)[0])
    assert np.linalg.norm(contaminated - base) > 3e6


@pytest.mark.parametrize("cutoff", [2, 4, 8, 100, None])
def test_omitted_even_ranges_against_exact_entrywise_control(cutoff):
    rng = np.random.default_rng(172)
    entries = rng.integers(-5, 6, size=(14, 14))
    comm = entries - entries.T
    actual, tail = momentum_profile(comm, cutoff)
    expected, omitted, full = [Fraction(0)] * 12, [Fraction(0)] * 12, [Fraction(0)] * 12
    # Enumerate matrix entries instead of the producer's midpoint/range loop.
    for i in range(14):
        for k in range(i + 2, 14, 2):
            r = k - i
            j = (i + k) // 2 - 1
            term = Fraction(r, 2) * (-1) ** ((r - 2) // 2) * int(comm[i, k])
            full[j] += term
            if cutoff is None or r <= cutoff:
                expected[j] += term
            else:
                omitted[j] += abs(term)
    np.testing.assert_array_equal(actual, np.array(expected, dtype=float))
    np.testing.assert_array_equal(tail, np.array(omitted, dtype=float))
    assert all(abs(a - b) <= t for a, b, t in zip(full, expected, omitted))


@pytest.mark.parametrize("cutoff", [-2, 0, 1, 3, True, 2.0, float("nan")])
def test_bad_cutoff_fails_before_fitting(cutoff):
    with pytest.raises(ValueError, match="rmax"):
        lie_closure_receipt(32, cutoff)
    with pytest.raises(ValueError, match="rmax"):
        momentum_profile(np.zeros((4, 4)), cutoff)


@pytest.mark.parametrize("rings", [(), (16, 16), (32, 16), (True,), (16.0,), (0,), (10,), (12,)])
def test_invalid_refinement_family_is_rejected(rings):
    with pytest.raises(ValueError):
        instrument_null_net(rings)


@pytest.mark.parametrize("m,g", [([], []), ([1], [1]), ([0, 0], [1, 1]),
                                ([1, 1], [0, 0]), ([float("nan"), 1], [1, 2]),
                                ([1, 2], [float("inf"), 1]), ([1j, 1], [1, 2])])
def test_unresolved_shape_data_cannot_pass(m, g):
    with pytest.raises(ValueError):
        _shape_fit(m, g)


def test_shape_fit_is_scale_stable_and_incorrect_shape_is_detected():
    for scale in (1e-300, 1, 1e300):
        r = _shape_fit(scale * np.array([1, 2, 3]), scale * np.array([2, 4, 6]))
        assert r["relative_residual"] < 1e-15
        assert r["normalization_alpha"] == pytest.approx(0.5)
    assert _shape_fit([1, -1, 1], [1, 1, 1])["relative_residual"] > 0.9


def test_shape_fit_cannot_return_zero_coefficient_with_perfect_nonzero_fit():
    with pytest.raises(ValueError, match="normalization"):
        _shape_fit([1e-300, 2e-300], [1e300, 2e300])
    with pytest.raises(ValueError, match="normalization"):
        _shape_fit([1e300, 2e300], [1e-300, 2e-300])
    # A genuine orthogonal signal still has the resolved coefficient zero.
    r = _shape_fit([1e300, -1e300], [1e-300, 1e-300])
    assert r["normalization_alpha"] == 0
    assert r["relative_residual"] == 1


@pytest.mark.parametrize("time", [0, -0.1, True, np.bool_(True), "0.1", float("nan"), float("inf"), 1j])
def test_invalid_flow_time_is_rejected(time):
    with pytest.raises(ValueError):
        modular_subspace_diagnostic(np.eye(3), np.eye(2), time, np.ones(2))


@pytest.mark.parametrize("packet", [[0, 0], [1], [float("nan"), 1], [float("inf"), 1], [True, False]])
def test_invalid_flow_packet_is_rejected(packet):
    with pytest.raises(ValueError):
        modular_subspace_diagnostic(np.eye(3), np.eye(2), 0.1, packet)


def test_nonhermitian_flow_and_invalid_embedding_are_rejected():
    with pytest.raises(ValueError, match="Hermitian"):
        modular_subspace_diagnostic(np.array([[1, 2], [0, 1]]), np.ones((1, 1)), 1, [1])
    with pytest.raises(ValueError, match="subspace"):
        modular_subspace_diagnostic(np.eye(2), np.eye(2), 1, [1, 1])
    for size, offset in ((1, 0), (4, -1), (3, 2), (True, 0)):
        with pytest.raises(ValueError):
            embed(np.eye(2), size, offset)
    y = np.array([[0, -1j], [1j, 0]])
    np.testing.assert_array_equal(embed(y, 4, 1)[1:3, 1:3], y)


def test_nonfinite_matrix_and_unresolved_spectrum_are_rejected(monkeypatch):
    import mpmath as mp
    import null_net_receipts as module
    for matrix in (np.array([[float("nan")]]), np.array([[float("inf")]]),
                   np.array([1, 2]), np.zeros((2, 3)), np.zeros((0, 0))):
        with pytest.raises(ValueError):
            embed(matrix, 4)
    monkeypatch.setattr(module.mp, "eigsy", lambda _: (mp.matrix([0, 1]), None))
    with pytest.raises(ValueError, match="unresolved"):
        separating_modulus_receipt(8)


def test_checked_report_replays_with_scoped_verdicts():
    expected = json.loads((HERE / "runs" / "null_net_receipt_report.json").read_text(encoding="utf-8"))
    actual = instrument_null_net()

    def compare(a, b):
        if isinstance(a, dict):
            assert a.keys() == b.keys()
            for key in a:
                compare(a[key], b[key])
        elif isinstance(a, list):
            assert len(a) == len(b)
            for x, y in zip(a, b):
                compare(x, y)
        elif type(a) is float:
            assert b == pytest.approx(a, abs=2e-12, rel=2e-10)
        else:
            assert type(a) is type(b) and a == b

    compare(expected, actual)
    assert actual["schema_version"] == 2
    assert actual["verdicts"]["lie_closure_residuals_decreasing"]
    assert not actual["verdicts"]["lie_closure_rate_certified"]
    assert all(r["all_even_ranges"] for r in actual["lie_closure"])
    json.dumps(actual, allow_nan=False)

"""Tests for the Z2 finite ground-state-transform diagnostic."""

from __future__ import annotations

import json
import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import z2_finite_transfer_receipt as z2  # noqa: E402
from verify_z2_finite_transfer_receipt import load_receipt, verify_receipt  # noqa: E402


@pytest.fixture(scope="module")
def orbits() -> z2.Z2GaugeOrbits:
    return z2.Z2GaugeOrbits(2)


def test_self_loop_lattice_is_rejected() -> None:
    with pytest.raises(ValueError, match="repeated self-loop incidences cancel"):
        z2.Z2GaugeOrbits(1)
    with pytest.raises(ValueError, match="requires L >= 2"):
        z2.star_masks(1)
    with pytest.raises(ValueError, match="requires L >= 2"):
        z2.plaquette_masks(1)


def test_orbit_count_and_free_action(orbits: z2.Z2GaugeOrbits) -> None:
    # 8 links, gauge group Z2^4 / global = 8 elements, free action
    assert orbits.n_links == 8
    assert len(orbits.gauge) == 8
    assert orbits.n_orbits == 32
    assert np.all(orbits.orbit_size == 8)


def test_free_control_receipt_is_exact(orbits: z2.Z2GaugeOrbits) -> None:
    """At beta_s = 0 the kinetic kernel factorises as
    prod_l (e^{beta_t} I + e^{-beta_t} X_l) = const * exp(b sum_l X_l) with
    b = log(coth beta_t) / 2, so the receipt holds exactly with every rate
    equal to log(coth beta_t)."""
    r = z2.evaluate(orbits, "wilson", beta_s=0.0, beta_t=0.5)
    fit = r["constant_rate_fit"]
    rate = math.log(1.0 / math.tanh(0.5))
    assert fit["relative_frobenius_residual"] < 1e-10
    assert fit["rate_min"] == pytest.approx(rate, abs=1e-9)
    assert fit["rate_max"] == pytest.approx(rate, abs=1e-9)
    assert r["spectral"]["gap_H"] == pytest.approx(2 * rate, abs=1e-9)
    assert r["fiber_dependent_rates"]["spread_max_over_min"] == pytest.approx(1.0, abs=1e-9)
    assert r["dobrushin"]["eta_star"] < 1e-9
    assert r["spectral"]["gap_unit_rate_heat_bath"] == pytest.approx(2.0, abs=1e-9)


def test_interacting_wilson_receipt_fails(orbits: z2.Z2GaugeOrbits) -> None:
    r = z2.evaluate(orbits, "wilson", beta_s=0.5, beta_t=0.5)
    assert r["doob_generator_rows_sum_zero"]
    # Row conservation alone is not the Markov sign condition; the Wilson
    # logarithm has positive off-diagonal entries at this interacting point.
    assert not r["doob_generator_offdiagonal_nonpositive"]
    assert r["constant_rate_fit"]["relative_frobenius_residual"] > 1e-2
    assert r["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] > 1.0
    assert r["fiber_dependent_rates"]["spread_max_over_min"] > 1.01


def test_kogut_susskind_single_flip_exact_but_fiber_dependent(orbits: z2.Z2GaugeOrbits) -> None:
    r = z2.evaluate(orbits, "kogut_susskind", lam=1.0)
    fib = r["fiber_dependent_rates"]
    assert r["doob_generator_rows_sum_zero"]
    assert r["doob_generator_offdiagonal_nonpositive"]
    assert fib["offdiagonal_mass_outside_single_flip"] < 1e-9
    assert fib["all_rates_positive"]
    assert fib["spread_max_over_min"] > 1.01
    # exact rate formula c = lam (r + 1/r) with r = Omega(o)/Omega(o') >= 2 lam
    assert fib["rate_min"] >= 2.0 - 1e-9
    floor = r["variable_rate_floor"]
    assert floor["lower_bound_value"] == pytest.approx(2.0)
    assert floor["numerical_min_respects_bound"] is True


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_rate_fit_resolves_each_coefficient_across_units(orbits, scale) -> None:
    projectors = z2.heat_bath_projectors(orbits, np.full(orbits.n_orbits, 1 / orbits.n_orbits))
    expected = np.arange(1, orbits.n_links + 1, dtype=float)
    generator = sum(c * (np.eye(orbits.n_orbits) - E)
                    for c, E in zip(expected, projectors))
    result = z2.constant_rate_fit(scale * generator, projectors)
    np.testing.assert_allclose(np.asarray(result["rates"]) / scale, expected, rtol=2e-14)
    assert result["relative_frobenius_residual"] < 2e-14
    # A real off-support perturbation must not disappear through a zero norm
    # or through a scale-independent absolute tolerance.
    generator[0, -1] += 0.3
    generator[0, 0] -= 0.3
    perturbed = z2.constant_rate_fit(scale * generator, projectors)
    assert perturbed["relative_frobenius_residual"] > 1e-3


@pytest.mark.parametrize("beta_s,beta_t", [
    (float("nan"), 0.5), (float("inf"), 0.5), (True, 0.5),
    (0.0, 0.0), (0.0, -1.0), (0.0, float("nan")),
    (0.0, float("inf")), (0.0, True), (0.0, 0.5j),
])
def test_wilson_refuses_nonphysical_or_nonfinite_parameters(orbits, beta_s, beta_t) -> None:
    with pytest.raises(ValueError):
        z2.evaluate(orbits, "wilson", beta_s=beta_s, beta_t=beta_t)


@pytest.mark.parametrize("matrix", [
    np.eye(2, dtype=complex), np.array([[1.0, 0.2], [0.3, 1.0]]),
    np.array([[float("nan"), 0.1], [0.1, 1.0]]), np.zeros((2, 2)),
    np.ones((2, 3)), np.ones((1, 1)), np.ones((2, 2), dtype=bool),
])
def test_generic_transfer_rejects_wrong_matrix_domain(matrix) -> None:
    with pytest.raises(ValueError):
        z2.symmetric_log_hamiltonian(matrix)


def test_unrepresentable_wilson_outputs_are_not_silent_zeros(orbits) -> None:
    with pytest.raises(ValueError, match="normalization.*range"):
        z2.evaluate(orbits, "wilson", beta_s=0.0, beta_t=100.0)
    with pytest.raises(RuntimeError, match="dual coupling.*range"):
        z2.dual_coupling(400.0)


def test_committed_receipt_matches_code() -> None:
    path = HERE / "receipts" / "z2_finite_transfer_receipt.json"
    receipt = load_receipt(path)
    assert receipt["schema"] == z2.SCHEMA
    assert receipt["physical_clay_receipt"] is False
    assert receipt["grid_scope"]["universal_no_go"] is False
    assert receipt["producer_sha256"] == hashlib.sha256(
        (HERE / "z2_finite_transfer_receipt.py").read_bytes()
    ).hexdigest()
    controls = [r for r in receipt["runs"] if r.get("control") == "beta_s_zero"]
    assert controls and all(
        r["constant_rate_fit"]["relative_frobenius_residual"] < 1e-10 for r in controls
    )
    interacting = [r for r in receipt["runs"] if r["transfer"] == "wilson" and "control" not in r]
    assert interacting and all(
        r["constant_rate_fit"]["relative_frobenius_residual"] > 1e-3 for r in interacting
    )
    local_runs = [r for r in receipt["runs"] if r["transfer"] == "kogut_susskind"]
    assert local_runs and all(
        r["variable_rate_floor"]["numerical_min_respects_bound"] for r in local_runs
    )
    fresh = z2.run([2, 3], [0.1, 0.3, 0.5, 0.7, 1.0], [0.5, 1.0, 2.0])
    # Floating eigensolver output need not be byte-identical across LAPACK/OS.
    # Verify every component with relative tolerances; hashes still bind the
    # exact retained bytes and the original analytic controls run separately.
    verify_receipt(receipt, fresh, HERE / "z2_finite_transfer_receipt.py")
    expected_runs_hash = hashlib.sha256(
        json.dumps(receipt["runs"], sort_keys=True).encode("utf-8")
    ).hexdigest()
    assert receipt["sha256_of_runs"] == expected_runs_hash

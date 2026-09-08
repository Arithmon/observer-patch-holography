#!/usr/bin/env python3
"""Finite profile checks and regression controls for modular instrumentation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from modular_clock_instrumentation import (  # noqa: E402
    arc_entanglement_hamiltonian,
    cft_profile,
    crossratio_receipt,
    instrument_tower,
    kms_profile_receipt,
    resummed_profile,
    ring_correlation_value,
    stage_data,
    truncation_tail_bound,
)
import modular_clock_instrumentation as clock  # noqa: E402


@pytest.fixture(scope="module")
def full_tower():
    # One real four-rung eigensolve run, including the cutoff that falsifies
    # the previous range-nine monotonicity claim. No source/data cache.
    return instrument_tower()


def test_correlations_are_half_filled_and_exact():
    # C(0) = 1/2 exactly at half filling
    assert abs(float(ring_correlation_value(16, 0)) - 0.5) < 1e-30


def test_extended_precision_resolves_the_spectrum():
    # float64 saturates interior EH bonds near |h| ~ 18; the extended-
    # precision computation must exceed that on the 64-ring, whose true
    # midpoint bond value is ~ beta_CFT/2 = 32
    h = arc_entanglement_hamiltonian(64, 32)
    interior = max(abs(h[j, j + 1]) for j in range(14, 18))
    assert interior > 25.0


def test_profile_matches_2pi_bw_and_wrong_normalization_separated():
    data = stage_data(32)
    receipt = kms_profile_receipt(data, 32)
    assert receipt["median_relative_residual_2pi"] < 5e-3
    assert receipt["separation_factor"] > 5.0


def test_nearest_neighbour_truncation_alone_fails():
    # The supplied free-chain diagnostic has a visible range-one defect;
    # this numerical comparison is not a physical KMS test.
    h = arc_entanglement_hamiltonian(32, 16)
    beta_r1 = resummed_profile(h, rmax=1)
    xs = np.arange(15) + 0.5
    beta_geo = cft_profile(32, xs, -0.5, 15.5)
    mid = slice(5, 10)
    rel = np.abs(beta_r1[mid] / beta_geo[mid] - 1.0)
    assert np.median(rel) > 0.05


def test_crossratio_error_improves_on_two_finite_cutoffs():
    r16 = crossratio_receipt(stage_data(16), 16)
    r32 = crossratio_receipt(stage_data(32), 32)
    assert r32["relative_error"] < r16["relative_error"]
    assert r32["relative_error"] < 0.05
    # the Moebius target is stage-independent for proportional quadruples
    assert abs(r16["cr_moebius"] - r32["cr_moebius"]) < 1e-9


def test_instrumented_tower_is_a_finite_diagnostic(full_tower):
    report = full_tower
    w = report["receipts_witnessed"]
    assert w["geometric_2pi_kms_boundary_collar"] is False
    assert w["modular_cross_ratio_boundary_collar"] is False
    assert all(report["finite_diagnostics"].values())
    assert report["rings"] == [16, 32, 64, 128]
    assert report["resummation_rmax_by_ring"] == [7, 15, 31, 63]
    assert report["error_scope"]["finite_range_truncation_removed"] is True
    assert report["error_scope"]["roundoff_interval_certified"] is False
    assert report["error_scope"]["continuum_profile_error_certified"] is False
    assert "boundary-collar" in report["scope"]
    assert len(report["receipts_pending"]) >= 4


def test_fixed_range_n128_falsifies_the_old_three_rung_trend(full_tower):
    fixed = full_tower["fixed_depth_control"]
    rows = fixed["kms_profile_receipt"]
    errors = [row["median_relative_residual_2pi"] for row in rows]
    assert errors[0] > errors[1] > errors[2]
    assert errors[3] > 2.9 * errors[2]
    assert fixed["verdicts"]["kms_residual_decreasing"] is False
    assert fixed["finite_diagnostics"]["profile_agreement_with_decreasing_finite_residual"] is False
    assert rows[3]["max_relative_residual_2pi"] > 0.002
    # A positive cross-ratio trend does not rescue the failed profile test.
    assert fixed["verdicts"]["crossratio_error_decreasing"] is True


def test_all_ranges_remove_the_finite_truncation_artifact(full_tower):
    rows = full_tower["kms_profile_receipt"]
    errors = [row["median_relative_residual_2pi"] for row in rows]
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert 7.0e-5 < errors[-1] < 7.2e-5
    assert rows[-1]["max_relative_residual_2pi"] < 1.01e-4
    assert full_tower["crossratio_receipt"][-1]["relative_error"] < 0.00068
    tails = full_tower["fixed_depth_control"]["omitted_range_diagnostics"]
    assert tails[0]["max_absolute_omitted_range_bound"] == 0
    assert tails[-1]["max_relative_omitted_range_bound"] > 0.002
    for row in tails:
        # Floating evaluation of an exact triangle bound is not an interval
        # certificate; allow ordinary roundoff only in this comparison.
        assert row["max_absolute_full_minus_fixed_profile"] <= row["max_absolute_omitted_range_bound"] + 1e-12


def test_full_sum_and_tail_against_independent_edge_enumeration():
    # Integer entries make all arithmetic exact in float64. Enumerate edges,
    # independently of the producer's loop over bond centers.
    h = np.arange(16 * 16, dtype=float).reshape(16, 16)
    full = np.zeros(15)
    retained = np.zeros(15)
    omitted = np.zeros(15)
    for i in range(16):
        for k in range(i + 1, 16):
            r = k - i
            if r % 2 == 0:
                continue
            j = (i + k - 1) // 2
            contribution = -2 * r * ((-1) ** ((r - 1) // 2)) * h[i, k]
            full[j] += contribution
            if r <= 9:
                retained[j] += contribution
            else:
                omitted[j] += abs(contribution)
    np.testing.assert_array_equal(resummed_profile(h), full)
    np.testing.assert_array_equal(resummed_profile(h, 9), retained)
    np.testing.assert_array_equal(truncation_tail_bound(h, 9), omitted)
    np.testing.assert_array_equal(truncation_tail_bound(h), np.zeros(15))
    assert np.all(abs(full - retained) <= omitted)


def test_omitted_long_range_cannot_silently_disappear():
    h = np.zeros((12, 12))
    h[0, 11] = h[11, 0] = 1
    assert resummed_profile(h)[5] == 22
    assert resummed_profile(h, 9)[5] == 0
    assert truncation_tail_bound(h, 9)[5] == 22


@pytest.mark.parametrize("rings", [(), (16,), (16, 32), (16, 16, 32), (32, 16, 64), (16, 32, True), (16, 32, 66), (16, 32, 68)])
def test_invalid_or_vacuous_towers_are_rejected(rings):
    with pytest.raises(ValueError):
        instrument_tower(rings)


@pytest.mark.parametrize("cutoff", [0, -1, 2, 9.0, True])
def test_invalid_depth_rejected(cutoff):
    with pytest.raises(ValueError):
        resummed_profile(np.eye(4), cutoff)


@pytest.mark.parametrize("h", [np.eye(4) * np.nan, np.ones((3, 4)), np.eye(4) * 1j, np.ones(4), np.zeros((1, 1))])
def test_invalid_hamiltonian_rejected(h):
    with pytest.raises(ValueError):
        resummed_profile(h)


@pytest.mark.parametrize("value", [0, -1, float("inf"), float("nan")])
def test_invalid_transport_profile_fails_closed(value):
    data = {"m": 8, "beta": np.ones(7), "u": -.5, "v": 7.5}
    data["beta"][3] = value
    with pytest.raises(ValueError):
        crossratio_receipt(data, 16)


def test_unresolved_eigenspectrum_fails_closed(monkeypatch):
    import mpmath as mp
    monkeypatch.setattr(clock.mp, "eigsy", lambda _: (mp.matrix([0, 1]), mp.eye(2)))
    with pytest.raises(ValueError, match="unresolved"):
        arc_entanglement_hamiltonian(8, 2)


def test_quarter_arc_readout_requires_integral_positions():
    with pytest.raises(ValueError, match="divisible by eight"):
        stage_data(12)


def test_eigensolve_precision_is_local_and_caller_precision_is_restored():
    import mpmath as mp
    reference = arc_entanglement_hamiltonian(16, 8)
    with mp.workdps(8):
        computed = arc_entanglement_hamiltonian(16, 8)
        assert mp.mp.dps == 8
    np.testing.assert_array_equal(computed, reference)

"""Focused tests for the near-hit attribution surface."""

from __future__ import annotations

import pytest

import derive_near_hit_attribution_surface as lane


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    out = tmp_path_factory.mktemp("nearhit") / "near_hit_attribution_surface.json"
    return lane.build(out)


def test_every_row_names_a_hypothesis_and_reduction(result):
    assert len(result["rows"]) >= 10
    for row in result["rows"]:
        hyp = row["missing_correction_hypothesis"]
        assert hyp["mechanism"]
        assert hyp["falsification_test"]
        assert hyp["reduces_to"]


def test_surface_is_compare_only(result):
    guards = result["guards"]
    assert guards["compare_only"] is True
    assert guards["public_promotion_allowed"] is False
    assert guards["changes_any_solve_path"] is False
    assert guards["new_axiom_introduced"] is False


def test_wz_chart_rows_fail_closed(result):
    assert {"MW_chart_gev", "MZ_chart_gev", "mt_pole_gev"} <= set(
        r["quantity"] for r in result["rows"]
    )
    charts = [
        row for row in result["rows"] if row["quantity"] in {"MW_chart_gev", "MZ_chart_gev"}
    ]
    for row in charts:
        assert row["physical_comparison_status"] == "NOT_EVALUABLE"
        assert row["physical_delta"] is None
        assert row["physical_pull"] is None
        assert "delta_over_sigma" not in row
        assert row["legacy_reference_experimental_sigma"] > 0.0


def test_alpha_rows_present_with_expected_signs(result):
    by_name = {r["quantity"]: r for r in result["rows"]}
    anchor = by_name["alpha_em_inv_MZ_anchor"]
    assert anchor["delta"] < 0.0, "anchor must sit below the physical on-shell value"
    endpoint = by_name["alpha_inv_thomson_endpoint"]
    assert endpoint["delta"] < 0.0, "endpoint must sit below CODATA"
    wz = by_name["alpha_em_eff_inv_on_shell_wz_lane"]
    assert wz["delta"] == pytest.approx(0.11796769814804975, rel=1e-9)


def test_lepton_rows_are_near_hits(result):
    lepton_rows = [r for r in result["rows"] if r["quantity"].startswith("m_")]
    assert len(lepton_rows) == 3
    for row in lepton_rows:
        assert abs(row["relative_delta"]) < 0.02
        lo, hi = row["oph_envelope"]
        assert lo < row["measured"] < hi


def test_mz_legacy_coordinate_diagnostic_has_no_physical_verdict(result):
    mz = next(r for r in result["rows"] if r["quantity"] == "MZ_chart_gev")
    executed = mz["missing_correction_hypothesis"]["legacy_coordinate_diagnostic"]
    assert executed["physical_comparison_status"] == "NOT_EVALUABLE"
    naive = executed["naive_coupling_injection"]
    assert naive["verdict"] == "LEGACY_DIAGNOSTIC_ONLY"
    assert all(delta < -50.0 for delta in naive["delta_mz_mev_over_gap_interval"])
    p_channel = executed["p_channel"]
    assert p_channel["verdict"] == "LEGACY_DIAGNOSTIC_ONLY"
    assert p_channel["p_calibration"] == 1.63094
    assert p_channel["p_residual_after_correction"] == pytest.approx(
        p_channel["p_corrected_endpoint_to_measured"] - p_channel["p_calibration"]
    )
    # The former near-zero residual compared to a label that did not generate the tuple.
    assert abs(p_channel["p_residual_after_correction"]) > 1e-5
    assert p_channel["mz_residual_after_correction_mev"] == pytest.approx(
        p_channel["p_residual_after_correction"] * p_channel["dmz_dp_gev_per_unit_p"] * 1e3
    )


def test_synthesis_names_single_dominant_reduction(result):
    synthesis = result["synthesis"]
    assert "#425" in synthesis["dominant_reduction"]
    assert "#545" in synthesis["dominant_reduction"]
    assert "#521" in synthesis["secondary_reduction"]

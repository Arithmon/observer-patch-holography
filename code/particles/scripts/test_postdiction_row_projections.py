"""Cheap projection checks; these never invoke the full structural-ledger build."""
import json

import build_postdiction_ledger as ledger


def test_lepton_projection_preserves_calibration_and_design_ancestry():
    surface = ledger._load("mass_surface")
    rows = ledger._lepton_rows(surface, ledger._load("kappa_rectangle"),
                              ledger._load("kappa_coherent"), ledger._load("koide_balance"),
                              ledger._load("koide_identity"))
    by_id = {row["id"]: row for row in rows}
    closure = by_id["charged_leptons_closure_target"]
    assert "empirical calibration" in closure["statement"]
    assert "measured-lepton information" in closure["statement"]
    mcpr = by_id["charged_leptons_mcpr_conditional"]
    parent = next(row for family in surface["families"] if family["family"] == "charged leptons"
                  for row in family["rows"] if row["lane"].startswith("MCPR"))
    assert mcpr["P"] == parent["P"]
    assert mcpr["historical_target_informed_architecture"] is True
    assert mcpr["scale_and_running_boundary"] == surface["display_note"]
    assert "three circulant" in mcpr["explanation"]
    assert not mcpr["dimensionless_ratio_comparison"]["m_mu_over_m_e"]["inside_componentwise_one_sigma_interval"]


def test_electroweak_projection_carries_actual_pixel_and_partial_scope():
    conditional, surface = ledger._load("conditional_ew"), ledger._load("mass_surface")
    rows = ledger._ew_rows(conditional, surface)
    assert len(rows) == 4
    for row in rows:
        assert row["P_basis"] == conditional["inputs"]["P_calibration"] == 1.63094
        assert row["historical_incorrect_P_label"] != row["P_basis"]
        assert row["sensitivity_scope"] == conditional["inputs"]["sensitivity_scope"]
        assert row["scale_and_running_boundary"] == surface["display_note"]
        if row["id"] in {"ew_mH_gev", "ew_mt_pole_gev"}:
            assert row["anchoring"] == "target_anchored_backsolve"
        else:
            assert row["physical_comparison_status"] == "NOT_EVALUABLE"


def test_principal_digest_keeps_calibration_ancestry_and_gauge_matter_premises():
    sections = json.loads(ledger.DEFAULT_OUT.read_text())["sections"]
    rows = {row["id"]: row for row in ledger._principal_results(sections)}
    gauge = rows["forced_gauge_structure"]["statement"]
    for premise in ("faithful twelve-dimensional compact", "commutator-closed",
                    "centre", "compact-simple classification", "spin lift",
                    "central trace redistribution", "nonvacuum even exterior sector"):
        assert premise in gauge
    assert "does not supply these response and centre hypotheses" in gauge
    closure = rows["lepton_closure_target"]["statement"]
    assert "measured-lepton information" in closure
    assert "empirical calibration" in closure
    assert "separately frozen inputs, matching and uncertainty" in closure
    assert "is a falsification target" not in closure

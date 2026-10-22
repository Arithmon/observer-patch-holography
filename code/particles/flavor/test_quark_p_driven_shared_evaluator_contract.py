#!/usr/bin/env python3
"""Guard the shared P-driven quark evaluator contract."""

from __future__ import annotations

import json
import math
from pathlib import Path

from p_driven_flavor_candidate import (
    ANCHOR_ALPHA_U,
    EXACT_SIGMA_TARGET_D,
    EXACT_SIGMA_TARGET_U,
    build_shared_p_driven_evaluator_contract,
    shared_candidate_quark_masses_from_alpha_u,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_p_driven_shared_evaluator_contract.json"


def test_shared_p_driven_contract_refuses_theorem_promotion_with_live_blockers() -> None:
    payload = json.loads(CONTRACT_JSON.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_p_driven_shared_evaluator_contract"
    assert payload["proof_status"] == "candidate_only"
    assert payload["runtime_status"] == "shared_candidate_evaluator"
    assert payload["public_promotion_allowed"] is False
    assert payload["theorem_grade_closure"] is False
    assert payload["input_statuses"]["pure_b_source_status"] == (
        "comparison_only_values_from_unpromotable_source_projection"
    )
    assert set(payload["promotion_blockers"]) == {
        "default_universe_anchor_not_removed",
        "edge_statistics_bridge_not_closed",
        "off_canonical_odd_response_not_closed",
        "pure_B_source_payload_not_closed",
    }


def test_shared_p_driven_contract_recovers_exact_default_quark_anchor() -> None:
    surface = shared_candidate_quark_masses_from_alpha_u(ANCHOR_ALPHA_U)
    assert math.isclose(surface["sigma_u_total_log_per_side"], EXACT_SIGMA_TARGET_U, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(surface["sigma_d_total_log_per_side"], EXACT_SIGMA_TARGET_D, rel_tol=0.0, abs_tol=1.0e-12)
    for row in surface["up_sector"] + surface["down_sector"]:
        assert math.isclose(row["mass_gev"], row["baseline_mass_gev"], rel_tol=0.0, abs_tol=1.0e-12)


def test_shared_p_driven_contract_moves_off_anchor_without_relabeling_status() -> None:
    surface = shared_candidate_quark_masses_from_alpha_u(ANCHOR_ALPHA_U * 1.05)
    assert any(
        not math.isclose(row["mass_gev"], row["baseline_mass_gev"], rel_tol=0.0, abs_tol=1.0e-12)
        for row in surface["up_sector"] + surface["down_sector"]
    )
    contract = build_shared_p_driven_evaluator_contract(
        edge_statistics_bridge_status="candidate_only",
        odd_response_proof_status="candidate_only",
        pure_b_source_status="source_values_shell_waiting_pure_B_payload_pair",
    )
    assert contract["theorem_grade_closure"] is False
    assert contract["public_promotion_allowed"] is False


def test_shared_p_driven_contract_distinguishes_selected_class_pure_b_from_off_canonical_family() -> None:
    contract = build_shared_p_driven_evaluator_contract(
        edge_statistics_bridge_status="candidate_only",
        odd_response_proof_status="candidate_only",
        pure_b_source_status="source_values_derived_from_source_emission",
    )
    assert "pure_B_source_payload_not_closed" not in contract["promotion_blockers"]
    assert "off_canonical_pure_B_payload_family_not_closed" in contract["promotion_blockers"]

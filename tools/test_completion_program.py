"""Fail-closed tests for completion scope, dependencies and retained foundations."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

import build_completion_program as program


@pytest.fixture
def fixture_root(tmp_path):
    source = json.loads((program.ROOT / program.PROGRAM).read_text())
    hypothesis = json.loads((program.ROOT / program.HYPOTHESIS).read_text())
    claims = json.loads((program.ROOT / "claims/claim_registry.yaml").read_text())["claims"]
    needed = {row["claim_id"] for row in source["foundations"]}
    source_claims = [{"claim_id": row["claim_id"], "status": row["status"]}
                     for row in claims if row["claim_id"] in needed]
    paths = {hypothesis["contract"]}
    for row in source["stages"] + source["foundations"]:
        paths.update(row["evidence"])
    bound_rows = [row for stage in source["stages"] for row in stage.get("progress_evidence", [])]
    bound_rows += [row["contract"] for row in source.get("component_campaigns", [])]
    bound_rows += source.get("next_wave", {}).get("evidence", [])
    for row in bound_rows:
        paths.update((row["path"], row["independent_check"]))
    for name in paths:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("Independent fixture evidence: " + name + "\n")
    for campaign in source.get("component_campaigns", []):
        name = campaign["contract"]["path"]
        (tmp_path / name).write_text((program.ROOT / name).read_text())
    for row in bound_rows:
        row["sha256"] = hashlib.sha256((tmp_path / row["path"]).read_bytes()).hexdigest()
    for name, value in ((program.PROGRAM, source), (program.HYPOTHESIS, hypothesis),
                        ("claims/claim_registry.yaml", {"claims": source_claims})):
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(value))
    return tmp_path


def change(root, name, mutate):
    path = root / name
    data = json.loads(path.read_text())
    mutate(data)
    path.write_text(json.dumps(data))


def test_retained_foundations_do_not_promote_integrated_gates(fixture_root):
    view = program.ledger_view(fixture_root)
    assert len(view["retained_foundations"]) == 5
    assert any("modular" in row["title"] for row in view["retained_foundations"])
    assert any("3+1" in row["title"] for row in view["retained_foundations"])
    assert all(row["canonical_claim_status"] for row in view["retained_foundations"])
    assert all(row["source_sha256"] for row in view["retained_foundations"])
    assert all(row["status"] in {"OPEN", "IN_PROGRESS"} for row in view["program"]["stages"])
    assert view["hypothesis"]["comparison_count_added"] == 0
    assert "not scientific gate acceptance" in view["verification"]


@pytest.mark.parametrize("damage", ["cycle", "missing_dependency", "duplicate_owner",
                                  "accepted_without_evidence", "premature_freeze",
                                  "whole_program_no_coverage", "unregistered_foundation"])
def test_damaged_program_cannot_pass(fixture_root, damage):
    def mutate(data):
        stages = {row["id"]: row for row in data["stages"]}
        if damage == "cycle":
            stages["CP-00"]["requires"] = ["CP-12"]
        elif damage == "missing_dependency":
            stages["CP-03"]["requires"].append("NONEXISTENT")
        elif damage == "duplicate_owner":
            stages["CP-01"]["issue"] = stages["CP-00"]["issue"]
        elif damage == "accepted_without_evidence":
            stages["CP-02"]["status"] = "ACCEPTED"
        elif damage == "premature_freeze":
            data["boundary_status"] = "FROZEN"
        elif damage == "whole_program_no_coverage":
            data["decision"] = "FALSIFIED_PROGRAM"
            data["whole_program_falsified"] = True
        else:
            data["foundations"][0]["claim_id"] = "INVENTED_CLAIM"
    change(fixture_root, program.PROGRAM, mutate)
    with pytest.raises(ValueError):
        program.load_program(fixture_root)


@pytest.mark.parametrize("flag", ["global_positive_mass_gap_claimed",
                                  "modular_flow_automatically_physical_time",
                                  "all_particles_stable_defects_established",
                                  "einstein_premises_realized"])
def test_unjustified_physical_promotion_fails(fixture_root, flag):
    change(fixture_root, program.HYPOTHESIS, lambda data: data.update({flag: True}))
    with pytest.raises(ValueError, match="unsupported hypothesis promotion"):
        program.load_program(fixture_root)


def test_hash_bound_artifact_does_not_replace_predecessor_acceptance(fixture_root):
    p = fixture_root / "evidence.txt"
    p.write_text("Bound scope witness\n")
    def mutate(data):
        stage = next(row for row in data["stages"] if row["id"] == "CP-02")
        stage["status"] = "ACCEPTED"
        stage["accepted_evidence"] = [{"path": p.name,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "independent_check": p.name, "scope": "fixture clock witness"}]
    change(fixture_root, program.PROGRAM, mutate)
    with pytest.raises(ValueError, match="before predecessors"):
        program.load_program(fixture_root)


def test_escape_or_missing_artifact_is_rejected(fixture_root):
    change(fixture_root, program.PROGRAM,
           lambda data: data["foundations"][0]["evidence"].append("../outside.txt"))
    with pytest.raises(ValueError, match="missing or escaping"):
        program.load_program(fixture_root)


def test_committed_issue_and_dag_view_is_current():
    data, hypothesis = program.load_program()
    assert (program.ROOT / program.OUTPUT).read_text() == program.render(data, hypothesis)
    assert all(row["issue"] for row in data["stages"])


@pytest.mark.parametrize("damage", ["custody", "unmapped", "promotion"])
def test_component_progress_cannot_change_scope_or_lose_custody(fixture_root, damage):
    def mutate(data):
        campaign = data["component_campaigns"][0]
        if damage == "custody":
            campaign["contract"]["sha256"] = "0" * 64
        elif damage == "promotion":
            campaign["whole_program_frozen"] = True
        else:
            next(row for row in data["stages"] if row.get("progress_evidence"))["progress_evidence"][0]["campaign_id"] = "UNKNOWN"
    change(fixture_root, program.PROGRAM, mutate)
    with pytest.raises(ValueError):
        program.load_program(fixture_root)


def test_component_progress_must_belong_to_campaign_supported_gate(fixture_root):
    def mutate(data):
        campaign = data["component_campaigns"][0]
        campaign["supports"].remove("CP-00")
    change(fixture_root, program.PROGRAM, mutate)
    with pytest.raises(ValueError, match="campaign for this gate"):
        program.load_program(fixture_root)


def test_next_wave_uses_four_existing_owners_without_acceptance_promotion(fixture_root):
    data, hypothesis = program.load_program(fixture_root)
    wave = data["next_wave"]
    assert {(row["gate"], row["issue"]) for row in wave["owners"]} == {
        ("CP-03", 1139), ("CP-04", 1140), ("CP-06", 1142), ("CP-13", 1137)}
    assert {row["gate"] for row in wave["embedded_reviews"]} == {
        "CP-00", "CP-01", "CP-02", "CP-05", "CP-07", "CP-08", "CP-09", "CP-12"}
    assert not wave["changes_frozen_campaign"]
    assert not wave["successor_class_frozen"]
    assert wave["physical_comparison_count_added"] == 0
    assert all(stage["status"] != "ACCEPTED" for stage in data["stages"])
    assert data["boundary_status"] == "NOT_FROZEN"
    rendered = program.render(data, hypothesis)
    assert "## Next four parallel workstreams" in rendered
    assert wave["acceptance_policy"] in rendered
    assert all(row["deliverable"] in rendered for row in wave["owners"])
    by_id = {stage["id"]: stage for stage in data["stages"]}
    for owner in wave["owners"]:
        issue = program.render_issue(by_id[owner["gate"]], data)
        assert owner["deliverable"] in issue
        assert owner["join_condition"] in issue


@pytest.mark.parametrize("damage", ["duplicate_owner", "wrong_issue", "missing_evidence",
                                  "wrong_baseline", "unilateral_primitive", "successor_freeze",
                                  "physical_comparison", "accepted_assignment", "review_owner",
                                  "empty_deliverable", "stale_evidence"])
def test_next_wave_cannot_lose_ownership_or_silently_widen_scope(fixture_root, damage):
    def mutate(data):
        wave = data["next_wave"]
        if damage == "duplicate_owner":
            wave["owners"][1] = deepcopy(wave["owners"][0])
        elif damage == "wrong_issue":
            wave["owners"][0]["issue"] = 1147
        elif damage == "missing_evidence":
            wave["owners"][0]["evidence"] = ["UNBOUND"]
        elif damage == "wrong_baseline":
            wave["shared_interface"]["baseline_contract"] = "other_campaign.json"
        elif damage == "unilateral_primitive":
            wave["shared_interface"]["unilateral_primitive_changes_allowed"] = True
        elif damage == "successor_freeze":
            wave["successor_class_frozen"] = True
        elif damage == "physical_comparison":
            wave["physical_comparison_count_added"] = 1
        elif damage == "accepted_assignment":
            wave["status"] = "ACCEPTED"
        elif damage == "review_owner":
            wave["embedded_reviews"][0]["issue"] = 9999
        elif damage == "empty_deliverable":
            wave["owners"][0]["deliverable"] = ""
        else:
            wave["evidence"][0]["sha256"] = "0" * 64
    change(fixture_root, program.PROGRAM, mutate)
    with pytest.raises(ValueError):
        program.load_program(fixture_root)


def test_readiness_distinguishes_work_underway_from_acceptance_dependencies(fixture_root):
    data, _ = program.load_program(fixture_root)
    by_id = {stage["id"]: stage for stage in data["stages"]}
    assert program.readiness(by_id["CP-00"], data) == "In progress"
    assert program.readiness(by_id["CP-02"], data) == "In progress; acceptance awaiting CP-01, CP-13"
    not_started = deepcopy(by_id["CP-00"])
    not_started["status"] = "OPEN"
    assert program.readiness(not_started, data) == "Ready to start"
    assert program.readiness(by_id["CP-04"], data) == "In progress; acceptance awaiting CP-02, CP-03"


def test_finite_diagnostics_register_progress_without_integrated_acceptance(fixture_root):
    data, _ = program.load_program(fixture_root)
    diagnostics = [row for row in data["component_campaigns"]
                   if row["scope"] == "FINITE_COMPONENT_DIAGNOSTICS"]
    assert len(diagnostics) >= 2
    by_id = {stage["id"]: stage for stage in data["stages"]}
    for diagnostic in diagnostics:
        linked = [by_id[gate] for gate in diagnostic["supports"] if any(
            row["campaign_id"] == diagnostic["id"]
            for row in by_id[gate].get("progress_evidence", []))]
        assert linked
        assert all(stage["status"] != "ACCEPTED" and not stage["accepted_evidence"]
                   for stage in linked)
    theorem_join = next(row for row in diagnostics
                        if row["id"] == "OPH-FINITE-THEOREM-JOINS-2026-10-28-v1")
    review = by_id["CP-13"]
    assert review["status"] == "IN_PROGRESS" and not review["accepted_evidence"]
    paths = {row["path"] for row in review["progress_evidence"]
             if row["campaign_id"] == theorem_join["id"]}
    assert {"model/clock_defect_energy/defect_history_results.json",
            "model/clock_defect_energy/mode_energy_results.json"} <= paths


@pytest.mark.parametrize("damage", ["unknown_scope", "contract_scope", "gate_promotion",
                                  "changed_family", "unfrozen_boundary"])
def test_diagnostic_contract_cannot_change_scope_even_with_updated_hash(fixture_root, damage):
    data = json.loads((fixture_root / program.PROGRAM).read_text())
    diagnostic = next(row for row in data["component_campaigns"]
                      if row["scope"] == "FINITE_COMPONENT_DIAGNOSTICS")
    path = fixture_root / diagnostic["contract"]["path"]
    contract = json.loads(path.read_text())
    if damage == "unknown_scope":
        diagnostic["scope"] = contract["scope"] = "PHYSICAL_PROGRAM"
    elif damage == "contract_scope":
        contract["scope"] = "FINITE_COMPONENT_FAMILY"
    elif damage == "gate_promotion":
        contract["integrated_gate_accepted"] = True
    elif damage == "changed_family":
        contract["native_family_changed"] = True
    else:
        contract["boundary_status"] = "OPEN"
    path.write_text(json.dumps(contract))
    # Refresh every occurrence so this tests semantic scope, not stale custody.
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    diagnostic["contract"]["sha256"] = digest
    for stage in data["stages"]:
        for row in stage.get("progress_evidence", []):
            if row["path"] == diagnostic["contract"]["path"]:
                row["sha256"] = digest
    for row in data["next_wave"]["evidence"]:
        if row["path"] == diagnostic["contract"]["path"]:
            row["sha256"] = digest
    (fixture_root / program.PROGRAM).write_text(json.dumps(data))
    with pytest.raises(ValueError):
        program.load_program(fixture_root)

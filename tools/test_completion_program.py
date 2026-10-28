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
    for name in paths:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("Independent fixture evidence: " + name + "\n")
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
    assert all(row["status"] == "OPEN" for row in view["program"]["stages"])
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

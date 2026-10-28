#!/usr/bin/env python3
"""Validate and render the executable completion DAG; no scientific replay."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO_URL = "https://github.com/FloatingPragma/observer-patch-holography"
PROGRAM = "claims/completion_program.json"
HYPOTHESIS = "claims/executable_model_hypothesis.json"
OUTPUT = "docs/COMPLETION_PROGRAM.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def local_file(root: Path, name: str) -> Path:
    require(isinstance(name, str) and bool(name), "missing local artifact")
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()) and path.is_file(),
            f"missing or escaping artifact: {name}")
    return path


def bound_evidence(root: Path, rows: list) -> None:
    require(isinstance(rows, list) and bool(rows), "missing bound decision evidence")
    for row in rows:
        require(isinstance(row, dict), "untyped decision evidence")
        path = local_file(root, row.get("path"))
        require(hashlib.sha256(path.read_bytes()).hexdigest() == row.get("sha256"),
                f"decision evidence custody changed: {row.get('path')}")
        local_file(root, row.get("independent_check"))
        require(bool(row.get("scope")), "missing evidence scope")


def load_program(root: Path = ROOT, *, allow_unlinked: bool = False) -> tuple[dict, dict]:
    root = root.resolve()
    program = json.loads(local_file(root, PROGRAM).read_text())
    hypothesis = json.loads(local_file(root, HYPOTHESIS).read_text())
    require(program.get("schema") == "oph.completion_program.v1", "wrong program schema")
    require(hypothesis.get("schema") == "oph.executable_model_hypothesis.v1", "wrong hypothesis schema")
    require(program.get("hypothesis") == HYPOTHESIS, "wrong hypothesis parent")
    require(hypothesis.get("claim_class") == "working_hypothesis"
            and hypothesis.get("status") == "PROPOSED_NOT_ESTABLISHED",
            "hypothesis cannot be silently promoted to established physics")
    for flag in ("global_positive_mass_gap_claimed", "modular_flow_automatically_physical_time",
                 "all_particles_stable_defects_established", "einstein_premises_realized"):
        require(hypothesis.get(flag) is False, f"unsupported hypothesis promotion: {flag}")
    require(hypothesis.get("common_process_required") is True
            and hypothesis.get("comparison_count_added") == 0, "wrong joint comparison scope")
    local_file(root, hypothesis.get("contract"))
    obligations = hypothesis.get("obligations", [])
    obligation_ids = {row["id"] for row in obligations}
    require(len(obligations) == 5 and obligation_ids == {f"EMH-{i}" for i in range(1, 6)},
            "missing or duplicate hypothesis obligation")
    require(all(row.get("status") == "NOT_ESTABLISHED" for row in obligations),
            "physical obligation promotion needs a separately reviewed contract revision")
    stages = program.get("stages", [])
    ids = [stage["id"] for stage in stages]
    require(len(ids) == len(set(ids)) and set(ids) == {f"CP-{i:02d}" for i in range(14)},
            "missing or duplicate completion gate")
    by_id = {stage["id"]: stage for stage in stages}
    claims = json.loads(local_file(root, "claims/claim_registry.yaml").read_text())["claims"]
    claim_ids = {row["claim_id"] for row in claims}
    foundations = program.get("foundations", [])
    require(bool(foundations) and len({row["id"] for row in foundations}) == len(foundations),
            "missing or duplicate retained foundations")
    for foundation in foundations:
        require(foundation.get("claim_id") in claim_ids, "foundation lacks a canonical claim")
        require(set(foundation.get("supports", [])) <= set(ids)
                and bool(foundation.get("achieved")) and bool(foundation.get("boundary")),
                "unscoped retained foundation")
        for path in foundation["evidence"]:
            local_file(root, path)
    issues = []
    for stage in stages:
        deps = stage.get("requires", [])
        require(isinstance(deps, list) and len(deps) == len(set(deps))
                and set(deps) <= set(ids) and stage["id"] not in deps,
                f"invalid dependencies: {stage['id']}")
        require(stage.get("status") in {"OPEN", "IN_PROGRESS", "ACCEPTED", "OBSTRUCTED"},
                f"invalid gate status: {stage['id']}")
        require(bool(stage.get("acceptance")) and bool(stage.get("obstruction"))
                and bool(stage.get("deliverable")) and bool(stage.get("evidence")),
                f"incomplete acceptance/obstruction contract: {stage['id']}")
        require(set(stage.get("obligations", [])) <= obligation_ids
                and bool(stage.get("obligations")), "unmapped hypothesis obligation")
        for path in stage["evidence"]:
            local_file(root, path)
        if stage["status"] == "ACCEPTED":
            bound_evidence(root, stage.get("accepted_evidence"))
            require(all(by_id[dep]["status"] == "ACCEPTED" for dep in deps),
                    f"gate accepted before predecessors: {stage['id']}")
        issue = stage.get("issue")
        if issue is not None:
            require(type(issue) is int and issue > 0, "invalid GitHub issue")
            issues.append(issue)
        else:
            require(allow_unlinked, f"unlinked completion issue: {stage['id']}")
    require(len(issues) == len(set(issues)), "duplicate completion issue ownership")
    visiting, visited = set(), set()

    def visit(id):
        require(id not in visiting, f"completion dependency cycle at {id}")
        if id in visited:
            return
        visiting.add(id)
        for dep in by_id[id]["requires"]:
            visit(dep)
        visiting.remove(id)
        visited.add(id)

    for id in ids:
        visit(id)
    umbrella = program.get("umbrella_issue")
    require((allow_unlinked and umbrella is None)
            or (type(umbrella) is int and umbrella > 0 and umbrella not in issues),
            "invalid umbrella issue")
    frozen = program.get("boundary_status") == "FROZEN"
    require(program.get("boundary_status") in {"NOT_FROZEN", "FROZEN"}, "invalid boundary state")
    require(frozen == (by_id["CP-00"]["status"] == "ACCEPTED"),
            "frozen class requires an accepted boundary contract")
    require(program.get("decision") in {"NOT_DECIDED", "FALSIFIED_PROGRAM", "STOPPED_RESOURCE", "COMPLETED_DOMAIN"},
            "invalid program decision")
    falsified = program.get("decision") == "FALSIFIED_PROGRAM"
    require(program.get("whole_program_falsified") is falsified, "inconsistent whole-program verdict")
    if program["decision"] != "NOT_DECIDED":
        bound_evidence(root, program.get("decision_evidence"))
    if falsified:
        require(frozen, "whole-program falsification requires a frozen quantified class")
        require(program.get("coverage") == "ALL_FROZEN_REALIZATIONS"
                and program.get("decision_kind") in {"INDISPENSABLE_CONTRADICTION", "UNIVERSAL_REALIZATION_NO_GO", "UNAVOIDABLE_EMPIRICAL_EXCLUSION"},
                "missing whole-program coverage or valid falsifier")
    if program["decision"] == "COMPLETED_DOMAIN":
        require(frozen and all(stage["status"] == "ACCEPTED" for stage in stages),
                "completion requires every accepted gate")
    return program, hypothesis


def issue_link(number: int | None, label: str | None = None) -> str:
    return f"[{label or '#' + str(number)}]({REPO_URL}/issues/{number})" if number else "Issue pending"


def readiness(stage: dict, program: dict) -> str:
    by_id = {row["id"]: row for row in program["stages"]}
    if stage["status"] in {"ACCEPTED", "OBSTRUCTED"}:
        return stage["status"]
    missing = [dep for dep in stage["requires"] if by_id[dep]["status"] != "ACCEPTED"]
    return "Awaiting " + ", ".join(missing) if missing else "Ready to start"


def mermaid(program: dict) -> str:
    lines = ["```mermaid", "flowchart TD"]
    for stage in program["stages"]:
        label = stage["id"] + ": " + stage["title"]
        lines.append(f"    {stage['id'].replace('-', '')}[{json.dumps(label)}]")
    for stage in program["stages"]:
        for dep in stage["requires"]:
            lines.append(f"    {dep.replace('-', '')} --> {stage['id'].replace('-', '')}")
    for stage in program["stages"]:
        if stage.get("issue"):
            lines.append(f"    click {stage['id'].replace('-', '')} \"{REPO_URL}/issues/{stage['issue']}\"")
    return "\n".join([*lines, "```"])


def render(program: dict, hypothesis: dict) -> str:
    lines = ["# Executable OPH completion program", "",
             f"Generated by `tools/build_completion_program.py` from [`{PROGRAM}`](../{PROGRAM}).",
             "", f"Coordination: {issue_link(program['umbrella_issue'])}. "
             "Physical target: [the executable hypothesis](OPH_EXECUTABLE_HYPOTHESIS.md).", "",
             "**One native repair/memory process must supply the observer clocks, "
             "defect/maintenance-energy mechanism, interacting field/particle description "
             "and realized Einstein premises together.** Component benchmarks retain their "
             "own contracts; they cannot be pooled into that realization.", "",
             f"Program: `{program['program_id']}`. Boundary: `{program['boundary_status']}`. "
             f"Decision: `{program['decision']}`. No integrated gate is accepted by opening an issue.", "",
             "## Existing foundations to extend", "",
             "The integrated gates remain open, but the program starts from retained "
             "positive constructions. The rows below link to existing canonical claims; "
             "they are reuse/navigation bindings, not new proof-discharge assertions. "
             "Preserve these capabilities while constructing the missing joins.", "",
             "| Foundation | Achieved component | Remaining join | Builds into |",
             "| --- | --- | --- | --- |"]
    for foundation in program["foundations"]:
        lines.append(f"| {foundation['title']} | {foundation['achieved']} | {foundation['boundary']} | {', '.join(foundation['supports'])} |")
    lines += ["", "In particular, retain the modular-clock candidate and the "
              "conditional 3+1 field construction based on the native response "
              "structure. The slow band's rank three and the native screen's "
              "spectral dimension two are distinct facts; the supplied radial "
              "extension remains part of that field construction's contract.", "",
             "## Next work and dependencies", "",
             "Start with CP-00. It must identify the law/resource/interpretation class "
             "whose complete exclusion would stop this program. Component investigations "
             "can continue provisionally; their results cannot discharge a dependent gate "
             "without the predecessor mappings. CP-12 audits every stage in parallel "
             "once the program boundary is fixed.", "", mermaid(program), "",
             "| Gate | Research task | Dependency readiness | Issue |", "| --- | --- | --- | --- |"]
    for stage in program["stages"]:
        lines.append(f"| {stage['id']} | {stage['title']} | {readiness(stage, program)} | {issue_link(stage['issue'])} |")
    lines += ["", "The arrows are research acceptance dependencies, not a new theorem "
              "entailment graph. The existing [claim graph](../claims/dependency_graph.json) "
              "continues to record theorem/claim dependencies.", "", "## Decisions and stopping", "",
              "| Verdict | Evidence needed | Consequence |", "| --- | --- | --- |",
              "| Candidate rejection | A replayable failure of the specified law/state/access/readout contract | Stop that candidate and retain its result |",
              "| Class obstruction | A proof covering every realization in an explicitly quantified subclass | Stop that class; list remaining admitted classes |",
              "| Whole-program falsification | Verified contradiction in indispensable premises, exclusion of every frozen realization of a required join, or an unavoidable independently excluded prediction | Stop the frozen program; any premise/class change starts a new version |",
              "| Inconclusive | Open maps, incomplete coverage or an unsuccessful finite search | No confirmation and no impossibility claim |",
              "| Resource/value stop | A declared budget or return threshold is reached | A research stop, not scientific falsification |", "",
              "The boundary is currently **not frozen**, so no whole-program exclusion "
              "is claimed. Broad closure alone does not guarantee this specific C/E "
              "candidate, and an indefinitely enlarged escape class cannot be declared "
              "exhausted. Freezing CP-00 makes the tested formulation and its possible "
              "rejection explicit. Changing it never erases an earlier verdict.", "",
              "Every obstruction records target revision, mandatory premises, quantified "
              "class, evidence and independent check, coverage, surviving escape routes "
              "and the exact consequence. Issue closure alone is not scientific acceptance. "
              "The registry validator checks structure and artifact custody, not the truth "
              "of a mathematical proof or a physical interpretation.", "",
              "## Gate contracts", ""]
    for stage in program["stages"]:
        lines += [f"### {stage['id']}: {stage['title']}", "",
                  f"Owner: {issue_link(stage['issue'])}. Status: `{stage['status']}`. "
                  f"Requires: {', '.join(stage['requires']) or 'none'}. "
                  f"Hypothesis obligations: {', '.join(stage['obligations'])}.", "",
                  stage["deliverable"], "", "Acceptance:", ""]
        lines += [f"- {row}" for row in stage["acceptance"]]
        lines += ["", "Obstruction/stopping test: " + stage["obstruction"], "", "Retained starting evidence:", ""]
        lines += [f"- [{Path(path).name}](../{path})" for path in stage["evidence"]]
        lines += ["", "Earlier related work: " + ", ".join(issue_link(number) for number in stage["legacy_issues"]) + ".", ""]
    lines += ["## Working organization", "",
              program["legacy_issue_policy"] + " Live observer code, "
              "tests and retained results stay in `model/`; the existing particle code "
              "and theorem library keep their scientific ownership. Do not move or "
              "rewrite frozen packets to make the navigation look simpler.", "",
              "New work identifies its CP gate, fixed contract, nearest existing result, "
              "actual model/theorem map, falsifier and preservation checks. Prefer the "
              "smallest exact counterexample or native witness before larger runs. "
              "A supplied benchmark is evidence only for its stated conditional task. "
              "Each accepted gate needs a retained, source-bound result and the explicit "
              "join to its accepted predecessors.", "",
              "The generated [postdiction ledger](POSTDICTION_LEDGER.md) carries the "
              "hypothesis and program custody without adding a physical comparison. "
              "The [falsification program](OPH_FALSIFICATION_PROGRAM.md) retains its "
              "claim/branch-specific rejection rules and the stronger scope required "
              "for whole-program decisions.", ""]
    return "\n".join(lines)


def render_issue(stage: dict, program: dict) -> str:
    by_id = {row["id"]: row for row in program["stages"]}
    lines = ["## Goal", "", stage["deliverable"], "",
             f"Completion owner: {issue_link(program['umbrella_issue'])}. "
             f"Gate `{stage['id']}` in program `{program['program_id']}`.", "",
             f"[Common hypothesis]({REPO_URL}/blob/main/docs/OPH_EXECUTABLE_HYPOTHESIS.md) · "
             f"[Dependency DAG and gate contracts]({REPO_URL}/blob/main/{OUTPUT})", "",
             "## Dependencies", ""]
    lines += ([f"- [ ] {issue_link(by_id[dep]['issue'])}: {by_id[dep]['title']}" for dep in stage["requires"]]
              or ["No predecessor. This issue freezes the program boundary before integrated acceptance."])
    lines += ["", "Existing component work can proceed provisionally; integrated gate "
              "acceptance requires the actual maps to accepted predecessors. The same "
              "law, state/reference, observer access, clocks, energy and refinement "
              "domain must compose. A similarly named benchmark is not that map.", "",
              "## Acceptance", ""]
    lines += [f"- [ ] {row}" for row in stage["acceptance"]]
    lines += ["", "## Obstruction and stopping", "", stage["obstruction"], "",
              "Route every exclusion to CP-12 with target version, mandatory premises, "
              "quantified class, evidence/independent checks, coverage, surviving "
              "escape routes and consequence. Candidate failure, class no-go, theorem "
              "error, empirical rejection, inconclusive search and resource stop remain "
              "distinct. A whole-program verdict requires coverage of every frozen "
              "admissible realization; a changed class is a new program version.", "",
              "## Retained starting evidence", ""]
    lines += [f"- [{path}]({REPO_URL}/blob/main/{path})" for path in stage["evidence"]]
    lines += ["", "Related earlier work: " + ", ".join(issue_link(number) for number in stage["legacy_issues"]) + ".", "",
              "Opening or closing this issue does not establish a physical result. "
              "Acceptance needs the source-bound evidence, independent checks and "
              "predecessor joins recorded in the program registry. No retained failure "
              "or empirical exposure is reset by this reorganization.", ""]
    return "\n".join(lines)


def ledger_view(root: Path = ROOT) -> dict:
    program, hypothesis = load_program(root)
    parents = {}
    for name in (PROGRAM, HYPOTHESIS, hypothesis["contract"], "claims/claim_registry.yaml"):
        parents[name] = hashlib.sha256(local_file(root, name).read_bytes()).hexdigest()
    claims = {row["claim_id"]: row for row in json.loads(local_file(root, "claims/claim_registry.yaml").read_text())["claims"]}
    foundations = []
    for foundation in program["foundations"]:
        row = dict(foundation)
        row["canonical_claim_status"] = claims[row["claim_id"]]["status"]
        row["source_sha256"] = {name: hashlib.sha256(local_file(root, name).read_bytes()).hexdigest()
                                for name in row["evidence"]}
        foundations.append(row)
    return dict(verification="Program structure and source custody only; not scientific gate acceptance or a new physical comparison",
                parents=parents, hypothesis=hypothesis, program=program, retained_foundations=foundations)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--allow-unlinked", action="store_true", help="local drafting only")
    args = parser.parse_args()
    program, hypothesis = load_program(allow_unlinked=args.allow_unlinked)
    expected = render(program, hypothesis)
    if args.check:
        require((ROOT / OUTPUT).read_text() == expected, "completion program view is stale")
        print("Completion program DAG and issue bindings verified; no physical gate promoted")
    else:
        (ROOT / OUTPUT).write_text(expected)
        print(OUTPUT)


if __name__ == "__main__":
    main()

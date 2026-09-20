#!/usr/bin/env python3
"""Independent verifier for the Stage-2 order-sensitive source inventory.

This module deliberately does not import the producer.  It checks custody,
classification semantics, the one-packet qualification rule, decisive source
facts, and the fixture-import firewall independently.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_INVENTORY = HERE / "manifests" / "source_current_order_sensitive_inventory.json"
PRODUCER = HERE / "source_current_order_sensitive_inventory.py"

SCHEMA = "oph.source_current_order_sensitive_inventory.v1"
VERDICT = "SOURCE_CURRENT_ORDER_SENSITIVE_OBJECT_NOT_PRESENT"

AUDITED_DIRECTORIES = (
    "code/source_feedback_transport",
    "code/source_scalar_instruments",
    "code/source_routing",
    "code/a5_closure",
    "code/consensus",
    "code/refinement",
    "code/angular_sprint",
    "Lean/Dynamics",
    "Lean/Time",
    "Lean/Screen",
    "Lean/Geometry",
    "Lean/InformationProjection",
    "Lean/Variational",
)

QUALIFICATION_FIELDS = (
    "source_native",
    "port_indexed",
    "reversible",
    "both_composition_orders_recorded",
    "raw_histories_serialized",
    "target_free",
    "same_twelve_port_carrier",
    "refinement_provenance_present",
    "twelve_reversible_perturbation_families",
)

EXPECTED_CLASSIFICATIONS = {
    "registered_adjacency_response": "COMMUTATIVE_ONLY",
    "source_seam_path_tomography": "IRREVERSIBLE_ONLY",
    "source_feedback_transport": "IRREVERSIBLE_ONLY",
    "source_scalar_sequential_instrument": "IRREVERSIBLE_ONLY",
    "record_counting_repair_628": "IRREVERSIBLE_ONLY",
    "directed_seam_repair": "IRREVERSIBLE_ONLY",
    "consensus_settling_packets": "IRREVERSIBLE_ONLY",
    "carrier_and_discrete_refinement_scaffold": "STATIC_ONLY",
    "issue_566_bracket_search_space": "MISSING_ORDER_INFORMATION",
    "b14_oriented_face_bracket": "STATIC_ONLY",
    "static_response_grammar_and_poles": "STATIC_ONLY",
    "spin_deck_transport": "STATIC_ONLY",
    "declared_port_current_fixture": "DOWNSTREAM_CONTAMINATED",
    "downstream_kinetic_form_selector": "DOWNSTREAM_CONTAMINATED",
    "generic_transport_from_gluing": "UNRELATED",
    "lean_source_history_abstractions": "PARTIAL",
    "lean_conditional_a2_words": "COMMUTATIVE_ONLY",
    "lean_reversible_flow_models": "MISSING_ORDER_INFORMATION",
    "lean_declared_fin12_transport_refinement": "PARTIAL",
    "lean_irreversible_provenance_histories": "IRREVERSIBLE_ONLY",
    "lean_transport_and_current_boundaries": "STATIC_ONLY",
}


def expected_properties(*true_fields: str) -> dict[str, bool]:
    return {field: field in true_fields for field in QUALIFICATION_FIELDS}


EXPECTED_PROPERTIES = {
    "registered_adjacency_response": expected_properties(
        "source_native",
        "port_indexed",
        "reversible",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "source_seam_path_tomography": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "source_feedback_transport": expected_properties(
        "source_native", "port_indexed", "raw_histories_serialized", "target_free"
    ),
    "source_scalar_sequential_instrument": expected_properties("target_free"),
    "record_counting_repair_628": expected_properties(
        "source_native",
        "port_indexed",
        "raw_histories_serialized",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "directed_seam_repair": expected_properties(
        "source_native", "target_free", "same_twelve_port_carrier"
    ),
    "consensus_settling_packets": expected_properties(
        "source_native", "raw_histories_serialized", "target_free"
    ),
    "carrier_and_discrete_refinement_scaffold": expected_properties(
        "source_native",
        "port_indexed",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "issue_566_bracket_search_space": expected_properties(
        "port_indexed", "target_free", "same_twelve_port_carrier"
    ),
    "b14_oriented_face_bracket": expected_properties(
        "source_native", "port_indexed", "target_free", "same_twelve_port_carrier"
    ),
    "static_response_grammar_and_poles": expected_properties(
        "reversible", "target_free", "same_twelve_port_carrier"
    ),
    "spin_deck_transport": expected_properties(
        "source_native",
        "reversible",
        "target_free",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "declared_port_current_fixture": expected_properties(
        "port_indexed",
        "reversible",
        "same_twelve_port_carrier",
        "refinement_provenance_present",
    ),
    "downstream_kinetic_form_selector": expected_properties(
        "same_twelve_port_carrier"
    ),
    "generic_transport_from_gluing": expected_properties("target_free"),
    "lean_source_history_abstractions": expected_properties("target_free"),
    "lean_conditional_a2_words": expected_properties(
        "port_indexed", "reversible", "target_free", "same_twelve_port_carrier"
    ),
    "lean_reversible_flow_models": expected_properties("reversible", "target_free"),
    "lean_declared_fin12_transport_refinement": expected_properties(
        "port_indexed", "reversible", "target_free", "same_twelve_port_carrier"
    ),
    "lean_irreversible_provenance_histories": expected_properties(
        "source_native",
        "raw_histories_serialized",
        "target_free",
        "refinement_provenance_present",
    ),
    "lean_transport_and_current_boundaries": expected_properties(
        "source_native", "target_free", "same_twelve_port_carrier"
    ),
}


class VerificationError(ValueError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def strict_load(path: Path) -> dict[str, Any]:
    def object_pairs(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            check(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise VerificationError(f"non-finite JSON constant: {value}")

    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=object_pairs,
        parse_constant=reject_constant,
    )
    check(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def independently_list_audited_files(repo_root: Path) -> list[str]:
    paths: set[str] = set()
    for relative in AUDITED_DIRECTORIES:
        directory = repo_root / relative
        check(directory.is_dir(), f"missing audited directory: {relative}")
        for path in directory.rglob("*"):
            if not path.is_file():
                continue
            repo_relative = path.relative_to(repo_root)
            if "__pycache__" in repo_relative.parts or path.suffix == ".pyc":
                continue
            paths.add(repo_relative.as_posix())
    return sorted(paths)


def verify_audited_file_snapshot(snapshot: Any, repo_root: Path) -> set[str]:
    check(isinstance(snapshot, Mapping), "audited file snapshot missing")
    check(snapshot.get("directories") == list(AUDITED_DIRECTORIES), "audit directory drift")
    paths = snapshot.get("paths")
    check(
        isinstance(paths, list)
        and all(isinstance(path, str) and path for path in paths),
        "audited file path list",
    )
    check(paths == sorted(set(paths)), "audited file paths not canonical")
    check(snapshot.get("path_count") == len(paths), "audited file path count")
    check(snapshot.get("paths_sha256") == canonical_sha256(paths), "audited file list hash")
    current_paths = independently_list_audited_files(repo_root)
    check(paths == current_paths, "audited file snapshot drift")
    return set(paths)


def verify_pin_rows(rows: Any, repo_root: Path, label: str) -> set[str]:
    check(isinstance(rows, list) and rows, f"{label} pins missing")
    observed: set[str] = set()
    for row in rows:
        check(isinstance(row, Mapping), f"{label} pin row")
        relative = row.get("path")
        check(isinstance(relative, str) and relative, f"{label} pin path")
        path = Path(relative)
        check(not path.is_absolute() and ".." not in path.parts, f"unsafe pin path: {relative}")
        check(relative not in observed, f"duplicate pin path: {relative}")
        observed.add(relative)
        absolute = repo_root / path
        check(absolute.is_file(), f"missing pinned file: {relative}")
        check(row.get("bytes") == absolute.stat().st_size, f"byte pin drift: {relative}")
        check(row.get("sha256") == file_sha256(absolute), f"SHA pin drift: {relative}")
    return observed


def verify_import_firewall(producer_path: Path) -> None:
    tree = ast.parse(producer_path.read_text(encoding="utf-8"), filename=str(producer_path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {
        "__future__",
        "argparse",
        "collections",
        "hashlib",
        "json",
        "pathlib",
        "typing",
    }
    check(imports <= allowed, f"producer import firewall failed: {sorted(imports - allowed)}")
    check(
        "port_current_inner_certificate" not in imports,
        "conditional current implementation imported by producer",
    )


def independently_check_sources(repo_root: Path) -> dict[str, int]:
    projection = strict_load(
        repo_root / "code/a5_closure/manifests/source_current_capability_projection.json"
    )
    capability = strict_load(
        repo_root / "code/a5_closure/receipts/source_current_capability.receipt.json"
    )
    recurrence = projection.get("registered_recurrence", {})
    algebra = capability.get("audit", {}).get("response_word_algebra", {})
    rechart = capability.get("audit", {}).get("recharting_audit", {})
    check(recurrence.get("raw_runtime_history_serialized") is False, "raw-history scope")
    check(recurrence.get("source_count") == 12, "source count")
    check(algebra.get("exact_dimension") == 4, "response algebra dimension")
    check(algebra.get("commutative") is True, "response algebra commutativity")
    check(algebra.get("commutator_nonzero_count") == 0, "response commutator count")
    check(rechart.get("proper_recharting_count") == 60, "proper recharting count")
    check(
        rechart.get("proper_rechartings_in_response_word_algebra") == 1,
        "response/recharting intersection",
    )

    routing = strict_load(repo_root / "code/source_routing/runtime/path_tomography_receipt.json")
    routing_ops = {
        event.get("op")
        for episode in routing.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    check("pair_mean" in routing_ops, "path-tomography pair mean")
    check(
        routing.get("scope", {}).get("full_metric_neighbor_refinement") is False,
        "path-tomography refinement scope",
    )

    feedback = strict_load(repo_root / "code/source_feedback_transport/transport_receipt.json")
    feedback_ops = {
        event.get("op")
        for episode in feedback.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    check({"mean", "reset", "export"} <= feedback_ops, "feedback destructive operations")

    scalar = strict_load(
        repo_root / "code/source_scalar_instruments/sequential_instrument_receipt.json"
    )
    check(
        scalar.get("scope", {}).get("sampled_quantum_outcome_histories") is False,
        "scalar sampled-history scope",
    )
    check(
        scalar.get("scope", {}).get("source_selected_quantum_operations") is False,
        "scalar source-selection scope",
    )

    counting = strict_load(
        repo_root / "code/a5_closure/manifests/record_counting_mechanism_reference.json"
    )
    theorem = counting.get("a2_clauses", {}).get("presentation_invariance", {}).get(
        "theorem", ""
    )
    check("addition is commutative" in theorem, "record-counting order erasure")
    check(
        "append-only" in counting.get("mechanism", {}).get("event_grammar", ""),
        "record-counting append-only grammar",
    )

    directed = strict_load(
        repo_root / "code/a5_closure/manifests/directed_seam_repair_reference.json"
    )
    check(
        directed.get("verdict", {}).get("refinement_semigroup_naturality") == "open",
        "directed repair refinement boundary",
    )

    bracket1 = strict_load(
        repo_root
        / "code/a5_closure/issue_566_bracket_space_stage1/"
        "a5_alternating_bracket_space_stage1.receipt.json"
    )
    bracket2 = strict_load(
        repo_root
        / "code/a5_closure/issue_566_bracket_space_stage2/"
        "a5_jacobi_stage2.receipt.json"
    )
    check(not any(bracket1.get("later_gates", {}).values()), "stage1 later gates")
    check(
        bracket2.get("not_attained", {}).get("preferred_bracket_selected") is False,
        "stage2 bracket selection boundary",
    )

    b14 = strict_load(
        repo_root / "code/b14_jacobi/oriented_face_bracket_selector.certificate.json"
    )
    check(b14.get("jacobi_failure", {}).get("nonzero_count") == 240, "B14 Jacobi count")
    check(b14.get("jacobi_failure", {}).get("squared_norm") == 240, "B14 Jacobi norm")

    fixture = strict_load(
        repo_root / "code/a5_closure/manifests/port_current_response_reference.json"
    )
    check(fixture.get("construction_model") == "charged_double_triplet", "fixture typing")
    kinetic = strict_load(
        repo_root / "code/angular_sprint/runtime/kinetic_form_selection_receipt.json"
    )
    check(
        kinetic.get("physical_selection_boundary", {}).get("current_lift_source_selected")
        is False,
        "kinetic current source boundary",
    )

    source_history_text = (
        repo_root / "Lean/Variational/SourceToHamiltonianComposed.lean"
    ).read_text(encoding="utf-8")
    check(
        "theorem here claims that the source selects" in source_history_text,
        "Lean source-selection boundary",
    )
    response_text = (repo_root / "Lean/Screen/A5ResponseWordAlgebra.lean").read_text(
        encoding="utf-8"
    )
    check("commutative" in response_text.lower(), "Lean response-word boundary")
    a2_words = (
        repo_root / "Lean/ObserverPatchHolography/A2EndpointCommutator.lean"
    ).read_text(encoding="utf-8")
    check("FifteenCommute" in a2_words, "Lean conditional A2 word boundary")
    provenance = (
        repo_root
        / "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean"
    ).read_text(encoding="utf-8")
    check("two serializations of the diamond" in provenance, "Lean serialization boundary")

    return {
        "response_algebra_dimension": int(algebra["exact_dimension"]),
        "response_commutator_nonzero_count": int(algebra["commutator_nonzero_count"]),
        "proper_recharting_count": int(rechart["proper_recharting_count"]),
        "b14_jacobi_nonzero_count": int(b14["jacobi_failure"]["nonzero_count"]),
    }


def verify(inventory_path: Path, repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    inventory = strict_load(inventory_path)
    check(inventory.get("schema") == SCHEMA, "inventory schema")
    body = {key: value for key, value in inventory.items() if key != "inventory_sha256"}
    check(inventory.get("inventory_sha256") == canonical_sha256(body), "inventory hash")
    check(inventory.get("verdict") == VERDICT, "inventory verdict")
    check(inventory.get("positive_stages_authorized") is False, "positive-stage stop")
    check(
        inventory.get("physical_current_source_bridge_attained") is False,
        "false physical source gate",
    )
    check(
        inventory.get("fixture_used_as_reconstruction_oracle") is False,
        "fixture oracle firewall",
    )
    check(inventory.get("target_labels_used") is False, "target-label firewall")
    check(
        inventory.get("qualification_fields") == list(QUALIFICATION_FIELDS),
        "qualification field schema",
    )

    snapshotted_paths = verify_audited_file_snapshot(
        inventory.get("audited_file_snapshot"), repo_root
    )

    source_paths = verify_pin_rows(inventory.get("source_pins"), repo_root, "source")
    implementation_paths = verify_pin_rows(
        inventory.get("implementation_pins"), repo_root, "implementation"
    )
    check(
        implementation_paths
        == {
            "code/a5_closure/source_current_order_sensitive_inventory.py",
            "code/a5_closure/verify_source_current_order_sensitive_inventory.py",
        },
        "implementation pin set",
    )
    check(
        {
            path
            for path in source_paths | implementation_paths
            if any(path == root or path.startswith(root + "/") for root in AUDITED_DIRECTORIES)
        }
        <= snapshotted_paths,
        "pinned audited file missing from snapshot",
    )
    verify_import_firewall(repo_root / "code/a5_closure/source_current_order_sensitive_inventory.py")

    candidates = inventory.get("candidates")
    check(isinstance(candidates, list), "candidate list")
    candidate_map: dict[str, Mapping[str, Any]] = {}
    for row in candidates:
        check(isinstance(row, Mapping), "candidate row")
        candidate_id = row.get("candidate_id")
        check(isinstance(candidate_id, str), "candidate id")
        check(candidate_id not in candidate_map, f"duplicate candidate: {candidate_id}")
        candidate_map[candidate_id] = row
    check(set(candidate_map) == set(EXPECTED_CLASSIFICATIONS), "candidate coverage")

    qualifying: list[str] = []
    for candidate_id, expected_classification in EXPECTED_CLASSIFICATIONS.items():
        row = candidate_map[candidate_id]
        check(row.get("classification") == expected_classification, f"classification: {candidate_id}")
        props = row.get("properties")
        check(props == EXPECTED_PROPERTIES[candidate_id], f"properties: {candidate_id}")
        failures = [field for field in QUALIFICATION_FIELDS if not props[field]]
        check(row.get("qualification_failures") == failures, f"failure list: {candidate_id}")
        check(
            isinstance(row.get("rejection_reason"), str) and row.get("rejection_reason"),
            f"rejection reason: {candidate_id}",
        )
        row_paths = row.get("source_paths")
        check(isinstance(row_paths, list) and row_paths, f"source paths: {candidate_id}")
        check(set(row_paths) <= source_paths, f"unpinned source path: {candidate_id}")
        if all(props[field] for field in QUALIFICATION_FIELDS):
            qualifying.append(candidate_id)
    check(not qualifying, f"unexpected qualifying candidates: {qualifying}")

    counts = Counter(EXPECTED_CLASSIFICATIONS.values())
    summary = inventory.get("summary", {})
    check(summary.get("candidate_count") == len(EXPECTED_CLASSIFICATIONS), "candidate count")
    check(summary.get("qualifying_candidate_count") == 0, "qualifying count")
    check(summary.get("classification_counts") == dict(sorted(counts.items())), "class counts")
    strongest = summary.get("strongest_registered_source_packet", {})
    check(strongest.get("candidate_id") == "registered_adjacency_response", "strongest packet")
    check(strongest.get("response_word_dimension") == 4, "reported response dimension")
    check(strongest.get("commutator_nonzero_count") == 0, "reported commutator count")
    check(strongest.get("proper_rechartings_total") == 60, "reported rechartings")
    check(strongest.get("proper_rechartings_in_response_words") == 1, "reported intersection")
    check(
        summary.get("richest_logged_near_candidate", {}).get("candidate_id")
        == "record_counting_repair_628",
        "richest near-candidate",
    )
    check("different primitives" in summary.get("cross_packet_noncomposition", ""), "no compositing rule")

    missing = inventory.get("minimal_missing_fields")
    check(isinstance(missing, list), "minimal missing fields")
    check(
        [row.get("field") for row in missing]
        == [
            "perturbation_families",
            "ordered_mixed_histories",
            "raw_history_custody",
            "common_source_binding",
            "reversibility_witness",
            "source_firewall",
        ],
        "minimal missing-field set",
    )
    boundary = inventory.get("claim_boundary", {})
    required_false_promotions = (
        "bounded_repository_inventory_only",
        "does_not_assert_a_no_go_for_future_source_data",
        "does_not_weaken_the_abstract_forced_lie_type_theorem",
        "does_not_reject_the_conditional_port_current_fixture",
        "does_not_make_a_physical_current_claim",
        "stages_3_through_6_not_run",
    )
    check(all(boundary.get(key) is True for key in required_false_promotions), "claim boundary")

    independent = independently_check_sources(repo_root)
    check(independent["response_algebra_dimension"] == 4, "independent dimension")
    return {
        "verdict": VERDICT,
        "candidate_count": len(candidate_map),
        "qualifying_candidate_count": len(qualifying),
        "classification_counts": dict(sorted(counts.items())),
        "response_algebra_dimension": independent["response_algebra_dimension"],
        "proper_recharting_count": independent["proper_recharting_count"],
        "b14_jacobi_nonzero_count": independent["b14_jacobi_nonzero_count"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()
    print(json.dumps(verify(args.inventory, args.repo_root), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

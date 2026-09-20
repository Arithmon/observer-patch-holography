#!/usr/bin/env python3
"""Stage-2 inventory for the source-current order-sensitive object.

This producer audits existing committed source surfaces.  It does not build a
new perturbation model and it never imports the conditional PORT-CURRENT-INNER
implementation.  A candidate qualifies only if one and the same source packet
contains twelve port-indexed reversible perturbation families, both mixed
composition orders, raw serialized histories, target-free provenance, and the
registered twelve-port refinement ancestry.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
INVENTORY_PATH = HERE / "manifests" / "source_current_order_sensitive_inventory.json"
INDEPENDENT_VERIFIER = HERE / "verify_source_current_order_sensitive_inventory.py"

SCHEMA = "oph.source_current_order_sensitive_inventory.v1"
VERDICT = "SOURCE_CURRENT_ORDER_SENSITIVE_OBJECT_NOT_PRESENT"
UPSTREAM_MAIN_SHA = "a197baac23c0dbec3770649f0f488d6df9a41541"

ALLOWED_CLASSIFICATIONS = {
    "QUALIFIES",
    "PARTIAL",
    "STATIC_ONLY",
    "COMMUTATIVE_ONLY",
    "IRREVERSIBLE_ONLY",
    "DOWNSTREAM_CONTAMINATED",
    "MISSING_ORDER_INFORMATION",
    "UNRELATED",
}

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

AUDITED_SOURCE_PATHS = (
    "code/a5_closure/manifests/echosahedral_federation_reference.json",
    "code/a5_closure/manifests/charged_response_semantic_artifact.json",
    "code/a5_closure/manifests/source_current_capability_projection.json",
    "code/a5_closure/receipts/source_current_capability.receipt.json",
    "code/source_routing/README.md",
    "code/source_routing/runtime/path_tomography_receipt.json",
    "code/source_feedback_transport/transport_receipt.json",
    "code/source_scalar_instruments/sequential_instrument_receipt.json",
    "code/a5_closure/manifests/record_counting_mechanism_reference.json",
    "code/a5_closure/manifests/directed_seam_repair_reference.json",
    "code/consensus/manifests/compiled_lattice_settling_reference.json",
    "code/refinement/runtime/discrete_refinement_theorem_receipt.json",
    "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
    "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json",
    "code/b14_jacobi/oriented_face_bracket_selector.certificate.json",
    "code/a5_closure/manifests/response_grammar_completeness_reference.json",
    "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json",
    "code/a5_closure/manifests/charged_response_pole_residue_artifact.json",
    "code/a5_closure/manifests/spin_statistics_semantic_artifact.json",
    "code/a5_closure/manifests/port_current_response_reference.json",
    "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
    "code/angular_sprint/runtime/kinetic_form_selection_receipt.json",
    "code/consensus/transport_from_gluing.py",
    "Lean/InformationProjection/SourceHistoryPacket.lean",
    "Lean/Variational/SourceToHamiltonianComposed.lean",
    "Lean/QFT/SourceOperatorGeneration.lean",
    "Lean/QFT/SourceCorrelationCapstone.lean",
    "Lean/QFT/SourceHistoryGNSDynamics.lean",
    "Lean/QFT/SourceDerivedEventPrecedence.lean",
    "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean",
    "Lean/ObserverPatchHolography/A2EndpointCommutator.lean",
    "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean",
    "Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean",
    "Lean/Geometry/SourceFeedbackTransport.lean",
    "Lean/Geometry/SourceSeamPathTomography.lean",
    "Lean/Geometry/CarrierDynamicsCompatibility.lean",
    "Lean/Geometry/WorldlineHopTransport.lean",
    "Lean/Screen/A2HolonomyBridge.lean",
    "Lean/Screen/A5ResponseWordAlgebra.lean",
    "Lean/Screen/CarrierEvolutionFlow.lean",
    "Lean/Screen/OrientedFaceBracketSelector.lean",
    "Lean/Screen/PrimitivePortTranslationBridge.lean",
    "Lean/Screen/RepairWordCarrierReadout.lean",
)


class InventoryError(ValueError):
    """Typed fail-closed inventory error."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise InventoryError(code, message)


def load_json(relative: str) -> dict[str, Any]:
    path = REPO_ROOT / relative
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), "JSON_OBJECT", f"{relative} is not an object")
    return value


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def canonical_sha256(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def source_pin(relative: str) -> dict[str, Any]:
    path = REPO_ROOT / relative
    require(path.is_file(), "SOURCE_PATH", f"missing audited source: {relative}")
    return {
        "path": relative,
        "bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def properties(**values: bool) -> dict[str, bool]:
    require(
        set(values) == set(QUALIFICATION_FIELDS),
        "PROPERTY_SCHEMA",
        "candidate property set drift",
    )
    require(
        all(isinstance(value, bool) for value in values.values()),
        "PROPERTY_TYPE",
        "candidate properties must be booleans",
    )
    return {key: values[key] for key in QUALIFICATION_FIELDS}


def candidate(
    candidate_id: str,
    classification: str,
    source_paths: Sequence[str],
    source_ancestry: str,
    history_custody: str,
    observed_properties: Mapping[str, bool],
    rejection_reason: str,
    evidence_facts: Sequence[str],
) -> dict[str, Any]:
    require(
        classification in ALLOWED_CLASSIFICATIONS,
        "CLASSIFICATION",
        f"invalid classification for {candidate_id}",
    )
    props = dict(observed_properties)
    qualifies = all(props[field] for field in QUALIFICATION_FIELDS)
    require(
        (classification == "QUALIFIES") == qualifies,
        "QUALIFICATION",
        f"classification/property mismatch for {candidate_id}",
    )
    return {
        "candidate_id": candidate_id,
        "classification": classification,
        "source_paths": list(source_paths),
        "source_ancestry": source_ancestry,
        "history_custody": history_custody,
        "properties": props,
        "qualification_failures": [
            field for field in QUALIFICATION_FIELDS if not props[field]
        ],
        "rejection_reason": rejection_reason,
        "evidence_facts": list(evidence_facts),
    }


def audit_source_semantics() -> dict[str, Any]:
    """Check the decisive source facts before any classification is emitted."""

    carrier = load_json("code/a5_closure/manifests/echosahedral_federation_reference.json")
    ports = carrier.get("carrier", {}).get("ports", [])
    require(len(ports) == 12 and len(set(ports)) == 12, "CARRIER", "twelve ports required")

    projection = load_json(
        "code/a5_closure/manifests/source_current_capability_projection.json"
    )
    capability = load_json(
        "code/a5_closure/receipts/source_current_capability.receipt.json"
    )
    recurrence = projection.get("registered_recurrence", {})
    algebra = capability.get("audit", {}).get("response_word_algebra", {})
    rechart = capability.get("audit", {}).get("recharting_audit", {})
    require(
        projection.get("schema") == "oph.source_current_capability_projection.v1"
        and recurrence.get("raw_runtime_history_serialized") is False
        and recurrence.get("source_count") == 12
        and recurrence.get("readback_count") == 12,
        "CAPABILITY_SCOPE",
        "registered recurrence scope drift",
    )
    require(
        algebra.get("exact_dimension") == 4
        and algebra.get("commutative") is True
        and algebra.get("commutator_nonzero_count") == 0
        and rechart.get("proper_recharting_count") == 60
        and rechart.get("proper_rechartings_in_response_word_algebra") == 1
        and capability.get("physical_current_source_bridge_attained") is False,
        "CAPABILITY_OBSTRUCTION",
        "bounded response-algebra obstruction drift",
    )

    path_tomography = load_json(
        "code/source_routing/runtime/path_tomography_receipt.json"
    )
    path_scope = path_tomography.get("scope", {})
    path_ops = {
        event.get("op")
        for episode in path_tomography.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    routing_readme = (REPO_ROOT / "code/source_routing/README.md").read_text(
        encoding="utf-8"
    )
    require(
        path_tomography.get("schema") == "oph.source_routing.path_tomography.v1"
        and path_scope.get("all_intermediate_events_retained") is True
        and path_scope.get("canonical_scalar_pair_means_only") is True
        and path_scope.get("full_metric_neighbor_refinement") is False
        and "pair_mean" in path_ops
        and "destructive readout construction" in routing_readme,
        "PATH_TOMOGRAPHY",
        "path tomography typing drift",
    )

    feedback = load_json("code/source_feedback_transport/transport_receipt.json")
    feedback_ops = {
        event.get("op")
        for episode in feedback.get("episodes", [])
        for event in episode.get("events", [])
        if isinstance(event, Mapping)
    }
    require(
        feedback.get("schema") == "oph.source_feedback_transport.receipt.v1"
        and "SUPPLIED_COPY_RESET_LAW" in str(feedback.get("scope"))
        and {"mean", "reset", "export"} <= feedback_ops,
        "FEEDBACK_TRANSPORT",
        "feedback transport operation typing drift",
    )

    scalar = load_json(
        "code/source_scalar_instruments/sequential_instrument_receipt.json"
    )
    scalar_scope = scalar.get("scope", {})
    require(
        scalar.get("schema") == "oph.source_scalar.centered_sequential_instrument.v1"
        and scalar_scope.get("sampled_quantum_outcome_histories") is False
        and scalar_scope.get("source_selected_quantum_operations") is False
        and scalar_scope.get("physical_clock_or_observed_outcomes") is False
        and len(scalar.get("regional_gadget", {}).get("detector_sites", [])) == 32,
        "SCALAR_INSTRUMENT",
        "sequential instrument scope drift",
    )

    counting = load_json(
        "code/a5_closure/manifests/record_counting_mechanism_reference.json"
    )
    counting_mechanism = counting.get("mechanism", {})
    presentation = counting.get("a2_clauses", {}).get("presentation_invariance", {})
    require(
        counting.get("schema") == "oph.record_counting_mechanism_certificate.v2"
        and "append-only" in str(counting_mechanism.get("event_grammar"))
        and "transfer" in str(counting_mechanism.get("repair_move"))
        and "addition is commutative" in str(presentation.get("theorem"))
        and presentation.get("permutations_checked") == 99,
        "RECORD_COUNTING",
        "record-counting irreversibility/order-erasure drift",
    )

    directed = load_json(
        "code/a5_closure/manifests/directed_seam_repair_reference.json"
    )
    directed_verdict = directed.get("verdict", {})
    require(
        directed.get("schema") == "oph.directed_seam_repair_certificate.v1"
        and directed_verdict.get("protected_record_and_checkpoint_instrument") == "open"
        and directed_verdict.get("refinement_semigroup_naturality") == "open"
        and directed.get("integer_macro_bridge", {}).get("expectation_equals_pair_average")
        is True,
        "DIRECTED_REPAIR",
        "directed repair boundary drift",
    )

    consensus = load_json(
        "code/consensus/manifests/compiled_lattice_settling_reference.json"
    )
    consensus_scope = consensus.get("scope", {})
    require(
        consensus.get("schema") == "oph.compiled_lattice_settling_certificate.v2"
        and "finite software patch dynamics" in str(consensus_scope.get("realized"))
        and "physical-hardware attachment" in str(consensus_scope.get("open")),
        "CONSENSUS_SCOPE",
        "consensus claim boundary drift",
    )

    refinement = load_json(
        "code/refinement/runtime/discrete_refinement_theorem_receipt.json"
    )
    require(
        refinement.get("schema") == "oph.discrete_refinement_theorem_packet.v1"
        and refinement.get("certified_exit")
        == "THEOREM_PACKET_AND_MESH_SCAFFOLD_ATTAINED__SOURCE_NATIVE_PHYSICAL_BIREFINEMENT_OPEN",
        "REFINEMENT_SCOPE",
        "refinement scaffold boundary drift",
    )

    bracket1 = load_json(
        "code/a5_closure/issue_566_bracket_space_stage1/"
        "a5_alternating_bracket_space_stage1.receipt.json"
    )
    bracket2 = load_json(
        "code/a5_closure/issue_566_bracket_space_stage2/"
        "a5_jacobi_stage2.receipt.json"
    )
    require(
        bracket1.get("status") == "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY"
        and not any(bracket1.get("later_gates", {}).values())
        and bracket2.get("status")
        == "EXACT_JACOBI_SYSTEM_AND_FIXED_LINE_REDUCTION__FULL_CLASSIFICATION_OPEN"
        and bracket2.get("not_attained", {}).get("preferred_bracket_selected") is False,
        "BRACKET_SPACE",
        "issue-566 search-space boundary drift",
    )

    b14 = load_json(
        "code/b14_jacobi/oriented_face_bracket_selector.certificate.json"
    )
    require(
        b14.get("schema") == "oph.b14.oriented_face_bracket_selector.v1"
        and b14.get("jacobi_failure", {}).get("nonzero_count") == 240
        and b14.get("jacobi_failure", {}).get("squared_norm") == 240,
        "B14_CONTROL",
        "B14 Jacobi control drift",
    )

    grammar = load_json(
        "code/a5_closure/manifests/response_grammar_completeness_reference.json"
    )
    grammar_receipt = load_json(
        "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json"
    )
    require(
        grammar.get("schema") == "oph.response_grammar_completeness_certificate.v1"
        and grammar_receipt.get("schema")
        == "oph.response_grammar_completeness_receipt.v1"
        and grammar_receipt.get("verdict", {}).get("bounded_exit")
        == "independence_limited"
        and "continuum dynamics"
        in grammar_receipt.get("claim_boundary", {}).get("does_not_close", []),
        "RESPONSE_GRAMMAR",
        "response grammar boundary drift",
    )

    spin = load_json(
        "code/a5_closure/manifests/spin_statistics_semantic_artifact.json"
    )
    require(
        spin.get("schema") == "oph.spin_statistics_semantic_artifact.v1"
        and spin.get("physical_source_gate", {}).get("passed") is True
        and spin.get("physical_source_gate", {}).get(
            "continuum_spin_statistics_theorem"
        )
        is False
        and spin.get("physical_source_gate", {}).get(
            "laboratory_exchange_measurement"
        )
        is False,
        "SPIN_TRANSPORT",
        "spin/deck source gate drift",
    )

    fixture = load_json(
        "code/a5_closure/manifests/port_current_response_reference.json"
    )
    require(
        fixture.get("schema") == "oph.port_current_response_manifest.v5"
        and fixture.get("construction_model") == "charged_double_triplet"
        and "fixture" in str(
            fixture.get("response_declaration_contract", {}).get("status")
        ),
        "FIXTURE_BOUNDARY",
        "conditional port-current fixture typing drift",
    )

    kinetic = load_json(
        "code/angular_sprint/runtime/kinetic_form_selection_receipt.json"
    )
    require(
        kinetic.get("schema") == "oph.kinetic_form_selection_receipt.v1"
        and kinetic.get("physical_selection_boundary", {}).get(
            "current_lift_source_selected"
        )
        is False
        and "charged-double-triplet" in str(kinetic.get("current_fixture")),
        "KINETIC_FIXTURE",
        "kinetic-form fixture boundary drift",
    )

    required_text = {
        "Lean/InformationProjection/SourceHistoryPacket.lean": (
            "sourceHistoryPacket_receipt",
        ),
        "Lean/Variational/SourceToHamiltonianComposed.lean": (
            "theorem here claims that the source selects",
        ),
        "Lean/Geometry/SourceFeedbackTransport.lean": (
            "supplied operations",
        ),
        "Lean/Geometry/SourceSeamPathTomography.lean": (
            "not a selected routing",
        ),
        "Lean/Screen/A5ResponseWordAlgebra.lean": (
            "commutative",
        ),
        "Lean/Screen/OrientedFaceBracketSelector.lean": (
            "Jacobi",
        ),
        "Lean/ObserverPatchHolography/A2EndpointCommutator.lean": (
            "fifteen strict-index commutators",
        ),
        "Lean/Screen/RepairWordCarrierReadout.lean": (
            "does not emit that complete word-action packet",
        ),
        "Lean/Screen/CarrierEvolutionFlow.lean": (
            "one-parameter group",
        ),
        "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean": (
            "two serializations of the diamond",
        ),
    }
    for relative, needles in required_text.items():
        text = (REPO_ROOT / relative).read_text(encoding="utf-8")
        require(
            all(needle.lower() in text.lower() for needle in needles),
            "LEAN_BOUNDARY",
            f"Lean boundary drift: {relative}",
        )

    return {
        "port_count": len(ports),
        "response_algebra_dimension": algebra["exact_dimension"],
        "response_commutator_nonzero_count": algebra["commutator_nonzero_count"],
        "proper_recharting_count": rechart["proper_recharting_count"],
        "path_episode_count": len(path_tomography.get("episodes", [])),
        "feedback_event_count": feedback.get("summary", {}).get("events"),
        "record_presentation_permutations": presentation["permutations_checked"],
        "b14_jacobi_nonzero_count": b14["jacobi_failure"]["nonzero_count"],
    }


def build_candidates(facts: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return the adjudicated list.  No candidate data are combined."""

    rows = [
        candidate(
            "registered_adjacency_response",
            "COMMUTATIVE_ONLY",
            (
                "code/a5_closure/manifests/charged_response_semantic_artifact.json",
                "code/a5_closure/manifests/source_current_capability_projection.json",
                "code/a5_closure/receipts/source_current_capability.receipt.json",
            ),
            "Pinned twelve-port carrier impulse readback and its exact adjacency recurrence.",
            "STATIC_RECONSTRUCTION",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "The only registered reversible channel is one functional calculus in A: exact dimension 4, zero nonzero commutators, and no raw runtime history.",
            (
                f"response-word dimension = {facts['response_algebra_dimension']}",
                f"nonzero commutators = {facts['response_commutator_nonzero_count']}",
                "registered_recurrence.raw_runtime_history_serialized = false",
            ),
        ),
        candidate(
            "source_seam_path_tomography",
            "IRREVERSIBLE_ONLY",
            (
                "code/source_routing/runtime/path_tomography_receipt.json",
                "code/source_routing/README.md",
            ),
            "Captured W12 level-three support with declared scalar pair-mean calibration episodes.",
            "RAW_DESTRUCTIVE_EVENT_LOG",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It serializes exact histories, but pair averaging is destructive, paths cross many carriers, one declared schedule is used, and full routing/refinement is false.",
            (
                f"serialized episodes = {facts['path_episode_count']}",
                "all intermediate events retained",
                "README types the construction as destructive",
            ),
        ),
        candidate(
            "source_feedback_transport",
            "IRREVERSIBLE_ONLY",
            ("code/source_feedback_transport/transport_receipt.json",),
            "Classical feedback episodes over the captured W12 routing support.",
            "RAW_DESTRUCTIVE_EVENT_LOG",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Mean, reset, and export overwrite state under one fixed serial schedule; the logs are transport histories, not reversible current perturbations.",
            (f"serialized events = {facts['feedback_event_count']}", "reset and mean operations present"),
        ),
        candidate(
            "source_scalar_sequential_instrument",
            "IRREVERSIBLE_ONLY",
            ("code/source_scalar_instruments/sequential_instrument_receipt.json",),
            "Supplied 64-site centered quantum instrument and analytic cross-time brackets.",
            "STATIC_INSTRUMENT_WITHOUT_SAMPLED_HISTORIES",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The full instrument includes reset/readout, is 64-site rather than W12, and explicitly has no sampled outcome histories or source-selected operations.",
            ("sampled_quantum_outcome_histories = false", "source_selected_quantum_operations = false"),
        ),
        candidate(
            "record_counting_repair_628",
            "IRREVERSIBLE_ONLY",
            ("code/a5_closure/manifests/record_counting_mechanism_reference.json",),
            "Operational signed record events and strict-descent repair on the pinned twelve-port carrier.",
            "RAW_APPEND_ONLY_RECORD_LOG",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "This is the richest near-candidate, but its append-only/strict-descent dynamics are irreversible and presentation order is explicitly quotiented by commutative addition.",
            (
                "event grammar is append-only",
                f"presentation permutations checked = {facts['record_presentation_permutations']}",
                "three coherent refinement maps are pinned",
            ),
        ),
        candidate(
            "directed_seam_repair",
            "IRREVERSIBLE_ONLY",
            ("code/a5_closure/manifests/directed_seam_repair_reference.json",),
            "Expected pair-average seam repair and representative integer microhistories.",
            "REPRESENTATIVE_REPAIR_HISTORY",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Directed placements are seam repairs, not per-port reversible perturbations; protected-record instrumentation and refinement naturality remain open.",
            ("expectation equals pair average", "refinement_semigroup_naturality = open"),
        ),
        candidate(
            "consensus_settling_packets",
            "IRREVERSIBLE_ONLY",
            ("code/consensus/manifests/compiled_lattice_settling_reference.json",),
            "Compiler scatter/project settling trajectory and its pinned upstream implementation.",
            "REPRESENTATIVE_SETTLING_TRAJECTORY",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Scatter/project settling erases state and its local arities are not the twelve-port current carrier.",
            ("settling trajectory is serialized", "scope excludes a physical current construction"),
        ),
        candidate(
            "carrier_and_discrete_refinement_scaffold",
            "STATIC_ONLY",
            (
                "code/a5_closure/manifests/echosahedral_federation_reference.json",
                "code/refinement/runtime/discrete_refinement_theorem_receipt.json",
            ),
            "Exact twelve-port carrier, port maps, and abstract mesh/refinement theorems.",
            "STATIC_COMBINATORIAL_SCAFFOLD",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "The carrier/refinement packet supplies ancestry and commuting maps, but no perturbation, inverse, or ordered source history.",
            ("twelve carrier ports", "source-native physical birefinement remains open"),
        ),
        candidate(
            "issue_566_bracket_search_space",
            "MISSING_ORDER_INFORMATION",
            (
                "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
                "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_stage2.receipt.json",
            ),
            "Target-free equivariant bracket coefficient space and exact Jacobi equations.",
            "STATIC_ALGEBRAIC_SEARCH_SPACE",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It parameterizes possible brackets but records no perturbation compositions and selects no bracket or physical current.",
            ("Stage 1 is search-space only", "Stage 2 full classification remains open"),
        ),
        candidate(
            "b14_oriented_face_bracket",
            "STATIC_ONLY",
            ("code/b14_jacobi/oriented_face_bracket_selector.certificate.json",),
            "Incidence-defined equal-weight bracket on twenty oriented faces.",
            "STATIC_INCIDENCE_CONSTRUCTION",
            properties(
                source_native=True,
                port_indexed=True,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It is not extracted from histories and its Jacobi tensor is exactly nonzero, so it is a negative control rather than a current.",
            (f"Jacobi nonzero coordinates = {facts['b14_jacobi_nonzero_count']}", "Jacobi squared norm = 240"),
        ),
        candidate(
            "static_response_grammar_and_poles",
            "STATIC_ONLY",
            (
                "code/a5_closure/manifests/response_grammar_completeness_reference.json",
                "code/a5_closure/receipts/response_grammar_completeness_reference.receipt.json",
                "code/a5_closure/manifests/charged_response_pole_residue_artifact.json",
            ),
            "Bounded involution grammar, projectors, and pole/residue response tables.",
            "STATIC_OPERATOR_TABLES",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Reversible static operators are not twelve source perturbation families and no two ordered histories are serialized.",
            ("grammar forcing is independence-limited", "operator grammar is bounded and static"),
        ),
        candidate(
            "spin_deck_transport",
            "STATIC_ONLY",
            ("code/a5_closure/manifests/spin_statistics_semantic_artifact.json",),
            "Finite deck-lift/group transport tables on the carrier and refinement maps.",
            "STATIC_GROUP_TABLE",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "Group elements are reversible algebraically, but there are no per-port perturbation histories and the physical source gate is false.",
            ("finite source-model gate passes", "continuum and laboratory rows remain false"),
        ),
        candidate(
            "declared_port_current_fixture",
            "DOWNSTREAM_CONTAMINATED",
            (
                "code/a5_closure/manifests/port_current_response_reference.json",
                "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
            ),
            "Conditional charged-double-triplet matrix-current fixture.",
            "DOWNSTREAM_MATRIX_FIXTURE",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=False,
                same_twelve_port_carrier=True,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "Its generators, bracket, block decomposition, and implementers are the forbidden comparison fixture, not reconstruction inputs.",
            ("construction_model = charged_double_triplet", "response status is declared fixture"),
        ),
        candidate(
            "downstream_kinetic_form_selector",
            "DOWNSTREAM_CONTAMINATED",
            ("code/angular_sprint/runtime/kinetic_form_selection_receipt.json",),
            "Kinetic-form analysis rebuilt from the conditional matrix-current receipt.",
            "DOWNSTREAM_FIXTURE_CONSUMER",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=False,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "It imports the charged-double-triplet current and explicitly leaves current-lift source selection open.",
            ("current_lift_source_selected = false", "current_fixture is charged-double-triplet"),
        ),
        candidate(
            "generic_transport_from_gluing",
            "UNRELATED",
            ("code/consensus/transport_from_gluing.py",),
            "Generic caller-supplied edge-matrix path products and holonomy obstruction.",
            "GENERIC_ABSTRACTION_NO_SOURCE_INSTANCE",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The caller supplies the matrices; there is no source-native twelve-port instance, history custody, or refinement ancestry.",
            ("generic path-product API", "no committed current packet"),
        ),
        candidate(
            "lean_source_history_abstractions",
            "PARTIAL",
            (
                "Lean/InformationProjection/SourceHistoryPacket.lean",
                "Lean/Variational/SourceToHamiltonianComposed.lean",
                "Lean/QFT/SourceOperatorGeneration.lean",
                "Lean/QFT/SourceCorrelationCapstone.lean",
                "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean",
            ),
            "Finite source-history packets, observed paths, and supplied Markov/instrument interfaces.",
            "FORMAL_INTERFACE_WITHOUT_RUNTIME_INSTANCE",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "These theorems type supplied or post-processed histories, but use Fin8/Fin13/abstract carriers and instantiate neither a W12 perturbation family nor mixed-order response data.",
            ("formal source-history structures exist", "source selection/runtime custody is explicitly not claimed"),
        ),
        candidate(
            "lean_conditional_a2_words",
            "COMMUTATIVE_ONLY",
            (
                "Lean/ObserverPatchHolography/A2EndpointCommutator.lean",
                "Lean/Screen/RepairWordCarrierReadout.lean",
            ),
            "Conditional six-axis/twelve-orientation step words and cumulative port-load readback.",
            "FORMAL_CONDITIONAL_WORD_ALGEBRA",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The files explicitly lack emitted step maps/diamonds and quotient the declared words by commuting relations; no runtime source packet instantiates them.",
            ("formal steps have inverses", "fifteen commutations force an abelian quotient", "source step maps are absent"),
        ),
        candidate(
            "lean_reversible_flow_models",
            "MISSING_ORDER_INFORMATION",
            (
                "Lean/QFT/SourceHistoryGNSDynamics.lean",
                "Lean/Screen/CarrierEvolutionFlow.lean",
            ),
            "Promoted unitary/GNS and one-parameter carrier-flow models.",
            "HISTORY_RECONSTRUCTED_FROM_DECLARED_GENERATOR",
            properties(
                source_native=False,
                port_indexed=False,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "The flows are mathematically reversible but reconstructed from supplied generators on M8/Fin5 or one aggregate Fin12 mode, with no twelve P_p families or mixed-order histories.",
            ("one-parameter/unitary group laws are formalized", "no per-port source perturbation packet"),
        ),
        candidate(
            "lean_declared_fin12_transport_refinement",
            "PARTIAL",
            (
                "Lean/Geometry/CarrierDynamicsCompatibility.lean",
                "Lean/Geometry/WorldlineHopTransport.lean",
                "Lean/Screen/PrimitivePortTranslationBridge.lean",
            ),
            "Declared Fin12 bijections, seam steps, translations, histories, and geometric transport/refinement interfaces.",
            "DECLARED_TRANSPORT_WITH_EXTERNAL_HISTORIES",
            properties(
                source_native=False,
                port_indexed=True,
                reversible=True,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Only a few global symmetries/declared steps are present; histories are arguments, and neither a run-pinned source producer nor all 132 mixed orders is supplied.",
            ("Fin12 reversible maps exist", "source/readout attachment remains open"),
        ),
        candidate(
            "lean_irreversible_provenance_histories",
            "IRREVERSIBLE_ONLY",
            (
                "Lean/QFT/SourceDerivedEventPrecedence.lean",
                "Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean",
                "Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean",
            ),
            "Read-from/Boolean commit histories with two serializations and a small provenance/refinement certificate.",
            "RAW_SMALL_CARRIER_IRREVERSIBLE_HISTORY",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=True,
                target_free=True,
                same_twelve_port_carrier=False,
                refinement_provenance_present=True,
                twelve_reversible_perturbation_families=False,
            ),
            "The Boolean commits/read-from edges are irreversible, order-independent where both serializations occur, and live on Fin2/Fin4-to-Fin3 rather than the twelve-port carrier.",
            ("two Boolean serializations are formalized", "refinement map is non-injective and supplied"),
        ),
        candidate(
            "lean_transport_and_current_boundaries",
            "STATIC_ONLY",
            (
                "Lean/Geometry/SourceFeedbackTransport.lean",
                "Lean/Geometry/SourceSeamPathTomography.lean",
                "Lean/Screen/A2HolonomyBridge.lean",
                "Lean/Screen/A5ResponseWordAlgebra.lean",
                "Lean/Screen/OrientedFaceBracketSelector.lean",
            ),
            "Kernel-checked finite statements mirroring destructive transports and static current obstructions.",
            "FORMALIZATION_OF_EXISTING_STATIC_OR_IRREVERSIBLE_OBJECTS",
            properties(
                source_native=True,
                port_indexed=False,
                reversible=False,
                both_composition_orders_recorded=False,
                raw_histories_serialized=False,
                target_free=True,
                same_twelve_port_carrier=True,
                refinement_provenance_present=False,
                twelve_reversible_perturbation_families=False,
            ),
            "Formal proofs do not create the absent runtime packet; they preserve the same destructive/static boundaries.",
            ("response-word algebra is commutative", "oriented-face bracket fails Jacobi"),
        ),
    ]
    require(
        len({row["candidate_id"] for row in rows}) == len(rows),
        "CANDIDATE_ID",
        "duplicate candidate id",
    )
    return rows


def build_inventory() -> dict[str, Any]:
    facts = audit_source_semantics()
    candidates = build_candidates(facts)
    qualifying = [
        row
        for row in candidates
        if all(row["properties"][field] for field in QUALIFICATION_FIELDS)
    ]
    require(not qualifying, "STOP_CONDITION", "a qualifying candidate requires Stage-2 review")
    counts = Counter(row["classification"] for row in candidates)
    inventory: dict[str, Any] = {
        "schema": SCHEMA,
        "program": "SOURCE-CURRENT-TOMOGRAPHY-0",
        "stage": 2,
        "audited_upstream_main_sha": UPSTREAM_MAIN_SHA,
        "verdict": VERDICT,
        "positive_stages_authorized": False,
        "physical_current_source_bridge_attained": False,
        "fixture_used_as_reconstruction_oracle": False,
        "target_labels_used": False,
        "audit_scope": {
            "rule": "one candidate packet must satisfy every qualification field; facts from distinct packets are not composited",
            "directories": [
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
            ],
            "repository_wide_search_terms": [
                "port-indexed perturbation",
                "reversible perturbation",
                "ordered response/history/composition",
                "mixed-order",
                "raw runtime history serialized",
                "sequential instrument",
                "commutator",
                "overlap word",
                "recharting",
                "refinement provenance",
            ],
        },
        "qualification_fields": list(QUALIFICATION_FIELDS),
        "source_pins": [source_pin(path) for path in AUDITED_SOURCE_PATHS],
        "implementation_pins": [
            source_pin("code/a5_closure/source_current_order_sensitive_inventory.py"),
            source_pin("code/a5_closure/verify_source_current_order_sensitive_inventory.py"),
        ],
        "candidates": candidates,
        "summary": {
            "candidate_count": len(candidates),
            "qualifying_candidate_count": len(qualifying),
            "classification_counts": dict(sorted(counts.items())),
            "strongest_registered_source_packet": {
                "candidate_id": "registered_adjacency_response",
                "response_word_dimension": facts["response_algebra_dimension"],
                "commutator_nonzero_count": facts[
                    "response_commutator_nonzero_count"
                ],
                "proper_rechartings_total": facts["proper_recharting_count"],
                "proper_rechartings_in_response_words": 1,
            },
            "richest_logged_near_candidate": {
                "candidate_id": "record_counting_repair_628",
                "fatal_failures": [
                    "irreversible strict-descent repair",
                    "presentation order erased by commutative addition",
                    "no twelve reversible perturbation families",
                ],
            },
            "cross_packet_noncomposition": (
                "The reversible four-dimensional response packet and the logged "
                "record-repair packet have different primitives and semantics. "
                "Combining their complementary fields would be a new producer, not "
                "an inventory finding."
            ),
        },
        "minimal_missing_fields": [
            {
                "field": "perturbation_families",
                "requirement": "exactly twelve reversible source-native families P_0,...,P_11, one per primitive port",
            },
            {
                "field": "ordered_mixed_histories",
                "requirement": "both p-then-q and q-then-p for all 66 unordered distinct-port pairs (132 ordered compositions)",
            },
            {
                "field": "raw_history_custody",
                "requirement": "canonical serialized raw histories with provenance and tamper-evident hashes, not reconstruction from one static generator",
            },
            {
                "field": "common_source_binding",
                "requirement": "all histories share the same twelve-port carrier, orientation, normalization, source family, and refinement stage/tower",
            },
            {
                "field": "reversibility_witness",
                "requirement": "an inverse/two-sided law for every perturbation family, excluding reset, averaging, measurement, and strict-descent repair",
            },
            {
                "field": "source_firewall",
                "requirement": "no current fixture, gauge/particle label, measurement target, or downstream selector enters the producer",
            },
        ],
        "claim_boundary": {
            "bounded_repository_inventory_only": True,
            "does_not_assert_a_no_go_for_future_source_data": True,
            "does_not_weaken_the_abstract_forced_lie_type_theorem": True,
            "does_not_reject_the_conditional_port_current_fixture": True,
            "does_not_make_a_physical_current_claim": True,
            "stages_3_through_6_not_run": True,
            "next_scientific_step": "add and independently review a new source producer containing the minimal missing fields; do not reinterpret existing static or irreversible packets",
        },
        "verifier_command": (
            "python3 code/a5_closure/verify_source_current_order_sensitive_inventory.py "
            "--inventory code/a5_closure/manifests/source_current_order_sensitive_inventory.json"
        ),
    }
    inventory["inventory_sha256"] = canonical_sha256(inventory)
    return inventory


def verify_inventory(inventory: Mapping[str, Any]) -> None:
    require(inventory.get("schema") == SCHEMA, "INVENTORY_SCHEMA", "wrong schema")
    body = {key: value for key, value in inventory.items() if key != "inventory_sha256"}
    require(
        inventory.get("inventory_sha256") == canonical_sha256(body),
        "INVENTORY_HASH",
        "inventory self-hash failed",
    )
    require(
        dict(inventory) == build_inventory(),
        "INVENTORY_REPLAY",
        "committed inventory does not replay exactly",
    )


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("all", "show", "verify"))
    parser.add_argument("--inventory", type=Path, default=INVENTORY_PATH)
    args = parser.parse_args()

    if args.command == "all":
        inventory = build_inventory()
        write_json(args.inventory, inventory)
        print(VERDICT)
        return 0
    if args.command == "show":
        print(json.dumps(build_inventory(), indent=2, sort_keys=True))
        return 0
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    require(isinstance(inventory, dict), "JSON_OBJECT", "inventory is not an object")
    verify_inventory(inventory)
    print("SOURCE_CURRENT_ORDER_SENSITIVE_INVENTORY_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

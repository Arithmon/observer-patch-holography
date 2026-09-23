#!/usr/bin/env python3
"""SOURCE-CURRENT-TOMOGRAPHY-0 Stage 1: exact order-sensitive gap contract.

Consumes the existing independent issue-566 source-current capability receipt.
It does not reconstruct a current and does not consume PORT-CURRENT-INNER as
an oracle.  Its purpose is to turn the bounded upstream obstruction into an
explicit minimal contract for the next source packet.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
CAPABILITY = HERE / "receipts" / "source_current_capability.receipt.json"

EXPECTED_SCHEMA = "oph.source_current_capability_receipt.v1"
EXPECTED_VERDICT = (
    "BOUNDED_REGISTERED_RESPONSE_ALGEBRA_INSUFFICIENT__"
    "ORDER_SENSITIVE_SOURCE_BRIDGE_OPEN"
)


class Stage1Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage1Error(message)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    require(isinstance(value, dict), f"{path} must contain a JSON object")
    return value


def payload() -> dict[str, Any]:
    receipt = load_json(CAPABILITY)
    require(receipt.get("schema") == EXPECTED_SCHEMA, "unexpected capability schema")
    require(receipt.get("verdict") == EXPECTED_VERDICT, "capability verdict drift")

    audit = receipt.get("audit", {})
    algebra = audit.get("response_word_algebra", {})
    rechart = audit.get("recharting_audit", {})
    unitary = audit.get("registered_unitary_channel_audit", {})
    refinement = audit.get("refinement_audit", {})
    classifications = audit.get("acceptance_classifications", {})
    next_source = receipt.get("next_source_object", {})

    current_rank = int(algebra.get("exact_dimension", -1))
    current_commutator_rank = int(algebra.get("commutator_nonzero_count", -1))
    current_rechartings = int(rechart.get("proper_rechartings_in_response_word_algebra", -1))
    target_rank = 12
    target_derived_rank = 11
    target_rechartings = int(rechart.get("proper_recharting_count", -1))

    require(current_rank == 4, "registered response-word algebra dimension drift")
    require(algebra.get("commutative") is True, "registered algebra no longer commutative")
    require(current_commutator_rank == 0, "registered commutator count drift")
    require(target_rechartings == 60, "proper recharting count drift")
    require(current_rechartings == 1, "response-word recharting membership drift")
    require(unitary.get("order_sensitive_port_indexed_tangent_available") is False,
            "order-sensitive tangent unexpectedly available")
    require(refinement.get("generator_bracket_implementer_intertwining_available") is False,
            "generator/bracket refinement unexpectedly available")

    deficits = {
        "algebraic_tangent_dimension_deficit": target_rank - current_rank,
        "native_tangent_rank_reconstructed": False,
        "derived_commutator_dimensions_missing": target_derived_rank,
        "nonidentity_rechartings_missing_from_response_words": target_rechartings - current_rechartings,
        "same_word_projective_implementers_missing": (
            classifications.get("same_word_projective_implementers") == "MISSING"
        ),
    }

    source_packet_contract = {
        "carrier_binding": "same pinned twelve-port carrier and orientation receipt",
        "perturbations": {
            "count": 12,
            "indexing": "one reversible perturbation family per primitive port",
            "target_labels_forbidden": True,
            "downstream_particle_or_gauge_labels_forbidden": True,
        },
        "ordered_histories": {
            "required": True,
            "both_composition_orders_per_unordered_port_pair": True,
            "same_source_packet": True,
            "unordered_port_pairs": 12 * 11 // 2,
            "ordered_compositions": 12 * 11,
            "reason": (
                "the registered single-generator functional calculus is commutative; "
                "order-sensitive mixed perturbations are the first source object capable "
                "of exposing a nonzero commutator"
            ),
        },
        "first_order_gate": {
            "tangent_generator_count": 12,
            "real_rank_exactly": 12,
            "skew_adjoint": True,
            "positive_pairing": "-Re tr(XY) positive definite on the source span",
        },
        "mixed_order_gate": {
            "commutators_reconstructed_from_order_difference": True,
            "all_commutators_close_in_source_span": True,
            "commutator_span_real_rank": 11,
            "center_dimension": 1,
            "center_source_characterization": (
                "constant linear combination of the twelve port generators"
            ),
        },
        "holonomy_gate": {
            "closed_overlap_words_source_derived": True,
            "proper_rechartings_covered": 60,
            "same_history_projective_implementers": True,
            "internal_modulo_pointwise_centralizer": True,
        },
        "refinement_gate": {
            "generator_intertwining": True,
            "bracket_intertwining": True,
            "implementer_intertwining": True,
            "same_registered_refinement_tower": True,
        },
    }

    falsifiers = [
        "first-order tangent rank < 12",
        "all mixed-order differences vanish",
        "commutator span rank != 11",
        "commutator leaves the twelve-dimensional source span",
        "center dimension != 1",
        "nonidentity proper recharting words are absent",
        "implementers require a matrix-current fixture not reconstructed from the same histories",
        "results change under allowed source relabeling/refinement without the declared intertwiner",
        "target/gauge/particle labels enter the source producer",
    ]

    exits = {
        "positive": (
            "source histories reconstruct a unique current up to the declared basis/gauge "
            "equivalence; comparison with PORT-CURRENT-INNER occurs only after reconstruction"
        ),
        "nonidentifiable": (
            "two inequivalent current realizations survive the identical source packet; "
            "emit explicit witnesses and the invariant separating their equivalence classes"
        ),
        "insufficient": (
            "the packet fails one or more first-order, mixed-order, holonomy, or refinement gates"
        ),
    }

    return {
        "schema": "arithmon.source_current_tomography_stage1_contract.v2",
        "upstream_capability_verdict": receipt["verdict"],
        "registered_response_word_algebra": {
            "basis": algebra.get("basis"),
            "dimension": current_rank,
            "commutative": algebra.get("commutative"),
            "commutator_nonzero_count": current_commutator_rank,
            "minimal_polynomial": algebra.get("minimal_polynomial"),
            "positive_skew_pairing": algebra.get("skew_pairing_positive_definite"),
        },
        "algebraic_capacity_gap": deficits,
        "next_source_object": next_source,
        "source_packet_contract": source_packet_contract,
        "falsifiers": falsifiers,
        "allowed_exits": exits,
        "claim_boundary": {
            "bounded_obstruction_only": True,
            "does_not_assert_no_go_for_order_sensitive_histories": True,
            "does_not_consume_conditional_current_fixture": True,
            "does_not_promote_physical_current": True,
            "abstract_forced_lie_type_theorem_preserved": receipt.get(
                "abstract_forced_lie_type_theorem_preserved"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("show", "verify"))
    args = parser.parse_args()

    result = payload()
    if args.command == "verify":
        gap = result["algebraic_capacity_gap"]
        require(gap["algebraic_tangent_dimension_deficit"] == 8, "algebraic capacity gap drift")
        require(gap["derived_commutator_dimensions_missing"] == 11, "derived-rank gap drift")
        require(gap["nonidentity_rechartings_missing_from_response_words"] == 59,
                "recharting gap drift")
        contract = result["source_packet_contract"]
        require(
            contract["ordered_histories"]["unordered_port_pairs"] == 66
            and contract["ordered_histories"]["ordered_compositions"] == 132,
            "ordered-pair contract drift",
        )
        require(result["claim_boundary"]["does_not_consume_conditional_current_fixture"],
                "fixture firewall failed")
        print("SOURCE-CURRENT-TOMOGRAPHY-0 stage-1 contract: verified")
        return 0

    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fail-closed source-information audit for SOURCE-CURRENT-TOMOGRAPHY-0.

This stage does not reconstruct a physical current. It asks whether the
currently committed source-bound artifacts contain the information required
to identify a Lie current from ordered response histories.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
SEMANTIC = HERE / "manifests" / "charged_response_semantic_artifact.json"
B14 = HERE.parent / "b14_jacobi" / "oriented_face_bracket_selector.certificate.json"


class AuditError(RuntimeError):
    pass


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise AuditError(f"{path} must contain a JSON object")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def audit_payload() -> dict[str, Any]:
    semantic = load_json(SEMANTIC)
    b14 = load_json(B14)

    require(
        semantic.get("schema") == "oph.charged_response_semantic_artifact.v3",
        "unexpected charged-response semantic schema",
    )
    require(
        b14.get("schema") == "oph.b14.oriented_face_bracket_selector.v1",
        "unexpected B14 selector schema",
    )

    lift = semantic.get("derived", {}).get("current_lift_status", {})
    response = semantic.get("source_response", {})
    bands = semantic.get("response_basis", {}).get("sector_dimensions", {})
    jacobi = b14.get("jacobi_failure", {})
    face = b14.get("source_face_bracket", {})
    discriminator = b14.get("orthogonal_compact_locus_discriminator", {})

    require(
        bands == {
            "unit_band": 1,
            "quintet_band": 5,
            "frame_band": 3,
            "kernel_band": 3,
        },
        "unexpected 1+5+3+3' sector dimensions",
    )
    require(
        response.get("target_labels_used") is None
        or response.get("impulse_readback_protocol", {}).get("target_labels_used") is False,
        "target labels entered the source response protocol",
    )
    require(
        response.get("impulse_readback_protocol", {}).get("downstream_labels_used") is False,
        "downstream labels entered the source response protocol",
    )

    source_bound_constraints = {
        "carrier_bound": bool(semantic.get("carrier_binding")),
        "target_blind_impulse_readback": bool(
            response.get("impulse_readback_response_executed")
        ),
        "response_sector_signs": response.get("sector_eigenvalues")
        == {
            "unit_band": -1,
            "quintet_band": -1,
            "frame_band": 1,
            "kernel_band": 1,
        },
        "orientation_convention": bool(semantic.get("orientation_convention")),
        "refinement_maps": bool(semantic.get("physical_refinement_maps")),
    }

    missing_tomography_data = {
        "ordered_two_sided_response_histories": False,
        "twelve_source_reconstructed_infinitesimal_generators": False,
        "source_reconstructed_positive_pairing": False,
        "source_reconstructed_commutator": bool(
            lift.get("commutator_reconstructed_from_ordered_response")
        ),
        "same_history_overlap_words": False,
        "same_current_projective_implementers": bool(
            lift.get("overlap_holonomy_internality_certified")
        ),
    }

    negative_control = {
        "name": "oriented_face_incidence_bracket",
        "source_face_definition_present": bool(face.get("definition")),
        "oriented_face_count": face.get("oriented_face_count"),
        "proper_action_order": face.get("proper_action_order"),
        "jacobi_nonzero_coordinate_count": jacobi.get("nonzero_count"),
        "jacobi_squared_norm": jacobi.get("squared_norm"),
        "fails_jacobi": bool(jacobi.get("nonzero_count", 0)),
        "nearest_compact_family_is_source_selection": (
            discriminator.get("conclusion_status")
            != "TARGET_CLEAN_CONDITIONAL_DISCRIMINATOR__NOT_SOURCE_SELECTION"
        ),
        "metric_is_source_derived": discriminator.get("metric_is_source_derived"),
        "repair_rule_is_source_derived": discriminator.get(
            "minimum_hs_or_jacobi_repair_is_source_derived"
        ),
    }

    enough_for_tomography = (
        all(source_bound_constraints.values())
        and all(missing_tomography_data.values())
    )

    verdict = (
        "SOURCE_CURRENT_IDENTIFIED"
        if enough_for_tomography
        else "INSUFFICIENT_SOURCE_DATA"
    )

    return {
        "schema": "arithmon.source_current_tomography_stage0.v1",
        "status": verdict,
        "source_bound_constraints": source_bound_constraints,
        "tomography_requirements": missing_tomography_data,
        "negative_control": negative_control,
        "conclusions": [
            "The committed source response fixes the carrier, target-blind antipodal response, four relative sector signs, orientation convention, and refinement persistence.",
            "The committed source response does not reconstruct the twelve infinitesimal current generators, their Lie bracket, or same-current overlap holonomy from ordered response histories.",
            "The oriented-face incidence bracket is an exact warning that incidence and A5-compatible source structure alone do not force a Lie bracket: its Jacobi tensor is nonzero.",
            "The B14 nearest-compact-family comparison is explicitly conditional on a non-source-derived metric and repair rule, so it cannot repair the missing source selector.",
        ],
        "claim_boundary": {
            "does_not_reject_conditional_port_current_algebra": True,
            "does_not_select_charged_double_triplet_fixture": True,
            "does_not_make_physical_current_claim": True,
            "next_positive_evidence_required": [
                "ordered two-sided port-response histories",
                "twelve reconstructed skew-adjoint tangent generators",
                "exact commutator reconstruction in their span",
                "closed overlap words and same-current implementers",
            ],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("audit", "verify"))
    args = parser.parse_args()

    payload = audit_payload()
    if args.command == "verify":
        require(payload["status"] == "INSUFFICIENT_SOURCE_DATA", "stage-0 verdict drift")
        require(payload["negative_control"]["fails_jacobi"], "B14 control no longer fails Jacobi")
        require(
            payload["negative_control"]["metric_is_source_derived"] is False,
            "B14 metric unexpectedly became source-derived",
        )
        require(
            payload["negative_control"]["repair_rule_is_source_derived"] is False,
            "B14 repair rule unexpectedly became source-derived",
        )
        print("SOURCE-CURRENT-TOMOGRAPHY-0 stage-0: verified fail-closed")
        return 0

    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

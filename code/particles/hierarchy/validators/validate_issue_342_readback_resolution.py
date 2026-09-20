#!/usr/bin/env python3
"""Validate the finite readback-resolution certificate in its current scope."""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from citation_resolution import scan_certificate  # noqa: E402


def decimal(value: str | None, default: str = "0") -> Decimal:
    return Decimal(value if value is not None else default)


def main(path: str = "certificates/R_readback_resolution_certificate.json") -> int:
    cert = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    normal = cert.get("normal_form", {})
    observer = cert.get("observer_sector", {})
    register = cert.get("capacity_register", {})
    selected = [row for row in register.get("atoms", []) if row.get("selected") is True]
    extractor = cert.get("extractor", {})
    readback = cert.get("readback_resolution", {})
    strict = cert.get("strict_source_capacity_check", {})
    refinement = cert.get("refinement_certificate", {})
    boundary = cert.get("claim_boundary", {})
    acceptance = cert.get("acceptance_criteria_status", {})
    ledger = cert.get("input_ledger", {})
    source_status = cert.get("source_status", {})
    unresolved = scan_certificate(cert, __file__)

    used = set(ledger.get("used_inputs", []))
    forbidden = set(ledger.get("forbidden_inputs", []))
    selected_capacity = decimal(selected[0].get("capacity")) if len(selected) == 1 else Decimal(0)
    tolerance = decimal(strict.get("relative_tolerance"), "1e-30")

    checks = {
        "issue_is_342": cert.get("issue") == 342,
        "accepted": cert.get("accepted") is True,
        "status_is_closed_finite_readback": (
            cert.get("status") == "closed_finite_readback_resolution_certificate"
        ),
        "normal_form_is_unique_at_fixed_cutoff": (
            normal.get("finite_state") is True
            and normal.get("unique") is True
            and normal.get("schedule_independent") is True
        ),
        "normal_form_cites_confluence_theorem": (
            "OBSERVERS_APPENDICES" in str(normal.get("source_artifact", ""))
        ),
        "observer_sector_is_stable_and_central": (
            observer.get("nonzero") is True
            and observer.get("stable_self_reading") is True
            and observer.get("central_projector") is True
            and observer.get("capacity_visible_on_quotient") is True
        ),
        "one_positive_sharp_capacity_atom": (
            len(selected) == 1
            and selected_capacity > 0
            and decimal(selected[0].get("probability")) == 1
            and decimal(register.get("selected_variance")) == 0
        ),
        "positive_root_extractor": (
            extractor.get("formula") == "rho_read = sqrt(pi / Cap_read)"
            and extractor.get("positive_root") is True
        ),
        "exact_capacity_readback_matches": (
            decimal(readback.get("cap_read")) == decimal(strict.get("N_CRC"))
            and decimal(readback.get("rho_read")) == decimal(strict.get("rho_star"))
            and abs(decimal(strict.get("capacity_relative_residual"))) <= tolerance
            and abs(decimal(strict.get("rho_relative_residual"))) <= tolerance
        ),
        "refinement_statement_is_explicit": (
            refinement.get("uniform_limit_declared") is True
            and refinement.get("positive_root_closure") is True
            and Decimal(0) <= decimal(refinement.get("contraction_kappa")) < Decimal(1)
            and refinement.get("kappa_matches_ew_lambda") is True
        ),
        "rounded_display_is_diagnostic_only": (
            readback.get("display_only") is True
            and cert.get("obstruction_record", {}).get("rounded_N_CRC_status")
            == "diagnostic_only_not_exact_bridge_witness"
        ),
        "input_ledger_is_source_separated": (
            not (used & forbidden)
            and "rounded_3p31e122_capacity_display_as_bridge_witness" in forbidden
        ),
        "scope_is_conjunction_limited": (
            "conjunction of the seven branches" in boundary.get("scope", "")
            and boundary.get("not_closed_here", []) == []
            and "diagnostic-only" in boundary.get("scope", "")
        ),
        "acceptance_checks_pass": bool(acceptance) and all(acceptance.values()),
        "internal_checks_pass": bool(cert.get("checks")) and all(cert["checks"].values()),
        "does_not_promote_full_hierarchy": (
            source_status.get("closes_gate") == "finite_readback_resolution"
            and source_status.get("does_not_promote_full_hierarchy_resonance") is True
        ),
        "every_citation_resolves_in_repository_or_is_marked_external": (
            unresolved == []
        ),
    }
    payload = {
        "checks": checks,
        "unresolved_citations": unresolved,
        "pass": all(checks.values()),
    }
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(
        main(
            sys.argv[1]
            if len(sys.argv) > 1
            else "certificates/R_readback_resolution_certificate.json"
        )
    )

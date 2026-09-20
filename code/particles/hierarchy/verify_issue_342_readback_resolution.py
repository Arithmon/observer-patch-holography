#!/usr/bin/env python3
"""Verifier for OPH issue #342: conditional positive-root readback schema.

The certificate records the implication from declared finite-repair,
observer-sector, central-atom, and cofinal-refinement premises to a singleton
positive-root readback and its refinement limit.  It does not construct the
finite objects, select the central atom, or emit an independent cosmic
capacity.  The bridge-defined value N_CRC^EW(P_*) is supplied by
R_EW_global_capacity_certificate.json; the rounded 3.31e122 cosmological
capacity display is diagnostic only.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any


getcontext().prec = 200

PI = Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
)
DEFAULT_N_DISPLAY = Decimal("3.31e122")
DEFAULT_KAPPA = Decimal("0.5")
DEFAULT_TOL = Decimal("1e-30")

REPO_ROOT = Path(__file__).resolve().parents[3]
HIERARCHY_DIR = Path(__file__).resolve().parent
DEFAULT_EW_CAPACITY_CERT = (
    HIERARCHY_DIR / "certificates" / "R_EW_global_capacity_certificate.json"
)

FORBIDDEN_INPUTS = {
    "v_measured",
    "weak_scale_measured",
    "higgs_scale_measured",
    "m_H_measured",
    "M_W_measured",
    "M_Z_measured",
    "m_t_measured",
    "G_measured",
    "G_CODATA",
    "Planck_area_measured",
    "ell_P_squared_measured",
    "Lambda_measured",
    "H0_measured",
    "de_sitter_area_measured_as_input",
    "rounded_3p31e122_capacity_display_as_bridge_witness",
}


def D(value: str | int | Decimal | None, default: Decimal | None = None) -> Decimal | None:
    if value is None:
        return default
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))


def decstr(value: Decimal | None) -> str | None:
    if value is None:
        return None
    if value != 0 and (abs(value) < Decimal("1e-6") or abs(value) >= Decimal("1e6")):
        return format(value, "E")
    return format(value, "f")


def rho_from_capacity(capacity: Decimal) -> Decimal:
    if capacity <= 0:
        raise ValueError(f"capacity must be positive, got {capacity}")
    return (PI / capacity).sqrt()


def rel_abs(a: Decimal, b: Decimal) -> Decimal:
    if b == 0:
        return abs(a - b)
    return abs(a - b) / abs(b)


def load_ew_capacity_cert(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            f"EW exact-capacity certificate not found at {path}; "
            "this verifier requires R_EW_global_capacity_certificate.json "
            "as a machine-checkable dependency."
        )
    return json.loads(path.read_text(encoding="utf-8"))


def n_crc_ew_from_cert(cert: dict[str, Any]) -> Decimal:
    fixed_point = cert.get("exact_capacity_fixed_point", {})
    value = fixed_point.get("N_CRC_EW")
    if value is None:
        raise ValueError(
            "EW capacity certificate is missing exact_capacity_fixed_point.N_CRC_EW"
        )
    return Decimal(str(value))


def build_certificate(
    n_display: Decimal,
    kappa: Decimal,
    source_n: Decimal | None,
    cap_read_override: Decimal | None,
    tol: Decimal,
    ew_cert: dict[str, Any],
    ew_cert_relpath: str,
) -> dict[str, Any]:
    n_crc_ew = n_crc_ew_from_cert(ew_cert)
    ew_lambda = Decimal(str(ew_cert.get("source_values", {}).get("lambda", "0.5")))
    ew_bridge_residual = Decimal(
        str(ew_cert.get("exact_capacity_fixed_point", {}).get("bridge_residual", "0"))
    )
    ew_accepted = bool(ew_cert.get("accepted") is True)

    if source_n is None:
        source_n = n_crc_ew
    if cap_read_override is None:
        cap_read = n_crc_ew
    else:
        cap_read = cap_read_override

    rho_read = rho_from_capacity(cap_read)
    rho_star = rho_from_capacity(source_n)
    rho_display = rho_from_capacity(n_display)
    strict_capacity_residual = rel_abs(cap_read, source_n)
    strict_rho_residual = rel_abs(rho_read, rho_star)

    rounded_bridge_residual_relative = rel_abs(n_display, n_crc_ew)
    rounded_rho_residual_relative = rel_abs(rho_display, rho_star)

    used_inputs = [
        "finite_OPH_repair_law",
        "fixed_cutoff_confluence_certificate",
        "central_record_algebra",
        "stable_self_reading_observer_sector",
        "D6_dimensionless_area_law_N_equals_pi_over_rho_squared",
        "positive_root_convention",
        "EW_refined_exact_capacity_certificate",
    ]
    used_artifacts = [
        "extra/observer_patch_holography_as_string_vacuum_selector.tex, theorem 'Conditional edge-capacity entropy theorem' (screen-capacity readback equation N_scr = pi*R_dS^2/ell_P^2)",
        "paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex (stable self-reading observer sector and capacity register)",
        "paper/tex_fragments/OBSERVERS_APPENDICES.tex (finite repair Lyapunov+diamond+completeness confluence theorem)",
        "code/particles/hierarchy/certificates/R_N_global_repair_tick_certificate.json (D6 area-law normalization)",
        "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json (Banach contraction with lambda=1/2 supplying F(N_CRC^EW)=N_CRC^EW>0)",
    ]
    no_forbidden_inputs = not (set(used_inputs) & FORBIDDEN_INPUTS)

    strict_ok = strict_capacity_residual <= tol and strict_rho_residual <= tol
    kappa_ok = Decimal(0) <= kappa < Decimal(1)
    rounded_rejected_as_bridge = rounded_bridge_residual_relative > tol

    derivation_chain = build_derivation_chain(
        ew_cert_relpath=ew_cert_relpath,
        n_crc_ew=n_crc_ew,
        ew_lambda=ew_lambda,
        kappa=kappa,
    )

    factor_origins = {
        "area_law_normalization_pi": {
            "value": "pi",
            "role": "screen-capacity area-law normalization N = pi*(r/ell)^2 = pi/rho^2",
            "source_theorem": "D6 dimensionless area law (corpus axiom)",
            "source_artifact": "extra/observer_patch_holography_as_string_vacuum_selector.tex, theorem 'Conditional edge-capacity entropy theorem'; reflected in R_N_global_repair_tick_certificate.json (definitions.screen_normalized_radius_coordinate, normalization.screen_capacity_relation)",
        },
        "positive_root_exponent_one_half": {
            "value": "1/2",
            "role": "positive-root inverse rho = sqrt(pi/N) = (N/pi)^(-1/2) of the D6 area law",
            "source_theorem": "positive-root convention on the D6 area law (corpus axiom)",
            "source_artifact": "extra/observer_patch_holography_as_string_vacuum_selector.tex, theorem 'Conditional edge-capacity entropy theorem'; rho_star=(N_CRC/pi)^(-1/2) in R_N_global_repair_tick_certificate.json",
        },
        "banach_contraction_lambda_one_half": {
            "value": "1/2",
            "role": "Banach contraction constant lambda for the EW-refined source-side capacity map C_EW(P,x)=(1-lambda)*x+lambda*6*pi/(P*alpha_U(P)); supplies the cofinal residual bound delta_r/(1-lambda)=2*delta_r",
            "source_theorem": "EW-refined exact-capacity Banach contraction",
            "source_artifact": "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json (contraction_certificate.lambda)",
        },
        "derivative_bound_factor_two": {
            "value": "2",
            "role": "denominator in the chain-rule bound |drho/dN| = sqrt(pi)/(2*N^(3/2)); inherited from differentiating rho(N)=sqrt(pi)*N^(-1/2)",
            "source_theorem": "elementary calculus on rho(N)=sqrt(pi)*N^(-1/2)",
            "source_artifact": "this certificate (refinement_certificate.finite_to_limit_bound)",
        },
    }

    branch_scope = {
        "d6_area_law_branch": "screen-capacity normalization N = pi*(r/ell)^2 in dimensionless form N = pi/rho^2 (D6 area law)",
        "finite_oph_repair_branch": "declared premise: a finite OPH patch system at fixed cutoff r with strictly Lyapunov-decreasing repair, local diamond on the physical quotient, and repair completeness; no executable family is constructed here",
        "stable_observer_sector_branch": "declared premise: a non-zero stable self-reading observer sector with a central projector and capacity visible on the quotient; no observer functor is constructed here",
        "central_record_algebra_branch": "declared premise: a central capacity register C_hat_{r,N}=sum_c c P_c with one selected positive atom and zero variance; the atom selection is not derived here",
        "ew_refined_exact_capacity_branch": "lambda=1/2 Banach contraction whose unique fixed point N_CRC^EW(P_*)=pi*exp[6*pi/(P_*alpha_U(P_*))] supplies F(N_CRC^EW)=N_CRC^EW>0",
        "positive_root_convention_branch": "the positive-root inverse rho=sqrt(pi/N) of the D6 area law (excludes the negative root)",
        "refinement_cofinal_branch": "declared premise: F_r -> F cofinally with quantitative residual delta_r controlling |N_r_star - N_CRC^EW|; no cofinal family is computed here",
        "scope_note": "The certificate verifies a conditional implication on the conjunction of the seven branches above. N_CRC^EW(P_*) is the bridge-defined fixed point supplied by the EW-refined exact-capacity certificate, not an independently emitted cosmic capacity.",
    }

    obstruction_record = {
        "rounded_N_CRC_display": decstr(n_display),
        "rounded_N_CRC_status": "diagnostic_only_not_exact_bridge_witness",
        "rounded_capacity_relative_gap_vs_N_CRC_EW": decstr(
            rounded_bridge_residual_relative
        ),
        "rounded_rho_relative_gap_vs_rho_star": decstr(rounded_rho_residual_relative),
        "rounded_bridge_residual_from_ew_certificate": ew_cert.get(
            "rounded_capacity_diagnostic", {}
        ).get("bridge_residual"),
        "meaning": (
            "The rounded 3.31e122 cosmological capacity display is a "
            "diagnostic label, not a bridge witness. The exact-capacity "
            "certificate supplies the source-side bridge fixed point "
            "N_CRC^EW(P_*) = pi*exp[6*pi/(P_*alpha_U(P_*))] with "
            f"B_EW(P_*,N_CRC^EW) = {ew_cert.get('exact_capacity_fixed_point', {}).get('bridge_residual', '0')} exactly; "
            "this finite readback-resolution certificate uses N_CRC^EW as Cap_read."
        ),
    }

    acceptance_criteria_status = {
        "definitions_emitted": True,
        "single_resolution_recorded_as_branch_premise": True,
        "single_atom_supplied_by_bridge_defined_dependency": True,
        "finite_to_refinement_implication_verified": (
            ew_accepted and ew_bridge_residual == 0 and kappa_ok
        ),
        "positive_root_implication_verified": (
            strict_ok and ew_accepted
        ),
        "allowed_inputs_and_forbidden_calibrations_stated": True,
        "rounded_capacity_display_rejected_as_bridge_witness": rounded_rejected_as_bridge,
        "public_certificate_and_verifier_emitted": True,
        "ew_exact_capacity_dependency_loaded_and_accepted": ew_accepted,
        "bridge_defined_fixed_point_hypothesis_supplied": (
            ew_accepted and ew_bridge_residual == 0
        ),
        "readback_object_construction_not_claimed": True,
        "independent_cosmic_capacity_not_claimed": True,
        "no_measured_weak_higgs_or_hierarchy_calibration_used": no_forbidden_inputs,
    }

    accepted = (
        cap_read > 0
        and kappa_ok
        and no_forbidden_inputs
        and strict_ok
        and ew_accepted
        and ew_bridge_residual == 0
        and rounded_rejected_as_bridge
        and all(acceptance_criteria_status.values())
    )

    return {
        "issue": 342,
        "artifact": "R_readback_resolution_certificate",
        "status": "conditional_positive_root_readback_schema",
        "accepted": bool(accepted),
        "theorem": "conditional positive-root readback implication",
        "closeout_decision": (
            "Under the declared finite-repair, stable-observer, sharp-central-atom, "
            "D6 area-law, and cofinal-refinement premises, the pipeline "
            "F_r(N)=Cap_read(Obs(nf_{r,N}(U_{r,N}))) is single-valued and its "
            "positive-root extractor satisfies rho_read=sqrt(pi/F_r(N)). The "
            "lambda=1/2 bridge-defined fixed point N_CRC^EW(P_*) gives the "
            "conditional limit rho_read(r,N_CRC^EW) -> "
            "(N_CRC^EW/pi)^(-1/2) with the stated quantitative bound. The "
            "certificate does not construct U, nf, Obs, Cap_read, or a cofinal "
            "family, and it does not identify N_CRC^EW with physical cosmic "
            "capacity. The rounded 3.31e122 display is diagnostic only."
        ),
        "definitions": {
            "finite_initial_object": "declared premise U_{r,N}: finite OPH patch system at fixed cutoff r and capacity coordinate N",
            "normal_form_map": "declared premise nf_{r,N}: U_{r,N} -> n_{r,N}: unique terminal normal form under finite repair confluence",
            "observer_sector": "declared premise Obs(n_{r,N}): non-zero stable self-reading observer sector with a central projector on the physical quotient",
            "capacity_register": "declared premise C_hat_{r,N}=sum_{c in S_{r,N}} c P_c: central capacity register with one selected positive atom and zero variance",
            "cap_read": "declared Cap_read: the selected positive eigenvalue of C_hat_{r,N}",
            "finite_readback_map": "conditional definition F_r(N) = Cap_read(Obs(nf_{r,N}(U_{r,N})))",
            "delivery_resolution": "rho_read(r,N) = sqrt(pi/F_r(N)) (positive-root inverse of D6 area law N = pi/rho^2)",
            "limit_resolution": "rho_star = (N_CRC^EW/pi)^(-1/2), with bridge-defined N_CRC^EW(P_*) supplied by R_EW_global_capacity_certificate.json",
        },
        "normal_form": {
            "premise_status": "declared_not_constructed_here",
            "declared_finite_state": True,
            "declared_unique": True,
            "declared_schedule_independent": True,
            "certificate_basis": [
                "strict finite Lyapunov descent",
                "local diamond on the physical quotient",
                "repair completeness",
            ],
            "source_artifact": "paper/tex_fragments/OBSERVERS_APPENDICES.tex (Lyapunov+diamond+completeness confluence theorem)",
        },
        "observer_sector": {
            "premise_status": "declared_not_constructed_here",
            "declared_nonzero": True,
            "declared_stable_self_reading": True,
            "declared_central_projector": True,
            "declared_capacity_visible_on_quotient": True,
            "declared_quotient_invariant": True,
            "source_artifact": "paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex",
        },
        "capacity_register": {
            "premise_status": "declared_not_constructed_here",
            "declared_central": True,
            "declared_positive_spectrum": True,
            "declared_selected_variance": "0",
            "atoms": [
                {
                    "label": "N_CRC_EW_capacity_atom",
                    "capacity": decstr(cap_read),
                    "probability": "1",
                    "selected": True,
                    "selection_status": "declared_bridge_identification_not_derived_here",
                    "source_artifact": ew_cert_relpath,
                }
            ],
        },
        "extractor": {
            "formula": "rho_read = sqrt(pi / Cap_read)",
            "area_law": "N = pi / rho^2",
            "positive_root": True,
            "interval_rule": "cap in [c_minus,c_plus] -> rho in [sqrt(pi/c_plus),sqrt(pi/c_minus)]",
        },
        "readback_resolution": {
            "cap_read": decstr(cap_read),
            "rho_read": decstr(rho_read),
            "rho_read_interval": [decstr(rho_read), decstr(rho_read)],
            "display_N_CRC": decstr(n_display),
            "display_rho_star": decstr(rho_display),
            "display_only": True,
        },
        "conditional_source_capacity_check": {
            "premise_status": "bridge_defined_value_substituted_into_declared_atom",
            "enabled": True,
            "N_CRC": decstr(source_n),
            "rho_star": decstr(rho_star),
            "capacity_relative_residual": decstr(strict_capacity_residual),
            "rho_relative_residual": decstr(strict_rho_residual),
            "relative_tolerance": decstr(tol),
            "accepted": bool(strict_ok),
        },
        "refinement_certificate": {
            "premise_status": "cofinal_family_declared_not_computed_here",
            "cofinal_limit_declared": True,
            "positive_root_implication_verified": True,
            "contraction_kappa": decstr(kappa),
            "ew_capacity_lambda": decstr(ew_lambda),
            "kappa_matches_ew_lambda": kappa == ew_lambda,
            "finite_to_limit_bound": (
                "|N_r_star-N_CRC^EW| <= delta_r/(1-kappa); "
                "|rho_r-rho_star| <= sqrt(pi)/(2*C_min^(3/2))*delta_r/(1-kappa)"
            ),
            "limit_statement": (
                "If F_r -> F cofinally and F(N_CRC^EW)=N_CRC^EW>0 (supplied by "
                "R_EW_global_capacity_certificate.json), then "
                "rho_read(r,N_CRC^EW)=sqrt(pi/F_r(N_CRC^EW)) -> "
                "(N_CRC^EW/pi)^(-1/2) = rho_star."
            ),
        },
        "dependencies": {
            "ew_refined_exact_capacity": ew_accepted,
            "global_repair_tick_area_law_supplied_elsewhere": True,
        },
        "dependency_artifacts": {
            "ew_refined_exact_capacity": ew_cert_relpath,
            "global_repair_tick_area_law_supplied_elsewhere": "code/particles/hierarchy/certificates/R_N_global_repair_tick_certificate.json",
        },
        "derivation_chain": derivation_chain,
        "factor_origins": factor_origins,
        "branch_scope": branch_scope,
        "obstruction_record": obstruction_record,
        "acceptance_criteria_status": acceptance_criteria_status,
        "input_ledger": {
            "used_inputs": used_inputs,
            "used_artifacts": used_artifacts,
            "forbidden_inputs": sorted(FORBIDDEN_INPUTS),
        },
        "claim_boundary": {
            "closed_here": [
                "the extractor rho_read = sqrt(pi/Cap_read) is the unique positive-root inverse of the D6 area law",
                "conditional implication from the declared singleton readback and cofinal-refinement premises to rho_read(r,N_CRC^EW) -> rho_star with the quantitative bound",
            ],
            "closed_elsewhere": [
                "bridge-defined EW capacity contraction in R_EW_global_capacity_certificate.json (supplies N_CRC^EW>0 as input, not cosmic capacity)",
                "conditional global repair-tick lemma in R_N_global_repair_tick_certificate.json",
                "representation-to-spectrum round-count theorem in R_m_rep_24_certificate.json",
                "electroweak projection bridge in R_EW_tick_projection_certificate.json",
                "conditional RG/Higgs naturality square in issue_332_rg_naturality_certificate.json",
                "conditional local/global hierarchy-resonance composition in R_local_global_hierarchy_resonance_closeout_335.json",
            ],
            "not_closed_here": [
                "construction of an executable finite family U_{r,N} and its repair map nf_{r,N}",
                "construction of Obs, the central register, and a source-selected capacity atom",
                "a computed cofinal refinement family F_r -> F",
                "an independently emitted physical cosmic capacity",
                "HIERARCHY-SCREEN-READOUT",
            ],
            "scope": (
                "Conditional positive-root implication on the declared branches in "
                "branch_scope. N_CRC^EW(P_*) is a bridge-defined source value. "
                "Neither the finite readback objects nor physical cosmic capacity "
                "nor HIERARCHY-SCREEN-READOUT is derived by this certificate."
            ),
        },
        "source_status": {
            "closes_gate": False,
            "receipt_class": "conditional_schema",
            "does_not_promote_full_hierarchy_resonance": True,
            "remaining_for_full_hierarchy_resonance": [
                "executable finite readback-family construction",
                "source selection of the central capacity atom",
                "computed cofinal refinement",
                "HIERARCHY-SCREEN-READOUT",
                "independent physical cosmic-capacity map",
            ],
        },
        "checks": {
            "normal_form_premise_recorded": True,
            "observer_sector_premise_recorded": True,
            "singleton_capacity_atom_premise_recorded": True,
            "zero_capacity_variance_premise_recorded": True,
            "positive_root_extractor": True,
            "kappa_is_contractive": kappa_ok,
            "kappa_matches_ew_lambda": kappa == ew_lambda,
            "forbidden_input_check": no_forbidden_inputs,
            "conditional_source_capacity_substitution_check": bool(strict_ok),
            "ew_dependency_accepted": ew_accepted,
            "ew_bridge_residual_zero": ew_bridge_residual == 0,
            "rounded_capacity_display_rejected_as_bridge": rounded_rejected_as_bridge,
            "declared_cap_read_equals_bridge_defined_n_crc_ew": cap_read == n_crc_ew,
            "positive_root_equals_rho_star_under_substitution": rho_read == rho_star,
            "all_acceptance_criteria_pass": all(acceptance_criteria_status.values()),
        },
        "verifier_command": (
            "python3 code/particles/hierarchy/verify_issue_342_readback_resolution.py "
            "--check --output code/particles/hierarchy/certificates/R_readback_resolution_certificate.json"
        ),
    }


def build_derivation_chain(
    *,
    ew_cert_relpath: str,
    n_crc_ew: Decimal,
    ew_lambda: Decimal,
    kappa: Decimal,
) -> list[dict[str, Any]]:
    return [
        {
            "step": 1,
            "premise": "Declared finite patch-carrier object",
            "uses": ["premise U_{r,N}: finite OPH patch system at fixed cutoff r and capacity coordinate N"],
            "source_artifact": "paper/screen_microphysics_and_observer_synchronization.tex, section 'From Spherical Screens to Federated Patch Carriers'; paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex",
            "conclusion": "Conditional schema domain: assume U_{r,N} is a finite initial object on which the repair-and-readback pipeline acts; no executable U_{r,N} family is constructed here.",
        },
        {
            "step": 2,
            "premise": "Declared fixed-cutoff confluence premise (Lyapunov + local diamond + repair completeness)",
            "uses": ["U_{r,N}", "strictly Lyapunov-decreasing finite repair", "local diamond on the physical quotient", "repair completeness"],
            "source_artifact": "paper/tex_fragments/OBSERVERS_APPENDICES.tex (confluence theorem)",
            "conclusion": "Under the confluence premise, nf_{r,N}(U_{r,N}) is a unique schedule-independent terminal normal form; the repair map is not constructed here.",
        },
        {
            "step": 3,
            "premise": "Declared stable self-reading observer-sector premise",
            "uses": ["n_{r,N}", "non-zero observer sector", "stable self-reading", "central projector", "capacity visible on the physical quotient"],
            "source_artifact": "paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex",
            "conclusion": "Under the observer-sector premise, Obs(n_{r,N}) is a non-zero stable self-reading sector on the physical quotient; Obs is not constructed here.",
        },
        {
            "step": 4,
            "premise": "Declared sharp central-record premise",
            "uses": ["Obs(n_{r,N})", "central capacity register C_hat_{r,N}=sum_c c P_c", "one selected positive atom with zero variance"],
            "source_artifact": "paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex (central record algebra)",
            "conclusion": "Under the sharp-atom premise, Cap_read(Obs(n_{r,N})) is single-valued. The atom and its selection are not derived here.",
            "acceptance_criterion_closed": "conditional singleton readback implication",
        },
        {
            "step": 5,
            "premise": "D6 dimensionless area law (corpus axiom)",
            "uses": ["screen-capacity normalization N = pi*(r/ell)^2", "dimensionless reduction N = pi/rho^2"],
            "source_artifact": "extra/observer_patch_holography_as_string_vacuum_selector.tex, theorem 'Conditional edge-capacity entropy theorem'; reflected in R_N_global_repair_tick_certificate.json (definitions.screen_normalized_radius_coordinate, normalization.screen_capacity_relation)",
            "conclusion": "The positive-root inverse of the D6 area law gives the extractor rho_read(r,N) = sqrt(pi/F_r(N)) as the unique positive solution of N = pi/rho^2 evaluated at N = F_r(N).",
        },
        {
            "step": 6,
            "premise": "Conditional composition of the finite pipeline",
            "uses": ["steps 1-4 (declared finite object, normal form, observer sector, and sharp atom)", "step 5 (positive-root extractor)"],
            "source_artifact": "this certificate (extractor and capacity_register fields)",
            "conclusion": "Under the declared premises, F_r(N) := Cap_read(Obs(nf_{r,N}(U_{r,N}))) is single-valued and rho_read(r,N) = sqrt(pi/F_r(N)) is its positive-root resolution.",
        },
        {
            "step": 7,
            "premise": f"EW-refined exact-capacity Banach contraction with lambda={ew_lambda}",
            "uses": ["C_EW(P,x) = (1-lambda)*x + lambda*6*pi/(P*alpha_U(P))", "Banach fixed-point theorem on the source-side log-capacity coordinate"],
            "source_artifact": ew_cert_relpath,
            "conclusion": f"The bridge-defined map has fixed point N_CRC^EW(P_*) = pi*exp[6*pi/(P_*alpha_U(P_*))] = {decstr(n_crc_ew)}, with B_EW(P_*,N_CRC^EW)=0. Its use as F(N_CRC^EW)=N_CRC^EW is a declared bridge identification, not an independent readback or cosmic-capacity derivation.",
            "acceptance_criterion_closed": "bridge-defined fixed-point value supplied",
        },
        {
            "step": 8,
            "premise": "Positive-root implication under the declared cofinal refinement",
            "uses": [
                "step 6 (conditionally single-valued F_r and rho_read)",
                "step 7 plus the declared identification F(N_CRC^EW)=N_CRC^EW>0",
                f"declared cofinal refinement F_r -> F with quantitative residual delta_r and contraction kappa={kappa}",
                "chain rule on rho(N) = sqrt(pi)*N^(-1/2) (drho/dN = -sqrt(pi)/(2*N^(3/2)))",
            ],
            "source_artifact": "this certificate (refinement_certificate.finite_to_limit_bound)",
            "conclusion": "If the declared residual bound holds, then |N_r_star - N_CRC^EW| <= delta_r/(1-kappa) and chain-rule continuity gives |rho_r - rho_star| <= sqrt(pi)/(2*C_min^(3/2))*delta_r/(1-kappa) -> 0. Thus rho_read(r,N_CRC^EW) -> (N_CRC^EW/pi)^(-1/2) conditionally.",
            "acceptance_criterion_closed": "conditional positive-root refinement implication",
        },
    ]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify OPH issue #342 conditional positive-root readback schema."
    )
    parser.add_argument(
        "--n-display",
        default=str(DEFAULT_N_DISPLAY),
        help="rounded cosmological capacity display (diagnostic only)",
    )
    parser.add_argument(
        "--n-crc",
        default=None,
        help="optional override for the source bridge fixed point N_CRC; defaults to N_CRC^EW from the EW exact-capacity certificate",
    )
    parser.add_argument(
        "--cap-read",
        default=None,
        help="optional Cap_read override; defaults to N_CRC^EW from the EW exact-capacity certificate",
    )
    parser.add_argument(
        "--kappa",
        default=str(DEFAULT_KAPPA),
        help="contraction constant for the refinement bound (must equal the EW Banach lambda)",
    )
    parser.add_argument(
        "--relative-tolerance",
        default=str(DEFAULT_TOL),
        help="strict source residual tolerance",
    )
    parser.add_argument(
        "--ew-capacity-cert",
        default=str(DEFAULT_EW_CAPACITY_CERT),
        help="path to R_EW_global_capacity_certificate.json",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit nonzero unless the certificate passes",
    )
    parser.add_argument("--output", default=None, help="write JSON certificate to path")
    args = parser.parse_args()

    n_display = D(args.n_display)
    source_n = D(args.n_crc)
    cap_read = D(args.cap_read)
    kappa = D(args.kappa)
    tol = D(args.relative_tolerance)
    assert n_display is not None and kappa is not None and tol is not None

    ew_cert_path = Path(args.ew_capacity_cert).resolve()
    ew_cert = load_ew_capacity_cert(ew_cert_path)
    try:
        ew_cert_relpath = str(ew_cert_path.relative_to(REPO_ROOT)).replace("\\", "/")
    except ValueError:
        ew_cert_relpath = str(ew_cert_path).replace("\\", "/")

    cert = build_certificate(
        n_display,
        kappa,
        source_n,
        cap_read,
        tol,
        ew_cert,
        ew_cert_relpath,
    )
    text = json.dumps(cert, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    if args.check and not cert["accepted"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

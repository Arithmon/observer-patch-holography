#!/usr/bin/env python3
"""Audit the conditional Majorana phases of the rejected weighted-cycle lane.

Chain role: retain the canonical Takagi phase pair as a candidate-only
diagnostic and enforce that the historical shared-basis transport can never
promote it.

Mathematics:
1. The repaired weighted-cycle artifact already emits one complex symmetric
   cycle matrix and one canonical Takagi PMNS surface.
2. A readout from PMNS columns alone is not supported because arbitrary
   computational column phases move `alpha_21` / `alpha_31` while leaving the
   oscillation observables unchanged.
3. The canonical congruence gauge is fixed by the Takagi condition
   `U_PMNS^T M_nu U_PMNS = diag(m_i)` with `diag(m_i)` real and non-negative.
4. In the electron-row gauge `U_e1 in R_{>0}`, the canonical pair is then read
   out as `alpha_21 = 2 arg(U_e2)` and
   `alpha_31 = 2 (arg(U_e3) + delta_PMNS)`.
5. The historical transport defined `U_nu_shared = U_e_left U_wc`, so
   `U_e_left^dagger U_nu_shared = U_wc` is an identity and cannot validate a
   physical charged basis.
6. A future physical result must use a distinct source-derived operator/basis
   artifact; this weighted-cycle transport object is permanently non-promotable.

Declared inputs: the rejected weighted-cycle candidate, the current open
charged-lepton basis artifact, and optionally its tautological shared-basis
transport audit.

Output: one candidate-only canonical phase pair plus the exact blocker for
public promotion. No physical Majorana phase is emitted.
"""

from __future__ import annotations

import argparse
import json
import sys
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.mixing import mixing_parameters as _standard_pmns_parameters
DEFAULT_WEIGHTED_CYCLE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_SHARED_CHARGED_LEFT = ROOT / "particles" / "runs" / "neutrino" / "shared_charged_lepton_left_basis.json"
DEFAULT_SHARED_BASIS_REPRESENTATION = (
    ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_shared_basis_representation.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_physical_majorana_phase_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _canonical_artifact_ref(path: Path | None) -> str | None:
    if path is None:
        return None
    if not path.is_absolute():
        return path.as_posix()
    try:
        rel = path.relative_to(ROOT)
    except ValueError:
        return path.as_posix()
    return f"code/{rel.as_posix()}"


def _complex_matrix(payload: dict[str, Any], real_key: str, imag_key: str) -> np.ndarray:
    return np.array(payload[real_key], dtype=float) + 1j * np.array(payload[imag_key], dtype=float)


def _canonicalize_takagi_column_signs(unitary: np.ndarray) -> np.ndarray:
    canonical = unitary.copy()
    for column in range(canonical.shape[1]):
        for row in range(canonical.shape[0]):
            entry = canonical[row, column]
            if abs(entry) <= 1.0e-12:
                continue
            if entry.real < -1.0e-12 or (abs(entry.real) <= 1.0e-12 and entry.imag < -1.0e-12):
                canonical[:, column] *= -1.0
            break
    return canonical


def _wrap_signed_phase(angle_rad: float) -> float:
    wrapped = float((angle_rad + math.pi) % (2.0 * math.pi) - math.pi)
    if wrapped <= -math.pi:
        return wrapped + 2.0 * math.pi
    return wrapped


def _row_gauge_pmns(unitary: np.ndarray) -> tuple[np.ndarray, float]:
    electron_row_phase = float(np.angle(unitary[0, 0]))
    row_gauged = np.diag([np.exp(-1j * electron_row_phase), 1.0, 1.0]) @ unitary
    u_e1 = row_gauged[0, 0]
    if abs(float(np.imag(u_e1))) > 1.0e-12:
        raise ValueError("electron-row Majorana readout gauge must make U_e1 real")
    if float(np.real(u_e1)) <= 0.0:
        raise ValueError("electron-row Majorana readout gauge must make U_e1 positive")
    return row_gauged, electron_row_phase


def _majorana_readout_details(unitary: np.ndarray, delta_rad: float) -> dict[str, Any]:
    row_gauged, electron_row_phase = _row_gauge_pmns(unitary)
    alpha21_signed = _wrap_signed_phase(2.0 * float(np.angle(row_gauged[0, 1])))
    alpha31_signed = _wrap_signed_phase(2.0 * (float(np.angle(row_gauged[0, 2])) + delta_rad))
    alpha21_mod = float(alpha21_signed % (2.0 * math.pi))
    alpha31_mod = float(alpha31_signed % (2.0 * math.pi))
    return {
        "row_gauged_pmns": row_gauged,
        "checks": {
            "row_gauged_u_e1_real": float(np.real(row_gauged[0, 0])),
            "row_gauged_u_e1_imag_abs": abs(float(np.imag(row_gauged[0, 0]))),
        },
        "parameters": {
            "electron_row_gauge_phase_rad": electron_row_phase,
            "electron_row_gauge_phase_deg": math.degrees(electron_row_phase),
            "alpha21_rad": alpha21_signed,
            "alpha21_deg": math.degrees(alpha21_signed),
            "alpha21_rad_0_to_2pi": alpha21_mod,
            "alpha21_deg_0_to_360": math.degrees(alpha21_mod),
            "alpha31_rad": alpha31_signed,
            "alpha31_deg": math.degrees(alpha31_signed),
            "alpha31_rad_0_to_2pi": alpha31_mod,
            "alpha31_deg_0_to_360": math.degrees(alpha31_mod),
        },
    }


def _majorana_pair_from_pmns(unitary: np.ndarray, delta_rad: float) -> dict[str, float]:
    return dict(_majorana_readout_details(unitary, delta_rad)["parameters"])


def _canonical_pmns_from_weighted_cycle_matrix(matrix: np.ndarray) -> dict[str, Any]:
    if np.max(np.abs(matrix - matrix.T)) > 1.0e-12:
        raise ValueError("weighted-cycle matrix must remain complex symmetric for Majorana readout")

    left_vectors, singular_values, vh = np.linalg.svd(matrix)
    order = np.argsort(singular_values)
    left_vectors = left_vectors[:, order]
    singular_values = singular_values[order]
    right_vectors = vh.conj().T[:, order]

    phase_relation = left_vectors.T @ right_vectors
    offdiag_relation = phase_relation - np.diag(np.diag(phase_relation))
    relation_offdiag_max = float(np.max(np.abs(offdiag_relation)))
    if relation_offdiag_max > 1.0e-8:
        raise ValueError("SVD phase relation is not diagonal enough to define a canonical Takagi gauge")

    congruence_half_angles = np.angle(np.diag(phase_relation))
    takagi_vectors = left_vectors @ np.diag(np.exp(-0.5j * congruence_half_angles))
    pmns = np.conjugate(takagi_vectors)
    pmns = _canonicalize_takagi_column_signs(pmns)
    diagonalized = pmns.T @ matrix @ pmns

    diagonalized_offdiag = diagonalized - np.diag(np.diag(diagonalized))
    diag_real = np.real(np.diag(diagonalized))
    if np.any(diag_real <= 0.0):
        raise ValueError("canonical Takagi readout must produce a positive diagonal mass row")

    return {
        "pmns": pmns,
        "singular_values": [float(x) for x in singular_values.tolist()],
        "congruence_half_angles_rad": [float(0.5 * angle) for angle in congruence_half_angles.tolist()],
        "congruence_half_angles_deg": [math.degrees(float(0.5 * angle)) for angle in congruence_half_angles.tolist()],
        "phase_relation_offdiag_max_abs": relation_offdiag_max,
        "diagonalized_real_masses": [float(x) for x in diag_real.tolist()],
        "diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(diagonalized))))),
        "diagonalized_offdiag_max_abs": float(np.max(np.abs(diagonalized_offdiag))),
    }


def _observable_match(observables: dict[str, float], weighted_cycle_observables: dict[str, Any]) -> dict[str, float]:
    return {
        "theta12_deg_abs_delta": abs(observables["theta12_deg"] - float(weighted_cycle_observables["theta12_deg"])),
        "theta23_deg_abs_delta": abs(observables["theta23_deg"] - float(weighted_cycle_observables["theta23_deg"])),
        "theta13_deg_abs_delta": abs(observables["theta13_deg"] - float(weighted_cycle_observables["theta13_deg"])),
        "delta_deg_abs_delta": abs(((observables["delta_deg"] - float(weighted_cycle_observables["delta_deg"]) + 180.0) % 360.0) - 180.0),
        "J_abs_delta": abs(observables["J"] - float(weighted_cycle_observables["J"])),
    }


def _shared_basis_transport_checks(
    shared_basis_matrix: np.ndarray,
    u_nu_shared: np.ndarray,
    pmns: np.ndarray,
    u_e_left: np.ndarray,
) -> dict[str, Any]:
    recovered_pmns = np.conjugate(u_e_left).T @ u_nu_shared
    shared_diagonalized = u_nu_shared.T @ shared_basis_matrix @ u_nu_shared
    shared_offdiag = shared_diagonalized - np.diag(np.diag(shared_diagonalized))
    shared_diag_real = np.real(np.diag(shared_diagonalized))
    if np.any(shared_diag_real <= 0.0):
        raise ValueError("shared-basis weighted-cycle representation must keep a positive Takagi diagonal")
    return {
        "shared_basis_symmetry_max_abs": float(np.max(np.abs(shared_basis_matrix - shared_basis_matrix.T))),
        "shared_basis_diagonalized_offdiag_max_abs": float(np.max(np.abs(shared_offdiag))),
        "shared_basis_diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(shared_diagonalized))))),
        "shared_basis_diagonalized_real_masses": [float(x) for x in shared_diag_real.tolist()],
        "pmns_recovery_max_abs": float(np.max(np.abs(recovered_pmns - pmns))),
    }


def _shared_basis_branch_identity_checks(
    weighted_cycle_matrix: np.ndarray,
    canonical_pmns: np.ndarray,
    charged_basis_matrix: np.ndarray,
    shared_basis_matrix: np.ndarray,
    u_nu_shared: np.ndarray,
    pmns: np.ndarray,
    u_e_left: np.ndarray,
) -> dict[str, float]:
    expected_shared_basis_matrix = u_e_left.conj() @ weighted_cycle_matrix @ u_e_left.conj().T
    expected_u_nu_shared = u_e_left @ canonical_pmns
    expected_pmns = np.conjugate(u_e_left).T @ expected_u_nu_shared
    return {
        "charged_basis_matrix_max_abs_delta": float(np.max(np.abs(charged_basis_matrix - weighted_cycle_matrix))),
        "shared_basis_matrix_max_abs_delta": float(np.max(np.abs(shared_basis_matrix - expected_shared_basis_matrix))),
        "u_nu_shared_max_abs_delta": float(np.max(np.abs(u_nu_shared - expected_u_nu_shared))),
        "pmns_matrix_max_abs_delta": float(np.max(np.abs(pmns - expected_pmns))),
    }


def build_payload(
    weighted_cycle: dict[str, Any],
    shared_charged_left: dict[str, Any],
    shared_basis_representation: dict[str, Any] | None = None,
    *,
    weighted_cycle_path: Path,
    shared_charged_left_path: Path,
    shared_basis_representation_path: Path | None = None,
) -> dict[str, Any]:
    basis_contract = dict(shared_charged_left.get("basis_contract", {}))
    if not isinstance(shared_charged_left.get("U_e_left"), dict):
        raise ValueError("shared charged-lepton artifact must carry a diagnostic U_e_left matrix")
    if not basis_contract.get("orientation_preserved", False):
        raise ValueError("shared charged-lepton basis must preserve orientation")
    if weighted_cycle.get("pmns_row_order_for_pdg") != ["e", "mu", "tau"]:
        raise ValueError("weighted-cycle branch must already be in PDG row order")

    charged_basis_eligible = (
        shared_charged_left.get("status") == "closed"
        and shared_charged_left.get("pmns_use_allowed") is True
        and basis_contract.get("physical_identification_closed") is True
    )
    input_contract_would_have_been_eligible = (
        weighted_cycle.get("source_only_prediction_eligible") is True
        and weighted_cycle.get("prediction_promotion_allowed") is True
        and weighted_cycle.get("historical_target_exposure") is False
        and charged_basis_eligible
    )
    # This artifact cannot close the physical-basis problem.  Its historical
    # "shared-basis" construction transports both M_nu and U_nu with the same
    # chosen U_e, so U_e^dagger U_nu_shared = U_weighted_cycle identically.
    # Treating that cancellation as evidence for the charged-lepton basis was
    # the circular step uncovered by the July 2026 audit.  A future physical
    # promotion must consume a distinct source-derived operator/basis artifact,
    # not revive this weighted-cycle transport object by toggling status flags.
    candidate_promotion_eligible = False

    matrix = _complex_matrix(weighted_cycle, "repaired_cycle_matrix_real", "repaired_cycle_matrix_imag")
    canonical = _canonical_pmns_from_weighted_cycle_matrix(matrix)
    pmns = canonical["pmns"]
    observables = _standard_pmns_parameters(pmns)
    majorana_readout = _majorana_readout_details(pmns, observables["delta_rad"])
    majorana_pair = dict(majorana_readout["parameters"])
    weighted_cycle_observables = dict(weighted_cycle["pmns_observables"])
    observable_match = _observable_match(observables, weighted_cycle_observables)
    if max(observable_match.values()) > 1.0e-8:
        raise ValueError("canonical Takagi PMNS readout must agree with the emitted weighted-cycle oscillation observables")

    public_promotion_status = "blocked_rejected_candidate_and_missing_basis_placement"
    public_promotion_blocker = (
        "The canonical Takagi readout is stable only conditional on the rejected target-informed weighted-cycle "
        "matrix. The corpus does not derive the placement of that f-labelled operator in the charged-lepton mass "
        "basis, and the historical shared-basis recovery is tautological."
    )
    shared_basis_use_status = "basis_placement_open_no_public_majorana_emission"
    shared_basis_representation_summary = None
    status = "candidate_only"
    proof_chain_role = "candidate_surface"
    public_surface_candidate_allowed = False
    theorem_object = "canonical_majorana_phase_pair_candidate"
    theorem_surface = "weighted_cycle_takagi_readout_candidate"
    statement = (
        "Canonical Takagi congruence fixes one Majorana-phase pair conditional on the rejected weighted-cycle "
        "matrix. The pair has comparison-only candidate status because the physical basis placement is not derived."
    )
    emitted_parameters: dict[str, Any] | None = None
    takagi_congruence_payload = {
        "singular_values": canonical["singular_values"],
        "congruence_half_angles_rad": canonical["congruence_half_angles_rad"],
        "congruence_half_angles_deg": canonical["congruence_half_angles_deg"],
        "phase_relation_offdiag_max_abs": canonical["phase_relation_offdiag_max_abs"],
        "diagonalized_real_masses": canonical["diagonalized_real_masses"],
        "diagonalized_imag_max_abs": canonical["diagonalized_imag_max_abs"],
        "diagonalized_offdiag_max_abs": canonical["diagonalized_offdiag_max_abs"],
    }
    readout_checks = dict(majorana_readout["checks"])
    row_gauged_pmns = majorana_readout["row_gauged_pmns"]
    shared_basis_identity_checks: dict[str, float] | None = None

    if shared_basis_representation is not None:
        if shared_basis_representation.get("artifact") != "oph_neutrino_weighted_cycle_shared_basis_representation":
            raise ValueError("unexpected shared-basis weighted-cycle artifact id")
        shared_basis_representation_summary = {
            "artifact": shared_basis_representation.get("artifact"),
            "status": shared_basis_representation.get("status"),
            "physical_branch_closed": bool(shared_basis_representation.get("physical_branch_closed", False)),
            "basis_placement_source_derived": bool(
                shared_basis_representation.get("basis_placement_source_derived", False)
            ),
            "statement": shared_basis_representation.get("statement"),
            "promotion_use": "blocked",
            "blocker": "tautological_common_transport_does_not_derive_physical_basis_placement",
        }
        if shared_basis_representation_path is not None:
            shared_basis_representation_summary["source_path"] = _canonical_artifact_ref(
                shared_basis_representation_path
            )

    return {
        "artifact": "oph_neutrino_physical_majorana_phase_theorem",
        "generated_utc": _timestamp(),
        "status": status,
        "proof_chain_role": proof_chain_role,
        "public_surface_candidate_allowed": public_surface_candidate_allowed,
        "candidate_promotion_eligible": candidate_promotion_eligible,
        "input_contract_would_have_been_eligible": input_contract_would_have_been_eligible,
        "charged_basis_eligible": charged_basis_eligible,
        "theorem_object": theorem_object,
        "theorem_surface": theorem_surface,
        "statement": statement,
        "readout_convention": {
            "stored_pmns_matrix": "pmns_matrix_* stores the canonical Takagi unitary on the repaired weighted-cycle branch before any charged-lepton row rephasing for Majorana readout. On the promoted path, the same matrix is revalidated against the explicit shared-basis transport data.",
            "takagi_condition": "U_PMNS^T M_nu U_PMNS = diag(m_i) with diag(m_i) in R_{>0}",
            "majorana_readout_row_gauge": "row_gauged_pmns_matrix_* = diag(exp(-i arg(U_e1)), 1, 1) * pmns_matrix_* so (row_gauged U)_{e1} in R_{>0}",
            "majorana_readout_formula": {
                "alpha21": "2 arg((row_gauged U)_{e2})",
                "alpha31": "2 (arg((row_gauged U)_{e3}) + delta_PMNS)",
            },
        },
        "public_promotion_status": public_promotion_status,
        "public_promotion_blocker": public_promotion_blocker,
        "shared_basis_use_status": shared_basis_use_status,
        "source_artifacts": {
            "weighted_cycle_branch": _canonical_artifact_ref(weighted_cycle_path),
            "shared_charged_left_basis": _canonical_artifact_ref(shared_charged_left_path),
            "shared_basis_representation": _canonical_artifact_ref(shared_basis_representation_path),
        },
        "depends_on": [
            dependency
            for dependency in [
                "oph_neutrino_weighted_cycle_repair",
                "oph_shared_charged_lepton_left_basis",
                "oph_neutrino_weighted_cycle_shared_basis_representation" if shared_basis_representation is not None else None,
            ]
            if dependency is not None
        ],
        "basis_contract": basis_contract,
        "takagi_congruence": takagi_congruence_payload,
        "pmns_matrix_real": np.real(pmns).tolist(),
        "pmns_matrix_imag": np.imag(pmns).tolist(),
        "row_gauged_pmns_matrix_real": np.real(row_gauged_pmns).tolist(),
        "row_gauged_pmns_matrix_imag": np.imag(row_gauged_pmns).tolist(),
        "pmns_observables": observables,
        "readout_checks": readout_checks,
        "weighted_cycle_observables_match": observable_match,
        "shared_basis_identity_checks": shared_basis_identity_checks,
        "shared_basis_representation": shared_basis_representation_summary,
        "candidate_parameters": majorana_pair,
        "emitted_parameters": emitted_parameters,
        "notes": [
            "A strict readout from PMNS columns alone would remain computational-gauge dependent, because arbitrary intermediate column phases leave the oscillation observables unchanged while moving alpha21 and alpha31.",
            "The readout fixes that ambiguity by first taking the canonical Takagi congruence of the emitted symmetric weighted-cycle matrix and then applying the readout-only electron-row gauge `U_e1 in R_{>0}`.",
            (
                "The resulting pair is promoted only after a source-eligible weighted-cycle branch and an independently derived shared-basis placement both close."
                if status == "theorem_grade_emitted"
                else "The resulting pair is a comparison-only coordinate of the rejected candidate. The shared-basis identity does not supply an independent physical-basis derivation."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Majorana-phase theorem artifact.")
    parser.add_argument("--weighted-cycle", default=str(DEFAULT_WEIGHTED_CYCLE))
    parser.add_argument("--shared-charged-left", default=str(DEFAULT_SHARED_CHARGED_LEFT))
    parser.add_argument("--shared-basis-representation", default=str(DEFAULT_SHARED_BASIS_REPRESENTATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    weighted_cycle_path = Path(args.weighted_cycle)
    shared_charged_left_path = Path(args.shared_charged_left)
    shared_basis_representation_path = Path(args.shared_basis_representation)
    payload = build_payload(
        _load_json(weighted_cycle_path),
        _load_json(shared_charged_left_path),
        _load_json(shared_basis_representation_path) if shared_basis_representation_path.exists() else None,
        weighted_cycle_path=weighted_cycle_path,
        shared_charged_left_path=shared_charged_left_path,
        shared_basis_representation_path=shared_basis_representation_path if shared_basis_representation_path.exists() else None,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Audit the weighted-cycle placement relative to the shared flavor basis.

The historical implementation called ``M_shared = U_e^* M_wc U_e^dagger`` a
shared-basis validation and then defined ``U_nu_shared = U_e U_wc``.  The
reported identity ``U_e^dagger U_nu_shared = U_wc`` follows algebraically from
that definition, so it supplies no independent basis information.

This builder preserves that change-of-coordinates calculation as a diagnostic
and also evaluates the literal source-basis interpretation.  For the latter it
reorders the declared ``[f3, f1, f2]`` cycle matrix into the independently
declared ``[f1, f2, f3]`` basis, diagonalizes it, and forms the actual mismatch
matrix ``U_PMNS = U_e^dagger U_nu``.  Neither interpretation is promoted because
the repository contains no source-derived map placing the weighted-cycle
operator in the charged-lepton mass basis.
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
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_shared_basis_representation.json"


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


def _majorana_pair_from_pmns(unitary: np.ndarray, delta_rad: float) -> dict[str, float]:
    electron_row_phase = float(np.angle(unitary[0, 0]))
    row_gauged = np.diag([np.exp(-1j * electron_row_phase), 1.0, 1.0]) @ unitary
    alpha21_signed = _wrap_signed_phase(2.0 * float(np.angle(row_gauged[0, 1])))
    alpha31_signed = _wrap_signed_phase(2.0 * (float(np.angle(row_gauged[0, 2])) + delta_rad))
    alpha21_mod = float(alpha21_signed % (2.0 * math.pi))
    alpha31_mod = float(alpha31_signed % (2.0 * math.pi))
    return {
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
    }


def _canonical_takagi_unitary(matrix: np.ndarray) -> dict[str, Any]:
    if np.max(np.abs(matrix - matrix.T)) > 1.0e-12:
        raise ValueError("weighted-cycle matrix must remain complex symmetric")

    left_vectors, singular_values, vh = np.linalg.svd(matrix)
    order = np.argsort(singular_values)
    left_vectors = left_vectors[:, order]
    singular_values = singular_values[order]
    right_vectors = vh.conj().T[:, order]

    phase_relation = left_vectors.T @ right_vectors
    offdiag_relation = phase_relation - np.diag(np.diag(phase_relation))
    relation_offdiag_max = float(np.max(np.abs(offdiag_relation)))
    if relation_offdiag_max > 1.0e-8:
        raise ValueError("SVD phase relation is not diagonal enough to define the canonical Takagi gauge")

    congruence_half_angles = np.angle(np.diag(phase_relation))
    takagi_vectors = left_vectors @ np.diag(np.exp(-0.5j * congruence_half_angles))
    unitary = np.conjugate(takagi_vectors)
    unitary = _canonicalize_takagi_column_signs(unitary)
    diagonalized = unitary.T @ matrix @ unitary
    diagonalized_offdiag = diagonalized - np.diag(np.diag(diagonalized))
    diagonal_masses = np.real(np.diag(diagonalized))
    if np.any(diagonal_masses <= 0.0):
        raise ValueError("canonical Takagi diagonal must stay positive")

    return {
        "unitary": unitary,
        "singular_values": [float(x) for x in singular_values.tolist()],
        "congruence_half_angles_rad": [float(0.5 * angle) for angle in congruence_half_angles.tolist()],
        "congruence_half_angles_deg": [math.degrees(float(0.5 * angle)) for angle in congruence_half_angles.tolist()],
        "phase_relation_offdiag_max_abs": relation_offdiag_max,
        "diagonalized_real_masses": [float(x) for x in diagonal_masses.tolist()],
        "diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(diagonalized))))),
        "diagonalized_offdiag_max_abs": float(np.max(np.abs(diagonalized_offdiag))),
    }


def build_payload(
    weighted_cycle: dict[str, Any],
    shared_charged_left: dict[str, Any],
    *,
    weighted_cycle_path: Path,
    shared_charged_left_path: Path,
) -> dict[str, Any]:
    basis_contract = dict(shared_charged_left.get("basis_contract", {}))
    if not isinstance(shared_charged_left.get("U_e_left"), dict):
        raise ValueError("shared charged-lepton artifact must carry a diagnostic U_e_left matrix")
    if not basis_contract.get("orientation_preserved", False):
        raise ValueError("shared charged-lepton basis must preserve orientation")
    weighted_cycle_matrix = _complex_matrix(weighted_cycle, "repaired_cycle_matrix_real", "repaired_cycle_matrix_imag")
    takagi = _canonical_takagi_unitary(weighted_cycle_matrix)
    declared_pmns = takagi["unitary"]
    weighted_cycle_observables = dict(weighted_cycle["pmns_observables"])

    u_e_left = _complex_matrix(shared_charged_left["U_e_left"], "real", "imag")
    # Historical construction.  This is a reversible change of coordinates of
    # an operator already assumed to be in the charged-lepton mass basis.
    shared_basis_matrix = u_e_left.conj() @ weighted_cycle_matrix @ u_e_left.conj().T
    u_nu_shared = u_e_left @ declared_pmns
    recovered_pmns = np.conjugate(u_e_left).T @ u_nu_shared

    shared_diagonalized = u_nu_shared.T @ shared_basis_matrix @ u_nu_shared
    shared_offdiag = shared_diagonalized - np.diag(np.diag(shared_diagonalized))
    shared_diag_real = np.real(np.diag(shared_diagonalized))
    if np.any(shared_diag_real <= 0.0):
        raise ValueError("shared-basis transport must preserve a positive Takagi diagonal")

    declared_observables = _standard_pmns_parameters(recovered_pmns)
    observable_match = {
        "theta12_deg_abs_delta": abs(declared_observables["theta12_deg"] - float(weighted_cycle_observables["theta12_deg"])),
        "theta23_deg_abs_delta": abs(declared_observables["theta23_deg"] - float(weighted_cycle_observables["theta23_deg"])),
        "theta13_deg_abs_delta": abs(declared_observables["theta13_deg"] - float(weighted_cycle_observables["theta13_deg"])),
        "delta_deg_abs_delta": abs(((declared_observables["delta_deg"] - float(weighted_cycle_observables["delta_deg"]) + 180.0) % 360.0) - 180.0),
        "J_abs_delta": abs(declared_observables["J"] - float(weighted_cycle_observables["J"])),
    }
    if max(observable_match.values()) > 1.0e-8:
        raise ValueError("shared-basis transported branch must recover the weighted-cycle PMNS observables exactly")

    declared_majorana_pair = _majorana_pair_from_pmns(recovered_pmns, declared_observables["delta_rad"])
    transport_checks = {
        "shared_basis_symmetry_max_abs": float(np.max(np.abs(shared_basis_matrix - shared_basis_matrix.T))),
        "shared_basis_diagonalized_offdiag_max_abs": float(np.max(np.abs(shared_offdiag))),
        "shared_basis_diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(shared_diagonalized))))),
        "shared_basis_diagonalized_real_masses": [float(x) for x in shared_diag_real.tolist()],
        "pmns_recovery_max_abs": float(np.max(np.abs(recovered_pmns - declared_pmns))),
    }

    cycle_basis_order = list(weighted_cycle.get("cycle_basis_order") or [])
    shared_basis_order = list(basis_contract.get("labels") or shared_charged_left.get("labels") or [])
    if sorted(cycle_basis_order) != sorted(shared_basis_order) or len(shared_basis_order) != 3:
        raise ValueError("cycle and shared basis labels must describe the same three labels")
    shared_from_cycle_indices = [cycle_basis_order.index(label) for label in shared_basis_order]
    source_basis_matrix = weighted_cycle_matrix[np.ix_(shared_from_cycle_indices, shared_from_cycle_indices)]
    source_takagi = _canonical_takagi_unitary(source_basis_matrix)
    source_u_nu = source_takagi["unitary"]
    source_pmns = np.conjugate(u_e_left).T @ source_u_nu
    source_observables = _standard_pmns_parameters(source_pmns)

    placement_contract = dict(weighted_cycle.get("basis_placement_contract") or {})
    placement_is_source_derived = (
        placement_contract.get("status") == "source_derived"
        and placement_contract.get("operator_input_basis") == shared_basis_order
        and placement_contract.get("physical_charged_basis_map") == "U_e_left"
    )

    return {
        "artifact": "oph_neutrino_weighted_cycle_shared_basis_representation",
        "generated_utc": _timestamp(),
        "status": "basis_placement_open_tautological_transport_audit",
        "proof_chain_role": "nonpromoting_basis_audit",
        "theorem_object": None,
        "theorem_surface": None,
        "public_surface_candidate_allowed": False,
        "physical_branch_closed": False,
        "prediction_promotion_allowed": False,
        "basis_placement_source_derived": placement_is_source_derived,
        "statement": (
            "The historical shared-basis recovery is an identity created by defining U_nu_shared = U_e_left U_wc. "
            "No source theorem places the f-labelled weighted-cycle operator in the charged-lepton mass basis, so the "
            "calculation cannot close a physical PMNS branch."
        ),
        "construction": {
            "historical_declared_charged_basis_pullback": "M_shared = U_e_left^* M_wc U_e_left^dagger",
            "historical_defined_neutrino_unitary": "U_nu_shared = U_e_left U_wc",
            "historical_identity": "U_e_left^dagger U_nu_shared = U_wc",
            "literal_source_basis_test": "reorder M_wc from cycle_basis_order to [f1,f2,f3], diagonalize independently, then form U_e_left^dagger U_nu",
        },
        "basis_contract": basis_contract,
        "basis_labels": list(shared_charged_left.get("labels") or []),
        "charged_basis_input_status": shared_charged_left.get("status"),
        "charged_basis_pmns_use_allowed": bool(shared_charged_left.get("pmns_use_allowed", False)),
        "source_artifacts": {
            "weighted_cycle_branch": _canonical_artifact_ref(weighted_cycle_path),
            "shared_charged_left_basis": _canonical_artifact_ref(shared_charged_left_path),
        },
        "depends_on": [
            "oph_neutrino_weighted_cycle_repair",
            "oph_shared_charged_lepton_left_basis",
        ],
        "weighted_cycle_takagi_congruence": {
            "singular_values": takagi["singular_values"],
            "congruence_half_angles_rad": takagi["congruence_half_angles_rad"],
            "congruence_half_angles_deg": takagi["congruence_half_angles_deg"],
            "phase_relation_offdiag_max_abs": takagi["phase_relation_offdiag_max_abs"],
            "diagonalized_real_masses": takagi["diagonalized_real_masses"],
            "diagonalized_imag_max_abs": takagi["diagonalized_imag_max_abs"],
            "diagonalized_offdiag_max_abs": takagi["diagonalized_offdiag_max_abs"],
        },
        "declared_charged_basis_matrix_real": np.real(weighted_cycle_matrix).tolist(),
        "declared_charged_basis_matrix_imag": np.imag(weighted_cycle_matrix).tolist(),
        "charged_basis_matrix_real": np.real(weighted_cycle_matrix).tolist(),
        "charged_basis_matrix_imag": np.imag(weighted_cycle_matrix).tolist(),
        "shared_basis_matrix_real": np.real(shared_basis_matrix).tolist(),
        "shared_basis_matrix_imag": np.imag(shared_basis_matrix).tolist(),
        "u_nu_shared_real": np.real(u_nu_shared).tolist(),
        "u_nu_shared_imag": np.imag(u_nu_shared).tolist(),
        "pmns_matrix_real": np.real(recovered_pmns).tolist(),
        "pmns_matrix_imag": np.imag(recovered_pmns).tolist(),
        "pmns_observables": declared_observables,
        "weighted_cycle_observables_match": observable_match,
        "transport_checks": transport_checks,
        "emitted_parameters": None,
        "candidate_parameters": declared_majorana_pair,
        "basis_audit": {
            "cycle_basis_order": cycle_basis_order,
            "shared_basis_order": shared_basis_order,
            "shared_from_cycle_indices": shared_from_cycle_indices,
            "historical_recovery_is_tautology": True,
            "historical_recovery_independent_empirical_content": False,
            "missing_source_object": "weighted_cycle_operator_basis_placement_and_physical_family_map",
            "literal_source_basis_interpretation_is_promotable": False,
            "literal_source_basis_matrix_real": np.real(source_basis_matrix).tolist(),
            "literal_source_basis_matrix_imag": np.imag(source_basis_matrix).tolist(),
            "literal_source_basis_u_nu_real": np.real(source_u_nu).tolist(),
            "literal_source_basis_u_nu_imag": np.imag(source_u_nu).tolist(),
            "literal_source_basis_pmns_real": np.real(source_pmns).tolist(),
            "literal_source_basis_pmns_imag": np.imag(source_pmns).tolist(),
            "literal_source_basis_pmns_abs": np.abs(source_pmns).tolist(),
            "literal_source_basis_pmns_observables": source_observables,
        },
        "notes": [
            "The exact congruence and PMNS recovery residuals verify algebra only; they do not select the input basis.",
            "The literal f-basis calculation is a diagnostic because the current corpus does not derive the physical family map or the weighted-cycle operator itself.",
            "A future promotion requires an independently source-emitted neutrino operator and charged-family basis map before oscillation data are consulted.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the shared-basis representation of the repaired weighted-cycle branch.")
    parser.add_argument("--weighted-cycle", default=str(DEFAULT_WEIGHTED_CYCLE))
    parser.add_argument("--shared-charged-left", default=str(DEFAULT_SHARED_CHARGED_LEFT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    weighted_cycle_path = Path(args.weighted_cycle)
    shared_charged_left_path = Path(args.shared_charged_left)
    payload = build_payload(
        _load_json(weighted_cycle_path),
        _load_json(shared_charged_left_path),
        weighted_cycle_path=weighted_cycle_path,
        shared_charged_left_path=shared_charged_left_path,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

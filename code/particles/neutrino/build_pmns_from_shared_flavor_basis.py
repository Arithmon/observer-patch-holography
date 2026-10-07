#!/usr/bin/env python3
"""Build PMNS from a shared charged-lepton basis on the same family labels."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.mixing import pmns_signed as _standard_pmns_parameters

from particles.artifact_paths import portable_json_dumps


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_complex_matrix(real_rows: list[list[float]], imag_rows: list[list[float]]) -> np.ndarray:
    return np.array(real_rows, dtype=float) + 1j * np.array(imag_rows, dtype=float)


def _basis_labels(payload: dict[str, Any]) -> list[str] | None:
    labels = payload.get("labels")
    if isinstance(labels, list) and labels:
        return [str(item) for item in labels]
    contract = payload.get("same_label_basis_contract")
    if isinstance(contract, dict) and isinstance(contract.get("labels"), list):
        return [str(item) for item in contract["labels"]]
    row_labels = payload.get("row_basis_labels")
    if isinstance(row_labels, list) and row_labels:
        return [str(item) for item in row_labels]
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Build PMNS from a shared charged-lepton basis if available.")
    ap.add_argument("--majorana", default="particles/runs/neutrino/intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json")
    ap.add_argument(
        "--charged-left",
        default="particles/runs/neutrino/shared_charged_lepton_left_basis.json",
        help="JSON artifact containing charged-lepton left diagonalizer",
    )
    ap.add_argument("--out", default="particles/runs/neutrino/pmns_from_shared_basis.json")
    args = ap.parse_args()

    majorana_path = pathlib.Path(args.majorana)
    if not majorana_path.is_absolute():
        majorana_path = ROOT / args.majorana
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    majorana = load_json(majorana_path)
    neutrino_source_closed = (
        majorana.get("source_only_physical_input_eligible") is True
        and (majorana.get("source_closure_status") or {}).get("closed") is True
    )
    if "U_nu_real" in majorana:
        u_nu = load_complex_matrix(majorana["U_nu_real"], majorana["U_nu_imag"])
    else:
        u_nu = load_complex_matrix(majorana["u_nu_real"], majorana["u_nu_imag"])
    payload: dict[str, Any] = {
        "status": "conditional_blocked",
        "majorana_artifact": str(majorana_path),
        "notes": [
            "Physical PMNS requires a forward-derived charged-lepton left basis.",
        ],
    }
    if args.charged_left:
        charged_path = pathlib.Path(args.charged_left)
        if not charged_path.is_absolute():
            charged_path = ROOT / args.charged_left
        if charged_path.exists():
            charged = load_json(charged_path)
        else:
            charged = {}
        basis_contract = dict(charged.get("basis_contract", {}))
        charged_labels = _basis_labels(charged)
        neutrino_labels = _basis_labels(majorana)
        charged_basis_closed = (
            charged.get("status") == "closed"
            and charged.get("pmns_use_allowed") is True
            and basis_contract.get("physical_identification_closed") is True
        )
        if neutrino_labels is None:
            neutrino_labels = charged_labels
        if charged_labels != neutrino_labels or not basis_contract.get("orientation_preserved", False):
            payload = {
                "status": "blocked_basis_mismatch",
                "majorana_artifact": str(majorana_path),
                "charged_left_artifact": str(charged_path),
                "notes": [
                    "Refusing PMNS because the charged-lepton basis contract is missing or mismatched.",
                ],
            }
        elif not charged_basis_closed:
            u_e = load_complex_matrix(charged["U_e_left"]["real"], charged["U_e_left"]["imag"])
            diagnostic_pmns = np.conjugate(u_e).T @ u_nu
            payload = {
                "status": "blocked_upstream_charged_basis_open",
                "majorana_artifact": str(majorana_path),
                "charged_left_artifact": str(charged_path),
                "basis_labels": neutrino_labels,
                "diagnostic_only": {
                    "status": "conditional_on_open_charged_basis",
                    "pmns_abs": np.abs(diagnostic_pmns).tolist(),
                    "pmns_real": np.real(diagnostic_pmns).tolist(),
                    "pmns_imag": np.imag(diagnostic_pmns).tolist(),
                    "standard_pmns_parameters": _standard_pmns_parameters(diagnostic_pmns),
                },
                "notes": [
                    "Refusing PMNS because the charged-lepton artifact does not close a stable physical left basis.",
                    (
                        "The neutrino scalar input is also not source-closed."
                        if not neutrino_source_closed
                        else "The neutrino scalar input passes its source-closure gate."
                    ),
                ],
            }
        elif not neutrino_source_closed:
            payload = {
                "status": "blocked_upstream_neutrino_source_open",
                "majorana_artifact": str(majorana_path),
                "charged_left_artifact": str(charged_path),
                "basis_labels": neutrino_labels,
                "notes": [
                    "Refusing PMNS because the intrinsic neutrino scalar input inherits template/candidate source artifacts.",
                ],
            }
        else:
            u_e = load_complex_matrix(charged["U_e_left"]["real"], charged["U_e_left"]["imag"])
            pmns = np.conjugate(u_e).T @ u_nu
            closed = bool(
                charged_basis_closed
                and neutrino_source_closed
                and (
                    majorana.get("source_scalar_certificate")
                    or majorana.get("completion_scope") == "intrinsic_mass_eigenstates_only"
                )
            )
            payload = {
                "status": "closed" if closed else "conditional_pmns",
                "majorana_artifact": str(majorana_path),
                "charged_left_artifact": str(charged_path),
                "basis_labels": neutrino_labels,
                "orientation_preserved": True,
                "pmns_abs": np.abs(pmns).tolist(),
                "pmns_real": np.real(pmns).tolist(),
                "pmns_imag": np.imag(pmns).tolist(),
                "standard_pmns_parameters": _standard_pmns_parameters(pmns),
                "notes": [
                    (
                        "This artifact is exact on the live intrinsic neutrino bundle and shared charged-lepton basis."
                        if closed
                        else "This artifact is conditional on the charged-lepton left basis being forward-derived."
                    ),
                ],
            }
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

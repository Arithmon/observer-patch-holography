#!/usr/bin/env python3
"""Enumerate the smallest already-local CKM-moving basis orbit.

This is diagnostic-only. It exposes the finite orbit obtained by replacing the
left singular basis with right or conjugate-right singular bases on the current
forward Yukawa surface. That orbit is useful as an exclusion result because it
shows that the missing physical object is not hidden in local linear-algebra
gauge freedom. It is not the theorem-grade `sigma_ud_orbit`.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.mixing import ckm_parameters as _ckm_tuple
FORWARD_YUKAWAS = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_local_basis_orbit_diagnostic.json"
TARGET_THETA_12 = 0.2256
TARGET_THETA_23 = 0.0438
TARGET_THETA_13 = 0.00347


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _svd_bases(matrix: np.ndarray) -> dict[str, np.ndarray]:
    left, _, right_h = np.linalg.svd(matrix)
    right = right_h.conj().T
    return {
        "L": left,
        "R": right,
        "Rbar": right.conj(),
    }


def _loss(ckm: dict[str, float]) -> float:
    return (
        (ckm["theta_12"] - TARGET_THETA_12) ** 2
        + (ckm["theta_23"] - TARGET_THETA_23) ** 2
        + (ckm["theta_13"] - TARGET_THETA_13) ** 2
    )


def build_artifact(forward: dict[str, Any]) -> dict[str, Any]:
    y_u = _complex_matrix(dict(forward["Y_u"]))
    y_d = _complex_matrix(dict(forward["Y_d"]))
    bases_u = _svd_bases(y_u)
    bases_d = _svd_bases(y_d)

    elements: list[dict[str, Any]] = []
    for basis_u_name, basis_u in bases_u.items():
        for basis_d_name, basis_d in bases_d.items():
            v_ckm = basis_u.conj().T @ basis_d
            ckm = _ckm_tuple(v_ckm)
            physical = basis_u_name == "L" and basis_d_name == "L"
            elements.append(
                {
                    "basis_u": basis_u_name,
                    "basis_d": basis_d_name,
                    "ckm": ckm,
                    "diagnostic_compare_shell_loss": _loss(ckm),
                    "physical_admissible": physical,
                    "exclusion_reason": None if physical else "CKM on the current corpus is defined on ordered same-label left eigenframes only",
                }
            )

    nonphysical = [item for item in elements if not item["physical_admissible"]]
    nonphysical.sort(key=lambda item: (item["diagnostic_compare_shell_loss"], item["basis_u"], item["basis_d"]))
    best_nonphysical = nonphysical[0] if nonphysical else None

    return {
        "artifact": "oph_quark_local_basis_orbit_diagnostic",
        "generated_utc": _timestamp(),
        "scope": "D12_current_reference_local_basis_only",
        "basis_set": ["L", "R", "Rbar"],
        "elements": elements,
        "best_nonphysical_candidate": best_nonphysical,
        "theorem_use": "diagnostic_exclusion_only",
        "notes": [
            "This is not the true sigma_ud orbit. It is the smallest finite local basis orbit already extractable from the current forward Yukawa surface.",
            "Only the L/L element remains physically admissible on the current theorem surface because CKM is defined on ordered same-label left eigenframes.",
            "A better nonphysical shell point here does not promote anything; it only proves that the remaining orbit must be left-handed and sheet-level.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Enumerate the local quark basis-orbit diagnostic.")
    parser.add_argument("--forward-yukawas", default=str(FORWARD_YUKAWAS))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    forward = _load_json(Path(args.forward_yukawas))
    payload = build_artifact(forward)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


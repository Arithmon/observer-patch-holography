#!/usr/bin/env python3
"""Finite-grid diagnostic of kernel admissibility and non-uniqueness.

For positive r,s, exact algebra gives positive gaps g21=s*r/(1+r),
g32=s/(1+r), ratio r and span s. Subtracting the mean and shifting the
spectrum positive realizes those gaps in a positive Hermitian matrix;
unitary conjugation preserves them. The corresponding identities are
rho_ord=3/(2+r) and x2=(r-1)/(r+1).

The numerical battery is a different statement. It imposes fixed positive
gap and squared-overlap floors. Two gaps strictly above GAP_FLOOR have
sum strictly above 2*GAP_FLOOR, so this battery cannot accept every positive
span. Shrinking frame drift can improve the defect inequality while making
cross-label overlaps fail their floor; nonzero overlap is not sufficient.
The bounded shrink loop below supplies no universal termination or joint
feasibility theorem. Its defect compares transport matrices with a descendant
gap as a declared numerical diagnostic, not a proved Riesz-projector bound.

This script evaluates twelve binary64 witnesses on a declared finite grid.
Passing distinct grid points is evidence of non-uniqueness within this
diagnostic, not surjectivity, a complete admissibility theorem or a physical
no-go. The historical grid contains rounded template-context coordinates;
it does not contain the exact on-disk operating point and is not an
independently source-selected prediction. No quark mass or fitted spread
is used. The stored kernel is read as declared context.

Run:
    python3 code/particles/flavor/derive_family_transport_kernel_admissibility_freedom_theorem.py
writes code/particles/runs/flavor/family_transport_kernel_admissibility_freedom_theorem.json.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone
from fractions import Fraction

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / "particles" / "runs" / "flavor"
DEFAULT_OUT = RUNS / "family_transport_kernel_admissibility_freedom_theorem.json"
KERNEL_PATH = RUNS / "family_transport_kernel.json"

GAP_FLOOR = 1.0e-12
AMPLITUDE_FLOOR = 1.0e-12  # Applied to squared eigenline overlaps.
MIXER = np.asarray([[0.0, 1.0, 1.0j],
                    [1.0, 0.0, 1.0],
                    [-1.0j, 1.0, 0.0]], dtype=complex)
EPS_FRAME_0 = 0.05
EPS_FRAME_DRIFT = 0.02
EIGENVALUE_DRIFT_REL = 0.003

R_GRID = (0.05, 0.317889, 1.0, 5.0)      # Includes a rounded context coordinate.
S_GRID = (0.1, 1.253553, 10.0)           # Includes a rounded context coordinate.


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _expm_hermitian_i(h: np.ndarray, eps: float) -> np.ndarray:
    evals, evecs = np.linalg.eigh(h)
    return evecs @ np.diag(np.exp(1j * eps * evals)) @ evecs.conj().T


def _centered_spectrum(r: float, s: float) -> np.ndarray:
    g21 = s * r / (1.0 + r)
    g32 = s / (1.0 + r)
    lam1 = -(2.0 * g21 + g32) / 3.0
    return np.asarray([lam1, lam1 + g21, lam1 + g21 + g32])


def _eigenlines(hermitian: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    evals, evecs = np.linalg.eigh(hermitian)
    return evals, evecs


def build_witness(r: float, s: float) -> dict:
    if not (math.isfinite(r) and math.isfinite(s) and r > 0 and s > 0):
        raise ValueError("positive finite ratio and span required")
    lam = _centered_spectrum(r, s)
    shift = 1.0 + abs(float(lam[0]))
    mu = lam + shift
    gap_min = float(min(lam[1] - lam[0], lam[2] - lam[1]))

    q0 = _expm_hermitian_i(MIXER, EPS_FRAME_0)
    t0 = q0 @ np.diag(np.sqrt(mu)) @ q0.conj().T

    shrink, attempts = 1.0, 0
    while True:
        attempts += 1
        eps1 = EPS_FRAME_0 + EPS_FRAME_DRIFT * shrink
        drift = EIGENVALUE_DRIFT_REL * shrink
        q1 = _expm_hermitian_i(MIXER, eps1)
        mu1 = mu * (1.0 + drift)
        t1 = q1 @ np.diag(np.sqrt(mu1)) @ q1.conj().T
        defect = float(np.linalg.norm(t1 - t0, ord=2))
        if defect < 0.5 * gap_min or attempts > 60:
            break
        shrink *= 0.5

    certificates = {}
    descendants = []
    projector_frames = []
    for level, t in enumerate((t0, t1)):
        desc = t @ t.conj().T
        evals, evecs = _eigenlines(desc)
        centered = desc - (np.trace(desc) / 3.0) * np.eye(3, dtype=complex)
        cevals = np.linalg.eigvalsh(centered)
        g21 = float(cevals[1] - cevals[0])
        g32 = float(cevals[2] - cevals[1])
        descendants.append({
            "level": level,
            "psd": bool(np.min(np.linalg.eigvalsh(desc)) > 0.0),
            "centered_eigenvalues": [float(x) for x in cevals],
            "g21": g21,
            "g32": g32,
        })
        projector_frames.append(evecs)
    certificates["psd_descendants"] = all(d["psd"] for d in descendants)
    certificates["three_cluster_gap_open_all_levels"] = all(
        min(d["g21"], d["g32"]) > 0.0 for d in descendants)
    certificates["simple_centered_spectrum"] = all(
        min(d["g21"], d["g32"]) > GAP_FLOOR for d in descendants)
    certificates["riesz_margin_defect_below_half_gap"] = bool(
        defect < 0.5 * gap_min)

    left, right = projector_frames
    overlaps = np.abs(left.conj().T @ right) ** 2
    same_label = [float(overlaps[i, i]) for i in range(3)]
    cross_label = [float(overlaps[i, j]) for i in range(3) for j in range(3)
                   if i != j]
    certificates["projector_labeling_persistent"] = bool(
        min(same_label) > max(cross_label))
    certificates["edge_amplitudes_above_floor"] = bool(
        min(cross_label) > AMPLITUDE_FLOOR)

    base = descendants[0]
    r_emitted = base["g21"] / base["g32"]
    s_emitted = base["g21"] + base["g32"]
    rho_emitted = 3.0 * base["g32"] / (2.0 * base["g32"] + base["g21"])
    x2_emitted = (r_emitted - 1.0) / (r_emitted + 1.0)

    return {
        "target": {"r": r, "s": s},
        "emitted": {
            "r": r_emitted,
            "s": s_emitted,
            "rho_ord": rho_emitted,
            "x2": x2_emitted,
        },
        "exactness": {
            "r_abs_error": abs(r_emitted - r),
            "s_abs_error": abs(s_emitted - s),
        },
        "riesz": {"defect_sup": defect, "half_gap": 0.5 * gap_min,
                  "shrink_attempts": attempts},
        "certificates": certificates,
        "all_certificates_pass": all(bool(v) for v in certificates.values()),
    }


def build() -> dict:
    kernel = json.loads(KERNEL_PATH.read_text(encoding="utf-8"))

    witnesses = []
    for r in R_GRID:
        for s in S_GRID:
            witnesses.append(build_witness(r, s))
    all_pass = all(w["all_certificates_pass"] for w in witnesses)
    max_r_err = max(w["exactness"]["r_abs_error"] for w in witnesses)
    max_s_err = max(w["exactness"]["s_abs_error"] for w in witnesses)
    if not all_pass:
        raise AssertionError("a freedom witness failed the certificate battery")
    if max(max_r_err, max_s_err) > 1.0e-9:
        raise AssertionError("a freedom witness missed its target invariants")

    span_r = (min(R_GRID), max(R_GRID))
    span_s = (min(S_GRID), max(S_GRID))
    context_level = max(kernel["refinements"], key=lambda row: row["level"])
    context_eigenvalues = sorted(float(x) for x in context_level["eigenvalues"])
    context_r = ((context_eigenvalues[1] - context_eigenvalues[0])
                 / (context_eigenvalues[2] - context_eigenvalues[1]))
    context_s = context_eigenvalues[2] - context_eigenvalues[0]

    return {
        "artifact": "oph_family_transport_kernel_admissibility_freedom_theorem",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "proof_status": "finite_grid_numerical_witnesses_only",
        "claim_tier": "declared_tolerance_battery_nonuniqueness_examples",
        "row_class": "numerical_diagnostic",
        "guards": {
            "quark_reference_values_consumed": False,
            "fitted_spreads_consumed": False,
            "kernel_template_consumed_for_context": True,
            "grid_is_source_selected": False,
            "universal_tolerance_battery_surjectivity": False,
            "physical_no_go": False,
            "public_promotion_allowed": False,
        },
        "statement": (
            "Twelve distinct declared (r,s) targets pass the binary64 "
            "diagnostic battery with target errors below 1e-9. These "
            "finite witnesses show non-uniqueness within that battery. "
            "They do not establish universal feasibility, exact numerical "
            "certification, source selection or a physical prediction."
        ),
        "proof_kind": "finite_binary64_witness_evaluation",
        "numeric_battery": {
            "strict_gap_floor": GAP_FLOOR,
            "strict_squared_overlap_floor": AMPLITUDE_FLOOR,
            "maximum_shrink_attempts": 61,
            "joint_feasibility_guaranteed": False,
            "riesz_margin_scope": "declared transport-defect/descendant-gap "
                                  "comparison; no projector perturbation theorem",
            "exact_span_obstruction": {
                "gap_floor": str(Fraction.from_float(GAP_FLOOR)),
                "excluded_span_upper_inclusive": str(
                    2 * Fraction.from_float(GAP_FLOOR)),
                "floor_semantics": "exact rational value of the binary64 "
                                   "threshold used by the numerical battery",
                "reason": "g21>f and g32>f imply s=g21+g32>2f",
                "scope": "fixed numerical acceptance floor, not a physical no-go",
            },
        },
        "construction": {
            "spectrum_placement": "trace-free target spectrum with gaps "
                                  "(s r/(1+r), s/(1+r)), shifted positive; "
                                  "T = Q diag(sqrt(mu)) Q_dagger so the "
                                  "descendant spectrum is exact under "
                                  "unitary conjugation",
            "frames": "Q_level = exp(i eps_level H) with a declared "
                      "Hermitian mixer; overlap floors and label dominance "
                      "are checked numerically at each grid target",
            "bounded_shrink_search": "frame and eigenvalue drifts shrink "
                                     "until the defect test passes or 61 "
                                     "attempts are used; other battery "
                                     "conditions are checked afterward",
        },
        "witness_grid": {
            "r_values": list(R_GRID),
            "s_values": list(S_GRID),
            "includes_on_disk_operating_point": bool(
                context_r in R_GRID and context_s in S_GRID),
            "selection": "declared grid with rounded kernel-template "
                         "context coordinates; not independent physics evidence",
            "r_span_orders_of_magnitude": math.log10(span_r[1] / span_r[0]),
            "s_span_orders_of_magnitude": math.log10(span_s[1] / span_s[0]),
            "all_certificates_pass": all_pass,
            "max_target_error": max(max_r_err, max_s_err),
        },
        "witnesses": witnesses,
        "identity_x2_of_r": {
            "statement": "x2 = (r - 1)/(r + 1) holds identically, so the "
                         "mean-law coordinate carries no information "
                         "beyond the gap ratio",
            "max_deviation_on_witnesses": max(
                abs(w["emitted"]["x2"]
                    - (w["emitted"]["r"] - 1.0) / (w["emitted"]["r"] + 1.0))
                for w in witnesses),
        },
        "on_disk_kernel_context": {
            "artifact": kernel.get("artifact"),
            "status": kernel.get("status"),
            "raw_gap_ratio_r": context_r,
            "spectral_span": context_s,
            "nearest_r_grid_distance": min(abs(r-context_r) for r in R_GRID),
            "nearest_s_grid_distance": min(abs(s-context_s) for s in S_GRID),
            "note": "declared template context, not an independently "
                    "source-selected or exactly sampled physical target",
        },
        "corollary_for_issue_377": (
            "The finite passing witnesses do not select a unique ratio "
            "or span. A source selection principle requires separate "
            "evidence; this diagnostic does not classify all admissible "
            "kernels or all physical quark interfaces."
        ),
        "relation_to_nonidentifiability": (
            "Exact positive-parameter spectrum placement is separate from "
            "the finite tolerance battery. The latter excludes sufficiently "
            "small spans and has multiple passing examples, without a "
            "universal nonidentifiability or physical no-go conclusion."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate a finite grid of kernel admissibility witnesses.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    report = build()
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    grid = report["witness_grid"]
    print(f"witness grid: r in {grid['r_values']}  s in {grid['s_values']}")
    print(f"all certificates pass: {grid['all_certificates_pass']}")
    print(f"max target error: {grid['max_target_error']:.2e}")
    print(f"x2 = (r-1)/(r+1) max deviation: "
          f"{report['identity_x2_of_r']['max_deviation_on_witnesses']:.2e}")
    for w in report["witnesses"]:
        t, e = w["target"], w["emitted"]
        print(f"  r={t['r']:<9g} s={t['s']:<9g} -> rho={e['rho_ord']:.6f} "
              f"x2={e['x2']:+.6f}  pass={w['all_certificates_pass']}")
    print(f"proof status: {report['proof_status']}")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

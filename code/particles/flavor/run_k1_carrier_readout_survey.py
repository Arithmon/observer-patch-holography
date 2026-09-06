#!/usr/bin/env python3
"""K1 carrier-readout survey: the finite grammar, swept with accounting (#377).

Purpose. Replace one-at-a-time carrier iteration with a single sweep of
every definable readout cell over the realized data on disk, under one
multiplicity-controlled protocol. The survey does not promote anything:
it publishes the full outcome table, the trials count, and a shortlist
under predeclared criteria; shortlist survivors earn only the right to a
preregistered replication on fresh run seeds.

Protocol (fixed before execution):

  grammar      four readout families over the realized data:
               A. response-channel second-moment Grams (modular response
                  cache, s3 class-delta channels) per run;
               B. class-fusion one-step composition law of the gauge
                  state (descendant M M^t) per run, with the analytic
                  Haar law as the family null;
               C. tower-compression sector kernels (majority-class cells,
                  conditional-expectation block scales 64 and 16) per
                  run with gauge plus freezeout exports;
               D. vacuum-inclusive population-transfer VAR kernels on the
                  dense timeline: series in {log, raw}, burn-in in
                  {0, quarter}, lag in {1, 2, 4}, with first-half versus
                  full-series stability.
               Basis variants are enumerated only where they change the
               spectrum; unitary basis changes of a Gram do not and are
               declared out.
  nulls        per family: label or time shuffles (two fixed seeds) and,
               for family B, the closed-form Haar law. A physical cell
               with no usable null cannot enter the shortlist.
  band         the quoted-convention band of the common-scheme scan; the
               template rho_ord is display context only.
  multiplicity the number of physical cells N is printed and stored; the
               expected chance count of in-band cells under a uniform
               prior on (0, 1.5) is N times the band fraction, and any
               shortlist must be read against it.
  shortlist    a physical cell qualifies only if (i) rho_ord is in band,
               (ii) its distance to every null rho_ord of its cell is at
               least 0.05, and (iii) where a level or half split exists,
               the two estimates differ by at most 0.05.
  replication  shortlist survivors are candidates only; the declared next
               step for any survivor is an identical rerun on fresh run
               seeds before harness gating and any further status.

No quark mass, fitted spread, or flavor template enters any step. Reruns
of already-executed attempts appear as table rows with their attempt
numbers for continuity with the issue ledger.

Run:
    <sim venv python> code/particles/flavor/run_k1_carrier_readout_survey.py
writes code/particles/runs/flavor/k1_carrier_readout_survey.json.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
RER_ROOT = HERE.parents[2]
SIM_ROOT = HERE.parents[3] / "oph-physics-sim"
RUNS = RER_ROOT / "code" / "particles" / "runs" / "flavor"
OUT_PATH = RUNS / "k1_carrier_readout_survey.json"
SCAN_PATH = RUNS / "quark_common_scheme_shape_law_scan.json"

sys.path.insert(0, str(HERE))
from derive_family_transport_kernel_from_sim_tower_compression import (  # noqa: E402
    load_state as tower_load_state,
    sector_kernel as tower_sector_kernel,
)
from derive_family_transport_kernel_from_sim_population_transfer import (  # noqa: E402
    load_population_series,
)

NULL_SEEDS = (11, 12)
NULL_DISTANCE_MIN = 0.05
STABILITY_MAX = 0.05

RESPONSE_RUNS = (
    "oph_universe_64k_final_audited_20260711",
    "oph_universe_4k_viewer_smoke_20260706_v2",
)
GAUGE_RUNS = (
    "e1_s3_bw_screen_64k_1783837538",
    "e1_s3_bw_screen_4k_1783837537",
    "k1_fusion_universe_4k_20260712",
    "k1_fusion_universe_64k_20260712",
)
TOWER_RUNS = (
    "k1_fusion_universe_4k_20260712",
    "k1_fusion_universe_64k_20260712",
)
TIMELINE_RUNS = (
    "k1_population_transfer_4k_dense_20260712",
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _sim_module(rel: str, name: str):
    spec = importlib.util.spec_from_file_location(name, SIM_ROOT / rel)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def readout(d: np.ndarray) -> dict:
    centered = d - (np.trace(d) / 3.0) * np.eye(3)
    evals = np.linalg.eigvalsh(centered)
    g21 = float(evals[1] - evals[0])
    g32 = float(evals[2] - evals[1])
    r = g21 / g32 if g32 > 1.0e-300 else float("inf")
    rho = 3.0 * g32 / (2.0 * g32 + g21) if (2.0 * g32 + g21) > 0 else float("nan")
    return {"rho_ord": rho, "raw_gap_ratio_r": r,
            "eigenvalues_centered": [float(x) for x in evals]}


# ---------------- family A: response-channel Grams ----------------

def family_a_cells() -> list[dict]:
    cells = []
    for run_id in RESPONSE_RUNS:
        run_dir = SIM_ROOT / "runs" / run_id
        cache_path = run_dir / "modular_response_kernel_cache.json"
        npz_path = run_dir / "modular_response_kernel_payload.npz"
        if not cache_path.exists() or not npz_path.exists():
            continue
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
        payload = np.load(npz_path)
        columns: dict[tuple[int, int], dict[int, int]] = {}
        fields = tuple(f"s3_sector_class:class_delta_{c}" for c in range(3))
        for row in cache["feature_rows"]:
            if row["field"] in fields:
                key = (int(row["cap_index"]), int(row["time_index"]))
                columns.setdefault(key, {})[fields.index(row["field"])] = int(
                    row["feature_index"])

        def slab(matrix: np.ndarray) -> np.ndarray:
            stacked = [matrix[:, [columns[k][c] for c in range(3)]]
                       for k in sorted(columns)]
            t = np.vstack(stacked)
            return t.T @ t

        matrix = np.asarray(payload["matrix"], dtype=float)
        phys = readout(slab(matrix))
        nulls = {}
        for seed in NULL_SEEDS:
            rng = np.random.default_rng(seed)
            shuffled = matrix[:, rng.permutation(matrix.shape[1])]
            nulls[f"column_shuffle_seed{seed}"] = readout(slab(shuffled))
        cells.append({
            "cell_id": f"A.response_gram.{run_id}",
            "family": "A_response_channel_gram",
            "attempt_ref": 1 if "audited" in run_id else None,
            "physical": phys, "nulls": nulls,
            "stability_pair": None,
        })
    return cells


# ---------------- family B: class-fusion law ----------------

def fusion_counts(left: np.ndarray, right: np.ndarray, labels: np.ndarray,
                  s3_class: np.ndarray, s3_mul: np.ndarray,
                  s3_inv: np.ndarray) -> np.ndarray:
    """Count classes on directed non-backtracking two-edge paths.

    A stored edge (i, j, g) represents i->j with g and j->i with g^-1.
    Reversing storage endpoints and inverting g leaves the path counts fixed.
    """
    head: dict[int, list[tuple[int, int]]] = {}
    for idx in range(left.shape[0]):
        i, j, g = int(left[idx]), int(right[idx]), int(labels[idx])
        head.setdefault(i, []).append((j, g))
        head.setdefault(j, []).append((i, int(s3_inv[g])))
    counts = np.zeros((3, 3), dtype=float)
    for idx in range(left.shape[0]):
        i, j, g1 = int(left[idx]), int(right[idx]), int(labels[idx])
        for k, g2 in head.get(j, ()):
            if k == i:
                continue
            counts[s3_class[g1], s3_class[s3_mul[g1, g2]]] += 1.0
        g1r = int(s3_inv[g1])
        for k, g2 in head.get(i, ()):
            if k == j:
                continue
            counts[s3_class[g1r], s3_class[s3_mul[g1r, g2]]] += 1.0
    return counts


def family_b_cells() -> list[dict]:
    holonomy = _sim_module("oph_fpe/defects/array_s3_holonomy.py", "s3holo")
    s3_class = np.asarray(holonomy.S3_CLASS, dtype=np.int64)
    s3_mul = np.asarray(holonomy.S3_MUL, dtype=np.int64)
    s3_inv = np.asarray(holonomy.S3_INV, dtype=np.int64)
    cells = []
    for run_id in GAUGE_RUNS:
        path = SIM_ROOT / "runs" / run_id / "s3_gauge_state.npz"
        if not path.exists():
            continue
        state = np.load(path)
        left = np.asarray(state["left"], dtype=np.int64)
        right = np.asarray(state["right"], dtype=np.int64)
        gauge = np.asarray(state["gauge"], dtype=np.int64)

        def fusion_matrix(labels: np.ndarray) -> np.ndarray:
            counts = fusion_counts(left, right, labels, s3_class, s3_mul, s3_inv)
            rows = counts.sum(axis=1, keepdims=True)
            m = counts / np.where(rows > 0, rows, 1.0)
            return m @ m.T

        phys = readout(fusion_matrix(gauge))
        nulls = {}
        for seed in NULL_SEEDS:
            rng = np.random.default_rng(seed)
            nulls[f"label_shuffle_seed{seed}"] = readout(
                fusion_matrix(rng.permutation(gauge)))
        haar_rows = np.tile(np.asarray([1 / 6, 1 / 2, 1 / 3]), (3, 1))
        nulls["haar_analytic"] = readout(haar_rows @ haar_rows.T)
        cells.append({
            "cell_id": f"B.class_fusion.{run_id}",
            "family": "B_class_fusion_law",
            "attempt_ref": 2 if run_id.startswith("e1_") else None,
            "physical": phys, "nulls": nulls,
            "stability_pair": None,
        })
    return cells


# ---------------- family C: tower compression ----------------

def family_c_cells() -> list[dict]:
    cells = []
    for run_id in TOWER_RUNS:
        run_dir = SIM_ROOT / "runs" / run_id
        if not ((run_dir / "s3_gauge_state.npz").exists()
                and (run_dir / "freezeout_fields.npz").exists()):
            continue
        state = tower_load_state(run_dir)
        for ratio in (64, 16):
            phys_block = tower_sector_kernel(state, ratio, state["mu"])
            phys = readout(np.asarray(phys_block["kernel"]))
            nulls = {}
            for seed in NULL_SEEDS:
                rng = np.random.default_rng(seed)
                shuffled = dict(state)
                shuffled["cell_class"] = rng.permutation(state["cell_class"])
                block = tower_sector_kernel(shuffled, ratio, state["mu"])
                nulls[f"label_shuffle_seed{seed}"] = readout(
                    np.asarray(block["kernel"]))
            cells.append({
                "cell_id": f"C.tower_ratio{ratio}.{run_id}",
                "family": "C_tower_compression",
                "attempt_ref": 4 if ratio == 16 else None,
                "physical": phys, "nulls": nulls,
                "stability_pair": None,
            })
    return cells


# ---------------- family D: population transfer ----------------

def var_transfer(series: np.ndarray, lag: int) -> np.ndarray | None:
    centered = series - series.mean(axis=0, keepdims=True)
    x0 = centered[:-lag]
    x1 = centered[lag:]
    c0 = x0.T @ x0
    if np.linalg.cond(c0) > 1.0e12:
        return None
    c1 = x1.T @ x0
    a = c1 @ np.linalg.inv(c0)
    return a @ a.T


def family_d_cells() -> list[dict]:
    cells = []
    for run_id in TIMELINE_RUNS:
        run_dir = SIM_ROOT / "runs" / run_id
        if not (run_dir / "defect_timeline_report.json").exists():
            continue
        series_raw = load_population_series(run_dir)["series"]
        variants = {"log": np.log(series_raw)} if (series_raw > 0).all() else {}
        variants["raw"] = series_raw
        for series_name, base in variants.items():
            for burn_label, burn in (("burn0", 0), ("burnQ", base.shape[0] // 4)):
                series = base[burn:]
                for lag in (1, 2, 4):
                    if series.shape[0] < 12 * lag:
                        continue
                    cell_id = (f"D.pop_var.{series_name}.{burn_label}"
                               f".lag{lag}.{run_id}")
                    d_full = var_transfer(series, lag)
                    if d_full is None:
                        cells.append({
                            "cell_id": cell_id,
                            "family": "D_population_transfer",
                            "attempt_ref": None,
                            "physical": {"rho_ord": float("nan"),
                                         "raw_gap_ratio_r": float("nan"),
                                         "degenerate_estimator": True},
                            "nulls": {},
                            "stability_pair": None,
                        })
                        continue
                    phys = readout(d_full)
                    d_half = var_transfer(series[: series.shape[0] // 2], lag)
                    half_rho = (readout(d_half)["rho_ord"]
                                if d_half is not None else float("nan"))
                    nulls = {}
                    for seed in NULL_SEEDS:
                        rng = np.random.default_rng(seed)
                        d_null = var_transfer(
                            series[rng.permutation(series.shape[0])], lag)
                        if d_null is not None:
                            nulls[f"time_shuffle_seed{seed}"] = readout(d_null)
                    cells.append({
                        "cell_id": cell_id,
                        "family": "D_population_transfer",
                        "attempt_ref": 5 if (series_name == "log" and burn == 0
                                             and lag == 1) else None,
                        "physical": phys, "nulls": nulls,
                        "stability_pair": {"first_half_rho": half_rho,
                                           "full_rho": phys["rho_ord"]},
                    })
    return cells


def require_declared_inputs() -> None:
    """Reject an incomplete grammar sweep before reading or replacing a report."""
    paths = {SCAN_PATH}
    for run_id in RESPONSE_RUNS:
        paths.update(SIM_ROOT / "runs" / run_id / name for name in (
            "modular_response_kernel_cache.json", "modular_response_kernel_payload.npz"))
    for run_id in GAUGE_RUNS:
        paths.add(SIM_ROOT / "runs" / run_id / "s3_gauge_state.npz")
    for run_id in TOWER_RUNS:
        paths.update(SIM_ROOT / "runs" / run_id / name for name in (
            "s3_gauge_state.npz", "freezeout_fields.npz"))
    for run_id in TIMELINE_RUNS:
        paths.add(SIM_ROOT / "runs" / run_id / "defect_timeline_report.json")
    missing = sorted(path.as_posix() for path in paths if not path.is_file())
    if missing:
        raise FileNotFoundError(
            "Declared survey inputs missing; refusing a partial survey: "
            + ", ".join(missing))


def build() -> dict:
    require_declared_inputs()
    scan = json.loads(SCAN_PATH.read_text(encoding="utf-8"))
    band = scan["findings"]["a2_band_by_convention"]["quoted_mixed"]
    template_rho = float(scan["template_rho_ord"])

    cells = (family_a_cells() + family_b_cells()
             + family_c_cells() + family_d_cells())

    shortlist = []
    for cell in cells:
        rho = cell["physical"]["rho_ord"]
        in_band = bool(band[0] <= rho <= band[1]) if np.isfinite(rho) else False
        null_rhos = [row["rho_ord"] for row in cell["nulls"].values()
                     if np.isfinite(row["rho_ord"])]
        null_dist = (min(abs(rho - nr) for nr in null_rhos)
                     if (null_rhos and np.isfinite(rho)) else None)
        stab = cell["stability_pair"]
        stable = (abs(stab["first_half_rho"] - stab["full_rho"]) <= STABILITY_MAX
                  if stab and all(np.isfinite(v) for v in stab.values())
                  else None)
        cell["evaluation"] = {
            "in_band": in_band,
            "min_null_distance": null_dist,
            "null_separated": bool(null_dist is not None
                                   and null_dist >= NULL_DISTANCE_MIN),
            "stable": stable,
        }
        qualifies = (in_band and cell["evaluation"]["null_separated"]
                     and stable is not False)
        cell["evaluation"]["shortlisted"] = bool(qualifies)
        if qualifies:
            shortlist.append(cell["cell_id"])

    n_cells = len(cells)
    band_fraction = (band[1] - band[0]) / 1.5
    return {
        "artifact": "oph_k1_carrier_readout_survey",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "row_class": "systematic_readout_survey",
        "guards": {
            "quark_reference_values_consumed": False,
            "fitted_spreads_consumed": False,
            "numerical_flavor_template_consumed": False,
            "public_promotion_allowed": False,
            "selection_policy": "full table published; shortlist criteria "
                                "predeclared; survivors require replication "
                                "on fresh run seeds before harness gating",
        },
        "band": band,
        "template_rho_ord_display_only": template_rho,
        "multiplicity_accounting": {
            "physical_cells": n_cells,
            "band_fraction_of_rho_range": band_fraction,
            "expected_chance_in_band_count": n_cells * band_fraction,
        },
        "shortlist_criteria": {
            "in_band": True,
            "min_null_distance": NULL_DISTANCE_MIN,
            "stability_max_delta_rho": STABILITY_MAX,
        },
        "cells": cells,
        "shortlist": shortlist,
        "declared_next_step_for_survivors": (
            "identical rerun on fresh run seeds (new sim runs of the same "
            "configs); only replicated survivors proceed to the acceptance "
            "harness"),
        "runs_missing_at_execution": [
            run_id for run_id in GAUGE_RUNS + TOWER_RUNS
            if not (SIM_ROOT / "runs" / run_id / "s3_gauge_state.npz").exists()],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Sweep the K1 carrier-readout grammar with nulls and "
                    "multiplicity accounting.")
    parser.add_argument("--output", default=str(OUT_PATH))
    args = parser.parse_args()

    report = build()
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    acc = report["multiplicity_accounting"]
    print(f"physical cells: {acc['physical_cells']}  "
          f"expected chance in-band: {acc['expected_chance_in_band_count']:.2f}")
    print(f"band: {report['band']}")
    for cell in report["cells"]:
        ev = cell["evaluation"]
        rho = cell["physical"]["rho_ord"]
        stab = cell["stability_pair"]
        stab_str = (f" half/full {stab['first_half_rho']:.3f}/{stab['full_rho']:.3f}"
                    if stab else "")
        print(f"  {cell['cell_id']:52s} rho={rho:8.4f} "
              f"in_band={str(ev['in_band']):5s} "
              f"null_dist={ev['min_null_distance'] if ev['min_null_distance'] is None else round(ev['min_null_distance'], 4)}"
              f"{stab_str}"
              f"{'  <-- SHORTLIST' if ev['shortlisted'] else ''}")
    print(f"shortlist: {report['shortlist'] or 'EMPTY'}")
    if report["runs_missing_at_execution"]:
        print(f"missing runs (rerun survey when they land): "
              f"{report['runs_missing_at_execution']}")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

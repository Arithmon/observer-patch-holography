"""Refine-and-settle ensemble on exact ancestral bins (issue 987), many short chains in parallel.

One chain: take a settled parent state at level m (readings q, q+1), refine it (every carrier's twelve
readings are inherited by its four children, the tower indexing is hierarchical), add fresh i.i.d. records
uniform on {0..4} (mean 2, which keeps the child's mean at a half-integer; see refine_settle_cocycle.py),
settle level m+1 with the canonical integer law, and measure how the inherited coarse field survives.

The measurement follows the audit of 2026-09-25 (plan/audits/SIM_COSMOLOGY_2026-09-25): equal-count
ancestral groups (the 4^g contiguous children of one level-(m+1-g) ancestor) instead of unequal cube
caps, one centering and weight convention, and a joint regression of the settled group means on the
inherited and the fresh group means:

    y = b x + c f + e,   <x, e> = <f, e> = 0,

with response b (the survival amplitude of the inherited pattern), fresh transfer c, the total power
ratio <y,y>/<x,x>, the innovation power <e,e>/<x,x> and the squared coherence <x,y>^2/(<x,x><y,y>).
The generator per refinement step is theta = -log2(b^2).  Controls: `shuffled` permutes the parent's
carriers before refinement (same multiset, no spatial structure; b must vanish), and `inherit_only`
(no fresh records) is the exact amplitude-one control, which needs no settlement and is recorded as such.

Chains run in parallel, one process each (the integer kernel is sequential); the ensemble spread over
parents and fresh seeds is the uncertainty.

usage: python3 refine_ensemble.py SIM_ROOT CACHE OUT_DIR --chains SPEC [SPEC ...] --processes P
  SPEC = parent_run_dir:parent_level:schedule_index:fresh_seed:mode   mode in {chain, shuffled}
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

P_STAR = 1.6309682094
GROUP_SIZES = [4 ** g for g in range(1, 9)]  # 4 .. 65536 cells


def settle(H, geo, loads: np.ndarray, seed: int, log=None) -> dict:
    lib = H.native_kernel()
    x = np.ascontiguousarray(loads.astype(np.int8).copy())
    tpl_a = np.ascontiguousarray(geo.template[:, 0], dtype=np.int8)
    tpl_b = np.ascontiguousarray(geo.template[:, 1], dtype=np.int8)
    v_min = int(H.expectation(geo, loads.astype(np.int64))["balanced_minimum"])
    v = np.array([int(np.dot(loads.astype(np.int64), loads.astype(np.int64)))], dtype=np.int64)
    counts = np.zeros(5, dtype=np.int64)
    rng = np.random.default_rng(seed)
    first = 0 if int(v[0]) == v_min else -1
    sweep = 0
    t0 = time.time()
    while first < 0 and sweep < 20000:
        seq, coin, _ = H.draw_sweep(rng, geo.seams, coins=True)
        first = int(lib.oph_huge_integer(H._ptr(x, ctypes.c_int8), H._ptr(tpl_a, ctypes.c_int8), H._ptr(tpl_b, ctypes.c_int8),
                                         H._ptr(geo.ia, ctypes.c_int32), H._ptr(geo.ib, ctypes.c_int32), geo.intra_count,
                                         H._ptr(seq, ctypes.c_int32), H._ptr(coin, ctypes.c_int8), geo.seams,
                                         H._ptr(counts, ctypes.c_int64), H._ptr(v, ctypes.c_int64), v_min))
        sweep += 1
        if log and sweep % 16 == 0:
            log(f"sweep {sweep} V-Vmin {int(v[0]) - v_min} {time.time() - t0:.0f}s")
    return {"state": x, "sweeps": sweep, "terminated": first >= 0, "seconds": round(time.time() - t0, 1),
            "counts": counts.tolist()}


def group_means(cell_field: np.ndarray, size: int) -> np.ndarray:
    n = cell_field.size
    return cell_field[: (n // size) * size].reshape(-1, size).mean(axis=1)


def regress(x: np.ndarray, f: np.ndarray | None, y: np.ndarray) -> dict:
    xc, yc = x - x.mean(), y - y.mean()
    xx, xy, yy = float(xc @ xc), float(xc @ yc), float(yc @ yc)
    row = {"groups": int(x.size), "response_b": xy / xx if xx > 0 else None,
           "total_power_ratio": yy / xx if xx > 0 else None,
           "squared_coherence": (xy * xy) / (xx * yy) if xx > 0 and yy > 0 else None}
    if row["response_b"] is not None:
        e = yc - row["response_b"] * xc
        row["innovation_power_ratio"] = float(e @ e) / xx
        b = row["response_b"]
        row["theta_per_step"] = float(-np.log2(b * b)) if b > 0 else None
    if f is not None and xx > 0:
        fc = f - f.mean()
        A = np.array([[xc @ xc, xc @ fc], [fc @ xc, fc @ fc]])
        rhs = np.array([xc @ yc, fc @ yc])
        try:
            bj, cj = np.linalg.solve(A, rhs)
            e2 = yc - bj * xc - cj * fc
            row.update({"joint_response_b": float(bj), "joint_fresh_c": float(cj), "joint_innovation_power_ratio": float(e2 @ e2) / xx,
                        "joint_theta_per_step": float(-np.log2(bj * bj)) if bj > 0 else None,
                        "fresh_power_over_inherited": float(fc @ fc) / xx})
        except np.linalg.LinAlgError:
            pass
    return row


def run_chain(task: dict) -> dict:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.path.insert(0, task["sim_root"])
    from oph_exact import federation_huge as H
    t0 = time.time()
    cache = Path(task["cache"])
    parent_level, child_level = task["parent_level"], task["parent_level"] + 1
    receipt = json.loads((Path(task["parent_run"]) / "receipt.json").read_text())
    e = receipt["integer_law"]["entries"][task["schedule_index"]]
    state = np.load(Path(task["parent_run"]) / f"integer_{e['seed']}" / e["terminal_state"]["path"]).astype(np.int64)
    assert hashlib.sha256(state.astype("<i1").tobytes()).hexdigest() == e["terminal_state"]["sha256"], "parent digest"
    bit = state - int(state.min())
    assert bit.max() <= 1
    blocks = bit.reshape(-1, 12)
    rng = np.random.default_rng(task["fresh_seed"])
    if task["mode"] == "shuffled":
        blocks = blocks[rng.permutation(blocks.shape[0])]
    inherited = np.repeat(blocks, 4, axis=0).reshape(-1)
    fresh = rng.integers(0, 5, size=inherited.size)
    loads = (inherited + fresh).astype(np.int64)
    geo = H.Geometry(cache, child_level)
    n = geo.carriers
    seed = H.schedule_seed(child_level, 50 + task["chain_index"])
    res = settle(H, geo, loads, seed)
    settled = res["state"].astype(np.int64)
    inh_cell = inherited.reshape(n, 12).mean(axis=1)
    fresh_cell = fresh.reshape(n, 12).mean(axis=1)
    set_cell = settled.reshape(n, 12).mean(axis=1)
    rows = []
    for size in GROUP_SIZES:
        if n // size < 8:
            break
        x, f, y = group_means(inh_cell, size), group_means(fresh_cell, size), group_means(set_cell, size)
        r = regress(x, f, y)
        r["cells_per_group"] = size
        rows.append(r)
    out = {"parent_run": task["parent_run"], "parent_level": parent_level, "child_level": child_level, "schedule_index": task["schedule_index"],
           "parent_seed": int(e["seed"]), "fresh_seed": task["fresh_seed"], "mode": task["mode"], "child_schedule_seed": seed,
           "carriers": n, "sweeps": res["sweeps"], "terminated": res["terminated"], "settle_seconds": res["seconds"], "move_counts": res["counts"],
           "child_terminal_sha256": hashlib.sha256(res["state"].astype("<i1").tobytes()).hexdigest(),
           "ancestral_groups": rows, "wall_seconds": round(time.time() - t0, 1)}
    Path(task["out"]).write_text(json.dumps(out, indent=1) + "\n")
    return out


def aggregate(results: list) -> dict:
    agg = {}
    for mode in sorted({r["mode"] for r in results}):
        for level in sorted({r["parent_level"] for r in results if r["mode"] == mode}):
            sel = [r for r in results if r["mode"] == mode and r["parent_level"] == level]
            rows = []
            for size in GROUP_SIZES:
                vals = {k: [] for k in ("response_b", "joint_response_b", "joint_fresh_c", "total_power_ratio", "innovation_power_ratio",
                                        "squared_coherence", "theta_per_step", "joint_theta_per_step")}
                for r in sel:
                    g = next((row for row in r["ancestral_groups"] if row["cells_per_group"] == size), None)
                    if g is None:
                        continue
                    for k in vals:
                        if g.get(k) is not None:
                            vals[k].append(g[k])
                if not vals["response_b"]:
                    continue
                row = {"cells_per_group": size, "chains": len(vals["response_b"])}
                for k, v in vals.items():
                    if v:
                        a = np.array(v)
                        row[k] = {"mean": float(a.mean()), "sd": float(a.std(ddof=1)) if len(a) > 1 else None,
                                  "se": float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else None}
                rows.append(row)
            agg[f"{mode}_L{level}_to_L{level + 1}"] = {"chains": len(sel), "sweeps": {"mean": float(np.mean([r["sweeps"] for r in sel])),
                                                                                         "min": int(min(r["sweeps"] for r in sel)), "max": int(max(r["sweeps"] for r in sel))},
                                                        "by_group_size": rows}
    return agg


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("sim_root"); ap.add_argument("cache"); ap.add_argument("out_dir")
    ap.add_argument("--chains", nargs="+", required=True); ap.add_argument("--processes", type=int, default=8)
    a = ap.parse_args()
    out_dir = Path(a.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for i, spec in enumerate(a.chains):
        run, lvl, idx, fseed, mode = spec.split(":")
        tasks.append({"sim_root": a.sim_root, "cache": a.cache, "parent_run": run, "parent_level": int(lvl), "schedule_index": int(idx),
                      "fresh_seed": int(fseed), "mode": mode, "chain_index": i, "out": str(out_dir / f"chain_{i:03d}_L{lvl}_s{idx}_f{fseed}_{mode}.json")})
    t0 = time.time()
    with mp.get_context("spawn").Pool(a.processes) as pool:
        results = []
        for r in pool.imap_unordered(run_chain, tasks):
            results.append(r)
            print(f"chain L{r['parent_level']}->L{r['child_level']} s{r['schedule_index']} f{r['fresh_seed']} {r['mode']}: sweeps {r['sweeps']}, "
                  f"b by group {[round(g['response_b'], 4) for g in r['ancestral_groups']]}, {time.time() - t0:.0f}s", flush=True)
    summary = {"schema": "oph.exploratory.refine-ensemble.v1", "targets": {"P_over_24": P_STAR / 24, "P_over_48": P_STAR / 48},
               "definition": __doc__.split("usage:")[0].strip(), "group_sizes": GROUP_SIZES, "chains": len(results),
               "aggregate": aggregate(results), "wall_seconds": round(time.time() - t0, 1)}
    (out_dir / "ensemble_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for key, block in summary["aggregate"].items():
        print(key, "chains", block["chains"], "sweeps", block["sweeps"])
        for row in block["by_group_size"]:
            print(f"   {row['cells_per_group']:6d} cells: b {row['response_b']['mean']:.4f} +- {row['response_b']['se'] or 0:.4f}, "
                  f"joint b {row.get('joint_response_b', {}).get('mean', float('nan')):.4f}, c {row.get('joint_fresh_c', {}).get('mean', float('nan')):.4f}, "
                  f"coherence {row['squared_coherence']['mean']:.4f}, theta {row['theta_per_step']['mean']:.4f} +- {row['theta_per_step']['se'] or 0:.4f}")


if __name__ == "__main__":
    main()

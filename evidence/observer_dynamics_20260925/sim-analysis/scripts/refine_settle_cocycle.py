"""Refine-and-settle cocycle: how much of a settled record's coarse covariance survives one refinement step.

The reserve-generator receipt asks for a covariance-survival cocycle under logarithmic refinement.
This experiment builds the most literal such object the tower offers: take a settled state at
level m, refine it (every carrier's twelve readings are inherited by its four children, the tower
indexing is hierarchical), add the fresh records of the finer level (i.i.d. loads on every port,
the same law as the initial loads), settle at level m+1 with the canonical integer law, and
measure at every scale the amplitude `b` with which the inherited pattern survives in the settled
field (regression of the settled cap excess on the inherited cap excess) and its covariance
survival `b^2`.  Repeating the step gives the cocycle test: the survival over two steps against
the product of the single-step survivals.  The generator per refinement step (`b = 2`) is
`theta = -log2(b^2)`; the certificate's targets are `P*/24 = 0.0680` (full collar) and
`P*/48 = 0.0340` (source-facing half).

usage: python3 refine_settle_cocycle.py SIM_ROOT CACHE START_RUN_DIR START_LEVEL STEPS OUT.json [--fresh uniform5|uniform6|none] [--inherit bit|full]

Fresh law.  The engine's own initial loads are uniform on {0..5} (mean 5/2).  A settled parent reads q or q+1 on every port
with the raised ports at fraction f; inheriting those readings and adding fresh uniform {0..5} loads puts the child's
mean at q + f + 3, an integer when f = 1/2, so the child's balanced minimum has (almost) no raised ports and every
fluctuation of the fresh loads must be transported across the whole sphere: the settlement becomes transport-limited
(level 7: 4,752 sweeps and still 80,064 above the minimum, against 76 sweeps for white loads).  The default fresh law
is therefore uniform on {0..4} (mean 2), which keeps the child's mean at a half-integer; `--inherit bit` inherits the
parent's raised-port indicator (its readings minus q), which is the same dynamics as inheriting the full reading
(the law depends on differences only) and keeps the loads inside the engine's declared range {0..5}.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

P_STAR = 1.6309682094


def caps(points: np.ndarray, res: int) -> tuple[np.ndarray, np.ndarray]:
    d = points / np.linalg.norm(points, axis=1, keepdims=True)
    key = np.minimum(np.floor((d + 1.0) * res / 2.0).astype(np.int64), res - 1)
    lab = key[:, 0] * res * res + key[:, 1] * res + key[:, 2]
    ids, inv = np.unique(lab, return_inverse=True)
    return inv, np.bincount(inv)


def cap_excess(field: np.ndarray, inv: np.ndarray, counts: np.ndarray) -> np.ndarray:
    keep = counts >= 4
    m = (np.bincount(inv, weights=field) / np.maximum(counts, 1))[keep]
    return m - float(field.mean())


def survival_by_scale(points: np.ndarray, inherited: np.ndarray, settled: np.ndarray, fresh: np.ndarray | None, resolutions=(3, 6, 12, 24, 48, 96)) -> list[dict]:
    rows = []
    for res in resolutions:
        inv, counts = caps(points, res)
        if int((counts >= 4).sum()) < 8:
            continue
        ei, es = cap_excess(inherited, inv, counts), cap_excess(settled, inv, counts)
        b = float(np.dot(ei, es) / np.dot(ei, ei)); c = float(np.corrcoef(ei, es)[0, 1])
        row = {"resolution": res, "cells_per_cap": round(float(counts[counts >= 4].mean()), 1), "amplitude_b": round(b, 6), "correlation": round(c, 6),
               "covariance_survival_b2": round(b * b, 6), "theta_per_step": round(float(-np.log2(b * b)) if b > 0 else float("nan"), 5)}
        if fresh is not None:
            ef = cap_excess(fresh, inv, counts)
            row["fresh_amplitude"] = round(float(np.dot(ef, es) / np.dot(ef, ef)), 6)
        rows.append(row)
    return rows


def settle(H, geo, loads: np.ndarray, seed: int) -> dict:
    """The engine's integer law (same kernel, draws and termination) on loads of any nonnegative range up to 127."""
    import ctypes
    lib = H.native_kernel()
    x = np.ascontiguousarray(loads.astype(np.int8).copy())
    tpl_a = np.ascontiguousarray(geo.template[:, 0], dtype=np.int8); tpl_b = np.ascontiguousarray(geo.template[:, 1], dtype=np.int8)
    v_min = int(H.expectation(geo, loads.astype(np.int64))["balanced_minimum"])
    v = np.array([int(np.dot(loads.astype(np.int64), loads.astype(np.int64)))], dtype=np.int64); counts = np.zeros(5, dtype=np.int64)
    rng = np.random.default_rng(seed); first = 0 if int(v[0]) == v_min else -1; sweep = 0; t0 = time.time()
    while first < 0 and sweep < 100000:
        seq, coin, _ = H.draw_sweep(rng, geo.seams, coins=True)
        first = int(lib.oph_huge_integer(H._ptr(x, ctypes.c_int8), H._ptr(tpl_a, ctypes.c_int8), H._ptr(tpl_b, ctypes.c_int8), H._ptr(geo.ia, ctypes.c_int32), H._ptr(geo.ib, ctypes.c_int32), geo.intra_count, H._ptr(seq, ctypes.c_int32), H._ptr(coin, ctypes.c_int8), geo.seams, H._ptr(counts, ctypes.c_int64), H._ptr(v, ctypes.c_int64), v_min))
        sweep += 1
        if sweep % 16 == 0:
            print(f"  L{geo.level} sweep {sweep} V-Vmin {int(v[0]) - v_min} {time.time() - t0:.0f}s", flush=True)
    return {"state": x, "sweeps": sweep, "terminated": first >= 0, "seconds": round(time.time() - t0, 1)}


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("sim_root"); ap.add_argument("cache"); ap.add_argument("start_run"); ap.add_argument("start_level", type=int); ap.add_argument("steps", type=int); ap.add_argument("out"); ap.add_argument("--fresh", default="uniform5", choices=["uniform5", "uniform6", "none"]); ap.add_argument("--inherit", default="bit", choices=["bit", "full"])
    a = ap.parse_args(); sys.path.insert(0, a.sim_root)
    from oph_exact import federation_huge as H
    cache = Path(a.cache); t0 = time.time()
    receipt = json.loads((Path(a.start_run) / "receipt.json").read_text()); e = receipt["integer_law"]["entries"][0]
    state = np.load(Path(a.start_run) / f"integer_{e['seed']}" / e["terminal_state"]["path"]).astype(np.int64)
    assert hashlib.sha256(state.astype("<i1").tobytes()).hexdigest() == e["terminal_state"]["sha256"]
    level = a.start_level
    origin_cell = state.reshape(-1, 12).mean(axis=1)  # the level-m0 settled field, to be tracked through every step
    out = {"schema": "oph.exploratory.refine-settle-cocycle.v1", "start_level": level, "start_seed": int(e["seed"]), "fresh_records": a.fresh, "inheritance": a.inherit,
           "fresh_law_note": "uniform5 = i.i.d. uniform on {0..4} per port (mean 2, keeps the child's mean at a half-integer); uniform6 = the engine's initial law {0..5}, "
                             "which puts the child's mean on an integer and makes the settlement transport-limited",
           "definition": "children inherit the parent's twelve readings; fresh i.i.d. loads {0..5} per port are added (or none); settle with the canonical integer law; "
                         "b = <settled cap excess, inherited cap excess> / <inherited, inherited>; covariance survival b^2; theta = -log2(b^2) per refinement step",
           "targets": {"P_over_24": round(P_STAR / 24, 6), "P_over_48": round(P_STAR / 48, 6)}, "steps": []}
    for step in range(a.steps):
        nxt = level + 1
        H.build_geometry(nxt, cache, log=lambda m: print(m, flush=True))
        geo = H.Geometry(cache, nxt); n = geo.carriers
        parent = state - int(state.min()) if a.inherit == "bit" else state
        assert parent.max() - parent.min() <= 1, "a settled parent reads two adjacent values"
        inherited = np.repeat(parent.reshape(-1, 12), 4, axis=0).reshape(-1)  # children 4c..4c+3 inherit carrier c's readings port by port
        rng = np.random.default_rng(31_000 + nxt)
        fresh = (rng.integers(0, 5, size=inherited.size) if a.fresh == "uniform5" else rng.integers(0, 6, size=inherited.size) if a.fresh == "uniform6"
                 else np.zeros_like(inherited))
        loads = (inherited + fresh).astype(np.int64)
        assert loads.max() <= 127
        seed = H.schedule_seed(nxt, 0)
        outdir = Path(a.out).parent / f"refine_L{nxt}"; outdir.mkdir(parents=True, exist_ok=True)
        res = settle(H, geo, loads, seed)
        np.save(outdir / "terminal_state.npy", res["state"])
        settled = res["state"].astype(np.int64)
        points = np.load(cache / f"L{nxt}" / "cell_points.npy")
        inh_cell = inherited.reshape(n, 12).mean(axis=1); set_cell = settled.reshape(n, 12).mean(axis=1); fresh_cell = fresh.reshape(n, 12).mean(axis=1)
        rows = survival_by_scale(points, inh_cell, set_cell, fresh_cell if a.fresh != "none" else None)
        # survival of the ORIGINAL level-m0 field through all steps so far (cocycle test)
        origin_here = np.repeat(origin_cell, 4 ** (step + 1))
        rows_origin = survival_by_scale(points, origin_here, set_cell, None)
        top = [r for r in rows if r["cells_per_cap"] >= 2000]
        entry = {"level": nxt, "carriers": n, "sweeps": res["sweeps"], "terminated": res["terminated"], "q": int(settled.min()), "seconds": res["seconds"],
                 "survival_by_scale": rows, "origin_survival_by_scale": rows_origin,
                 "super_horizon_theta_per_step": [r["theta_per_step"] for r in top], "super_horizon_b": [r["amplitude_b"] for r in top]}
        out["steps"].append(entry)
        print(f"step L{level}->L{nxt}: sweeps {res['sweeps']}, super-horizon b {entry['super_horizon_b']}, theta/step {entry['super_horizon_theta_per_step']}, "
              f"origin survival b {[r['amplitude_b'] for r in rows_origin if r['cells_per_cap'] >= 2000]}, {time.time() - t0:.0f}s", flush=True)
        state = settled; level = nxt
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    # cocycle summary: product of single-step super-horizon b^2 against the origin survival after all steps
    if len(out["steps"]) >= 2:
        prod = float(np.prod([np.mean(s["super_horizon_b"]) ** 2 for s in out["steps"]]))
        last = out["steps"][-1]["origin_survival_by_scale"]; origin_b2 = float(np.mean([r["amplitude_b"] ** 2 for r in last if r["cells_per_cap"] >= 2000]))
        out["cocycle"] = {"product_of_single_step_survivals": round(prod, 6), "origin_survival_after_all_steps": round(origin_b2, 6),
                          "ratio": round(origin_b2 / prod, 6) if prod > 0 else None, "theta_per_step_from_origin": round(float(-np.log2(origin_b2) / len(out["steps"])), 5) if origin_b2 > 0 else None}
        print("cocycle:", out["cocycle"], flush=True)
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

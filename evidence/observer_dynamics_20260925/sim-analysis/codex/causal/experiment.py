"""Frozen boost-clock experiment; run on the designated AWS host, one process.

The source graph is the unchanged source_net.py at commit 14d1699, vendored
with its sole local import. No existing simulation directories are mutated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "vendor"))
from oph_exact import source_net as S


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def lens_volume(r: float, s: float, d: float) -> float:
    """Intersection volume of two Euclidean 3-balls, from spherical caps."""
    if min(r, s) <= 0 or d >= r + s:
        return 0.0
    if d <= abs(r - s):
        return 4 * math.pi * min(r, s) ** 3 / 3
    return math.pi * (r + s - d) ** 2 * (d * d + 2 * d * (r + s) - 3 * (r - s) ** 2) / (12 * d)


def lens_volume_cross_section(r: float, s: float, d: float) -> float:
    """Independent 1D integration of circular cross sections along tip axis."""
    lo, hi = max(-r, d - s), min(r, d + s)
    if hi <= lo:
        return 0.0
    points = [] if d == 0 else [(r * r - s * s + d * d) / (2 * d)]
    points = [p for p in points if lo < p < hi]
    return quad(lambda x: math.pi * max(0, min(r * r - x * x, s * s - (x - d) ** 2)),
                lo, hi, points=points, epsabs=1e-13, epsrel=1e-12)[0]


def geometry_control(k: int, delta: float, displacement: float) -> dict:
    t = k * delta
    if not 0 <= displacement < t:
        raise ValueError("timelike endpoints required")
    exact = math.pi * (t * t - displacement * displacement) ** 2 / 24
    integrated = quad(lambda u: lens_volume(u, t - u, displacement), 0, t,
                      points=[(t - displacement) / 2, (t + displacement) / 2],
                      epsabs=1e-13, epsrel=1e-12)[0]
    mismatch = abs(integrated / exact - 1)
    if mismatch > 1e-10:
        raise ValueError("lens integration fails Lorentzian volume calibration")
    layer_volumes = [lens_volume(j * delta, (k - j) * delta, displacement) for j in range(k + 1)]
    sampled = delta * sum(layer_volumes)
    return {"continuum_volume": exact, "integrated_lens_volume": integrated,
            "calibration_relative_error": mismatch, "finite_layer_volume": sampled,
            "finite_layer_over_continuum_volume": sampled / exact,
            "finite_layer_slice_volumes": layer_volumes,
            "proper_time_ratio": math.sqrt(1 - (displacement / t) ** 2)}


def axis_graph(indptr, indices, sites):
    """Same population, deleting only reads that change >1 coordinates."""
    pieces, pointers = [], [0]
    for u in range(len(sites)):
        row = indices[indptr[u]:indptr[u + 1]]
        kept = row[np.count_nonzero(sites[row] != sites[u], axis=1) <= 1]
        pieces.append(kept)
        pointers.append(pointers[-1] + len(kept))
    return np.array(pointers, dtype=np.int64), np.concatenate(pieces)


def longest_read_history(indptr, indices, start: int, k: int) -> np.ndarray:
    """Dynamic programming on actual one-layer reads, not a formula from K."""
    score = np.full(len(indptr) - 1, -100000, dtype=np.int32)
    score[start] = 0
    for _ in range(k):
        score = np.maximum.reduceat(score[indices], indptr[:-1]) + 1
    return score


def witnessed_history(indptr, indices, x: int, y: int, to_y, k: int) -> list[dict]:
    remaining = int(to_y[x])
    if not 0 <= remaining <= k:
        return []
    sites = [x] * (k - remaining + 1)
    u = x
    while u != y:
        neighbors = indices[indptr[u]:indptr[u + 1]]
        candidates = neighbors[to_y[neighbors] == remaining - 1]
        if len(candidates) == 0:
            raise ValueError("BFS path cannot be witnessed")
        u = int(candidates.min())
        sites.append(u)
        remaining -= 1
    if len(sites) != k + 1:
        raise ValueError("wrong history length")
    for u, v in zip(sites, sites[1:]):
        if v not in indices[indptr[u]:indptr[u + 1]]:
            raise ValueError("history uses a nonexistent read")
    return [{"layer": j, "site": u} for j, u in enumerate(sites)]


def point_label(point, values, q):
    coords = [int(np.argmin(abs(values - x))) for x in point]
    return coords[0] * q * q + coords[1] * q + coords[2]


def interval_counts(distance_x, distance_y, k):
    return [int(np.count_nonzero((distance_x >= 0) & (distance_x <= j)
                                & (distance_y >= 0) & (distance_y <= k - j))) for j in range(k + 1)]


def run(spec: dict, output: Path) -> dict:
    start = time.monotonic()
    rows, graph_rows = [], []
    source_pins = {}
    for name, expected in spec["source_pins"].items():
        got = digest((ROOT / "vendor" / name).read_bytes())
        if got != expected:
            raise ValueError("vendored source differs from frozen pin")
        source_pins[name] = got
    result = {"schema": "oph.codex.causal-observer-clock-experiment.receipt.v1",
              "spec_sha256": digest((ROOT / "SPEC.json").read_bytes()),
              "script_sha256": digest(Path(__file__).read_bytes()),
              "source_pins_checked": source_pins, "source_commit": spec["source_commit"],
              "python": sys.version.split()[0], "numpy": np.__version__,
              "thread_env": {name: os.environ.get(name) for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
              "complete": False, "graphs": graph_rows, "rows": rows}

    def save():
        result["elapsed_seconds"] = time.monotonic() - start
        output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")

    for q in spec["q_values"]:
        if time.monotonic() - start > 600:
            break
        k, delta = math.floor(.9 * math.sqrt(q)), 1 / math.sqrt(q)
        t = k * delta
        tick = time.monotonic()
        indptr, indices, values, _, _, site_labels = S.build_site_graph(q, 3)
        coords_1d = np.array([S.phi_float(x) for x in values])
        positions = coords_1d[site_labels]
        main_indptr, main_indices = indptr, indices
        rest = point_label(np.full(3, .5), coords_1d, q)
        rest_control = geometry_control(k, delta, 0)
        for law in ("euclidean", "axis_only_control"):
            if law == "axis_only_control":
                indptr, indices = axis_graph(main_indptr, main_indices, site_labels)
            graph_rows.append({"q": q, "law": law, "sites": q ** 3, "K": k, "delta": delta,
                               "neighbor_entries": int(len(indices)),
                               "same_site_wait_present_at_every_site": all(u in indices[indptr[u]:indptr[u + 1]] for u in range(q ** 3))})
            distances, scores = {}, {}

            def fetch(u):
                if u not in distances:
                    distances[u] = S.bfs(indptr, indices, u, q ** 3, cap=k)
                    scores[u] = longest_read_history(indptr, indices, u, k)
                return distances[u]

            rest_distance = fetch(rest)
            reference_count = sum(interval_counts(rest_distance, rest_distance, k))
            reference_rank = int(scores[rest][rest])
            for direction in spec["directions"]:
                direction_unit = np.array(direction) / np.linalg.norm(direction)
                for beta in spec["beta_targets"]:
                    if time.monotonic() - start > 600:
                        save()
                        return result
                    x = point_label(np.full(3, .5) - beta * t * direction_unit / 2, coords_1d, q)
                    y = point_label(np.full(3, .5) + beta * t * direction_unit / 2, coords_1d, q)
                    vector = positions[y] - positions[x]
                    separation = float(np.linalg.norm(vector))
                    actual_beta = separation / t
                    row = {"q": q, "law": law, "direction": direction, "beta_target": beta,
                           "K": k, "delta": delta, "source": x, "target": y,
                           "source_position": positions[x].tolist(), "target_position": positions[y].tolist(),
                           "actual_beta": actual_beta, "reference_count": reference_count, "reference_rank_edges": reference_rank}
                    if actual_beta >= 1:
                        row["status"] = "not_timelike_after_tip_rounding"
                        rows.append(row)
                        continue
                    control = geometry_control(k, delta, separation)
                    midpoint = (positions[x] + positions[y]) / 2
                    half_width = .5 * np.sqrt(t * t - separation * separation + vector * vector)
                    clearance = float(min((midpoint - half_width).min(), (1 - midpoint - half_width).min()))
                    dx, dy = fetch(x), fetch(y)
                    counts = interval_counts(dx, dy, k)
                    total = sum(counts)
                    related = 0 <= int(dx[y]) <= k
                    row.update({"control": control, "counts_by_layer": counts, "event_count": total,
                                "continuum_clearance_from_cube": clearance, "graph_distance": int(dx[y]),
                                "actual_rank_edges": int(scores[x][y]) if related else None,
                                "status": "timelike_and_graph_related" if related else "timelike_but_graph_unrelated"})
                    if clearance < 0:
                        row["status"] = "continuum_diamond_clipped"
                    if related:
                        history = witnessed_history(indptr, indices, x, y, dy, k)
                        row["witness_history"] = history
                        row["history_has_K_distinct_layer_events_plus_one"] = len(history) == k + 1 and len({(e['layer'], e['site']) for e in history}) == k + 1
                        row["rank_clock_ratio"] = row["actual_rank_edges"] / reference_rank
                        raw = (total / reference_count) ** .25
                        correction = (control["finite_layer_over_continuum_volume"] / rest_control["finite_layer_over_continuum_volume"]) ** .25
                        corrected = raw / correction if correction > 0 else None
                        proper = control["proper_time_ratio"]
                        row.update({"volume_clock_ratio": raw, "finite_layer_correction": correction,
                                    "corrected_volume_clock_ratio": corrected,
                                    "rank_relative_error": row["rank_clock_ratio"] / proper - 1,
                                    "volume_relative_error": raw / proper - 1,
                                    "corrected_volume_relative_error": None if corrected is None else corrected / proper - 1})
                    rows.append(row)
            save()
        print(f"q={q}: {time.monotonic()-tick:.2f}s; {len(rows)} rows retained", flush=True)
    result["complete"] = len(rows) == len(spec["q_values"]) * 2 * len(spec["directions"]) * len(spec["beta_targets"])
    main = [r for r in rows if r["q"] == 34 and r["beta_target"] in [.25, .5, .75] and r["status"] == "timelike_and_graph_related"]
    summaries = {}
    for law in ("euclidean", "axis_only_control"):
        subset = [r for r in main if r["law"] == law]
        summaries[law] = {"rows": len(subset), "median_abs_rank_relative_error": float(np.median([abs(r["rank_relative_error"]) for r in subset])),
                          "median_abs_volume_relative_error": float(np.median([abs(r["volume_relative_error"]) for r in subset])),
                          "median_abs_corrected_volume_relative_error": float(np.median([abs(r["corrected_volume_relative_error"]) for r in subset if r["corrected_volume_relative_error"] is not None]))}
    valid = [r for r in rows if r["status"] == "timelike_and_graph_related"]
    result["summary"] = {"q34_primary_comparison": summaries,
                         "uniform_rank_clock_identification_rejected_at_5pct": any(abs(r["rank_relative_error"]) > .05 for r in valid if r["law"] == "euclidean" and r["beta_target"] >= .5),
                         "all_actual_rank_edges_equal_K": all(r["actual_rank_edges"] == r["K"] for r in valid),
                         "all_history_witnesses_valid": all(r["history_has_K_distinct_layer_events_plus_one"] for r in valid),
                         "calibration_max_relative_error": max(r["control"]["calibration_relative_error"] for r in rows if "control" in r)}
    save()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "receipt.json")
    args = parser.parse_args()
    spec_bytes = (ROOT / "SPEC.json").read_bytes()
    if digest(spec_bytes) != (ROOT / "SPEC.sha256").read_text().split()[0]:
        raise SystemExit("frozen specification hash mismatch")
    result = run(json.loads(spec_bytes), args.output)
    print(json.dumps(result.get("summary", {"partial": True}), indent=2))


if __name__ == "__main__":
    main()

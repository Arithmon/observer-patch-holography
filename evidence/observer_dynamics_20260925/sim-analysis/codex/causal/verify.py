"""Independent bounded verification using KD-tree edges and SciPy BFS.

Does not import the experiment or its source graph implementation. Runs on AWS;
the q34 graph is rebuilt in memory, then discarded. No raw arrays are saved.
"""
from __future__ import annotations
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b, label, tolerance=2e-10):
    if not math.isclose(a, b, rel_tol=tolerance, abs_tol=1e-12):
        raise ValueError(f"{label}: {a} != {b}")


def intersection(r, s, d):
    lo, hi = max(-r, d - s), min(r, d + s)
    if hi <= lo:
        return 0.0
    switch = (r * r - s * s + d * d) / (2 * d) if d else 0
    return quad(lambda x: math.pi * max(0, min(r * r - x * x, s * s - (x - d) ** 2)), lo, hi,
                points=[switch] if lo < switch < hi else None, epsabs=1e-13, epsrel=1e-12)[0]


def main():
    start = time.monotonic()
    spec = json.loads((ROOT / "SPEC.json").read_text())
    receipt = json.loads((ROOT / "receipt.json").read_text())
    assert sha(ROOT / "SPEC.json") == receipt["spec_sha256"] == (ROOT / "SPEC.sha256").read_text().split()[0]
    assert sha(ROOT / "experiment.py") == receipt["script_sha256"]
    assert receipt["complete"]
    keys = {(r["q"], r["law"], tuple(r["direction"]), r["beta_target"]) for r in receipt["rows"]}
    expected = set(itertools.product(spec["q_values"], ["euclidean", "axis_only_control"], map(tuple, spec["directions"]), spec["beta_targets"]))
    assert len(keys) == len(receipt["rows"]) == len(expected) and keys == expected
    checks = []
    for q in spec["q_values"]:
        values = np.array([b * (1 + math.sqrt(5)) / 2 % 1 for b in range(q)])
        sites = np.array(list(itertools.product(range(q), repeat=3)), dtype=np.int32)
        points = values[sites]
        tree = cKDTree(points)
        radius = 1 / math.sqrt(q)
        low = tree.query_ball_point(points, radius * (1 - 1e-10), return_length=True, workers=1)
        high = tree.query_ball_point(points, radius * (1 + 1e-10), return_length=True, workers=1)
        assert np.array_equal(low, high), "floating edge decisions too close to boundary for independent check"
        neighbors = tree.query_ball_point(points, radius, return_sorted=True, workers=1)
        k = math.floor(.9 * math.sqrt(q))
        t = k * radius
        for law in ["euclidean", "axis_only_control"]:
            rows = [r for r in receipt["rows"] if r["q"] == q and r["law"] == law]
            if law == "axis_only_control":
                neighbor_rows = [np.array(row, dtype=np.int32)[np.count_nonzero(sites[row] != sites[u], axis=1) <= 1] for u, row in enumerate(neighbors)]
            else:
                neighbor_rows = neighbors
            ptr = np.r_[0, np.cumsum([len(row) for row in neighbor_rows])].astype(np.int64)
            indices = np.concatenate(neighbor_rows).astype(np.int32)
            graph = csr_matrix((np.ones(len(indices), dtype=np.int8), indices, ptr), shape=(q**3, q**3))
            advertised = next(g for g in receipt["graphs"] if g["q"] == q and g["law"] == law)
            assert len(indices) == advertised["neighbor_entries"]
            cache = {}

            def distance(u):
                if u not in cache:
                    order, parents = breadth_first_order(graph, u, directed=True, return_predecessors=True)
                    dist = np.full(q**3, -1, dtype=np.int16)
                    dist[u] = 0
                    for v in order[1:]:
                        dist[v] = dist[parents[v]] + 1
                    cache[u] = dist
                return cache[u]

            center_label = int(np.argmin(abs(values - .5)))
            rest = center_label * (q * q + q + 1)
            drest = distance(rest)
            rest_count = sum(int(np.count_nonzero((drest >= 0) & (drest <= min(j, k-j)))) for j in range(k+1))
            rest_volume = radius * sum(intersection(j*radius, (k-j)*radius, 0) for j in range(k+1))
            for row in rows:
                x, y = row["source"], row["target"]
                close(float(np.linalg.norm(points[y] - points[x])) / t, row["actual_beta"], "actual beta")
                direction = np.array(row["direction"], dtype=float)
                direction /= np.linalg.norm(direction)
                target_positions = [np.full(3, .5) + sign * row["beta_target"] * t * direction/2 for sign in [-1, 1]]
                labels = []
                for point in target_positions:
                    cs = [int(np.argmin(abs(values-z))) for z in point]
                    labels.append(cs[0]*q*q+cs[1]*q+cs[2])
                assert labels == [x, y], "tip selection changed"
                if row["actual_beta"] >= 1:
                    assert row["status"] == "not_timelike_after_tip_rounding"
                    continue
                dx, dy = distance(x), distance(y)
                counts = [int(np.count_nonzero((dx >= 0)&(dx <= j)&(dy >= 0)&(dy <= k-j))) for j in range(k+1)]
                assert counts == row["counts_by_layer"] and sum(counts) == row["event_count"]
                assert rest_count == row["reference_count"]
                graph_related = 0 <= dx[y] <= k
                if not graph_related:
                    assert row["status"] == "timelike_but_graph_unrelated"
                    assert sum(counts) == 0
                    continue
                history = row["witness_history"]
                assert len(history) == k+1
                assert [e["layer"] for e in history] == list(range(k+1))
                assert history[0]["site"] == x and history[-1]["site"] == y
                for a,b in zip(history,history[1:]):
                    assert graph[a["site"],b["site"]] == 1
                assert row["actual_rank_edges"] == k == row["reference_rank_edges"]
                # At most one strict event per layer; witness attains the bound.
                close(row["rank_clock_ratio"], 1.0, "rank ratio")
                d = float(np.linalg.norm(points[y]-points[x]))
                tau = math.sqrt(1-(d/t)**2)
                volume = math.pi*(t*t-d*d)**2/24
                sampled = radius*sum(intersection(j*radius,(k-j)*radius,d) for j in range(k+1))
                correction = ((sampled/volume)/(rest_volume/(math.pi*t**4/24)))**.25
                raw = (sum(counts)/rest_count)**.25
                close(row["control"]["continuum_volume"], volume,"continuum volume")
                close(row["control"]["finite_layer_volume"],sampled,"independent layer volume")
                close(row["volume_clock_ratio"],raw,"count clock")
                close(row["rank_relative_error"],1/tau-1,"rank error")
                close(row["corrected_volume_relative_error"],raw/correction/tau-1,"corrected error")
            checks.append({"q":q,"law":law,"rows_verified":len(rows),"independent_neighbor_entries":len(indices),
                           "no_float_edge_ambiguity_in_relative_radius_band":1e-10})
        print(f"independent verification q={q} complete",flush=True)
    output = {"schema":"oph.codex.causal-observer-clock-independent-verification.v1", "verified":True,
              "receipt_sha256":sha(ROOT/'receipt.json'),"spec_sha256":sha(ROOT/'SPEC.json'),
              "verifier_sha256":sha(Path(__file__)), "method":"independent KD-tree graph, scipy BFS, exact layer-count upper bound plus edge witnesses, numerical cross-section integration",
              "checks":checks,"rows_verified":len(receipt['rows']),"elapsed_seconds":time.monotonic()-start,
              "scope":"finite deterministic instrument verification, not physical identification or continuum proof"}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()

"""Census of carrier states in a settled federation.

For every carrier of a terminal state the twelve readings are in {q, q+1} (balanced class); the
carrier's state is the subset of ports holding q+1.  Two carriers are in the same *class* when
their subsets are related by one of the sixty rotations of the icosahedral port graph (A5), so
the census counts the settled world's carriers by A5-orbit: how many ports are raised (k) and
which arrangement, up to symmetry.  Also reported: the raised-port count distribution against
the binomial expectation of independent ports at the same mean, the distribution at the twelve
pentagonal-vertex cells against the rest, and the slow-band norm by class.

usage: python3 carrier_census.py SIM_ROOT RUN_DIR LEVEL OUT.json [--terminals N]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from math import comb
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sim_root")
    parser.add_argument("run_dir")
    parser.add_argument("level", type=int)
    parser.add_argument("out")
    parser.add_argument("--terminals", type=int, default=4)
    parser.add_argument("--cache", default=None, help="geometry cache (for the vertex cells); default RUN_DIR/../../geo or the receipt's cache_dir")
    args = parser.parse_args()
    sys.path.insert(0, args.sim_root)
    from oph_exact import carrier
    from oph_exact.federation_archive import vertex_cells

    rots = np.asarray(carrier.rotations(), dtype=np.int64)  # (60, 12) port permutations
    assert rots.shape == (60, 12)
    # canonical orbit id of every 12-bit pattern: the minimum image over the sixty rotations
    patterns = np.arange(4096, dtype=np.int64)
    bits = ((patterns[:, None] >> np.arange(12)) & 1).astype(np.int64)  # (4096, 12)
    images = np.zeros((60, 4096), dtype=np.int64)
    for r, perm in enumerate(rots):
        moved = bits[:, perm]  # port p of the image reads bit perm[p]
        images[r] = (moved << np.arange(12)).sum(axis=1)
    canon = images.min(axis=0)
    orbit_ids, orbit_index = np.unique(canon, return_inverse=True)
    orbit_size = np.bincount(orbit_index, minlength=orbit_ids.size)
    k_of_orbit = np.array([bin(int(o)).count("1") for o in orbit_ids])
    p_slow = carrier.slow_band_projector()

    run = Path(args.run_dir)
    receipt = json.loads((run / "receipt.json").read_text())
    n = int(receipt["carriers"])
    cache = Path(args.cache) if args.cache else Path(receipt["geometry"]["cache_dir"]).parent
    vtx = np.asarray(vertex_cells(args.level)).ravel()
    vtx = np.unique(vtx[vtx >= 0])
    is_vertex = np.zeros(n, dtype=bool)
    is_vertex[vtx] = True
    out = {"level": args.level, "carriers": n, "orbits_total": int(orbit_ids.size), "orbits_by_k": {str(k): int(np.sum(k_of_orbit == k)) for k in range(13)}, "states": []}
    for e in receipt["integer_law"]["entries"][: args.terminals]:
        x = np.load(run / f"integer_{e['seed']}" / e["terminal_state"]["path"]).astype(np.int64)
        if hashlib.sha256(x.astype("<i1").tobytes()).hexdigest() != e["terminal_state"]["sha256"]:
            raise ValueError("terminal state digest")
        blocks = x.reshape(n, 12)
        q = int(blocks.min())
        raised = (blocks == q + 1).astype(np.int64)
        assert np.all((blocks == q) | (blocks == q + 1)), "not in the balanced class"
        code = (raised << np.arange(12)).sum(axis=1)
        k = raised.sum(axis=1)
        orb = orbit_index[code]
        counts = np.bincount(orb, minlength=orbit_ids.size)
        p = float(raised.mean())
        k_hist = np.bincount(k, minlength=13)
        binom = np.array([comb(12, j) * p**j * (1 - p) ** (12 - j) for j in range(13)]) * n
        # class census: top classes by count, with their expected count under independent ports (orbit size x p^k (1-p)^(12-k) x n)
        expected = orbit_size * (p ** k_of_orbit) * ((1 - p) ** (12 - k_of_orbit)) * n
        top = np.argsort(counts)[::-1][:20]
        slow = np.sqrt(np.sum(((raised - raised.mean(axis=1, keepdims=True)) @ p_slow) ** 2, axis=1))
        out["states"].append({
            "seed": int(e["seed"]), "q": q, "raised_fraction": round(p, 6),
            "raised_count_histogram": k_hist.tolist(), "binomial_expectation": [round(float(v), 1) for v in binom],
            "raised_count_chi2_vs_binomial": round(float(np.sum((k_hist - binom) ** 2 / np.where(binom > 0, binom, 1))), 2),
            "raised_fraction_vertex_cells": round(float(raised[is_vertex].mean()), 5), "raised_fraction_other_cells": round(float(raised[~is_vertex].mean()), 5),
            "vertex_cells": int(is_vertex.sum()),
            "occupied_orbits": int(np.count_nonzero(counts)),
            "top_classes": [{"orbit_pattern": int(orbit_ids[o]), "k": int(k_of_orbit[o]), "orbit_size": int(orbit_size[o]), "count": int(counts[o]), "expected_independent": round(float(expected[o]), 1),
                             "slow_norm": round(float(slow[orb == o].mean()), 4) if counts[o] else None} for o in top],
            "orbit_count_chi2_vs_independent_over_occupied": round(float(np.sum(((counts - expected) ** 2 / np.where(expected > 0, expected, 1))[expected > 5])), 1),
            "slow_norm_mean": round(float(slow.mean()), 5), "slow_norm_zero_fraction": round(float(np.mean(slow < 1e-9)), 5),
        })
        print(f"L{args.level} seed {e['seed']}: q={q}, raised fraction {p:.4f}, k histogram {k_hist.tolist()}, chi2 vs binomial {out['states'][-1]['raised_count_chi2_vs_binomial']}, "
              f"occupied orbits {out['states'][-1]['occupied_orbits']}/{orbit_ids.size}, vertex raised {out['states'][-1]['raised_fraction_vertex_cells']} vs other {out['states'][-1]['raised_fraction_other_cells']}", flush=True)
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

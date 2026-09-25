"""Independent checks of the sampled chain-count and interval-spectrum receipt.

Standard library and numpy only for the arithmetic; the small-q rows are rebuilt from the
definition of the order (graph distance on the neighbour graph, one event per site and layer)
with a brute-force relation matrix, the continuum references are regenerated from their seeds,
and every derived field of every region (falling factorials, flat coefficients, inverted
dimensions, chain ratios, spectrum distances, closest dimension) is recomputed from the stored
counts.  The sampled draws themselves are not replayed above the brute-force levels; their
standard errors are the receipt's own.

usage: python3 -m oph_exact.verify_source_net_manifold_sampled_independent [RECEIPT] [--brute-force-levels 5 6]
"""

from __future__ import annotations

import argparse
import json
import sys
from math import gamma as Gamma
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "data/exact/source_net_manifold_sampled_receipt.json"


def chi(k: int, d: float) -> float:
    return (1.0 / k) * (Gamma(d + 1) / 2.0) ** (k - 1) * Gamma(d / 2.0) * Gamma(d) / (Gamma(k * d / 2.0) * Gamma((k + 1) * d / 2.0))


def invert(k: int, value: float):
    if not value > 0.0:
        return None
    lo, hi = 1.05, 12.0
    if value >= chi(k, lo) or value <= chi(k, hi):
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if chi(k, mid) > value:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def falling(N: int, k: int) -> int:
    out = 1
    for i in range(k):
        out *= N - i
    return out


def close(a, b, tol=1e-9) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(a)), abs(float(b)))


def check_region(row: dict, dim: int, references: dict, grid, failures: list, where: str) -> None:
    N = int(row["event_count"])
    d = dim + 1
    for k in ("2", "3", "4"):
        ch = row["chains"][k]
        est = ch["count"]["estimate"]
        ff = falling(N, int(k))
        ratio = est / ff if ff > 0 else None
        if not close(ch["count_over_falling_factorial"], ratio, 1e-9):
            failures.append(f"{where}: C{k} ratio")
        if not close(ch["flat_chi"], chi(int(k), d), 1e-9):
            failures.append(f"{where}: chi_{k}")
        if ratio is not None and not close(ch["relative_deviation_from_flat"], ratio / chi(int(k), d) - 1.0, 1e-8):
            failures.append(f"{where}: C{k} deviation")
        inv = None if ratio is None else invert(int(k), ratio)
        if not close(ch["inverted_dimension"], inv, 1e-7):
            failures.append(f"{where}: inverted dimension C{k}")
    c2, c3, c4 = (row["chains"][k]["count"]["estimate"] for k in ("2", "3", "4"))
    r = row["chain_ratios"]
    if c2 > 0 and not (close(r["C3_N_over_C2_squared"], c3 * N / c2 ** 2, 1e-8) and close(r["C4_N2_over_C2_cubed"], c4 * N * N / c2 ** 3, 1e-8)):
        failures.append(f"{where}: chain ratios")
    if not (close(r["flat_C3_N_over_C2_squared"], chi(3, d) / chi(2, d) ** 2, 1e-9) and close(r["flat_C4_N2_over_C2_cubed"], chi(4, d) / chi(2, d) ** 3, 1e-9)):
        failures.append(f"{where}: flat chain ratios")
    sp = row["interval_spectrum"]
    if not close(sp["flat_mean_chi3_over_chi2"], chi(3, d) / chi(2, d), 1e-9):
        failures.append(f"{where}: flat spectrum mean")
    if not close(sp["relative_deviation_of_mean_from_flat"], sp["weighted_mean"] / (chi(3, d) / chi(2, d)) - 1.0, 1e-8):
        failures.append(f"{where}: spectrum mean deviation")
    cdf = sp["cdf_on_grid"]
    if any(a > b + 1e-12 for a, b in zip(cdf, cdf[1:])) or len(cdf) != len(grid):
        failures.append(f"{where}: spectrum cdf not monotone")
    dist = {}
    for key, ref in references.items():
        ks = max(abs(a - b) for a, b in zip(cdf, ref["cdf_on_grid"]))
        dist[key] = ks
        if not close(sp["distance_to_flat_references"][key]["max_cdf_difference_on_grid"], ks, 1e-9):
            failures.append(f"{where}: distance to d={key}")
    closest = min(dist, key=lambda k_: dist[k_])
    if int(closest) != int(sp["closest_flat_spacetime_dimension"]):
        failures.append(f"{where}: closest dimension")


def brute_force_rows(n: int, dim: int) -> list:
    """Exact chain counts of the centre (and moving) diamond from the order's definition."""
    from oph_exact import source_net_manifold_sampled as M  # geometry only: sites, graph, tips
    from oph_exact.source_net import bfs
    fam = M.Family(n, dim)
    K = fam.K
    out = []
    regions = [("centre", fam.centre, fam.centre)]
    moving = fam.moving_tips() if dim >= 2 else None
    if moving is not None:
        regions.append(("moving", moving["x"], moving["y"]))
    D = np.stack([bfs(fam.indptr, fam.indices, s, fam.count) for s in range(fam.count)]).astype(np.int64)
    for name, x, y in regions:
        A = D[x]
        G = D[y]
        events = [(j, s) for s in range(fam.count) for j in range(int(A[s]), K - int(G[s]) + 1)]
        J = np.array([e[0] for e in events])
        Sx = np.array([e[1] for e in events])
        R = ((J[None, :] > J[:, None]) & (D[Sx[:, None], Sx[None, :]] <= (J[None, :] - J[:, None]))).astype(np.int64)
        R2 = R @ R
        out.append({"region": name, "N": len(events), "C2": int(R.sum()), "C3": int(R2.sum()), "C4": int((R2 @ R).sum()),
                    "intervals": np.sort(R2[R.astype(bool)])})
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", nargs="?", type=Path, default=DEFAULT)
    parser.add_argument("--brute-force-levels", type=int, nargs="*", default=[5, 6])
    parser.add_argument("--reference-pairs", type=int, default=None, help="regenerate the continuum references (slow at 10^6 pairs)")
    args = parser.parse_args(argv)
    receipt = json.loads(args.receipt.read_text())
    failures: list = []
    grid = receipt["spectrum_grid"]
    refs = receipt["continuum_references"]
    for key, ref in refs.items():
        d = int(key)
        if not close(ref["ordering_fraction_exact"], 2.0 * chi(2, d), 1e-9) or not close(ref["mean_exact_chi3_over_chi2"], chi(3, d) / chi(2, d), 1e-9):
            failures.append(f"reference d={key}: closed forms")
        # the mean of (tau/T)^d has a heavy tail in high d: 5% on the mean, 2% on the ordering fraction, at 10^5 pairs or more
        if abs(ref["mean"] / (chi(3, d) / chi(2, d)) - 1.0) > 0.05 or abs(ref["ordering_fraction_monte_carlo"] / (2 * chi(2, d)) - 1.0) > 0.02:
            failures.append(f"reference d={key}: Monte Carlo off its closed form")
    if args.reference_pairs:
        from oph_exact import source_net_manifold_sampled as M
        for key, ref in refs.items():
            again = M.continuum_spectrum(int(key), ref["pairs"] if args.reference_pairs < 0 else args.reference_pairs, ref["seed"], grid)
            if args.reference_pairs < 0 and again["cdf_on_grid"] != ref["cdf_on_grid"]:
                failures.append(f"reference d={key}: not reproduced from its seed")
    for lv in receipt["levels"]:
        for fam in lv["families"]:
            dim = fam["dimension"]
            rows = list(fam["regions"]) + (fam["homogeneity"]["diamonds"] if fam.get("homogeneity") else [])
            for row in rows:
                check_region(row, dim, refs, grid, failures, f"q={fam['q']} dim={dim} {row['region']}")
            if lv["fibonacci_index"] in args.brute_force_levels:
                bf = {r["region"]: r for r in brute_force_rows(lv["fibonacci_index"], dim)}
                for row in fam["regions"]:
                    b = bf.get(row["region"])
                    if b is None:
                        continue
                    if int(row["event_count"]) != b["N"]:
                        failures.append(f"q={fam['q']} dim={dim} {row['region']}: event count")
                    ex = row.get("exact")
                    if ex is not None and (ex["C2"], ex["C3"], ex["C4"]) != (b["C2"], b["C3"], b["C4"]):
                        failures.append(f"q={fam['q']} dim={dim} {row['region']}: exact chain counts")
                    for k, key in ((2, "2"), (3, "3"), (4, "4")):
                        est, se = row["chains"][key]["count"]["estimate"], row["chains"][key]["count"]["standard_error"]
                        if abs(est - b["C" + key]) > 5.0 * se + 1e-9:
                            failures.append(f"q={fam['q']} dim={dim} {row['region']}: sampled C{k} {est} vs exact {b['C' + key]} (SE {se})")
    for line in failures:
        print("FAIL", line)
    print("MANIFOLD_SAMPLED_VERIFY", "PASS" if not failures else f"FAIL ({len(failures)})")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())

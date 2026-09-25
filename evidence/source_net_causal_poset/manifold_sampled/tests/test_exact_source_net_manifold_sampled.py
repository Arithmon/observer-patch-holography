"""Tests of the sampled chain-count and interval-spectrum lane against brute force at small q."""

from __future__ import annotations

import tempfile
from math import isclose
from pathlib import Path

import numpy as np
import pytest

from oph_exact import source_net_manifold_sampled as M
from oph_exact.source_net import bfs


def brute_force(fam: M.Family, K: int, x: int, y: int) -> dict:
    """Strict relation matrix of the diamond from the definition; chains by matrix products."""
    n = fam.count
    A = fam.distances_from(x)
    G = A if y == x else fam.distances_from(y)
    D = np.stack([bfs(fam.indptr, fam.indices, s, n) for s in range(n)]).astype(np.int64)
    events = [(j, s) for s in range(n) for j in range(int(A[s]), K - int(G[s]) + 1)]
    N = len(events)
    J = np.array([e[0] for e in events]); Sx = np.array([e[1] for e in events])
    R = (J[None, :] > J[:, None]) & (D[Sx[:, None], Sx[None, :]] <= (J[None, :] - J[:, None]))
    R = R.astype(np.int64)
    R2 = R @ R
    return {"N": N, "C2": int(R.sum()), "C3": int(R2.sum()), "C4": int((R2 @ R).sum()),
            "intervals": np.sort(R2[R.astype(bool)])}


def test_chain_coefficients_and_inversion():
    assert isclose(M.chi(2, 4), 0.05, rel_tol=1e-12)
    assert isclose(2 * M.chi(2, 3), 8 / 35, rel_tol=1e-12)
    assert isclose(M.chi(2, 2), 0.25, rel_tol=1e-12)
    assert isclose(M.chi(3, 2), 1 / 36, rel_tol=1e-12)
    assert isclose(M.chi(3, 4) / M.chi(2, 4), 1 / 105, rel_tol=1e-12)
    assert isclose(M.invert_chi(2, 0.05), 4.0, abs_tol=1e-9)
    assert isclose(M.invert_chi(3, M.chi(3, 3)), 3.0, abs_tol=1e-9)
    assert M.invert_chi(2, 0.0) is None
    assert M.falling(10, 3) == 720


def test_continuum_spectrum_means_match_chain_ratios():
    ref4 = M.continuum_spectrum(4, 200_000, 7, M.SPECTRUM_GRID)
    assert abs(ref4["ordering_fraction_monte_carlo"] - 0.1) < 0.003
    assert abs(ref4["mean"] / (1 / 105) - 1.0) < 0.04
    ref2 = M.continuum_spectrum(2, 100_000, 8, M.SPECTRUM_GRID)
    assert abs(ref2["ordering_fraction_monte_carlo"] - 0.5) < 0.01
    assert abs(ref2["mean"] / (1 / 9) - 1.0) < 0.03
    cdf = ref4["cdf_on_grid"]
    assert all(a <= b + 1e-12 for a, b in zip(cdf, cdf[1:])) and cdf[-1] == 1.0


@pytest.mark.parametrize("n,dim", [(5, 3), (6, 2), (6, 3)])
def test_exact_enumeration_equals_brute_force(n, dim):
    fam = M.Family(n, dim)
    K = fam.K
    A = fam.distances_from(fam.centre)
    ex = M.exact_region(fam, K, A, A)
    bf = brute_force(fam, K, fam.centre, fam.centre)
    assert ex["event_count"] == bf["N"]
    assert (ex["C2"], ex["C3"], ex["C4"]) == (bf["C2"], bf["C3"], bf["C4"])
    assert np.array_equal(np.sort(ex["intervals"]), bf["intervals"])
    moving = fam.moving_tips()
    if moving is not None:
        Am, Gm = fam.distances_from(moving["x"]), fam.distances_from(moving["y"])
        ex = M.exact_region(fam, K, Am, Gm)
        bf = brute_force(fam, K, moving["x"], moving["y"])
        assert (ex["event_count"], ex["C2"], ex["C3"], ex["C4"]) == (bf["N"], bf["C2"], bf["C3"], bf["C4"])


def test_sampled_estimators_agree_with_exact_counts():
    fam = M.Family(6, 3)  # q = 8, 512 sites, K = 3
    K = fam.K
    with tempfile.TemporaryDirectory() as tmp:
        fam.open_pool(1, Path(tmp))
        refs = {str(d): M.continuum_spectrum(d, 20_000, 100 + d, M.SPECTRUM_GRID) for d in M.CONTINUUM_DIMENSIONS}
        row = M.run_region(fam, "centre", K, fam.centre, fam.centre, 1500, 11, Path(tmp), "test", refs, M.SPECTRUM_GRID, True)
        fam.close_pool()
    ex = row["exact"]
    for k in ("2", "3", "4"):
        est, se = row["chains"][k]["count"]["estimate"], row["chains"][k]["count"]["standard_error"]
        assert abs(est - ex["C" + k]) <= 4.0 * se + 1e-9, (k, est, se, ex["C" + k])
    # the interval sum estimates C_3 as well, and the weighted spectrum mean is C_3 / (C_2 N)
    N = row["event_count"]
    assert abs(row["interval_spectrum"]["weighted_mean"] - ex["C3"] / (ex["C2"] * N)) < 0.15 * ex["C3"] / (ex["C2"] * N)
    assert row["chains"]["2"]["inverted_dimension"] is not None


def test_region_events_and_draws():
    fam = M.Family(5, 3)
    A = fam.distances_from(fam.centre)
    layers, N = M.region_events(fam.K, A, A)
    assert int(N) == sum(max(fam.K - 2 * int(a) + 1, 0) for a in A)
    ev = M.draw_events(layers, A, 200, np.random.default_rng(0))
    assert np.all(ev[:, 0] >= A[ev[:, 1]]) and np.all(ev[:, 0] <= fam.K - A[ev[:, 1]])

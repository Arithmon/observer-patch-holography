"""Chain counts and interval spectrum of the source-net order at large q, from sampled pairs.

The causal-limit lane (``oph_exact/source_net.py``) counts the strict pairs of the centre diamond
exactly and reads one number per diamond from them: the ordering fraction, its Myrheim-Meyer
inversion and the count clock.  The manifold-observations lane
(``oph_exact/manifold_observations.py``) adds interval-abundance profiles, the Benincasa-Dowker
action, boosts, off-centre tips and link directions on the same order, but it needs dense distance
tables, which caps it at ``q = 34``.  This lane reads the statistics of the order that take three
and four events at a time, at every ``q`` whose neighbour graph fits in memory (``q = 55, 89, 144``
on one box), from breadth-first searches alone:

1. the chain counts ``C_2, C_3, C_4`` of a diamond (numbers of related pairs, of chains
   ``x < y < z`` and of chains of four events), each with its standard error, against the flat
   ``(dim+1)``-dimensional expectations ``C_k = N^(k) chi_k(d)`` (Meyer 1988; ``N^(k)`` the
   falling factorial, ``2 chi_2`` the Myrheim-Meyer ordering fraction) and the dimension each of
   them inverts to;
2. the interval spectrum: the distribution of ``|I(y, z)| / N`` over the related pairs ``(y, z)`` of
   the diamond, against the continuum distribution of ``(tau / T)^(dim+1)`` over the related pairs
   of the flat diamond of proper duration ``T`` (its mean is ``chi_3 / chi_2``, ``1/105`` in
   ``3+1`` dimensions), with the flat references in every spacetime dimension from two to six;
3. the same readouts on the moving-tip diamond of the causal-limit lane (boost) and on interior
   diamonds at off-centre tips (homogeneity);
4. the two- and one-dimensional golden control populations under the same read law.

Estimators.  Events ``y = (j, s)`` are drawn uniformly from the diamond.  One breadth-first search
from ``s`` with the causal-limit lane's exact cone pruning gives ``|past(y)|`` and ``|fut(y)|``
inside the diamond; a partner ``z = (m, u)`` is drawn uniformly from ``fut(y)``; a second search from
``u`` gives ``|fut(z)|``, ``|past(z)|`` and, together with the first, the interval ``|I(y, z)|``.
Then ``C_2 = N E|fut(y)|``, ``C_3 = N E[|past(y)| |fut(y)|]``, ``C_4 = N E[|fut(y)| |past(y)|
|fut(z)|]`` (the pair ``(y, z)`` has probability ``1 / (N |fut(y)|)``), ``sum |I(y, z)|`` over
related pairs is a second estimate of ``C_3``, and the interval spectrum is the ``|fut(y)|``-weighted
distribution of ``|I(y, z)| / N``.  With ``--exact`` every event and every partner is enumerated from
an all-pairs distance table, which is the exact count (the tests compare it with a brute-force
relation matrix at ``q = 5`` and ``q = 8``).

The population, the read law, the layer duration and one event per site and layer are supplied,
as in the theory.  The flat-diamond expectations are the comparison values of the declared limit.
Nothing here selects the population or the read law from native repairs, identifies a physical
clock or spacetime, or demonstrates the asymptotic limit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import sys
import tempfile
import time
from math import gamma as Gamma
from math import log, sqrt
from pathlib import Path

import numpy as np

from oph_exact import source_net as S
from oph_exact.source_net import (
    bfs,
    build_site_graph,
    canonical,
    ceil_sqrt,
    cone_masks,
    diamond_volume,
    digest,
    ellipsoid_buffer,
    fibonacci,
    file_sha256,
    metric_tables,
    orbit,
    phi_float,
    phi_scale,
    phi_sign,
    phi_square,
    phi_sub,
    rounded,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/exact/source_net_manifold_sampled_receipt.json"
SCHEMA = "oph.exact.source-net-manifold-sampled.v1"
LEVELS = (10, 11, 12)  # q = 55, 89, 144
DIMENSIONS = (3, 2, 1)
EVENT_SAMPLES = {"centre": 4000, "moving": 3000, "tip": 800}
TIP_DIRECTIONS = 14  # six axis directions and eight body diagonals
SAMPLE_SEED = 20260925
CONTINUUM_PAIRS = 1_000_000
CONTINUUM_SEED = 4040
CONTINUUM_DIMENSIONS = (2, 3, 4, 5, 6)
SPECTRUM_GRID = [float(v) for v in np.logspace(-6, 0, 25)]
EXACT_MAX_SITES = 10_000
POOL_MINIMUM_SITES = 4000
BATCH_SIZE = 16
PINS = (
    "oph_exact/source_net_manifold_sampled.py",
    "oph_exact/verify_source_net_manifold_sampled_independent.py",
    "tests/test_exact_source_net_manifold_sampled.py",
    "oph_exact/source_net.py",
)
SCOPE = {
    "population_supplied": True,
    "read_law_supplied": True,
    "tick_supplied": True,
    "one_event_per_site_and_layer_supplied": True,
    "exact_gram_metric_edge_decisions": True,
    "exact_graph_distance_order": True,
    "readouts_from_order_and_cardinality_only": True,
    "pair_sampling_declared_seeds": True,
    "flat_diamond_expectations_are_the_declared_limit": True,
    "native_repair_selected": False,
    "physical_clock_or_spacetime_identified": False,
    "finite_runs_demonstrate_asymptotic_limit": False,
    "dimension_statistic_used_as_acceptance": False,
    "poisson_sprinkling": False,
}
CLAIM_BOUNDARY = (
    "Sampled readouts of the declared source-record family (q^3 golden sites on the rank-three "
    "source Gram metric, complete-neighbour reads inside a_q = L/sqrt(q), one event per site and "
    "layer, layer duration a_q/c) at q <= 144 from the event order and event cardinalities alone: "
    "chain counts of two, three and four events with their standard errors, the dimension each "
    "inverts to, and the interval spectrum of the centre, moving-tip and off-centre diamonds, "
    "against the flat-diamond expectations of the declared limit and the two- and one-dimensional "
    "golden controls. The readouts are finite diagnostics of the supplied law; they do not select "
    "the population or the read law from native repairs, identify a physical clock or spacetime, "
    "or demonstrate the asymptotic limit of the theory's propositions."
)
CONVENTIONS = {
    "event_order": "(j,s) <= (j',t) iff graph distance d(s,t) <= j'-j on the exact site graph, waiting included",
    "finite_schedule": "K_q = ceil(sqrt q); T_q = K_q L/sqrt q; density rho_q = q^(dim+1/2)/L^(dim+1) with c = 1",
    "diamond": "inclusive interval between (0,x) and (K',y): events (j,s) with d(x,s) <= j <= K' - d(s,y)",
    "chains": "C_k = number of chains of k distinct events of the diamond under the strict order; "
              "C_2 is the strict pair count of the causal-limit lane",
    "flat_expectation": "C_k = N (N-1) ... (N-k+1) chi_k(d) with chi_k(d) = (1/k) (Gamma(d+1)/2)^(k-1) "
                        "Gamma(d/2) Gamma(d) / (Gamma(k d/2) Gamma((k+1) d/2)) for d = dim+1 (Meyer 1988); "
                        "2 chi_2(d) is the Myrheim-Meyer ordering fraction",
    "interval_spectrum": "distribution of |I(y,z)|/N over related pairs y < z of the diamond, |I| the number of "
                         "events strictly between; continuum reference: (tau/T)^d for two uniform points of the "
                         "flat d-dimensional diamond of duration T conditioned on being related; grid = SPECTRUM_GRID",
    "estimators": "events y uniform on the diamond; z uniform on fut(y); C_2 = N E fut(y); C_3 = N E[past(y) fut(y)]; "
                  "C_4 = N E[fut(y) past(y) fut(z)]; sum of intervals = N E[fut(y) |I(y,z)|] = C_3; spectrum weights fut(y); "
                  "standard errors from the sample variance of the per-draw terms",
    "seed_rule": "SAMPLE_SEED + 1000 dim + q + 10 region_index (numpy PCG64)",
}


# --------------------------------------------------------------------------
# Flat-diamond expectations
# --------------------------------------------------------------------------


def chi(k: int, d: float) -> float:
    """Meyer's chain coefficient: E C_k = chi_k (rho V)^k for a Poisson sprinkling of a flat d-diamond."""
    return (1.0 / k) * (Gamma(d + 1) / 2.0) ** (k - 1) * Gamma(d / 2.0) * Gamma(d) / (Gamma(k * d / 2.0) * Gamma((k + 1) * d / 2.0))


def invert_chi(k: int, value: float) -> float | None:
    """The spacetime dimension d in [1.05, 12] with chi_k(d) = value (chi_k decreases in d), or None."""
    if not (value > 0.0):
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


def continuum_spectrum(d: int, pairs: int, seed: int, grid) -> dict:
    """(tau/T)^d over related pairs of uniform points in the flat d-diamond of unit duration, by Monte Carlo."""
    rng = np.random.default_rng(seed)
    got, total, related = [], 0, 0
    while sum(len(g) for g in got) < pairs:
        B = 400_000
        u = rng.random((2, B))
        t = np.where(u < 0.5, (2.0 * u) ** (1.0 / d) / 2.0, 1.0 - (2.0 * (1.0 - u)) ** (1.0 / d) / 2.0)
        r = np.minimum(t, 1.0 - t)
        if d - 1 == 1:
            x = (rng.random((2, B, 1)) * 2.0 - 1.0) * r[..., None]
        else:
            g = rng.standard_normal((2, B, d - 1))
            g /= np.linalg.norm(g, axis=-1, keepdims=True)
            x = g * (r * rng.random((2, B)) ** (1.0 / (d - 1)))[..., None]
        dt = np.abs(t[1] - t[0])
        dx = np.linalg.norm(x[1] - x[0], axis=-1)
        rel = dx <= dt
        total += B
        related += int(rel.sum())
        got.append((dt[rel] ** 2 - dx[rel] ** 2) ** (d / 2.0))
    v = np.concatenate(got)[:pairs]
    qs = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    return {"spacetime_dimension": d, "pairs": int(pairs), "seed": seed,
            "ordering_fraction_monte_carlo": rounded(related / total), "ordering_fraction_exact": rounded(2.0 * chi(2, d)),
            "mean": rounded(float(v.mean())), "mean_exact_chi3_over_chi2": rounded(chi(3, d) / chi(2, d)),
            "second_moment": rounded(float((v * v).mean())),
            "quantiles": {str(p): rounded(float(np.quantile(v, p))) for p in qs},
            "cdf_on_grid": [rounded(float(np.mean(v <= g))) for g in grid]}


# --------------------------------------------------------------------------
# Family geometry (as in the causal-limit lane)
# --------------------------------------------------------------------------


class Family:
    def __init__(self, n: int, dim: int):
        self.n, self.dim = n, dim
        self.q, self.p = fibonacci(n)
        self.K = ceil_sqrt(self.q)
        q, dim = self.q, self.dim
        self.indptr, self.indices, self.values, self.A, self.B, self.sites = build_site_graph(q, dim)
        self.count = q ** dim
        degrees = np.diff(self.indptr)
        if int(degrees.min()) < 1:
            raise ValueError("a site lost its same-site read")
        self.degrees = degrees
        self.perm = [(b * self.p) % q for b in range(q)]
        self.ordered = sorted(range(q), key=lambda b: self.perm[b])
        values = self.values
        radii = [phi_sub((1, 0), values[self.ordered[-1]])]
        radii += [phi_scale(S.Fraction(1, 2), phi_sub(values[b], values[a])) for a, b in zip(self.ordered, self.ordered[1:])]
        radius = radii[0]
        for r in radii[1:]:
            if phi_sign(phi_sub(r, radius)) > 0:
                radius = r
        self.h_squared = phi_scale(dim, phi_square(radius))
        centre_axis = 0
        for b in range(1, q):
            if phi_sign(phi_sub(phi_square(phi_sub(values[b], (S.Fraction(1, 2), 0))),
                                phi_square(phi_sub(values[centre_axis], (S.Fraction(1, 2), 0))))) < 0:
                centre_axis = b
        self.centre_axis = centre_axis
        self.centre = self.site_of([centre_axis] * dim)
        self.a_q = 1.0 / sqrt(q)
        self.H_q = 2.0 * sqrt(3.0) / q if dim == 3 else 0.0
        self.density = q ** (dim + 0.5)
        self.xi = np.array([phi_float(v) for v in values])
        self.h_over_a = sqrt(q * phi_float(self.h_squared))
        self.pool = None
        self._dist_cache: dict = {}

    def site_of(self, labels) -> int:
        return int(sum(int(b) * self.q ** (self.dim - 1 - i) for i, b in enumerate(labels)))

    def labels_of(self, site: int) -> list[int]:
        return [int(v) for v in self.sites[site]]

    def distances_from(self, site: int) -> np.ndarray:
        if site not in self._dist_cache:
            d = bfs(self.indptr, self.indices, site, self.count)
            if int(d.min()) < 0:
                raise ValueError("disconnected site graph")
            self._dist_cache[site] = d
        return self._dist_cache[site]

    def cones(self, site: int, K: int) -> np.ndarray:
        At, Bt = metric_tables(self.q, self.A, self.B, self.sites, site)
        return cone_masks(self.q, At, Bt, K)

    def open_pool(self, processes: int, workdir: Path) -> dict:
        payload = {"K": self.K, "n": self.count}
        tag = f"{self.dim}_{self.q}"
        if processes > 1 and self.count >= POOL_MINIMUM_SITES:
            np.save(workdir / f"indptr_{tag}.npy", self.indptr)
            np.save(workdir / f"indices_{tag}.npy", self.indices)
            payload["indptr_path"] = str(workdir / f"indptr_{tag}.npy")
            payload["indices_path"] = str(workdir / f"indices_{tag}.npy")
        else:
            payload["indptr"], payload["indices"] = self.indptr, self.indices
        S._worker_init(payload)
        if "indptr_path" in payload:
            self.pool = mp.get_context("spawn").Pool(processes, initializer=S._worker_init, initargs=(payload,))
        return payload

    def close_pool(self) -> None:
        if self.pool is not None:
            self.pool.close()
            self.pool.join()
            self.pool = None

    def admissible_shift(self, K_prime: int) -> int:
        """Largest rank shift s such that every tip of the 14 directions keeps the K'-diamond's ball inside the cube."""
        need = K_prime * self.a_q / 2.0 + self.H_q
        rank_centre = self.perm[self.centre_axis]
        best = 0
        for s in range(1, self.q):
            lo, hi = rank_centre - s, rank_centre + s
            if lo < 0 or hi >= self.q:
                break
            ok = True
            for label in (self.ordered[lo], self.ordered[hi]):
                x = self.xi[label]
                if min(x, 1.0 - x) < need:
                    ok = False
            if ok:
                best = s
            else:
                break
        return best

    def tip_sites(self, shift: int) -> list[tuple[int, tuple[int, ...]]]:
        rank_centre = self.perm[self.centre_axis]
        dirs = []
        for axis in range(self.dim):
            for sign in (1, -1):
                dirs.append(tuple(sign if i == axis else 0 for i in range(self.dim)))
        if self.dim >= 2:
            for corner in range(2 ** self.dim):
                dirs.append(tuple(1 if (corner >> i) & 1 else -1 for i in range(self.dim)))
        out = []
        for u in dirs:
            labels = [self.ordered[rank_centre + shift * c] for c in u]
            out.append((self.site_of(labels), u))
        return out

    def moving_tips(self) -> dict | None:
        """The causal-limit lane's moving-tip rule: first rank shift with q ell^2 < K^2 and buffer >= H_q (then > 0)."""
        q, K, dim = self.q, self.K, self.dim
        rank_centre = self.perm[self.centre_axis]
        rules = [("buffer_at_least_H_q", self.H_q), ("buffer_positive", 0.0)] if dim == 3 else [("buffer_positive", 0.0)]
        for rule, floor_ in rules:
            for s in range(1, q):
                lo, hi = rank_centre - s, rank_centre + s
                if lo < 0 or hi >= q:
                    break
                xl, yl = self.ordered[lo], self.ordered[hi]
                ell2 = phi_scale(dim, (int(self.A[xl][yl]), int(self.B[xl][yl])))
                if phi_sign(phi_sub(phi_scale(q, ell2), (K * K, 0))) >= 0:
                    continue
                buffer = ellipsoid_buffer([self.xi[xl]] * dim, [self.xi[yl]] * dim, K * self.a_q, dim)
                if buffer > 0.0 and buffer >= floor_:
                    return {"rule": rule, "shift": s, "x_axis": xl, "y_axis": yl,
                            "x": self.site_of([xl] * dim), "y": self.site_of([yl] * dim),
                            "ell2": phi_float(ell2), "buffer": buffer}
        return None


# --------------------------------------------------------------------------
# Regions and the per-event searches (workers)
# --------------------------------------------------------------------------


_CTX: dict = {}


def _ctx(key: str) -> dict:
    """One cached region context per worker (A, G, cones, K, n)."""
    if _CTX.get("key") != key:
        with np.load(key, allow_pickle=False) as z:
            _CTX.clear()
            _CTX.update({k: z[k] for k in z.files})
            _CTX["key"] = key
    return _CTX


def save_region(path: Path, K: int, n: int, A, G, cone_x, cone_y) -> str:
    np.savez(path, K=np.int64(K), n=np.int64(n), A=np.asarray(A, dtype=np.int16), G=np.asarray(G, dtype=np.int16),
             cone_x=np.asarray(cone_x, dtype=bool), cone_y=np.asarray(cone_y, dtype=bool))
    return str(path)


def region_events(K: int, A: np.ndarray, G: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Number of layers per site and the event count of the region: (j, s) with A[s] <= j <= K - G[s]."""
    layers = np.maximum(K - G.astype(np.int64) - A.astype(np.int64) + 1, 0)
    return layers, layers.sum()


def _search(c: dict, s: int) -> np.ndarray:
    """Pruned breadth-first search from ``s`` exact on the past and the future of every event of the region at ``s``."""
    K, n = int(c["K"]), int(c["n"])
    A, G, cx, cy = c["A"], c["G"], c["cone_x"], c["cone_y"]
    a_s, g_s = int(A[s]), int(G[s])
    cap = K - min(a_s, g_s)
    allowed = []
    for m in range(cap):
        ix, iy = K - g_s - m, K - a_s - m
        mask = None
        if ix >= 0:
            mask = cx[ix].copy()
        if iy >= 0:
            mask = cy[iy].copy() if mask is None else (mask | cy[iy])
        if mask is None:
            mask = np.zeros(n, dtype=bool)
        allowed.append(mask)
    return bfs(S._G["indptr"], S._G["indices"], s, n, cap=cap, allowed=allowed)


def _cumulative(dist: np.ndarray, other: np.ndarray, cap: int, K: int) -> np.ndarray:
    """C[delta, g] = #{sites: 0 <= dist <= delta, other <= g}."""
    ok = (dist >= 0) & (other <= K)  # sites beyond K layers from a tip belong to no event of the region
    h = np.zeros((cap + 1, K + 1), dtype=np.int64)
    np.add.at(h, (dist[ok].astype(np.int64), other[ok].astype(np.int64)), 1)
    return h.cumsum(axis=0).cumsum(axis=1)


def _future_by_layer(CG: np.ndarray, j: int, K: int, cap: int) -> np.ndarray:
    """n_m = |{(m, u) in fut(j, s)}| for m = 0..K (zero for m <= j)."""
    out = np.zeros(K + 1, dtype=np.int64)
    for m in range(j + 1, K + 1):
        out[m] = CG[min(m - j, cap), K - m]
    return out


def _past_total(CA: np.ndarray, j: int, cap: int) -> int:
    total = 0
    for m in range(0, j):
        total += int(CA[min(j - m, cap), m])
    return total


def _interval(dist_s: np.ndarray, dist_u: np.ndarray, j: int, m: int) -> int:
    """|I((j,s),(m,u))| = sum over j < m' < m of #{w: d(s,w) <= m'-j, d(u,w) <= m-m'}."""
    if m - j < 2:
        return 0
    ok = (dist_s >= 0) & (dist_u >= 0)
    ds, du = dist_s[ok].astype(np.int64), dist_u[ok].astype(np.int64)
    span = m - j
    keep = ds + du <= span  # waiting reads count: d(s, w) = 0 or d(u, w) = 0 are allowed
    ds, du = ds[keep], du[keep]
    h = np.zeros((span + 1, span + 1), dtype=np.int64)
    np.add.at(h, (ds, du), 1)
    C = h.cumsum(axis=0).cumsum(axis=1)
    total = 0
    for mp_ in range(j + 1, m):
        total += int(C[mp_ - j, m - mp_])
    return total


def _sample_task(task) -> dict:
    key, events, seed = task
    c = _ctx(key)
    K, n = int(c["K"]), int(c["n"])
    A, G = c["A"], c["G"]
    rng = np.random.default_rng(seed)
    fut_y = np.zeros(len(events), dtype=np.int64)
    past_y = np.zeros(len(events), dtype=np.int64)
    fut_z = np.zeros(len(events), dtype=np.int64)
    past_z = np.zeros(len(events), dtype=np.int64)
    interval = np.zeros(len(events), dtype=np.int64)
    partners = np.zeros((len(events), 2), dtype=np.int64)
    for i, (j, s) in enumerate(events):
        j, s = int(j), int(s)
        dist_s = _search(c, s)
        cap_s = K - min(int(A[s]), int(G[s]))
        CG = _cumulative(dist_s, G, cap_s, K)
        CA = _cumulative(dist_s, A, cap_s, K)
        by_layer = _future_by_layer(CG, j, K, cap_s)
        fut_y[i] = int(by_layer.sum())
        past_y[i] = _past_total(CA, j, cap_s)
        if fut_y[i] == 0:
            partners[i] = (-1, -1)
            continue
        m = int(rng.choice(K + 1, p=by_layer / by_layer.sum()))
        cand = np.flatnonzero((dist_s >= 0) & (dist_s <= m - j) & (G <= K - m))
        u = int(cand[rng.integers(len(cand))])
        partners[i] = (m, u)
        dist_u = _search(c, u)
        cap_u = K - min(int(A[u]), int(G[u]))
        CGu = _cumulative(dist_u, G, cap_u, K)
        CAu = _cumulative(dist_u, A, cap_u, K)
        fut_z[i] = int(_future_by_layer(CGu, m, K, cap_u).sum())
        past_z[i] = _past_total(CAu, m, cap_u)
        interval[i] = _interval(dist_s, dist_u, j, m)
    return {"events": np.asarray(events, dtype=np.int64), "fut_y": fut_y, "past_y": past_y, "fut_z": fut_z,
            "past_z": past_z, "interval": interval, "partners": partners}


# --------------------------------------------------------------------------
# Exact enumeration (small q)
# --------------------------------------------------------------------------


def exact_region(fam: Family, K: int, A: np.ndarray, G: np.ndarray) -> dict:
    """Exact C_2, C_3, C_4 and interval spectrum from an all-pairs distance table (q^dim <= EXACT_MAX_SITES)."""
    n = fam.count
    if n > EXACT_MAX_SITES:
        raise ValueError("exact enumeration is limited to small populations")
    D = np.stack([bfs(fam.indptr, fam.indices, s, n, cap=K) for s in range(n)]).astype(np.int64)  # -1 beyond K
    D[D < 0] = K + 1
    A64, G64 = A.astype(np.int64), G.astype(np.int64)
    layers, N = region_events(K, A, G)
    events = [(j, s) for s in range(n) for j in range(int(A64[s]), K - int(G64[s]) + 1)]
    idx = {e: i for i, e in enumerate(events)}
    N = len(events)
    fut = np.zeros(N, dtype=np.int64)
    past = np.zeros(N, dtype=np.int64)
    fut_sets = []
    for i, (j, s) in enumerate(events):
        members = []
        for m in range(j + 1, K + 1):
            us = np.flatnonzero((D[s] <= m - j) & (G64 <= K - m))
            members.extend((m, int(u)) for u in us)
        fut_sets.append(members)
        fut[i] = len(members)
        p = 0
        for m in range(0, j):
            p += int(np.count_nonzero((D[s] <= j - m) & (A64 <= m)))
        past[i] = p
    C2 = int(fut.sum())
    C3 = int((past * fut).sum())
    C4 = 0
    intervals = []
    for i, (j, s) in enumerate(events):
        for (m, u) in fut_sets[i]:
            C4 += int(past[i]) * int(fut[idx[(m, u)]])
            intervals.append(_interval(D[s], D[u], j, m))
    intervals = np.asarray(intervals, dtype=np.int64)
    if int(intervals.sum()) != C3:
        raise ValueError("interval sum does not equal the three-chain count")
    return {"event_count": N, "C2": C2, "C3": C3, "C4": C4, "intervals": intervals}


# --------------------------------------------------------------------------
# Region readouts
# --------------------------------------------------------------------------


def draw_events(layers: np.ndarray, A: np.ndarray, size: int, rng) -> np.ndarray:
    """``size`` events uniform on the region: site with probability proportional to its layer count, then a layer."""
    p = layers / layers.sum()
    sites = rng.choice(len(layers), size=size, p=p)
    js = A[sites].astype(np.int64) + rng.integers(0, layers[sites])
    return np.stack([js, sites], axis=1)


def summarise_region(N: int, dim: int, out: dict, reference: dict, grid, exact: dict | None) -> dict:
    d = dim + 1
    fut_y, past_y, fut_z = out["fut_y"].astype(float), out["past_y"].astype(float), out["fut_z"].astype(float)
    interval = out["interval"].astype(float)
    M = len(fut_y)

    def est(terms: np.ndarray) -> dict:
        mean = float(terms.mean())
        se = float(terms.std(ddof=1) / sqrt(len(terms))) if len(terms) > 1 else 0.0
        return {"estimate": rounded(N * mean), "standard_error": rounded(N * se)}

    c2 = est(fut_y)
    c3 = est(past_y * fut_y)
    c3_from_intervals = est(fut_y * interval)
    c4 = est(fut_y * past_y * fut_z)
    chains = {}
    for k, e in ((2, c2), (3, c3), (4, c4)):
        ff = falling(N, k)
        ratio = e["estimate"] / ff if ff > 0 else None
        chains[str(k)] = {"count": e, "count_over_falling_factorial": None if ratio is None else rounded(ratio),
                          "flat_chi": rounded(chi(k, d)),
                          "relative_deviation_from_flat": None if ratio is None else rounded(ratio / chi(k, d) - 1.0),
                          "inverted_dimension": None if ratio is None else (None if invert_chi(k, ratio) is None else rounded(invert_chi(k, ratio))),
                          "relative_standard_error": rounded(e["standard_error"] / e["estimate"]) if e["estimate"] > 0 else None}
    ratios = {"C3_N_over_C2_squared": rounded(c3["estimate"] * N / c2["estimate"] ** 2) if c2["estimate"] > 0 else None,
              "flat_C3_N_over_C2_squared": rounded(chi(3, d) / chi(2, d) ** 2),
              "C4_N2_over_C2_cubed": rounded(c4["estimate"] * N * N / c2["estimate"] ** 3) if c2["estimate"] > 0 else None,
              "flat_C4_N2_over_C2_cubed": rounded(chi(4, d) / chi(2, d) ** 3)}
    # interval spectrum, weighted by fut(y) so that related pairs are uniform
    w = fut_y / fut_y.sum() if fut_y.sum() > 0 else np.full(M, 1.0 / M)
    v = interval / N
    order = np.argsort(v)
    vs, ws = v[order], w[order]
    cw = np.cumsum(ws)
    cdf = [rounded(float(cw[np.searchsorted(vs, g, side="right") - 1])) if np.any(vs <= g) else 0.0 for g in grid]
    mean_v = float((w * v).sum())
    spectrum = {"weighted_mean": rounded(mean_v), "flat_mean_chi3_over_chi2": rounded(chi(3, d) / chi(2, d)),
                "relative_deviation_of_mean_from_flat": rounded(mean_v / (chi(3, d) / chi(2, d)) - 1.0),
                "weighted_second_moment": rounded(float((w * v * v).sum())),
                "empty_interval_weight": rounded(float(w[interval == 0].sum())),
                "cdf_on_grid": cdf}
    dist = {}
    for key, ref in reference.items():
        ks = max(abs(a - b) for a, b in zip(cdf, ref["cdf_on_grid"]))
        dist[key] = {"max_cdf_difference_on_grid": rounded(ks),
                     "mean_ratio": rounded(mean_v / ref["mean"]) if ref["mean"] > 0 else None}
    closest = min(dist, key=lambda k_: dist[k_]["max_cdf_difference_on_grid"])
    spectrum["distance_to_flat_references"] = dist
    spectrum["closest_flat_spacetime_dimension"] = int(closest)
    row = {"event_count": N, "draws": M, "chains": chains, "chain_ratios": ratios, "interval_spectrum": spectrum,
           "sampled_events_sha256": digest(out["events"].tolist()), "partners_sha256": digest(out["partners"].tolist())}
    if exact is not None:
        row["exact"] = {"C2": exact["C2"], "C3": exact["C3"], "C4": exact["C4"],
                        "C3_over_falling_factorial": rounded(exact["C3"] / falling(N, 3)) if N >= 3 else None,
                        "C4_over_falling_factorial": rounded(exact["C4"] / falling(N, 4)) if N >= 4 else None,
                        "interval_mean_over_N": rounded(float(exact["intervals"].mean()) / N) if len(exact["intervals"]) else None,
                        "interval_cdf_on_grid": [rounded(float(np.mean(exact["intervals"] / N <= g))) for g in grid]}
    return row


def run_region(fam: Family, name: str, K: int, x: int, y: int, samples: int, seed: int, workdir: Path, tag: str,
               references: dict, grid, exact: bool) -> dict:
    A = fam.distances_from(x)
    G = A if y == x else fam.distances_from(y)
    layers, N = region_events(K, A, G)
    cone_x = fam.cones(x, K)
    cone_y = cone_x if y == x else fam.cones(y, K)
    key = save_region(workdir / f"region_{tag}.npz", K, fam.count, A, G, cone_x, cone_y)
    rng = np.random.default_rng(seed)
    events = draw_events(layers, A, samples, rng)
    tasks = [(key, events[i:i + BATCH_SIZE].tolist(), seed * 1000 + i) for i in range(0, len(events), BATCH_SIZE)]
    parts = list(fam.pool.imap(_sample_task, tasks) if fam.pool is not None else map(_sample_task, tasks))
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    ex = exact_region(fam, K, A, G) if exact else None
    row = summarise_region(int(N), fam.dim, out, references, grid, ex)
    row.update({"region": name, "layers": K, "lower_tip_site": x, "upper_tip_site": y, "seed": seed,
                "support_site_count": int(np.count_nonzero(layers))})
    return row


def build_family(n: int, dim: int, processes: int, workdir: Path, references: dict, scale: float, exact: bool) -> dict:
    t0 = time.time()
    fam = Family(n, dim)
    q, K = fam.q, fam.K
    refs = {str(d): references[str(d)] for d in CONTINUUM_DIMENSIONS}
    seed_base = SAMPLE_SEED + 1000 * dim + q
    fam.open_pool(processes, workdir)
    regions = []
    try:
        samples = max(64, int(EVENT_SAMPLES["centre"] * scale))
        regions.append(run_region(fam, "centre", K, fam.centre, fam.centre, samples, seed_base, workdir, f"c_{dim}_{q}", refs, SPECTRUM_GRID, exact))
        centre = regions[-1]
        centre["count_volume_coefficient"] = rounded(centre["event_count"] / (fam.density * diamond_volume(dim, K * fam.a_q)))
        S._log(f"level n={n} q={q} dim={dim}: centre N={centre['event_count']} draws={samples} in {time.time() - t0:.1f}s")
        moving = fam.moving_tips() if dim >= 2 else None
        if moving is not None:
            t1 = time.time()
            samples = max(64, int(EVENT_SAMPLES["moving"] * scale))
            row = run_region(fam, "moving", K, moving["x"], moving["y"], samples, seed_base + 10, workdir, f"m_{dim}_{q}", refs, SPECTRUM_GRID, exact)
            T = K * fam.a_q
            tau = sqrt(T * T - moving["ell2"])
            row.update({"tip_rule": moving["rule"], "rank_shift": moving["shift"], "tip_separation_over_T": rounded(sqrt(moving["ell2"]) / T),
                        "proper_duration_over_L": rounded(tau), "spatial_buffer_over_L": rounded(moving["buffer"]),
                        "count_volume_coefficient": rounded(row["event_count"] / (fam.density * diamond_volume(dim, tau)))})
            regions.append(row)
            S._log(f"  moving N={row['event_count']} in {time.time() - t1:.1f}s")
        t2 = time.time()
        tip_rows = []
        chosen = None
        minimum_shift = -(-q // 8)  # displace the tips by at least L/8 when a diamond of at least two layers allows it
        for want in (minimum_shift, 1):
            for K_prime in range(K, 1, -1):
                shift = fam.admissible_shift(K_prime)
                if shift >= want:
                    chosen = (K_prime, shift)
                    break
            if chosen is not None:
                break
        if chosen is not None and not exact:
            K_prime, shift = chosen
            samples = max(64, int(EVENT_SAMPLES["tip"] * scale))
            for i, (site, direction) in enumerate(fam.tip_sites(shift)):
                row = run_region(fam, "tip", K_prime, site, site, samples, seed_base + 20 + i, workdir, f"t_{dim}_{q}_{i}", refs, SPECTRUM_GRID, False)
                row.update({"direction": list(direction), "rank_shift": shift,
                            "count_volume_coefficient": rounded(row["event_count"] / (fam.density * diamond_volume(dim, K_prime * fam.a_q)))})
                tip_rows.append(row)
            centre_prime = run_region(fam, "centre_reference_for_tips", K_prime, fam.centre, fam.centre, samples, seed_base + 19, workdir,
                                      f"cr_{dim}_{q}", refs, SPECTRUM_GRID, False)
            centre_prime["count_volume_coefficient"] = rounded(centre_prime["event_count"] / (fam.density * diamond_volume(dim, K_prime * fam.a_q)))
            tip_rows.append(centre_prime)
            S._log(f"  {len(tip_rows) - 1} tips at {K_prime} layers, shift {shift}, in {time.time() - t2:.1f}s")
    finally:
        fam.close_pool()
        for f in workdir.glob(f"region_*_{dim}_{q}*.npz"):
            f.unlink()

    def spread(rows: list, path) -> dict | None:
        vals = [path(r) for r in rows if path(r) is not None]
        if not vals:
            return None
        return {"count": len(vals), "mean": rounded(float(np.mean(vals))), "sd": rounded(float(np.std(vals))),
                "min": rounded(float(min(vals))), "max": rounded(float(max(vals)))}

    homogeneity = None
    if tip_rows:
        homogeneity = {"layers": tip_rows[0]["layers"], "rank_shift": tip_rows[0].get("rank_shift"),
                       "tip_rule": "tips at the 14 axis and body-diagonal directions at the largest rank shift whose K'-diamond ball "
                                   "(radius K' a_q/2 + H_q) stays inside the cube; K' the largest layer count admitting a shift of at "
                                   "least ceil(q/8) ranks (about L/8), else the largest admitting any shift",
                       "diamonds": tip_rows,
                       "inverted_dimension_C3_spread": spread(tip_rows, lambda r: r["chains"]["3"]["inverted_dimension"]),
                       "inverted_dimension_C2_spread": spread(tip_rows, lambda r: r["chains"]["2"]["inverted_dimension"]),
                       "interval_mean_spread": spread(tip_rows, lambda r: r["interval_spectrum"]["weighted_mean"]),
                       "count_volume_coefficient_spread": spread(tip_rows, lambda r: r["count_volume_coefficient"])}
    return {"dimension": dim, "fibonacci_index": n, "q": q, "p": fam.p, "layer_steps": K, "site_count": fam.count,
            "centre_axis": fam.centre_axis, "centre_site": fam.centre, "edge_radius_a_over_L": rounded(fam.a_q),
            "model_time_T_over_L": rounded(K * fam.a_q), "h_over_a": rounded(fam.h_over_a),
            "event_density_times_L_to_dim_plus_one": rounded(fam.density),
            "minimum_neighbor_count_including_wait": int(fam.degrees.min()), "maximum_neighbor_count_including_wait": int(fam.degrees.max()),
            "regions": regions, "homogeneity": homogeneity, "wall_seconds": rounded(time.time() - t0)}


# --------------------------------------------------------------------------
# Receipt
# --------------------------------------------------------------------------


def references_table(pairs: int = CONTINUUM_PAIRS) -> dict:
    return {str(d): continuum_spectrum(d, pairs, CONTINUUM_SEED + d, SPECTRUM_GRID) for d in CONTINUUM_DIMENSIONS}


def trend(levels: list) -> dict:
    rows = []
    for lv in levels:
        for fam in lv["families"]:
            c = fam["regions"][0]
            mv = next((r for r in fam["regions"] if r["region"] == "moving"), None)
            row = {"dimension": fam["dimension"], "q": fam["q"], "layers": fam["layer_steps"], "centre_event_count": c["event_count"],
                   "centre_C2_over_falling": c["chains"]["2"]["count_over_falling_factorial"],
                   "centre_C3_over_falling": c["chains"]["3"]["count_over_falling_factorial"],
                   "centre_C4_over_falling": c["chains"]["4"]["count_over_falling_factorial"],
                   "flat_chi_2_3_4": [c["chains"][k]["flat_chi"] for k in ("2", "3", "4")],
                   "centre_inverted_dimension_C2_C3_C4": [c["chains"][k]["inverted_dimension"] for k in ("2", "3", "4")],
                   "centre_interval_mean_over_flat": c["interval_spectrum"]["relative_deviation_of_mean_from_flat"],
                   "centre_closest_flat_dimension": c["interval_spectrum"]["closest_flat_spacetime_dimension"],
                   "centre_max_cdf_difference_to_own_flat": c["interval_spectrum"]["distance_to_flat_references"][str(fam["dimension"] + 1)]["max_cdf_difference_on_grid"],
                   "moving_inverted_dimension_C2_C3_C4": None if mv is None else [mv["chains"][k]["inverted_dimension"] for k in ("2", "3", "4")],
                   "moving_interval_mean_over_flat": None if mv is None else mv["interval_spectrum"]["relative_deviation_of_mean_from_flat"],
                   "homogeneity_C3_dimension_spread": None if fam["homogeneity"] is None else fam["homogeneity"]["inverted_dimension_C3_spread"]}
            rows.append(row)
    return {"rows": rows, "statement": "per-level readouts of one nested construction; the finite runs do not demonstrate the asymptotic limit"}


def build(levels=LEVELS, dimensions=DIMENSIONS, processes: int | None = None, scale: float = 1.0, exact: bool = False,
          workdir: Path | None = None, reference_pairs: int = CONTINUUM_PAIRS) -> dict:
    processes = processes or max(1, min(6, os.cpu_count() or 1))
    t0 = time.time()
    references = references_table(reference_pairs)
    S._log(f"continuum references in {time.time() - t0:.1f}s")
    rows = []
    with tempfile.TemporaryDirectory(prefix="oph_manifold_sampled_", dir=None if workdir is None else str(workdir)) as tmp:
        for n in levels:
            families = [build_family(n, dim, processes, Path(tmp), references, scale, exact) for dim in dimensions]
            rows.append({"fibonacci_index": n, "q": fibonacci(n)[0], "families": families})
    pins = {p: file_sha256(ROOT / p) for p in PINS if (ROOT / p).exists()}
    return {"schema": SCHEMA, "scope": SCOPE, "claim_boundary": CLAIM_BOUNDARY, "conventions": CONVENTIONS,
            "sampling": {"event_samples": EVENT_SAMPLES, "scale": scale, "seed_base": SAMPLE_SEED, "batch_size": BATCH_SIZE,
                         "exact_enumeration": exact},
            "spectrum_grid": SPECTRUM_GRID, "continuum_references": references, "levels": rows,
            "trend": trend(rows), "source_pins": pins, "wall_seconds": rounded(time.time() - t0)}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--out", type=Path, default=OUTPUT)
    parser.add_argument("--processes", type=int, default=None)
    parser.add_argument("--levels", type=int, nargs="*", default=None)
    parser.add_argument("--dimensions", type=int, nargs="*", default=None)
    parser.add_argument("--scale", type=float, default=1.0, help="multiplier on the event sample sizes")
    parser.add_argument("--exact", action="store_true", help="also enumerate every event and partner (small q only)")
    parser.add_argument("--workdir", type=Path, default=None, help="directory for the memory-mapped graph (large q)")
    parser.add_argument("--reference-pairs", type=int, default=CONTINUUM_PAIRS)
    args = parser.parse_args(argv)
    t0 = time.time()
    data = canonical(build(levels=tuple(args.levels) if args.levels else LEVELS,
                           dimensions=tuple(args.dimensions) if args.dimensions else DIMENSIONS,
                           processes=args.processes, scale=args.scale, exact=args.exact, workdir=args.workdir,
                           reference_pairs=args.reference_pairs))
    S._log(f"build wall time {time.time() - t0:.1f}s")
    if args.write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_bytes(data)
    if args.check and args.out.read_bytes() != data:
        raise ValueError("manifold-sampled receipt is stale")
    print("MANIFOLD_SAMPLED_REPLAYED", len(data), hashlib.sha256(data).hexdigest())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Rebuild the expected OPH causal poset at one level from its definitions.

Standard library and numpy only; no simulator code is imported.  The script
restates the declared source-record family and the carrier realization of
the two mirrored receipts and compares every rebuilt quantity with the
receipt values.

The construction at level ``q`` (a Fibonacci number):

* golden orbit ``xi_b = b*phi - floor(b*phi)`` for ``0 <= b < q``, held as
  the integer pair ``(m, b)`` meaning ``m + b*phi`` with
  ``m = -floor(b*phi)`` and ``phi = (1 + sqrt5)/2``; all metric decisions are
  exact in ``Q(phi)`` through the sign rule of the simulator's
  ``source_net.py`` (for ``x = a + b*phi`` put ``c = 2a + b``; the sign of
  ``c^2 - 5 b^2`` and the signs of ``c`` and ``b`` decide);
* population ``[0, q)^3`` in product order, record
  ``z(b) = (b2 - m1, b2 + m1, b3 - m2, b3 + m2, b1 - m3, b1 + m3)``, its
  integer current section and primitive word length;
* one twelve-port icosahedral carrier per site holding ``z(b)`` as port
  loads on the six antipodal pairs, reading its position through the
  rank-three projector ``x = 2 P_slow N``; the readback metric equals the
  paper's source metric ``||s(b) - s(b')||`` with scale exactly one, which is
  checked in float on every pair and exactly in integers through
  ``10 S2 + 2 X sqrt5 = 20 A + (8 B - 4 A) sqrt5``;
* the complete-neighbour read law inside ``a_q = L/sqrt(q)`` with waiting,
  decided from the source metric (``q (A + B phi) <= 1``) and, independently,
  from the carrier readbacks (``sign((5 q S2 - 10) + (q X + 2) sqrt5) <= 0``);
  the two relations are compared edge for edge;
* rounds as layers, ``K_q = ceil(sqrt q)``, one event per site and round; the
  layered order ``(j, s) <= (j', t)`` iff ``d(s, t) <= j' - j`` with ``d`` the
  graph distance; the centre ``K``-diamond, its exact strict pair count, the
  ordering fraction, the inverted Myrheim-Meyer dimension and the count clock;
* the read law ``q_i(0) = i + 1``, ``q_i(j) = 1 + sum`` of the neighbours'
  previous values with the hash-chained audit trace, the ``+1`` intervention
  at the centre, and, when the stored event log is present, the provenance
  order generated from the log's read-after-write records alone.

Usage (from the package directory or anywhere):

    python3 build_causal_poset.py                 # q = 5 and q = 8, all checks
    python3 build_causal_poset.py --q 13          # one level
    python3 build_causal_poset.py --q 21 --no-pairs --no-reads
    python3 build_causal_poset.py --dim 2 --q 8   # a control family

Exit status 0 when every comparison agrees, 1 otherwise.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import time
from fractions import Fraction
from math import isqrt, lgamma, log, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE_NET_RECEIPT = HERE / "source_net_causal_limit_receipt.json"
CARRIER_RECEIPT = HERE / "carrier_source_net_receipt.json"
LOG_DIR = HERE / "carrier_source_net_logs"

FIBONACCI_INDEX = {3: 4, 5: 5, 8: 6, 13: 7, 21: 8, 34: 9, 55: 10}
REFERENCE_FRACTION = {1: Fraction(1, 2), 2: Fraction(8, 35), 3: Fraction(1, 10)}
PORT_COUNT = 12
POSITIVE_PORTS = (0, 1, 4, 5, 8, 9)
ID_BYTES = 4
VERSION_BYTES = 2
# The oriented faces of the committed twelve-port carrier (the simulator's
# ORIENTED_BASE_FACES); the thirty seams are their edges.
ORIENTED_FACES = (
    (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
    (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
    (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
    (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
)


# --------------------------------------------------------------------------
# Canonical JSON and digests
# --------------------------------------------------------------------------


def canonical(x) -> bytes:
    return (json.dumps(x, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def digest(x) -> str:
    return hashlib.sha256(canonical(x)).hexdigest()


def rounded(x: float) -> float:
    """Twelve significant digits, the receipt convention for derived floats."""
    return float(f"{float(x):.12g}")


# --------------------------------------------------------------------------
# Q(phi): (a, b) means a + b*phi
# --------------------------------------------------------------------------


def phi_sign(x) -> int:
    a, b = x
    c = 2 * a + b
    if b == 0:
        return (c > 0) - (c < 0)
    if c >= 0 and b > 0:
        return 1
    if c <= 0 and b < 0:
        return -1
    d = c * c - 5 * b * b
    return ((d > 0) - (d < 0)) * (1 if c > 0 else -1)


def phi_sign_array(a, b) -> np.ndarray:
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    c = 2 * a + b
    d = c * c - 5 * b * b
    mixed = np.sign(d) * np.where(c > 0, 1, -1)
    return np.where(b == 0, np.sign(c),
                    np.where((c >= 0) & (b > 0), 1,
                             np.where((c <= 0) & (b < 0), -1, mixed))).astype(np.int8)


def phi_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def phi_scale(c, x):
    return (c * x[0], c * x[1])


def phi_square(x):
    a, b = x
    return (a * a + b * b, 2 * a * b + b * b)


def phi_float(x) -> float:
    return float(x[0]) + float(x[1]) * (1.0 + sqrt(5.0)) / 2.0


def sqrt5_sign_array(a, b) -> np.ndarray:
    """Exact sign of a + b*sqrt5 (int64 arrays)."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    d = a * a - 5 * b * b
    mixed = np.sign(d) * np.where(a > 0, 1, -1)
    return np.where(b == 0, np.sign(a),
                    np.where((a >= 0) & (b > 0), 1,
                             np.where((a <= 0) & (b < 0), -1, mixed))).astype(np.int8)


def root_upper(x, denominator: int = 10 ** 6) -> Fraction:
    """Rational upper bound for sqrt(a + b*phi) on a fixed grid, by squaring."""
    if phi_sign(x) < 0:
        raise ValueError("negative radicand")
    if phi_sign(x) == 0:
        return Fraction(0)
    lo, hi = 0, denominator
    while phi_sign(phi_sub((Fraction(hi, denominator) ** 2, 0), x)) < 0:
        hi *= 2
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if phi_sign(phi_sub((Fraction(mid, denominator) ** 2, 0), x)) >= 0:
            hi = mid
        else:
            lo = mid
    return Fraction(hi, denominator)


# --------------------------------------------------------------------------
# Golden orbit, population, records
# --------------------------------------------------------------------------


def golden_floor(b: int) -> int:
    """floor(b*phi) exactly: (b + isqrt(5 b^2)) // 2."""
    return (b + isqrt(5 * b * b)) // 2


def ceil_sqrt(q: int) -> int:
    return isqrt(q) + (isqrt(q) ** 2 < q)


def orbit(q: int) -> list[tuple[int, int]]:
    return [(-golden_floor(b), b) for b in range(q)]


def axis_tables(values):
    """(xi_i - xi_j)^2 = A + B*phi in units of L^2, as two int64 tables."""
    m = np.array([v[0] for v in values], dtype=np.int64)
    b = np.array([v[1] for v in values], dtype=np.int64)
    da = m[:, None] - m[None, :]
    db = b[:, None] - b[None, :]
    return da * da + db * db, 2 * da * db + db * db


def site_coordinates(q: int, dim: int) -> np.ndarray:
    return np.indices((q,) * dim).reshape(dim, -1).T.copy()


def source_records(q: int, sites: np.ndarray):
    """Records z(b), current sections, word lengths, and RER's record list."""
    m = np.array([-golden_floor(b) for b in range(q)], dtype=np.int64)
    b = sites.astype(np.int64)
    a = m[b]
    z = np.stack([b[:, 1] - a[:, 0], b[:, 1] + a[:, 0], b[:, 2] - a[:, 1],
                  b[:, 2] + a[:, 1], b[:, 0] - a[:, 2], b[:, 0] + a[:, 2]], axis=1)
    total = z.sum(axis=1)
    if np.any(total % 2):
        raise ValueError("a source record left the even-sum sublattice")
    h = total // 2
    currents = np.stack([-z[:, 5], h - z[:, 0], h - z[:, 2] - z[:, 3],
                         h - z[:, 1] - z[:, 4] - z[:, 5], h - z[:, 2] - z[:, 3] - z[:, 4],
                         z[:, 3]], axis=1)
    lengths = np.abs(currents).sum(axis=1)
    records = [[bb, zz, cc] for bb, zz, cc in zip(b.tolist(), z.tolist(), currents.tolist())]
    return z, currents, lengths, records


def centre_axis_label(values) -> int:
    """The orbit label nearest to 1/2 (smallest label on ties), exactly."""
    half = (Fraction(1, 2), 0)
    centre = 0
    for b in range(1, len(values)):
        if phi_sign(phi_sub(phi_square(phi_sub(values[b], half)),
                            phi_square(phi_sub(values[centre], half)))) < 0:
            centre = b
    return centre


# --------------------------------------------------------------------------
# The twelve-port carrier: seams, slow band, antipode, readback
# --------------------------------------------------------------------------


def carrier_tables() -> dict:
    seams = sorted({tuple(sorted((f[i], f[(i + 1) % 3]))) for f in ORIENTED_FACES for i in range(3)})
    if len(seams) != 30:
        raise ValueError("the oriented faces do not define thirty seams")
    adjacency = np.zeros((PORT_COUNT, PORT_COUNT), dtype=np.int64)
    for i, j in seams:
        adjacency[i, j] = adjacency[j, i] = 1
    if not np.all(adjacency.sum(axis=1) == 5):
        raise ValueError("the carrier is not five-regular")
    laplacian = 5 * np.eye(PORT_COUNT, dtype=np.int64) - adjacency
    # Graph distances by repeated squaring of the reachability.
    dist = np.full((PORT_COUNT, PORT_COUNT), -1, dtype=np.int64)
    np.fill_diagonal(dist, 0)
    reach = np.eye(PORT_COUNT, dtype=bool)
    for step in range(1, PORT_COUNT):
        nxt = reach | ((reach.astype(np.int64) @ adjacency) > 0)
        dist[nxt & ~reach] = step
        reach = nxt
        if reach.all():
            break
    antipode = []
    for p in range(PORT_COUNT):
        far = [r for r in range(PORT_COUNT) if dist[p, r] == 3]
        if len(far) != 1:
            raise ValueError("the distance-three partner is not unique")
        antipode.append(far[0])
    if any(antipode[antipode[p]] != p or antipode[p] == p for p in range(PORT_COUNT)):
        raise ValueError("the antipode is not a fixed-point-free involution")
    # The slow band: the eigenspace of the seam Laplacian at 5 - sqrt5.
    w, v = np.linalg.eigh(laplacian.astype(float))
    slow = np.abs(w - (5.0 - sqrt(5.0))) < 1e-9
    if int(slow.sum()) != 3:
        raise ValueError("the slow band is not three-dimensional")
    p_slow = v[:, slow] @ v[:, slow].T
    gram = 4.0 * p_slow
    expected = {0: 1.0, 1: 1.0 / sqrt(5.0), 2: -1.0 / sqrt(5.0), 3: -1.0}
    for p in range(PORT_COUNT):
        for r in range(PORT_COUNT):
            if abs(gram[p, r] - expected[int(dist[p, r])]) > 1e-12:
                raise ValueError("the Gram entries by port distance are not {1, 1/sqrt5, -1/sqrt5, -1}")
    sigma = np.zeros((6, 6), dtype=np.int64)
    for k in range(6):
        for l in range(6):
            if k != l:
                d = int(dist[POSITIVE_PORTS[k], POSITIVE_PORTS[l]])
                if d not in (1, 2):
                    raise ValueError("positive ports are not pairwise non-antipodal")
                sigma[k, l] = 1 if d == 1 else -1
    return {"seams": seams, "laplacian": laplacian, "distance": dist, "antipode": antipode,
            "p_slow": p_slow, "gram": gram, "sigma": sigma}


def place_loads(z: np.ndarray, antipode) -> np.ndarray:
    """N_{p_k} = max(z_k, 0), N_{-p_k} = max(-z_k, 0) on the six antipodal pairs."""
    loads = np.zeros((len(z), PORT_COUNT), dtype=np.int64)
    for k, p in enumerate(POSITIVE_PORTS):
        loads[:, p] = np.maximum(z[:, k], 0)
        loads[:, antipode[p]] = np.maximum(-z[:, k], 0)
    return loads


def readback(loads: np.ndarray, p_slow: np.ndarray) -> np.ndarray:
    """x = 2 P_slow N for every carrier (rows)."""
    return loads.astype(float) @ (2.0 * p_slow).T


def pair_blocks(n: int, rows: int):
    for lo in range(0, n, rows):
        hi = min(lo + rows, n)
        i = np.repeat(np.arange(lo, hi, dtype=np.int64), n)
        j = np.tile(np.arange(n, dtype=np.int64), hi - lo)
        keep = j > i
        yield i[keep], j[keep]


def metric_identity(q: int, sites: np.ndarray, z: np.ndarray, x: np.ndarray, sigma: np.ndarray,
                    A: np.ndarray, B: np.ndarray) -> dict:
    """Readback metric against the source metric: float on every pair, exact in integers."""
    n = len(z)
    l2 = 2.0 - 2.0 / sqrt(5.0)
    rows = max(1, 400_000 // n)
    float_max = 0.0
    exact = True
    pairs = 0
    for i, j in pair_blocks(n, rows):
        dz = z[i] - z[j]
        s2 = (dz * dz).sum(axis=1)
        xq = np.einsum("mk,kl,ml->m", dz, sigma, dz)
        At = np.zeros(len(i), dtype=np.int64)
        Bt = np.zeros(len(i), dtype=np.int64)
        for axis in range(3):
            At += A[sites[i, axis], sites[j, axis]]
            Bt += B[sites[i, axis], sites[j, axis]]
        # 10 ||dx||^2 = 10 S2 + 2 X sqrt5 and 10 L^2 |dxi|^2 = 20 A + (8 B - 4 A) sqrt5.
        exact = exact and bool(np.array_equal(10 * s2, 20 * At) and np.array_equal(2 * xq, 8 * Bt - 4 * At))
        dx = x[i] - x[j]
        lhs = (dx * dx).sum(axis=1)
        rhs = l2 * (At + Bt * (1.0 + sqrt(5.0)) / 2.0)
        float_max = max(float_max, float(np.max(np.abs(lhs - rhs))))
        pairs += len(i)
    return {"pairs": pairs, "exact_integer_identity": exact, "float_max_deviation": float_max,
            "float_identity": float_max < 1e-9, "scale_to_paper_position": 1}


# --------------------------------------------------------------------------
# Neighbour graph: the read law
# --------------------------------------------------------------------------


def site_graph(q: int, dim: int, values, A, B, sites: np.ndarray):
    """CSR neighbour lists (same-site read included, rows ascending), exact decisions."""
    within = phi_sign_array(q * A - 1, q * B) <= 0
    cand = [np.flatnonzero(within[x]) for x in range(q)]
    n = q ** dim
    powers = [q ** (dim - 1 - i) for i in range(dim)]
    rows_out = []
    for s in range(n):
        coords = sites[s]
        lists = [cand[int(coords[axis])] for axis in range(dim)]
        grids = np.meshgrid(*lists, indexing="ij")
        Asum = np.zeros(grids[0].shape, dtype=np.int64)
        Bsum = np.zeros(grids[0].shape, dtype=np.int64)
        nb = np.zeros(grids[0].shape, dtype=np.int64)
        for axis in range(dim):
            Asum += A[int(coords[axis])][grids[axis]]
            Bsum += B[int(coords[axis])][grids[axis]]
            nb += grids[axis] * powers[axis]
        ok = phi_sign_array(q * Asum - 1, q * Bsum) <= 0
        rows_out.append(np.sort(nb[ok]))
    indptr = np.zeros(n + 1, dtype=np.int64)
    indptr[1:] = np.cumsum([len(r) for r in rows_out])
    indices = np.concatenate(rows_out).astype(np.int32)
    return indptr, indices


def carrier_neighbour_relation(q: int, z: np.ndarray, sigma: np.ndarray) -> set:
    """All pairs decided from the readbacks alone: sign((5 q S2 - 10) + (q X + 2) sqrt5) <= 0."""
    n = len(z)
    rows = max(1, 400_000 // n)
    edges = []
    for i, j in pair_blocks(n, rows):
        dz = z[i] - z[j]
        s2 = (dz * dz).sum(axis=1)
        xq = np.einsum("mk,kl,ml->m", dz, sigma, dz)
        ok = sqrt5_sign_array(5 * q * s2 - 10, q * xq + 2) <= 0
        edges.append(np.stack([i[ok], j[ok]], axis=1))
    e = np.concatenate(edges) if edges else np.zeros((0, 2), dtype=np.int64)
    return set(map(tuple, e.tolist()))


def csr_edges(indptr, indices) -> set:
    out = set()
    for i in range(len(indptr) - 1):
        for j in indices[indptr[i]:indptr[i + 1]].tolist():
            if j > i:
                out.add((i, j))
    return out


def neighbour_digest(indptr, indices) -> str:
    n = len(indptr) - 1
    h = hashlib.sha256()
    h.update(b"[")
    for i in range(n):
        if i:
            h.update(b",")
        h.update(("[" + ",".join(map(str, indices[indptr[i]:indptr[i + 1]].tolist())) + "]").encode("ascii"))
    h.update(b"]\n")
    return h.hexdigest()


# --------------------------------------------------------------------------
# Distances, cones, the layered order
# --------------------------------------------------------------------------


def bfs(indptr, indices, start: int, n: int, cap: int | None = None) -> np.ndarray:
    dist = np.full(n, -1, dtype=np.int64)
    dist[start] = 0
    frontier = np.array([start], dtype=np.int64)
    layers = n if cap is None else cap
    for m in range(layers):
        if frontier.size == 0:
            break
        pieces = [indices[indptr[u]:indptr[u + 1]] for u in frontier.tolist()]
        nxt = np.unique(np.concatenate(pieces)).astype(np.int64)
        nxt = nxt[dist[nxt] < 0]
        if nxt.size == 0:
            break
        dist[nxt] = m + 1
        frontier = nxt
    return dist


def metric_tables(A, B, sites: np.ndarray, site: int):
    dim = sites.shape[1]
    At = np.zeros(len(sites), dtype=np.int64)
    Bt = np.zeros(len(sites), dtype=np.int64)
    for axis in range(dim):
        At += A[sites[site, axis]][sites[:, axis]]
        Bt += B[sites[site, axis]][sites[:, axis]]
    return At, Bt


def myrheim_meyer_fraction(d: float) -> float:
    return float(np.exp(lgamma(d + 1.0) + lgamma(d / 2.0) - log(2.0) - lgamma(3.0 * d / 2.0)))


def invert_myrheim_meyer(f: float):
    if not 0.0 < f < 1.0:
        return None
    lo, hi = 1.01, 20.0
    if not (myrheim_meyer_fraction(hi) <= f <= myrheim_meyer_fraction(lo)):
        return None
    for _ in range(96):
        mid = (lo + hi) / 2.0
        if myrheim_meyer_fraction(mid) > f:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def diamond_counts(alpha: np.ndarray, k: int) -> list[int]:
    return [int(np.count_nonzero(alpha <= min(j, k - j))) for j in range(k + 1)]


def diamond_event_digest(alpha: np.ndarray, k: int) -> str:
    events = []
    for j in range(k + 1):
        events.extend([j, s] for s in np.flatnonzero(alpha <= min(j, k - j)).tolist())
    return digest(events)


def strict_pair_counts(indptr, indices, alpha: np.ndarray, K: int, n: int) -> dict[int, int]:
    """Exact strict ordered-pair counts of the centre k-diamonds, k = 1..K.

    For sites ``s, t`` at graph distances ``alpha_s, alpha_t`` from the centre
    and ``d = d(s, t)``, the number of ordered event pairs
    ``(j, s) < (j', t)`` inside the k-diamond is the number of layer pairs
    with ``alpha_s <= j <= k - alpha_s``, ``alpha_t <= j' <= k - alpha_t`` and
    ``j' - j >= max(d, 1)``.  The search from ``s`` stops at distance
    ``K - alpha_s``: an ordered pair needs ``d <= j' - j <= K - alpha_s``, so
    larger distances cannot form one inside any diamond.
    """
    support = np.flatnonzero(alpha <= K // 2)
    bins = K + 2
    hist = np.zeros((K + 1, K + 1, bins), dtype=np.int64)
    a_sup = alpha[support]
    for s in support.tolist():
        a_s = int(alpha[s])
        dist = bfs(indptr, indices, s, n, cap=K - a_s)
        d = dist[support]
        d = np.where(d < 0, K + 1, d)
        hist[a_s] += np.bincount(a_sup * bins + d, minlength=(K + 1) * bins).reshape(K + 1, bins)
    counts = {}
    for k in range(1, K + 1):
        C = 0
        for a_s in range(k // 2 + 1):
            for a_t in range(k // 2 + 1):
                for d in range(k + 1):
                    h = int(hist[a_s, a_t, d])
                    if h == 0:
                        continue
                    pairs = 0
                    for j in range(a_s, k - a_s + 1):
                        for jp in range(a_t, k - a_t + 1):
                            if jp - j >= max(d, 1):
                                pairs += 1
                    C += h * pairs
        counts[k] = C
    return counts


# --------------------------------------------------------------------------
# The read law with the hash-chained audit trace
# --------------------------------------------------------------------------


def run_reads(indptr, indices, rounds: int, intervention: int | None = None) -> dict:
    """q_i(0) = i + 1 (+1 at the intervention site); q_i(j) = 1 + sum of version-j reads."""
    n = len(indptr) - 1
    rows = [indices[indptr[i]:indptr[i + 1]].tolist() for i in range(n)]
    chain = bytes(32)
    layer_hashes, layer_sums, values = [], [], []
    prev = None
    for j in range(rounds + 1):
        if j == 0:
            cur = [1 + i + (1 if i == intervention else 0) for i in range(n)]
            for i in range(n):
                chain = hashlib.sha256(chain + canonical([[0, i], [], [i, 1, [0, i], cur[i]]])).digest()
        else:
            cur = [1 + sum(prev[r] for r in row) for row in rows]
            for i in range(n):
                reads = [[r, j, [j - 1, r], prev[r]] for r in rows[i]]
                chain = hashlib.sha256(chain + canonical([[j, i], reads, [i, j + 1, [j, i], cur[i]]])).digest()
        layer_hashes.append(digest(cur))
        layer_sums.append(sum(cur))
        values.append(cur)
        prev = cur
    return {"audit_trace_sha256": chain.hex(), "layer_value_sha256": layer_hashes,
            "layer_value_sums": layer_sums, "event_count": (rounds + 1) * n,
            "authenticated_read_count": rounds * int(indptr[-1]), "values": values}


def byte_length(v: int) -> int:
    return max(1, (int(v).bit_length() + 7) // 8)


def operation_costs(indptr, values: list[list[int]]) -> dict:
    n = len(indptr) - 1
    degrees = np.diff(indptr).astype(np.int64)
    total_reads = total_read_bytes = total_write_bytes = 0
    per_round = []
    for j, cur in enumerate(values):
        write_bytes = sum(ID_BYTES + VERSION_BYTES + byte_length(v) for v in cur)
        if j == 0:
            reads, read_bytes = 0, 0
        else:
            prev_len = np.array([byte_length(v) for v in values[j - 1]], dtype=np.int64)
            reads = int(degrees.sum())
            read_bytes = int((degrees * (2 * ID_BYTES + VERSION_BYTES + prev_len)).sum())
        per_round.append({"round": j, "reads": reads, "writes": n, "read_bytes": read_bytes,
                          "write_bytes": int(write_bytes), "total_bytes": read_bytes + int(write_bytes),
                          "maximum_value_bytes": max(byte_length(v) for v in cur)})
        total_reads += reads
        total_read_bytes += read_bytes
        total_write_bytes += int(write_bytes)
    return {"per_round": per_round, "total_reads": total_reads, "total_writes": n * len(values),
            "total_read_bytes": total_read_bytes, "total_write_bytes": total_write_bytes,
            "total_bytes": total_read_bytes + total_write_bytes}


# --------------------------------------------------------------------------
# Provenance order from the stored event log alone
# --------------------------------------------------------------------------


def load_log(q: int) -> dict | None:
    path = LOG_DIR / f"q{q}_event_log.json.gz"
    if not path.is_file():
        return None
    data = gzip.open(path).read()
    return {"payload": json.loads(data), "uncompressed_sha256": hashlib.sha256(data).hexdigest(),
            "uncompressed_bytes": len(data)}


def provenance_from_log(log: dict, centre: int) -> dict:
    """Generate the order from the log by the read-after-write rule; no layer labels are used."""
    n = int(log["site_count"])
    rounds = int(log["rounds"])
    events = log["events"]
    writer_of = {}
    for event_id, ev in enumerate(events):
        register, version, _value = ev["write"]
        if (register, version) in writer_of:
            raise ValueError("a register version has two writers")
        writer_of[(register, version)] = event_id
    parents = []
    rank = [0] * len(events)
    chronological = True
    rank_equals_round = True
    edge_count = 0
    read_relation = {}
    relation_identical = True
    for event_id, ev in enumerate(events):
        ps = []
        for register, version in ev["reads"]:
            if (register, version) not in writer_of:
                raise ValueError("a read names an unwritten register version")
            p = writer_of[(register, version)]
            if p == event_id:
                raise ValueError("an event reads its own committed version")
            chronological = chronological and p < event_id
            ps.append(p)
        parents.append(ps)
        edge_count += len(ps)
        if ps:
            rank[event_id] = 1 + max(rank[p] for p in ps)
        j, i = ev["event"]
        rank_equals_round = rank_equals_round and rank[event_id] == j
        if j >= 1:
            regs = [r for r, _v in ev["reads"]]
            if j == 1:
                read_relation[i] = regs
            else:
                relation_identical = relation_identical and read_relation.get(i) == regs
    # Cones from the parent relation: future of the centre's seed, past of its last event.
    children = [[] for _ in events]
    for e, ps in enumerate(parents):
        for p in ps:
            children[p].append(e)
    index = {tuple(ev["event"]): e for e, ev in enumerate(events)}
    seed = index[(0, centre)]
    last = index[(rounds, centre)]
    future = {seed}
    stack = [seed]
    while stack:
        e = stack.pop()
        for c in children[e]:
            if c not in future:
                future.add(c)
                stack.append(c)
    past = {last}
    stack = [last]
    while stack:
        e = stack.pop()
        for p in parents[e]:
            if p not in past:
                past.add(p)
                stack.append(p)
    interval = sorted(tuple(events[e]["event"]) for e in future & past)
    counts = [sum(1 for j, _s in interval if j == jj) for jj in range(rounds + 1)]
    future_counts = [sum(1 for e in future if events[e]["event"][0] == jj) for jj in range(rounds + 1)]
    future_digests = [digest(sorted(events[e]["event"][1] for e in future if events[e]["event"][0] == jj))
                      for jj in range(rounds + 1)]
    rows = [read_relation[i] for i in range(n)]
    rel_indptr = np.zeros(n + 1, dtype=np.int64)
    rel_indptr[1:] = np.cumsum([len(r) for r in rows])
    rel_indices = np.concatenate([np.asarray(r, dtype=np.int32) for r in rows])
    return {"edge_count": edge_count, "single_writer": True, "all_reads_resolved": True,
            "parents_precede_children": chronological, "derived_rank_equals_round": rank_equals_round,
            "read_relation_identical_across_rounds": relation_identical,
            "read_relation_sha256": neighbour_digest(rel_indptr, rel_indices),
            "interval_event_set_sha256": digest([list(e) for e in interval]),
            "interval_counts_by_round": counts, "interval_event_count": len(interval),
            "future_cone_counts": future_counts, "future_cone_sha256": future_digests}


# --------------------------------------------------------------------------
# Comparison helpers
# --------------------------------------------------------------------------


class Report:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def check(self, name: str, ours, theirs) -> bool:
        ok = canonical(ours) == canonical(theirs)
        self.rows.append({"check": name, "agree": ok, "rebuilt": ours, "receipt": theirs})
        return ok

    def note(self, name: str, value) -> None:
        self.rows.append({"check": name, "agree": True, "rebuilt": value, "receipt": None})

    @property
    def ok(self) -> bool:
        return all(r["agree"] for r in self.rows)


def load_receipts() -> tuple[dict, dict]:
    return (json.loads(SOURCE_NET_RECEIPT.read_bytes()), json.loads(CARRIER_RECEIPT.read_bytes()))


def family_row(receipt: dict, q: int, dim: int) -> dict | None:
    for level in receipt["levels"]:
        if level["q"] == q:
            for fam in level["families"]:
                if fam["dimension"] == dim:
                    return fam
    return None


def carrier_row(receipt: dict, q: int) -> dict | None:
    for level in receipt["levels"]:
        if level["q"] == q:
            return level
    return None


# --------------------------------------------------------------------------
# One level
# --------------------------------------------------------------------------


def build_level(q: int, dim: int = 3, pairs: bool | None = None, reads: bool | None = None,
                log_out=sys.stderr) -> Report:
    t0 = time.time()
    if q not in FIBONACCI_INDEX:
        raise ValueError("q must be a Fibonacci number of the receipt family")
    if pairs is None:
        pairs = q <= 13 or dim < 3
    if reads is None:
        reads = q <= 13 and dim == 3
    source_receipt, carrier_receipt = load_receipts()
    fam = family_row(source_receipt, q, dim)
    crow = carrier_row(carrier_receipt, q) if dim == 3 else None
    R = Report()
    R.note("level", {"q": q, "dimension": dim, "fibonacci_index": FIBONACCI_INDEX[q]})
    if fam is None:
        R.check("family present in the source-net receipt", True, False)
        return R
    K = ceil_sqrt(q)
    n = q ** dim
    values = orbit(q)
    A, B = axis_tables(values)
    sites = site_coordinates(q, dim)
    R.check("layer_steps K = ceil(sqrt q)", K, fam["layer_steps"])
    R.check("site_count q^dim", n, fam["site_count"])
    R.check("orbit (m, b) pairs", [[str(v[0]), str(v[1])] for v in values], fam["orbit_Qphi"])
    caxis = centre_axis_label(values)
    centre = sum(caxis * q ** (dim - 1 - i) for i in range(dim))
    R.check("centre axis label", caxis, fam["centre_axis"])
    R.check("centre site", centre, fam["intervention_source_id"])

    # Records and the carrier realization (three-dimensional family only).
    carrier = None
    if dim == 3:
        z, currents, lengths, records = source_records(q, sites)
        R.check("source records digest", digest(records), fam["source_records_sha256"])
        R.check("maximum word length", int(lengths.max()), fam["maximum_word_length"])
        R.check("sum of word lengths", int(lengths.sum()), fam["sum_word_lengths"])
        R.check("word length bound 27(q-1)", 27 * (q - 1), fam["word_length_bound"])
        R.check("word lengths within bound", bool(lengths.max() <= 27 * (q - 1)), True)
        carrier = carrier_tables()
        loads = place_loads(z, carrier["antipode"])
        quotient = np.stack([loads[:, p] - loads[:, carrier["antipode"][p]] for p in POSITIVE_PORTS], axis=1)
        R.check("load placement round trip (antipodal-odd quotient returns z)", bool(np.array_equal(quotient, z)), True)
        R.check("loads nonnegative with total |z|_1", bool(loads.min() >= 0 and loads.sum() == np.abs(z).sum()), True)
        x = readback(loads, carrier["p_slow"])
        gens = np.stack([2.0 * carrier["p_slow"][:, p] for p in POSITIVE_PORTS], axis=0)
        R.check("readback equals the signed generator sum sum_k z_k v_{p_k}",
                bool(np.max(np.abs(x - z.astype(float) @ gens)) < 1e-9), True)
        ident = metric_identity(q, sites, z, x, carrier["sigma"], A, B)
        R.check("readback metric equals the source metric, scale one (float, all pairs)", ident["float_identity"], True)
        R.check("readback metric identity exact in integers (all pairs)", ident["exact_integer_identity"], True)
        R.note("metric identity pairs checked", ident["pairs"])
        if crow is not None:
            R.check("carrier receipt: antipode", carrier["antipode"], carrier_receipt["carrier"]["antipode"])
            R.check("carrier receipt: records digest", digest(records), crow["records"]["source_records_sha256"])
            R.check("carrier receipt: load total", int(loads.sum()), crow["load_placement"]["total_load"])
            R.check("carrier receipt: maximum port load", int(loads.max()), crow["load_placement"]["maximum_port_load"])
            R.check("carrier receipt: metric scale to the paper position", "1",
                    carrier_receipt["readback_metric"]["scale_to_paper_position_s"])
    _log(log_out, f"q={q} dim={dim}: records and readbacks {time.time() - t0:.1f}s")

    # Neighbour graph from the source metric, and from the carrier readbacks.
    t1 = time.time()
    indptr, indices = site_graph(q, dim, values, A, B, sites)
    degrees = np.diff(indptr)
    nb_sha = neighbour_digest(indptr, indices)
    R.check("neighbour digest (waiting included)", nb_sha, fam["neighbors_including_wait_sha256"])
    R.check("undirected spatial edges", int((int(degrees.sum()) - n) // 2), fam["undirected_spatial_edges"])
    R.check("minimum neighbour count", int(degrees.min()), fam["minimum_neighbor_count_including_wait"])
    R.check("maximum neighbour count", int(degrees.max()), fam["maximum_neighbor_count_including_wait"])
    if dim == 3:
        carrier_edges = carrier_neighbour_relation(q, z, carrier["sigma"])
        R.check("carrier-decided relation equals the source-decided relation (edge sets)",
                carrier_edges == csr_edges(indptr, indices), True)
        if crow is not None:
            R.check("carrier receipt: neighbour digest", nb_sha, crow["neighbours"]["neighbors_including_wait_sha256"])
            R.check("carrier receipt: accepted pairs", len(carrier_edges), crow["neighbours"]["accepted_pairs"])
    _log(log_out, f"q={q} dim={dim}: neighbour graph {time.time() - t1:.1f}s (edges {(int(degrees.sum()) - n) // 2})")

    # Cone sandwich probes and the covering data.
    p_next = fibonacci_next(q)
    perm = [(b * p_next) % q for b in range(q)]
    ordered = sorted(range(q), key=lambda b: perm[b])
    radii = [phi_sub((1, 0), values[ordered[-1]])]
    radii += [phi_scale(Fraction(1, 2), phi_sub(values[b], values[a])) for a, b in zip(ordered, ordered[1:])]
    radius = radii[0]
    for r in radii[1:]:
        if phi_sign(phi_sub(r, radius)) > 0:
            radius = r
    h_squared = phi_scale(dim, phi_square(radius))
    ratio_upper = root_upper(phi_scale(q, h_squared))
    inner_speed = max(Fraction(0), 1 - 2 * ratio_upper)
    R.check("grid permutation b p mod q", perm, fam["grid_permutation"])
    R.check("one-dimensional fill radius r_q", [str(Fraction(radius[0])), str(Fraction(radius[1]))],
            fam["one_dimensional_fill_over_L_Qphi"])
    R.check("h_q^2 = dim r_q^2 (units L^2)", [str(Fraction(h_squared[0])), str(Fraction(h_squared[1]))],
            fam["whole_cube_fill_h_squared_over_L2_Qphi"])
    R.check("h/a upper bound", str(ratio_upper), fam["h_over_a_upper"])
    R.check("certified inner speed lower bound", str(inner_speed), fam["certified_inner_speed_lower"])
    R.check("positive inner cone 4 q h^2 < 1", phi_sign(phi_sub((1, 0), phi_scale(4 * q, h_squared))) > 0,
            fam["positive_inner_cone"])
    distance_from = {}
    probes = []
    for start in sorted({0, centre, n - 1}):
        dist = bfs(indptr, indices, start, n)
        if int(dist.min()) < 0:
            raise ValueError("the site graph is disconnected")
        distance_from[start] = dist
        At, Bt = metric_tables(A, B, sites, start)
        for k in range(1, K + 1):
            cone = phi_sign_array(q * At - k * k, q * Bt) <= 0
            reached = dist <= k
            outer = int(np.count_nonzero(reached & ~cone))
            missing = np.flatnonzero(cone & ~reached)
            inner_misses = 0
            if inner_speed > 0:
                bound = k * k * inner_speed * inner_speed
                for i in missing.tolist():
                    if phi_sign(phi_sub(phi_scale(q, (int(At[i]), int(Bt[i]))), (bound, 0))) <= 0:
                        inner_misses += 1
            probes.append({"start": start, "layers": k, "reachable_count": int(reached.sum()),
                           "reachable_ids_sha256": digest(np.flatnonzero(reached).tolist()),
                           "outer_cone_violations": outer, "certified_inner_cone_misses": inner_misses,
                           "exact_finite_cone_missing_count": int(len(missing)),
                           "missing_ids_sha256": digest(missing.tolist())})
    R.check("reachability probes (cone sandwich, reachable and missing digests)", probes, fam["reachability_probes"])
    alpha = distance_from[centre]
    R.check("shell counts from the centre", np.bincount(alpha, minlength=K + 1)[:K + 1].tolist(),
            fam["shell_counts_from_centre"])

    # The layered order: centre diamonds, counts, exact strict pairs, fraction, dimension, clock.
    t2 = time.time()
    rows = {r["layers"]: r for r in fam["vertical_intervals"]}
    counts = {k: diamond_counts(alpha, k) for k in range(1, K + 1)}
    for k in range(1, K + 1):
        R.check(f"k={k} diamond counts by layer", counts[k], rows[k]["counts_by_layer"])
        R.check(f"k={k} diamond event count", sum(counts[k]), rows[k]["inclusive_event_count"])
    R.check("exact width = site count", n, fam["exact_width"])
    R.check("exact height in events = K + 1", K + 1, fam["exact_height_in_events"])
    if crow is not None:
        ci = crow["provenance"]["centre_interval"]
        R.check("carrier receipt: layered-order event-set digest of the K-diamond",
                diamond_event_digest(alpha, K), ci["layered_order_event_set_sha256"])
        R.check("carrier receipt: provenance interval equals the layered order",
                diamond_event_digest(alpha, K), ci["event_set_sha256"])
        R.check("carrier receipt: shell counts", np.bincount(alpha, minlength=K + 1)[:K + 1].tolist(),
                crow["provenance"]["shell_counts_from_centre"])
    if pairs:
        C = strict_pair_counts(indptr, indices, alpha, K, n)
        for k in range(1, K + 1):
            row = rows[k]
            N = sum(counts[k])
            if row.get("pair_counting") == "exact_all_pairs":
                R.check(f"k={k} strict pair count", C[k], row["strict_pair_count"])
            else:
                R.note(f"k={k} strict pair count (receipt row is a stratified estimate)", C[k])
            if N >= 2:
                f = Fraction(2 * C[k], N * (N - 1))
                if row.get("pair_counting") == "exact_all_pairs":
                    R.check(f"k={k} ordering fraction 2C/(N(N-1))", str(f), row["ordering_fraction"])
                    R.check(f"k={k} ordering fraction (12 digits)", rounded(float(f)), row["ordering_fraction_float"])
                    R.check(f"k={k} distance to the reference fraction",
                            rounded(float(f - REFERENCE_FRACTION[dim])), row["distance_to_reference"])
                    est = invert_myrheim_meyer(float(f))
                    R.check(f"k={k} inverted Myrheim-Meyer dimension",
                            None if est is None else rounded(est), row["myrheim_meyer_dimension"])
            if crow is not None:
                mrow = {r["layers"]: r for r in crow["manifold"]["intervals"]}[k]
                R.check(f"carrier receipt: k={k} strict pair count on the provenance order", C[k], mrow["strict_pair_count"])
        _log(log_out, f"q={q} dim={dim}: strict pairs {time.time() - t2:.1f}s (C_K = {C[K]})")
    Kp = K // 2
    NI, NJ = sum(counts[K]), sum(counts[Kp])
    clock = fam["count_clock"]
    R.check("count clock (N_K/N_K')^(1/4)", rounded((NI / NJ) ** 0.25), clock["count_clock"])
    R.check("count clock model-time ratio K/K'", rounded(K / Kp), clock["model_time_ratio"])
    R.check("count clock relative deviation", rounded((NI / NJ) ** 0.25 / (K / Kp) - 1.0), clock["relative_deviation"])
    R.check("count clock interval and reference counts", [NI, NJ], [clock["interval_count"], clock["reference_count"]])

    # Reads, audit chains, intervention, costs, and the stored log.
    if reads and crow is not None:
        t3 = time.time()
        forward = run_reads(indptr, indices, K)
        R.check("forward audit trace", forward["audit_trace_sha256"], crow["reads"]["forward"]["audit_trace_sha256"])
        R.check("forward layer digests", forward["layer_value_sha256"], crow["reads"]["forward"]["layer_value_sha256"])
        R.check("forward layer sums", forward["layer_value_sums"], crow["reads"]["forward"]["layer_value_sums"])
        R.check("forward event and read counts", [forward["event_count"], forward["authenticated_read_count"]],
                [crow["reads"]["forward"]["event_count"], crow["reads"]["forward"]["authenticated_read_count"]])
        interv = run_reads(indptr, indices, K, intervention=centre)
        R.check("intervention audit trace", interv["audit_trace_sha256"], crow["reads"]["intervention"]["audit_trace_sha256"])
        R.check("intervention layer digests", interv["layer_value_sha256"], crow["reads"]["intervention"]["layer_value_sha256"])
        support_rows = []
        probe_by_layer = {p["layers"]: p for p in probes if p["start"] == centre}
        for j in range(K + 1):
            deltas = [b - a for a, b in zip(forward["values"][j], interv["values"][j])]
            changed = [i for i, d in enumerate(deltas) if d != 0]
            row = {"round": j, "support_count": len(changed), "support_ids_sha256": digest(changed),
                   "positive_integer_delta_sum": sum(deltas), "positive_integer_delta_maximum": max(deltas),
                   "all_deltas_nonnegative": min(deltas) >= 0}
            theirs = crow["intervention"]["rounds"][j]
            R.check(f"intervention round {j} support and deltas", row,
                    {k_: theirs[k_] for k_ in row})
            if j >= 1:
                R.check(f"intervention round {j} support equals the future cone probe",
                        row["support_ids_sha256"], probe_by_layer[j]["reachable_ids_sha256"])
            support_rows.append(row)
        costs = operation_costs(indptr, forward["values"])
        theirs = crow["operation_costs"]
        R.check("operation costs (byte model)", {k_: costs[k_] for k_ in ("per_round", "total_reads", "total_writes",
                                                                        "total_read_bytes", "total_write_bytes", "total_bytes")},
                {k_: theirs[k_] for k_ in ("per_round", "total_reads", "total_writes",
                                           "total_read_bytes", "total_write_bytes", "total_bytes")})
        R.check("reads per round = 2E + n", all(r["reads"] == int(degrees.sum()) for r in costs["per_round"][1:]), True)
        _log(log_out, f"q={q}: reads, chains, intervention, costs {time.time() - t3:.1f}s")
        stored = load_log(q)
        if stored is not None and crow.get("stored_log"):
            t4 = time.time()
            R.check("stored log uncompressed digest", stored["uncompressed_sha256"], crow["stored_log"]["uncompressed_sha256"])
            R.check("stored log uncompressed bytes", stored["uncompressed_bytes"], crow["stored_log"]["uncompressed_bytes"])
            log = stored["payload"]
            R.check("stored log event count", len(log["events"]), crow["stored_log"]["event_count"])
            written = [[int(ev["write"][2]) for ev in log["events"] if ev["event"][0] == j] for j in range(K + 1)]
            R.check("stored log write values equal the rebuilt read-law values", written, forward["values"])
            prov = provenance_from_log(log, centre)
            theirs = crow["provenance"]
            for key in ("edge_count", "single_writer", "all_reads_resolved", "parents_precede_children",
                        "derived_rank_equals_round", "read_relation_identical_across_rounds",
                        "read_relation_sha256", "future_cone_counts", "future_cone_sha256"):
                R.check(f"provenance from the log: {key}", prov[key], theirs[key])
            R.check("provenance from the log: read relation equals the neighbour digest", prov["read_relation_sha256"], nb_sha)
            R.check("provenance from the log: centre interval event set", prov["interval_event_set_sha256"],
                    theirs["centre_interval"]["event_set_sha256"])
            R.check("provenance from the log: centre interval counts by round", prov["interval_counts_by_round"],
                    theirs["centre_interval"]["counts_by_round"])
            R.check("provenance from the log: interval equals the layered order rebuilt here",
                    prov["interval_event_set_sha256"], diamond_event_digest(alpha, K))
            _log(log_out, f"q={q}: stored log replay {time.time() - t4:.1f}s")
    _log(log_out, f"q={q} dim={dim}: {sum(1 for r in R.rows if r['receipt'] is not None)} comparisons, "
                  f"{'all agree' if R.ok else 'DISAGREEMENT'}, {time.time() - t0:.1f}s")
    return R


def fibonacci_next(q: int) -> int:
    a, b = 0, 1
    while b < q:
        a, b = b, a + b
    if b != q:
        raise ValueError("q is not a Fibonacci number")
    return a + b


def _log(stream, msg: str) -> None:
    if stream is not None:
        print(msg, file=stream, flush=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--q", type=int, nargs="*", default=[5, 8], help="levels to rebuild (default 5 8)")
    parser.add_argument("--dim", type=int, default=3, choices=(1, 2, 3), help="family (3 = source net; 2, 1 = controls)")
    parser.add_argument("--no-pairs", action="store_true", help="skip the exact strict pair count")
    parser.add_argument("--pairs", action="store_true", help="force the exact strict pair count at every requested level")
    parser.add_argument("--no-reads", action="store_true", help="skip the reads, audit chains and log replay")
    parser.add_argument("--reads", action="store_true", help="force the reads at every requested level")
    parser.add_argument("--json", type=Path, default=None, help="write the comparison rows to this file")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    pairs = True if args.pairs else (False if args.no_pairs else None)
    reads = True if args.reads else (False if args.no_reads else None)
    reports = []
    ok = True
    for q in args.q:
        R = build_level(q, dim=args.dim, pairs=pairs, reads=reads, log_out=None if args.quiet else sys.stderr)
        reports.append({"q": q, "dimension": args.dim, "agree": R.ok,
                        "comparisons": sum(1 for r in R.rows if r["receipt"] is not None),
                        "disagreements": [r for r in R.rows if not r["agree"]], "rows": R.rows})
        ok = ok and R.ok
    summary = {"status": "CAUSAL_POSET_REBUILT_AND_EQUAL_TO_RECEIPTS" if ok else "CAUSAL_POSET_DISAGREES_WITH_RECEIPTS",
               "levels": [{k: v for k, v in r.items() if k != "rows"} for r in reports]}
    if args.json is not None:
        args.json.write_text(json.dumps(reports, indent=1, default=str) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

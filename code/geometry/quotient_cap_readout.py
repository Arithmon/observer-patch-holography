#!/usr/bin/env python3
"""Machine receipts for the quotient-intrinsic geometry producer (GitHub #523).

Implements the objects of Theorem 4.3c (quotient-intrinsic geometry producer
and dimension-selection boundary) in `paper/tex_fragments/PAPER.tex`, the
technical fragment included by The spacetime and Einstein paper; the clause
names below are descriptive and carry no LaTeX label in the .tex sources:

* a finite transactional quotient repair system whose record layer binds one
  record token per 2-cell of an abstract patch adjacency (the adjacency is
  never consumed directly by the readout: the incidence complex is recomputed
  from the repaired record tokens of the normal form);
* the support-visible incidence complex `K(W)`, with the invariance
  checks (schedule, gauge, refinement);
* the computable spherical-incidence receipt (i) of Theorem 4.3c
  (closed combinatorial surface, connected, Euler characteristic 2, coherent
  orientation) and the wrong-beta scalar KMS receipt clause;
* the topology-production step (surface classification by orientability +
  Euler characteristic);
* the conformal production step: modular cross-ratios against a fixed gauge
  triple reconstruct the celestial embedding, cap normals `n_C` are produced
  from boundary circles by the round-cap-normal formula, and the residuals
  shrink under refinement;
* the underdetermination countermodels of Theorem 4.3c: six repair systems
  with one rewrite signature `(3, 1)` (patches in the seed record, repair
  steps) over S^2, T^2 (Csaszar torus), RP^2 (hemi-icosahedron), the Klein
  bottle, the 2-skeleton of the boundary of the 4-simplex, and a wedge of
  two spheres, distinguished only by the receipts.  They show that
  confluence alone does not select topology, orientability, dimension, or
  modular normalization.  On a connected closed surface the orientability
  clause of the spherical-incidence receipt follows from chi = 2 by the
  classification of surfaces, so both nonorientable systems also fail that
  receipt at its Euler clause;
* an orientation-framing countermodel: the two coherent orientations of the
  S^2 system share every confluence-level datum, and a mirror pair of
  embeddings carries complex-conjugate cross-ratios, so the framing enters
  only through its separately typed receipt.  No MGNS-state countermodel is
  built here, and the receipt set is one sufficient certificate whose
  semantic minimality is unproved.
"""

from __future__ import annotations

import itertools
import random
from dataclasses import dataclass, field
from numbers import Integral, Real

import numpy as np

if __package__:
    from .conformal_readout import (
        cross_ratio, cross_ratio_receipts, fit_cap, mobius_normalize,
        produced_cap_normal, reconstruct_from_cross_ratios, stereographic,
    )
    from .finite_incidence import (
        IncidenceComplex, certify_midpoint_subdivision, complexes_equal,
        complexes_isomorphic_under, dimension, euler_characteristic,
        incidence_complex, refinement_is_simplicial, validate_complex, validate_records,
    )
else:
    from conformal_readout import (
        cross_ratio, cross_ratio_receipts, fit_cap, mobius_normalize,
        produced_cap_normal, reconstruct_from_cross_ratios, stereographic,
    )
    from finite_incidence import (
        IncidenceComplex, certify_midpoint_subdivision, complexes_equal,
        complexes_isomorphic_under, dimension, euler_characteristic,
        incidence_complex, refinement_is_simplicial, validate_complex, validate_records,
    )

TOL = 1e-9


# ---------------------------------------------------------------------------
# abstract 2-complexes (record layers)
# ---------------------------------------------------------------------------

def icosahedron() -> tuple[list[tuple[int, int, int]], np.ndarray]:
    """Return (triangles, vertex coordinates on the unit sphere)."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    verts = []
    for a, b in itertools.product((-1.0, 1.0), (-phi, phi)):
        verts += [(0.0, a, b), (a, b, 0.0), (b, 0.0, a)]
    coords = np.array(verts)
    coords /= np.linalg.norm(coords, axis=1, keepdims=True)
    # triangles: triples of mutually nearest vertices (edge length = 2/sqrt(phi^2+1))
    edge2 = 4.0 / (phi * phi + 1.0)
    n = len(coords)
    tris = []
    for i, j, k in itertools.combinations(range(n), 3):
        d = (
            np.sum((coords[i] - coords[j]) ** 2),
            np.sum((coords[j] - coords[k]) ** 2),
            np.sum((coords[i] - coords[k]) ** 2),
        )
        if all(abs(x - edge2) < 1e-6 for x in d):
            tris.append((i, j, k))
    assert len(tris) == 20
    return tris, coords


def csaszar_torus() -> list[tuple[int, int, int]]:
    """The Moebius 7-vertex triangulation of the torus (complete graph K7):
    triangles {i, i+1, i+3} and {i, i+2, i+3} mod 7."""
    tris = [tuple(sorted(((i) % 7, (i + 1) % 7, (i + 3) % 7))) for i in range(7)]
    tris += [tuple(sorted(((i) % 7, (i + 2) % 7, (i + 3) % 7))) for i in range(7)]
    assert len(set(tris)) == 14
    return tris


def boundary_4_simplex_2_skeleton() -> list[tuple[int, int, int]]:
    """All ten triangles on five vertices: the 2-skeleton of boundary(Delta^4).

    Locally three-dimensional (every edge lies in three triangles), so the
    closed-surface clause of the spherical incidence receipt must fail.
    """
    return list(itertools.combinations(range(5), 3))


def wedge_of_two_spheres() -> list[tuple[int, int, int]]:
    """Two icosahedra glued at one vertex: connected, chi = 3, nonmanifold."""
    tris, _ = icosahedron()
    shift = 12

    def relabel(v: int) -> int:
        # vertex 0 of the second copy is identified with vertex 0 of the first
        return 0 if v == 0 else v + shift - 1

    second = [tuple(relabel(v) for v in t) for t in tris]
    return tris + second


def real_projective_plane() -> list[tuple[int, int, int]]:
    """The hemi-icosahedron: the antipodal quotient of `icosahedron()`.

    Antipodal vertex pairs of the icosahedron become six vertices and
    antipodal face pairs become ten triangles. The edge graph is the
    complete graph K6 and chi = 1: the vertex-minimal triangulation of the
    real projective plane, a closed nonorientable surface.
    """
    tris, coords = icosahedron()
    n = len(coords)
    antipode = [int(np.argmin(np.linalg.norm(coords + coords[i], axis=1))) for i in range(n)]
    classes = sorted({min(i, antipode[i]) for i in range(n)})
    label = {i: classes.index(min(i, antipode[i])) for i in range(n)}
    quotient = sorted({tuple(sorted(label[v] for v in t)) for t in tris})
    assert len(classes) == 6 and len(quotient) == 10
    return quotient


def klein_bottle() -> list[tuple[int, int, int]]:
    """A nine-vertex Klein bottle from a 3 x 3 grid with a glide identification.

    The plane modulo the translation (x, y) -> (x, y + 3) and the glide
    reflection (x, y) -> (x + 3, -y) is a Klein bottle. Each unit square is
    cut by one diagonal whose direction alternates with the column; the
    glide has an odd shift and reflects y, so it maps the diagonal pattern to
    itself and the triangulation descends to the quotient: 9 vertices,
    27 edges, 18 triangles, chi = 0. The vertex-minimal Klein bottle has 8
    vertices; this closed form is used because it can be checked by hand.
    """
    m = n = 3

    def vertex(x: int, y: int) -> int:
        sheet, column = divmod(x, m)
        return column * n + ((-1) ** sheet * y) % n

    tris = []
    for x, y in itertools.product(range(m), range(n)):
        a, b = vertex(x, y), vertex(x + 1, y)
        c, d = vertex(x, y + 1), vertex(x + 1, y + 1)
        if x % 2 == 0:
            tris += [(a, b, d), (a, c, d)]
        else:
            tris += [(a, b, c), (b, c, d)]
    tris = [tuple(sorted(t)) for t in tris]
    assert len(set(tris)) == 18
    return tris


# ---------------------------------------------------------------------------
# finite transactional quotient repair system
# ---------------------------------------------------------------------------

@dataclass
class RepairSystem:
    """Patch bits with equality constraints along record tokens.

    Records demand that all member patch bits agree. The admitted state has
    at most one bit equal to one; every record has at least two members.
    Majority repair with a zero-valued tie break erases that defect in one
    effective update, under every complete record schedule. This is a finite
    one-defect fixture, not a confluence theorem for arbitrary majority states.
    """

    records: list[tuple[int, ...]]
    seed_record: int = 0
    state: dict[int, int] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        self._validate_records()
        patches = sorted({p for rec in self.records for p in rec})
        self.state = {p: 0 for p in patches}
        # seed the conflict: flip the highest-label patch of the seed record
        bad = max(self.records[self.seed_record])
        self.state[bad] = 1

    def _validate_records(self) -> None:
        validate_records(self.records)
        if (not self.records or any(len(rec) < 2 for rec in self.records)
                or not isinstance(self.seed_record, Integral)
                or isinstance(self.seed_record, (bool, np.bool_))
                or not 0 <= self.seed_record < len(self.records)):
            raise ValueError("repair fixture requires nonsingleton records and a valid seed")

    def conflicted_records(self) -> list[int]:
        return [
            k for k, rec in enumerate(self.records)
            if len({self.state[p] for p in rec}) > 1
        ]

    def repair(self, schedule: list[int] | None = None) -> dict[int, int]:
        """Run transactional repair; return the quotient normal form."""
        self._validate_records()
        order = list(schedule) if schedule is not None else list(range(len(self.records)))
        if (len(order) != len(self.records)
                or any(not isinstance(k, Integral) or isinstance(k, (bool, np.bool_)) for k in order)
                or set(order) != set(range(len(self.records)))):
            raise ValueError("schedule must contain every record index exactly once")
        patches = {p for rec in self.records for p in rec}
        if (set(self.state) != patches
                or any(not isinstance(k, Integral) or isinstance(k, (bool, np.bool_))
                       for k in self.state)
                or any(not isinstance(v, Integral) or isinstance(v, (bool, np.bool_))
                       or v not in (0, 1) for v in self.state.values())
                or sum(int(v) for v in self.state.values()) > 1):
            raise ValueError("repair fixture only admits zero or one seeded defect")
        working = dict(self.state)
        progress = True
        while progress:
            progress = False
            for k in order:
                rec = self.records[k]
                vals = [working[p] for p in rec]
                if len(set(vals)) > 1:
                    # Supplied fixture law: majority with a zero-valued tie break.
                    counts = {v: vals.count(v) for v in set(vals)}
                    best = min(sorted(counts), key=lambda v: (-counts[v], v))
                    for p in rec:
                        working[p] = best
                    progress = True
        return working

    def rewrite_signature(self) -> tuple[int, int]:
        """Confluence-level invariant: (#patches touched by conflict, #steps).

        Isomorphic across all six countermodel systems by construction: the
        conflict component is a single record with three patches everywhere.
        """
        rec = self.records[self.seed_record]
        return (len(rec), 1)


def normal_form_records(system: RepairSystem, normal_form: dict[int, int]) -> list[tuple[int, ...]]:
    """The record tokens carried by the repaired normal form.

    A record is quotient-visible exactly when its constraint is satisfied on
    the normal form (all tokens after successful repair). The incidence
    readout below consumes ONLY this list, i.e. normal-form data.
    """
    return [rec for rec in system.records
            if len({normal_form[p] for p in rec}) == 1]


# ---------------------------------------------------------------------------
# support-visible incidence complex and receipts
# ---------------------------------------------------------------------------

def is_connected(K: IncidenceComplex) -> bool:
    validate_complex(K)
    if not K.vertices:
        return False
    adj: dict[int, set[int]] = {v: set() for v in K.vertices}
    for e in K.edges:
        a, b = tuple(e)
        adj[a].add(b)
        adj[b].add(a)
    seen = {K.vertices[0]}
    stack = [K.vertices[0]]
    while stack:
        for w in adj[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(K.vertices)


def is_closed_surface(K: IncidenceComplex) -> bool:
    """Every edge in exactly two triangles and every vertex link a single cycle."""
    validate_complex(K)
    if not K.triangles or K.higher_simplices:
        return False
    for e in K.edges:
        cofaces = [t for t in K.triangles if e < t]
        if len(cofaces) != 2:
            return False
    for v in K.vertices:
        star = [t for t in K.triangles if v in t]
        if not star:
            return False
        link_edges = [tuple(sorted(t - {v})) for t in star]
        nodes = sorted({x for le in link_edges for x in le})
        deg = {x: sum(1 for le in link_edges if x in le) for x in nodes}
        if any(d != 2 for d in deg.values()):
            return False
        # single cycle: connected link
        adj = {x: set() for x in nodes}
        for a, b in link_edges:
            adj[a].add(b)
            adj[b].add(a)
        seen, stack = {nodes[0]}, [nodes[0]]
        while stack:
            for w in adj[stack.pop()]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != len(nodes):
            return False
    return True


def orient(K: IncidenceComplex) -> list[tuple[int, int, int]] | None:
    """Return coherently oriented triangles, or None if nonorientable."""
    validate_complex(K)
    if K.higher_simplices:
        return None
    tris = sorted(tuple(sorted(t)) for t in K.triangles)
    if not tris:
        return None
    oriented: dict[tuple, tuple] = {}
    first = tris[0]
    oriented[first] = first
    stack = [first]
    remaining = set(tris[1:])
    while stack:
        t = stack.pop()
        a, b, c = oriented[t]
        # induced edge orientations of t
        induced = {(a, b), (b, c), (c, a)}
        for u in list(remaining):
            shared = set(t) & set(u)
            if len(shared) == 2:
                x, y = tuple(shared)
                (z,) = set(u) - shared
                # neighbor must induce the opposite orientation on the shared edge
                if (x, y) in induced:
                    ori = (y, x, z)
                else:
                    ori = (x, y, z)
                oriented[u] = ori
                remaining.discard(u)
                stack.append(u)
    if remaining:
        return None  # disconnected; caller checks connectivity separately
    # verify global coherence
    edge_orientations: dict[frozenset, list] = {}
    for t, (a, b, c) in oriented.items():
        for x, y in ((a, b), (b, c), (c, a)):
            edge_orientations.setdefault(frozenset((x, y)), []).append((x, y))
    for e, oris in edge_orientations.items():
        if len(oris) == 2 and oris[0] == oris[1]:
            return None  # incoherent: nonorientable
    return list(oriented.values())


def spherical_incidence_receipt(K: IncidenceComplex) -> bool:
    """SphInc: connected closed orientable combinatorial surface with chi = 2."""
    return (
        is_connected(K)
        and is_closed_surface(K)
        and euler_characteristic(K) == 2
        and orient(K) is not None
    )


def classify_surface(K: IncidenceComplex) -> str:
    """Topology-production step of Theorem 4.3c: classification by receipts.

    Orientability and Euler characteristic fix a connected closed surface up
    to homeomorphism: the orientable surface of genus g has chi = 2 - 2g, and
    the nonorientable surface with k crosscaps has chi = 2 - k (k = 1 is the
    real projective plane, k = 2 the Klein bottle).
    """
    if not is_connected(K):
        return "DISCONNECTED"
    if not is_closed_surface(K):
        return "NOT_A_CLOSED_SURFACE"
    chi = euler_characteristic(K)
    orientable = orient(K) is not None
    if orientable and chi == 2:
        return "S2"
    if orientable and chi == 0:
        return "T2"
    if orientable:
        return f"GENUS_{(2 - chi) // 2}"
    if chi == 1:
        return "RP2"
    if chi == 0:
        return "KLEIN_BOTTLE"
    return f"NONORIENTABLE_GENUS_{2 - chi}"


# ---------------------------------------------------------------------------
# orientation framing (separately typed receipt of Theorem 4.3c)
# ---------------------------------------------------------------------------

def reverse_orientation(oriented: list[tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    """The opposite framing: every oriented triangle with its cycle reversed."""
    return [(b, a, c) for a, b, c in oriented]


def orientation_is_coherent(oriented: list[tuple[int, int, int]]) -> bool:
    """Every edge is traversed exactly once in each direction."""
    if not oriented or any(len(set(t)) != 3 for t in oriented):
        return False
    if len({frozenset(t) for t in oriented}) != len(oriented):
        return False
    directed: dict[tuple[int, int], int] = {}
    for a, b, c in oriented:
        for x, y in ((a, b), (b, c), (c, a)):
            directed[(x, y)] = directed.get((x, y), 0) + 1
    return all(
        count == 1 and directed.get((y, x)) == 1
        for (x, y), count in directed.items()
    )


def is_outward_framing(oriented: list[tuple[int, int, int]], coords: np.ndarray) -> bool:
    """The framing agrees with the outward normals of an embedding star-shaped about the origin."""
    if not orientation_is_coherent(oriented):
        return False
    if coords.ndim != 2 or coords.shape[1] != 3 or not np.isfinite(coords).all():
        return False
    if any(v < 0 or v >= len(coords) for t in oriented for v in t):
        return False
    for a, b, c in oriented:
        normal = np.cross(coords[b] - coords[a], coords[c] - coords[a])
        signed_volume = float(np.dot(normal, coords[a] + coords[b] + coords[c]))
        if not np.isfinite(signed_volume) or signed_volume <= 0.0:
            return False
    return True


# ---------------------------------------------------------------------------
# modular cross-ratio production (Theorem thm:conformal-cap-production)
# ---------------------------------------------------------------------------

def minkowski(u: np.ndarray, v: np.ndarray) -> float:
    return float(-u[0] * v[0] + np.dot(u[1:], v[1:]))


def kms_receipt(clock_scale: float, beta_target: float = 2.0 * np.pi, tol: float = 1e-9) -> bool:
    """Wrong-normalization separation clause: the independently normalized
    geometric comparison certifies exactly the declared modular temperature."""
    for value in (clock_scale, beta_target, tol):
        if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
                or not np.isfinite(value) or value <= 0):
            raise ValueError("clock, target and tolerance must be finite positive reals")
    if tol > 1e-6:
        raise ValueError("comparison tolerance exceeds the finite diagnostic's limit")
    return abs(clock_scale - beta_target) < tol


# ---------------------------------------------------------------------------
# invariance drivers (Lemma lem:incidence-invariance)
# ---------------------------------------------------------------------------

def readout_from_system(system: RepairSystem, schedule: list[int] | None = None) -> IncidenceComplex:
    nf = system.repair(schedule)
    return incidence_complex(normal_form_records(system, nf))


def gauge_relabel(records: list[tuple[int, ...]], perm: dict[int, int]) -> list[tuple[int, ...]]:
    return [tuple(sorted(perm[p] for p in rec)) for rec in records]


def refinement_subdivide(records: list[tuple[int, ...]]) -> tuple[list[tuple[int, ...]], dict[int, int]]:
    """Midpoint subdivision of a triangular record layer (one refinement
    stage) together with the coarse-graining projection on patch labels."""
    K = incidence_complex(records)
    if not K.triangles or any(len(rec) != 3 for rec in records):
        raise ValueError("subdivision producer requires nonempty triangular records")
    next_label = int(max(K.vertices)) + 1
    midpoint: dict[frozenset, int] = {}
    projection: dict[int, int] = {v: v for v in K.vertices}
    for e in sorted(K.edges, key=lambda e: tuple(sorted(e))):
        midpoint[e] = next_label
        a, _b = tuple(sorted(e))
        projection[next_label] = a  # midpoints project to an endpoint
        next_label += 1
    refined = []
    for rec in records:
        a, b, c = sorted(rec)
        mab = midpoint[frozenset((a, b))]
        mbc = midpoint[frozenset((b, c))]
        mac = midpoint[frozenset((a, c))]
        refined += [
            (a, mab, mac), (b, mab, mbc), (c, mbc, mac), (mab, mbc, mac),
        ]
    if not certify_midpoint_subdivision(K, incidence_complex(refined), midpoint):
        raise RuntimeError("subdivision failed its complete carrier check")
    return refined, projection

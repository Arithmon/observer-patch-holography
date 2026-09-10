#!/usr/bin/env python3
"""Exact certificate: the A1 twelve-port carrier as the coset geometry of a
supplied binary icosahedral source.

The source manifest declares one finite group, SL(2, F_5), as the exact finite
representative of the binary icosahedral group 2I.  The identification
SL(2, 5) = 2I is classical and is cited here without formalization.  The
verifier:

* constructs SL(2, F_5) and checks its order, centre {+I, -I}, unique
  involution, perfectness, and binary icosahedral element-order profile;
* forms the central quotient G of order sixty;
* enumerates every placement (C5, C3, C2) of cyclic subgroups of G, 900 in
  all, and builds the coset geometry with vertices G/C5, faces G/C3, edges
  G/C2, and incidence by nonempty coset intersection;
* classifies every placement by exact checks on all three incidences;
* compares every compatible placement with the committed A1 packet, read
  from the Lean sources, through an explicit relabelling.

The local counts are identical for all 900 placements.  Exactly 120
placements close into the committed oriented carrier.  They are the
placements admitting generators x, y, z of orders 2, 3, 5 with xyz = 1, and
they correspond bijectively to the solutions of a^2 = b^3 = c^5 = abc in
SL(2, F_5).  They form two orbits under inner automorphisms and one orbit
under all automorphisms of the supplied group.

The reconstruction is conditional on the supplied group and factors through
its central quotient.  It makes no physical selection of the source.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

SCHEMA = "oph.coset_carrier_source.v1"
RECEIPT_SCHEMA = "oph.coset_carrier_receipt.v1"
NEGATIVE_SCHEMA = "oph.coset_carrier_negative_controls.v1"
GENERATED_BY = "code/a5_closure/coset_carrier_certificate.py"
REPO_ROOT = Path(__file__).resolve().parents[2]

FIXTURE_DECLARATIONS = {
    "adjacency": ("Lean/Screen/PortFrameGram.lean", "neighbors"),
    "port_action": ("Lean/Screen/A5PortAction.lean", "perms"),
    "oriented_faces": ("Lean/ObserverPatchHolography/CoreAxioms.lean", "orientedFaces"),
}

# Element orders of the binary icosahedral group: order -> number of elements.
BINARY_ICOSAHEDRAL_ORDER_PROFILE = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}

# Words that name target structure.  None of them may occur in the source block.
SOURCE_FIREWALL_TOKENS = ("icosa", "dodeca", "a5", "port", "vertex", "edge", "face", "neighbor", "orient", "carrier")

Matrix = tuple[int, int, int, int]
Placement = tuple[frozenset, frozenset, frozenset]


class CertificateError(ValueError):
    """Fail-closed manifest, construction, or receipt error carrying a stable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def require_exact_keys(value: Mapping[str, Any], expected: set[str], path: str) -> None:
    actual = set(value)
    require(
        actual == expected,
        "SCHEMA_FIELDS",
        f"{path} fields must be exactly {sorted(expected)}; "
        f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}",
    )


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def load_json(path: Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("JSON_READ", f"cannot read {path}: {exc}") from exc


def write_json(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Source manifest
# ---------------------------------------------------------------------------


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(isinstance(manifest, Mapping), "SCHEMA", "the manifest must be a JSON object")
    require_exact_keys(manifest, {"schema", "source", "comparison_fixture"}, "$")
    require(manifest["schema"] == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    source = manifest["source"]
    require(isinstance(source, Mapping), "SCHEMA", "source must be an object")
    require_exact_keys(source, {"group", "prime"}, "$.source")
    for key, value in source.items():
        if isinstance(value, str):
            lowered = value.lower()
            hit = [token for token in SOURCE_FIREWALL_TOKENS if token in lowered]
            require(not hit, "SOURCE_FIREWALL", f"$.source.{key} names target structure {hit}")
    require(source["group"] == "SL2", "SOURCE_GROUP", "the declared source family must be SL2")
    prime = source["prime"]
    require(isinstance(prime, int) and not isinstance(prime, bool) and prime >= 2, "SOURCE_PRIME", "prime must be an integer >= 2")
    fixture = manifest["comparison_fixture"]
    require(isinstance(fixture, Mapping), "SCHEMA", "comparison_fixture must be an object")
    require_exact_keys(fixture, set(FIXTURE_DECLARATIONS), "$.comparison_fixture")
    for role, (module, declaration) in FIXTURE_DECLARATIONS.items():
        entry = fixture[role]
        require(isinstance(entry, Mapping), "SCHEMA", f"comparison_fixture.{role} must be an object")
        require_exact_keys(entry, {"module", "declaration"}, f"$.comparison_fixture.{role}")
        require(
            (entry["module"], entry["declaration"]) == (module, declaration),
            "FIXTURE_DECLARATION",
            f"comparison_fixture.{role} must name {module} `{declaration}`",
        )


# ---------------------------------------------------------------------------
# The supplied group and its central quotient
# ---------------------------------------------------------------------------


def element_order(mul: Callable[[Any, Any], Any], identity: Any, x: Any) -> int:
    k, y = 1, x
    while y != identity:
        y = mul(x, y)
        k += 1
    return k


def cyclic_subgroup(mul: Callable[[Any, Any], Any], identity: Any, x: Any) -> frozenset:
    out, y = {identity}, x
    while y != identity:
        out.add(y)
        y = mul(x, y)
    return frozenset(out)


@dataclass(frozen=True)
class Source:
    p: int
    elements: tuple[Matrix, ...]

    @property
    def identity(self) -> Matrix:
        return (1, 0, 0, 1)

    @property
    def minus_identity(self) -> Matrix:
        return (self.p - 1, 0, 0, self.p - 1)

    def mul(self, x: Matrix, y: Matrix) -> Matrix:
        a, b, c, d = x
        e, f, g, h = y
        p = self.p
        return ((a * e + b * g) % p, (a * f + b * h) % p, (c * e + d * g) % p, (c * f + d * h) % p)

    def inv(self, x: Matrix) -> Matrix:
        a, b, c, d = x
        p = self.p
        return (d % p, (-b) % p, (-c) % p, a % p)

    def neg(self, x: Matrix) -> Matrix:
        return tuple((-t) % self.p for t in x)


def build_source(manifest: Mapping[str, Any]) -> Source:
    validate_manifest(manifest)
    p = manifest["source"]["prime"]
    elements = tuple(m for m in itertools.product(range(p), repeat=4) if (m[0] * m[3] - m[1] * m[2]) % p == 1)
    source = Source(p, elements)
    require(len(elements) == 120, "SOURCE_ORDER", f"SL(2, F_{p}) has order {len(elements)}, expected 120")
    centre = sorted(z for z in elements if all(source.mul(z, g) == source.mul(g, z) for g in elements))
    require(centre == sorted([source.identity, source.minus_identity]), "SOURCE_CENTRE", f"centre is {centre}")
    profile = Counter(element_order(source.mul, source.identity, g) for g in elements)
    require(dict(profile) == BINARY_ICOSAHEDRAL_ORDER_PROFILE, "SOURCE_ORDER_PROFILE", f"element-order profile {dict(profile)}")
    derived = {source.mul(source.mul(x, y), source.mul(source.inv(x), source.inv(y))) for x in elements for y in elements}
    while True:
        grown = derived | {source.mul(a, b) for a in derived for b in derived}
        if grown == derived:
            break
        derived = grown
    require(len(derived) == 120, "SOURCE_NOT_PERFECT", f"derived subgroup has order {len(derived)}")
    return source


@dataclass(frozen=True)
class Quotient:
    """G = SL(2, F_p)/{+I, -I}, each coset stored by its lexicographically least matrix."""

    source: Source
    elements: tuple[Matrix, ...]

    @property
    def identity(self) -> Matrix:
        return self.canon(self.source.identity)

    def canon(self, x: Matrix) -> Matrix:
        return min(x, self.source.neg(x))

    def mul(self, x: Matrix, y: Matrix) -> Matrix:
        return self.canon(self.source.mul(x, y))

    def inv(self, x: Matrix) -> Matrix:
        return self.canon(self.source.inv(x))


def central_quotient(source: Source) -> tuple[Quotient, dict[str, Any]]:
    quotient = Quotient(source, tuple(sorted({min(g, source.neg(g)) for g in source.elements})))
    G = quotient.elements
    require(len(G) == 60, "QUOTIENT_ORDER", f"central quotient has order {len(G)}")
    classes: list[frozenset] = []
    seen: set[Matrix] = set()
    for g in G:
        if g in seen:
            continue
        cls = frozenset(quotient.mul(quotient.mul(h, g), quotient.inv(h)) for h in G)
        seen |= cls
        classes.append(cls)
    sizes = sorted(len(c) for c in classes)
    require(sizes == [1, 12, 12, 15, 20], "QUOTIENT_CLASSES", f"conjugacy class sizes {sizes}")
    nontrivial = [len(c) for c in classes if quotient.identity not in c]
    normal_orders = sorted({1 + sum(pick) for r in range(len(nontrivial) + 1) for pick in itertools.combinations(nontrivial, r)
                            if 60 % (1 + sum(pick)) == 0})
    require(normal_orders == [1, 60], "QUOTIENT_NOT_SIMPLE", f"class unions of admissible order {normal_orders}")
    return quotient, {"order": 60, "conjugacy_class_sizes": sizes, "simple": True}


# ---------------------------------------------------------------------------
# Placements and their coset geometry
# ---------------------------------------------------------------------------


def placement_key(placement: Placement) -> tuple:
    return tuple(tuple(sorted(H)) for H in placement)


def cyclic_subgroups(quotient: Quotient, order: int) -> list[frozenset]:
    found = {
        cyclic_subgroup(quotient.mul, quotient.identity, g)
        for g in quotient.elements
        if element_order(quotient.mul, quotient.identity, g) == order
    }
    return sorted(found, key=sorted)


def subgroup_data(quotient: Quotient) -> tuple[dict[int, list[frozenset]], dict[str, Any]]:
    subs = {n: cyclic_subgroups(quotient, n) for n in (5, 3, 2)}
    counts = {n: len(subs[n]) for n in subs}
    require(counts == {5: 6, 3: 10, 2: 15}, "SUBGROUP_COUNT", f"cyclic subgroup counts {counts}")
    source = quotient.source
    for n in (5, 3, 2):
        for H in subs[n]:
            preimage = [g for g in source.elements if quotient.canon(g) in H]
            orders = {element_order(source.mul, source.identity, g) for g in preimage}
            require(
                len(preimage) == 2 * n and 2 * n in orders,
                "PREIMAGE_NOT_CYCLIC",
                f"the preimage of a C{n} is not cyclic of order {2 * n}",
            )
    return subs, {"C5": 6, "C3": 10, "C2": 15, "placements": 900, "preimages_cyclic_of_orders": [10, 6, 4]}


@dataclass(frozen=True)
class Incidence:
    """Coset geometry of one placement.  Incidence is nonempty coset intersection,
    so every relation is {(coset of g, coset of g) : g in G}."""

    edge_vertices: tuple[frozenset, ...]
    face_vertices: tuple[frozenset, ...]
    edge_faces: tuple[frozenset, ...]
    chambers: tuple[tuple[int, int, int], ...]
    vertex_representatives: tuple[Matrix, ...]

    @property
    def n_vertices(self) -> int:
        return len(self.vertex_representatives)


def coset_index(quotient: Quotient, H: frozenset) -> tuple[dict[Matrix, int], list[Matrix]]:
    """Left cosets gH; each is labelled in order of its least element."""
    index: dict[Matrix, int] = {}
    representatives: list[Matrix] = []
    for g in quotient.elements:
        if g not in index:
            for h in H:
                index[quotient.mul(g, h)] = len(representatives)
            representatives.append(g)
    return index, representatives


def build_incidence(quotient: Quotient, placement: Placement) -> Incidence:
    H5, H3, H2 = placement
    iV, vertex_reps = coset_index(quotient, H5)
    iF, face_reps = coset_index(quotient, H3)
    iE, edge_reps = coset_index(quotient, H2)
    counts = (len(vertex_reps), len(edge_reps), len(face_reps))
    require(counts == (12, 30, 20), "COSET_COUNTS", f"coset counts {counts}")
    ev = [set() for _ in edge_reps]
    fv = [set() for _ in face_reps]
    ef = [set() for _ in edge_reps]
    chambers = set()
    for g in quotient.elements:
        v, e, f = iV[g], iE[g], iF[g]
        ev[e].add(v)
        fv[f].add(v)
        ef[e].add(f)
        chambers.add((v, e, f))
    return Incidence(
        tuple(frozenset(x) for x in ev),
        tuple(frozenset(x) for x in fv),
        tuple(frozenset(x) for x in ef),
        tuple(sorted(chambers)),
        tuple(vertex_reps),
    )


def local_counts(inc: Incidence) -> tuple[tuple[int, ...], ...]:
    """(faces per vertex, vertices per face, edges per vertex, vertices per edge,
    faces per edge, edges per face), each as its set of observed values."""
    faces_at = Counter(v for verts in inc.face_vertices for v in verts)
    edges_at = Counter(v for pair in inc.edge_vertices for v in pair)
    edges_on = Counter(f for faces in inc.edge_faces for f in faces)
    return (
        tuple(sorted(set(faces_at.values()))),
        tuple(sorted({len(x) for x in inc.face_vertices})),
        tuple(sorted(set(edges_at.values()))),
        tuple(sorted({len(x) for x in inc.edge_vertices})),
        tuple(sorted({len(x) for x in inc.edge_faces})),
        tuple(sorted(set(edges_on.values()))),
    )


def adjacency(inc: Incidence) -> tuple[frozenset, ...]:
    adj: list[set[int]] = [set() for _ in range(inc.n_vertices)]
    for pair in inc.edge_vertices:
        if len(pair) == 2:
            a, b = sorted(pair)
            adj[a].add(b)
            adj[b].add(a)
    return tuple(frozenset(x) for x in adj)


def oriented_faces(inc: Incidence) -> tuple[tuple[int, int, int], ...]:
    """Orient each face by the chambers that share one coset representative.

    A chamber (xC5, xC2, xC3) contributes the directed edge from xC5 to the other
    end of xC2 on the face xC3.  Each face is written from its least vertex."""
    darts: dict[int, set[tuple[int, int]]] = {}
    for v, e, f in inc.chambers:
        others = inc.edge_vertices[e] - {v}
        require(len(others) == 1, "ORIENTATION_INCOHERENT", f"chamber edge {e} has no second vertex")
        darts.setdefault(f, set()).add((v, next(iter(others))))
    faces = []
    for f, verts in enumerate(inc.face_vertices):
        succ = dict(darts.get(f, ()))
        require(len(succ) == 3 and set(succ) == set(verts), "ORIENTATION_INCOHERENT", f"face {f} is not cyclically oriented")
        a = min(succ)
        require(succ[succ[succ[a]]] == a, "ORIENTATION_INCOHERENT", f"face {f} darts do not close")
        faces.append((a, succ[a], succ[succ[a]]))
    flat = [d for ds in darts.values() for d in ds]
    require(len(flat) == len(set(flat)) == 2 * len(inc.edge_vertices), "ORIENTATION_INCOHERENT", "a directed edge repeats")
    return tuple(sorted(faces))


def classify(inc: Incidence) -> str:
    """Run the checks in a fixed order and raise the first failure; CARRIER otherwise."""
    pairs = inc.edge_vertices
    require(all(len(p) == 2 for p in pairs) and len(set(pairs)) == len(pairs), "COLLAPSED_EDGES", "coset edges repeat a vertex pair")
    adj = adjacency(inc)
    triangles = {
        frozenset((a, b, c))
        for a in range(inc.n_vertices)
        for b in adj[a]
        for c in adj[a]
        if a < b < c and c in adj[b]
    }
    require(set(inc.face_vertices) == triangles, "FACES_NOT_TRIANGLES", "coset faces are not the triangles of the coset edge graph")
    for e, faces in enumerate(inc.edge_faces):
        for f in faces:
            require(
                pairs[e] <= inc.face_vertices[f],
                "FACE_BOUNDARY_MISMATCH",
                f"coset edge {sorted(pairs[e])} is coset-incident to face {sorted(inc.face_vertices[f])}",
            )
    for f, verts in enumerate(inc.face_vertices):
        sides = {pairs[e] for e, faces in enumerate(inc.edge_faces) if f in faces}
        require(sides == {frozenset(s) for s in itertools.combinations(sorted(verts), 2)}, "FACE_SIDES", f"face {sorted(verts)} has wrong sides")
    for v in range(inc.n_vertices):
        star_edges = [e for e, pair in enumerate(pairs) if v in pair]
        star_faces = {f for f, verts in enumerate(inc.face_vertices) if v in verts}
        seen, stack = {("e", star_edges[0])}, [("e", star_edges[0])]
        while stack:
            kind, x = stack.pop()
            if kind == "e":
                step = [("f", f) for f in inc.edge_faces[x]]
            else:
                step = [("e", e) for e in star_edges if x in inc.edge_faces[e]]
            for y in step:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        require(len(seen) == len(star_edges) + len(star_faces), "VERTEX_LINK", f"the link of vertex {v} is not one cycle")
    oriented_faces(inc)
    return "CARRIER"


def verdict(inc: Incidence) -> str:
    try:
        return classify(inc)
    except CertificateError as exc:
        return exc.code


def triangle_generators(quotient: Quotient, placement: Placement) -> list[tuple[Matrix, Matrix, Matrix]]:
    """Generators x in C2, y in C3, z in C5 with xyz = 1."""
    H5, H3, H2 = placement
    e = quotient.identity
    return [
        (x, y, z)
        for x in H2 if x != e
        for y in H3 if y != e
        for z in H5 if z != e
        if quotient.mul(x, quotient.mul(y, z)) == e
    ]


def normalizes(quotient: Quotient, g: Matrix, H: frozenset) -> bool:
    return frozenset(quotient.mul(quotient.mul(g, h), quotient.inv(g)) for h in H) == H


def classify_placements(quotient: Quotient, subs: Mapping[int, list[frozenset]]) -> tuple[dict[Placement, str], dict[str, Any]]:
    verdicts: dict[Placement, str] = {}
    local: set[tuple] = set()
    splits: set[tuple] = set()
    carrier_iff_triangle = unique_triangle = collapse_iff_normalizer = True
    e = quotient.identity
    for H5 in subs[5]:
        for H3 in subs[3]:
            split: Counter = Counter()
            for H2 in subs[2]:
                placement = (H5, H3, H2)
                inc = build_incidence(quotient, placement)
                local.add(local_counts(inc))
                v = verdict(inc)
                verdicts[placement] = v
                split[v] += 1
                generators = triangle_generators(quotient, placement)
                carrier_iff_triangle &= (v == "CARRIER") == bool(generators)
                if v == "CARRIER":
                    unique_triangle &= len(generators) == 1
                t = next(x for x in H2 if x != e)
                collapse_iff_normalizer &= (v == "COLLAPSED_EDGES") == normalizes(quotient, t, H5)
            splits.add(tuple(sorted(split.items())))
    counts = dict(Counter(verdicts.values()))
    expected = {"CARRIER": 120, "FACE_BOUNDARY_MISMATCH": 180, "FACES_NOT_TRIANGLES": 300, "COLLAPSED_EDGES": 300}
    require(counts == expected, "PLACEMENT_CLASSES", f"placement classes {counts}")
    require(local == {((5,), (3,), (5,), (2,), (2,), (3,))}, "LOCAL_COUNTS", f"local counts {local}")
    require(len(splits) == 1, "PAIR_SPLIT", f"the (C5, C3) pairs split the involutions differently: {splits}")
    require(carrier_iff_triangle, "TRIANGLE_RELATION", "carrier placements differ from the xyz = 1 placements")
    require(unique_triangle, "TRIANGLE_RELATION", "a carrier placement admits more than one generator choice with xyz = 1")
    require(collapse_iff_normalizer, "COLLAPSE_CRITERION", "collapsed placements differ from C2 inside N(C5)")
    summary = {
        "placement_classes": counts,
        "local_counts_placement_independent": {
            "faces_per_vertex": 5, "vertices_per_face": 3, "edges_per_vertex": 5,
            "vertices_per_edge": 2, "faces_per_edge": 2, "edges_per_face": 3,
        },
        "per_C5_C3_pair_split_of_the_15_involutions": dict(next(iter(splits))),
        "characterizations": {
            "carrier_iff_generators_with_xyz_equal_1": True,
            "xyz_equal_1_generator_choice_unique_on_carriers": True,
            "collapsed_iff_C2_inside_normalizer_of_C5": True,
        },
    }
    return verdicts, summary


def gl2_conjugators(quotient: Quotient) -> tuple[list[Callable[[Matrix], Matrix]], int]:
    """Conjugation by GL(2, F_p) on G, with the number of distinct automorphisms it induces."""
    p = quotient.source.p
    maps, images = [], set()
    for M in itertools.product(range(p), repeat=4):
        det = (M[0] * M[3] - M[1] * M[2]) % p
        if det == 0:
            continue
        s = pow(det, p - 2, p)
        Minv = ((M[3] * s) % p, (-M[1] * s) % p, (-M[2] * s) % p, (M[0] * s) % p)

        def conj(h: Matrix, M: Matrix = M, Minv: Matrix = Minv) -> Matrix:
            return quotient.canon(quotient.source.mul(quotient.source.mul(M, h), Minv))

        maps.append(conj)
        images.add(tuple(conj(g) for g in quotient.elements))
    return maps, len(images)


def orbit_sizes(placements: Sequence[Placement], conjugators: Sequence[Callable[[Matrix], Matrix]]) -> list[int]:
    remaining = set(placements)
    sizes = []
    for t in sorted(placements, key=placement_key):
        if t not in remaining:
            continue
        orbit = {tuple(frozenset(c(h) for h in H) for H in t) for c in conjugators}
        require(orbit <= set(placements), "ORBIT_CLOSURE", "a conjugate of a compatible placement is incompatible")
        remaining -= orbit
        sizes.append(len(orbit))
    return sorted(sizes)


def presentation_placements(quotient: Quotient) -> tuple[int, set[Placement]]:
    """Solutions of a^2 = b^3 = c^5 = abc (= -I) in the source, mapped to placements."""
    source = quotient.source
    minus = source.minus_identity
    solutions = []
    for a in source.elements:
        if source.mul(a, a) != minus:
            continue
        for b in source.elements:
            if source.mul(b, source.mul(b, b)) != minus:
                continue
            c = source.mul(source.inv(source.mul(a, b)), minus)
            c5 = c
            for _ in range(4):
                c5 = source.mul(c, c5)
            if c5 == minus:
                solutions.append((a, b, c))
    image = {
        tuple(frozenset(quotient.canon(g) for g in cyclic_subgroup(source.mul, source.identity, w)) for w in (c, b, a))
        for a, b, c in solutions
    }
    return len(solutions), image


# ---------------------------------------------------------------------------
# The committed A1 packet, read from the Lean sources
# ---------------------------------------------------------------------------


COMMITTED_ANTIPODE_DECLARATION = "def antipode (i : Fin 12) : Fin 12 := 11 - i"


@dataclass(frozen=True)
class Fixture:
    adjacency: tuple[frozenset, ...]
    port_action: frozenset
    oriented_faces: tuple[tuple[int, int, int], ...]
    antipode: tuple[int, ...] = tuple(11 - i for i in range(12))


def read_lean(module: str) -> str:
    try:
        return (REPO_ROOT / module).read_text(encoding="utf-8")
    except OSError as exc:
        raise CertificateError("LEAN_FIXTURE_READ", f"cannot read {module}: {exc}") from exc


def lean_block(text: str, marker: str, terminator: str, module: str) -> str:
    start = text.find(marker)
    require(start >= 0, "LEAN_FIXTURE_PARSE", f"{module}: `{marker}` is missing")
    end = text.find(terminator, start + len(marker))
    require(end >= 0, "LEAN_FIXTURE_PARSE", f"{module}: `{marker}` block is unterminated")
    return text[start + len(marker) : end + len(terminator)]


def load_fixture() -> Fixture:
    module = FIXTURE_DECLARATIONS["adjacency"][0]
    text = read_lean(module)
    block = lean_block(text, "def neighbors : Fin 12 → List (Fin 12)", "\n\n", module)
    rows = {int(k): tuple(int(t) for t in v.split(",")) for k, v in re.findall(r"\|\s*(\d+)\s*=>\s*\[([0-9,\s]+)\]", block)}
    require(sorted(rows) == list(range(12)) and all(len(r) == 5 for r in rows.values()), "LEAN_FIXTURE_PARSE", f"{module}: expected twelve five-entry rows")
    adj = tuple(frozenset(rows[i]) for i in range(12))
    require(all(i in adj[j] for i in range(12) for j in adj[i]), "LEAN_FIXTURE_PARSE", f"{module}: the neighbor table is not symmetric")
    require(COMMITTED_ANTIPODE_DECLARATION in text, "LEAN_FIXTURE_PARSE", f"{module}: the committed antipode declaration is missing")
    antipode = tuple(11 - i for i in range(12))
    require(distance_three_partner(adj) == list(antipode), "LEAN_FIXTURE_PARSE", f"{module}: the committed antipode is not the distance-three partner")

    module = FIXTURE_DECLARATIONS["port_action"][0]
    block = lean_block(read_lean(module), "def perms : List (List Nat) := [", "]]", module)
    perms = [tuple(int(t) for t in row.split(",")) for row in re.findall(r"\[([0-9,\s]+)\]", block)]
    require(
        len(perms) == 60 and len(set(perms)) == 60 and all(sorted(p) == list(range(12)) for p in perms),
        "LEAN_FIXTURE_PARSE",
        f"{module}: expected sixty distinct permutations of twelve ports",
    )

    module = FIXTURE_DECLARATIONS["oriented_faces"][0]
    block = lean_block(read_lean(module), "def orientedFaces : List (Fin 12 × Fin 12 × Fin 12) :=", "]", module)
    faces = tuple(tuple(int(t) for t in m) for m in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", block))
    require(len(faces) == 20, "LEAN_FIXTURE_PARSE", f"{module}: expected twenty oriented faces, got {len(faces)}")
    require(
        all(b in adj[a] and c in adj[b] and a in adj[c] for a, b, c in faces),
        "LEAN_FIXTURE_PARSE",
        f"{module}: an oriented face is not a triangle of the committed adjacency",
    )
    return Fixture(adj, frozenset(perms), faces, antipode)


def fixture_hashes(fixture: Fixture) -> dict[str, str]:
    return {
        "adjacency_sha256": sha256_json([sorted(row) for row in fixture.adjacency]),
        "port_action_sha256": sha256_json(sorted(list(p) for p in fixture.port_action)),
        "oriented_faces_sha256": sha256_json([list(f) for f in fixture.oriented_faces]),
    }


def breadth_first_order(adj: Sequence[frozenset], root: int = 0) -> list[int]:
    order, seen = [root], {root}
    for v in order:
        for w in sorted(adj[v]):
            if w not in seen:
                seen.add(w)
                order.append(w)
    return order + [v for v in range(len(adj)) if v not in seen]


def graph_isomorphisms(source: Sequence[frozenset], target: Sequence[frozenset]) -> list[tuple[int, ...]]:
    """Every bijection phi with u ~ v iff phi(u) ~ phi(v), as the tuple (phi(0), ..., phi(n-1))."""
    n = len(source)
    order = breadth_first_order(source)
    phi: dict[int, int] = {}
    used: set[int] = set()
    found: list[tuple[int, ...]] = []

    def extend(i: int) -> None:
        if i == n:
            found.append(tuple(phi[v] for v in range(n)))
            return
        v = order[i]
        for t in range(n):
            if t in used or len(target[t]) != len(source[v]):
                continue
            if all((phi[u] in target[t]) == (u in source[v]) for u in phi):
                phi[v] = t
                used.add(t)
                extend(i + 1)
                del phi[v]
                used.discard(t)

    extend(0)
    return found


def relabel_faces(faces: Sequence[tuple[int, int, int]], phi: Sequence[int]) -> tuple[tuple[int, int, int], ...]:
    out = []
    for face in faces:
        a, b, c = (phi[x] for x in face)
        out.append(min((a, b, c), (b, c, a), (c, a, b)))
    return tuple(sorted(out))


def distance_three_partner(adj: Sequence[frozenset]) -> list[int]:
    partner = []
    for start in range(len(adj)):
        dist, frontier = {start: 0}, [start]
        while frontier:
            step = []
            for u in frontier:
                for w in adj[u]:
                    if w not in dist:
                        dist[w] = dist[u] + 1
                        step.append(w)
            frontier = step
        far = [v for v, d in dist.items() if d == 3]
        require(len(far) == 1, "ANTIPODE_MISMATCH", f"vertex {start} has {len(far)} vertices at distance three")
        partner.append(far[0])
    return partner


def relabel(quotient: Quotient, placement: Placement, inc: Incidence, fixture: Fixture) -> tuple[tuple[int, ...], dict[str, int]]:
    """Explicit relabelling of one placement's coset geometry into the committed packet."""
    adj = adjacency(inc)
    isos = graph_isomorphisms(adj, fixture.adjacency)
    require(bool(isos), "RELABEL_ADJACENCY", "the coset edge graph is not isomorphic to the committed adjacency")
    faces = oriented_faces(inc)
    committed = relabel_faces(fixture.oriented_faces, range(12))
    matching = [phi for phi in isos if relabel_faces(faces, phi) == committed]
    require(bool(matching), "RELABEL_ORIENTATION", "no relabelling carries the coset orientation onto the committed oriented faces")
    reversed_packet = relabel_faces([(a, c, b) for a, b, c in fixture.oriented_faces], range(12))
    reversing = {phi for phi in isos if relabel_faces(faces, phi) == reversed_packet}
    require(
        reversing == {tuple(fixture.antipode[port] for port in phi) for phi in matching} and len(reversing) + len(matching) == len(isos),
        "ORIENTATION_ANTIPODE",
        "the relabellings onto the reversed packet are not the antipode compositions of the matching ones",
    )
    phi = matching[0]
    iV, reps = coset_index(quotient, placement[0])
    port_to_vertex = {port: v for v, port in enumerate(phi)}
    transported = frozenset(
        tuple(phi[iV[quotient.mul(g, reps[port_to_vertex[k]])]] for k in range(12)) for g in quotient.elements
    )
    require(transported == fixture.port_action, "PORT_ACTION_MISMATCH", "the transported action differs from the committed port rotations")
    partner = distance_three_partner(adj)
    require(all(phi[partner[v]] == fixture.antipode[phi[v]] for v in range(12)), "ANTIPODE_MISMATCH", "the distance-three partner is not the committed antipode")
    return phi, {"graph_isomorphisms": len(isos), "orientation_matching_relabellings": len(matching)}


# ---------------------------------------------------------------------------
# Receipt
# ---------------------------------------------------------------------------

CHOICE_ACCOUNTING = {
    "subgroup_orders": (
        "The orders five, three, and two fix the coset sizes 12, 20, and 30 and every local "
        "incidence count, identically for all 900 placements. They do not determine the global incidence."
    ),
    "relative_placement": (
        "Exactly 120 of the 900 placements close into the committed carrier. They are the placements "
        "with generators x, y, z of orders 2, 3, 5 satisfying xyz = 1; they form two orbits of 60 under "
        "inner automorphisms and one orbit under all automorphisms of the source. The placement is "
        "therefore unique up to automorphisms of the supplied group, and the presentation "
        "a^2 = b^3 = c^5 = abc supplies a compatible placement with no further choice."
    ),
    "duality": (
        "Twelve ports select the stabilizer of order five: G/C5 has twelve points, G/C3 has twenty, "
        "and G/C2 has thirty."
    ),
    "orientation": (
        "The construction outputs a coherent orientation, read from the chambers whose three cosets "
        "share one representative. For every compatible placement, 60 of the 120 graph relabellings "
        "carry it onto the committed orientedFaces, and a fixture with every face reversed is matched "
        "through the antipode."
    ),
    "residual_premise": (
        "The supplied group. The certificate takes SL(2, F_5) as given and makes no selection among "
        "candidate source groups."
    ),
}

CLAIM_BOUNDARY = (
    "Conditional on the supplied group SL(2, F_5), taken as the exact finite representative of the "
    "binary icosahedral group 2I; that identification is classical and is cited without formalization. "
    "The reconstruction factors through the central quotient, so it cannot distinguish the source from "
    "its quotient. It does not select the source group, identify ports with physical objects, or alter "
    "the declared status of the A1 boundary packet."
)


@dataclass
class Construction:
    quotient: Quotient
    subs: dict[int, list[frozenset]]
    verdicts: dict[Placement, str]
    carriers: list[Placement]
    fixture: Fixture


def construct(manifest: Mapping[str, Any]) -> tuple[Construction, dict[str, Any]]:
    source = build_source(manifest)
    quotient, quotient_summary = central_quotient(source)
    subs, subgroup_summary = subgroup_data(quotient)
    verdicts, classification = classify_placements(quotient, subs)
    carriers = sorted((t for t, v in verdicts.items() if v == "CARRIER"), key=placement_key)
    inner = [lambda h, g=g: quotient.mul(quotient.mul(g, h), quotient.inv(g)) for g in quotient.elements]
    automorphisms, automorphism_count = gl2_conjugators(quotient)
    require(automorphism_count == 120, "AUTOMORPHISM_GROUP", f"GL(2, F_5) induces {automorphism_count} automorphisms")
    inner_sizes = orbit_sizes(carriers, inner)
    all_sizes = orbit_sizes(carriers, automorphisms)
    require(inner_sizes == [60, 60] and all_sizes == [120], "ORBITS", f"orbit sizes {inner_sizes}, {all_sizes}")
    solutions, image = presentation_placements(quotient)
    require(solutions == 120 and image == set(carriers), "PRESENTATION", f"{solutions} solutions, image matches: {image == set(carriers)}")
    fixture = load_fixture()
    summary = {
        "source": {
            "group": "SL(2, F_5)",
            "order": 120,
            "centre": [list(source.identity), list(source.minus_identity)],
            "involutions": 1,
            "perfect": True,
            "element_order_profile": {str(k): v for k, v in sorted(BINARY_ICOSAHEDRAL_ORDER_PROFILE.items())},
        },
        "central_quotient": quotient_summary,
        "subgroups": subgroup_summary,
        "classification": classification,
        "orbits_of_compatible_placements": {
            "inner_automorphism_group_order": 60,
            "automorphism_group_order": automorphism_count,
            "orbit_sizes_under_inner_automorphisms": inner_sizes,
            "orbit_sizes_under_all_automorphisms": all_sizes,
        },
        "presentation": {
            "relation": "a^2 = b^3 = c^5 = abc",
            "solutions_in_source": solutions,
            "bijective_with_compatible_placements": True,
        },
    }
    return Construction(quotient, subs, verdicts, carriers, fixture), summary


def certificate_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    built, summary = construct(manifest)
    quotient, fixture = built.quotient, built.fixture
    iso_counts, match_counts = set(), set()
    canonical = None
    for t in built.carriers:
        inc = build_incidence(quotient, t)
        phi, stats = relabel(quotient, t, inc, fixture)
        iso_counts.add(stats["graph_isomorphisms"])
        match_counts.add(stats["orientation_matching_relabellings"])
        if canonical is None:
            canonical = (t, inc, phi)
    require(iso_counts == {120} and match_counts == {60}, "RELABEL_COUNTS", f"relabelling counts {iso_counts}, {match_counts}")
    t, inc, phi = canonical
    ((x, y, z),) = triangle_generators(quotient, t)
    faces = relabel_faces(oriented_faces(inc), phi)
    return {
        "schema": RECEIPT_SCHEMA,
        "generated_by": GENERATED_BY,
        "manifest_sha256": sha256_json(manifest),
        **summary,
        "comparison_fixture": {
            "modules": {role: {"module": m, "declaration": d} for role, (m, d) in FIXTURE_DECLARATIONS.items()},
            **fixture_hashes(fixture),
        },
        "reconstruction": {
            "compatible_placements_relabelled": len(built.carriers),
            "graph_relabellings_per_placement": 120,
            "orientation_matching_relabellings_per_placement": 60,
            "oriented_faces_equal_committed": True,
            "transported_action_equals_committed_port_rotations": True,
            "distance_three_partner_equals_committed_antipode": True,
            "relabellings_onto_the_reversed_packet_are_antipode_compositions": True,
        },
        "canonical_placement": {
            "generators_with_xyz_equal_1": {"x": list(x), "y": list(y), "z": list(z)},
            "relabelling": [
                {"port": phi[v], "vertex_coset_least_element": list(inc.vertex_representatives[v])}
                for v in sorted(range(12), key=lambda v: phi[v])
            ],
            "oriented_faces_in_port_labels": [list(f) for f in faces],
        },
        "choice_accounting": CHOICE_ACCOUNTING,
        "claim_boundary": CLAIM_BOUNDARY,
    }


def verify_receipt(manifest: Mapping[str, Any], receipt: Mapping[str, Any]) -> None:
    require(receipt == certificate_payload(manifest), "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def mutated(manifest: Mapping[str, Any], edit: Callable[[dict[str, Any]], None]) -> dict[str, Any]:
    copy = json.loads(json.dumps(manifest))
    edit(copy)
    return copy


def first_placement(built: Construction, code: str) -> Placement:
    return min((t for t, v in built.verdicts.items() if v == code), key=placement_key)


def with_edge_faces(inc: Incidence, edge: int, faces: frozenset) -> Incidence:
    ef = list(inc.edge_faces)
    ef[edge] = faces
    return Incidence(inc.edge_vertices, inc.face_vertices, tuple(ef), inc.chambers, inc.vertex_representatives)


def degree_preserving_swap(adj: Sequence[frozenset]) -> tuple[tuple[frozenset, ...], tuple]:
    edges = sorted(tuple(sorted((a, b))) for a in range(12) for b in adj[a] if a < b)
    for (a, b), (c, d) in itertools.combinations(edges, 2):
        if len({a, b, c, d}) == 4 and d not in adj[a] and b not in adj[c]:
            rows = [set(r) for r in adj]
            for u, w in ((a, b), (c, d)):
                rows[u].discard(w)
                rows[w].discard(u)
            for u, w in ((a, d), (c, b)):
                rows[u].add(w)
                rows[w].add(u)
            return tuple(frozenset(r) for r in rows), ((a, b), (c, d))
    raise CertificateError("CONTROL_SETUP", "no degree-preserving edge swap exists")


def negative_control_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    built, _ = construct(manifest)
    quotient, fixture = built.quotient, built.fixture
    parent = built.carriers[0]
    parent_inc = build_incidence(quotient, parent)
    require(classify(parent_inc) == "CARRIER", "CONTROL_SETUP", "the parent placement is not a carrier")
    relabel(quotient, parent, parent_inc, fixture)
    edge0_true = frozenset(f for f, verts in enumerate(parent_inc.face_vertices) if parent_inc.edge_vertices[0] <= verts)
    edge0_false = next(f for f, verts in enumerate(parent_inc.face_vertices) if not parent_inc.edge_vertices[0] <= verts)
    swapped_adj, swap = degree_preserving_swap(fixture.adjacency)
    bad_row = list(sorted(fixture.port_action)[1])
    bad_row[0], bad_row[1] = bad_row[1], bad_row[0]
    a, b, c = fixture.oriented_faces[0]

    def run(thunk: Callable[[], Any]) -> str:
        try:
            thunk()
        except CertificateError as exc:
            return exc.code
        return "ACCEPTED"

    cases: list[tuple[str, Callable[[], Any], str]] = [
        ("declared_prime_3", lambda: build_source(mutated(manifest, lambda m: m["source"].update(prime=3))), "SOURCE_ORDER"),
        ("declared_prime_7", lambda: build_source(mutated(manifest, lambda m: m["source"].update(prime=7))), "SOURCE_ORDER"),
        ("declared_family_GL2", lambda: build_source(mutated(manifest, lambda m: m["source"].update(group="GL2"))), "SOURCE_GROUP"),
        ("source_names_target_structure", lambda: build_source(mutated(manifest, lambda m: m["source"].update(group="SL2 icosahedral"))), "SOURCE_FIREWALL"),
        ("source_extra_field", lambda: build_source(mutated(manifest, lambda m: m["source"].update(order=120))), "SCHEMA_FIELDS"),
        (
            "fixture_declaration_swapped",
            lambda: validate_manifest(mutated(manifest, lambda m: m["comparison_fixture"]["adjacency"].update(declaration="perms"))),
            "FIXTURE_DECLARATION",
        ),
        ("placement_collapsed_edges", lambda: classify(build_incidence(quotient, first_placement(built, "COLLAPSED_EDGES"))), "COLLAPSED_EDGES"),
        ("placement_faces_not_triangles", lambda: classify(build_incidence(quotient, first_placement(built, "FACES_NOT_TRIANGLES"))), "FACES_NOT_TRIANGLES"),
        ("placement_face_boundary_mismatch", lambda: classify(build_incidence(quotient, first_placement(built, "FACE_BOUNDARY_MISMATCH"))), "FACE_BOUNDARY_MISMATCH"),
        (
            "carrier_edge_face_incidence_corrupted",
            lambda: classify(with_edge_faces(parent_inc, 0, frozenset([min(edge0_true), edge0_false]))),
            "FACE_BOUNDARY_MISMATCH",
        ),
        ("carrier_edge_dropped_from_a_face", lambda: classify(with_edge_faces(parent_inc, 0, frozenset([min(edge0_true)]))), "FACE_SIDES"),
        (
            "carrier_chamber_moved_to_another_face",
            lambda: classify(Incidence(
                parent_inc.edge_vertices, parent_inc.face_vertices, parent_inc.edge_faces,
                ((parent_inc.chambers[0][0], parent_inc.chambers[0][1], (parent_inc.chambers[0][2] + 1) % 20),) + parent_inc.chambers[1:],
                parent_inc.vertex_representatives,
            )),
            "ORIENTATION_INCOHERENT",
        ),
        (
            "committed_adjacency_double_edge_swap",
            lambda: relabel(quotient, parent, parent_inc, Fixture(swapped_adj, fixture.port_action, fixture.oriented_faces)),
            "RELABEL_ADJACENCY",
        ),
        (
            "committed_face_reversed",
            lambda: relabel(quotient, parent, parent_inc, Fixture(fixture.adjacency, fixture.port_action, ((a, c, b),) + fixture.oriented_faces[1:])),
            "RELABEL_ORIENTATION",
        ),
        (
            "committed_port_row_perturbed",
            lambda: relabel(quotient, parent, parent_inc, Fixture(fixture.adjacency, (fixture.port_action - {sorted(fixture.port_action)[1]}) | {tuple(bad_row)}, fixture.oriented_faces)),
            "PORT_ACTION_MISMATCH",
        ),
    ]
    results = []
    for name, thunk, expected in cases:
        actual = run(thunk)
        require(actual == expected, "NEGATIVE_CONTROL_FAILED", f"{name}: expected {expected}, got {actual}")
        results.append({"name": name, "expected_error": expected, "actual_error": actual, "passed": True})

    reversed_fixture = Fixture(fixture.adjacency, fixture.port_action, tuple((p, r, q) for p, q, r in fixture.oriented_faces))
    reversed_ok = run(lambda: relabel(quotient, parent, parent_inc, reversed_fixture)) == "ACCEPTED"
    require(reversed_ok, "CONTROL_SETUP", "a fixture with every face reversed is not matched")

    mismatch = first_placement(built, "FACE_BOUNDARY_MISMATCH")
    minc = build_incidence(quotient, mismatch)
    phi = graph_isomorphisms(adjacency(minc), fixture.adjacency)[0]
    bad_edge = next(e for e, fs in enumerate(minc.edge_faces) if any(not minc.edge_vertices[e] <= minc.face_vertices[f] for f in fs))
    ports = lambda verts: sorted(phi[v] for v in verts)  # noqa: E731
    face_sets_agree = {frozenset(phi[v] for v in f) for f in minc.face_vertices} == {frozenset(f) for f in fixture.oriented_faces}

    collapsed = build_incidence(quotient, first_placement(built, "COLLAPSED_EDGES"))
    multiplicities = sorted(Counter(collapsed.edge_vertices).values())
    require(multiplicities == [5] * 6, "CONTROL_SETUP", f"collapsed edge multiplicities {multiplicities}")
    nontriangles = build_incidence(quotient, first_placement(built, "FACES_NOT_TRIANGLES"))
    nadj = adjacency(nontriangles)
    triangle_faces = sum(1 for f in nontriangles.face_vertices if all(v in nadj[u] for u, v in itertools.combinations(sorted(f), 2)))

    counts = built_counts = Counter(built.verdicts.values())
    return {
        "schema": NEGATIVE_SCHEMA,
        "discussion": 769,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "right_graph_and_faces_wrong_boundaries": {
                "vertex_graph_isomorphic_to_committed": True,
                "face_vertex_sets_equal_committed_under_that_relabelling": face_sets_agree,
                "edge_in_port_labels": ports(minc.edge_vertices[bad_edge]),
                "coset_incident_faces_in_port_labels": sorted(ports(minc.face_vertices[f]) for f in minc.edge_faces[bad_edge]),
                "faces_containing_that_edge": sorted(ports(v) for v in minc.face_vertices if minc.edge_vertices[bad_edge] <= v),
            },
            "collapsed_edges": {
                "distinct_vertex_pairs": len(set(collapsed.edge_vertices)),
                "edges_per_pair": 5,
                "mechanism": "the involution of C2 normalizes C5, so every coset edge joins gC5 to gtC5",
            },
            "faces_not_triangles": {"coset_faces_that_are_triangles_of_the_coset_graph": triangle_faces},
            "every_committed_face_reversed": {"relabelling_exists": reversed_ok, "absorbed_by": "the antipode i -> 11 - i"},
        },
        "classified_larger_families": {
            "subgroup_placements": 900,
            "placement_classes": dict(sorted(counts.items())),
            "placements_passing_vertex_level_checks_only": built_counts["CARRIER"] + built_counts["FACE_BOUNDARY_MISMATCH"],
            "placements_passing_all_three_incidences": built_counts["CARRIER"],
            "degree_preserving_swap_used": [list(e) for e in swap],
        },
    }


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def default_paths() -> tuple[Path, Path, Path]:
    here = Path(__file__).resolve().parent
    return (
        here / "manifests" / "coset_carrier_reference.json",
        here / "receipts" / "coset_carrier_reference.receipt.json",
        here / "negative_controls" / "coset_carrier_negative_controls.json",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify", help="create the deterministic exact receipt")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify", help="recompute and compare a receipt")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls", help="run and write the finite control bundle")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    everything = sub.add_parser("all", help="regenerate receipt and negative controls at the default paths")
    everything.add_argument("--manifest", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "certify":
        receipt = certificate_payload(load_json(args.manifest))
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output), "sha256": sha256_json(receipt)}, indent=2))
    elif args.command == "verify":
        verify_receipt(load_json(args.manifest), load_json(args.receipt))
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        write_json(args.output, negative_control_payload(load_json(args.manifest)))
        print(json.dumps({"status": "PASS", "negative_controls": str(args.output)}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest = load_json(args.manifest or default_manifest)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(json.dumps({"status": "PASS", "receipt": str(default_receipt), "negative_controls": str(default_negative)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

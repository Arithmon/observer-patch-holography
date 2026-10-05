"""Finite support complexes, including faces above dimension two.

Combinatorial predicates check the supplied complex, not its two-skeleton.
This module contains no quantum channel or physical geometry certificate.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from numbers import Integral


def _label(x: object) -> bool:
    return isinstance(x, Integral) and not isinstance(x, bool)


@dataclass
class IncidenceComplex:
    vertices: list[int]
    edges: set[frozenset]
    triangles: set[frozenset]
    higher_simplices: set[frozenset] = field(default_factory=set)


def validate_complex(k: IncidenceComplex) -> None:
    """Reject malformed face data, including omitted boundaries and aliases."""
    if not isinstance(k, IncidenceComplex):
        raise ValueError("an IncidenceComplex is required")
    if (not isinstance(k.vertices, list) or not all(_label(v) for v in k.vertices)
            or len(set(k.vertices)) != len(k.vertices)):
        raise ValueError("vertices must be distinct integer labels")
    vertices = set(k.vertices)
    for faces, size in ((k.edges, 2), (k.triangles, 3), (k.higher_simplices, None)):
        if not isinstance(faces, set):
            raise ValueError("face collections must be sets")
        for face in faces:
            if (not isinstance(face, frozenset) or not all(_label(v) for v in face)
                    or not face <= vertices
                    or (len(face) != size if size else len(face) < 4)):
                raise ValueError("invalid simplex or undeclared vertex")
    all_faces = ({frozenset((v,)) for v in vertices}
                 | k.edges | k.triangles | k.higher_simplices)
    for face in k.edges | k.triangles | k.higher_simplices:
        if any(face - {v} not in all_faces for v in face):
            raise ValueError("complex is not closed under taking faces")


def simplices(k: IncidenceComplex) -> set[frozenset]:
    validate_complex(k)
    return ({frozenset((v,)) for v in k.vertices}
            | k.edges | k.triangles | k.higher_simplices)


def validate_records(records: list[tuple[int, ...]]) -> None:
    """Validate support labels without allocating their simplicial closure."""
    if not isinstance(records, (list, tuple)):
        raise ValueError("records must be a finite sequence")
    for record in records:
        if (not isinstance(record, (list, tuple)) or not record
                or not all(_label(v) for v in record)
                or len(set(record)) != len(record)):
            raise ValueError("each record needs distinct integer patch labels")


def incidence_complex(records: list[tuple[int, ...]]) -> IncidenceComplex:
    """Every nonempty subset of a record's support is a simplex."""
    validate_records(records)
    faces: set[frozenset] = set()
    for record in records:
        for size in range(1, len(record) + 1):
            faces.update(frozenset(f) for f in itertools.combinations(record, size))
    return IncidenceComplex(
        sorted(next(iter(f)) for f in faces if len(f) == 1),
        {f for f in faces if len(f) == 2},
        {f for f in faces if len(f) == 3},
        {f for f in faces if len(f) >= 4},
    )


def euler_characteristic(k: IncidenceComplex) -> int:
    return sum((-1) ** (len(f) - 1) for f in simplices(k))


def dimension(k: IncidenceComplex) -> int:
    return max((len(f) - 1 for f in simplices(k)), default=-1)


def complexes_equal(left: IncidenceComplex, right: IncidenceComplex) -> bool:
    return simplices(left) == simplices(right)


def refinement_is_simplicial(fine: IncidenceComplex, coarse: IncidenceComplex,
                             projection: dict[int, int]) -> bool:
    """Check all faces and a total vertex map; no topology or CPTP inference."""
    fine_faces, coarse_faces = simplices(fine), simplices(coarse)
    if (not isinstance(projection, dict) or set(projection) != set(fine.vertices)
            or not all(_label(v) for v in projection)
            or not all(_label(v) and v in coarse.vertices for v in projection.values())):
        return False
    return all(frozenset(projection[v] for v in face) in coarse_faces
               for face in fine_faces)


def complexes_isomorphic_under(left: IncidenceComplex, right: IncidenceComplex,
                                permutation: dict[int, int]) -> bool:
    if not refinement_is_simplicial(left, right, permutation):
        return False
    if (len(set(permutation.values())) != len(left.vertices)
            or set(permutation.values()) != set(right.vertices)):
        return False
    return {frozenset(permutation[v] for v in f) for f in simplices(left)} == simplices(right)


def certify_midpoint_subdivision(coarse: IncidenceComplex, fine: IncidenceComplex,
                                 edge_midpoints: dict[frozenset, int]) -> bool:
    """Recognize the entire 1-to-4 subdivision, including isolated vertices.

    The certificate supplies one distinct fresh vertex per coarse edge. The
    realization sending it to that edge's midpoint is a PL homeomorphism.
    An endpoint-valued simplicial projection is a different map.
    """
    validate_complex(coarse)
    validate_complex(fine)
    if coarse.higher_simplices or fine.higher_simplices:
        return False
    if (not isinstance(edge_midpoints, dict) or set(edge_midpoints) != coarse.edges
            or not all(isinstance(e, frozenset) and all(_label(v) for v in e)
                       for e in edge_midpoints)
            or not all(_label(v) for v in edge_midpoints.values())):
        return False
    midpoints = set(edge_midpoints.values())
    if len(midpoints) != len(coarse.edges) or midpoints & set(coarse.vertices):
        return False
    records = [(v,) for v in coarse.vertices]
    for edge, midpoint in edge_midpoints.items():
        records.extend((v, midpoint) for v in edge)
    for triangle in coarse.triangles:
        a, b, c = sorted(triangle)
        ab, bc, ac = (edge_midpoints[frozenset(e)] for e in ((a, b), (b, c), (a, c)))
        records.extend(((a, ab, ac), (b, ab, bc), (c, ac, bc), (ab, ac, bc)))
    return complexes_equal(incidence_complex(records), fine)

#!/usr/bin/env python3
"""Exact icosahedral sharp-code and D-optimal tomography certificate.

The geometric checks needed for the Cohn--Kumar universal-optimality theorem
are verified symbolically:
  * exactly three inner products among distinct vertices;
  * spherical moments through degree five equal the uniform S^2 moments.
Because the configuration is antipodal, odd moments vanish.  Degree 2 and 4
are checked exactly, which is sufficient for the spherical 5-design claim.

External theorem: a sharp configuration with m distances and strength 2m-1
is universally optimal; for the icosahedron the strictly completely monotonic
minimum is unique up to O(3).

The internal D-optimal check proves the vector and quadrupole frame equalities,
the Seidel relation, and the unique five-cycle switching class.
"""
from __future__ import annotations

import itertools
import json
import sympy as sp

SQRT5 = sp.sqrt(5)
PHI = (1 + SQRT5) / 2
NORM2 = sp.simplify(1 + PHI**2)


class CertificateError(RuntimeError):
    """Fail-closed certificate error with a stable code."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")
        self.code = code
        self.detail = detail


def require(condition: bool, code: str, detail: str) -> None:
    if not condition:
        raise CertificateError(code, detail)


def vertices() -> list[sp.Matrix]:
    out: list[sp.Matrix] = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            out.append(sp.Matrix([0, s1, s2 * PHI]) / sp.sqrt(NORM2))
            out.append(sp.Matrix([s1, s2 * PHI, 0]) / sp.sqrt(NORM2))
            out.append(sp.Matrix([s1 * PHI, 0, s2]) / sp.sqrt(NORM2))
    return out


def axes() -> list[sp.Matrix]:
    n = sp.sqrt(NORM2)
    return [
        sp.Matrix(v) / n
        for v in [
            (0, 1, PHI),
            (0, 1, -PHI),
            (1, PHI, 0),
            (1, -PHI, 0),
            (PHI, 0, 1),
            (PHI, 0, -1),
        ]
    ]


def sym0_basis() -> list[sp.Matrix]:
    return [
        sp.diag(1, -1, 0) / sp.sqrt(2),
        sp.diag(1, 1, -2) / sp.sqrt(6),
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / sp.sqrt(2),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sp.sqrt(2),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / sp.sqrt(2),
    ]


def delta(i: int, j: int) -> int:
    return 1 if i == j else 0


def payload() -> dict:
    vs = vertices()
    require(len(vs) == 12, "VERTEX_COUNT", f"the orbit has {len(vs)} vertices, not twelve")
    norms = [sp.simplify(v.dot(v)) for v in vs]
    require(
        all(v == 1 for v in norms),
        "UNIT_NORMS",
        "a vertex does not lie on the unit sphere",
    )

    ips = sorted({sp.simplify(vs[i].dot(vs[j])) for i in range(12) for j in range(i + 1, 12)}, key=float)
    expected_ips = [-sp.Integer(1), -1 / SQRT5, 1 / SQRT5]
    require(
        len(ips) == len(expected_ips)
        and all(sp.simplify(a - b) == 0 for a, b in zip(ips, expected_ips, strict=True)),
        "DISTINCT_INNER_PRODUCTS",
        f"the distinct inner products are {ips}, not {expected_ips}",
    )

    mean = sp.simplify(sum(vs, sp.zeros(3, 1)) / 12)
    require(
        mean == sp.zeros(3, 1),
        "FIRST_MOMENT",
        "the vertex configuration is not centred",
    )

    m2 = sp.simplify(sum((v * v.T for v in vs), sp.zeros(3, 3)) / 12)
    require(
        m2 == sp.eye(3) / 3,
        "SECOND_MOMENT",
        "the second moment does not equal the uniform S^2 second moment",
    )

    fourth_bad = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    actual = sp.simplify(sum(v[i] * v[j] * v[k] * v[l] for v in vs) / 12)
                    expected = sp.Rational(1, 15) * (
                        delta(i, j) * delta(k, l)
                        + delta(i, k) * delta(j, l)
                        + delta(i, l) * delta(j, k)
                    )
                    if sp.simplify(actual - expected) != 0:
                        fourth_bad.append((i, j, k, l, actual, expected))
    require(
        not fourth_bad,
        "FOURTH_MOMENT",
        f"{len(fourth_bad)} fourth-moment components miss the uniform S^2 value",
    )

    # Antipodality makes every odd moment vanish exactly.
    antipodal = all(any(sp.simplify(v + w) == sp.zeros(3, 1) for w in vs) for v in vs)
    require(
        antipodal,
        "ANTIPODALITY",
        "the configuration is not antipodal, so odd moments need not vanish",
    )

    us = axes()
    ps = [sp.simplify(u * u.T) for u in us]
    qs = [sp.simplify(p - sp.eye(3) / 3) for p in ps]
    f1 = sp.simplify(sum(ps, sp.zeros(3)))
    basis = sym0_basis()
    qcoords = sp.Matrix(5, 6, lambda a, i: sp.simplify(sp.trace(basis[a] * qs[i])))
    f2 = sp.simplify(qcoords * qcoords.T)
    qgram = sp.simplify(qcoords.T * qcoords)
    require(
        f1 == 2 * sp.eye(3),
        "VECTOR_FISHER_FRAME",
        "the axis projectors do not sum to twice the identity",
    )
    require(
        f2 == sp.Rational(4, 5) * sp.eye(5),
        "QUADRUPOLE_FISHER_FRAME",
        "the quadrupole coordinates do not form a tight frame",
    )
    require(
        qgram == sp.Rational(4, 5) * sp.eye(6) - sp.Rational(2, 15) * sp.ones(6),
        "QUADRUPOLE_GRAM",
        "the quadrupole Gram matrix does not have the equiangular form",
    )

    dots = sp.Matrix(6, 6, lambda i, j: sp.simplify(us[i].dot(us[j])))
    seidel = sp.Matrix(
        6, 6, lambda i, j: 0 if i == j else sp.simplify(SQRT5 * dots[i, j])
    )
    require(
        seidel * seidel == 5 * sp.eye(6),
        "SEIDEL_RELATION",
        "the Seidel matrix does not satisfy S^2 = 5 I_6",
    )

    edges = list(itertools.combinations(range(5), 2))
    switched_solutions = 0
    solution_degree_profiles: set[tuple[int, ...]] = set()
    for mask in range(1 << len(edges)):
        s = sp.zeros(6)
        for j in range(1, 6):
            s[0, j] = s[j, 0] = 1
        for k, (a, b) in enumerate(edges):
            s[a + 1, b + 1] = s[b + 1, a + 1] = -1 if mask >> k & 1 else 1
        if s * s == 5 * sp.eye(6):
            switched_solutions += 1
            degrees = tuple(
                sum(int(s[a + 1, b + 1] == -1) for b in range(5) if a != b)
                for a in range(5)
            )
            solution_degree_profiles.add(degrees)
    require(
        switched_solutions == 12,
        "SEIDEL_SWITCHING_COUNT",
        f"{switched_solutions} switched Seidel solutions, not twelve",
    )
    require(
        solution_degree_profiles == {(2, 2, 2, 2, 2)},
        "SEIDEL_SWITCHING_CLASS",
        "the switched Seidel solutions are not a single five-cycle class",
    )

    return {
        "schema": "A5 icosahedral selection certificate v1",
        "n_vertices": 12,
        "unit_norms": True,
        "distinct_inner_products": [str(sp.simplify(x)) for x in ips],
        "number_of_distances": 3,
        "mean_zero": True,
        "second_moment": [[str(x) for x in row] for row in m2.tolist()],
        "fourth_moment_matches_uniform_S2": True,
        "odd_moments_through_degree_5_vanish_by_antipodality": True,
        "spherical_design_strength": 5,
        "sharp_condition": "m=3 distances and t=2m-1=5",
        "d_optimal_vector_fisher": [[str(x) for x in row] for row in f1.tolist()],
        "d_optimal_quadrupole_fisher": [[str(x) for x in row] for row in f2.tolist()],
        "d_optimal_determinants": [str(f1.det()), str(f2.det())],
        "equiangular_axis_squared_inner_product": "1/5",
        "seidel_relation": "S^2 = 5 I_6",
        "switched_seidel_solutions": switched_solutions,
        "switched_seidel_isomorphism_classes": 1,
        "normalized_negative_graph": "C5",
        "external_theorem": "Cohn--Kumar: sharp configurations are universally optimal; strict completely monotonic energy gives the icosahedron uniquely up to O(3)",
        "oph_consequence": "Either the D-optimal tomography objective or a source-derived strictly completely monotonic pair potential selects the icosahedral A5 orbit.",
        "open_physical_gate": "source production of one selector without an icosahedral template",
    }


if __name__ == "__main__":
    print(json.dumps(payload(), indent=2))

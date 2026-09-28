#!/usr/bin/env python3
"""Bounded alternative-carrier and response-law specificity ensemble.

Issue #644 asks how specific the registered structural hits are: run the
same structural probes on a declared bounded ensemble of alternative finite
carriers and publish the full hit table. The probes are source-blind
structural computations in exact integer and rational arithmetic:

* twelve ports; degree-five regularity; exactly four spectral bands
  (minimal-polynomial degree by exact rational linear algebra); a
  quadratic-irrational Galois band pair (rational-root stripping of the
  minimal polynomial); the (1, 3, 3, 5) band multiplicity pattern
  (exhaustion against exact power traces); a unique-antipode involutive
  automorphism; twenty triangles; and the declared response polynomial
  ``J = (A^3 - 4A^2 - 5A + 10I)/10`` squaring to the identity.

The response-law half enumerates the complete involutive-response family on
the icosahedral spectrum: a degree-three polynomial with ``J^2 = I`` takes a
sign at each of the four bands, giving sixteen laws by exact Lagrange
interpolation in the quadratic field; exactly the Galois-symmetric sign
patterns have rational coefficients, and the declared law is recovered
coefficient-exactly from its sign pattern.

The orientation stratum uses the one orientation-changing operation available
on a carrier with a unique antipode involution: identifying antipodes. For the
icosahedron the quotient is the hemi-icosahedron (edge graph K6, the minimal
triangulation of the real projective plane), the carrier obtained by dropping
orientation; it joins the menu as a member. Exact rational band projectors
give the antipode parity of every band, and the published compact-current
dimension argument (current dimension equal to the port count, one fixed line
of the transitive A5 action, centre of dimension at most one, inner A5 action)
is replayed on the twelve ports of the icosahedron and on the six ports of the
quotient.

The ensemble is a declared bounded menu, the score is a structural
specificity calibration, and nothing here is a physical forecast. No public
measurement is read and no comparison is permitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Any, Callable

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
OUTPUT_PATH = RUNTIME / "carrier_specificity_receipt.json"

SCHEMA = "oph.carrier_specificity_receipt.v2"
STATUS = "NAMED_MENU_FINGERPRINT__EXHAUSTIVE_CALIBRATION_OPEN"


class SpecificityError(ValueError):
    """The specificity ensemble refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SpecificityError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Declared carrier ensemble: explicit exact constructions
# ---------------------------------------------------------------------------


def _cycle_edges(names: list[str]) -> list[tuple[str, str]]:
    return [(names[i], names[(i + 1) % len(names)]) for i in range(len(names))]


def icosahedron() -> list[tuple[str, str]]:
    upper = [f"u{i}" for i in range(5)]
    lower = [f"l{i}" for i in range(5)]
    edges = []
    for i in range(5):
        edges.append(("n", upper[i]))
        edges.append(("s", lower[i]))
        edges.append((upper[i], upper[(i + 1) % 5]))
        edges.append((lower[i], lower[(i + 1) % 5]))
        edges.append((upper[i], lower[i]))
        edges.append((upper[i], lower[(i + 1) % 5]))
    return edges


def cuboctahedron() -> list[tuple[str, str]]:
    top = [f"t{i}" for i in range(4)]
    middle = [f"m{i}" for i in range(4)]
    bottom = [f"b{i}" for i in range(4)]
    edges = _cycle_edges(top) + _cycle_edges(bottom)
    for i in range(4):
        edges.append((top[i], middle[i]))
        edges.append((top[i], middle[(i - 1) % 4]))
        edges.append((bottom[i], middle[i]))
        edges.append((bottom[i], middle[(i - 1) % 4]))
    return edges


def truncated_tetrahedron() -> list[tuple[str, str]]:
    triangles = [[f"{corner}{i}" for i in range(3)] for corner in "abcd"]
    edges = []
    for triangle in triangles:
        edges.extend(_cycle_edges(triangle))
    edges.extend(
        [
            ("a0", "b0"),
            ("a1", "c0"),
            ("a2", "d0"),
            ("b1", "c2"),
            ("b2", "d1"),
            ("c1", "d2"),
        ]
    )
    return edges


def hexagonal_prism() -> list[tuple[str, str]]:
    top = [f"p{i}" for i in range(6)]
    bottom = [f"q{i}" for i in range(6)]
    edges = _cycle_edges(top) + _cycle_edges(bottom)
    edges.extend((top[i], bottom[i]) for i in range(6))
    return edges


def octahedron() -> list[tuple[str, str]]:
    ring = [f"r{i}" for i in range(4)]
    edges = _cycle_edges(ring)
    for pole in ("n", "s"):
        edges.extend((pole, vertex) for vertex in ring)
    return edges


def cube() -> list[tuple[str, str]]:
    top = [f"t{i}" for i in range(4)]
    bottom = [f"b{i}" for i in range(4)]
    edges = _cycle_edges(top) + _cycle_edges(bottom)
    edges.extend((top[i], bottom[i]) for i in range(4))
    return edges


def dodecahedron() -> list[tuple[str, str]]:
    outer = [f"o{i}" for i in range(5)]
    upper = [f"u{i}" for i in range(5)]
    lower = [f"l{i}" for i in range(5)]
    inner = [f"i{i}" for i in range(5)]
    edges = _cycle_edges(outer) + _cycle_edges(inner)
    for i in range(5):
        edges.append((outer[i], upper[i]))
        edges.append((upper[i], lower[i]))
        edges.append((upper[i], lower[(i - 1) % 5]))
        edges.append((lower[i], inner[i]))
    return edges


def tetrahedron() -> list[tuple[str, str]]:
    names = [f"v{i}" for i in range(4)]
    return [
        (names[i], names[j]) for i in range(4) for j in range(i + 1, 4)
    ]


def petersen() -> list[tuple[str, str]]:
    outer = [f"o{i}" for i in range(5)]
    inner = [f"i{i}" for i in range(5)]
    edges = _cycle_edges(outer)
    edges.extend((inner[i], inner[(i + 2) % 5]) for i in range(5))
    edges.extend((outer[i], inner[i]) for i in range(5))
    return edges


def complete_k12() -> list[tuple[str, str]]:
    names = [f"k{i}" for i in range(12)]
    return [
        (names[i], names[j]) for i in range(12) for j in range(i + 1, 12)
    ]


def cycle_c12() -> list[tuple[str, str]]:
    return _cycle_edges([f"c{i}" for i in range(12)])


def complete_bipartite_k66() -> list[tuple[str, str]]:
    left = [f"x{i}" for i in range(6)]
    right = [f"y{i}" for i in range(6)]
    return [(a, b) for a in left for b in right]


def antipodal_quotient(edges: list[tuple[str, str]], prefix: str) -> list[tuple[str, str]]:
    """Identify every vertex with its unique antipode and keep the edge image."""

    names = sorted({name for edge in edges for name in edge})
    index = {name: position for position, name in enumerate(names)}
    pairing = antipode_pairing(adjacency_matrix(edges))
    require(pairing is not None, "carrier has no unique antipode involution")
    classes = sorted({min(v, pairing[v]) for v in range(len(names))})
    label = {
        v: f"{prefix}{classes.index(min(v, pairing[v]))}" for v in range(len(names))
    }
    quotient = set()
    for left, right in edges:
        a, b = label[index[left]], label[index[right]]
        require(a != b, "antipodal vertices are adjacent")
        quotient.add((min(a, b), max(a, b)))
    return sorted(quotient)


def hemi_icosahedron() -> list[tuple[str, str]]:
    """Antipodal quotient of the icosahedron: six port classes, edge graph K6.

    Each of the fifteen quotient edges is the image of one antipodal pair of
    icosahedral edges. A5 acts on the six classes through its rotation action
    on the six fivefold axes, so the carrier keeps the icosahedral rotation
    symmetry and loses the orientation (the quotient surface is the real
    projective plane).
    """

    return antipodal_quotient(icosahedron(), "h")


ENSEMBLE: dict[str, Callable[[], list[tuple[str, str]]]] = {
    "icosahedron": icosahedron,
    "hemi_icosahedron": hemi_icosahedron,
    "cuboctahedron": cuboctahedron,
    "truncated_tetrahedron": truncated_tetrahedron,
    "hexagonal_prism": hexagonal_prism,
    "octahedron": octahedron,
    "cube": cube,
    "dodecahedron": dodecahedron,
    "tetrahedron": tetrahedron,
    "petersen": petersen,
    "complete_k12": complete_k12,
    "cycle_c12": cycle_c12,
    "complete_bipartite_k66": complete_bipartite_k66,
}


# ---------------------------------------------------------------------------
# Exact structural probes
# ---------------------------------------------------------------------------


def adjacency_matrix(edges: list[tuple[str, str]]) -> list[list[int]]:
    names = sorted({name for edge in edges for name in edge})
    index = {name: position for position, name in enumerate(names)}
    size = len(names)
    matrix = [[0] * size for _ in range(size)]
    for left, right in edges:
        require(left != right, "self loop in declared carrier")
        require(matrix[index[left]][index[right]] == 0, "duplicate edge")
        matrix[index[left]][index[right]] = 1
        matrix[index[right]][index[left]] = 1
    return matrix


def matmul(x: list[list[int]], y: list[list[int]]) -> list[list[int]]:
    size = len(x)
    return [
        [sum(x[row][k] * y[k][column] for k in range(size)) for column in range(size)]
        for row in range(size)
    ]


def matrix_power_traces(matrix: list[list[int]], count: int) -> list[int]:
    size = len(matrix)
    traces = []
    power = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    for _ in range(count):
        power = matmul(power, matrix)
        traces.append(sum(power[i][i] for i in range(size)))
    return traces


def minimal_polynomial(matrix: list[list[int]]) -> list[Fraction]:
    """Monic minimal polynomial coefficients (low to high) by exact algebra."""

    size = len(matrix)
    flat_powers: list[list[Fraction]] = []
    power = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    for _degree in range(size + 1):
        flat_powers.append(
            [Fraction(power[i][j]) for i in range(size) for j in range(size)]
        )
        power = matmul(power, matrix)
    for degree in range(1, size + 1):
        rows = flat_powers[: degree + 1]
        width = size * size
        system = [
            [rows[k][column] for k in range(degree)] + [rows[degree][column]]
            for column in range(width)
        ]
        solution = _solve_least_exact(system, degree)
        if solution is not None:
            coefficients = [-value for value in solution] + [Fraction(1)]
            return coefficients
    raise SpecificityError("minimal polynomial search failed")


def _solve_least_exact(
    system: list[list[Fraction]], unknowns: int
) -> list[Fraction] | None:
    rows = [row[:] for row in system]
    pivots: list[int] = []
    row_index = 0
    for column in range(unknowns):
        pivot_row = next(
            (
                candidate
                for candidate in range(row_index, len(rows))
                if rows[candidate][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            continue
        rows[row_index], rows[pivot_row] = rows[pivot_row], rows[row_index]
        pivot_value = rows[row_index][column]
        rows[row_index] = [value / pivot_value for value in rows[row_index]]
        for other in range(len(rows)):
            if other != row_index and rows[other][column] != 0:
                factor = rows[other][column]
                rows[other] = [
                    value - factor * pivot
                    for value, pivot in zip(rows[other], rows[row_index])
                ]
        pivots.append(column)
        row_index += 1
    for row in rows[row_index:]:
        if row[unknowns] != 0:
            return None
    solution = [Fraction(0)] * unknowns
    for position, column in enumerate(pivots):
        solution[column] = rows[position][unknowns]
    return solution


def polynomial_rational_roots(coefficients: list[Fraction]) -> list[Fraction]:
    """Strip rational roots; return them and mutate nothing."""

    working = coefficients[:]
    roots: list[Fraction] = []
    changed = True
    while changed and len(working) > 1:
        changed = False
        scale = 1
        for value in working:
            scale = scale * value.denominator // _gcd(scale, value.denominator)
        integers = [int(value * scale) for value in working]
        lead = integers[-1]
        constant = integers[0]
        if constant == 0:
            root = Fraction(0)
        else:
            root = None
            for p in _divisors(abs(constant)):
                for q in _divisors(abs(lead)):
                    for sign in (1, -1):
                        candidate = Fraction(sign * p, q)
                        if _evaluate(working, candidate) == 0:
                            root = candidate
                            break
                    if root is not None:
                        break
                if root is not None:
                    break
        if root is not None:
            roots.append(root)
            working = _deflate(working, root)
            changed = True
    return roots


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a) if a else 1


def _divisors(value: int) -> list[int]:
    return [d for d in range(1, value + 1) if value % d == 0] or [1]


def _evaluate(coefficients: list[Fraction], point: Fraction) -> Fraction:
    total = Fraction(0)
    for coefficient in reversed(coefficients):
        total = total * point + coefficient
    return total


def _deflate(coefficients: list[Fraction], root: Fraction) -> list[Fraction]:
    reversed_coefficients = list(reversed(coefficients))
    quotient: list[Fraction] = []
    carry = Fraction(0)
    for coefficient in reversed_coefficients[:-1]:
        carry = coefficient + carry * root
        quotient.append(carry)
    return list(reversed(quotient))


def spectral_structure(matrix: list[list[int]]) -> dict[str, Any]:
    minimal = minimal_polynomial(matrix)
    band_count = len(minimal) - 1
    rational_roots = polynomial_rational_roots(minimal)
    residual_degree = band_count - len(rational_roots)
    quadratic_pair = None
    if residual_degree == 2:
        remaining = minimal[:]
        for root in rational_roots:
            remaining = _deflate(remaining, root)
        c0, c1, _c2 = remaining[0], remaining[1], remaining[2]
        quadratic_pair = {"linear": str(c1), "constant": str(c0)}
    return {
        "band_count": band_count,
        "rational_bands": sorted(str(root) for root in rational_roots),
        "irrational_residual_degree": residual_degree,
        "quadratic_pair": quadratic_pair,
    }


def multiplicity_pattern(
    matrix: list[list[int]], structure: dict[str, Any]
) -> list[int] | None:
    """Exact band multiplicities via trace exhaustion, when solvable."""

    size = len(matrix)
    if structure["irrational_residual_degree"] not in (0, 2):
        return None
    rational = [Fraction(text) for text in structure["rational_bands"]]
    has_pair = structure["irrational_residual_degree"] == 2
    if has_pair:
        pair = structure["quadratic_pair"]
        if pair is None or Fraction(pair["linear"]) != 0:
            return None
        radicand = -Fraction(pair["constant"])
    traces = matrix_power_traces(matrix, 4)
    solutions = []
    slots = len(rational) + (1 if has_pair else 0)
    for combo in _compositions(size, slots):
        if has_pair and combo[-1] % 2 == 1:
            continue
        moments_match = True
        for power in range(1, 5):
            rational_sum = sum(
                Fraction(m) * root ** power
                for m, root in zip(combo, rational)
            )
            if has_pair:
                half = combo[-1] // 2
                if power % 2 == 0:
                    rational_sum += 2 * half * radicand ** (power // 2)
            if rational_sum != traces[power - 1]:
                moments_match = False
                break
        if moments_match:
            solutions.append(list(combo))
    if len(solutions) != 1:
        return None
    solution = solutions[0]
    multiplicities = solution[: len(rational)]
    if has_pair:
        half = solution[-1] // 2
        multiplicities.extend([half, half])
    return sorted(multiplicities)


def _compositions(total: int, slots: int):
    if slots == 1:
        yield (total,)
        return
    for head in range(total + 1):
        for tail in _compositions(total - head, slots - 1):
            yield (head, *tail)


def antipode_pairing(matrix: list[list[int]]) -> dict[int, int] | None:
    """The unique-farthest-vertex pairing, when it is an involutive automorphism."""

    size = len(matrix)
    distances = _all_distances(matrix)
    eccentricity = [max(row) for row in distances]
    diameter = max(eccentricity)
    pairing = {}
    for vertex in range(size):
        far = [
            other
            for other in range(size)
            if distances[vertex][other] == diameter
        ]
        if len(far) != 1:
            return None
        pairing[vertex] = far[0]
    if any(pairing[pairing[v]] != v or pairing[v] == v for v in range(size)):
        return None
    if not all(
        matrix[a][b] == matrix[pairing[a]][pairing[b]]
        for a in range(size)
        for b in range(size)
    ):
        return None
    return pairing


def antipodal_involution(matrix: list[list[int]]) -> bool:
    return antipode_pairing(matrix) is not None


def _all_distances(matrix: list[list[int]]) -> list[list[int]]:
    size = len(matrix)
    infinity = size + 1
    distances = [
        [0 if i == j else (1 if matrix[i][j] else infinity) for j in range(size)]
        for i in range(size)
    ]
    for k in range(size):
        for i in range(size):
            for j in range(size):
                through = distances[i][k] + distances[k][j]
                if through < distances[i][j]:
                    distances[i][j] = through
    require(
        all(value <= size for row in distances for value in row),
        "carrier is disconnected",
    )
    return distances


def response_involution(matrix: list[list[int]]) -> bool:
    if len(matrix) != 12:
        return False
    a2 = matmul(matrix, matrix)
    a3 = matmul(a2, matrix)
    size = 12
    j_numerator = [
        [
            a3[i][j] - 4 * a2[i][j] - 5 * matrix[i][j]
            + (10 if i == j else 0)
            for j in range(size)
        ]
        for i in range(size)
    ]
    j_squared = matmul(j_numerator, j_numerator)
    return all(
        j_squared[i][j] == (100 if i == j else 0)
        for i in range(size)
        for j in range(size)
    )


def probe_member(name: str) -> dict[str, Any]:
    matrix = adjacency_matrix(ENSEMBLE[name]())
    size = len(matrix)
    degrees = {sum(row) for row in matrix}
    structure = spectral_structure(matrix)
    multiplicities = multiplicity_pattern(matrix, structure)
    traces = matrix_power_traces(matrix, 3)
    probes = {
        "twelve_ports": size == 12,
        "degree_five_regular": degrees == {5},
        "four_bands": structure["band_count"] == 4,
        "galois_band_pair": structure["irrational_residual_degree"] == 2,
        "multiplicities_1_3_3_5": multiplicities == [1, 3, 3, 5],
        "unique_antipode_involution": antipodal_involution(matrix),
        "twenty_triangles": traces[2] // 6 == 20,
        "response_involution": response_involution(matrix),
    }
    cells = {}
    for probe, value in probes.items():
        if probe == "response_involution" and size != 12:
            cells[probe] = "failed_precondition_non_twelve_port"
        elif probe == "multiplicities_1_3_3_5" and multiplicities is None:
            cells[probe] = "non_identifiable_at_declared_trace_depth"
        else:
            cells[probe] = "hit" if value else "miss"
    gauge_stratum = (
        "reached"
        if all(probes.values())
        else (
            "failed_precondition_response_involution"
            if size == 12
            else "failed_precondition_non_twelve_port"
        )
    )
    return {
        "member": name,
        "vertices": size,
        "edges": sum(sum(row) for row in matrix) // 2,
        "band_count": structure["band_count"],
        "rational_bands": structure["rational_bands"],
        "band_multiplicities": multiplicities,
        "probes": probes,
        "cell_typing": cells,
        "gauge_stratum": gauge_stratum,
        "full_hit": all(probes.values()),
    }


# ---------------------------------------------------------------------------
# Response-law ensemble: the sixteen involutive sign laws
# ---------------------------------------------------------------------------


def _q5(a, b):
    return (Fraction(a), Fraction(b))


def _q5_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def _q5_mul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def _q5_scale(x, factor: Fraction):
    return (x[0] * factor, x[1] * factor)


def _q5_div(x, y):
    norm = y[0] * y[0] - 5 * y[1] * y[1]
    numerator = _q5_mul(x, (y[0], -y[1]))
    return (numerator[0] / norm, numerator[1] / norm)


def response_law_family() -> dict[str, Any]:
    """All sixteen involutive response laws on the icosahedral spectrum.

    A degree-three polynomial with values in {+1, -1} at the four bands is
    determined by its sign pattern through exact Lagrange interpolation over
    the quadratic field. The Galois-symmetric patterns (equal signs on the
    conjugate band pair) are exactly the rational laws, and the declared
    pattern reproduces the declared coefficients exactly.
    """

    nodes = [_q5(5, 0), _q5(0, 1), _q5(0, -1), _q5(-1, 0)]
    laws = []
    rational_count = 0
    declared_recovered = False
    for signs in product((1, -1), repeat=4):
        coefficients = [_q5(0, 0)] * 4
        for node_index, node in enumerate(nodes):
            basis = [_q5(1, 0)]
            denominator = _q5(1, 0)
            for other_index, other in enumerate(nodes):
                if other_index == node_index:
                    continue
                new_basis = [_q5(0, 0)] * (len(basis) + 1)
                for degree, coefficient in enumerate(basis):
                    new_basis[degree + 1] = _q5_add(
                        new_basis[degree + 1], coefficient
                    )
                    new_basis[degree] = _q5_add(
                        new_basis[degree],
                        _q5_mul(coefficient, (-other[0], -other[1])),
                    )
                basis = new_basis
                denominator = _q5_mul(
                    denominator,
                    _q5_add(node, (-other[0], -other[1])),
                )
            weight = _q5_div(_q5(signs[node_index], 0), denominator)
            padded = basis + [_q5(0, 0)] * (4 - len(basis))
            coefficients = [
                _q5_add(coefficient, _q5_mul(term, weight))
                for coefficient, term in zip(coefficients, padded)
            ]
        rational = all(part[1] == 0 for part in coefficients)
        if rational:
            rational_count += 1
        is_declared = signs == (1, -1, -1, 1)
        if is_declared:
            declared_recovered = [part[0] for part in coefficients] == [
                Fraction(1),
                Fraction(-1, 2),
                Fraction(-2, 5),
                Fraction(1, 10),
            ]
        laws.append(
            {
                "sign_pattern": list(signs),
                "rational_coefficients": rational,
                "coefficients_low_to_high": [
                    f"{part[0]}+{part[1]}*sqrt5" for part in coefficients
                ],
                "declared_law": is_declared,
            }
        )
    return {
        "band_nodes": ["5", "sqrt5", "-sqrt5", "-1"],
        "law_count": len(laws),
        "rational_law_count": rational_count,
        "declared_pattern": [1, -1, -1, 1],
        "declared_coefficients_recovered": bool(declared_recovered),
        "laws": laws,
    }


# ---------------------------------------------------------------------------
# Orientation stratum: antipode parity and the antipodal quotient
# ---------------------------------------------------------------------------

# Imported classification inputs of the replayed dimension argument. The
# compact simple Lie algebras of dimension at most twelve are su(2), su(3) and
# so(5). An inner A5 action on su(2) or su(3) is trivial or factors through an
# icosahedral SO(3) representation; the adjoint restrictions are 3 and 3 + 5,
# so the fixed space has dimension 0 or the full factor dimension.
COMPACT_SIMPLE_DIMENSIONS: dict[str, int] = {"su(2)": 3, "su(3)": 8, "so(5)": 10}
INNER_A5_FIXED_DIMENSIONS: dict[str, tuple[int, ...]] = {
    "su(2)": (0, 3),
    "su(3)": (0, 8),
}


def _polynomial_at_matrix(
    coefficients: list[Fraction], matrix: list[list[int]]
) -> list[list[Fraction]]:
    """Horner evaluation of a polynomial (low to high) at an integer matrix."""

    size = len(matrix)
    result = [[Fraction(0)] * size for _ in range(size)]
    for coefficient in reversed(coefficients):
        result = [
            [
                sum(result[i][k] * matrix[k][j] for k in range(size))
                + (coefficient if i == j else 0)
                for j in range(size)
            ]
            for i in range(size)
        ]
    return result


def band_projectors(
    matrix: list[list[int]],
) -> list[tuple[str, list[list[Fraction]]]]:
    """Exact spectral projectors: one per rational band, one for the rest.

    For a rational root r of the squarefree minimal polynomial m, the band
    projector is q(A)/q(r) with q = m/(x - r). The irrational bands together
    carry the identity minus the rational projectors, which is rational as
    well; it is labelled ``galois_pair`` when the irrational residual has
    degree two.
    """

    size = len(matrix)
    minimal = minimal_polynomial(matrix)
    structure = spectral_structure(matrix)
    projectors: list[tuple[str, list[list[Fraction]]]] = []
    total = [[Fraction(0)] * size for _ in range(size)]
    for root in sorted(polynomial_rational_roots(minimal)):
        quotient = _deflate(minimal, root)
        scale = _evaluate(quotient, root)
        require(scale != 0, "repeated root in a minimal polynomial")
        values = _polynomial_at_matrix(quotient, matrix)
        projector = [[value / scale for value in row] for row in values]
        projectors.append((str(root), projector))
        total = [
            [a + b for a, b in zip(row_total, row_projector)]
            for row_total, row_projector in zip(total, projector)
        ]
    residual = [
        [(1 if i == j else 0) - total[i][j] for j in range(size)]
        for i in range(size)
    ]
    if any(value != 0 for row in residual for value in row):
        label = (
            "galois_pair"
            if structure["irrational_residual_degree"] == 2
            else "irrational_bands"
        )
        projectors.append((label, residual))
    return projectors


def _exact_integer(value: Fraction) -> int:
    require(value.denominator == 1, "projector trace is not an integer")
    return int(value)


def band_dimensions(matrix: list[list[int]]) -> dict[str, int]:
    size = len(matrix)
    return {
        label: _exact_integer(sum(projector[i][i] for i in range(size)))
        for label, projector in band_projectors(matrix)
    }


def antipode_band_parities(matrix: list[list[int]]) -> list[dict[str, Any]] | None:
    """Trace of the antipode permutation on every band, or None without one."""

    pairing = antipode_pairing(matrix)
    if pairing is None:
        return None
    size = len(matrix)
    rows = []
    for label, projector in band_projectors(matrix):
        dimension = _exact_integer(sum(projector[i][i] for i in range(size)))
        trace = _exact_integer(sum(projector[pairing[i]][i] for i in range(size)))
        if trace == dimension:
            parity = "even"
        elif trace == -dimension:
            parity = "odd"
        else:
            parity = "mixed"
        rows.append(
            {
                "band": label,
                "dimension": dimension,
                "antipode_trace": trace,
                "parity": parity,
            }
        )
    return rows


def _simple_factor_multisets(total: int, names: list[str]) -> list[tuple[str, ...]]:
    if total == 0:
        return [()]
    multisets = []
    for position, name in enumerate(names):
        size = COMPACT_SIMPLE_DIMENSIONS[name]
        if size <= total:
            for rest in _simple_factor_multisets(total - size, names[position:]):
                multisets.append((name, *rest))
    return multisets


def replay_dimension_argument(ports: int) -> dict[str, Any]:
    """Replay the published compact-current dimension argument on ``ports``.

    Premises, as in the gauge-structure paper: the faithful response gives
    a compact reductive current algebra g with dim g equal to the port
    count; the transitive A5 action on the ports leaves one fixed line, so
    dim g^A5 = 1; the inner A5 action fixes the centre pointwise, so
    dim Z(g) <= 1. A split into a centre and simple factors is admissible
    when some choice of inner A5 actions on its factors gives a fixed space
    of dimension exactly one.
    """

    require(ports <= 12, "the imported simple-dimension menu covers twelve ports")
    candidates: list[str] = []
    admissible: list[str] = []
    for centre in (0, 1):
        for factors in _simple_factor_multisets(
            ports - centre, list(COMPACT_SIMPLE_DIMENSIONS)
        ):
            fixed = {centre}
            for factor in factors:
                require(
                    factor in INNER_A5_FIXED_DIMENSIONS,
                    f"no imported fixed-dimension menu for {factor}",
                )
                fixed = {a + b for a in fixed for b in INNER_A5_FIXED_DIMENSIONS[factor]}
            label = "+".join((["u(1)"] if centre else []) + list(factors))
            candidates.append(label)
            if 1 in fixed:
                admissible.append(label)
    return {"ports": ports, "candidates": candidates, "admissible": admissible}


def antipode_parity_table() -> dict[str, list[dict[str, Any]]]:
    table = {}
    for name in sorted(ENSEMBLE):
        rows = antipode_band_parities(adjacency_matrix(ENSEMBLE[name]()))
        if rows is not None:
            table[name] = rows
    return table


def orientation_stratum() -> dict[str, Any]:
    carrier = adjacency_matrix(icosahedron())
    parities = antipode_band_parities(carrier)
    require(parities is not None, "icosahedron antipode involution drift")
    even = {row["band"]: row["dimension"] for row in parities if row["parity"] == "even"}
    odd = {row["band"]: row["dimension"] for row in parities if row["parity"] == "odd"}
    require(
        all(row["parity"] != "mixed" for row in parities),
        "an icosahedral band mixes antipode parities",
    )
    quotient = adjacency_matrix(hemi_icosahedron())
    quotient_bands = band_dimensions(quotient)
    degree = sum(carrier[0])
    require(
        all(sum(row) == degree for row in carrier),
        "icosahedron is not regular",
    )
    charge_mode_band = str(degree)
    return {
        "carrier": "icosahedron",
        "antipodal_quotient_member": "hemi_icosahedron",
        "band_parities": parities,
        "even_sector_bands": even,
        "odd_sector_bands": odd,
        "quotient_bands": quotient_bands,
        "quotient_bands_equal_even_sector": quotient_bands == even,
        "galois_pair_parity": next(
            row["parity"] for row in parities if row["band"] == "galois_pair"
        ),
        "charge_mode_band": charge_mode_band,
        "charge_mode_parity": next(
            row["parity"] for row in parities if row["band"] == charge_mode_band
        ),
        "replayed_dimension_argument": {
            "icosahedron": replay_dimension_argument(len(carrier)),
            "hemi_icosahedron": replay_dimension_argument(len(quotient)),
        },
        "imported_inputs": {
            "compact_simple_dimensions": COMPACT_SIMPLE_DIMENSIONS,
            "inner_a5_fixed_dimensions": {
                name: list(values) for name, values in INNER_A5_FIXED_DIMENSIONS.items()
            },
            "a5_transitive_on_ports_and_on_antipodal_classes": True,
        },
        "claim_boundary": (
            "exact finite statements on the declared carriers: antipode "
            "parities of the adjacency bands, equality of the quotient "
            "spectrum with the antipode-even sector, and a replay of the "
            "published compact-current dimension argument on the imported "
            "inputs listed here; the hemi-icosahedron is a counterfactual "
            "ablation outside A1 (it is nonorientable), not an OPH branch; "
            "reading the antipode-odd sector as charge-conjugation-twisted "
            "configurations, whose total charge vanishes because the charge "
            "mode is antipode-even, consumes the orientation-reversal and "
            "charge-conjugation identification of the sector category and a "
            "physical charge attachment, neither of which is supplied here"
        ),
    }


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_receipt() -> dict[str, Any]:
    table = [probe_member(name) for name in sorted(ENSEMBLE)]
    invariants = {
        (row["vertices"], row["edges"], tuple(sorted(row["rational_bands"])))
        for row in table
    }
    require(
        len(invariants) == len(table),
        "ensemble members are not pairwise distinguished by invariants",
    )
    full_hits = [row["member"] for row in table if row["full_hit"]]
    probe_names = list(table[0]["probes"])
    per_probe = {
        probe: sum(1 for row in table if row["probes"][probe])
        for probe in probe_names
    }
    laws = response_law_family()
    require(laws["declared_coefficients_recovered"], "declared law drift")
    require(laws["rational_law_count"] == 8, "rational law count drift")
    receipt = {
        "schema": SCHEMA,
        "issue": 644,
        "status": STATUS,
        "ensemble_size": len(table),
        "carrier_hit_table": table,
        "per_probe_pass_counts": per_probe,
        "full_hit_members": full_hits,
        "unique_full_hit": full_hits == ["icosahedron"],
        "response_law_family": laws,
        "orientation_stratum": orientation_stratum(),
        "antipode_parity_table": antipode_parity_table(),
        "gauge_stratum_note": (
            "the current-algebra, center-quotient, charge-lattice, and "
            "matter-menu producers consume the twelve-port response "
            "involution; every member except the icosahedron fails that "
            "precondition, so the gauge stratum is typed failed_precondition "
            "rather than run, and only the icosahedron reaches it on the "
            "declared pipeline"
        ),
        "claim_boundary": (
            "the ensemble is a declared bounded menu of thirteen carriers and "
            "sixteen response laws; the score calibrates the "
            "orientation-blind spectral and response-law stratum of the "
            "registered hits and is not a physical forecast; the orientation "
            "stratum classifies antipode parity on the menu and the antipodal "
            "quotient of the icosahedron; the oriented structure beyond that "
            "quotient, the gauge-side producers beyond their typed "
            "precondition rows, and carriers outside the menu are not "
            "classified"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def write_runtime() -> Path:
    RUNTIME.mkdir(exist_ok=True)
    OUTPUT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    return OUTPUT_PATH


def verify_runtime() -> None:
    if OUTPUT_PATH.read_bytes() != canonical_json_bytes(build_receipt()):
        raise SystemExit("carrier specificity receipt is stale")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        print(write_runtime())
    if args.verify:
        verify_runtime()
        print("CARRIER_SPECIFICITY_VALID")
    if not args.write and not args.verify:
        receipt = build_receipt()
        print(
            json.dumps(
                {
                    "full_hit_members": receipt["full_hit_members"],
                    "per_probe_pass_counts": receipt["per_probe_pass_counts"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

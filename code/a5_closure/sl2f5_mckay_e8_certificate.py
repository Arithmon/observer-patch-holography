#!/usr/bin/env python3
"""Exact McKay certificate for the certified SL(2,F5) spin doublet.

This certificate starts only after the explicit bridge in
``sl2f5_port_spin_bridge_certificate.py`` has identified the canonical
``SL(2,F5)`` source with the executable ``PORT-SPIN-LIFT``.  It then derives
finite character data directly from that faithful exact two-dimensional
representation over Q(sqrt(5), i).

No character table is hard-coded.  The nine irreducible characters are
recovered from exact symmetric-power characters of the certified doublet and
orthogonal subtraction in the finite-group character inner product.  Tensoring
all recovered irreducibles by the certified doublet gives an exact graph, which
is compared with a separately encoded affine-E8 graph only after the fusion
matrix has been derived.

A Galois control applies sqrt(5) -> -sqrt(5) to the faithful doublet.  The
conjugate doublet is checked to remain faithful and to produce an isomorphic
affine-E8 McKay graph.  Therefore the affine-E8 graph type alone does not select
one of the two real embeddings of Q(sqrt(5)); in particular this certificate
does not select phi, state a mass law, or make a physical identification.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import sl2f5_port_spin_bridge_certificate as gl3

MODULE_DIR = Path(__file__).resolve().parent
SCHEMA = "oph.sl2f5_mckay_e8.v1"

c784 = gl3.c784
m314 = gl3.m314
Matrix = gl3.Matrix
SpinMatrix = gl3.SpinMatrix
Character = dict[Matrix, m314.C5]

CZ = m314.C5(m314.ZERO, m314.ZERO)
CO = m314.C5(m314.ONE, m314.ZERO)


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise c784.CertificateError(code, message)


def c5_integer(n: int) -> m314.C5:
    return m314.C5(m314.F5(n), m314.ZERO)


def c5_rational(q: Fraction) -> m314.C5:
    return m314.C5(m314.F5(q), m314.ZERO)


def c5_galois(z: m314.C5) -> m314.C5:
    """Apply sqrt(5) -> -sqrt(5), fixing i."""

    return m314.C5(z.re.conj(), z.im.conj())


def matrix_galois(matrix: Sequence[Sequence[m314.C5]]) -> SpinMatrix:
    return [[c5_galois(entry) for entry in row] for row in matrix]


def c5_key(z: m314.C5) -> tuple:
    return (z.re.a, z.re.b, z.im.a, z.im.b)


def c5_text(z: m314.C5) -> str:
    return z.text()


def character_zero(elements: Iterable[Matrix]) -> Character:
    return {g: CZ for g in elements}


def character_one(elements: Iterable[Matrix]) -> Character:
    return {g: CO for g in elements}


def character_add(left: Character, right: Character) -> Character:
    return {g: left[g] + right[g] for g in left}


def character_sub(left: Character, right: Character) -> Character:
    return {g: left[g] - right[g] for g in left}


def character_mul(left: Character, right: Character) -> Character:
    return {g: left[g] * right[g] for g in left}


def character_scale(character: Character, n: int) -> Character:
    scalar = c5_integer(n)
    return {g: scalar * value for g, value in character.items()}


def character_galois(character: Character) -> Character:
    return {g: c5_galois(value) for g, value in character.items()}


def character_equal(left: Character, right: Character) -> bool:
    return all(left[g] == right[g] for g in left)


def character_inner(left: Character, right: Character) -> m314.C5:
    total = CZ
    for g in left:
        total = total + left[g] * right[g].conj()
    return total * c5_rational(Fraction(1, len(left)))


def exact_integer(value: m314.C5, code: str, context: str) -> int:
    require(value.im.is_zero(), code, f"{context}: inner product has nonzero imaginary part {value.text()}")
    require(value.re.b == 0, code, f"{context}: inner product is not rational {value.text()}")
    require(value.re.a.denominator == 1, code, f"{context}: inner product is not integral {value.text()}")
    return int(value.re.a)


def inner_integer(left: Character, right: Character, context: str) -> int:
    return exact_integer(character_inner(left, right), "CHARACTER_INNER_PRODUCT", context)


def character_dimension(character: Character, identity: Matrix) -> int:
    return exact_integer(character[identity], "CHARACTER_DIMENSION", "character dimension")


def character_is_real(character: Character) -> bool:
    return all(value.im.is_zero() for value in character.values())


def symmetric_power_characters(spin: Character, identity: Matrix, count: int) -> list[Character]:
    """Characters of Sym^n(spin), using det(spin)=1: S_(n+1)=spin*S_n-S_(n-1)."""

    elements = list(spin)
    out = [character_one(elements), dict(spin)]
    for _n in range(1, count):
        out.append(character_sub(character_mul(spin, out[-1]), out[-2]))
    for n, chi in enumerate(out):
        require(
            character_dimension(chi, identity) == n + 1,
            "SYMMETRIC_POWER_DIMENSION",
            f"Sym^{n} character has wrong dimension",
        )
    return out


def subtract_known(character: Character, known: Mapping[str, Character]) -> tuple[Character, dict[str, int]]:
    residual = dict(character)
    multiplicities: dict[str, int] = {}
    for name, irrep in known.items():
        multiplicity = inner_integer(residual, irrep, f"projection onto {name}")
        require(multiplicity >= 0, "NEGATIVE_MULTIPLICITY", f"negative multiplicity for {name}")
        multiplicities[name] = multiplicity
        if multiplicity:
            residual = character_sub(residual, character_scale(irrep, multiplicity))
    return residual, multiplicities


def verify_irreducible(name: str, character: Character, identity: Matrix) -> None:
    dim = character_dimension(character, identity)
    require(dim > 0, "IRREP_DIMENSION", f"{name} has nonpositive dimension")
    require(inner_integer(character, character, f"norm of {name}") == 1, "IRREP_NORM", f"{name} is not irreducible")


def verify_orthonormal(characters: Mapping[str, Character]) -> None:
    names = list(characters)
    for i, left_name in enumerate(names):
        for j, right_name in enumerate(names):
            value = inner_integer(characters[left_name], characters[right_name], f"{left_name},{right_name}")
            expected = 1 if i == j else 0
            require(value == expected, "IRREP_ORTHOGONALITY", f"<{left_name},{right_name}>={value}, expected {expected}")


def reconstruct_bridge_data() -> tuple[c784.Source, dict[Matrix, SpinMatrix]]:
    """Rebuild the explicit GL3-A source-to-spin isomorphism, without trusting its JSON summary."""

    coset_manifest = c784.load_json(MODULE_DIR / "manifests" / "coset_carrier_reference.json")
    source = c784.build_source(coset_manifest)
    quotient, _summary = c784.central_quotient(source)
    fixture = c784.load_fixture()

    solutions = gl3.presentation_solutions(source)
    require(len(solutions) == 120, "PRESENTATION", f"expected 120 presentation solutions, got {len(solutions)}")
    source_generators = min(solutions)
    placement = gl3.placement_of_solution(quotient, source_generators)
    incidence = c784.build_incidence(quotient, placement)
    require(c784.classify(incidence) == "CARRIER", "PRESENTATION", "chosen presentation solution is not a carrier")
    phi, _relabel_summary = c784.relabel(quotient, placement, incidence, fixture)
    source_to_port = gl3.port_action_for_placement(quotient, placement, phi)

    matter_manifest = m314.load_json(MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json")
    upstream = m314.load_upstream(matter_manifest, MODULE_DIR)
    algebra = m314.CurrentAlgebra(upstream["current_manifest"], MODULE_DIR)
    spin = m314.spin_lift_certificate(algebra)
    raw_lifts: Mapping[gl3.PortRow, SpinMatrix] = spin["lifts"]
    current_phi, _current_relabel_summary = gl3.current_carrier_relabel(algebra, fixture)
    lifts = {
        gl3.relabel_port_row(tuple(raw_row), current_phi): lift
        for raw_row, lift in raw_lifts.items()
    }

    generator_rows = tuple(source_to_port(g) for g in source_generators)
    base_lifts = tuple(lifts[row] for row in generator_rows)
    choices = gl3.signed_generator_choices(base_lifts)
    require(len(choices) == 1, "GENERATOR_SIGN_CHOICE", f"expected one sign choice, got {len(choices)}")
    _signs, spin_generators = choices[0]
    source_to_spin = gl3.extend_generator_map(source, source_generators, spin_generators)

    identity2 = m314.cidentity(2)
    minus_identity2 = gl3.spin_neg(identity2)
    require(gl3.spin_equal(source_to_spin[source.identity], identity2), "CENTER_MAP", "+I does not map to +I2")
    require(gl3.spin_equal(source_to_spin[source.minus_identity], minus_identity2), "CENTER_MAP", "-I does not map to -I2")
    return source, source_to_spin


def spin_character(source: c784.Source, source_to_spin: Mapping[Matrix, SpinMatrix]) -> Character:
    chi = {g: m314.ctrace(source_to_spin[g]) for g in source.elements}
    require(character_dimension(chi, source.identity) == 2, "SPIN_CHARACTER", "spin character does not have dimension two")
    require(character_is_real(chi), "SPIN_CHARACTER", "certified spin character is not real-valued")
    require(inner_integer(chi, chi, "spin character norm") == 1, "SPIN_CHARACTER", "spin doublet is not irreducible")
    return chi


def derive_irreducibles(source: c784.Source, spin: Character) -> tuple[dict[str, Character], dict[str, Any]]:
    """Recover all nine irreducibles from Sym^n(spin), without a character table."""

    sym = symmetric_power_characters(spin, source.identity, 7)
    names = ["one", "spin2", "triplet3", "spin4", "quintet5", "spin6"]
    irreps: dict[str, Character] = {}
    for n, name in enumerate(names):
        verify_irreducible(name, sym[n], source.identity)
        irreps[name] = sym[n]
    verify_orthonormal(irreps)

    sym6_residual, sym6_known = subtract_known(sym[6], irreps)
    require(character_dimension(sym6_residual, source.identity) == 7, "SYMMETRIC_POWER_SPLIT", "unexpected known constituent in Sym^6")
    require(inner_integer(sym6_residual, sym6_residual, "Sym^6 residual norm") == 2, "SYMMETRIC_POWER_SPLIT", "Sym^6 residual should contain two new irreducibles")

    sym7_residual, sym7_known = subtract_known(sym[7], irreps)
    verify_irreducible("spin2_galois", sym7_residual, source.identity)
    require(character_dimension(sym7_residual, source.identity) == 2, "SYMMETRIC_POWER_SPLIT", "Sym^7 residual is not two-dimensional")
    irreps["spin2_galois"] = sym7_residual

    tensor_candidate = character_mul(spin, irreps["spin2_galois"])
    quotient4, tensor_known = subtract_known(tensor_candidate, irreps)
    verify_irreducible("quotient4", quotient4, source.identity)
    require(character_dimension(quotient4, source.identity) == 4, "SYMMETRIC_POWER_SPLIT", "spin2_galois tensor spin2 does not leave a four-dimensional irrep")
    irreps["quotient4"] = quotient4

    triplet3_galois, sym6_after = subtract_known(sym6_residual, {"quotient4": quotient4})
    verify_irreducible("triplet3_galois", triplet3_galois, source.identity)
    require(character_dimension(triplet3_galois, source.identity) == 3, "SYMMETRIC_POWER_SPLIT", "remaining Sym^6 constituent is not three-dimensional")
    irreps["triplet3_galois"] = triplet3_galois

    verify_orthonormal(irreps)
    dimensions = {name: character_dimension(chi, source.identity) for name, chi in irreps.items()}
    require(len(irreps) == 9, "IRREP_COMPLETENESS", f"expected nine irreducibles, got {len(irreps)}")
    require(sum(dim * dim for dim in dimensions.values()) == 120, "IRREP_COMPLETENESS", "sum of squared irreducible dimensions is not 120")

    require(
        character_equal(character_galois(irreps["spin2"]), irreps["spin2_galois"]),
        "GALOIS_IRREP",
        "derived second doublet is not the sqrt(5)-Galois conjugate of the certified spin doublet",
    )
    require(
        character_equal(character_galois(irreps["triplet3"]), irreps["triplet3_galois"]),
        "GALOIS_IRREP",
        "derived second triplet is not the sqrt(5)-Galois conjugate of the first triplet",
    )

    derivation = {
        "sym0_through_sym5_irreducible_dimensions": [dimensions[name] for name in names],
        "sym6_known_projections": sym6_known,
        "sym6_residual_dimension": character_dimension(sym6_residual, source.identity),
        "sym6_residual_norm": inner_integer(sym6_residual, sym6_residual, "Sym^6 residual norm report"),
        "sym7_known_projections": sym7_known,
        "sym7_residual_dimension": dimensions["spin2_galois"],
        "spin2_times_spin2_galois_known_projections": tensor_known,
        "sym6_final_projection": sym6_after,
    }
    return irreps, derivation


def conjugacy_classes(source: c784.Source) -> list[list[Matrix]]:
    unseen = set(source.elements)
    classes: list[list[Matrix]] = []
    while unseen:
        g = min(unseen)
        cls = {
            source.mul(source.mul(h, g), source.inv(h))
            for h in source.elements
        }
        require(cls <= unseen | (set(source.elements) - unseen), "CONJUGACY_CLASS", "internal class error")
        unseen -= cls
        classes.append(sorted(cls))
    classes.sort(key=lambda cls: (c784.element_order(source.mul, source.identity, cls[0]), len(cls), cls[0]))
    require(sum(len(cls) for cls in classes) == len(source.elements), "CONJUGACY_CLASS", "conjugacy classes do not partition the group")
    require(len(classes) == 9, "CONJUGACY_CLASS", f"expected nine conjugacy classes, got {len(classes)}")
    return classes


def verify_class_function(character: Character, classes: Sequence[Sequence[Matrix]], name: str) -> None:
    for cls in classes:
        value = character[cls[0]]
        require(all(character[g] == value for g in cls), "CLASS_FUNCTION", f"{name} is not constant on a conjugacy class")


def center_parity(source: c784.Source, character: Character) -> int:
    dim = character_dimension(character, source.identity)
    central = exact_integer(character[source.minus_identity], "CENTER_CHARACTER", "central character value")
    require(central in (-dim, dim), "CENTER_CHARACTER", "central -I does not act by a scalar sign")
    return 1 if central == dim else -1


def fusion_matrix(
    source: c784.Source,
    tensor_character: Character,
    irreps: Mapping[str, Character],
) -> tuple[list[str], list[list[int]]]:
    names = list(irreps)
    dimensions = [character_dimension(irreps[name], source.identity) for name in names]
    rows: list[list[int]] = []
    for i, name in enumerate(names):
        product = character_mul(tensor_character, irreps[name])
        row = [inner_integer(product, irreps[target], f"{name} tensor doublet -> {target}") for target in names]
        require(all(m >= 0 for m in row), "FUSION_MULTIPLICITY", f"negative fusion multiplicity in row {name}")
        require(
            sum(m * dimensions[j] for j, m in enumerate(row)) == 2 * dimensions[i],
            "FUSION_DIMENSION",
            f"fusion row {name} has wrong total dimension",
        )
        rows.append(row)
    return names, rows


def graph_edges(names: Sequence[str], matrix: Sequence[Sequence[int]]) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            multiplicity = matrix[i][j]
            require(multiplicity in (0, 1), "MCKAY_MULTIPLICITY", "McKay graph is not simply-laced")
            require(i != j or multiplicity == 0, "MCKAY_LOOP", "McKay graph has a loop")
            require(matrix[i][j] == matrix[j][i], "MCKAY_SYMMETRY", "McKay fusion matrix is not symmetric")
            if i < j and multiplicity:
                edges.add((left, right))
    return edges


def connected(names: Sequence[str], edges: set[tuple[str, str]]) -> bool:
    if not names:
        return True
    adjacency = {name: set() for name in names}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen = {names[0]}
    queue = [names[0]]
    for node in queue:
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return len(seen) == len(names)


E8_NODES = ("e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8")
E8_DIMENSIONS = {
    "e0": 1,
    "e1": 2,
    "e2": 3,
    "e3": 4,
    "e4": 5,
    "e5": 6,
    "e6": 4,
    "e7": 2,
    "e8": 3,
}
E8_EDGES = {
    ("e0", "e1"),
    ("e1", "e2"),
    ("e2", "e3"),
    ("e3", "e4"),
    ("e4", "e5"),
    ("e5", "e6"),
    ("e6", "e7"),
    ("e5", "e8"),
}


def normalized_edge(edge: tuple[str, str]) -> tuple[str, str]:
    left, right = edge
    return (left, right) if left < right else (right, left)


def affine_e8_isomorphisms(
    names: Sequence[str],
    dimensions: Mapping[str, int],
    edges: set[tuple[str, str]],
) -> list[dict[str, str]]:
    source_by_dim: dict[int, list[str]] = {}
    target_by_dim: dict[int, list[str]] = {}
    for name in names:
        source_by_dim.setdefault(dimensions[name], []).append(name)
    for node in E8_NODES:
        target_by_dim.setdefault(E8_DIMENSIONS[node], []).append(node)
    require(set(source_by_dim) == set(target_by_dim), "E8_DIMENSIONS", "dimension multisets differ from affine E8 marks")
    for dim in source_by_dim:
        require(len(source_by_dim[dim]) == len(target_by_dim[dim]), "E8_DIMENSIONS", f"multiplicity mismatch at dimension {dim}")

    groups = []
    for dim in sorted(source_by_dim):
        source_group = sorted(source_by_dim[dim])
        target_group = sorted(target_by_dim[dim])
        groups.append((source_group, list(itertools.permutations(target_group))))

    target_edges = {normalized_edge(edge) for edge in E8_EDGES}
    out: list[dict[str, str]] = []
    for choices in itertools.product(*(permutations for _source, permutations in groups)):
        mapping: dict[str, str] = {}
        for (source_group, _permutations), target_order in zip(groups, choices):
            mapping.update(dict(zip(source_group, target_order)))
        mapped_edges = {
            normalized_edge((mapping[left], mapping[right]))
            for left, right in edges
        }
        if mapped_edges == target_edges:
            out.append(mapping)
    return out


def perron_mark_check(
    names: Sequence[str],
    matrix: Sequence[Sequence[int]],
    dimensions: Mapping[str, int],
) -> bool:
    index = {name: i for i, name in enumerate(names)}
    for name in names:
        i = index[name]
        if sum(matrix[i][j] * dimensions[names[j]] for j in range(len(names))) != 2 * dimensions[name]:
            return False
    return True


def galois_faithful(
    source: c784.Source,
    source_to_spin: Mapping[Matrix, SpinMatrix],
) -> tuple[bool, list[Matrix]]:
    identity2 = m314.cidentity(2)
    kernel = [
        g for g in source.elements
        if gl3.spin_equal(matrix_galois(source_to_spin[g]), identity2)
    ]
    return kernel == [source.identity], kernel


def class_table_payload(
    source: c784.Source,
    classes: Sequence[Sequence[Matrix]],
    irreps: Mapping[str, Character],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, cls in enumerate(classes):
        representative = cls[0]
        rows.append(
            {
                "class": index,
                "size": len(cls),
                "element_order": c784.element_order(source.mul, source.identity, representative),
                "representative": list(representative),
                "characters": {name: c5_text(character[representative]) for name, character in irreps.items()},
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    source, source_to_spin = reconstruct_bridge_data()
    spin = spin_character(source, source_to_spin)
    irreps, derivation = derive_irreducibles(source, spin)
    classes = conjugacy_classes(source)
    for name, character in irreps.items():
        verify_class_function(character, classes, name)

    dimensions = {name: character_dimension(chi, source.identity) for name, chi in irreps.items()}
    parities = {name: center_parity(source, chi) for name, chi in irreps.items()}
    require(sum(1 for sign in parities.values() if sign == 1) == 5, "CENTER_PARITY", "expected five center-trivial irreducibles")
    require(sum(1 for sign in parities.values() if sign == -1) == 4, "CENTER_PARITY", "expected four spinorial irreducibles")

    names, fusion = fusion_matrix(source, spin, irreps)
    edges = graph_edges(names, fusion)
    require(len(edges) == 8, "MCKAY_GRAPH", f"expected eight edges, got {len(edges)}")
    require(connected(names, edges), "MCKAY_GRAPH", "McKay graph is disconnected")
    require(len(edges) == len(names) - 1, "MCKAY_GRAPH", "McKay graph is not a tree")
    require(perron_mark_check(names, fusion, dimensions), "MCKAY_PERRON", "dimension vector does not satisfy A d = 2 d")
    e8_maps = affine_e8_isomorphisms(names, dimensions, edges)
    require(e8_maps, "MCKAY_E8", "derived McKay graph is not isomorphic to affine E8")

    galois_spin = character_galois(spin)
    galois_faithful_ok, galois_kernel = galois_faithful(source, source_to_spin)
    require(galois_faithful_ok, "GALOIS_FAITHFUL", f"Galois-conjugate doublet has kernel of size {len(galois_kernel)}")
    decomposition = {name: inner_integer(galois_spin, chi, f"Galois spin against {name}") for name, chi in irreps.items()}
    require(sum(decomposition.values()) == 1, "GALOIS_DOUBLET", "Galois-conjugate spin character is not one recovered irreducible")
    galois_name = next(name for name, multiplicity in decomposition.items() if multiplicity == 1)
    require(dimensions[galois_name] == 2, "GALOIS_DOUBLET", "Galois-conjugate spin character is not two-dimensional")
    require(galois_name != "spin2", "GALOIS_DOUBLET", "Galois conjugation did not exchange the faithful doublet")

    g_names, g_fusion = fusion_matrix(source, galois_spin, irreps)
    require(g_names == names, "GALOIS_MCKAY", "Galois fusion order changed unexpectedly")
    g_edges = graph_edges(g_names, g_fusion)
    g_e8_maps = affine_e8_isomorphisms(g_names, dimensions, g_edges)
    require(g_e8_maps, "GALOIS_MCKAY", "Galois-conjugate doublet does not produce affine E8")

    two_trivial = character_scale(irreps["one"], 2)
    _bad_names, bad_fusion = fusion_matrix(source, two_trivial, irreps)
    bad_edges = set()
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            if i < j and bad_fusion[i][j]:
                bad_edges.add((left, right))
    bad_e8_maps = affine_e8_isomorphisms(names, dimensions, bad_edges) if bad_edges else []
    require(not bad_e8_maps, "NEGATIVE_CONTROL", "reducible two-trivial representation falsely reproduces affine E8")

    spin_trace_values = sorted({c5_text(value) for value in spin.values()})
    galois_trace_values = sorted({c5_text(value) for value in galois_spin.values()})
    irrational_spin_values = sorted({c5_text(value) for value in spin.values() if value.re.b != 0 or value.im.b != 0})
    require(irrational_spin_values, "GOLDEN_FIELD", "spin character contains no sqrt(5)-sensitive trace values")
    require(not character_equal(spin, galois_spin), "GOLDEN_FIELD", "Galois action fixes the spin character")

    return {
        "schema": SCHEMA,
        "source": {
            "group": "SL(2,F5)",
            "order": len(source.elements),
            "conjugacy_class_count": len(classes),
            "certified_doublet_field": "Q(sqrt(5), i)",
        },
        "irreducible_recovery": {
            "method": "exact symmetric-power recurrence plus finite-group character inner products",
            "hard_coded_character_table": False,
            "irreducible_count": len(irreps),
            "dimensions": dimensions,
            "sum_squared_dimensions": sum(dim * dim for dim in dimensions.values()),
            "center_signs": parities,
            "derivation": derivation,
        },
        "character_table": class_table_payload(source, classes, irreps),
        "mckay": {
            "tensor_doublet": "spin2",
            "node_order": names,
            "fusion_matrix": fusion,
            "edges": [list(edge) for edge in sorted(edges)],
            "affine_e8_isomorphism_count": len(e8_maps),
            "affine_e8_isomorphism": e8_maps[0],
            "dimension_vector_is_two_eigenvector": True,
            "connected_tree": True,
        },
        "galois_control": {
            "automorphism": "sqrt(5) -> -sqrt(5), i fixed",
            "conjugate_doublet": galois_name,
            "faithful": galois_faithful_ok,
            "kernel_size": len(galois_kernel),
            "fusion_matrix": g_fusion,
            "edges": [list(edge) for edge in sorted(g_edges)],
            "affine_e8_isomorphism_count": len(g_e8_maps),
            "same_affine_e8_graph_type": True,
            "same_labeled_graph": g_edges == edges,
            "spin_trace_values": spin_trace_values,
            "galois_trace_values": galois_trace_values,
            "sqrt5_sensitive_spin_trace_values": irrational_spin_values,
            "affine_e8_graph_type_selects_galois_embedding": False,
        },
        "negative_control": {
            "representation": "1 direct-sum 1",
            "dimension": 2,
            "affine_e8_isomorphism_count": len(bad_e8_maps),
            "passes": not bad_e8_maps,
        },
        "claim_boundary": {
            "proves": (
                "the McKay graph of the certified faithful two-dimensional SL(2,F5) spin representation "
                "is exactly affine E8, with all nine irreducible characters recovered from exact character arithmetic"
            ),
            "also_proves": (
                "the sqrt(5)-Galois-conjugate faithful doublet produces the same affine-E8 graph type, "
                "so graph type alone does not select one real embedding of Q(sqrt(5))"
            ),
            "does_not_prove": [
                "selection or derivation of phi rather than its Galois conjugate",
                "27^phi or any mass relation",
                "physical identification of the finite spin action",
                "source selection of the charged-double-triplet current fixture",
                "a typed Lean McKay or affine-E8 theorem",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("verify", "print"), nargs="?", default="verify")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    certificate = build_certificate()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8", newline="\n")
    if args.command == "print" or args.output is None:
        print(encoded, end="")


if __name__ == "__main__":
    main()

"""Exact incidence, orientation, covering and deformation controls."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product

import pytest

from . import statistics_topology as t
from ..clock_defect_energy import defect_history as history
from ..observer import native_loop_readout


def point(weights):
    return tuple(F(weights.get(i, 0)) for i in range(12))


def project(cell):
    left, right = cell
    if (len(left), left) < (len(right), right):
        return cell, 1
    return (right, left), (-1)**((len(left)-1)*(len(right)-1))


def independent_boundary(cell, ordered):
    """Tensor boundary before projection; no production boundary helper."""
    left, right = cell
    terms = Counter()
    for sign, fixed, moving, first in ((1, right, left, True),
                                      ((-1)**(len(left)-1), left, right, False)):
        if len(moving) == 1:
            continue
        for i in range(len(moving)):
            face = tuple(v for j, v in enumerate(moving) if i != j)
            child = (face, fixed) if first else (fixed, face)
            orientation = 1
            if not ordered:
                child, orientation = project(child)
            terms[child] += sign*((-1)**i)*orientation
    return {cell: value for cell, value in terms.items() if value}


def test_exact_surface_and_full_cell_census_matches_existing_two_skeleton():
    result = t.certificate()
    assert result["surface"]["vertex_link_lengths"] == (5,)*12
    assert result["unordered_cell_counts"] == (66, 300, 495, 360, 100)
    assert result["ordered_cover"]["ordered_cell_counts"] == (132, 600, 990, 720, 200)
    assert result["full_deleted_product_euler_characteristic"] == 1
    assert result["native_two_skeleton_euler_characteristic"] == 261
    assert result["native_fundamental_group"] == "Z/2"
    assert result["ordered_cover_fundamental_group"] == "trivial"
    assert result["unitary_irreducible_characters"] == (1, -1)
    assert result["native_quantum_character_selected"] is False
    assert result["native_two_skeleton"] == {
        "vertices": 66, "edges": 300, "two_cells": 495,
        "all_oriented_boundaries_match_up_to_cell_orientation": True}
    kinds = Counter(tuple(sorted((len(a)-1, len(b)-1)))
                    for layer in t.deleted_cells() for a, b in layer)
    assert kinds == {(0, 0): 66, (0, 1): 300, (0, 2): 180,
                     (1, 1): 315, (1, 2): 360, (2, 2): 100}


@pytest.mark.parametrize("ordered", (False, True))
def test_every_integral_boundary_and_boundary_squared_independently(ordered):
    layers = t.deleted_cells(ordered=ordered)
    for dimension, layer in enumerate(layers):
        for cell in layer:
            expected = independent_boundary(cell, ordered)
            assert dict(t.cell_boundary(cell, ordered=ordered)) == expected
            assert all(face in layers[dimension-1] for face in expected)
            twice = Counter()
            for face, coefficient in expected.items():
                for subface, second in independent_boundary(face, ordered).items():
                    twice[subface] += coefficient*second
            assert not any(twice.values())


def test_deck_swap_is_free_and_has_necessary_orientation_sign():
    for layer in t.deleted_cells(ordered=True):
        for cell in layer:
            left, right = cell
            swapped = (right, left)
            assert swapped != cell and swapped in layer
            sign = (-1)**((len(left)-1)*(len(right)-1))
            moved_boundary = Counter()
            for (a, b), coefficient in t.cell_boundary(cell, ordered=True):
                moved_boundary[(b, a)] += coefficient*((-1)**((len(a)-1)*(len(b)-1)))
            assert dict(moved_boundary) == {
                face: sign*coefficient for face, coefficient in t.cell_boundary(swapped, ordered=True)}
    square = next(cell for cell in t.deleted_cells()[2] if len(cell[0]) == len(cell[1]) == 2)
    assert t.quotient_cell(square[::-1]) == (square, -1)
    # Forgetting the Koszul sign reverses the square's boundary.
    assert independent_boundary(square[::-1], False) == {
        face: -coefficient for face, coefficient in independent_boundary(square, False).items()}


def test_cover_exchange_lift_and_all_admitted_faces():
    pairs, _, cells, _, _, _ = history.geometry()
    _, path = history.exchange_loop()
    ordered = pairs[path[0]]
    for start, finish in zip(path, path[1:]):
        old = next(iter(set(pairs[start])-set(pairs[finish])))
        new = next(iter(set(pairs[finish])-set(pairs[start])))
        ordered = tuple(new if v == old else v for v in ordered)
    assert ordered == pairs[path[0]][::-1]
    for _, vertices in cells:
        assert history.transport(vertices+vertices[:1]) == 0
    # A proposed extra filling of the single exchange would not lift closed.
    assert history.transport(path) == 1


def test_removing_declared_triangles_changes_flat_character_domain():
    result = t.graph_only_certificate()
    assert result["square_boundary_rank_mod2"] == 215
    assert result["first_cohomology_mod2_dimension"] == 20
    assert result["graph_only_fundamental_group_computed"] is False
    assert result["changes_declared_deformation_contract"] is True
    witness = result["continuous_character_witness"]
    origin, destination = witness["directed_one_particle_edge"]
    pairs, edges, cells, _, _, _ = history.geometry()

    def exponent(left, right):
        old, new = set(left)-set(right), set(right)-set(left)
        return int(old == {origin} and new == {destination})-int(
            old == {destination} and new == {origin})

    # This integer identity verifies every angle at once, not only a chosen
    # root of unity.  The removed triangle is an explicit violation of the
    # ORIGINAL flatness contract for generic theta.
    for kind, vertices in cells:
        value = sum(exponent(pairs[a], pairs[b])
                    for a, b in zip(vertices, vertices[1:]+vertices[:1]))
        if kind == "disjoint_square":
            assert value == 0
    path = witness["removed_triangle_occupied_pairs"]
    assert sum(exponent(a, b) for a, b in zip(path, path[1:])) == 1
    # Independently eliminate F_2 boundary columns in the original edge
    # basis, without using the production rank or boundary functions.
    positions = {edge: i for i, edge in enumerate(edges)}
    pivots = {}
    for kind, vertices in cells:
        if kind != "disjoint_square":
            continue
        support = {positions[tuple(sorted((a, b)))]
                   for a, b in zip(vertices, vertices[1:]+vertices[:1])}
        while support and min(support) in pivots:
            support ^= pivots[min(support)]
        if support:
            pivots[min(support)] = support
    assert len(pivots) == 215


def test_surface_validator_rejects_boundary_duplicate_and_wrong_euler_manifold():
    edges, faces = native_loop_readout.carrier()
    vertices = tuple(range(12))
    with pytest.raises(ValueError, match="two incident"):
        t.sphere_triangulation_certificate(vertices, edges, faces[:-1])
    with pytest.raises(ValueError, match="duplicate triangle"):
        t.sphere_triangulation_certificate(vertices, edges, faces+(faces[0][::-1],))
    with pytest.raises(ValueError, match="boundary edge absent"):
        t.sphere_triangulation_certificate(vertices, edges[1:], faces)
    # A 4x4 periodic triangulation is a connected oriented closed torus,
    # with circular links: Euler two is a separate necessary premise.
    at = lambda x, y: 4*(x % 4)+(y % 4)
    torus_faces = []
    for x, y in product(range(4), repeat=2):
        a, b, c, d = at(x, y), at(x+1, y), at(x+1, y+1), at(x, y+1)
        torus_faces.extend(((a, b, c), (a, c, d)))
    torus_edges = tuple(sorted({tuple(sorted(edge)) for face in torus_faces
                                for edge in combinations(face, 2)}))
    with pytest.raises(ValueError, match="Euler characteristic two"):
        t.sphere_triangulation_certificate(tuple(range(16)), torus_edges, tuple(torus_faces))


@pytest.mark.parametrize("time", (F(0), F(1, 7), F(1, 2), F(1)))
def test_barycentric_deformation_exact_invariants_on_all_face_pairs(time):
    faces = t.simplices()[2]
    for left_face, right_face in product(faces, repeat=2):
        left = point(dict(zip(left_face, (F(1, 2), F(1, 3), F(1, 6)))))
        right = point(dict(zip(right_face, (F(1, 5), F(2, 5), F(2, 5)))))
        out_left, out_right = t.deleted_product_deformation(left, right, time)
        assert out_left != out_right
        assert sum(out_left) == sum(out_right) == 1
        assert all(x >= 0 for x in out_left+out_right)
        assert all(x == 0 or v in left_face for v, x in enumerate(out_left))
        assert all(x == 0 or v in right_face for v, x in enumerate(out_right))
        assert t.deleted_product_deformation(right, left, time) == (out_right, out_left)
        common = sum(min(a, b) for a, b in zip(left, right))
        assert tuple(a-b for a, b in zip(out_left, out_right)) == tuple(
            (a-b)/(1-time*common) for a, b in zip(left, right))
        if time == 0:
            assert (out_left, out_right) == (left, right)
        if time == 1:
            assert all(not (a and b) for a, b in zip(out_left, out_right))
        if set(left_face).isdisjoint(right_face):
            assert (out_left, out_right) == (left, right)


def test_barycentric_boundary_and_near_diagonal_do_not_need_subdivision():
    a, b, c = t.simplices()[2][0]
    left = point({a: F(1, 2), b: F(1, 2)})
    right = point({a: F(1, 2)-F(1, 10**20), b: F(1, 2), c: F(1, 10**20)})
    assert t.deleted_product_deformation(left, right, 1) == (point({a: 1}), point({c: 1}))
    for cell in t.deleted_cells()[4]:
        a, b = (point({v: F(1, len(simplex)) for v in simplex}) for simplex in cell)
        assert t.deleted_product_deformation(a, b, F(3, 7)) == (a, b)


@pytest.mark.parametrize("time", (F(0), F(1, 9), F(1, 2), F(1)))
def test_exact_antipodal_deformation_normalized_numerators(time):
    points = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, 0, -1),
              (F(3, 5), F(4, 5), 0), (F(2, 3), F(-2, 3), F(1, 3)))
    for left, right in product(points, repeat=2):
        if left == right:
            continue
        a, b, norm = t.antipodal_deformation(left, right, time)
        assert norm > 0 and norm == sum(x*x for x in a) == sum(x*x for x in b)
        assert a != b
        assert t.antipodal_deformation(right, left, time) == (b, a, norm)
        assert tuple(x-y for x, y in zip(a, b)) == tuple(F(x)-F(y) for x, y in zip(left, right))
        if time == 0 or tuple(-F(x) for x in left) == right:
            assert (a, b, norm) == (left, right, 1)
        if time == 1:
            assert a == tuple(-x for x in b)


def test_deformation_and_cell_domain_guards_fail_after_warm_calls():
    face = t.simplices()[2][0]
    left, right = point({face[0]: 1}), point({face[1]: 1})
    t.deleted_product_deformation(left, right, 1)
    for bad in (True, 1.0, F(-1), F(2)):
        with pytest.raises(ValueError):
            t.deleted_product_deformation(left, right, bad)
    with pytest.raises(ValueError, match="distinct"):
        t.deleted_product_deformation(left, left, F(1, 2))
    nonedge = next(pair for pair in combinations(range(12), 2) if pair not in t.simplices()[1])
    with pytest.raises(ValueError, match="native simplex"):
        t.deleted_product_deformation(point({v: F(1, 2) for v in nonedge}), right, 1)
    with pytest.raises(ValueError, match="unit-sphere"):
        t.antipodal_deformation((1, 1, 0), (0, 1, 0), 1)
    with pytest.raises(ValueError, match="distinct"):
        t.antipodal_deformation((1, 0, 0), (1, 0, 0), 0)
    cell = t.deleted_cells()[2][0]
    t.cell_boundary(cell)
    with pytest.raises(ValueError, match="canonical"):
        t.cell_boundary(cell[::-1])
    with pytest.raises(ValueError, match="vertex-disjoint"):
        t.cell_boundary((face, (face[0],)), ordered=True)
    with pytest.raises(ValueError, match="boolean"):
        t.deleted_cells(ordered=1)
    with pytest.raises(ValueError, match="boolean"):
        t.cell_boundary(cell, ordered=0)

"""Independent face, homology and hostile-input controls for finite incidence."""

import itertools
import sys
from pathlib import Path

import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from finite_incidence import (
    IncidenceComplex, certify_midpoint_subdivision, complexes_equal,
    complexes_isomorphic_under, dimension, euler_characteristic,
    incidence_complex, refinement_is_simplicial, validate_complex,
)
from quotient_cap_readout import (
    RepairSystem, classify_surface, is_closed_surface, orient, refinement_subdivide,
    spherical_incidence_receipt,
)


TETRA_BOUNDARY = list(itertools.combinations(range(4), 3))


def _gf2_rank(columns):
    """Bit-column Gaussian elimination, independent of the implementation."""
    pivots = {}
    for column in columns:
        while column:
            lead = column.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = column
                break
            column ^= pivots[lead]
    return len(pivots)


def _betti(k):
    faces = [[frozenset((v,)) for v in k.vertices], sorted(k.edges, key=sorted),
             sorted(k.triangles, key=sorted)]
    for size in range(4, max((len(f) for f in k.higher_simplices), default=3) + 1):
        faces.append(sorted((f for f in k.higher_simplices if len(f) == size), key=sorted))
    ranks = [0]
    for lower, upper in zip(faces, faces[1:]):
        ranks.append(_gf2_rank([sum(1 << lower.index(f - {v}) for v in f) for f in upper]))
    ranks.append(0)
    return [len(f) - ranks[i] - ranks[i + 1] for i, f in enumerate(faces)]


def test_joint_four_patch_record_is_a_ball_not_its_spherical_boundary():
    solid = incidence_complex([(0, 1, 2, 3)])
    boundary = incidence_complex(TETRA_BOUNDARY)
    assert solid.edges == boundary.edges and solid.triangles == boundary.triangles
    assert not complexes_equal(solid, boundary)
    assert dimension(solid) == 3 and dimension(boundary) == 2
    assert euler_characteristic(solid) == 1 and euler_characteristic(boundary) == 2
    assert _betti(solid) == [1, 0, 0, 0]
    assert _betti(boundary) == [1, 0, 1]
    assert classify_surface(solid) == 'NOT_A_CLOSED_SURFACE'
    assert not spherical_incidence_receipt(solid)
    assert spherical_incidence_receipt(boundary)


@pytest.mark.parametrize('size', range(1, 9))
def test_full_support_simplex_has_all_faces_and_contractible_homology(size):
    k = incidence_complex([tuple(range(size))])
    assert len(k.vertices) + len(k.edges) + len(k.triangles) + len(k.higher_simplices) == 2**size - 1
    assert dimension(k) == size - 1
    assert euler_characteristic(k) == 1
    assert _betti(k) == [1] + [0] * (max(3, size) - 1)


@pytest.mark.parametrize('records', [[()], [(0, 0, 1)], [(True, 2)], [(0, 1.0)],
                                     [(0, 'a')], [None], '012'])
def test_rejects_malformed_record_support(records):
    with pytest.raises(ValueError):
        incidence_complex(records)


@pytest.mark.parametrize('bad', [
    IncidenceComplex([0, 0], set(), set()),
    IncidenceComplex([True], set(), set()),
    IncidenceComplex([0, 1], {frozenset((0,))}, set()),
    IncidenceComplex([0, 1], {frozenset((0, 2))}, set()),
    IncidenceComplex([0, 1, 2], set(), {frozenset((0, 1, 2))}),
    IncidenceComplex(list(range(4)), set(), set(), {frozenset(range(4))}),
    IncidenceComplex([0], [], set()),
])
def test_malformed_complex_never_receives_a_topological_verdict(bad):
    for operation in (validate_complex, spherical_incidence_receipt, classify_surface):
        with pytest.raises(ValueError):
            operation(bad)


def test_empty_and_isolated_vertices_are_not_closed_surfaces():
    for records in ([], [(0,)], TETRA_BOUNDARY + [(10,)]):
        assert not is_closed_surface(incidence_complex(records))
        assert not spherical_incidence_receipt(incidence_complex(records))


def test_simplicial_check_includes_zero_one_and_higher_dimensional_faces():
    edge = incidence_complex([(0, 1)])
    vertices = incidence_complex([(0,), (1,)])
    assert not refinement_is_simplicial(edge, vertices, {0: 0, 1: 1})
    assert not refinement_is_simplicial(vertices, vertices, {0: 9, 1: 9})
    assert not refinement_is_simplicial(vertices, vertices, {0: 0})
    assert not refinement_is_simplicial(vertices, vertices, {0: 0, 1: 1, 2: 1})
    assert not refinement_is_simplicial(vertices, vertices, {0: False, 1: 1})
    sphere = incidence_complex(TETRA_BOUNDARY)
    solid = incidence_complex([tuple(range(4))])
    identity = {i: i for i in range(4)}
    assert not refinement_is_simplicial(solid, sphere, identity)
    assert refinement_is_simplicial(sphere, solid, identity)


def test_simplicial_and_classical_channel_properties_do_not_preserve_topology():
    sphere, point = incidence_complex(TETRA_BOUNDARY), incidence_complex([(7,)])
    collapse = {v: 7 for v in sphere.vertices}
    assert refinement_is_simplicial(sphere, point, collapse)
    # Deterministic classical pushforward: nonnegative columns, each summing to one.
    stochastic_matrix = [[int(collapse[v] == 7) for v in sphere.vertices]]
    assert all(sum(column) == 1 for column in zip(*stochastic_matrix))
    assert _betti(sphere)[2] == 1 and _betti(point)[2] == 0
    assert not complexes_isomorphic_under(sphere, point, collapse)


def _midpoints(coarse, projection):
    # Producer uses sorted edges; certificate is exposed independently to mutations.
    start = max(coarse.vertices) + 1
    return {e: start + i for i, e in enumerate(sorted(coarse.edges, key=lambda e: tuple(sorted(e))))}


def test_subdivision_full_certificate_and_independent_homology():
    coarse = incidence_complex(TETRA_BOUNDARY)
    records, projection = refinement_subdivide(TETRA_BOUNDARY)
    fine = incidence_complex(records)
    midpoints = _midpoints(coarse, projection)
    assert certify_midpoint_subdivision(coarse, fine, midpoints)
    assert refinement_is_simplicial(fine, coarse, projection)
    assert _betti(fine) == _betti(coarse) == [1, 0, 1]
    assert len(fine.triangles) == 4 * len(coarse.triangles)
    # Endpoint projection collapses fine edges; it is not itself a homeomorphism.
    assert any(len({projection[v] for v in e}) == 1 for e in fine.edges)
    assert not certify_midpoint_subdivision(coarse, coarse, midpoints)
    assert not certify_midpoint_subdivision(coarse, incidence_complex(records[:-1]), midpoints)
    assert not certify_midpoint_subdivision(coarse, incidence_complex(records + [(99,)]), midpoints)
    assert not certify_midpoint_subdivision(coarse, fine, {})
    wrong = dict(midpoints)
    first, second = list(wrong)[:2]
    wrong[first] = wrong[second]
    assert not certify_midpoint_subdivision(coarse, fine, wrong)
    wrong[first] = coarse.vertices[0]
    assert not certify_midpoint_subdivision(coarse, fine, wrong)


def test_midpoint_realization_covers_each_triangle_with_four_equal_areas():
    a, b, c = sp.Matrix([0, 0]), sp.Matrix([1, 0]), sp.Matrix([0, 1])
    ab, bc, ac = (a+b)/2, (b+c)/2, (a+c)/2
    areas = []
    for p, q, r in ((a, ab, ac), (b, bc, ab), (c, ac, bc), (ab, bc, ac)):
        areas.append(sp.det(sp.Matrix.hstack(q-p, r-p)))
    assert areas == [sp.Rational(1, 4)] * 4
    assert sum(areas) == 1


def test_endpoint_projection_carries_the_fundamental_class_with_unit_degree():
    coarse = incidence_complex(TETRA_BOUNDARY)
    records, projection = refinement_subdivide(TETRA_BOUNDARY)
    fine = incidence_complex(records)
    def sign(triangle):
        return (-1)**sum(triangle[i] > triangle[j] for i in range(3) for j in range(i+1, 3))
    image = {face: 0 for face in coarse.triangles}
    for t in orient(fine):
        projected = tuple(projection[v] for v in t)
        if len(set(projected)) == 3:
            image[frozenset(projected)] += sign(projected)
    degrees = {image[frozenset(t)] * sign(t) for t in orient(coarse)}
    assert degrees in ({1}, {-1})


def test_certificate_handles_lower_dimensional_components():
    coarse = incidence_complex([(0, 1), (2,)])
    fine = incidence_complex([(0, 3), (3, 1), (2,)])
    assert certify_midpoint_subdivision(coarse, fine, {frozenset((0, 1)): 3})
    assert not certify_midpoint_subdivision(coarse, incidence_complex([(0, 3), (3, 1)]),
                                           {frozenset((0, 1)): 3})


def test_relabeling_preserves_high_dimensional_incidence():
    k = incidence_complex([(0, 1, 2, 3), (1, 2, 8)])
    perm = {v: 100 - v for v in k.vertices}
    other = incidence_complex([(100, 99, 98, 97), (99, 98, 92)])
    assert complexes_isomorphic_under(k, other, perm)
    assert euler_characteristic(k) == euler_characteristic(other)


def test_one_defect_repair_all_small_schedules_and_support_dimensions():
    for records in ([(0, 1)], TETRA_BOUNDARY, [(0, 1, 2, 3)]):
        system = RepairSystem(records)
        expected = {v: 0 for v in system.state}
        for order in itertools.permutations(range(len(records))):
            assert system.repair(order) == expected
        if len(records[0]) == 4:
            assert not spherical_incidence_receipt(incidence_complex(records))
        system.state = expected
        assert system.repair() == expected


@pytest.mark.parametrize('order', [[], [0], [0, 0, 1, 2], [0, 1, 2, 4], [True, 0, 2, 3]])
def test_repair_rejects_incomplete_or_aliased_schedules(order):
    with pytest.raises(ValueError, match='schedule'):
        RepairSystem(TETRA_BOUNDARY).repair(order)


def test_repair_does_not_claim_confluence_outside_one_defect_family():
    system = RepairSystem(TETRA_BOUNDARY)
    system.state = {0: 1, 1: 1, 2: 0, 3: 0}
    with pytest.raises(ValueError, match='one seeded defect'):
        system.repair()
    for records, seed in (([], 0), ([(0,)], 0), (TETRA_BOUNDARY, True), (TETRA_BOUNDARY, -1)):
        with pytest.raises(ValueError):
            RepairSystem(records, seed_record=seed)

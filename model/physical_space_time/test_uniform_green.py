"""Independent source, graph, observation and noncollapse controls."""
from fractions import Fraction as F
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp
from scipy.sparse.linalg import splu

from . import uniform_green as u


def test_exact_base_geometry_and_rational_comparison_constants():
    data = u.exact_base_geometry()
    r2 = sp.sympify(data['squared_inradius'])
    a2 = sp.sympify(data['squared_unit_side'])
    assert sp.simplify(8*r2/a2-(7+3*sp.sqrt(5))/3) == 0
    assert r2 > sp.Rational(9, 16)
    # Explicit conservative chain from spherical isoperimetry through degree3.
    cheeger_square = F(3)*F(5, 3)*F(3, 4)**3/20
    assert cheeger_square == F(27, 256)
    assert cheeger_square/6 == F(9, 512)
    assert 24/(cheeger_square/6) == u.GREEN_BOUND
    assert u.GREEN_BOUND/2 < 27**2
    assert F(1+2*27, 132) == F(5, 12)


@pytest.mark.parametrize('level', range(4))
def test_independent_affine_stencil_matches_actual_spherical_source(level):
    assert u.check_affine_geometry(level)['cells'] == 20*4**level
    assert u.source_binding(level)['actual_faces_cross_seams_and_all_ancestries_checked']


def test_subdivision_local_affine_identity_without_a_special_parent_triangle():
    # Arbitrary independent 6-vector coordinates exercise the affine stencil;
    # equality is coefficientwise, rather than a floating geometric test.
    vertices = ((1, 2, 3, 4, 5, 6), (7, 11, 13, 17, 19, 23), (29, 31, 37, 41, 43, 47))
    refined, faces = u.refine(vertices, ((0, 1, 2),))
    assert len(refined) == 6 and len(faces) == 4
    total = tuple(sum(refined[v][j] for face in faces for v in face) for j in range(6))
    assert total == tuple(8*sum(v[j] for v in vertices) for j in range(6))
    # Two glued triangles share one midpoint, including reversed edge order.
    glued, cells = u.refine((*vertices, (53, 59, 61, 67, 71, 73)), ((0, 1, 2), (1, 0, 3)))
    assert len(glued) == 9 and len(cells) == 8
    assert cells[0][1] == cells[4][1]
    # The barycentric stencil divides arbitrary parent area exactly by four.
    corners = ((0, 0), (2, 0), (0, 2), (1, 0), (1, 1), (0, 1))
    for face in faces:
        a, b, c = (corners[i] for i in face)
        assert (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]) == 1


def test_omitted_face_and_false_midpoint_are_rejected(monkeypatch):
    with pytest.raises(ArithmeticError, match='two cells'):
        u.dual_edges(u.FACES[:-1])
    vertices, faces = u.affine_prefix(1)
    changed = list(vertices)
    changed[-1] = tuple(x+int(i == 0) for i, x in enumerate(changed[-1]))
    monkeypatch.setattr(u, 'affine_prefix', lambda level: (tuple(changed), faces))
    with pytest.raises(ArithmeticError, match='dyadic length'):
        u.check_affine_geometry(1)


def test_a_relabelled_source_dictionary_is_not_silently_accepted(monkeypatch):
    builder = u.ancestry.elimination.source._source()
    tower = builder.build_geodesic_icosahedral_tower(1)
    faces = np.array(tower.levels[1].faces, copy=True)
    faces[[0, 1]] = faces[[1, 0]]
    fake = SimpleNamespace(build_geodesic_icosahedral_tower=lambda level:
        SimpleNamespace(levels=(tower.levels[0], SimpleNamespace(faces=faces))))
    monkeypatch.setattr(u.ancestry.elimination.source, '_source', lambda: fake)
    with pytest.raises(ArithmeticError, match='stencil'):
        u.source_binding(1)


@pytest.mark.parametrize('coarse', range(3))
def test_enlarging_readout_dictionary_exact_Gram_and_refinement(coarse):
    assert u.dictionary_control(3, coarse)['coarse_channels'] == 20*4**coarse


def test_fine_and_coarse_complete_native_clock_covariances_obey_same_isometry():
    # Solve all 80 cell observations, independently of the retained20-column
    # certificate; compare Green, attempt and incident outputs after projection.
    raw, laplacian, _, incident, _, _ = u.ancestry.problem(1)
    z = (np.eye(80)-np.ones((80, 80))/80)
    w = np.repeat(z/np.sqrt(12), 12, axis=0)
    solution = np.zeros_like(w)
    solution[1:] = splu(laplacian[1:, 1:].astype(float)).solve(w[1:])
    solution -= solution.mean(axis=0)
    assert np.linalg.norm(laplacian@solution-w) < 2e-11
    green = w.T@solution/4
    pairing = sum(np.outer(w[a]-w[b], solution[a]-solution[b]) for a, b in incident)
    attempt = green-z/(16*raw.seams)
    clock = green+(z-pairing-pairing.T)/(132*4)
    j = np.repeat(np.eye(20)/2, 4, axis=0)
    assert np.array_equal(j.T@j, np.eye(20))
    coarse = u.ancestry.covariance(1)
    for name, matrix in (('green', green), ('attempt', attempt), ('incident', clock)):
        assert np.linalg.norm(j.T@matrix@j-np.asarray(coarse[name], dtype=float), 2) < float(coarse[name+'_error'])+2e-11


@pytest.mark.parametrize('permutation', ((0, 1, 2), (11, 7, 3)))
def test_Green_and_noncollapse_bounds_do_not_require_geometric_port_assignment(permutation):
    # A different valid matching retains the same source graph and native rule.
    local, _ = u.ancestry.elimination.passage.recovered_carrier()
    edges = [(12*c+a, 12*c+b) for c in range(20) for a, b in local]
    used = [0]*20
    for a, b in u.dual_edges(u.FACES):
        edges.append((12*a+permutation[used[a]], 12*b+permutation[used[b]]))
        used[a] += 1; used[b] += 1
    laplacian = np.zeros((240, 240))
    for a, b in edges:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    w = np.repeat((np.eye(20)-np.ones((20, 20))/20)/np.sqrt(12), 12, axis=0)
    green = w.T@np.linalg.solve(laplacian+np.ones((240, 240))/240, w)
    assert np.linalg.eigvalsh(green)[-1] < float(u.GREEN_BOUND)
    root = np.sqrt(5.)
    vertices = np.asarray([[(v[i]+v[i+1]*root)/2 for i in (0, 2, 4)] for v in u.VERTICES])
    centroid = np.asarray([vertices[list(face)].mean(axis=0) for face in u.FACES])
    basis, _ = np.linalg.qr(centroid)
    assert np.linalg.eigvalsh(basis.T@green@basis)[0] > (7+3*np.sqrt(5))/3


def test_degree_and_census_alone_do_not_supply_the_spherical_gap_bound():
    # A3regular prism with the same census has gapO(N^-2), not the source's
    # proved O(N^-1) lower bound. The source refinement hypothesis matters.
    level = 10
    cycle = 10*4**level
    # 2-2cos(2pi/cycle) <= (2pi/cycle)^2 <40/cycle^2.
    assert F(40, cycle**2) < u.bounds(level)['cell_gap_lower']


@pytest.mark.parametrize('level', (0, 3, 20, 100))
def test_all_level_bounds_do_not_require_allocating_a_finite_prefix(level):
    value = u.bounds(level)
    assert value['Green_norm_upper'] == F(4096, 3)
    assert value['incident_Green_difference_upper'] == F(5, 12*4**level)
    assert value['attempt_Green_difference_upper'] == F(1, 2520*16**level)
    assert u.bounds(level+1)['incident_Green_difference_upper'] == value['incident_Green_difference_upper']/4


@pytest.mark.parametrize('value', (True, -1, 1., '2'))
def test_noninteger_level_is_rejected(value):
    with pytest.raises(ValueError):
        u.bounds(value)


def test_incompatible_coarse_dictionary_rejected():
    with pytest.raises(ValueError):
        u.bounds(0, 1)
    with pytest.raises(ValueError):
        u.dictionary_control(1, 1)
    with pytest.raises(ValueError, match='overflow'):
        u.dictionary_control(20, 19)


def test_cached_base_geometry_is_immutable():
    with pytest.raises(TypeError):
        u.exact_base_geometry()['squared_inradius'] = '1'

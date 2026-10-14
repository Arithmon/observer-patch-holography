"""Inverse tests use independently generated hidden quadratic systems."""
import ast
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import sinm, sqrtm

from model.observer_reconstruction import recovery as r


def cube(mass=1.3, coupling=.7):
    sites = np.array(list(product((0, 1), repeat=3)))
    adjacency = (np.sum(abs(sites[:, None]-sites[None, :]), axis=-1) == 1).astype(int)
    return mass**2*np.eye(8)+coupling*(3*np.eye(8)-adjacency), adjacency, sites


def independent_response(matrix, time):
    root = sqrtm(matrix)
    return np.linalg.solve(root, sinm(time*root))


def probabilities(response, gamma=.3, eta=.4):
    # Different detector visibilities are nuisance variables: the inversion
    # must use phases, not a fitted signal amplitude or a supplied vacuum.
    visibility = np.linspace(.3, .9, len(response))[None, :]
    phase = eta*gamma*response.T
    chi = np.array([visibility*np.exp(1j*phase), visibility*np.exp(-1j*phase)])
    return np.stack(((1+chi.real)/2, (1-chi.imag)/2), axis=-1)


def clock(time):
    return [(1+np.cos(time))/2, (1-np.sin(time))/2]


def test_hidden_permutation_recovers_generator_graph_and_geometry():
    matrix, adjacency, sites = cube()
    permutation = np.array([5, 2, 0, 7, 3, 1, 6, 4])
    matrix = matrix[np.ix_(permutation, permutation)]
    adjacency = adjacency[np.ix_(permutation, permutation)]
    first, second = [probabilities(independent_response(matrix, t)) for t in (.2, .4)]
    result = r.recover(first, second, clock(.2), clock(.4), gamma=.3, eta=.4,
                       frequency_upper_bound=3.)
    assert np.allclose(result['stiffness'], matrix, atol=1e-11)
    assert np.array_equal(result['adjacency'], adjacency)
    assert result['mass'] == pytest.approx(1.3, abs=1e-11)
    assert result['coupling'] == pytest.approx(.7, abs=1e-11)
    positions = sites[permutation]-.5
    assert np.allclose(result['gram'], positions@positions.T)
    assert np.allclose(result['coordinates']@result['coordinates'].T, result['gram'])
    assert result['spatial_dimension'] == 3


def test_inverse_is_covariant_under_nonspatial_orthogonal_basis_change():
    matrix, _, _ = cube()
    rotation, _ = np.linalg.qr(np.random.default_rng(15).normal(size=(8, 8)))
    transformed = rotation@matrix@rotation.T
    result = r.infer_stiffness(independent_response(transformed, .2),
                               independent_response(transformed, .4), .2, .4,
                               frequency_upper_bound=3.)
    assert np.allclose(result['stiffness'], transformed, atol=1e-11)
    # A general basis change preserves frequencies, but not local site access.
    with pytest.raises(ValueError, match='cube|couplings'):
        r.cube_geometry(result['stiffness'])


def test_sign_order_and_ramsey_phase_convention():
    response = np.array([[.2, -.03], [-.03, .17]])
    records = probabilities(response)
    assert np.allclose(r.response_from_records(records, gamma=.3, eta=.4), response)
    assert np.allclose(r.response_from_records(records[::-1], gamma=.3, eta=.4), -response)
    assert r.clock_time(clock(.37)) == pytest.approx(.37)
    with pytest.raises(ValueError, match='positive window'):
        r.clock_time(clock(-.37))


def test_predicts_unfitted_time_and_vacuum_wightman():
    matrix, _, _ = cube()
    first, second = [independent_response(matrix, t) for t in (.2, .4)]
    inferred = r.infer_stiffness(first, second, .2, .4, frequency_upper_bound=3.)['stiffness']
    root = sqrtm(matrix)
    covariance = np.linalg.inv(root)/2
    response = independent_response(matrix, .73)
    expected = np.exp(-.4**2*covariance[5, 5]/2+1j*.4*.3*response[5, 2])
    assert r.predict_chi(inferred, .73, 2, 5, gamma=.3, eta=.4) == pytest.approx(expected, abs=1e-13)
    wightman = r.free_wightman(inferred, .73)
    assert np.allclose(wightman.imag, -response/2, atol=1e-12)
    assert np.allclose(r.free_wightman(inferred, 0), covariance, atol=1e-12)
    assert abs(r.predict_chi(inferred, .73, 2, 5, gamma=.3, eta=.4, temperature=1)) < abs(expected)


def test_global_response_scale_cannot_fake_another_canonical_field():
    matrix, _, _ = cube()
    with pytest.raises(ValueError, match='amplitude'):
        r.infer_stiffness(1.1*independent_response(matrix, .2),
                         1.1*independent_response(matrix, .4), .2, .4,
                         frequency_upper_bound=3.)


def test_nonquadratic_or_inconsistent_response_is_rejected():
    matrix, _, _ = cube()
    first = independent_response(matrix, .2)
    second = independent_response(matrix, .4)
    second[0, 0] += .003
    with pytest.raises(ValueError, match='generator|amplitude'):
        r.infer_stiffness(first, second, .2, .4, frequency_upper_bound=3.)


@pytest.mark.parametrize('first,second,t1,t2,upper,match', [
    (np.zeros((2, 2)), np.eye(2), .2, .4, 3., 'singular'),
    (np.eye(2), np.eye(2), .2, .5, 3., 'ratio two'),
    (np.eye(2), np.eye(2), 1., 2., 4., 'nonaliasing'),
    (np.eye(2), 2.1*np.eye(2), .2, .4, 3., 'spectrum'),
    (np.eye(2), np.eye(2), .2, .4, .1, 'frequency band'),
    (np.ones((2, 3)), np.ones((2, 3)), .2, .4, 3., 'square'),
])
def test_rejects_nonidentifiable_or_inadmissible_inputs(first, second, t1, t2, upper, match):
    with pytest.raises(ValueError, match=match):
        r.infer_stiffness(first, second, t1, t2, frequency_upper_bound=upper)


@pytest.mark.parametrize('probabilities', [
    [1., 0.], [-.1, .5], [np.nan, .5], [1+0j, .5+0j], [True, False], [[.5]],
])
def test_invalid_ramsey_probabilities(probabilities):
    with pytest.raises(ValueError):
        r.characteristic(probabilities)


def test_endpoint_roundoff_admission_preserves_raw_records_and_rejects_larger_errors():
    probabilities = np.array([1+5e-13, .5])
    raw = probabilities.copy()
    assert r.characteristic(probabilities) == 1
    assert np.array_equal(probabilities, raw)
    assert r.characteristic([-5e-13, .5]) == -1
    for invalid in ([1+2e-12, .5], [-2e-12, .5]):
        with pytest.raises(ValueError, match='endpoint roundoff'):
            r.characteristic(invalid)


def test_phase_resolution_and_clock_visibility_are_checked():
    with pytest.raises(ValueError, match='visibility'):
        r.response_from_records(np.full((2, 2, 2, 2), .5), gamma=.3, eta=.4)
    with pytest.raises(ValueError, match='visibility'):
        r.clock_time([.6, .4])
    with pytest.raises(ValueError, match='unwrapped'):
        r.response_from_records(np.tile([0., .5], (2, 2, 2, 1)), gamma=.3, eta=.4)


def test_cube_fit_does_not_install_equal_edges_or_mass():
    matrix, _, _ = cube()
    matrix[0, 0] += .05
    with pytest.raises(ValueError, match='cube-Laplacian'):
        r.cube_geometry(matrix)


def test_inverse_module_has_no_executor_or_hamiltonian_dependency():
    tree = ast.parse(Path(r.__file__).read_text())
    modules = {node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)}
    modules.update(alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
                   for alias in node.names)
    assert modules <= {'__future__', 'math', 'numbers', 'numpy'}

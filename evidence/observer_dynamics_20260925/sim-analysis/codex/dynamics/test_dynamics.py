"""Small independent controls and adversarial checks for the dynamics experiment."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("codex_dynamics_run", PATH)
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


def path_laplacian(n):
    adjacency = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    return np.diag(adjacency.sum(axis=1)) - adjacency


def test_heat_trace_and_derivative_tail_bounds_against_full_spectrum():
    eigen = la.eigvalsh(path_laplacian(12))
    rows = RUN.heat_rows(eigen[:4], 12, 4., [.0001, .1, 1., 10., 100.])
    for row in rows:
        s = row["sweeps"] / 2
        ex = np.exp(-s * eigen)
        exact_return = ex.mean()
        dimension = 2 * s * np.dot(eigen, ex) / ex.sum()
        assert row["return_probability_lower"] <= exact_return + 1e-14
        assert row["return_probability_upper"] >= exact_return - 1e-14
        assert row["dimension_lower"] <= dimension + 1e-13
        assert row["dimension_upper"] >= dimension - 1e-13
    assert not rows[0]["truncation_accepted"]  # A few low modes are insufficient at short times.
    assert rows[-1]["truncation_accepted"]


def test_no_missing_modes_collapses_heat_bounds():
    eigen = RUN.torus_eigenvalues(8, 2)
    rows = RUN.heat_rows(eigen, len(eigen), 8., [.2, 3., 10.])
    for row in rows:
        assert row["return_probability_lower"] == row["return_probability_upper"]
        assert row["dimension_lower"] == row["dimension_upper"]
        assert row["tail_over_known_trace_bound"] == 0


def test_torus_spectrum_matches_direct_graph_and_cartesian_product():
    n = 9
    adjacency = np.zeros((n, n))
    for k in range(n):
        adjacency[k, (k + 1) % n] = adjacency[(k + 1) % n, k] = 1
    lap = 2 * np.eye(n) - adjacency
    assert np.allclose(la.eigvalsh(lap), RUN.torus_eigenvalues(n, 1), atol=1e-13)
    product = np.kron(lap, np.eye(n)) + np.kron(np.eye(n), lap)
    assert np.allclose(la.eigvalsh(product), RUN.torus_eigenvalues(n, 2), atol=1e-13)


def test_dimension_controls_distinguish_one_two_and_three():
    controls = RUN.controls()
    at_eight = [next(r for r in c["heat"] if r["sweeps"] == 8.)["dimension_lower"] for c in controls]
    assert at_eight == pytest.approx([1., 2., 3.], abs=.1)
    # Finite volume eventually reduces the dimension: do not extrapolate t→∞ as bulk d.
    assert controls[-1]["heat"][-1]["dimension_lower"] < 1


def test_shift_invert_returns_lowest_modes_not_highest_and_has_small_residuals():
    lap = path_laplacian(40)
    expected = la.eigvalsh(lap)[:5]
    values, vectors = sla.eigsh(sp.csr_matrix(lap), k=5, sigma=-1e-6, which="LM", tol=1e-11)
    assert np.allclose(values, expected, atol=1e-12)
    assert np.max(np.linalg.norm(lap @ vectors - vectors * values, axis=0)) < 1e-10
    # The high spectrum would also have small residuals, so residuals alone do not fix index.
    wrong = la.eigvalsh(lap)[-5:]
    assert not np.allclose(wrong, expected)


def test_static_elimination_includes_coupling_and_preserves_passivity():
    # Two boundary nodes separated by an interior node: effective conductance 1/2.
    lap = path_laplacian(3)
    boundary = [0, 2]
    a = lap[np.ix_(boundary, boundary)]
    b = lap[np.ix_(boundary, [1])]
    c = lap[np.ix_([1], [1])]
    schur = a - b @ la.solve(c, b.T)
    assert np.allclose(schur, [[.5, -.5], [-.5, .5]])
    assert np.allclose(schur @ np.ones(2), 0.)
    assert la.eigvalsh(schur)[0] >= -1e-14
    assert not np.allclose(a, schur)  # Freezing hidden state is not its elimination.


def test_monotone_scalar_clocks_do_not_create_mode_sign_reversal():
    time = np.linspace(0, 10, 201)
    for clock in (time, time ** 2, np.log1p(time)):
        relaxation = np.exp(-clock)
        assert np.all(relaxation > 0)
        assert np.all(np.diff(relaxation) <= 0)
    wave = np.cos(time)
    assert np.any(wave < 0) and np.count_nonzero(np.diff(np.sign(wave))) >= 3
    # Backwards/nonmonotone clocks violate the stated assumption rather than the theorem.
    assert np.any(np.diff(np.exp(-np.sin(time))) > 0)


def test_positive_local_clock_matrix_remains_similar_to_symmetric_generator():
    lap = path_laplacian(7)
    rates = np.array([.3, .6, 1., 2., 4., 1.5, .7])
    generator = -np.diag(rates) @ lap
    symmetric = -np.diag(np.sqrt(rates)) @ lap @ np.diag(np.sqrt(rates))
    assert np.allclose(np.sort(la.eigvals(generator).real), la.eigvalsh(symmetric), atol=1e-12)
    assert np.max(abs(la.eigvals(generator).imag)) < 1e-13
    assert la.eigvalsh(symmetric)[-1] <= 1e-13


def test_reversible_autocorrelation_is_positive_mixture_of_relaxations():
    lap = path_laplacian(7)
    eigen, vectors = la.eigh(lap)
    observable = np.array([1., -2., 3., -1., .5, -1., -.5])
    observable -= observable.mean()
    amplitudes = vectors.T @ observable
    time = np.linspace(0, 20, 31)
    correlation = (amplitudes ** 2)[None, :] * np.exp(-time[:, None] * np.maximum(eigen, 0)[None, :])
    correlation = correlation.sum(axis=1)
    assert np.all(correlation >= 0)
    assert np.all(np.diff(correlation) <= 1e-13)
    assert correlation[0] == pytest.approx(np.dot(observable, observable))


def geometry_module():
    spec = importlib.util.spec_from_file_location("codex_geometry_control", PATH.with_name("geometry_control.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_cotangent_fem_equilateral_triangle_exact_stiffness_and_area():
    geometry = geometry_module()
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [.5, np.sqrt(3) / 2, 0.]])
    faces = np.array([[0, 1, 2]])
    stiffness, mass = geometry.cotangent_fem(vertices, faces)
    expected = (3 * np.eye(3) - np.ones((3, 3))) / (2 * np.sqrt(3))
    assert np.allclose(stiffness.toarray(), expected)
    assert mass.diagonal() == pytest.approx(np.full(3, np.sqrt(3) / 12))
    assert np.allclose(stiffness @ np.ones(3), 0.)
    assert la.eigvalsh(stiffness.toarray(), mass.toarray()) == pytest.approx([0., 6., 6.], abs=1e-13)


def test_fem_metric_scaling_changes_diffusion_spectrum_by_inverse_length_squared():
    geometry = geometry_module()
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [.2, .7, 0.]])
    faces = np.array([[0, 1, 2]])
    stiffness, mass = geometry.cotangent_fem(vertices, faces)
    scaled_stiffness, scaled_mass = geometry.cotangent_fem(vertices * 3, faces)
    assert np.allclose(scaled_stiffness.toarray(), stiffness.toarray())
    assert np.allclose(scaled_mass.toarray(), 9 * mass.toarray())
    first = la.eigvalsh(stiffness.toarray(), mass.toarray())
    second = la.eigvalsh(scaled_stiffness.toarray(), scaled_mass.toarray())
    assert np.allclose(first, 9 * second)


def test_equilateral_angle_deficit_uses_metric_assumption_and_euler_identity():
    # Tetrahedron: four degree3 vertices each carry pi, total4pi. The same
    # total does not imply all triangulations have twelve defects.
    tetra_valence = np.full(4, 3)
    assert np.sum((6 - tetra_valence) * np.pi / 3) == pytest.approx(4 * np.pi)
    # The particular icosahedral refinement has twelve degree5 vertices,
    # all newly introduced vertices degree6.
    ico_refined_valence = np.r_[np.full(12, 5), np.full(150, 6)]
    deficits = (6 - ico_refined_valence) * np.pi / 3
    assert np.count_nonzero(deficits) == 12
    assert deficits.sum() == pytest.approx(4 * np.pi)


def test_exact_lazy_balanced_chain_and_nonlazy_sign_changing_mutation():
    spec = importlib.util.spec_from_file_location("codex_stationary_control", PATH.with_name("stationary_control.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.build()
    assert all(result["exact_checks"].values())
    assert result["nonlazy_adversarial_autocorrelation"] == [str((-1) ** k) for k in range(9)]
    assert result["autocorrelation_float"][0] > result["autocorrelation_float"][-1] > 0

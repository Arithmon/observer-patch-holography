"""Independent quadrature, oscillator-boundary and tensor-order regressions."""
from itertools import product
from math import pi, prod, sqrt

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from .lattice import (
    GalerkinParameters, SevenModeLattice, mode_functions,
    quartic_polynomial, stiffness,
)


def sampled_modes(grid):
    """Evaluate the declared physical smearings without implementation helpers."""
    points = np.asarray(list(product(range(grid), repeat=3)))*2*pi/grid
    columns = [np.ones(len(points))]
    for axis in range(3):
        columns.extend((sqrt(2)*np.cos(points[:, axis]),
                        sqrt(2)*np.sin(points[:, axis])))
    return points, np.column_stack(columns)


def oscillator(d, omega):
    """Extra levels retain all four virtual ladder steps before projection."""
    annihilation = np.diag(np.sqrt(np.arange(1, d+6)), 1)
    q = (annihilation+annihilation.T)/sqrt(2*omega)
    p = 1j*sqrt(omega/2)*(annihilation.T-annihilation)
    return {"q": q[:d, :d], "q2": (q@q)[:d, :d],
            "p2": (p@p)[:d, :d], "q4": np.linalg.matrix_power(q, 4)[:d, :d]}


def embed_by_kron(dimensions, factors):
    result = np.ones((1, 1))
    for mode in reversed(range(7)):
        result = np.kron(result, factors.get(mode, np.eye(dimensions[mode])))
    return result


@pytest.mark.parametrize("grid", [5, 6, 9, 12])
def test_fixed_smearing_quartic_integral_is_exact_without_aliasing(grid):
    points, modes = sampled_modes(grid)
    np.testing.assert_allclose(mode_functions(points), modes/sqrt((2*pi)**3),
                               atol=2e-14, rtol=0)
    np.testing.assert_allclose(modes.T@modes/len(modes), np.eye(7), atol=2e-14)
    rng = np.random.default_rng(713)
    for q in rng.normal(size=(5, 7)):
        expected = np.mean((modes@q)**4)
        assert quartic_polynomial(q) == pytest.approx(expected, rel=2e-14, abs=2e-13)


def test_four_point_grid_aliases_quartic_even_though_modes_are_orthogonal():
    _, modes = sampled_modes(4)
    np.testing.assert_allclose(modes.T@modes/len(modes), np.eye(7), atol=1e-14)
    cosine = np.eye(7)[1]
    assert np.mean((modes@cosine)**4) == pytest.approx(2.)
    assert quartic_polynomial(cosine) == pytest.approx(1.5)
    with pytest.raises(ValueError):
        GalerkinParameters(grid=4)


def test_spatial_refinement_holds_box_mass_coupling_and_band_fixed():
    for grid in (5, 7, 11, None):
        parameters = GalerkinParameters(grid=grid, mass=.7, coupling=.3, omega=1.2)
        assert parameters.mass == .7
        assert parameters.volume == pytest.approx((2*pi)**3)
        assert parameters.physical_coupling == pytest.approx(.3*(2*pi)**3)
        assert parameters.omega == 1.2
        expected = 1. if grid is None else (grid*np.sin(pi/grid)/pi)**2
        assert stiffness(grid) == pytest.approx(expected)
        assert parameters.kappa == pytest.approx(expected)
    # This is convergence of one stencil coefficient on seven fixed modes.
    # It makes no assertion that the retained Fourier band grows to a field.
    errors = [1-stiffness(n) for n in (20, 40)]
    assert 3.98 < errors[0]/errors[1] < 4.01


@pytest.mark.parametrize("dimension", [2, 4])
def test_projected_oscillator_powers_keep_virtual_cutoff_excursions(dimension):
    lattice = SevenModeLattice((dimension,)*7)
    expected = oscillator(dimension, lattice.parameters.omega)
    for name, matrix in expected.items():
        np.testing.assert_allclose(lattice.operators[0][name], matrix, atol=3e-14)
    op = lattice.operators[0]
    assert np.linalg.norm(op["q4"]-np.linalg.matrix_power(op["q"], 4)) > .1
    assert np.linalg.norm(op["q4"]-op["q2"]@op["q2"]) > .1
    defect = op["q2"]-op["q"]@op["q"]
    assert defect[-1, -1] == pytest.approx(dimension/(2*lattice.parameters.omega))


def test_sparse_hamiltonian_matches_independent_spatial_quadrature():
    parameters = GalerkinParameters(grid=7, mass=.8, coupling=.6, omega=1.3)
    dimensions = (4, 2, 2, 2, 2, 2, 2)
    lattice = SevenModeLattice(dimensions, parameters)
    _, modes = sampled_modes(9)
    operators = [oscillator(d, parameters.omega) for d in dimensions]
    expected = np.zeros((prod(dimensions),)*2, dtype=complex)
    kappa = (7*np.sin(pi/7)/pi)**2
    for i, op in enumerate(operators):
        frequency_squared = parameters.mass**2+(kappa if i else 0)
        local = (op["p2"]+frequency_squared*op["q2"])/2
        local += parameters.coupling*np.mean(modes[:, i]**4)*op["q4"]/24
        expected += embed_by_kron(dimensions, {i: local})
    for i in range(7):
        for j in range(i+1, 7):
            coefficient = parameters.coupling*np.mean(modes[:, i]**2*modes[:, j]**2)/4
            assert lattice.pair_coefficient(i, j) == pytest.approx(coefficient)
            expected += coefficient*embed_by_kron(
                dimensions, {i: operators[i]["q2"], j: operators[j]["q2"]})
    assert sparse.isspmatrix_csr(lattice.hamiltonian)
    np.testing.assert_allclose(lattice.hamiltonian.toarray(), expected, atol=4e-13)
    np.testing.assert_allclose(expected, expected.conj().T, atol=2e-14)
    assert np.linalg.eigvalsh(expected)[0] > 0


def test_nonadjacent_reversed_tensor_embedding_preserves_reference_columns():
    lattice = SevenModeLattice((4, 2, 2, 2, 2, 2, 2))
    rng = np.random.default_rng(39)
    left = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
    right = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    # First site supplied is the least-significant local tensor factor.
    local = np.kron(right, left)
    expected = embed_by_kron(lattice.dimensions, {4: left, 0: right})
    state = .03*(rng.normal(size=(lattice.dimension, 3))
                 +1j*rng.normal(size=(lattice.dimension, 3)))
    np.testing.assert_allclose(lattice.embed_local(local, (4, 0)).toarray(), expected)
    np.testing.assert_allclose(lattice.apply_local(state, (4, 0), local), expected@state)
    np.testing.assert_allclose(lattice.apply_local(state, (4, 0), .25*np.eye(8)), .25*state)


def test_evolution_preserves_unnormalized_branch_and_absolute_phase():
    lattice = SevenModeLattice()
    state = .31*lattice.vacuum_state()
    duration = .037
    expected = expm(-1j*duration*lattice.hamiltonian.toarray())@state
    actual = lattice.evolve(state, duration)
    np.testing.assert_allclose(actual, expected, atol=2e-14)
    assert np.linalg.norm(actual) == pytest.approx(.31)
    np.testing.assert_allclose(lattice.evolve(actual, -duration), state, atol=2e-14)


@pytest.mark.parametrize("parameters", [
    {"grid": 4}, {"grid": 5.5}, {"grid": True}, {"mass": 0},
    {"mass": float("nan")}, {"coupling": -1}, {"coupling": True},
    {"omega": 0}, {"omega": float("inf")},
])
def test_invalid_parameters_rejected(parameters):
    with pytest.raises(ValueError):
        GalerkinParameters(**parameters)


@pytest.mark.parametrize("dimensions", [(2,)*6, (2,)*8, (True,)+(2,)*6,
                                          (0,)+(2,)*6, (2.5,)+(2,)*6])
def test_invalid_occupation_codes_rejected(dimensions):
    with pytest.raises(ValueError):
        SevenModeLattice(dimensions)

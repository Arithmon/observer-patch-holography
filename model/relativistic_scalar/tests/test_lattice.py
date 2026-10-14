"""Independent oscillator-boundary, tensor-order and free-action checks."""
from itertools import product
from math import prod, sqrt

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from model.relativistic_scalar.lattice import (
    CUBE_EDGES, CUBE_SITES, ScalarLattice, ScalarParameters,
    free_dispersion, local_operators,
)


def analytic_moment(power, row, column, omega):
    """Closed oscillator matrix elements, independent of truncated products."""
    lo, hi = sorted((row, column))
    if power == 1:
        return sqrt(hi/(2*omega)) if hi-lo == 1 else 0.
    if power == 2:
        if hi == lo:
            return (2*lo+1)/(2*omega)
        return sqrt((lo+1)*(lo+2))/(2*omega) if hi-lo == 2 else 0.
    if power == 4:
        if hi == lo:
            return 3*(2*lo*lo+2*lo+1)/(4*omega**2)
        if hi-lo == 2:
            return (4*lo+6)*sqrt((lo+1)*(lo+2))/(4*omega**2)
        return sqrt(prod(range(lo+1, lo+5)))/(4*omega**2) if hi-lo == 4 else 0.
    raise AssertionError("unsupported moment")


def dense_oracle(lattice):
    """Enumerate occupation transitions using analytic matrix elements."""
    p = lattice.parameters
    states = list(product(*(range(d) for d in reversed(lattice.dimensions))))
    states = [tuple(reversed(state)) for state in states]
    out = np.zeros((lattice.dimension, lattice.dimension))
    for row, left in enumerate(states):
        for column, right in enumerate(states):
            changed = {i for i in range(8) if left[i] != right[i]}
            for site in range(8):
                if changed <= {site}:
                    q2 = analytic_moment(2, left[site], right[site], p.omega)
                    p2 = p.omega**2*q2*(1 if left[site] == right[site] else -1)
                    q4 = analytic_moment(4, left[site], right[site], p.omega)
                    degree = sum(site in edge for edge in CUBE_EDGES)
                    out[row, column] += .5*p2+.5*(p.mass**2+p.kappa*degree)*q2+p.g*q4
            for x, y in CUBE_EDGES:
                if changed <= {x, y}:
                    out[row, column] -= p.kappa*analytic_moment(1, left[x], right[x], p.omega)*analytic_moment(1, left[y], right[y], p.omega)
    return out


def test_open_cube_and_little_endian_occupation_order():
    lattice = ScalarLattice((3, 2, 2, 1, 1, 1, 1, 1))
    assert lattice.sites == CUBE_SITES
    assert lattice.origin == 0 and lattice.sites[0] == (0, 0, 0)
    assert len(CUBE_EDGES) == 12 and np.array_equal(lattice.degree, np.full(8, 3))
    assert {tuple(np.abs(np.subtract(CUBE_SITES[x], CUBE_SITES[y]))) for x, y in CUBE_EDGES} == {(1, 0, 0), (0, 1, 0), (0, 0, 1)}
    state = lattice.basis_state((2, 1, 1, 0, 0, 0, 0, 0))
    assert np.flatnonzero(state).tolist() == [2+3+6]
    assert np.array_equal(lattice.vacuum_state(), lattice.basis_state((0,)*8))


@pytest.mark.parametrize("dimension", [1, 2, 3, 4, 8])
@pytest.mark.parametrize("omega", [.7, 2.])
def test_polynomial_compression_matches_infinite_oscillator_elements(dimension, omega):
    op = local_operators(dimension, omega)
    for power in (1, 2, 4):
        expected = np.array([[analytic_moment(power, m, n, omega)
                              for n in range(dimension)] for m in range(dimension)])
        assert np.allclose(op["q" if power == 1 else f"q{power}"], expected, atol=2e-13)
    boundary = np.zeros((dimension, dimension))
    boundary[-1, -1] = dimension/(2*omega)
    assert np.allclose(op["q2"]-op["q"]@op["q"], boundary)
    assert np.allclose(op["p2"]-op["p"]@op["p"], omega**2*boundary)
    assert not np.allclose(op["q4"], np.linalg.matrix_power(op["q"], 4))
    assert not np.allclose(op["q4"], op["q2"]@op["q2"])
    # No finite-dimensional CCR is asserted: the top state records its defect.
    commutator = op["q"]@op["p"]-op["p"]@op["q"]
    expected = 1j*np.eye(dimension)
    expected[-1, -1] -= 1j*dimension
    assert np.allclose(commutator, expected)


def test_sparse_hamiltonian_matches_independent_dense_moments():
    p = ScalarParameters(mass=.8, speed=.6, spacing=1.3, coupling=.9, omega=1.1)
    lattice = ScalarLattice((4, 3, 1, 1, 1, 1, 1, 1), p)
    assert sparse.isspmatrix_csr(lattice.hamiltonian)
    assert np.allclose(lattice.hamiltonian.toarray(), dense_oracle(lattice), atol=3e-14)


def test_hamiltonian_hermiticity_positivity_and_default_onsite_identity():
    lattice = ScalarLattice((2,)*8)
    dense = lattice.hamiltonian.toarray()
    assert np.max(abs(dense-dense.conj().T)) < 2e-15
    assert np.linalg.eigvalsh(dense)[0] > 0
    assert np.allclose(lattice.onsite(0)-lattice.parameters.g*lattice.operators[0]["q4"],
                       np.diag([1., 3.]))
    # Repulsive interaction is exactly the sum of compressed q^4 terms.
    free = ScalarLattice((2,)*8, ScalarParameters(coupling=0))
    interaction = (lattice.hamiltonian-free.hamiltonian).toarray()
    assert np.linalg.eigvalsh(interaction)[0] > 0


@pytest.mark.parametrize("sites", [(0,), (1, 0), (2, 0), (2, 1, 0)])
def test_local_tensor_application_matches_independent_occupation_enumeration(sites):
    lattice = ScalarLattice((3, 2, 2, 1, 1, 1, 1, 1))
    rng = np.random.default_rng(451)
    size = prod(lattice.dimensions[site] for site in sites)
    gate = .07*(rng.normal(size=(size, size))+1j*rng.normal(size=(size, size)))
    state = rng.normal(size=lattice.dimension)+1j*rng.normal(size=lattice.dimension)
    occupations = [tuple(reversed(x)) for x in product(*(range(d) for d in reversed(lattice.dimensions)))]
    oracle = np.zeros_like(state)
    for row, left in enumerate(occupations):
        for column, right in enumerate(occupations):
            if any(left[site] != right[site] for site in range(8) if site not in sites):
                continue
            r = c = 0
            stride = 1
            for site in sites:
                r += left[site]*stride
                c += right[site]*stride
                stride *= lattice.dimensions[site]
            oracle[row] += gate[r, c]*state[column]
    assert np.allclose(lattice.embed_local(gate, sites)@state, oracle)
    assert np.allclose(lattice.apply_local(state, sites, gate), oracle)
    assert np.allclose(lattice.apply_local(state, sites, sparse.csr_matrix(gate)), oracle)
    columns = np.column_stack((state, .2j*state))
    assert np.allclose(lattice.apply_local(columns, sites, gate), np.column_stack((oracle, .2j*oracle)))
    assert np.allclose(lattice.apply_local(state, sites, .25*np.eye(size)), .25*state)


def test_model_time_evolution_preserves_unnormalized_branch_and_matches_dense():
    lattice = ScalarLattice((3, 2, 1, 1, 1, 1, 1, 1))
    state = .3*lattice.vacuum_state()
    expected = expm(-.17j*dense_oracle(lattice))@state
    actual = lattice.evolve(state, .17)
    assert np.linalg.norm(actual-expected) < 2e-15
    assert np.linalg.norm(actual) == pytest.approx(.3)
    assert np.linalg.norm(lattice.evolve(actual, -.17)-state) < 2e-15
    assert np.array_equal(lattice.evolve(state, 0), state)


def test_open_free_normal_frequencies_and_infinite_lattice_dispersion():
    p = ScalarParameters(mass=.7, speed=1.2, spacing=.8)
    lattice = ScalarLattice((1,)*8, p)
    expected = sorted(sqrt(p.mass**2+2*p.kappa*sum(mode)) for mode in product((0, 1), repeat=3))
    assert np.allclose(lattice.free_frequencies(), expected)
    # Independent periodic nearest-neighbor stencil on a five-site-per-axis
    # torus: a plane wave is an eigenvector of the untruncated free stiffness.
    side = 5
    sites = np.array(list(product(range(side), repeat=3)))
    k = 2*np.pi*np.array([1, 2, -1])/(side*p.spacing)
    field = np.exp(1j*p.spacing*(sites@k))
    action = p.mass**2*field.copy()
    for axis in range(3):
        for direction in (-1, 1):
            shifted = sites.copy()
            shifted[:, axis] = (shifted[:, axis]+direction) % side
            action += p.kappa*(field-np.exp(1j*p.spacing*(shifted@k)))
    omega = free_dispersion(k, p.mass, p.speed, p.spacing)
    assert np.max(abs(action-omega**2*field)) < 3e-13
    assert np.allclose(free_dispersion([k, -k], p.mass, p.speed, p.spacing), [omega, omega])
    assert free_dispersion([0, 0, 0], p.mass, p.speed, p.spacing) == p.mass
    momentum = np.array([.2, -.1, .3])
    continuum = sqrt(p.mass**2+p.speed**2*np.dot(momentum, momentum))
    error1 = continuum-free_dispersion(momentum, p.mass, p.speed, .1)
    error2 = continuum-free_dispersion(momentum, p.mass, p.speed, .05)
    assert 3.99 < error1/error2 < 4.01


@pytest.mark.parametrize("name,value", [
    ("mass", 0), ("mass", -1), ("speed", 0), ("spacing", 0),
    ("omega", 0), ("coupling", -1), ("coupling", True),
    ("omega", np.inf), ("mass", np.nan), ("mass", "1"),
    ("mass", 1e308), ("spacing", 1e-308),
])
def test_invalid_model_parameters_rejected(name, value):
    with pytest.raises(ValueError):
        ScalarParameters(**{name: value})


@pytest.mark.parametrize("dimensions", [
    2, (2,)*7, (2,)*9, (0,)+(2,)*7, (-1,)+(2,)*7,
    (2.5,)+(2,)*7, (True,)+(2,)*7, (np.nan,)+(2,)*7,
])
def test_invalid_cutoff_shapes_rejected(dimensions):
    with pytest.raises(ValueError):
        ScalarLattice(dimensions)


def test_invalid_local_operations_states_and_momenta_rejected():
    lattice = ScalarLattice((2,)*8)
    state = lattice.vacuum_state()
    for sites, matrix in [((0, 0), np.eye(4)), ((8,), np.eye(2)),
                          ((False,), np.eye(2)), ((), np.eye(1)),
                          ((0,), np.eye(3)), ((0,), np.full((2, 2), np.nan))]:
        with pytest.raises(ValueError):
            lattice.apply_local(state, sites, matrix)
    for bad_state in [np.zeros(3), np.full(lattice.dimension, np.inf),
                      np.ma.array(state, mask=False), np.ones(lattice.dimension, dtype=bool)]:
        with pytest.raises(ValueError):
            lattice.evolve(bad_state, .1)
    for bad_time in [True, np.inf, 1j]:
        with pytest.raises(ValueError):
            lattice.evolve(state, bad_time)
    for k in [1, [1, 2], [1, 2, np.nan], [1, 2, 3j], [True, False, True]]:
        with pytest.raises(ValueError):
            free_dispersion(k)
    with pytest.raises(ValueError):
        lattice.basis_state((2,)+(0,)*7)
    with pytest.raises(ValueError):
        ScalarLattice(parameters={"mass": 1})

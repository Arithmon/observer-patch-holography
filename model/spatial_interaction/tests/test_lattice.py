"""Independent particle-coordinate and local-code checks for spatial fields."""
import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.spatial_interaction.lattice import CubicBosonLattice, local_sector_leakage


def particle_coordinate_oracle(lattice, edges=None, interaction=None):
    """Build two labeled particles, then independently impose Bose symmetry."""
    volume = lattice.volume
    one_particle = np.zeros((volume, volume))
    edges = lattice.edges if edges is None else edges
    for x, y in edges:
        difference = np.zeros(volume)
        difference[x], difference[y] = 1, -1
        one_particle += lattice.hopping*np.outer(difference, difference)
    kinetic = sparse.csr_matrix(one_particle)
    ordered = sparse.kron(kinetic, sparse.eye(volume)) + sparse.kron(sparse.eye(volume), kinetic)
    collision = np.zeros(volume*volume)
    collision[np.arange(volume)*(volume+1)] = lattice.interaction if interaction is None else interaction
    ordered += sparse.diags(collision)
    rows, columns, values = [], [], []
    for column, (x, y) in enumerate(lattice.basis):
        rows.append(x*volume+y)
        columns.append(column)
        values.append(1 if x == y else 1/np.sqrt(2))
        if x != y:
            rows.append(y*volume+x)
            columns.append(column)
            values.append(1/np.sqrt(2))
    symmetric = sparse.coo_matrix((values, (rows, columns)),
                                 shape=(volume*volume, lattice.dimension)).tocsr()
    return (symmetric.T @ ordered @ symmetric).tocsr()


def local_code_oracle(lattice, state, edge, gate):
    """Enumerate occupations and spectators without the runtime's blocks."""
    x, y = edge
    output = np.zeros_like(state, dtype=complex)
    for index, pair in enumerate(lattice.basis):
        nx, ny = pair.count(x), pair.count(y)
        spectator = tuple(site for site in pair if site not in edge)
        for target_x in range(nx+ny+1):
            target_y = nx+ny-target_x
            target = tuple(sorted(spectator+(x,)*target_x+(y,)*target_y))
            output[lattice.index[target]] += gate[target_x+4*target_y, nx+4*ny]*state[index]
    return output


@pytest.fixture(scope="module")
def lattice():
    return CubicBosonLattice(1, hopping=.7, interaction=1.3)


def test_genuinely_three_dimensional_boxes_and_consistent_matchings(lattice):
    assert lattice.volume == 27
    assert lattice.dimension == 27*28//2
    assert len(lattice.edges) == 3*2*3**2
    assert len(set(lattice.edges)) == len(lattice.edges)
    for group, edges in enumerate(lattice.edge_groups):
        assert len({site for edge in edges for site in edge}) == 2*len(edges)
        for x, y in edges:
            displacement = np.subtract(lattice.sites[y], lattice.sites[x])
            assert np.array_equal(displacement, np.eye(3, dtype=int)[group//2])
            assert lattice.sites[x][group//2] % 2 == group % 2
    larger = CubicBosonLattice(2)
    for group in range(6):
        inner = {(lattice.sites[x], lattice.sites[y]) for x, y in lattice.edge_groups[group]}
        outer = {(larger.sites[x], larger.sites[y]) for x, y in larger.edge_groups[group]}
        assert inner <= outer


def test_hamiltonian_matches_independent_symmetric_particle_construction(lattice):
    reference = particle_coordinate_oracle(lattice)
    delta = (reference-lattice.hamiltonian).toarray()
    assert np.max(abs(delta)) < 4e-15
    assert np.max(abs((lattice.hamiltonian-lattice.hamiltonian.T).data), initial=0) == 0
    # The interaction is exactly the doublon projector, not a fitted shift.
    free = particle_coordinate_oracle(lattice, interaction=0)
    assert np.allclose((lattice.hamiltonian-free).toarray(),
                       np.diag(lattice.interaction*lattice.collision_mask), atol=4e-15)


def test_local_gate_indexing_and_vacuum_factor_for_kraus_maps(lattice):
    rng = np.random.default_rng(773)
    state = rng.normal(size=lattice.dimension)+1j*rng.normal(size=lattice.dimension)
    state /= np.linalg.norm(state)
    gate = np.zeros((16, 16), complex)
    gate[0, 0] = .3+.2j
    for sector in ((1, 4), (2, 5, 8)):
        trial = rng.normal(size=(len(sector), len(sector)))+1j*rng.normal(size=(len(sector), len(sector)))
        unitary, _ = np.linalg.qr(trial)
        gate[np.ix_(sector, sector)] = .7*unitary
    assert local_sector_leakage(gate) == 0
    expected = local_code_oracle(lattice, state, lattice.ramsey_edge, gate)
    actual = lattice.apply_local(state, lattice.ramsey_edge, gate)
    assert np.allclose(actual, expected, atol=2e-16)
    # Reversing edge orientation reverses which occupation is least significant.
    reverse = tuple(reversed(lattice.ramsey_edge))
    assert np.allclose(lattice.apply_local(state, reverse, gate),
                       local_code_oracle(lattice, state, reverse, gate), atol=2e-16)
    columns = np.column_stack((state, state.conj()))
    assert np.allclose(lattice.apply_local(columns, lattice.ramsey_edge, gate),
                       np.column_stack((actual, local_code_oracle(lattice, state.conj(),
                                                                  lattice.ramsey_edge, gate))))


def test_one_local_evolution_matches_independent_sparse_exponential(lattice):
    lowering = np.diag(np.sqrt(np.arange(1, 4)), 1)
    ax, ay = np.kron(np.eye(4), lowering), np.kron(lowering, np.eye(4))
    local = lattice.hopping*(ax-ay).conj().T @ (ax-ay)
    gate = expm(-.19j*local)
    rng = np.random.default_rng(811)
    state = rng.normal(size=lattice.dimension)+1j*rng.normal(size=lattice.dimension)
    state /= np.linalg.norm(state)
    reference = particle_coordinate_oracle(lattice, edges=(lattice.ramsey_edge,), interaction=0)
    expected = expm_multiply(-.19j*reference, state)
    assert np.linalg.norm(lattice.apply_local(state, lattice.ramsey_edge, gate)-expected) < 8e-15


def test_matching_matches_separate_gates_and_commutes_on_disjoint_edges(lattice):
    lowering = np.diag(np.sqrt(np.arange(1, 4)), 1)
    ax, ay = np.kron(np.eye(4), lowering), np.kron(lowering, np.eye(4))
    gate = expm(-.037j*(ax-ay).conj().T @ (ax-ay))
    expected = lattice.initial_state()
    for edge in lattice.edge_groups[0]:
        expected = lattice.apply_local(expected, edge, gate)
    assert np.allclose(lattice.apply_matching(lattice.initial_state(), 0, gate), expected)
    assert np.allclose(lattice.apply_matching(lattice.initial_state(),
                                              list(reversed(lattice.edge_groups[0])), gate), expected)
    assert np.linalg.norm(expected) == pytest.approx(1, abs=3e-15)


@pytest.mark.parametrize("vacuum", (.8*np.exp(.17j), .3+.1j, 0j))
def test_matching_retains_arbitrary_kraus_vacuum_amplitude(lattice, vacuum):
    gate = np.diag(np.linspace(.1, .8, 16)).astype(complex)
    gate[0, 0] = vacuum
    rng = np.random.default_rng(901)
    state = rng.normal(size=lattice.dimension)+1j*rng.normal(size=lattice.dimension)
    state /= np.linalg.norm(state)
    expected = state.copy()
    for edge in lattice.edge_groups[0]:
        expected = local_code_oracle(lattice, expected, edge, gate)
    actual = lattice.apply_matching(state, 0, gate)
    assert np.linalg.norm(actual-expected) < 3e-16


def test_observer_operations_keep_all_branch_weight(lattice):
    state = lattice.initial_state()
    assert np.linalg.norm(state) == pytest.approx(1)
    assert abs(state[lattice.state_index((0, 0, 0), (1, 0, 0))]) == pytest.approx(1/np.sqrt(2))
    mask = lattice.origin_occupation_mask(threshold=2)
    assert mask.sum() == 1
    first = [np.where(mask == outcome, state, 0) for outcome in (False, True)]
    assert [np.vdot(branch, branch).real for branch in first] == pytest.approx([.5, .5])
    outcomes = []
    for bit, branch in enumerate(first):
        evolved = lattice.reference_evolve(branch, .12)
        fed = lattice.phase_feedback(evolved, bit*.7)
        for last in lattice.ramsey_branches(fed):
            outcomes.append(np.vdot(last, last).real)
    assert sum(outcomes) == pytest.approx(1, abs=3e-14)
    # The Ramsey probability uses the declared positive-i quadrature.
    no, yes = lattice.ramsey_branches(state)
    assert np.vdot(yes, yes).real == pytest.approx(.5)
    assert np.vdot(no, yes) == pytest.approx(0, abs=2e-16)


def test_positive_kinetic_and_uniform_condensate_zero_mode():
    lattice = CubicBosonLattice(1, hopping=1, interaction=0)
    uniform = np.where(lattice.collision_mask, 1., np.sqrt(2.))/lattice.volume
    assert np.linalg.norm(uniform) == pytest.approx(1)
    assert np.linalg.norm(lattice.hamiltonian @ uniform) < 5e-15
    rng = np.random.default_rng(933)
    for _ in range(4):
        state = rng.normal(size=lattice.dimension)+1j*rng.normal(size=lattice.dimension)
        assert np.vdot(state, lattice.hamiltonian @ state).real >= 0


def test_invalid_domains_and_number_changing_local_maps_are_rejected(lattice):
    for radius in (0, -1, 1.5, True):
        with pytest.raises(ValueError):
            CubicBosonLattice(radius)
    for kwargs in ({"hopping": -1}, {"interaction": float("nan")}, {"hopping": True}):
        with pytest.raises(ValueError):
            CubicBosonLattice(1, **kwargs)
    gate = np.eye(16)
    gate[0, 1] = .001
    with pytest.raises(ValueError, match="particle-number"):
        lattice.apply_local(lattice.initial_state(), lattice.ramsey_edge, gate)
    with pytest.raises(ValueError, match="disjoint"):
        lattice.apply_matching(lattice.initial_state(), [lattice.ramsey_edge]*2, np.eye(16))
    with pytest.raises(ValueError):
        lattice.origin_occupation_mask(0)

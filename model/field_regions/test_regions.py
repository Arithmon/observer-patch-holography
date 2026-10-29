"""Independent full-Fock, complete-history and hostile restriction controls."""
from itertools import combinations
import numpy as np
import pytest
from scipy.linalg import block_diag
from ..interacting_field.lattice import TwoFermions, WeylLattice
from ..interacting_field import source
from .regions import (
    regional_basis, restrict_pure, restrict_density, restrict_regional,
    lift_number_observable, conditional_expectation, embed_two_particles,
)


def annihilators(modes):
    result = []
    for p in range(modes):
        op = np.zeros((2**modes, 2**modes), complex)
        for j in range(2**modes):
            if j & (1 << p):
                op[j ^ (1 << p), j] = (-1)**((j & ((1 << p)-1)).bit_count())
        result.append(op)
    return result


def full_fock_split(modes, kept):
    """Independent oracle: construct every reordered Fock basis by CAR creation."""
    a = annihilators(modes)
    outside = [p for p in range(modes) if p not in kept]
    unitary = np.zeros((2**modes, 2**modes), complex)
    for inside_bits in range(2**len(kept)):
        for outside_bits in range(2**len(outside)):
            occupied = [p for j, p in enumerate(kept) if inside_bits & (1 << j)]
            occupied += [p for j, p in enumerate(outside) if outside_bits & (1 << j)]
            vector = np.zeros(2**modes, complex); vector[0] = 1
            for p in reversed(occupied):
                vector = a[p].conj().T @ vector
            unitary[:, inside_bits*2**len(outside)+outside_bits] = vector
    return unitary


def oracle_restrict(space, state, kept):
    if state.ndim == 1:
        state = state[:, None]
    full = np.zeros((2**space.modes, state.shape[1]), complex)
    for j, (p, q) in enumerate(zip(space.first, space.second, strict=True)):
        full[(1 << p)+(1 << q)] = state[j]
    split = full_fock_split(space.modes, kept).conj().T @ full
    split = split.reshape(2**len(kept), 2**(space.modes-len(kept)), state.shape[1])
    density = np.einsum('ier,jes->irjs', split, split.conj()).reshape(
        2**len(kept)*state.shape[1], 2**len(kept)*state.shape[1])
    positions = {p: j for j, p in enumerate(kept)}
    masks = [sum(1 << positions[p] for p in occupied) for occupied in regional_basis(kept)]
    selected = [mask*state.shape[1]+r for mask in masks for r in range(state.shape[1])]
    return density[np.ix_(selected, selected)]


def random_state(space, reference=1, seed=17):
    rng = np.random.default_rng(seed)
    state = rng.normal(size=(space.dimension, reference))+1j*rng.normal(size=(space.dimension, reference))
    return state/np.linalg.norm(state)


@pytest.mark.parametrize('kept', [(), (0,), (1, 4), (0, 2, 5), (0, 1, 2, 3, 4, 5)])
def test_signed_restriction_agrees_with_independent_full_fock_and_complex_reference(kept):
    space = TwoFermions(6); state = random_state(space, 2)
    actual = restrict_pure(space, state, kept)
    oracle = oracle_restrict(space, state, kept)
    np.testing.assert_allclose(actual, oracle, atol=3e-15)
    flat = state.reshape(-1)
    np.testing.assert_allclose(restrict_density(space, np.outer(flat, flat.conj()), kept, 2), oracle, atol=3e-15)
    assert abs(np.trace(actual)-1) < 1e-14
    assert np.linalg.eigvalsh(actual)[0] > -1e-14


def test_nonconsecutive_split_sign_is_observable_in_local_coherence():
    space = TwoFermions(4)
    state = (space.basis_state(0, 1)+space.basis_state(1, 3))/np.sqrt(2)
    reduced = restrict_pure(space, state, (0, 3))
    assert reduced[1, 2] == pytest.approx(-.5)
    # An unsigned occupation split would produce +1/2 and a different read law.
    assert reduced[1, 2].real < 0


def test_region_needs_zero_and_one_particle_sectors_without_postselection():
    space = TwoFermions(4)
    state = (space.basis_state(0, 1)+space.basis_state(0, 2)+space.basis_state(2, 3))/np.sqrt(3)
    reduced = restrict_pure(space, state, (0, 1))
    np.testing.assert_allclose(np.diag(reduced), [1/3, 1/3, 0, 1/3], atol=1e-15)
    assert np.trace(reduced[3:, 3:]).real == pytest.approx(1/3)
    assert np.trace(reduced).real == pytest.approx(1)


def test_nested_trace_maps_compose_with_reference_and_cross_number_coherence():
    space = TwoFermions(6); state = random_state(space, 3)
    large = (0, 1, 3, 4); small = (1, 4)
    first = restrict_pure(space, state, large)
    np.testing.assert_allclose(restrict_regional(first, large, small, 3),
                               restrict_pure(space, state, small), atol=1e-14)
    # The generic regional map also preserves admitted vacuum/single coherences.
    basis = regional_basis(large)
    vector = np.zeros(len(basis)*2, complex)
    vector[0] = 1/np.sqrt(2); vector[2*basis.index((1,))+1] = 1j/np.sqrt(2)
    density = np.outer(vector, vector.conj())
    via = restrict_regional(restrict_regional(density, large, (0, 1, 4), 2), (0, 1, 4), small, 2)
    direct = restrict_regional(density, large, small, 2)
    np.testing.assert_allclose(via, direct, atol=1e-15)
    assert abs(direct[0, 3]) > .49


def number_operator(kept, mode):
    return np.diag([int(mode in occupied) for occupied in regional_basis(kept)])


def test_observable_duality_and_number_preserving_star_algebra_inclusion():
    space = TwoFermions(6); kept = (1, 4, 5); rng = np.random.default_rng(31)
    blocks = []
    for d in (1, 3, 3):
        blocks.append(rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)))
    x = block_diag(*blocks); y = x.conj().T
    lx = lift_number_observable(space, x, kept).toarray()
    ly = lift_number_observable(space, y, kept).toarray()
    np.testing.assert_allclose(lift_number_observable(space, x@y, kept).toarray(), lx@ly, atol=1e-14)
    np.testing.assert_allclose(lift_number_observable(space, y, kept).toarray(), lx.conj().T, atol=1e-14)
    np.testing.assert_array_equal(lift_number_observable(space, np.eye(7), kept).toarray(), np.eye(space.dimension))
    state = random_state(space).reshape(-1); rho = np.outer(state, state.conj())
    assert np.trace(rho@lx) == pytest.approx(np.trace(restrict_density(space, rho, kept)@x))


def test_disjoint_number_preserving_algebras_commute_and_local_read_keeps_remote_state():
    space = TwoFermions(6); left = (0, 3); right = (1, 5)
    x = np.zeros((4, 4), complex); x[1, 2] = 1; x[2, 1] = 1j
    y = np.zeros((4, 4), complex); y[1, 2] = .3j; y[2, 1] = 1
    a = lift_number_observable(space, x, left).toarray()
    b = lift_number_observable(space, y, right).toarray()
    np.testing.assert_array_equal(a@b, b@a)
    state = random_state(space, 2)
    branches = source.read_ramsey(space, state, *left)
    after = sum(restrict_pure(space, branch, right) for branch in branches)
    np.testing.assert_allclose(after, restrict_pure(space, state, right), atol=5e-13)
    assert abs(sum(np.linalg.norm(branch)**2 for branch in branches)-1) < 1e-12


def test_trace_preserving_ambient_expectation_matches_kraus_projection_and_tower():
    space = TwoFermions(5); rng = np.random.default_rng(61)
    matrix = rng.normal(size=(10, 10))+1j*rng.normal(size=(10, 10))
    kept = (1, 3); actual = conditional_expectation(space, matrix, kept)
    # Independent Fock creation oracle identifies each local/environment factor.
    split = full_fock_split(5, kept)
    sector = np.zeros((32, 10))
    for j, (p, q) in enumerate(zip(space.first, space.second, strict=True)):
        sector[(1 << p)+(1 << q), j] = 1
    kraus = []
    for inside_number in range(3):
        imasks = [k for k in range(4) if k.bit_count() == inside_number]
        emasks = [k for k in range(8) if k.bit_count() == 2-inside_number]
        for a in emasks:
            for b in emasks:
                op = np.zeros((32, 32), complex)
                for i in imasks:
                    op[8*i+a, 8*i+b] = 1/np.sqrt(len(emasks))
                kraus.append(sector.T @ split @ op @ split.conj().T @ sector)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(10), atol=1e-15)
    np.testing.assert_allclose(actual, sum(k@matrix@k.conj().T for k in kraus), atol=2e-15)
    np.testing.assert_allclose(conditional_expectation(space, actual, kept), actual, atol=1e-15)
    assert np.trace(actual) == pytest.approx(np.trace(matrix))
    np.testing.assert_allclose(conditional_expectation(space, np.eye(10), kept), np.eye(10), atol=1e-15)
    via = conditional_expectation(space, conditional_expectation(space, matrix, (0, 1, 3)), kept)
    np.testing.assert_allclose(via, actual, atol=1e-15)
    # Complete positivity with an arbitrary complex two-dimensional reference.
    state = random_state(space, 2); rho = np.outer(state.reshape(-1), state.reshape(-1).conj())
    amplified = sum(np.kron(k, np.eye(2))@rho@np.kron(k.conj().T, np.eye(2)) for k in kraus)
    assert np.linalg.eigvalsh(amplified)[0] > -1e-14


def test_conditional_expectation_is_not_a_homomorphic_restriction():
    space = TwoFermions(4)
    h = np.zeros((4, 4)); h[2, 3] = h[3, 2] = 1
    x = space.lift_one_body(h).toarray()
    ex = conditional_expectation(space, x, (0, 1))
    np.testing.assert_array_equal(ex, np.zeros((6, 6)))
    defect = conditional_expectation(space, x@x, (0, 1))-ex@ex
    assert np.linalg.norm(defect, 2) == 1
    assert np.trace(defect) == 4


def test_exterior_square_embedding_composes_and_retains_reference():
    a, b, c = TwoFermions(4), TwoFermions(6), TwoFermions(7)
    first = (4, 0, 3, 1); second = (6, 1, 4, 2, 5, 0)
    state = random_state(a, 2)
    via = embed_two_particles(b, c, second, embed_two_particles(a, b, first, state))
    direct = embed_two_particles(a, c, tuple(second[j] for j in first), state)
    np.testing.assert_array_equal(via, direct)
    np.testing.assert_allclose(via.conj().T@via, state.conj().T@state, atol=1e-15)


@pytest.mark.parametrize('bad', [(1, 0), (0, 0), (False,), (-1,), (4,)])
def test_invalid_regions_rejected(bad):
    space = TwoFermions(4)
    with pytest.raises(ValueError):
        restrict_pure(space, space.basis_state(0, 1), bad)


def test_number_changing_compression_cannot_be_mislabeled_as_a_homomorphism():
    space = TwoFermions(4); pair = np.zeros((4, 4)); pair[3, 0] = 1
    with pytest.raises(ValueError, match='number preservation'):
        lift_number_observable(space, pair, (0, 1))
    # Full CAR: P a†_0 a†_1 P=0 but P a†_0 a†_1 a_1 a_0 P != 0.
    a = annihilators(4); j = np.zeros((16, 6))
    for i, (p, q) in enumerate(zip(space.first, space.second, strict=True)):
        j[(1 << p)+(1 << q), i] = 1
    creation = a[0].conj().T@a[1].conj().T
    np.testing.assert_array_equal(j.T@creation@j, np.zeros((6, 6)))
    assert np.linalg.norm(j.T@creation@creation.conj().T@j, 2) == 1

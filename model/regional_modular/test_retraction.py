"""Independent Kraus, fixed-N CAR and exact boundary-current controls."""
from fractions import Fraction as F
from itertools import combinations
import numpy as np
import pytest
from scipy.linalg import block_diag, expm
from scipy import sparse
from ..interacting_field.lattice import TwoFermions, WeylLattice
from ..field_regions.regions import conditional_expectation
from .retraction import block_expectation, crossing_certificate


def kraus_operators(dimensions, states):
    """Independent spectral Kraus oracle; no partial-trace contraction."""
    dimension = sum(a*b for a, b in dimensions)
    operators = []
    start = 0
    for (a, b), tau in zip(dimensions, states, strict=True):
        values, vectors = np.linalg.eigh(tau)
        for i, value in enumerate(values):
            for j in range(b):
                exterior = np.outer(vectors[:, i], np.eye(b)[j])*np.sqrt(value)
                operator = np.zeros((dimension, dimension), complex)
                operator[start:start+a*b, start:start+a*b] = np.kron(np.eye(a), exterior)
                operators.append(operator)
        start += a*b
    return operators


def random_matrix(dimension, seed=27):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))


def test_exact_actual_field_boundary_current_and_interaction_independence():
    interacting = crossing_certificate()
    free = crossing_certificate(WeylLattice(1, interaction=0))
    assert interacting == free
    assert interacting['modes'] == 54
    assert interacting['N2_dimension'] == 1431
    assert interacting['ordered_crossing_entries'] == 36
    assert F(interacting['one_body_commutator_HS_squared']) == F(9, 64)
    assert F(interacting['full_N2_commutator_HS_squared']) == F(117, 16)
    assert interacting['Gibbs_preserving_regional_retraction_excluded_for_beta_nonzero']
    assert not interacting['ordinary_regional_CP_restriction_excluded']

    lattice = WeylLattice(1)
    kept = set(interacting['kept_modes'])
    numbers = np.array([len(set(pair) & kept)
                        for pair in combinations(range(lattice.modes), 2)])
    number = sparse.diags(numbers, dtype=np.int64)
    matrix = 16*(number@lattice.hamiltonian-lattice.hamiltonian@number)
    assert np.array_equal(matrix.data.real, np.rint(matrix.data.real))
    assert np.array_equal(matrix.data.imag, np.rint(matrix.data.imag))
    independent = sum(int(x.real)**2+int(x.imag)**2 for x in matrix.data)
    assert F(independent, 256) == F(117, 16)
    assert not crossing_certificate(lattice, ())[
        'Gibbs_preserving_regional_retraction_excluded_for_beta_nonzero']
    assert not crossing_certificate(lattice, tuple(range(54)))[
        'Gibbs_preserving_regional_retraction_excluded_for_beta_nonzero']


def test_sector_weighted_product_reference_has_complete_UCP_retraction():
    dimensions = ((1, 1), (2, 2), (1, 1))
    tau = np.array([[.5, .25j], [-.25j, .5]])
    local = np.array([[.5, .125+.125j], [.125-.125j, .5]])
    states = (np.ones((1, 1)), tau, np.ones((1, 1)))
    # Arbitrary number-sector weights; no full-Fock global product is assumed.
    rho = block_diag(np.array([[.25]]), .5*np.kron(local, tau), np.array([[.25]]))
    operators = kraus_operators(dimensions, states)
    actual = lambda x: block_expectation(x, dimensions, states)
    x = random_matrix(6)
    np.testing.assert_allclose(actual(x), sum(k.conj().T@x@k for k in operators), atol=2e-15)
    np.testing.assert_allclose(sum(k.conj().T@k for k in operators), np.eye(6), atol=2e-15)
    np.testing.assert_allclose(sum(k@rho@k.conj().T for k in operators), rho, atol=2e-15)
    np.testing.assert_allclose(actual(actual(x)), actual(x), atol=2e-15)
    assert abs(np.trace(rho@actual(x))-np.trace(rho@x)) < 2e-15
    regional = block_diag(np.array([[2j]]), np.kron(random_matrix(2), np.eye(2)), np.array([[3]]))
    np.testing.assert_allclose(actual(regional), regional, atol=2e-15)
    np.testing.assert_allclose(actual(regional@x), regional@actual(x), atol=2e-15)
    # A nontracial exterior state need not preserve ordinary ambient trace.
    probe = block_diag(np.zeros((1, 1)), np.kron(np.eye(2), tau), np.zeros((1, 1)))
    assert abs(np.trace(actual(probe))-np.trace(probe)) > .1
    # Complete CP amplification, including a complex untouched reference.
    vector = np.arange(18)+1j*np.arange(18)[::-1]
    joint = np.outer(vector, vector.conj())
    out = sum(np.kron(k.conj().T, np.eye(3))@joint@np.kron(k, np.eye(3)) for k in operators)
    assert np.linalg.eigvalsh(out)[0] > -1e-11


def test_tracial_block_map_matches_actual_nonconsecutive_CAR_expectation():
    space = TwoFermions(6)
    kept, outside = (0, 2, 5), (1, 3, 4)
    basis = tuple(combinations(range(6), 2))
    order = []
    dimensions = []
    for n in range(3):
        local = tuple(combinations(kept, n))
        exterior = tuple(combinations(outside, 2-n))
        dimensions.append((len(local), len(exterior)))
        for a in local:
            for b in exterior:
                words = a+b
                sign = (-1)**sum(words[i] > words[j]
                                  for i in range(2) for j in range(i+1, 2))
                order.append((basis.index(tuple(sorted(words))), sign))
    split = np.zeros((15, 15), int)
    for row, (column, sign) in enumerate(order):
        split[row, column] = sign
    x = random_matrix(15)
    states = tuple(np.eye(b)/b for _, b in dimensions)
    expected = split.T@block_expectation(split@x@split.T, dimensions, states)@split
    np.testing.assert_allclose(conditional_expectation(space, x, kept), expected, atol=2e-15)


def test_cut_action_positive_control_and_boundary_hopping_negative_control():
    # Canonical N2 block order: empty region, one in each factor, full region.
    dimensions = ((1, 1), (2, 2), (1, 1))
    z = np.diag([1., -1.]); x = np.array([[0., 1.], [1., 0.]])
    h = block_diag(np.array([[.25]]), np.kron(z, np.eye(2))+np.kron(np.eye(2), x)/2,
                   np.array([[.75]]))
    rho = expm(-h); rho /= np.trace(rho)
    tau = expm(-x/2); tau /= np.trace(tau)
    states = (np.ones((1, 1)), tau, np.ones((1, 1)))
    e = lambda a: block_expectation(a, dimensions, states)
    a = random_matrix(6)
    assert abs(np.trace(rho@e(a))-np.trace(rho@a)) < 2e-15
    np.testing.assert_allclose(e(h@a-a@h), h@e(a)-e(a)@h, atol=2e-15)
    # A crossing hopping matrix element changes local number.
    crossing = h.copy(); crossing[0, 2] = crossing[2, 0] = 1/8
    changed = expm(-crossing); changed /= np.trace(changed)
    number = np.diag([0., 1., 1., 1., 1., 2.])
    assert np.linalg.norm(crossing@number-number@crossing) > .1
    assert np.linalg.norm(changed@number-number@changed) > .01
    cross_observable = np.zeros((6, 6)); cross_observable[0, 2] = cross_observable[2, 0] = 1
    assert np.count_nonzero(e(cross_observable)) == 0
    assert abs(np.trace(changed@cross_observable)) > .01


def test_central_commutation_is_necessary_but_not_sufficient():
    # Correlation inside a single factor block: the center is automatically fixed.
    rho = np.diag([3., 1., 1., 3.])/8
    tau = np.eye(2)/2
    product = np.eye(4)/4
    assert np.linalg.det(np.diag(rho).reshape(2, 2)) != 0
    zz = np.diag([1., -1., -1., 1.])
    assert np.trace(rho@zz) == .5
    assert np.trace(rho@block_expectation(zz, ((2, 2),), (tau,))) == 0
    operators = kraus_operators(((2, 2),), (tau,))
    np.testing.assert_allclose(sum(k@rho@k.conj().T for k in operators), product, atol=2e-15)


def test_tracial_generator_covariance_does_not_require_Gibbs_preservation():
    local = np.array([[0., 1.], [1., 0.]])
    exterior = np.diag([1., -1.])
    h = np.kron(local, np.eye(2))+np.kron(np.eye(2), exterior)/2
    rho = expm(-h); rho /= np.trace(rho)
    e = lambda a: block_expectation(a, ((2, 2),), (np.eye(2)/2,))
    a = random_matrix(4)
    np.testing.assert_allclose(e(h@a-a@h), h@e(a)-e(a)@h, atol=2e-15)
    probe = np.kron(np.eye(2), exterior)
    assert abs(np.trace(rho@e(probe))-np.trace(rho@probe)) > .4


def test_zero_beta_is_an_explicit_exception_to_the_Gibbs_obstruction():
    rho = np.eye(4)/4
    x = random_matrix(4)
    actual = block_expectation(x, ((2, 2),), (np.eye(2)/2,))
    assert abs(np.trace(rho@actual)-np.trace(rho@x)) < 2e-15
    # Infinite-temperature preservation imposes no condition on a chosen H.
    h = np.kron(np.array([[0., 1.], [1., 0.]]), np.diag([1., -1.]))
    e = lambda a: block_expectation(a, ((2, 2),), (np.eye(2)/2,))
    assert np.linalg.norm(e(h@x-x@h)-(h@e(x)-e(x)@h)) > 1


def test_common_reference_gives_weighted_tower_and_mismatched_reference_does_not():
    tau_b = np.diag([.75, .25]); tau_c = np.diag([.625, .375])
    smaller = lambda x: block_expectation(x, ((2, 4),), (np.kron(tau_b, tau_c),))
    larger = lambda x: block_expectation(x, ((4, 2),), (tau_c,))
    x = random_matrix(8)
    np.testing.assert_allclose(smaller(larger(x)), smaller(x), atol=2e-15)
    np.testing.assert_allclose(larger(smaller(x)), smaller(x), atol=2e-15)
    inconsistent = lambda x: block_expectation(x, ((4, 2),), (np.eye(2)/2,))
    probe = np.kron(np.eye(4), np.diag([1., -1.]))
    assert np.linalg.norm(smaller(inconsistent(probe))-smaller(probe)) > .5


def test_nonfaithful_operator_system_projection_explains_the_faithfulness_hypothesis():
    def e(x):
        return block_diag(x[:2, :2], np.array([[np.trace(x[:2, :2])/2]]))
    operators = [np.diag([1., 1., 0.])]
    for i in (0, 1):
        k = np.zeros((3, 3)); k[i, 2] = 1/np.sqrt(2); operators.append(k)
    rho = np.diag([.5, .5, 0.])
    matrix = random_matrix(3)
    np.testing.assert_allclose(e(matrix), sum(k.conj().T@matrix@k for k in operators), atol=1e-15)
    np.testing.assert_allclose(e(e(matrix)), e(matrix), atol=1e-15)
    assert np.array_equal(e(np.eye(3)), np.eye(3))
    assert np.trace(rho@e(matrix)) == np.trace(rho@matrix)
    fixed = np.diag([1., -1., 0.])
    assert np.array_equal(e(fixed), fixed)
    assert not np.array_equal(e(fixed@fixed), fixed@fixed)


@pytest.mark.parametrize('damage', ['unfaithful', 'unnormalized', 'nonhermitian', 'masked', 'missing'])
def test_declared_block_states_reject_invalid_inputs(damage):
    tau = np.eye(2)/2
    if damage == 'unfaithful': tau = np.diag([1., 0.])
    elif damage == 'unnormalized': tau *= 2
    elif damage == 'nonhermitian': tau[0, 1] = .125
    elif damage == 'masked': tau = np.ma.array(tau, mask=False)
    states = () if damage == 'missing' else (tau,)
    with pytest.raises(ValueError): block_expectation(np.eye(4), ((2, 2),), states)

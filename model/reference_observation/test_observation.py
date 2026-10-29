import numpy as np
import pytest
from scipy.linalg import expm

from .observation import ReferenceObservation, car_restriction
from ..field_regions.regions import restrict_density, restrict_regional
from ..interacting_field.lattice import TwoFermions, WeylLattice


def tensor_trace(local, exterior):
    return tuple(np.kron(np.eye(local), np.eye(exterior)[j:j+1])
                 for j in range(exterior))


def random_density(size, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
    rho = a@a.conj().T+np.eye(size)
    return rho/np.trace(rho)


def test_general_channel_tower_and_reference_adjoint_without_factorization():
    rho = random_density(8, 721)
    ksb, kba = tensor_trace(4, 2), tensor_trace(2, 2)
    direct = ReferenceObservation(rho, tuple(b@s for b in kba for s in ksb))
    first = ReferenceObservation(rho, ksb)
    second = ReferenceObservation(first.target_reference, kba)
    x = random_density(16, 331)
    np.testing.assert_allclose(second.observe(first.observe(x, 2), 2),
                               direct.observe(x, 2), atol=3e-13)
    np.testing.assert_allclose(first.recover(second.recover(random_density(4, 88), 2), 2),
                               direct.recover(random_density(4, 88), 2), atol=3e-13)
    np.testing.assert_allclose(first.recover(first.target_reference), rho, atol=3e-14)
    wrong = ReferenceObservation(np.diag([1, 2, 3, 4])/10, kba)
    assert np.linalg.norm(wrong.observe(first.observe(x, 2), 2)-direct.observe(x, 2)) > 1e-4


def test_complete_reference_valued_povm_and_pairing():
    rho = random_density(6, 43)
    channel = ReferenceObservation(rho, tensor_trace(2, 3))
    effect = random_density(12, 45)
    a, b = channel.observe(effect, 2), channel.observe(np.eye(12)-effect, 2)
    assert np.linalg.eigvalsh(a)[0] > 0 and np.linalg.eigvalsh(b)[0] > 0
    np.testing.assert_allclose(a+b, np.eye(4), atol=3e-14)
    ref = random_density(2, 44)
    assert abs(np.trace(np.kron(channel.target_reference, ref)@a)
               -np.trace(np.kron(rho, ref)@effect)) < 3e-14
    recovered = channel.recover(random_density(4, 44), 2)
    assert abs(np.trace(recovered)-1) < 3e-14
    assert np.linalg.eigvalsh(recovered)[0] > 0


def test_correlated_reference_is_not_a_retraction_or_modular_intertwiner():
    rho = np.diag([3, 1, 1, 3])/8
    channel = ReferenceObservation(rho, tensor_trace(2, 2))
    local = np.array([[0, 1], [1, 0]], complex)
    lifted = np.kron(local, np.eye(2))
    observed = channel.observe(lifted)
    np.testing.assert_allclose(observed, np.sqrt(3)/2*local, atol=1e-14)
    assert np.linalg.norm(observed-local) > .1
    u = np.diag(np.exp(.7j*np.log(np.diag(rho))))
    assert np.linalg.norm(channel.observe(u@lifted@u.conj().T)-observed) > .1


def test_product_reference_recovers_weighted_expectation_and_modular_flow():
    a, b = random_density(2, 60), random_density(3, 61)
    rho = np.kron(a, b)
    channel = ReferenceObservation(rho, tensor_trace(2, 3))
    x = random_density(6, 62)
    weighted = np.einsum('aibj,ji->ab', x.reshape(2, 3, 2, 3), b)
    np.testing.assert_allclose(channel.observe(x), weighted, atol=2e-14)
    np.testing.assert_allclose(channel.observe(np.kron(a, np.eye(3))), a, atol=2e-14)
    def flow(r):
        w, v = np.linalg.eigh(r)
        return (v*np.exp(.9j*np.log(w)))@v.conj().T
    u, v = flow(rho), flow(a)
    np.testing.assert_allclose(channel.observe(u@x@u.conj().T),
                               v@weighted@v.conj().T, atol=3e-14)


def test_same_weyl_contact_action_signed_car_nested_tower():
    parent = WeylLattice(1)
    sites = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    modes = [2*parent.site_index[x]+spin for x in sites for spin in (0, 1)]
    space = TwoFermions(8)
    h = space.lift_one_body(parent.single_h[np.ix_(modes, modes)]).toarray()
    for p in range(0, 8, 2):
        h[space.index(p, p+1), space.index(p, p+1)] += 1
    rho = expm(-.7*h)
    rho /= np.trace(rho)
    b, a = (0, 1, 4, 5), (0, 1)
    first = ReferenceObservation(rho, car_restriction(tuple(range(8)), b))
    second = ReferenceObservation(first.target_reference, car_restriction(b, a, regional_input=True))
    direct = ReferenceObservation(rho, car_restriction(tuple(range(8)), a))
    np.testing.assert_allclose(first.target_reference, restrict_density(space, rho, b), atol=2e-15)
    np.testing.assert_allclose(second.target_reference,
                               restrict_regional(first.target_reference, b, a), atol=2e-15)
    x = random_density(56, 811)
    np.testing.assert_allclose(second.observe(first.observe(x, 2), 2),
                               direct.observe(x, 2), atol=2e-13)
    assert np.linalg.norm(rho-np.diag(np.diag(rho))) > .005
    assert first.source_dimension == 28 and first.target_dimension == 11


@pytest.mark.parametrize('reference,restriction', [
    (np.diag([1., 0.]), [np.eye(2)]),
    (np.eye(2), [np.eye(2)]),
    (np.eye(2)/2, [np.eye(2)/2]),
    (np.eye(2)/2, [np.array([[1., 0.], [0., np.nan]])]),
    (np.eye(2)/2, [np.array([[1., 0.], [0., 1.], [0., 0.]])]),
    (np.ones((2, 2), bool), [np.eye(2)]),
    (np.eye(2)/2, []),
])
def test_invalid_reference_support_or_incomplete_channel_fails(reference, restriction):
    with pytest.raises(ValueError):
        ReferenceObservation(reference, restriction)


def test_immutable_arrays_and_exact_arbitrary_operator_domain():
    rho = np.eye(4)/4
    channel = ReferenceObservation(rho, tensor_trace(2, 2))
    rho[0, 0] = 3
    with pytest.raises(ValueError):
        channel.kraus[0].setflags(write=True)
    x = np.zeros((4, 4), complex)
    x[0, 2] = -3+2j
    expected = np.array([[0, (-3+2j)/2], [0, 0]])
    np.testing.assert_allclose(channel.observe(x), expected)
    with pytest.raises(ValueError):
        channel.observe(x, True)
    with pytest.raises(ValueError):
        channel.recover(np.eye(2), 0)


def test_zero_support_car_region_requires_explicit_support_choice():
    with pytest.raises(ValueError, match='faithful'):
        ReferenceObservation(np.eye(6)/6, car_restriction((0, 1, 2, 3), (0, 1, 2, 3)))


def test_small_eigenvalue_does_not_impose_an_unstated_support_cutoff():
    rho = np.diag([1e-24, 1-1e-24])
    channel = ReferenceObservation(rho, [np.eye(2)])
    np.testing.assert_allclose(channel.observe(np.eye(2)), np.eye(2), atol=1e-14)
    np.testing.assert_allclose(channel.recover(rho), rho, atol=1e-14)


def test_ill_conditioned_rotated_reference_cannot_silently_lose_unitality():
    rotation = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                         [1, 1, -1, -1], [1, -1, -1, 1.]])/2
    rho = (rotation*np.array([1e-15, .2, .3, .5-1e-15]))@rotation.T
    # The exact identity observation is unital. Its floating whitening has
    # amplified cancellation error and must fail admission, not clip rho.
    with pytest.raises(ArithmeticError, match='residual budget'):
        ReferenceObservation(rho, [np.eye(4)])

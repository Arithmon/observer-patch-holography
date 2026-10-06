"""Independent controls for support, channels, and observable recovery errors."""

from fractions import Fraction
import itertools

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from quantum_information.recovery import (
    matrix_inv_sqrt_psd, matrix_sqrt_psd, petz_recovery, petz_recovery_kraus,
    state_fidelity, trace_distance,
)
from quantum_information.states import conditional_mutual_information, partial_trace


def apply(kraus, value):
    return sum(k @ value @ k.conj().T for k in kraus)


@pytest.mark.parametrize("dims", ((2, 2), (2, 3), (3, 2)))
@pytest.mark.parametrize("seed", range(4))
def test_full_channel_and_reference_identities_for_noncommuting_data(dims, seed):
    b, c = dims
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(b*c, b*c)) + 1j*rng.normal(size=(b*c, b*c))
    tau = x @ x.conj().T + np.eye(b*c)
    tau /= np.trace(tau)
    kraus = petz_recovery_kraus(tau, dims)
    np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(b), atol=8e-15)
    marginal = partial_trace(tau, dims, [0])
    np.testing.assert_allclose(apply(kraus, marginal), tau, atol=8e-15)
    # Assemble the full Choi matrix by applying the map to every matrix unit.
    units = [np.eye(b)[:, i:i+1] @ np.eye(b)[j:j+1] for i in range(b) for j in range(b)]
    choi = sum(np.kron(unit, apply(kraus, unit)) for unit in units)
    assert np.linalg.eigvalsh(choi)[0] > -2e-15
    np.testing.assert_allclose(np.trace(choi.reshape(b,b*c,b,b*c), axis1=1, axis2=3),
                               np.eye(b), atol=8e-15)
    # Independent high-precision matrix-function evaluation of the Petz formula.
    probe = np.zeros((b,b), complex)
    probe[0,1], probe[1,0] = 1j, -1j
    with mp.workdps(65):
        root = mp.sqrtm(mp.matrix(tau.tolist()))
        inverse = mp.sqrtm(mp.matrix(marginal.tolist()))**-1
        small = inverse * mp.matrix(probe.tolist()) * inverse
        lifted = mp.matrix(b*c)
        for i, j, k in itertools.product(range(b), range(b), range(c)):
            lifted[i*c+k, j*c+k] = small[i,j]
        expected = np.array((root*lifted*root).tolist(), complex)
    np.testing.assert_allclose(apply(kraus, probe), expected, atol=6e-15, rtol=2e-13)


def test_singular_reference_has_explicit_trace_preserving_completion():
    tau = np.diag([.75, .25, 0, 0])
    kraus = petz_recovery_kraus(tau, (2,2))
    # The usual support-restricted Petz expression sends |1><1| to zero.
    # The declared full channel instead sends it to tau, without normalization.
    for state in (np.diag([1.,0]), np.diag([0.,1]), np.array([[.5,.5j],[-.5j,.5]])):
        np.testing.assert_allclose(apply(kraus, state), tau, atol=3e-16)
    mixture = .3*np.diag([1.,0])+.7*np.diag([0.,1])
    np.testing.assert_allclose(apply(kraus, mixture),
                               .3*apply(kraus, np.diag([1.,0]))+.7*apply(kraus, np.diag([0.,1])),
                               atol=3e-16)
    np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(2), atol=3e-16)


def test_symbolic_bell_reference_channel_is_replacement_not_an_identity():
    # Exact K_j = sqrt(2) P_Bell (I tensor |j>), so R(X)=Tr(X)P_Bell.
    psi = sp.Matrix([1,0,0,1])/sp.sqrt(2)
    tau = psi*psi.T
    kraus = [sp.sqrt(2)*tau*sp.kronecker_product(sp.eye(2), sp.eye(2)[:,j]) for j in range(2)]
    assert sum((k.T*k for k in kraus), sp.zeros(2)) == sp.eye(2)
    x = sp.Matrix(2,2,sp.symbols("x:4"))
    assert sp.simplify(sum((k*x*k.T for k in kraus), sp.zeros(4))-sp.trace(x)*tau) == sp.zeros(4)
    actual = petz_recovery_kraus(np.array(tau).astype(complex), (2,2))
    # Keeping an entangled reference can necessarily disturb the input B state.
    np.testing.assert_allclose(apply(actual, np.diag([1.,0])), np.array(tau).astype(complex), atol=5e-16)


def test_classical_recovery_and_information_identity_on_all_255_supports():
    labels = list(itertools.product(range(2), repeat=3))
    for mask in range(1,256):
        weights = [(i+1) if mask & (1<<i) else 0 for i in range(8)]
        p = [Fraction(w,sum(weights)) for w in weights]
        ab = {(a,b): sum(p[i] for i,t in enumerate(labels) if t[:2] == (a,b))
              for a,b in itertools.product(range(2),repeat=2)}
        bc = {(b,c): sum(p[i] for i,t in enumerate(labels) if t[1:] == (b,c))
              for b,c in itertools.product(range(2),repeat=2)}
        pb = {b: sum(p[i] for i,t in enumerate(labels) if t[1] == b) for b in range(2)}
        q = [ab[a,b]*bc[b,c]/pb[b] if pb[b] else Fraction(0) for a,b,c in labels]
        assert sum(q) == 1
        recovered = petz_recovery(np.diag([float(x) for x in p]))
        np.testing.assert_allclose(recovered, np.diag([float(x) for x in q]), atol=3e-16, rtol=3e-14)
        with mp.workdps(65):
            exact = lambda v: mp.mpf(v.numerator)/v.denominator
            entropy = lambda vs: -mp.fsum(exact(v)*mp.log(exact(v)) for v in vs if v)
            divergence = mp.fsum(exact(pi)*mp.log(exact(pi)/exact(qi)) for pi,qi in zip(p,q) if pi)
            cmi = entropy(ab.values())+entropy(bc.values())-entropy(pb.values())-entropy(p)
            assert abs(cmi-divergence) < mp.mpf("1e-60")


@pytest.mark.parametrize("seed", range(8))
def test_quantum_recovery_preserves_a_and_bc_but_is_not_a_perfect_repair(seed):
    rng = np.random.default_rng(seed+300)
    x = rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    rho = x @ x.conj().T+np.eye(8)
    rho /= np.trace(rho)
    recovered = petz_recovery(rho)
    np.testing.assert_allclose(partial_trace(recovered,[2,2,2],[0]), partial_trace(rho,[2,2,2],[0]), atol=3e-15)
    np.testing.assert_allclose(partial_trace(recovered,[2,2,2],[1,2]), partial_trace(rho,[2,2,2],[1,2]), atol=3e-15)
    assert trace_distance(rho,recovered) > .01
    # Do not claim that quantum Petz preserves AB, or makes every input Markov.
    assert np.linalg.norm(partial_trace(recovered,[2,2,2],[0,1])-partial_trace(rho,[2,2,2],[0,1])) > 1e-4
    assert conditional_mutual_information(recovered,[2,2,2],[0],[1],[2]) > 1e-6
    assert trace_distance(petz_recovery(recovered),recovered) > 1e-4


@pytest.mark.parametrize("seed", range(12))
def test_squared_fidelity_matches_independent_qubit_determinant_law(seed):
    rng = np.random.default_rng(seed)
    states = []
    for _ in range(2):
        x = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        a = x @ x.conj().T+np.eye(2)
        states.append(a/np.trace(a))
    a,b = states
    expected = np.trace(a@b).real+2*np.sqrt(np.linalg.det(a).real*np.linalg.det(b).real)
    assert state_fidelity(a,b) == pytest.approx(expected, abs=2e-15)
    assert state_fidelity(a,b) == pytest.approx(state_fidelity(b,a), abs=2e-15)


@pytest.mark.parametrize("mass", (2.**-40,1e-100,1e-310,np.nextafter(0.,1.)))
def test_root_and_inverse_retain_coordinate_support(mass):
    state = np.diag([1-mass,mass])
    root = matrix_sqrt_psd(state)
    inverse = matrix_inv_sqrt_psd(state)
    assert root[1,1].real > 0
    assert inverse[1,1].real > 1
    np.testing.assert_allclose(root @ inverse, np.eye(2), rtol=3e-16, atol=0)


def test_unresolved_dense_support_is_refused_instead_of_thresholded():
    epsilon = 2.**-50
    rho = np.array([[.5,.5-epsilon],[.5-epsilon,.5]])
    with pytest.raises(ValueError, match="support is numerically unresolved"):
        matrix_inv_sqrt_psd(rho)
    with pytest.raises(ValueError, match="support is numerically unresolved"):
        petz_recovery_kraus(np.kron(rho,np.eye(2)/2), (2,2))


@pytest.mark.parametrize("cutoff", (1e-10, -1., True, [0.], float("nan")))
def test_support_cutoff_cannot_disable_a_positive_direction(cutoff):
    with pytest.raises(ValueError):
        matrix_inv_sqrt_psd(np.eye(2)/2, cutoff=cutoff)


@pytest.mark.parametrize("bad", (np.zeros((2,2)), np.eye(2), np.diag([1.1,-.1]),
                                [[True,0.],[0.,0.]], np.ma.array(np.eye(2)/2,mask=True)))
def test_shared_recovery_entry_points_reject_invalid_evidence(bad):
    for function in (matrix_sqrt_psd, matrix_inv_sqrt_psd):
        with pytest.raises(ValueError):
            function(bad)
    for function in (state_fidelity, trace_distance):
        with pytest.raises(ValueError):
            function(bad,np.eye(2)/2)
        with pytest.raises(ValueError):
            function(np.eye(2)/2,bad)
    with pytest.raises(ValueError):
        petz_recovery_kraus(bad,(2,1))

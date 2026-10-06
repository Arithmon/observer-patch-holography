"""Independent exact spectra, product-state witnesses and Gibbs stability."""

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm, logm

from maxent.information_projection import gibbs_state


I = np.eye(2)
PAULI = (np.array([[0., 1.], [1., 0.]]),
         np.array([[0., -1j], [1j, 0.]]), np.diag([1., -1.]))
H = sum(np.kron(a, a) for a in PAULI)


@pytest.mark.parametrize('beta', (-2., -.5, 0., .1, .25, .3, .5, 1., 2., 4.))
def test_complete_thermal_family_has_the_correct_entanglement_and_response(beta):
    rho, log_z = gibbs_state([H]*3, [2.**60, beta, -2.**60])
    p = 1/(1+3*np.exp(-4*beta))
    expected_spectrum = np.sort([p]+[(1-p)/3]*3)
    np.testing.assert_allclose(np.linalg.eigvalsh(rho), expected_spectrum, atol=2e-15, rtol=0)
    pt = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    pt_spectrum = np.sort([.5-p]+[(1+2*p)/6]*3)
    np.testing.assert_allclose(np.linalg.eigvalsh(pt), pt_spectrum, atol=2e-15, rtol=0)
    assert log_z == pytest.approx(3*beta+np.log1p(3*np.exp(-4*beta)), abs=3e-15)
    for a in PAULI:
        assert np.trace(rho@np.kron(a, a)).real == pytest.approx((1-4*p)/3, abs=3e-15)
    # Both one-site states stay maximally mixed although correlations vary.
    np.testing.assert_allclose(np.trace(rho.reshape(2,2,2,2), axis1=1, axis2=3), I/2, atol=2e-15)


def test_exact_projectors_prove_the_entire_separability_boundary():
    identity = sp.eye(2)
    pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
    h = sum((sp.kronecker_product(a, a) for a in pauli), sp.zeros(4))
    singlet = (sp.eye(4)-h)/4
    assert singlet*singlet == singlet and sp.trace(singlet) == 1
    assert h.eigenvals() == {sp.Integer(-3): 1, sp.Integer(1): 3}
    aligned, opposite = sp.zeros(4), sp.zeros(4)
    for a in pauli:
        for sign in (-1, 1):
            first, second = (identity+sign*a)/2, (identity-sign*a)/2
            assert first*first == first and sp.trace(first) == 1
            aligned += sp.kronecker_product(first, first)/6
            opposite += sp.kronecker_product(first, second)/6
    p = sp.Symbol('p', real=True)
    state = p*singlet+(1-p)*(sp.eye(4)-singlet)/3
    # For 0<=p<=1/2, this is a convex mixture of twelve product states.
    assert sp.simplify(state-((1-2*p)*aligned+2*p*opposite)) == sp.zeros(4)
    pt = sp.Matrix(4, 4, lambda r,c: state[2*(r//2)+c%2, 2*(c//2)+r%2])
    t = sp.Symbol('t')
    assert sp.expand(pt.charpoly(t).as_expr()-(t-(sp.Rational(1,2)-p))*(t-(1+2*p)/6)**3) == 0


@pytest.mark.parametrize('seed', range(6))
def test_noncommuting_gibbs_stability_bound_is_independent_of_energy_origin(seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    b = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    h = (a+a.conj().T)/3
    perturbation = (b+b.conj().T)/20
    k = h+perturbation
    def thermal(matrix):
        value = expm(-matrix)
        return value/np.trace(value)
    rho, sigma = thermal(h), thermal(k)
    trace_distance = sum(abs(np.linalg.eigvalsh(rho-sigma)))
    energies = np.linalg.eigvalsh(perturbation)
    delta = (energies[-1]-energies[0])/2
    symmetric_entropy = np.trace((rho-sigma)@(logm(rho)-logm(sigma))).real
    assert trace_distance <= delta+2e-14
    assert trace_distance**2 <= symmetric_entropy+2e-14
    assert symmetric_entropy <= delta*trace_distance+2e-14
    shifted, _ = gibbs_state([k, np.eye(3)], [1., 7.])
    np.testing.assert_allclose(shifted, sigma, atol=2e-15, rtol=0)

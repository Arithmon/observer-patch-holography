"""Independent controls for the Hamiltonian supplied to finite A3 numerics."""

import itertools
from fractions import Fraction

import numpy as np
import pytest
from scipy.linalg import expm

from maxent.information_projection import constrained_hamiltonian, gibbs_state


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])
HEISENBERG = sum(np.kron(a, a) for a in (X, Y, Z))


@pytest.mark.parametrize('order', tuple(itertools.permutations(range(3))))
def test_bookkeeping_cannot_erase_an_entangled_thermal_state(order):
    coefficients = [1e20, 1., -1e20]
    multipliers = [coefficients[j] for j in order]
    # This exact binary64 sum is one in every order.
    assert sum(map(Fraction, multipliers)) == 1
    actual_h = constrained_hamiltonian([HEISENBERG]*3, multipliers)
    np.testing.assert_array_equal(actual_h, HEISENBERG)
    rho, log_z = gibbs_state([HEISENBERG]*3, multipliers)
    # H=I-4|singlet><singlet|: diagonalize the known projectors, not H.
    singlet = np.array([0., 1., -1., 0.])/np.sqrt(2)
    projector = np.outer(singlet, singlet)
    probability = 1/(1+3*np.exp(-4))
    expected = probability*projector+(1-probability)*(np.eye(4)-projector)/3
    np.testing.assert_allclose(rho, expected, atol=2e-15, rtol=0)
    assert log_z == pytest.approx(3+np.log1p(3*np.exp(-4)), abs=2e-15)
    partial_transpose = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    assert np.linalg.eigvalsh(partial_transpose)[0] == pytest.approx(.5-probability, abs=2e-15)


@pytest.mark.parametrize('observable', (X, Y, Z, HEISENBERG))
def test_products_are_accumulated_before_their_roundoff_can_cancel(observable):
    # (1+2^-27)(1-2^-27)-1 = -2^-54 exactly. Rounding the product
    # first gives zero even with a compensated sum of rounded products.
    a, b = 1+2.**-27, 1-2.**-27
    expected = -2.**-54*observable
    actual = constrained_hamiltonian([b*observable, observable], [a, -1.])
    np.testing.assert_array_equal(actual, expected)


def test_scalar_partition_energy_survives_product_cancellation():
    a, b = 1+2.**-27, 1-2.**-27
    # Three terms have exact scalar energy -1 after the two large products
    # cancel. This affects log Z even though it cannot change the state.
    rho, log_z = gibbs_state([b*I, I, X], [a*2.**54, -2.**54, .25])
    np.testing.assert_allclose(rho, (I-np.tanh(.25)*X)/2, atol=2e-15, rtol=0)
    assert log_z == pytest.approx(1+np.log(2*np.cosh(.25)), abs=2e-15)


def test_large_cancelled_terms_do_not_hide_a_noncommuting_remainder():
    small = .3*X+.2*Y+.4*Z
    rho, _ = gibbs_state([X, Y, small, X, Y], [1e20, -1e20, 1., -1e20, 1e20])
    exponential = expm(-small)
    np.testing.assert_allclose(rho, exponential/np.trace(exponential), atol=2e-15, rtol=0)

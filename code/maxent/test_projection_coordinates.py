"""Independent controls for observable-coordinate invariance (issue #1033)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from maxent_closure_acceptance import (
    duhamel_covariance, gibbs_state, independent_operator_count, i_projection,
)

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


@pytest.mark.parametrize("scale", [1e-100, 1e-5, 1., 1e5, 1e100, -1e-5])
def test_diagonal_projection_does_not_depend_on_units(scale):
    # Analytic Bernoulli solution, without the Gibbs producer as the oracle.
    target = np.diag([.7, .3])
    assert independent_operator_count([scale*Z]) == 1
    lam, _ = i_projection(target, [scale*Z])
    assert scale*lam[0] == pytest.approx(-.5*np.log(7/3), abs=2e-11)
    rho, _ = gibbs_state([scale*Z], lam)
    np.testing.assert_allclose(rho, target, atol=2e-11, rtol=0)


def test_identity_shift_cannot_destroy_constraint_rank():
    ops = [Z + 1e8*I]
    assert independent_operator_count(ops) == 1
    lam, _ = i_projection(np.diag([.7, .3]), ops)
    assert lam[0] == pytest.approx(-.5*np.log(7/3), abs=2e-11)


def test_noncommuting_rechart_preserves_the_analytic_bloch_solution():
    target = (I + .3*X + .2*Y + .4*Z)/2
    ops = [Z, Z + 1e-5*X]
    assert independent_operator_count(ops) == 2
    lam, _ = i_projection(target, ops)
    rho, _ = gibbs_state(ops, lam)
    # MaxEnt removes only the unconstrained Y component.
    np.testing.assert_allclose(rho, (I + .3*X + .4*Z)/2, atol=2e-10, rtol=0)


def test_scalar_energy_origin_does_not_erase_resolved_offdiagonal_physics():
    # The literal spectrum of 1e20 I + X loses the splitting in float64.
    # That scalar has no effect on the exact normalized Gibbs state.
    rho, _ = gibbs_state([1e20*I + X], [.4])
    np.testing.assert_allclose(rho, (I - np.tanh(.4)*X)/2, atol=2e-14)
    cov = duhamel_covariance([1e20*I + X], [.4])
    np.testing.assert_allclose(cov, [[1/np.cosh(.4)**2]], atol=2e-13)

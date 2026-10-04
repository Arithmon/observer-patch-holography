#!/usr/bin/env python3
"""Regression tests for GitHub issue #539's MaxEnt multiplier-family acceptance test."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from maxent_closure_acceptance import (  # noqa: E402
    cell_densities,
    decimate,
    duhamel_covariance,
    gibbs_state,
    global_sum_constraints,
    i_projection,
    independent_operator_count,
    relative_entropy,
    run_acceptance,
    run_lattice_pair,
    trace_norm,
)


def test_multiplier_count_is_cutoff_independent_for_global_sums() -> None:
    for n_sites in (3, 4, 6, 8):
        constraints = global_sum_constraints(n_sites)
        assert independent_operator_count(constraints) == 2
        assert sum(len(ops) for ops in cell_densities(n_sites)) == 2 * n_sites


def test_i_projection_moment_matches_and_is_unique() -> None:
    fine = global_sum_constraints(6)
    coarse = global_sum_constraints(3)
    omega, _ = gibbs_state(fine, np.array([0.7, 0.4]))
    sigma = decimate(omega, 6)
    lam_star, residual = i_projection(sigma, coarse)
    assert residual < 1e-9
    hessian_floor = np.linalg.eigvalsh(duhamel_covariance(coarse, lam_star)).min()
    assert hessian_floor > 1e-9


def test_pinsker_residual_bound() -> None:
    run = run_lattice_pair(6, np.array([0.7, 0.4]))
    assert run["trace_norm_residual"] <= run["pinsker_residual_bound"] + 1e-9
    assert run["closure_defect_nats"] > 1e-6


def test_transverse_field_subfamily_is_exactly_closed() -> None:
    run = run_lattice_pair(6, np.array([0.0, 0.4]))
    assert run["closure_defect_nats"] < 1e-10
    assert abs(run["induced_map_R_multipliers"][0]) < 1e-8
    assert abs(run["induced_map_R_multipliers"][1] - 0.4) < 1e-8


def test_relative_entropy_is_zero_only_on_the_diagonal() -> None:
    coarse = global_sum_constraints(3)
    rho, _ = gibbs_state(coarse, np.array([0.3, 0.5]))
    tau, _ = gibbs_state(coarse, np.array([0.4, 0.5]))
    assert relative_entropy(rho, rho) < 1e-12
    gap = relative_entropy(rho, tau)
    assert gap > 0
    assert trace_norm(rho - tau) <= (2 * gap) ** 0.5 + 1e-9


def test_full_acceptance_run_passes() -> None:
    results = run_acceptance()
    assert results["all_checks_pass"], results["checks"]


@pytest.mark.parametrize("lam", ([.5], [.5,.2,.1], [np.nan,.2], [True,False]))
def test_gibbs_cannot_silently_drop_or_corrupt_constraints(lam):
    constraints=global_sum_constraints(3)
    for operation in (gibbs_state,duhamel_covariance):
        with pytest.raises(ValueError):
            operation(constraints,np.asarray(lam))


def test_projection_rejects_nonunique_parameters_and_invalid_target():
    identity=np.eye(2)
    with pytest.raises(ValueError,match="independent"):
        i_projection(identity/2,[identity,identity])
    with pytest.raises(ValueError):
        i_projection(np.diag([1.1,-.1]),[np.diag([1.,-1.])])


def test_projection_budget_does_not_return_an_unconverged_minimizer():
    constraints=global_sum_constraints(3)
    sigma,_=gibbs_state(constraints,np.array([.7,.4]))
    with pytest.raises(RuntimeError,match="converge"):
        i_projection(sigma,constraints,max_iter=1)
    with pytest.raises(ValueError):
        i_projection(sigma,constraints,max_iter=0)


def test_gibbs_underflow_is_not_silently_regularized():
    for operation in (gibbs_state,duhamel_covariance):
        with pytest.raises(ValueError,match="underflow"):
            operation([np.diag([0.,1.])],np.array([1000.]))


def test_legacy_relative_entropy_keyword_order_is_preserved():
    state,reference=np.diag([.8,.2]),np.diag([.3,.7])
    expected=.8*np.log(.8/.3)+.2*np.log(.2/.7)
    assert relative_entropy(sigma=state,rho=reference) == pytest.approx(expected)
    assert abs(relative_entropy(sigma=reference,rho=state)-expected) > .01


def test_finite_inputs_cannot_overflow_into_a_nan_state():
    with pytest.raises(ValueError,match="numerical range"):
        gibbs_state([np.diag([1e308,0.])],np.array([1e308]))


def test_gibbs_normalization_cannot_erase_a_subnormal_positive_weight():
    with pytest.raises(ValueError,match="underflow"):
        gibbs_state([np.diag([0.,0.,0.,745.])],np.array([1.]))


def test_gibbs_reconstruction_cannot_erase_a_faithful_direction():
    # Neither exp(-60) nor its normalization underflows. Adding it to the
    # dominant rotated eigenspace loses it when the dense state is formed.
    x=np.array([[0.,1.],[1.,0.]])
    for operation in (gibbs_state,duhamel_covariance):
        with pytest.raises(ValueError,match="faithful"):
            operation([x],np.array([30.]))
    rho,_=gibbs_state([np.diag([0.,60.])],np.array([1.]))
    assert rho[1,1] == pytest.approx(np.exp(-60),rel=1e-14,abs=0)


@pytest.mark.parametrize("bad", (1j,1+0j,np.bool_(True),True,np.nan,np.inf,[1.]))
def test_projection_rejects_nonreal_or_nonscalar_tolerance(bad):
    with pytest.raises(ValueError,match="finite real scalar"):
        i_projection(np.eye(2)/2,[np.diag([1.,-1.])],tol=bad)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))

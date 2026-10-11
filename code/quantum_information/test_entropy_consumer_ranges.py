"""Public entropy receipts must retain the original finite comparison."""

import json

import mpmath
import numpy as np
import pytest

from maxent.information_projection import duhamel_covariance
from maxent.maxent_closure_acceptance import relative_entropy as closure_entropy
from quantum_information.entropy_response import (
    finite_entropy_balance, gibbs_entropy_response,
)


def _coherent_pair(power):
    p, t = 2.**-40, 2.**-power
    rho = np.array([[p, t], [t, 1-p]])
    sigma = np.diag([p, 1-p])
    # Both traces are exactly one, and det(rho)=p(1-p)-t^2>0.
    # Characteristic roots give the entropy without a matrix eigensolver.
    ctx = mpmath.mp.clone()
    ctx.dps = 500
    p, t = ctx.mpf(p), ctx.mpf(t)
    high = (1+ctx.sqrt((1-2*p)**2+4*t*t))/2
    low = (p*(1-p)-t*t)/high
    divergence = ctx.fsum([low*ctx.log(low), high*ctx.log(high),
                           -p*ctx.log(p), -(1-p)*ctx.log1p(-p)])
    return rho, sigma, float(divergence)


@pytest.mark.parametrize("power", [535, 538, 539])
def test_closure_and_serialized_balance_keep_a_positive_coherence_cost(power):
    rho, sigma, expected = _coherent_pair(power)
    assert expected > 0
    assert closure_entropy(sigma=rho, rho=sigma) == expected
    receipt = json.loads(json.dumps(finite_entropy_balance(rho, sigma), allow_nan=False))
    assert receipt["trace_difference"] == receipt["modular_change"] == 0
    assert receipt["relative_entropy"] == receipt["bregman_remainder"] == expected
    assert receipt["entropy_change"] == -expected


def test_unrepresentable_coherence_cost_is_unavailable_in_both_consumers():
    rho, sigma, rounded = _coherent_pair(545)
    assert rounded == 0  # The exact positive cost lies below output range.
    for evaluate in (lambda: closure_entropy(sigma=rho, rho=sigma),
                     lambda: finite_entropy_balance(rho, sigma)):
        with pytest.raises(ValueError, match="range|precision|unresolved"):
            evaluate()


@pytest.mark.parametrize("power", [40, 60, 100])
def test_entropy_response_and_maxent_hessian_describe_the_same_source(power):
    energies = [0., 2.**-power, 1.]
    beta = 700.
    observable = np.diag(energies)
    ctx = mpmath.mp.clone()
    ctx.dps = 400
    exact_energies = [ctx.mpf(e) for e in energies]
    # Independent partition-function derivative; no state producer or
    # covariance implementation supplies this answer.
    def log_partition(b):
        return ctx.log(ctx.fsum(ctx.exp(-b*e) for e in exact_energies))
    expected = float(ctx.diff(log_partition, ctx.mpf(beta), 2))
    response = gibbs_entropy_response(observable, beta)
    hessian = duhamel_covariance([observable], [beta])[0, 0]
    assert response["energy_variance"] == pytest.approx(expected, rel=2e-13, abs=0)
    assert hessian == pytest.approx(expected, rel=2e-13, abs=0)
    assert response["dt_dlambda"] == pytest.approx(-expected, rel=2e-13, abs=0)
    assert response["ds_dlambda"] == pytest.approx(-beta*expected, rel=2e-13, abs=0)

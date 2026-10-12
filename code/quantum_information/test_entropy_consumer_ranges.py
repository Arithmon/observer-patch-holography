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


@pytest.mark.parametrize("phase", [1., 1j])
def test_balance_preserves_original_hermitian_half_subnormal_blocks(phase):
    tiny = np.nextafter(0., 1.)
    reference = np.diag([1., *([tiny]*6)]).astype(complex)
    source = reference.copy()
    for i in (1, 3, 5):
        source[i, i+1] = phase*tiny
    # Exact Hermitian averaging gives three blocks with eigenvalues
    # 3t/2,t/2. Their relative entropy adds before final output rounding.
    ctx = mpmath.mp.clone()
    ctx.dps = 450
    expected = float(3*ctx.mpf(tiny)*(ctx.mpf('1.5')*ctx.log(ctx.mpf('1.5'))
                                     +ctx.mpf('.5')*ctx.log(ctx.mpf('.5'))))
    assert expected == tiny
    report = json.loads(json.dumps(finite_entropy_balance(source, reference), allow_nan=False))
    assert report["relative_entropy"] == report["bregman_remainder"] == expected
    assert report["entropy_change"] == -expected
    assert report["trace_difference"] == report["modular_change"] == 0


def test_balance_keeps_symmetric_subnormal_trace_offsets():
    tiny = np.nextafter(0., 1.)
    source, reference = np.diag([1., 2*tiny]), np.diag([1., 4*tiny])
    report = finite_entropy_balance(source, reference)
    ctx = mpmath.mp.clone()
    ctx.dps = 450
    a, b = 2*ctx.mpf(tiny), 4*ctx.mpf(tiny)
    assert report["relative_entropy"] == float(a*ctx.log(a/b))
    assert report["modular_change"] == float(-(a-b)*ctx.log(b))
    assert report["trace_difference"] == -2*tiny
    assert report["bregman_remainder"] > 0


def test_balance_refuses_a_reference_changed_by_rounded_averaging():
    tiny = np.nextafter(0., 1.)
    reference = np.diag([1., tiny, tiny])
    reference[1, 2] = tiny
    with pytest.raises(ValueError, match="reference Hermitian averaging"):
        finite_entropy_balance(np.diag([1., tiny, tiny]), reference)


def test_balance_refuses_indefinite_source_hidden_by_rounded_averaging():
    source = np.diag([1., 0., 0.])
    source[1, 2] = np.nextafter(0., 1.)
    # The exact Hermitian representative has eigenvalues 1,+t/2,-t/2.
    with pytest.raises(ValueError, match="positive semidefinite"):
        finite_entropy_balance(source, np.eye(3)/3)


def test_balance_retains_ordinary_hermitian_averaging_roundoff():
    reference = np.array([[.4, .1], [np.nextafter(.1, 1.), .6]])
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    coherence = (ctx.mpf(float(reference[0, 1]))+ctx.mpf(float(reference[1, 0])))/2
    determinant = ctx.mpf(.4)*ctx.mpf(.6)-coherence**2
    # For source I/2, Tr(source log(reference))=log(det(reference))/2.
    expected = float(-ctx.log(2)-ctx.log(determinant)/2)
    report = finite_entropy_balance(np.eye(2)/2, reference)
    assert report["relative_entropy"] == pytest.approx(expected, rel=2e-14, abs=0)
    assert report["bregman_remainder"] > 0


def test_balance_pairs_the_original_difference_before_rounding():
    x = .25-2.**-20
    reference = np.array([[.25, x, 0], [x, .5, x], [0, x, .25]], dtype=complex)
    source = reference.copy()
    source[0, 1] += 1j*2.**-30
    source[1, 0] -= 1j*2.**-30
    tiny = np.nextafter(0., 1.)
    source[0, 2] = tiny
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    logarithm = ctx.logm(ctx.matrix(reference.tolist()))
    # The imaginary perturbation pairs to zero with the real logarithm.
    # Exact Hermitian averaging puts t/2 in each of the two corner entries.
    expected = float(-ctx.mpf(tiny)*ctx.re(logarithm[0, 2]))
    assert expected == 4*tiny
    report = finite_entropy_balance(source, reference)
    assert report["modular_change"] == expected
    assert report["relative_entropy"] > 1e-18


@pytest.mark.parametrize("steps", [-11, -3, 3, 11])
@pytest.mark.parametrize("power", [48, 53])
def test_bregman_remainder_precedes_rounding_of_divergence(steps, power):
    value = .5+steps*2.**-power
    source, reference = np.eye(2)*value, np.eye(2)/2
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    x = ctx.mpf(value)
    # Direct scalar Bregman expression, keeping the accepted trace offset.
    expected = float(2*(x*ctx.log(2*x)-x+ctx.mpf('.5')))
    assert expected > 0
    report = finite_entropy_balance(source, reference)
    assert report["trace_difference"] == 2*(value-.5)
    assert report["bregman_remainder"] == pytest.approx(expected, rel=2e-14, abs=0)

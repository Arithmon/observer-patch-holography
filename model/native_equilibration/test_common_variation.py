"""Independent quantum, tangent, region and source-extension controls."""
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import common_variation as v, cut_completion as cut
from ..native_clock import regional_preparation as incident
from ..native_modular_reference import common_caps as common, local_preparation as local
from ..observer import native


def entropy(rho):
    eigenvalues = np.linalg.eigvalsh(rho)
    assert min(eigenvalues) > -1e-12
    positive = eigenvalues[eigenvalues > 1e-15]
    return -float(positive@np.log(positive))


def relative(rho, weights):
    return -entropy(rho)-float(np.diag(rho).real@np.log(weights))


def apply(op, rho, ref=1):
    data = rho.reshape(66, ref, 66, ref).transpose(0, 2, 1, 3).reshape(4356, ref*ref)
    return (op@data).reshape(66, 66, ref, ref).transpose(0, 2, 1, 3).reshape(66*ref, 66*ref)


def state(ref=1):
    rng = np.random.default_rng(9601)
    a = rng.normal(size=(66*ref, 4))+1j*rng.normal(size=(66*ref, 4))
    rho = a@a.conj().T
    return .7*rho/np.trace(rho)+.3*np.eye(66*ref)/(66*ref)


def independent_center(rho, ref=1):
    # Conditional reference states must remain attached to their charge.
    data = rho.reshape(66, ref, 66, ref)
    out = np.zeros_like(data)
    for q in range(3):
        indices = [i for i, x in enumerate(native.states(binary=True)) if x[0]+x[1] == q]
        conditional = sum(data[i, :, i, :] for i in indices)/len(indices)
        for i in indices:
            out[i, :, i, :] = conditional
    return out.reshape(rho.shape)


def test_retained_exact_source_and_all_region_receipt():
    report = v.run()
    assert json.loads(v.OUTPUT.read_text()) == report
    for law in ("C", "E"):
        assert report["native_source"]["laws"][law]["matrix_units"] == 4356
    assert report["arbitrary_input_preparation"]["native_work_attempts"] == 18601
    for row in report["selected_references"].values():
        assert row["region_envelope"]["regions_checked"] == 1024


@pytest.mark.parametrize("law", ["C", "E"])
def test_actual_quantum_source_converges_to_two_moment_maxent(law):
    rho = state(ref=2)
    z, zs = v.center_counts()
    target = independent_center(rho, ref=2)
    np.testing.assert_allclose(apply(z/zs, rho, 2), target, atol=1e-16)
    p, ps = common.mixed_counts(law)
    core = incident.seam_counts(v.CORE, law)/2
    current = rho.copy()
    for _ in range(420):
        current = apply(p/ps, current, 2)
    current = apply(core, current, 2)
    assert np.linalg.norm(current-target) < 2e-7
    e, es = cut.completion_counts(v.CORE)
    np.testing.assert_allclose(apply(core, apply(e/es, rho, 2), 2), target, atol=2e-16)
    # The reference remains correlated with charge; no product convergence.
    source = np.trace(target.reshape(66, 2, 66, 2), axis1=1, axis2=3)
    reference = np.trace(target.reshape(66, 2, 66, 2), axis1=0, axis2=2)
    assert np.linalg.norm(target-np.kron(source, reference)) > .003
    assert entropy(target)-entropy(rho) > .1


@pytest.mark.parametrize("reference", [common.reference, local.reference])
def test_coherent_first_law_and_strict_full_matrix_stationarity(reference):
    w = np.array(list(map(float, reference())))
    sigma = np.diag(w).astype(complex)
    rho = state()
    q = np.array(v.charges())
    pair = q*(q-1)/2
    cert = v.family_certificate(*v.moments(reference()))
    h1 = np.log(float(cert["linear_multiplier_log_argument"]))
    h2 = np.log(float(cert["pair_multiplier_log_argument"]))
    linear = h1*(q@(np.diag(rho).real-w))+h2*(pair@(np.diag(rho).real-w))
    assert abs(entropy(rho)-entropy(sigma)-linear+relative(rho, w)) < 3e-14
    x = rho-independent_center(rho)
    x /= np.linalg.norm(x)
    assert abs(np.trace(x)) < 1e-15
    assert abs(q@np.diag(x).real)+abs(pair@np.diag(x).real) < 1e-15
    # The entropy gradient on normalized tangents is defined modulo identity.
    # Remove that constant before contraction: otherwise the tiny floating
    # trace residual above is amplified by log(w[0]), producing a false
    # stationarity failure on a different BLAS implementation.
    assert abs(np.diag(x).real@np.log(w/w[0])) < 1e-15
    divided = np.empty((66, 66))
    for i, a in enumerate(w):
        for j, b in enumerate(w):
            divided[i, j] = 1/a if a == b else (np.log(a)-np.log(b))/(a-b)
    curvature = -float(np.sum(abs(x)**2*divided))
    assert -1/min(w) <= curvature <= -1/max(w)
    step = 2e-6
    difference = (entropy(sigma+step*x)+entropy(sigma-step*x)-2*entropy(sigma))/step**2
    assert abs(difference-curvature) < .003
    for signed in (-1, 1):
        varied = sigma+signed*min(w)/3*x
        assert abs(entropy(varied)-entropy(sigma)+relative(varied, w)) < 2e-14
        assert entropy(varied) < entropy(sigma)


@pytest.mark.parametrize("reference", [common.reference, local.reference])
def test_stationarity_readout_removes_only_the_identity_component(reference):
    w = np.array(list(map(float, reference())))
    q = np.array(v.charges())
    within = np.flatnonzero(q == 0)
    other = np.flatnonzero(q == 1)[0]
    diagonal = np.zeros(66)
    diagonal[within[0]] = .5
    diagonal[within[1]] = -.5+2*np.finfo(float).eps
    # This fixture satisfies the existing trace/moment precision gates, yet
    # its irrelevant constant contribution exceeds the old slope threshold.
    assert 0 < abs(diagonal.sum()) < 1e-15
    assert abs(diagonal@np.log(w)) > 1e-15
    assert diagonal@np.log(w/w[0]) == 0
    # A genuine moment-changing direction still has a nonzero slope.
    diagonal[:] = 0
    diagonal[within[0]], diagonal[other] = .5, -.5
    assert abs(diagonal@np.log(w/w[0])) > 1e-3


@pytest.mark.parametrize("reference", [common.reference, local.reference])
def test_explicit_envelope_derivatives_and_two_parameter_strict_hessian(reference):
    mean, pair = v.moments(reference())
    cert = v.family_certificate(mean, pair)
    h = np.array([[float(x) for x in row] for row in cert["entropy_hessian"]])
    assert np.linalg.eigvalsh(h).max() < 0
    assert abs(np.linalg.det(h)-float(cert["entropy_hessian_determinant"])) < 2e-10
    delta = F(1, 10**6)
    for dm, dc, key in ((1, 0, "linear_multiplier_log_argument"), (0, 1, "pair_multiplier_log_argument")):
        plus = np.diag(list(map(float, v.family(mean+dm*delta, pair+dc*delta))))
        minus = np.diag(list(map(float, v.family(mean-dm*delta, pair-dc*delta))))
        derivative = (entropy(plus)-entropy(minus))/(2*float(delta))
        assert abs(derivative-np.log(float(cert[key]))) < 2e-9


def test_charge_mean_alone_does_not_select_the_common_entropy_maximum():
    reference = common.reference()
    mean, pair = v.moments(reference)
    assert mean == F(1, 3) == v.moments((F(1, 66),)*66)[0]
    assert pair == F(1, 75) != F(1, 66)
    cert = v.family_certificate(mean, pair)
    assert cert["pair_multiplier_log_argument"] > 1
    moved = np.diag(list(map(float, v.family(mean, pair+F(1, 10000)))))
    assert entropy(moved) > entropy(np.diag(list(map(float, reference))))
    assert entropy(np.eye(66)/66) > entropy(moved)


@pytest.mark.parametrize("law", ["C", "E"])
def test_core_attempt_removes_a_real_fixed_coherence_and_keeps_both_moments(law):
    sigma = np.diag(list(map(float, common.reference()))).astype(complex)
    # Exact bare-core |10><01| tensored with identity on its ten spectators.
    basis = native.states(binary=True)
    x = np.zeros((66, 66), complex)
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            if a[2:] == b[2:] and a[:2] == (1, 0) and b[:2] == (0, 1):
                x[i, j] = x[j, i] = .001
    rho = sigma+x
    p, ps = common.mixed_counts(law)
    np.testing.assert_allclose(apply(p/ps, rho), rho, atol=1e-17)
    np.testing.assert_allclose(apply(incident.seam_counts(v.CORE, law)/2, rho), sigma, atol=1e-17)
    assert entropy(rho) < entropy(sigma)


def test_region_changes_keep_one_envelope_but_unprotected_changes_do_not():
    reference = local.reference()
    small = v.region_entropy(v.CORE, reference)
    big = v.region_entropy(common.CAPS[0], reference)
    full = v.region_entropy(tuple(range(12)), reference)
    assert small["counted_edge"] != big["counted_edge"]
    assert full["counted_edge"] == {}
    assert small["completion"] == big["completion"] == full["completion"]
    outside = v.region_entropy((2, 3, 4, 6), reference)
    assert outside["completion"] != outside["source"]
    sigma = np.diag(list(map(float, reference)))
    e, es = cut.completion_counts((2, 3, 4, 6))
    out = apply(e/es, sigma)
    assert entropy(out) > entropy(sigma)+.01
    assert abs(np.diag(out)@v.charges()-np.diag(sigma)@v.charges()) > .01


@pytest.mark.parametrize("bits", [1, 24, 1000])
def test_exact_entropy_budget_and_both_finite_selected_reference_tangents(bits):
    b = v.preparation_budget(bits)
    epsilon = b["reference_uniform_halftrace_upper"]
    assert epsilon < F(1, 2**bits)
    assert b["native_work_attempts"] == b["mixed_complement_attempts"]+1
    for kind in ("global", "local"):
        cert = v.stationarity_budget(kind)
        assert 0 < cert["entropy_tangent_error_per_trace_norm_upper"] < F(1, 5000)
        assert cert["simultaneously_admitted_caps"] == 1024


@pytest.mark.parametrize("args", [(True, F(1, 4)), (.5, .1), (F(1, 3), 0), (1, F(1, 2)), (2, F(1, 4)), (0, 0)])
def test_nonfaithful_or_inexact_family_rejected(args):
    with pytest.raises(ValueError):
        v.family(*args)


@pytest.mark.parametrize("cap", [[0, 1], (1, 0), (True, 1), (0, 0), (0, 12)])
def test_malformed_region_rejected(cap):
    with pytest.raises(ValueError):
        v.region_entropy(cap, local.reference())


def test_nonfamily_populations_do_not_get_a_region_stationarity_certificate():
    weights = list(local.reference())
    weights[0] += F(1, 100000); weights[1] -= F(1, 100000)
    with pytest.raises(ValueError):
        v.region_entropy(v.CORE, tuple(weights))
    with pytest.raises(ValueError):
        v.stationarity_budget("fitted")

"""All-column checks against independent oscillator and matrix-function oracles."""
from math import pi, sqrt

import numpy as np
import pytest
from scipy.linalg import cosm, expm, sinm

from . import source
from .lattice import GalerkinParameters


def oscillator(d, omega):
    extended = d+6
    lowering = np.zeros((extended, extended))
    for n in range(1, extended):
        lowering[n-1, n] = sqrt(n)
    q = (lowering+lowering.T)/sqrt(2*omega)
    p = 1j*sqrt(omega/2)*(lowering.T-lowering)
    return q[:d, :d], (q@q)[:d, :d], (p@p)[:d, :d], (q@q@q@q)[:d, :d]


def check_execution(metadata):
    resources = metadata["resources"]
    assert resources["native_pulses"] > 0
    assert resources["proper_code_load_events"] == resources["native_pulses"]
    assert resources["proper_code_export_events"] == resources["native_pulses"]
    assert resources["proper_code_service_events"] == 2*resources["native_pulses"]
    for key in source.COUNT_KEYS:
        assert resources[key] == sum(stage["resources"][key]
                                     for stage in metadata["source_stages"])
    assert metadata["full_phase_target_residual"] < 2e-10
    assert metadata["service_failure_outcomes_retained"]
    assert not metadata["physical_service_error_certified"]


@pytest.mark.parametrize("d,mode,grid", [(2, 0, None), (4, 0, 5), (4, 3, 9)])
def test_onsite_full_phase_matches_independent_projected_oscillator(d, mode, grid):
    p = GalerkinParameters(grid=grid, mass=.8, coupling=.6, omega=1.4)
    q, q2, p2, q4 = oscillator(d, p.omega)
    spatial = 1 if grid is None else (grid*np.sin(pi/grid)/pi)**2
    frequency = p.mass**2+(spatial if mode else 0)
    quartic = p.coupling*(1.5 if mode else 1)/24
    h = p2/2+frequency*q2/2+quartic*q4
    gate, metadata = source.onsite(d, mode, .13, p)
    np.testing.assert_allclose(gate, expm(-.13j*h), atol=7e-12, rtol=0)
    np.testing.assert_allclose(gate.conj().T@gate, np.eye(d), atol=8e-13, rtol=0)
    assert np.linalg.norm(gate-expm(-.13j*(h-h[0, 0]*np.eye(d))), 2) > .01
    wrong = p2/2+frequency*q2/2+quartic*np.linalg.matrix_power(q, 4)
    assert np.linalg.norm(gate-expm(-.13j*wrong), 2) > 1e-3
    check_execution(metadata)


@pytest.mark.parametrize("dimensions", [(2, 4), (4, 2), (4, 4)])
def test_pair_gate_uses_projected_second_moments_and_unequal_tensor_order(dimensions):
    left, right = dimensions
    omega, coefficient, duration = 1.6, .23, .09
    qleft, q2left, _, _ = oscillator(left, omega)
    qright, q2right, _, _ = oscillator(right, omega)
    expected = expm(-1j*duration*coefficient*np.kron(q2right, q2left))
    gate, metadata = source.pair(left, right, coefficient, duration, omega=omega)
    np.testing.assert_allclose(gate, expected, atol=7e-12, rtol=0)
    wrong = expm(-1j*duration*coefficient*np.kron(qright@qright, qleft@qleft))
    assert np.linalg.norm(gate-wrong, 2) > .01
    inverse, _ = source.pair(left, right, coefficient, -duration, omega=omega)
    np.testing.assert_allclose(inverse@gate, np.eye(left*right), atol=9e-12, rtol=0)
    check_execution(metadata)


@pytest.mark.parametrize("d,eta,theta", [(2, 1., pi/2), (4, .7, 0.), (4, -.8, .23)])
def test_smooth_instrument_kraus_completeness_and_reference_backaction(d, eta, theta):
    omega = 1.3
    q = oscillator(d, omega)[0]
    angle = (eta*q+theta*np.eye(d))/2
    expected = (cosm(angle), sinm(angle))
    kraus, metadata = source.instrument(d, eta=eta, theta=theta, omega=omega)
    for actual, target in zip(kraus, expected, strict=True):
        np.testing.assert_allclose(actual, target, atol=7e-12, rtol=0)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(d),
                               atol=8e-13, rtol=0)
    np.testing.assert_allclose(kraus[0].conj().T@kraus[0],
                               (np.eye(d)+cosm(eta*q+theta*np.eye(d)))/2,
                               atol=8e-13, rtol=0)
    pointer = np.block([[kraus[0], -kraus[1]], [kraus[1], kraus[0]]])
    np.testing.assert_allclose(pointer.conj().T@pointer, np.eye(2*d),
                               atol=1e-12, rtol=0)
    rng = np.random.default_rng(701)
    reference = rng.normal(size=(d, 3))+1j*rng.normal(size=(d, 3))
    reference *= .27/np.linalg.norm(reference)
    for actual, target in zip(kraus, expected, strict=True):
        np.testing.assert_allclose(actual@reference, target@reference, atol=2e-12)
    assert sum(np.linalg.norm(k@reference)**2 for k in kraus) == pytest.approx(.27**2)
    check_execution(metadata)


def test_zero_exposure_is_state_independent_and_feedback_has_declared_sign():
    kraus, _ = source.instrument(4, eta=0., theta=.7)
    np.testing.assert_allclose(kraus[0], np.cos(.35)*np.eye(4), atol=5e-12)
    np.testing.assert_allclose(kraus[1], np.sin(.35)*np.eye(4), atol=5e-12)
    q = oscillator(4, 1.7)[0]
    feedback, metadata = source.feedback(4, kick=.21, omega=1.7)
    np.testing.assert_allclose(feedback, expm(.21j*q), atol=6e-12, rtol=0)
    assert np.linalg.norm(feedback-expm(-.21j*q), 2) > .1
    control, _ = source.feedback(4, kick=0., omega=1.7)
    np.testing.assert_allclose(control, np.eye(4), atol=6e-12, rtol=0)
    check_execution(metadata)


def test_finite_detector_is_not_projection_of_unrestricted_detector():
    dimension, omega, eta, theta = 4, 2., 1., pi/2
    kraus, metadata = source.instrument(dimension, eta=eta, theta=theta, omega=omega)
    # A larger independent oscillator exposes leakage across the retained
    # occupation boundary.  Merely projecting its Kraus matrices does not
    # reproduce the declared complete finite-code instrument.
    large_q = oscillator(24, omega)[0]
    large_angle = (eta*large_q+theta*np.eye(24))/2
    projected = [function(large_angle)[:dimension, :dimension]
                 for function in (cosm, sinm)]
    assert max(np.linalg.norm(actual-target, 2)
               for actual, target in zip(kraus, projected, strict=True)) > .01
    missing = np.eye(dimension)-sum(k.conj().T@k for k in projected)
    assert np.linalg.norm(missing, 2) > .01
    assert metadata["unrestricted_instrument_compression_claimed"] is False


def test_source_execution_fault_cannot_be_hidden_by_returning_target(monkeypatch):
    original = source.c.execute

    def erroneous_phase(*args, **kwargs):
        return np.exp(.001j)*original(*args, **kwargs)

    monkeypatch.setattr(source.c, "execute", erroneous_phase)
    # Unique duration prevents a prior cached circuit from bypassing execution.
    with pytest.raises(ArithmeticError, match="full.phase"):
        source.onsite(2, 0, .01713971, GalerkinParameters())


def test_returned_source_buffers_and_metadata_do_not_mutate_cached_execution():
    matrix, metadata = source.feedback(2, kick=.173)
    original = matrix.copy()
    matrix[:] = 0
    metadata["resources"]["native_pulses"] = -1
    again, ledger = source.feedback(2, kick=.173)
    np.testing.assert_allclose(again, original)
    assert ledger["resources"]["native_pulses"] > 0


@pytest.mark.parametrize("call", [
    lambda: source.onsite(3, 0, .1),
    lambda: source.onsite(2, True, .1),
    lambda: source.onsite(2, 7, .1),
    lambda: source.onsite(2, 0, float("nan")),
    lambda: source.pair(2, 4, .2, float("inf")),
    lambda: source.instrument(2, eta=float("nan")),
    lambda: source.instrument(2, theta=True),
    lambda: source.feedback(4, kick=.1, omega=0),
])
def test_invalid_source_inputs_rejected(call):
    with pytest.raises(ValueError):
        call()

"""Independent finite-operator, Lindblad and selection checks on executed pulses."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.anharmonic_observations import source as s
from model.scalar_observations.experimental_readout import executable_source as native
from model.observer.source_scalar import compiler as c

DELTA = 2*math.pi*500
OMEGA = math.pi/150e-6
TAU = 50e-6
SETTINGS = [{'u': 0., 'v': .4, 'theta': 0.},
            {'u': 0., 'v': .4, 'theta': -math.pi/2}]


def independent_hamiltonian(dimension, sign, omega=OMEGA):
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    return DELTA*np.kron(np.eye(2), np.diag(np.arange(dimension))) - sign*omega/2*np.kron(np.array([[0., 1.], [1., 0.]]), a+a.T)


def independent_initial(dimension, alpha=.7, thermal=.1):
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    displacement = expm(alpha*(a.T-a))
    ratio = thermal/(1+thermal)
    weights = (1-ratio)*ratio**np.arange(dimension)
    weights[-1] += ratio**dimension
    rho = np.zeros((2*dimension, 2*dimension), complex)
    rho[:dimension, :dimension] = (displacement*weights) @ displacement.conj().T
    return rho


def independent_evolution(dimension, steps, gamma=0., reset=False):
    size = 2*dimension
    number = np.kron(np.eye(2), np.diag(np.arange(dimension)))
    identity = np.eye(size)
    propagators = {}
    for sign, duration in [(-1, TAU), (1, 2*TAU)]:
        h = independent_hamiltonian(dimension, sign)
        generator = -1j*(np.kron(identity, h)-np.kron(h.T, identity))
        generator += gamma*(np.kron(number.conj(), number)-.5*np.kron(identity, number@number)-.5*np.kron((number@number).T, identity))
        propagators[sign] = expm(duration*generator)
    spin = np.kron(np.diag(np.exp(-.4j*np.array([1., -1.]))), np.eye(dimension))
    rho = independent_initial(dimension)
    for _ in range(steps):
        for sign, rotation in [(-1, True), (1, True), (-1, False)]:
            rho = (propagators[sign] @ rho.reshape(-1, order='F')).reshape((size, size), order='F')
            if rotation:
                rho = spin @ rho @ spin.conj().T
        if reset:
            rho[dimension:, :] = 0
            rho[:, dimension:] = 0
            rho /= np.trace(rho)
    return rho


@pytest.mark.parametrize('dimension', [2, 4, 8])
@pytest.mark.parametrize('sign', [-1, 1])
def test_source_pulse_matches_independent_finite_continuous_detuning(dimension, sign):
    actual, metadata = s.pulse_unitary(dimension, DELTA, OMEGA, TAU/3, sign)
    target = expm(-1j*TAU/3*independent_hamiltonian(dimension, sign))
    assert np.linalg.norm(actual-target, 2) < 3e-12
    assert np.linalg.norm(actual.conj().T @ actual-np.eye(2*dimension), 2) < 3e-12
    assert metadata['resources']['native_pulses'] > 0
    assert not metadata['cutoff_CCR_identity_used']
    assert all(stage['matrix_origin'] == 'contracted_native_source_pulses' for stage in metadata['stages'])


def test_complete_spin_evolution_and_conditional_readout_match_dense_reference():
    dimension = 4
    result = s.run(dimension, .7, .1, [0, 1, 3], SETTINGS, dephasing_rate=0., subdivisions=2)
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    for record in result['records']:
        rho = independent_evolution(dimension, record['steps'])
        accepted = rho[:dimension, :dimension]
        acceptance = np.trace(accepted).real
        assert record['spin_probabilities'] == pytest.approx([acceptance, 1-acceptance], abs=2e-12)
        characteristic = np.trace(accepted @ expm(.4j*(a.T+a)))/acceptance
        assert [read['selected_mean'] for read in record['reads']] == pytest.approx([characteristic.real, characteristic.imag], abs=2e-12)
        for read in record['reads']:
            joint = np.asarray(read['spin_detector_joint_probabilities'])
            assert joint.sum() == pytest.approx(1., abs=3e-12)
            assert joint.sum(axis=1) == pytest.approx(record['spin_probabilities'], abs=3e-12)
            assert min(joint.ravel()) > -2e-13
    assert result['both_final_spin_branches_retained']
    assert not result['intermediate_spin_projection_or_reset']


def test_intermediate_spin_projection_is_measurably_different():
    rho = independent_evolution(4, 6)
    reset = independent_evolution(4, 6, reset=True)
    conditional = rho[:4, :4]/np.trace(rho[:4, :4])
    assert np.linalg.norm(conditional-reset[:4, :4]) > .02
    result = s.run(4, .7, .1, [6], SETTINGS, dephasing_rate=0.)
    assert result['records'][0]['spin_probabilities'][1] > .02


def test_supplied_diffusion_preserves_spin_and_has_declared_number_rate():
    rho = np.ones((8, 8), complex)/8
    actual = s.phase_diffusion(rho, 4, 18., .002)
    assert actual[0, 4] == rho[0, 4]
    assert actual[0, 3] == pytest.approx(rho[0, 3]*math.exp(-18*.002*9/2))
    assert np.trace(actual) == pytest.approx(1.)
    assert np.min(np.linalg.eigvalsh(actual)) > -2e-15


def test_source_noise_splitting_converges_to_independent_lindblad_channel():
    rho = independent_evolution(4, 4, gamma=18.)
    expected = np.trace(rho[:4, :4]).real
    errors = []
    for subdivisions in (1, 2, 4):
        result = s.run(4, .7, .1, [4], SETTINGS, subdivisions=subdivisions)
        record = result['records'][0]
        errors.append(abs(record['spin_probabilities'][0]-expected))
        assert record['total_probability_mass'] == pytest.approx(1., abs=3e-12)
    assert errors[1] < errors[0]/3
    assert errors[2] < errors[1]/3
    assert errors[-1] < 4e-5


def test_zero_force_is_executed_harmonic_control():
    result = s.run(4, .7, .1, [3], SETTINGS, omega=0., dephasing_rate=0.)
    rho = independent_initial(4)[:4, :4]
    rotation = np.diag(np.exp(-1j*DELTA*(12*TAU)*np.arange(4)))
    rho = rotation @ rho @ rotation.conj().T
    a = np.diag(np.sqrt(np.arange(1, 4)), 1)
    characteristic = np.trace(rho @ expm(.4j*(a.T+a)))
    record = result['records'][0]
    assert record['spin_probabilities'] == pytest.approx([1., 0.], abs=2e-12)
    assert [read['selected_mean'] for read in record['reads']] == pytest.approx([characteristic.real, characteristic.imag], abs=3e-12)


def test_record_snapshots_do_not_change_continuing_state():
    sparse = s.run(4, .7, .1, [4], SETTINGS, subdivisions=2)
    dense = s.run(4, .7, .1, [0, 1, 2, 3, 4], SETTINGS, subdivisions=2)
    assert sparse['records'][-1] == dense['records'][-1]


def test_every_subpulse_application_is_in_resource_ledger():
    result = s.run(2, .7, .1, [1], SETTINGS, subdivisions=3)
    operations = result['operations']
    for key in native.COUNT_KEYS:
        expected = 6*sum(operations[name]['resources'][key] for name in ('positive_subpulse', 'negative_subpulse'))
        expected += 2*sum(operations[name]['resources'][key] for name in ('spin_z', 'spin_x', 'spin_x_inverse'))
        assert result['resources_per_evolution_step'][key] == expected


def test_failed_native_execution_cannot_be_replaced_with_hamiltonian_target(monkeypatch):
    s.spin_rotation.cache_clear()
    monkeypatch.setattr(c, 'execute', lambda word, qubits: np.zeros((2**qubits, 2**qubits), complex))
    with pytest.raises(ArithmeticError):
        s.spin_rotation('Z', .81234567)
    s.spin_rotation.cache_clear()


@pytest.mark.parametrize('steps', [[], [1, 1], [2, 1], [-1], [True], [1.2], None])
def test_invalid_step_schedules_rejected(steps):
    with pytest.raises(ValueError):
        s.run(2, .7, .1, steps, SETTINGS)


@pytest.mark.parametrize('settings', [[], [{'u': 0., 'v': .4}], [{'u': 0., 'v': .4, 'theta': 0., 'observed': .2}], [{'u': True, 'v': .4, 'theta': 0.}]])
def test_observed_values_and_invalid_settings_rejected(settings):
    with pytest.raises(ValueError):
        s.run(2, .7, .1, [0], settings)


@pytest.mark.parametrize('kwargs', [{'dephasing_rate': -1.}, {'subdivisions': 0}, {'subdivisions': True}, {'pulse_duration': 0.}, {'theta': float('nan')}, {'phi': True}, {'omega': -1.}, {'delta': 0.}])
def test_invalid_apparatus_rejected(kwargs):
    with pytest.raises(ValueError):
        s.run(2, .7, .1, [0], SETTINGS, **kwargs)


@pytest.mark.parametrize('alpha0', [.1, .3])
def test_small_angle_step_has_positive_x_cosine_potential_before_any_data(alpha0):
    """The stated target fixes the phase/sign convention independently of a fit."""
    from scipy.linalg import cosm
    dimension = 4
    # A tiny positive delta approaches the frozen-force identity while retaining
    # the same exact finite Hamiltonian source path; subtract the theta=0 drift.
    positive, _ = s.pulse_unitary(dimension, 1e-5, 2*alpha0/TAU, TAU, 1)
    negative, _ = s.pulse_unitary(dimension, 1e-5, 2*alpha0/TAU, TAU, -1)
    angle = 1e-6
    spin, _ = s.spin_rotation('Z', angle)
    spin = np.kron(spin, np.eye(dimension))
    step = negative @ spin @ positive @ positive @ spin @ negative
    baseline = negative @ positive @ positive @ negative
    derivative = (step[:dimension, :dimension]-baseline[:dimension, :dimension])/angle
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    expected = -1j*cosm(2*alpha0*(ladder+ladder.T))
    assert np.linalg.norm(derivative-expected, 2) < 2e-6

"""Independent pulse, gap, readout and complete-history execution checks."""
import math

import numpy as np
import pytest
from scipy.linalg import cosm, expm

from model.anharmonic_observations import reference
from model.anharmonic_observations import source as old
from model.anharmonic_transfer import source as s
from model.scalar_observations.experimental_readout import executable_source as native


DELTA = 2*math.pi*500
OMEGA = math.pi/150e-6
CARRIER = math.pi/35e-6
TAU = 50e-6
SETTINGS = [{'u': 0., 'v': .4, 'theta': -math.pi/2},
            {'u': .3, 'v': .2, 'theta': 0.}]


def independent_step(rho, *, phi=.2, gamma=18., finite=True):
    """Unsplit DOP853 sidebands and an independently constructed gap map."""
    d = len(rho)//2
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    number = np.kron(np.eye(2), np.diag(np.arange(d)))
    zrotation = expm(-.4j*np.kron(z, np.eye(d)))
    levels = np.tile(np.arange(d), 2)

    def gap(state, angle):
        duration = abs(angle)/CARRIER if finite else 0.
        h_times_duration = DELTA*duration*number+angle*np.kron(x, np.eye(d))/2
        u = expm(-1j*h_times_duration)
        output = u @ state @ u.conj().T
        return output*np.exp(-gamma*duration*(levels[:, None]-levels[None, :])**2/2)

    for sign, duration, rotations in [(-1, TAU, (-phi, phi)),
                                      (1, 2*TAU, (phi, -phi)),
                                      (-1, TAU, None)]:
        rho = reference.lindblad_pulse(rho, delta=DELTA, omega=OMEGA,
                                       duration=duration, sign=sign,
                                       dephasing_rate=gamma)
        if rotations:
            rho = gap(rho, rotations[0])
            rho = zrotation @ rho @ zrotation.conj().T
            rho = gap(rho, rotations[1])
    return rho


@pytest.mark.parametrize('dimension', [2, 4, 8])
@pytest.mark.parametrize('angle', [-.3, 0., .3])
def test_finite_gap_is_full_spin_and_number_exponential(dimension, angle):
    actual, meta = s.spin_gap(dimension, angle, DELTA, CARRIER)
    duration = abs(angle)/CARRIER
    hdt = (DELTA*duration*np.kron(np.eye(2), np.diag(np.arange(dimension)))
           + angle*np.kron(np.array([[0., 1.], [1., 0.]]), np.eye(dimension))/2)
    assert np.linalg.norm(actual-expm(-1j*hdt), 2) < 2e-12
    assert meta['duration_seconds'] == duration
    assert meta['matrix_origin'] == 'contracted_native_source_pulses'
    assert all(stage['word_sha256'] for stage in meta['stages'])


@pytest.mark.parametrize('finite', [False, True])
def test_complete_native_history_matches_independent_unsplit_master_equation(finite):
    result = s.run(4, .7, .1, [0, 1, 3], SETTINGS, phi=.2,
                   subdivisions=16, finite_spin_duration=finite)
    rho = reference.preparation(4, .7, .1)
    wanted = {record['steps']: record for record in result['records']}
    for step in range(4):
        if step in wanted:
            expected = reference.readout(rho, SETTINGS)
            actual = wanted[step]
            budget = actual['execution_error_budget']
            assert actual['spin_probabilities'] == pytest.approx(
                expected['spin_probabilities'], abs=budget['complete_trace_distance_upper']+4e-12)
            for got, target in zip(actual['reads'], expected['reads'], strict=True):
                p, q = (np.asarray(row['spin_detector_joint_probabilities']) for row in (got, target))
                assert np.sum(np.abs(p-q))/2 <= budget['complete_trace_distance_upper']+4e-12
                assert abs(got['selected_mean']-target['selected_mean']) <= budget['selected_readout']['conditional_spin_mean_error_upper']+6e-12
                assert np.sum(p) == pytest.approx(1., abs=4e-12)
                assert np.sum(p, axis=1) == pytest.approx(actual['spin_probabilities'], abs=4e-12)
            assert actual['drive_elapsed_duration_seconds'] == pytest.approx(step*(4*TAU+4*abs(.2)/CARRIER*finite))
            assert actual['simulated_duration_seconds'] == pytest.approx(step*4*TAU)
        if step < 3:
            rho = independent_step(rho, finite=finite)


def test_literal_beta_readout_and_two_point_estimate_are_not_exact_position():
    result = s.run(8, -.8, .1, [0], subdivisions=1)
    actual = result['records'][0]
    rho = reference.preparation(8, -.8, .1)
    expected = reference.readout(rho, [s.position_setting()])
    assert actual['reads'][0]['selected_mean'] == pytest.approx(expected['reads'][0]['selected_mean'], abs=3e-12)
    assert actual['reads'][0]['v'] == .4
    estimated = actual['reads'][0]['selected_mean']/(math.sqrt(2)*.4)
    assert abs(estimated-actual['selected_position_mean']) > .1


def test_zero_phase_retains_frozen_symmetric_predictions():
    kwargs = dict(phi=0., subdivisions=4)
    before = old.run(4, .7, .1, [0, 1, 3], SETTINGS, **kwargs)
    after = s.run(4, .7, .1, [0, 1, 3], SETTINGS, **kwargs)
    for a, b in zip(before['records'], after['records'], strict=True):
        assert a['spin_probabilities'] == pytest.approx(b['spin_probabilities'], abs=3e-12)
        assert [x['selected_mean'] for x in a['reads']] == pytest.approx([x['selected_mean'] for x in b['reads']], abs=3e-12)
        assert b['drive_elapsed_duration_seconds'] == b['simulated_duration_seconds']


def test_finite_spin_duration_changes_dynamics_and_has_exact_declared_clock():
    primary = s.run(4, .7, .1, [5], SETTINGS, phi=math.pi/10, subdivisions=4)
    instantaneous = s.run(4, .7, .1, [5], SETTINGS, phi=math.pi/10,
                          subdivisions=4, finite_spin_duration=False)
    assert primary['settings']['spin_gap_duration_seconds'] == pytest.approx(3.5e-6)
    a, b = primary['records'][0], instantaneous['records'][0]
    assert a['drive_elapsed_duration_seconds'] == pytest.approx(5*214e-6)
    assert b['drive_elapsed_duration_seconds'] == pytest.approx(5*200e-6)
    assert abs(a['reads'][0]['selected_mean']-b['reads'][0]['selected_mean']) > .005


def test_snapshots_do_not_reset_or_project_the_continuing_spin():
    kwargs = dict(phi=.2, subdivisions=4)
    sparse = s.run(4, .7, .1, [5], SETTINGS, **kwargs)
    dense = s.run(4, .7, .1, [0, 1, 2, 3, 4, 5], SETTINGS, **kwargs)
    assert sparse['records'][-1] == dense['records'][-1]
    assert sparse['records'][-1]['spin_probabilities'][1] > .005
    reset = reference.preparation(4, .7, .1)
    for _ in range(5):
        reset = independent_step(reset)
        reset[4:, :] = reset[:, 4:] = 0
        reset /= np.trace(reset)
    wrong = reference.readout(reset, SETTINGS)
    assert abs(wrong['reads'][0]['selected_mean']-sparse['records'][0]['reads'][0]['selected_mean']) > .001
    assert not sparse['intermediate_spin_projection_or_reset']


def test_phase_convention_is_fixed_by_native_small_angle_generator():
    dimension, theta, phi, alpha0 = 4, 1e-6, .2, .3
    delta = 1e-12
    omega = 2*alpha0/TAU
    negative, _ = old.pulse_unitary(dimension, delta, omega, TAU, -1)
    positive, _ = old.pulse_unitary(dimension, delta, omega, 2*TAU, 1)
    plus, _ = s.spin_gap(dimension, phi, delta, CARRIER, False)
    minus, _ = s.spin_gap(dimension, -phi, delta, CARRIER, False)
    z, _ = old.spin_rotation('Z', theta)
    z = np.kron(z, np.eye(dimension))
    unitary = negative @ minus @ z @ plus @ positive @ plus @ z @ minus @ negative
    generator = 1j*(unitary[:dimension, :dimension]-np.eye(dimension))/theta
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    argument = 2*alpha0*(a+a.T)
    expected = cosm(argument-phi*np.eye(dimension))
    reversed_phase = cosm(argument+phi*np.eye(dimension))
    assert np.linalg.norm(generator-expected) < 3e-6
    assert np.linalg.norm(generator-reversed_phase) > .2


def test_every_source_environment_draw_and_native_word_is_counted():
    result = s.run(2, .7, .1, [0, 3], phi=.2, subdivisions=3, quadrature_order=4)
    operations = result['operations']
    assert result['independent_environment_draws_per_step'] == 28
    assert result['records'][-1]['independent_environment_draws'] == 84
    for key in native.COUNT_KEYS:
        unitary = 6*sum(operations[name]['resources'][key] for name in ('positive_subpulse', 'negative_subpulse'))
        unitary += 2*sum(operations[name]['resources'][key] for name in ('spin_z', 'positive_spin_gap', 'negative_spin_gap'))
        env = 24*operations['sideband_half_environment']['resources_per_realization'][key]
        env += 4*operations['spin_gap_environment']['resources_per_realization'][key]
        assert result['resources_per_evolution_step'][key] == unitary+env
        assert result['environment_resources_for_all_branch_contractions_per_step'][key] == 4*env
        assert result['records'][-1]['evolution_resources_per_realization'][key] == 3*(unitary+env)
    budget = result['records'][-1]['execution_error_budget']
    assert [term['invocations'] for term in budget['environment']['terms']] == [72, 12]
    assert not result['analytic_diffusion_channel_inserted']
    assert result['all_environment_phase_unitaries_source_executed']


def test_coupled_spin_reference_coherences_survive_gap_execution():
    u, _ = s.spin_gap(4, .2, DELTA, CARRIER)
    state = np.zeros(8, complex)
    state[0] = state[7] = 1/math.sqrt(2)
    actual = u @ np.outer(state, state.conj()) @ u.conj().T
    assert np.trace(actual @ actual).real == pytest.approx(1., abs=2e-12)
    assert np.linalg.norm(actual[:4, 4:]) > .49


def test_asymmetric_finite_gap_split_converges_to_unsplit_reference():
    rho = reference.preparation(4, .7, .1)
    for _ in range(3):
        rho = independent_step(rho)
    expected = reference.readout(rho, SETTINGS)['spin_probabilities'][0]
    errors = []
    for subdivisions in (4, 8, 16):
        actual = s.run(4, .7, .1, [3], SETTINGS, phi=.2,
                       subdivisions=subdivisions)['records'][0]
        errors.append(abs(actual['spin_probabilities'][0]-expected))
    assert errors[1] < errors[0]/3.5
    assert errors[2] < errors[1]/3.5
    assert errors[-1] < 1.3e-6


@pytest.mark.parametrize('kwargs', [{'finite_spin_duration': 1}, {'omega_carrier': 0},
    {'omega_carrier': float('inf')}, {'phi': True}, {'subdivisions': 0},
    {'quadrature_order': 0}, {'dephasing_rate': -1}, {'pulse_duration': 0}])
def test_invalid_apparatus_rejected(kwargs):
    with pytest.raises(ValueError):
        s.run(2, .7, .1, [0], **kwargs)


@pytest.mark.parametrize('steps', [[], [1, 1], [2, 1], [-1], [True], None])
def test_invalid_step_schedules_rejected(steps):
    with pytest.raises(ValueError):
        s.run(2, .7, .1, steps)


def test_observed_values_cannot_enter_execution_settings():
    with pytest.raises(ValueError):
        s.run(2, .7, .1, [0], [{**s.position_setting(), 'observed': .2}])

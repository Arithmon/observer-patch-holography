"""Independent circuit and intervention checks for executed field probes."""
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm, sqrtm

from model.observer_reconstruction import recovery, source


def _dense_two_mode_probability(states, controls, phase):
    """Explicit full tensor circuit, with mode zero least significant."""
    d0, d1 = map(len, states)
    hadamard = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
    state = np.kron(hadamard[:, 0], np.kron(states[1], states[0]))
    for index, control in enumerate(controls):
        dimension = (d0, d1)[index]
        blocks = control.reshape(2, dimension, 2, dimension)
        full = np.block([[np.kron(np.eye(d1), blocks[a, :, b, :]) if index == 0
                          else np.kron(blocks[a, :, b, :], np.eye(d0))
                          for b in range(2)] for a in range(2)])
        state = full@state
    read = hadamard@np.diag([1, np.exp(1j*phase)])
    state = np.kron(read, np.eye(d0*d1))@state
    return np.sum(abs(state.reshape(2, d0*d1))**2, axis=1)


@pytest.mark.parametrize('phase', [0., np.pi/2, -.37])
def test_complete_pointer_transfer_matches_dense_circuit_with_off_block_mixing(phase):
    # A mode interacts once with the common pointer and may then be traced.
    # Deliberate pointer-flip terms prevent a diagonal-block shortcut passing.
    rng = np.random.default_rng(91)
    states, controls = [], []
    for dimension in (2, 3):
        state = rng.normal(size=dimension)+1j*rng.normal(size=dimension)
        states.append(state/np.linalg.norm(state))
        generator = rng.normal(size=(2*dimension, 2*dimension))
        generator = generator+generator.T
        controls.append(expm(.27j*generator))
    actual = source.pointer_probabilities(states, controls, phase)
    expected = _dense_two_mode_probability(states, controls, phase)
    assert np.allclose(actual, expected, atol=2e-13)
    assert actual.sum() == pytest.approx(1., abs=1e-13)
    diagonal_only = []
    for control, state in zip(controls, states, strict=True):
        d = len(state)
        incomplete = control.copy()
        incomplete[:d, d:] = 0
        incomplete[d:, :d] = 0
        diagonal_only.append(incomplete)
    assert np.linalg.norm(_dense_two_mode_probability(states, diagonal_only, phase)-expected) > .01


def test_executed_two_mode_controls_match_independent_joint_pointer():
    states, controls = [], []
    for omega, kick, angle in ((1.2, .17, .31), (1.8, -.21, -.27)):
        prepared = source.coordinate(8, omega, kick)[0][:, 0]
        states.append(source.evolution(8, omega, .39)[0]@prepared)
        controls.append(source.coordinate(8, omega, angle, controlled=True)[0])
    for phase in (0., np.pi/2):
        assert np.allclose(source.pointer_probabilities(states, controls, phase),
                           _dense_two_mode_probability(states, controls, phase), atol=2e-13)


def test_walsh_modes_diagonalize_independently_assembled_hidden_cube():
    positions = np.array(list(product((0, 1), repeat=3)))[list(source.PERMUTATION)]
    adjacency = (np.sum(abs(positions[:, None]-positions[None, :]), axis=-1) == 1).astype(int)
    target = source.MASS**2*np.eye(8)+source.KAPPA*(3*np.eye(8)-adjacency)
    frequencies, basis = source.modes()
    assert np.allclose(basis.T@basis, np.eye(8), atol=1e-14)
    assert np.allclose((basis*frequencies**2)@basis.T, target, atol=1e-14)
    assert np.allclose(source.supplied_target(), target, atol=1e-14)
    assert np.allclose(np.unique(frequencies**2),
                       source.MASS**2+source.KAPPA*np.array([0., 2., 4., 6.]))


@pytest.mark.parametrize('dimension,duration,kick_site,read_site,gamma,eta', [
    (8, .37, 0, 0, .3, .4),
    (8, .61, 2, 5, -.47, .23),
    (8, .23, 7, 1, .19, -.37),
    (16, .73, 3, 6, .45, .55),
])
def test_source_probability_matches_independent_vacuum_characteristic(
        dimension, duration, kick_site, read_site, gamma, eta):
    probabilities, metadata = source.field_readout(dimension, duration, kick_site,
                                                   read_site, gamma, eta)
    matrix = source.supplied_target()
    root = sqrtm(matrix)
    inverse_root = np.linalg.inv(root)
    # Independent dense matrix functions, without the backend's Walsh modes.
    response = (inverse_root@(expm(1j*duration*root)-expm(-1j*duration*root))/(2j)).real
    variance = inverse_root[read_site, read_site]/2
    expected = np.exp(-eta**2*variance/2+1j*eta*gamma*response[read_site, kick_site])
    actual = recovery.characteristic(probabilities[:, 0])
    assert actual == pytest.approx(expected, abs=2e-10)
    assert recovery.predict_chi(matrix, duration, kick_site, read_site, gamma=gamma, eta=eta) == pytest.approx(expected, abs=1e-13)
    assert np.allclose(probabilities.sum(axis=1), 1, atol=2e-13)
    assert metadata['maximum_local_target_residual'] < 2e-10
    resources = metadata['resources_for_two_fresh_phase_trials']
    assert resources['native_pulses'] > 0
    assert resources['proper_code_service_events'] == 2*resources['native_pulses']


@pytest.mark.parametrize('time', [.13, .37, 1.1])
def test_executed_clock_has_correct_orientation_and_can_be_inverted(time):
    probabilities, _ = source.clock_readout(time)
    assert np.allclose(probabilities[:, 0],
                       [(1+np.cos(time))/2, (1-np.sin(time))/2], atol=2e-13)
    assert recovery.clock_time(probabilities[:, 0]) == pytest.approx(time, abs=2e-13)
    assert np.allclose(probabilities.sum(axis=1), 1, atol=2e-13)


def test_weyl_loop_detects_order_and_commuting_control():
    forward = source.weyl_loop(8, a=.3, b=.4)
    backward = source.weyl_loop(8, reverse=True, a=.3, b=.4)
    commuting = source.weyl_loop(8, a=0., b=.4)
    assert recovery.characteristic(forward[:, 0]) == pytest.approx(np.exp(.12j), abs=1e-10)
    assert recovery.characteristic(backward[:, 0]) == pytest.approx(np.exp(-.12j), abs=1e-10)
    assert recovery.characteristic(commuting[:, 0]) == pytest.approx(1., abs=1e-12)
    assert abs(forward[1, 0]-backward[1, 0]) > .1
    assert np.allclose(forward.sum(axis=1), 1, atol=1e-12)


def test_source_corruption_cannot_be_hidden_by_analytic_targets(monkeypatch):
    original = source.c.execute

    def corrupted(word, bits):
        matrix = original(word, bits)
        # Left multiplication by a nontrivial row phase keeps unitarity but
        # changes the full source operation; target checks must catch it.
        matrix[0] *= -1
        return matrix

    source._coordinate.cache_clear()
    monkeypatch.setattr(source.c, 'execute', corrupted)
    try:
        with pytest.raises(ArithmeticError, match='full-phase target'):
            source.field_readout(8, .37, 0, 0, .311, .411)
    finally:
        source._coordinate.cache_clear()


def test_audit_target_is_not_a_probability_oracle(monkeypatch):
    def forbidden():
        raise AssertionError('audit target used as a probability oracle')

    monkeypatch.setattr(source, 'supplied_target', forbidden)
    probabilities, _ = source.field_readout(8, .31, 1, 4, .17, .23)
    assert np.allclose(probabilities.sum(axis=1), 1, atol=1e-12)


def test_positive_and_negative_kicks_retain_conjugate_characteristics():
    positive, _ = source.field_readout(8, .43, 0, 0, .39, .41)
    negative, _ = source.field_readout(8, .43, 0, 0, -.39, .41)
    left, right = map(lambda p: recovery.characteristic(p[:, 0]), (positive, negative))
    assert left == pytest.approx(np.conj(right), abs=1e-12)
    assert left.imag > .01


@pytest.mark.parametrize('dimension,omega,angle,controlled', [
    (4, 1., .1, False), (True, 1., .1, False), (8, 0., .1, False),
    (8, 1., np.nan, False), (8, 1., .1, 1),
])
def test_invalid_coordinate_controls_are_rejected(dimension, omega, angle, controlled):
    with pytest.raises(ValueError):
        source.coordinate(dimension, omega, angle, controlled=controlled)

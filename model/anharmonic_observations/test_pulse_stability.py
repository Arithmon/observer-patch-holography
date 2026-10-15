"""Inverse-word cancellation, literal parity phase and unchanged raw guards."""
from copy import deepcopy
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import source as s
from ..scalar_observations.experimental_readout import executable_source as native


DELTA, OMEGA = 2*math.pi*500, math.pi/.000150


@pytest.mark.parametrize('dimension', [2, 4, 8, 32, 64])
def test_production_controlled_parity_is_two_pulses_with_full_phase(dimension):
    actual, metadata = s._controlled_parity(dimension)
    wanted = np.diag(np.r_[np.ones(dimension), (-1.)**np.arange(dimension)])
    assert np.linalg.norm(actual-wanted, 2) < 2e-14
    assert np.linalg.norm(actual@actual-np.eye(2*dimension), 2) < 3e-14
    assert metadata['resources']['native_pulses'] == 2
    assert metadata['resources']['pair_code_pulses'] == 1


@pytest.fixture(autouse=True)
def separate_pulse_cache():
    s.pulse_unitary.cache_clear()
    yield
    s.pulse_unitary.cache_clear()


def _perturb_inverse(monkeypatch, phase):
    original = s._forced_basis
    def nearby(*args):
        values, basis, inverse, forward_meta, inverse_meta = original(*args)
        return values, basis, inverse*np.exp(1j*phase), forward_meta, inverse_meta
    monkeypatch.setattr(s, '_forced_basis', nearby)


@pytest.mark.parametrize('sign', [-1, 1])
def test_guarded_inverse_contraction_roundoff_does_not_accumulate_at_zero_time(monkeypatch, sign):
    # A deterministic nearby contraction models the measured platform drift;
    # the exact compiled word and its separately executed inverse are unchanged.
    _perturb_inverse(monkeypatch, 2e-11)
    actual, metadata = s.pulse_unitary(4, DELTA, OMEGA, 0., sign)
    assert np.linalg.norm(actual-np.eye(8), 2) < 5e-14
    diagnostics = metadata['inverse_word_cancellation']
    assert 1e-11 < diagnostics['raw_full_phase_target_residual'] < 3e-11
    assert 1e-11 < diagnostics['representation_change'] < 3e-11
    assert not diagnostics['outward_roundoff_or_physical_error_certificate']
    assert not diagnostics['physical_inverse_service_errors_cancelled']


def test_tiny_pulses_keep_their_actual_generator_after_inverse_cancellation(monkeypatch):
    _perturb_inverse(monkeypatch, 2e-11)
    duration, repetitions, dimension = 1e-12, 10000, 4
    actual, _ = s.pulse_unitary(dimension, DELTA, OMEGA, duration, -1)
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    hamiltonian = DELTA*np.kron(np.eye(2), np.diag(np.arange(dimension)))
    hamiltonian += OMEGA/2*np.kron(np.array([[0., 1.], [1., 0.]]), ladder+ladder.T)
    target = expm(-1j*duration*repetitions*hamiltonian)
    assert np.linalg.norm(np.linalg.matrix_power(actual, repetitions)-target, 2) < 5e-11
    assert np.linalg.norm(target-np.eye(8), 2) > 1e-4


def _scale_inverse(monkeypatch, excess):
    original = s._forced_basis
    def nearby(*args):
        values, basis, inverse, forward_meta, inverse_meta = original(*args)
        return values, basis, (1+excess)*inverse, forward_meta, inverse_meta
    monkeypatch.setattr(s, '_forced_basis', nearby)


def test_guarded_polar_representation_removes_accumulating_nonunitary_drift(monkeypatch):
    # Independently contracted inverse words can leave a small Hermitian
    # defect after cancellation even when both source target guards pass.
    # Correct the pulse representation once, without normalizing any state.
    _, baseline = s.pulse_unitary(4, DELTA, OMEGA, .00005, 1)
    s.pulse_unitary.cache_clear()
    _scale_inverse(monkeypatch, 2e-11)
    actual, metadata = s.pulse_unitary(4, DELTA, OMEGA, .00005, 1)
    diagnostics = metadata['inverse_word_cancellation']
    guard = diagnostics['post_cancellation_unitary_stabilization']
    assert 1e-11 < guard['raw_unitarity_defect'] < 1e-10
    assert 1e-12 < guard['projection_change'] <= guard['maximum_projection_change']
    assert guard['restored_unitarity_defect'] < 1e-13
    assert diagnostics['complete_representation_change'] <= diagnostics['representation_change_limit']
    assert not guard['physical_error_certificate']
    repeated = np.linalg.matrix_power(actual, 10000)
    assert np.linalg.norm(repeated.conj().T@repeated-np.eye(8), 2) < 1e-10
    assert metadata['resources'] == baseline['resources']
    assert [stage['word_sha256'] for stage in metadata['stages']] == [
        stage['word_sha256'] for stage in baseline['stages']]


def test_polar_representation_does_not_hide_excessive_nonunitary_drift(monkeypatch):
    # This still passes the raw 3e-10 target guard, but violates the existing
    # fixed 1e-10 polar-input defect limit. The limit is not enlarged.
    _scale_inverse(monkeypatch, 1e-10)
    with pytest.raises(ArithmeticError, match='executed-unitary defect'):
        s.pulse_unitary(4, DELTA, OMEGA, .00005, 1)


def test_representation_change_has_its_own_fixed_guard(monkeypatch):
    _perturb_inverse(monkeypatch, 2.5e-10)
    with pytest.raises(ArithmeticError, match='source pulse.*representation'):
        s.pulse_unitary(4, DELTA, OMEGA, 0., 1)


@pytest.mark.parametrize('sign', [-1, 1])
def test_repeated_inverse_increment_keeps_complete_probability_mass(monkeypatch, sign):
    # The perturbed inverse is still unitary: its phase is harmless in the
    # straight product but makes I+B(D-I)Binverse slightly nonunitary.
    _perturb_inverse(monkeypatch, 2e-11)
    dimension, duration, repetitions = 4, .00005/8, 1024
    actual, metadata = s.pulse_unitary(dimension, DELTA, OMEGA, duration, sign)
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    hamiltonian = DELTA*np.kron(np.eye(2), np.diag(np.arange(dimension)))
    hamiltonian -= sign*OMEGA/2*np.kron(np.array([[0., 1.], [1., 0.]]), ladder+ladder.T)
    _, eigenvectors = np.linalg.eigh(hamiltonian)
    initial = eigenvectors[:, -1]
    continued = np.linalg.matrix_power(actual, repetitions) @ initial
    assert abs(float(np.vdot(continued, continued).real)-1.) < 2e-10
    assert np.linalg.norm(actual-expm(-1j*duration*hamiltonian), 2) < 3e-11
    assert metadata['inverse_word_cancellation']['representation_change'] < 2e-10


@pytest.mark.parametrize('corruption', ['wrong_phase', 'nonfinite', 'over_total_change'])
def test_final_unitary_restoration_is_checked_against_actual_source(monkeypatch, corruption):
    # Warm all primitive factors so only the new complete-pulse polar call
    # receives an adversarial return. The per-word source checks stay live.
    s.pulse_unitary(4, DELTA, OMEGA, .00005/8, 1)
    s.pulse_unitary.cache_clear()
    original = s.c.stabilize_executed_unitary
    seen = []
    def manufactured(value):
        corrected, metadata = original(value)
        if value.shape == (8, 8):
            seen.append(True)
            if corruption == 'nonfinite':
                corrected = corrected.copy()
                corrected[0, 0] = float('nan')
            else:
                phase = 1e-6 if corruption == 'wrong_phase' else 2.5e-10
                corrected = corrected*np.exp(1j*phase)
        return corrected, metadata
    monkeypatch.setattr(s.c, 'stabilize_executed_unitary', manufactured)
    with pytest.raises(ArithmeticError, match='source pulse'):
        s.pulse_unitary(4, DELTA, OMEGA, .00005/8, 1)
    assert seen


@pytest.mark.parametrize('sign', [-1, 1])
def test_literal_parity_square_is_not_replaced_with_identity(monkeypatch, sign):
    parity, meta = s._controlled_parity(4)
    phase = 2e-11
    changed = np.exp(1j*phase)*parity
    monkeypatch.setattr(s, '_controlled_parity', lambda _: (changed, deepcopy(meta)))
    actual, metadata = s.pulse_unitary(4, DELTA, OMEGA, 0., sign)
    # The applied parity is repeated, not followed by its inverse word.
    expected = np.exp(2j*phase)*np.eye(8)
    assert np.linalg.norm(actual-expected, 2) < 5e-14
    assert np.linalg.norm(actual-np.eye(8), 2) > 3e-11
    assert metadata['inverse_word_cancellation']['literal_parity_square_retained']


@pytest.mark.parametrize('corruption', ['inverse', 'parity', 'axis', 'nonfinite'])
def test_stable_identity_cannot_hide_bad_source_contractions(monkeypatch, corruption):
    if corruption == 'inverse':
        _perturb_inverse(monkeypatch, 1e-6)
    elif corruption == 'nonfinite':
        original = s._forced_basis
        def corrupted(*args):
            values, basis, inverse, forward_meta, inverse_meta = original(*args)
            inverse = inverse.copy()
            inverse[0, 0] = float('nan')
            return values, basis, inverse, forward_meta, inverse_meta
        monkeypatch.setattr(s, '_forced_basis', corrupted)
    elif corruption == 'parity':
        parity, meta = s._controlled_parity(4)
        monkeypatch.setattr(s, '_controlled_parity', lambda _: (np.exp(1e-6j)*parity, meta))
    else:
        original = s.spin_rotation
        def corrupted(axis, angle):
            value, metadata = original(axis, angle)
            return (np.exp(1e-6j)*value if axis == 'Y' and angle > 0 else value), metadata
        monkeypatch.setattr(s, 'spin_rotation', corrupted)
    with pytest.raises(ArithmeticError, match='source pulse'):
        s.pulse_unitary(4, DELTA, OMEGA, 0., 1)


def test_stabilization_keeps_every_source_stage_and_resource():
    _, metadata = s.pulse_unitary(4, DELTA, OMEGA, .00005/256, -1)
    assert [stage['name'] for stage in metadata['stages']] == [
        'spin_rotation_Z', 'spin_rotation_Y', 'spin_controlled_oscillator_parity',
        'forced_oscillator_inverse_basis', 'forced_eigenphases', 'forced_oscillator_basis',
        'spin_controlled_oscillator_parity', 'spin_rotation_Y', 'spin_rotation_Z']
    for key in native.COUNT_KEYS:
        assert metadata['resources'][key] == sum(stage['resources'][key] for stage in metadata['stages'])
    assert all(len(stage['word_sha256']) == 64 for stage in metadata['stages'])

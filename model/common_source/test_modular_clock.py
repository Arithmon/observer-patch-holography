"""Actual-history clock instruments, sign checks and failed clock controls."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import experiment as e
from . import modular_clock as m
from . import operational as o
from ..observer.scalar_response.experiment import canonical, compare_receipt


@pytest.fixture(scope='module')
def history():
    return e.history(2, 16)


def test_clock_rate_is_the_actual_reduced_state_gap(history):
    rho, mass = o.scalar_state(history)
    clock = m.parameters(history)
    eigenvalues = np.linalg.eigvalsh(rho)
    assert clock['modular_frequency'] == pytest.approx(math.log(eigenvalues[1]/eigenvalues[0]), abs=2e-13)
    assert clock['contrast'] == pytest.approx(eigenvalues[1]-eigenvalues[0], abs=2e-13)
    assert clock['mass'] == pytest.approx(mass, abs=2e-13)
    assert clock['p_minus'] == pytest.approx(history['second_record_marginal'][0], abs=2e-13)
    assert clock['p_plus'] == pytest.approx(history['second_record_marginal'][1], abs=2e-13)


def test_clock_receipt_preserves_scalar_types_through_json(history):
    result = m.receipt(history)
    compare_receipt(json.loads(canonical(result)), result)


def test_source_instrument_against_independent_modular_and_pointer_matrices(history):
    rho, _ = o.scalar_state(history)
    k = o.modular_hamiltonian(rho)
    z = np.diag([1., -1.])
    n = np.diag([0., 1.])
    s = .73
    preparation = expm(1j*math.pi*z/4)
    u = expm(1j*s*k)@preparation
    # The Schrödinger unitary exp(+is K) is the dual of the declared sigma_s.
    assert np.trace(u@rho@u.conj().T@n) == pytest.approx(
        np.trace(preparation@rho@preparation.conj().T@o.modular_flow(rho, n, s)), abs=2e-13)
    full = np.kron(np.kron(np.eye(2), u), np.eye(4))
    vectors = o.terminal_vectors(history)
    row = m.sample_setting(history, s)
    for branch in row['branches']:
        bit = branch['clock_bit']
        projector = np.diag([int((i >> 2 & 1) == bit) for i in range(16)])
        target = projector@full@vectors[tuple(branch['old_record'])]
        values = np.asarray(branch['unnormalized_terminal_state'])
        actual = values[:, 0]+1j*values[:, 1]
        # The source omits only the scalar phase of K: compare branch densities.
        np.testing.assert_allclose(np.outer(actual, actual.conj()), np.outer(target, target.conj()), atol=2e-13)
    assert row['clock_bit_marginal'][1] == pytest.approx(np.trace(u@rho@u.conj().T@n).real, abs=2e-13)


def test_all_old_records_and_quantum_branches_retained_with_charged_resources(history):
    row = m.sample_setting(history, .7)
    assert len(row['branches']) == 8
    vectors = o.terminal_vectors(history)
    for record, vector in vectors.items():
        subtotal = sum(x['probability'] for x in row['branches'] if tuple(x['old_record']) == record)
        assert subtotal == pytest.approx(np.vdot(vector, vector).real, abs=2e-13)
    appended = row['appended_source_resources_per_realized_history']
    pointer = row['execution']['clock_pointer']['resources']
    total = row['total_source_resources_per_realized_history']
    assert appended['native_pulses'] == 2+pointer['native_pulses']
    assert appended['proper_code_service_events'] == 2*appended['native_pulses']
    assert total['native_pulses'] == history['resources_per_realized_history']['native_pulses']+appended['native_pulses']
    assert total['fresh_pointer_blanks'] == total['committed_classical_bits'] == 3
    # A zero-angle pulse is still an executed request, not a free service.
    zero = m.sample_setting(history, 0., prepared=False)
    assert zero['appended_source_resources_per_realized_history']['native_pulses'] == appended['native_pulses']


def test_monotone_clock_inversion_and_genuine_phase_ambiguity(history):
    receipt = m.receipt(history)
    assert receipt['all_pass']
    clock = receipt['clock_parameters']
    for row in receipt['settings']:
        assert m.inferred_parameter(row['clock_bit_marginal'][1], clock) == pytest.approx(row['modular_parameter'], abs=1e-9)
    first = receipt['settings'][2]
    for control in (receipt['wrapped_parameter_control'], receipt['outside_window_mirror_control']):
        assert control['modular_parameter'] > clock['monotone_parameter_window'][1]
        assert control['clock_bit_marginal'] == pytest.approx(first['clock_bit_marginal'], abs=2e-13)
    assert all(abs(x['clock_bit_marginal'][1]-clock['mass']/2) < 2e-13
               for x in receipt['unprepared_stationary_control'])
    assert not receipt['scope']['one_binary_outcome_is_an_exact_time_reading']
    assert not receipt['scope']['physical_time_calibration_established']
    assert not receipt['scope']['original_hamiltonian_model_time_identification_established']


def test_actual_feedback_changes_state_dependent_clock_rate_and_contrast(history):
    changed = e.history(2, 16, feedback_enabled=False)
    old, new = m.parameters(history), m.parameters(changed)
    assert abs(old['modular_frequency']-new['modular_frequency']) > .01
    assert abs(old['contrast']-new['contrast']) > .01
    # Compare the same modular parameter, without refitting a desired period.
    first, second = m.sample_setting(history, .7), m.sample_setting(changed, .7)
    assert abs(first['clock_bit_marginal'][1]-second['clock_bit_marginal'][1]) > .005


def test_disconnected_actual_history_has_no_modular_clock_contrast():
    disconnected = e.history(2, 16, coupling=0.)
    rho, _ = o.scalar_state(disconnected)
    np.testing.assert_allclose(rho, np.eye(2)/2, atol=2e-13)
    with pytest.raises(ValueError, match='no nontrivial modular clock contrast'):
        m.receipt(disconnected)


@pytest.mark.parametrize('parameter', (True, float('nan'), float('inf'), 1j, '1'))
def test_rejects_invalid_modular_parameters(history, parameter):
    with pytest.raises(ValueError, match='finite real'):
        m.sample_setting(history, parameter)


def test_rejects_foreign_or_incomplete_clock_histories_and_invalid_window(history):
    incomplete = deepcopy(history)
    incomplete['branches'].pop()
    with pytest.raises(ValueError, match='complete'):
        m.receipt(incomplete)
    foreign = deepcopy(history)
    foreign['dimension'] = 4
    with pytest.raises(ValueError, match='d=2'):
        m.receipt(foreign)
    with pytest.raises(ValueError, match='Boolean'):
        m.sample_setting(history, 1, prepared=1)
    clock = m.parameters(history)
    for probability in (0., 1.):
        with pytest.raises(ValueError, match='outside'):
            m.inferred_parameter(probability, clock)

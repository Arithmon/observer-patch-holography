"""Complete recorded instruments, physical countercontrols and shared outputs."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import experiment as e
from .test_scalar_coupling import independent_terms


def decode(state):
    values = np.asarray(state)
    return values[:, 0]+1j*values[:, 1]


@pytest.fixture(scope='module')
def histories():
    return {(d, current): e.history(d, 16, current=current)
            for d in (2, 4) for current in (-1, 0, 1)}


@pytest.mark.parametrize('d', (2, 4))
def test_pointer_is_complete_qnd_instrument_on_entangled_reference(d):
    # Q assembled from physical oscillator matrix elements, independently of
    # the compiler and its chosen eigenvectors.
    q = np.zeros((d, d))
    for n in range(d-1):
        q[n, n+1] = q[n+1, n] = np.sqrt((n+1)/2)
    values, basis = np.linalg.eigh(q)
    read, _ = e.scalar_pointer(d)
    rng = np.random.default_rng(411+d)
    reference = rng.normal(size=(4*d*d, 3))+1j*rng.normal(size=(4*d*d, 3))
    reference /= np.linalg.norm(reference)
    for bit, k in enumerate(read):
        mask = values > 0 if bit else values < 0
        expected = np.kron(np.kron(np.eye(d), basis[:, mask]@basis[:, mask].T), np.eye(4))
        np.testing.assert_allclose(k@reference, expected@reference, atol=2e-13)
        np.testing.assert_allclose(k@k, k, atol=2e-13)
    np.testing.assert_allclose(sum(k.conj().T@k for k in read), np.eye(4*d*d), atol=2e-13)


def test_full_record_law_and_terminal_states_against_continuous_hamiltonian(histories):
    # An independent continuous target tests more than replaying a split word.
    # Charge both recorded interventions against the exact temporal bound.
    for (d, current), result in histories.items():
        size = 4*d*d
        h = sum(independent_terms(d).values())
        scalar = np.zeros(d*d)
        scalar[0] = 1
        if d == 4:
            scalar[2] = scalar[2*d] = -1/(18*np.sqrt(2))
            scalar[1+d] = 1/18
        scalar /= np.linalg.norm(scalar)
        initial = np.kron(scalar, [1, 0, 0, 0])
        x = np.zeros((size, size))
        for i in range(size):
            x[i ^ 1, i] = 1
        initial = expm(-1j*current*e.FEEDBACK_ANGLE*x)@initial
        first, second = expm(-1j*h), expm(-.5j*h)
        total_trace_distance = 0.
        for row in result['branches']:
            a, b = row['record']
            pa = np.diag([int((i & 1) == a) for i in range(size)])
            feedback = expm(-1j*a*e.FEEDBACK_ANGLE*x)
            target = e.scalar_projectors(d)[b]@second@feedback@pa@first@initial
            actual = decode(row['unnormalized_terminal_state'])
            difference = np.outer(actual, actual.conj())-np.outer(target, target.conj())
            total_trace_distance += np.sum(abs(np.linalg.eigvalsh(difference)))/2
            assert row['probability'] == pytest.approx(float(np.vdot(actual, actual).real), abs=1e-14)
        certificate = F(result['ideal_finite_continuous_comparison']['complete_instrument_half_diamond_error_upper'])
        assert 0 < total_trace_distance < float(certificate)
        assert sum(row['probability'] for row in result['branches']) == pytest.approx(1, abs=2e-11)


@pytest.mark.parametrize('d', (2, 4))
def test_interaction_transmits_feedback_but_disconnected_scalar_cannot(d, histories):
    enabled = histories[d, 0]
    disabled = e.history(d, 16, feedback_enabled=False)
    assert enabled['second_record_marginal'][1]-disabled['second_record_marginal'][1] < -.02
    np.testing.assert_allclose(enabled['first_record_marginal'], disabled['first_record_marginal'], atol=2e-12)
    decoupled = [e.history(d, 16, coupling=0., feedback_enabled=flag) for flag in (False, True)]
    np.testing.assert_allclose(decoupled[0]['second_record_marginal'], decoupled[1]['second_record_marginal'], atol=2e-12)
    # Coarse feedback reaches an interacting scalar even though its first
    # committed matter record cannot change retrospectively.
    assert abs(enabled['terminal_gauge_scalar_covariance']) > .001


def test_resources_include_zero_controls_and_every_service(histories):
    for (d, _), result in histories.items():
        execution = result['execution']
        first, second = execution['feedback']
        assert first['angle'] == 0 and second['angle'] != 0
        assert first['resources']['native_pulses'] == second['resources']['native_pulses'] == 1
        rows = [execution[key]['resources'] for key in (
            'preparation', 'native_current_preparation', 'first_evolution',
            'second_evolution', 'first_pointer', 'second_pointer')]+[first['resources']]
        resources = result['resources_per_realized_history']
        for key in e.COUNT_KEYS:
            assert resources[key] == sum(row[key] for row in rows)
        assert resources['proper_code_service_events'] == 2*resources['native_pulses']
        assert resources['fresh_pointer_blanks'] == 2
        assert resources['system_qubits'] == 2+2*(d.bit_length()-1)


def test_temporal_bound_covers_whole_instrument_and_refines_quadratically():
    a = e.history(2, 8)
    b = e.history(2, 16)
    # Use uncapped segment bounds here since these deliberately coarse words
    # can have a vacuous, capped full-instrument certificate.
    for old, new in zip(a['ideal_finite_continuous_comparison']['segment_certificates'],
                        b['ideal_finite_continuous_comparison']['segment_certificates']):
        assert F(old['unclipped_operator_error_upper']) == 4*F(new['unclipped_operator_error_upper'])


@pytest.mark.parametrize('kwargs', ({'current': True}, {'current': 2}, {'steps': 0},
                                   {'steps': 1.5}, {'feedback_enabled': 1}))
def test_protocol_rejects_undeclared_controls(kwargs):
    with pytest.raises(ValueError):
        e.history(2, **kwargs)

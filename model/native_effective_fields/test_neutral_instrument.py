"""Complete source-coded phase/spatial detector and its genuine boundaries."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import neutral_instrument as n


@pytest.mark.parametrize('law', ['C', 'E'])
def test_native_two_record_herald_is_pure_and_keeps_both_signs(law):
    weights, outputs = n.herald_audit(law)
    assert weights == (F(131, 132), F(1, 264), F(1, 264))
    # Arbitrary complex superposition entangled with a reference; every
    # successful branch has only its declared one-dimensional work support.
    v = np.arange(132)+1j*np.arange(132)[::-1]
    v /= np.linalg.norm(v)
    rho = np.outer(v, v.conj())
    for flag, group in enumerate(n.herald(law)):
        result = sum((np.kron(k, np.eye(2))@rho@np.kron(k.conj().T, np.eye(2)) for _, k in group),
                     np.zeros((132, 132), complex))
        if flag:
            allowed = [2*outputs[flag], 2*outputs[flag]+1]
            p = np.zeros((132, 132)); p[allowed, allowed] = 1
            assert np.linalg.norm(result-p@result@p) < 1e-14


def test_source_copy_words_use_actual_pointer_and_full_raw_certificate():
    a = n.copy_matrix()
    wanted = np.zeros((4, 4))
    for j in range(4):
        wanted[j ^ (2 if j & 1 else 0), j] = 1
    assert np.linalg.norm(a-wanted) < 1e-13
    assert n.copy_certificate()['half_diamond_upper'] < F(1, 10**12)
    assert sum(k.conj().T@k for k in n.pointer_kraus()) == pytest.approx(np.eye(2), abs=1e-13)
    # Execute emitted words on nonlocal bits with spectator/reference columns.
    for bit in (0, 4, 8):
        state = np.zeros((1024, 4), complex)
        rows = [0, 1 << bit, 512, 512+(1 << bit)]
        state[rows, np.arange(4)] = [1, 1j, -1, -.5j]
        actual = n.compiler.execute(n.copy_word(bit), 10, state)
        expected = np.zeros_like(state)
        for col, row in enumerate(rows):
            expected[row ^ (512 if row >> bit & 1 else 0), col] = state[row, col]
        assert np.linalg.norm(actual-expected) < 1e-13


def test_raw_copy_certificate_detects_all_column_errors_and_no_polar_repair():
    damaged = np.array(n.copy_matrix())
    damaged[3, 3] += .01j
    assert n.copy_certificate(damaged)['operator_error_upper'] >= F(1, 101)
    assert n.copy_certificate(2*np.array(n.copy_matrix()))['channel_norm_upper'] >= 4
    class Hostile(np.ndarray):
        pass
    with pytest.raises(ValueError):
        n.copy_certificate(np.array(n.copy_matrix()).view(Hostile))


def test_address_detector_matches_povm_but_not_complement_lueders():
    audit = n.analyzer_audit()
    assert audit['same_complete_21_outcome_POVM_residual_diagnostic'] < 1e-12
    assert audit['complement_update_trace_distance_countercontrol'] == pytest.approx(.5)
    assert audit['raw_pointer_readout_count'] == 512
    assert not audit['full_Lueders_instrument_claimed']
    # Every nonbinary or padding address is a retained complement outcome.
    active = [code for code in range(512) if n.coarse_address(code) < 20]
    assert len(active) == 20
    assert len({n.coarse_address(code) for code in active}) == 20
    assert n.coarse_address(511) == 20


def test_full_original_pointer_arrays_include_both_signs_and_inactive_states():
    for sign in (-1, 1):
        state = n.stopped.dense(n.neutral.preparation(n.PAIR, n.COMMON, sign))
        raw = n.raw_readout_probabilities(state)
        index = n.neutral.addresses(n.PAIR, n.COMMON)[0 if sign == 1 else 1]
        code = n.circuit.encode(n.quantum.STATES[index])
        assert len(raw) == 512 and raw[code] == pytest.approx(1, abs=1e-13)
        assert sum(raw) == pytest.approx(1, abs=1e-13)
    rho = np.eye(66)/66
    raw = n.raw_readout_probabilities(rho)
    assert sum(raw) == pytest.approx(1, abs=1e-13)
    assert sum(p for code, p in enumerate(raw) if n.coarse_address(code) == 20) == pytest.approx(F(46, 66))


@pytest.mark.parametrize('law', ['C', 'E'])
def test_capped_complete_execution_and_exact_postmeasurement_states(law):
    row = n.execution(law)
    assert row['complete_coarse_record_count'] == 127
    assert row['failed_herald_probability'] == '131/132'
    assert row['success_and_timeout_probability'] == '25/4752'
    assert row['runtime_probability_residual_diagnostic'] < 1e-12
    assert row['checkpoint_complete_distance_diagnostic'] < 1e-12
    assert row['raw_source_analyzer_L1_residual_diagnostic'] < 1e-12
    n.compare_receipt(n.strict_json.loads(n.canonical(row)), row)


def test_discrimination_keeps_herald_failure_and_finite_timeout():
    c, e = (n.exact_law(law) for law in ('C', 'E'))
    distance = sum(abs(p-e[key]) for key, p in c.items())/2
    complete_only = sum(abs(p-e[key]) for key, p in c.items() if key[0] and key[1])/2
    assert distance == F(29, 39600)
    assert complete_only == F(1, 23760)
    assert distance != F(4, 69)
    # A first-attempt moving channel is identical in C and E; only later
    # zero-current coherence can make completed-movement records differ.
    assert all(p == e[key] for key, p in c.items() if key[1] == 1)
    assert any(p != e[key] for key, p in c.items() if key[1] == 2)
    assert all(p >= 0 for p in c.values()) and sum(c.values()) == sum(e.values()) == 1


def test_error_budget_uses_complete_output_not_rare_herald_conditioning():
    budget = n.source_budget()
    assert budget['raw_source_half_diamond_error_upper'] < F(1, 10**12)
    assert budget['complete_output_error_upper'] < F(1, 10**6)
    assert budget['maximum_total_native_attempts'] == 373
    assert budget['maximum_source_pulses_per_trial'] == 95
    assert budget['maximum_proper_code_service_events'] == 190
    assert budget['pointer_blank_resets_per_success'] == budget['pointer_bit_reads_per_success'] == 9
    assert not budget['physical_control_reset_read_error_included']
    audit = n.finite_warmup_audit()
    assert audit['actual_full_flagged_law_error_diagnostic'] <= audit['actual_prior_error_diagnostic']


@pytest.mark.parametrize('bad', [True, -1, 9, .5, '1'])
def test_invalid_copy_bit(bad):
    with pytest.raises(ValueError):
        n.copy_word(bad)


@pytest.mark.parametrize('bad', [True, -1, 512, .5, '1'])
def test_invalid_full_address(bad):
    with pytest.raises(ValueError):
        n.coarse_address(bad)

"""Adversarial controls for the canonical protected observer integration."""
from dataclasses import replace
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import encoded_memory as memory, encoded_observer as observer
from . import quantum_native as q


@pytest.fixture(scope='module')
def histories():
    return {mode: q.execute('C', 2, 8, native_mode=mode)
            for mode in ('protected_qubit', 'passive_qubit')}


def _changed(history, index, **updates):
    branches = list(history.branches)
    branches[index] = replace(branches[index], **updates)
    return replace(history, branches=tuple(branches))


def _active(history):
    return next(i for i, b in enumerate(history.branches)
                if not b.abort and np.linalg.norm(b.native_map) > .1
                and np.linalg.norm(b.work_state) > .1)


def test_complete_memory_receipt_checks_both_modes(histories):
    for history in histories.values():
        assert q.audit(history)['complete_stacked_map_frobenius_error_upper'] < 2e-9
    result = observer.memory_witness(histories['protected_qubit'], histories['passive_qubit'])
    assert result['active_encoded_bell_fidelity_diagnostic'] > 1-2e-9
    assert result['evaluated_passive_decoder_error_lower_diagnostic'] > 1/16
    assert result['terminal_read_is_destructive']


def test_relabeled_passive_execution_cannot_certify_active_memory(histories):
    forged = replace(histories['passive_qubit'], native_mode='protected_qubit')
    with pytest.raises(ArithmeticError, match='independent targets'):
        q.audit(forged)
    with pytest.raises(ArithmeticError, match='logical channel'):
        observer.memory_witness(forged, histories['passive_qubit'])


@pytest.mark.parametrize('mutation', ['logical_phase', 'padding', 'lost_abort', 'norm_preserving_wrong_work'])
def test_complete_audit_rejects_probability_and_menu_blind_spots(histories, mutation):
    history = histories['protected_qubit']
    i = _active(history)
    if mutation == 'lost_abort':
        changed = replace(history, branches=tuple(b for b in history.branches if b.record[2]))
        error = ValueError
    elif mutation == 'norm_preserving_wrong_work':
        changed = _changed(history, i, work_state=np.roll(history.branches[i].work_state, 1, axis=0))
        error = ArithmeticError
    else:
        matrix = history.branches[i].native_map.copy()
        if mutation == 'logical_phase':
            matrix[:, 1] *= 1j
            np.testing.assert_array_equal(abs(matrix)**2, abs(history.branches[i].native_map)**2)
        else:
            matrix[511, 0] += .01
        changed = _changed(history, i, native_map=matrix)
        error = ArithmeticError
    with pytest.raises(error):
        q.audit(changed)


def test_bypassing_source_recovery_is_detected(monkeypatch):
    q._native_tree.cache_clear()
    monkeypatch.setattr(memory, 'recovery_unitary', lambda *args: np.eye(512))
    try:
        wrong = q.execute('C', 2, 8, native_mode='protected_qubit')
        with pytest.raises(ArithmeticError, match='independent targets'):
            q.audit(wrong)
    finally:
        q._native_tree.cache_clear()


def test_wrong_source_terminal_decoder_is_detected(histories, monkeypatch):
    monkeypatch.setattr(memory, 'decoding', lambda value: value)
    with pytest.raises(ArithmeticError, match='decoded logical read'):
        observer.memory_witness(histories['protected_qubit'], histories['passive_qubit'])


def test_executor_seam_change_cannot_redefine_its_oracle(monkeypatch):
    q._native_tree.cache_clear()
    q._ideal_encoded_tree.cache_clear()
    monkeypatch.setattr(q, 'final_seam', lambda mode, bit: (0, 1))
    try:
        wrong = q.execute('C', 2, 8, native_mode='passive_qubit')
        with pytest.raises(ArithmeticError, match='independent targets'):
            q.audit(wrong)
    finally:
        q._native_tree.cache_clear()
        q._ideal_encoded_tree.cache_clear()


@pytest.mark.parametrize('mode', ['protected_qubit', 'passive_qubit'])
def test_encoded_preparation_rejects_tiny_escape_into_another_valid_native_code(mode, monkeypatch):
    wrong = np.zeros(512, complex)
    wrong[0], wrong[2] = math.sqrt(1-1e-12), 1e-6
    assert 2 in q.CODES  # The old 78-dimensional support test missed this.
    assert abs(np.vdot(wrong, wrong)-1) < 1e-15
    monkeypatch.setattr(q.compiler, 'execute', lambda *args: wrong)
    with pytest.raises(ArithmeticError, match='escaped'):
        q.prepared_native('phase_0', native_mode=mode)


@pytest.mark.parametrize('bad', [True, None, 'protected', [], 2])
def test_native_mode_must_be_explicitly_admitted(bad):
    with pytest.raises(ValueError):
        q.execute(native_mode=bad)


@pytest.mark.parametrize('bad', [[1., False], [1e308+1e308j, 0.], [np.nan, 0.],
                                np.eye(78)[:, 0], [0., 0.], np.empty((2, 0))])
def test_encoded_input_is_not_truncated_coerced_or_overflow_normalized(histories, bad):
    with pytest.raises(ValueError):
        q.output(histories['protected_qubit'], bad)


@pytest.mark.parametrize('mode', ['protected_qubit', 'passive_qubit'])
def test_path_budget_charges_causal_source_schedule_and_abort_branches(mode):
    costs = observer.resources('C', 2, 8, native_mode=mode)
    record = costs['maximizing_scheduled_record']
    work = next(b for b in q._programs(2, 8, 'source')[
        q._record_current('C', record['first_coin'], record['first_raw'])].branches
        if list(b.record) == record['work_record'])
    seam = (2, 3) if work.record[2] == 0 else (0, 1)
    assert list(seam) == record['final_seam']
    expected = (3+memory.encoding_resources()['native_pulses']
                +q.routing.routed_repair_resources('C', (0, 1))['native_pulses']
                +dict(work.resources)['native_pulses']
                +q.routing.routed_repair_resources('C', seam)['native_pulses'])
    if mode == 'protected_qubit':
        expected += memory.recovery_resources('C', (0, 1), record['first_coin'], record['first_raw'])['native_pulses']
        expected += memory.recovery_resources('C', seam, record['last_coin'], record['last_raw'])['native_pulses']
    if not record['last_abort']:
        expected += memory.decoding_resources()['native_pulses']+1
    assert expected == costs['maximum_native_pulses_per_attempt']
    assert costs['maximum_proper_code_service_events_per_attempt'] == 2*expected
    assert costs['initial_blank_native_code_qubits'] == 8
    assert costs['blank_quantum_qubits_added_during_recovery'] == 0
    assert costs['archive_record_rewrites_during_recovery'] == 0


def test_new_native_services_remain_explicit_unknowns():
    program = q.protocol.continuation_program(2, 8)
    keys = observer.required_service_keys(program)
    result = q.bounds.service_error_budget(F(1, 100), required_stage_keys=keys)
    assert result['conditional_complete_error_upper'] is None
    assert set(keys) == set(result['missing_service_stage_keys'])
    assert {'native_encoder', 'first_native_recovery', 'last_native_recovery',
            'native_terminal_decoder_rotation_and_readout',
            'native_logical_storage_through_work'} <= set(keys)


@pytest.mark.parametrize('bad', [True, 2., 8, 0])
def test_exact_finite_obstruction_is_not_extrapolated(bad):
    with pytest.raises(ValueError):
        observer.passive_obstruction(bad)

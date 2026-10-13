"""Complete-channel rejection controls for the quantum native observer."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import quantum_native as q


@pytest.fixture(scope='module')
def history():
    return q.execute('C', 2, 16)


def changed_branch(history, index, **updates):
    branches = list(history.branches)
    branches[index] = replace(branches[index], **updates)
    return replace(history, branches=tuple(branches))


def active_index(history):
    return next(i for i, b in enumerate(history.branches)
                if not b.abort and np.linalg.norm(b.native_map) > .1
                and np.linalg.norm(b.work_state) > .1)


def test_full_reference_stable_map_audit(history):
    report = q.audit(history)
    assert report['full_outcome_count_including_abort'] == 770
    assert report['first_abort_outcome_count'] == 2
    assert report['final_abort_outcome_count'] == 192
    assert report['complete_stacked_map_frobenius_error_upper'] < 2e-8
    assert report['prepared_native_input_completeness_operator_error'] < 2e-9


@pytest.mark.parametrize('mutation', ['omit_zero_abort', 'duplicate', 'wrong_abort', 'bool_record',
                                     'missing_record', 'nan_native', 'nan_work', 'nan_prepared'])
def test_complete_menu_rejects_trash(history, mutation):
    i = active_index(history)
    b = history.branches[i]
    if mutation == 'omit_zero_abort':
        changed = replace(history, branches=tuple(x for x in history.branches if x.record[2]))
    elif mutation == 'duplicate':
        changed = replace(history, branches=(*history.branches, history.branches[0]))
    elif mutation == 'wrong_abort':
        changed = changed_branch(history, i, abort=True)
    elif mutation == 'bool_record':
        changed = changed_branch(history, i, record=(bool(b.record[0]), *b.record[1:]))
    elif mutation == 'missing_record':
        changed = changed_branch(history, i, record=(0, 0, (), None, None))
    else:
        field = {'nan_native': 'native_map', 'nan_work': 'work_map',
                 'nan_prepared': 'work_state'}[mutation]
        value = getattr(b, field).copy()
        value.flat[0] = np.nan
        changed = changed_branch(history, i, **{field: value})
    with pytest.raises(ValueError):
        q.audit(changed)


def test_audit_rejects_probability_invisible_native_relative_phase(history):
    # Preserve every population for the maximally mixed native input while
    # changing a surviving off-diagonal matrix unit.
    i = next(i for i, b in enumerate(history.branches)
             if not b.abort and np.count_nonzero(np.linalg.norm(b.native_map, axis=0) > .1) > 1
             and np.linalg.norm(b.work_state) > .1)
    b = history.branches[i]
    value = b.native_map.copy()
    column = np.flatnonzero(np.linalg.norm(value, axis=0) > .1)[0]
    value[:, column] *= 1j
    assert np.allclose(abs(value)**2, abs(b.native_map)**2, rtol=0, atol=0)
    with pytest.raises(ArithmeticError, match='complete source-native'):
        q.audit(changed_branch(history, i, native_map=value))


def test_prepared_output_is_checked_separately_from_kraus_map(history):
    i = active_index(history)
    value = history.branches[i].work_state.copy()
    # Exactly norm-preserving wrong work output, with the original map intact.
    value = np.roll(value, 1, axis=0)
    with pytest.raises(ArithmeticError, match='complete source-native'):
        q.audit(changed_branch(history, i, work_state=value))


def test_retained_padding_is_checked_instead_of_projected_away(history):
    i = active_index(history)
    value = history.branches[i].native_map.copy()
    padding = next(k for k in range(512) if k not in q.CODES)
    value[padding, 0] += .01
    with pytest.raises(ArithmeticError, match='complete source-native'):
        q.audit(changed_branch(history, i, native_map=value))


@pytest.mark.parametrize('bad', [True, np.ones(78, dtype=bool), [True]+[0.]*77,
                                np.zeros(78), np.ones((78, 0)), np.zeros((77, 1)),
                                [np.nan]+[0.]*77, [np.inf]+[0.]*77,
                                np.ma.array(np.eye(78)[:, 0], mask=False)])
def test_native_input_rejects_invalid_original_data(history, bad):
    with pytest.raises(ValueError):
        q.output(history, bad)


@pytest.mark.parametrize('law,dimension,steps', [('Q', 2, 16), (True, 2, 16),
                                               ('C', True, 16), ('C', 8, 16),
                                               ('C', 2, 0), ('C', 2, True)])
def test_program_domain_is_not_silently_coerced(law, dimension, steps):
    with pytest.raises(ValueError):
        q.execute(law, dimension, steps)


def test_complex_reference_is_retained_through_actual_terminal_read(history):
    original = np.zeros((78, 2), complex)
    original[0, 0], original[1, 1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    marginal = np.zeros((2, 2), complex)
    outcomes = 0
    for _, code, native, work in q.readout_outputs(history, original):
        assert native.shape == ((512, 2) if code is None else (1, 2))
        marginal += float(np.vdot(work, work).real)*(native.conj().T@native)
        outcomes += 1
    assert outcomes == 576*512+194
    np.testing.assert_allclose(marginal, original.conj().T@original, atol=2e-10, rtol=0)


def test_empty_or_unnormalized_joint_reference_is_rejected(history):
    for value in (np.zeros((78, 16, 0)), np.zeros((78, 16, 1)),
                  np.zeros((78, 15, 1)), np.zeros((78, 16))):
        with pytest.raises(ValueError):
            q.output_on_joint_input(history, value)


def test_resource_count_covers_routing_and_all_code_services():
    costs = q.resources('C', 2, 16)
    assert costs['maximum_native_pulses_per_attempt'] == (
        costs['field_path_maxima']['native_pulses']+5371+2198007+4)
    assert costs['maximum_proper_code_service_events_per_attempt'] == 2*costs['maximum_native_pulses_per_attempt']
    assert costs['fixed_width_terminal_record_bits_including_two_status_bits'] == 21
    assert costs['peak_code_qubits_excluding_six_level_processor'] == 15
    assert costs['field_preparations'] == 1


def test_service_unknowns_do_not_become_zero_error():
    program = q.protocol.continuation_program()
    keys = (*q.bounds.required_service_keys(program.operations), 'native_repair')
    result = q.bounds.service_error_budget(F(1, 100), required_stage_keys=keys)
    assert result['conditional_complete_error_upper'] is None
    assert 'native_repair' in result['missing_service_stage_keys']


@pytest.mark.parametrize('mutation', ['wrong_count', 'wrong_flag', 'nonfinite', 'erased_error',
                                     'unknown_field', 'missing_field'])
def test_receipt_comparison_rejects_changed_science(mutation):
    actual = {'count': 770, 'full_quantum': True, 'error': 0.01, 'coverage': None}
    recorded = deepcopy(actual)
    if mutation == 'wrong_count':
        recorded['count'] = 769
    elif mutation == 'wrong_flag':
        recorded['full_quantum'] = 1
    elif mutation == 'nonfinite':
        recorded['error'] = float('nan')
    elif mutation == 'erased_error':
        recorded['error'] = 0.
    elif mutation == 'unknown_field':
        recorded['all_pass'] = True
    else:
        del recorded['coverage']
    with pytest.raises(ValueError):
        q.compare_receipt(recorded, actual)


def test_duplicate_receipt_keys_are_not_accepted():
    with pytest.raises(ValueError):
        q.strict_json.loads('{"count":770,"count":1}')


def test_canonical_receipt_export_disallows_nan():
    with pytest.raises(ValueError):
        q.canonical({'diagnostic': float('nan')})

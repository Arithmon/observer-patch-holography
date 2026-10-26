"""Complete source-history, record, reference and custody controls."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import math
from unittest.mock import patch
import zipfile

import numpy as np
import pytest
import sympy as s

from . import native_field_history as h


@pytest.fixture(scope='module')
def execution():
    return h.execute()


def test_exact_uniform_native_onsite_and_exchange_commute():
    p = h.field.logical_parameters()
    n = s.diag(0, 1, 1, 2)
    onsite = 2*p['onsite_baselines'][0]*s.eye(4)+p['positive_occupation_gaps'][0]*n
    exchange = p['exchange_strength']*(s.eye(4)-h.field.carrier_swap(0, 1, 2))
    assert onsite*exchange == exchange*onsite
    assert p['exchange_strength'] == s.Rational(125, 1024)


def test_composed_native_record_law_against_independent_trigonometry(execution):
    report, arrays = execution
    j = 125/1024
    a, b = math.sin(3*j)**2, math.sin(5*j)**2
    expected = np.array([[a*(1-b), a*b], [(1-a)*b, (1-a)*(1-b)]])
    np.testing.assert_allclose(report['execution']['detector_joint_probabilities'], expected,
                               atol=h.SOURCE_ADMISSION, rtol=0)
    assert report['execution']['terminal_records'] == 319
    assert report['execution']['complete_input_columns'] == 4
    assert report['execution']['all_record_dilation_Frobenius_error_upper_diagnostic'] < h.SOURCE_ADMISSION
    assert report['execution']['first_read_backaction_final_marginal_L1_control'] > .6
    np.testing.assert_array_equal(arrays['labels'], h.record_menu())
    assert set(arrays['labels'][:, 0]) == {0, 1, 2, 3}


def test_number_history_does_not_identify_uniform_native_onsite_gap(execution):
    report, _ = execution
    control = h.number_identifiability_control()
    assert control['fixed_charge'] == 1
    assert control['number_law_depends_only_on_J_and_declared_intervention_parameters']
    assert control['comparison_uniform_onsite_coefficients'] == ['0', '1/256', '1']
    assert control['complete_four_input_dilations_still_contain_relative_charge_phases']
    assert not control['reported_number_history_identifies_native_mass_gap']
    assert not report['scope']['reported_number_history_identifies_native_mass_gap']


def test_native_rows_abort_records_and_failed_ancilla_are_retained(execution):
    report, arrays = execution
    assert arrays['A'].shape[1:] == (512, 4)
    assert arrays['B'].shape[1:] == (512, 4, 4)
    assert len(arrays['exchange_output_A']) == 6
    labels = arrays['labels']
    assert np.sum(labels[:, 0] == 1) == 3
    assert np.sum(labels[:, 0] == 2) == 4
    assert np.sum(labels[:, 0] == 3) == 72
    assert report['execution']['retained_halt_probability_diagnostic'] > 0
    assert report['execution']['all512_row_padding_probability_diagnostic'] > 0
    assert not report['scope']['complete_native_hardware_protection_load_export_history_executed']


def test_four_input_coherence_gate_detects_probability_invisible_phase(execution):
    _, arrays = execution
    changed = deepcopy(arrays)
    # A phase on one input column preserves every individual input's record
    # probabilities, but changes its coherence with an untouched reference.
    changed['B'] *= np.array([1., 1., 1., 1j])[None, None, None, :]
    with pytest.raises(ValueError, match='input/reference'):
        h.metrics(changed)


def test_tiny_halt_branch_is_not_conditionally_normalised(execution):
    _, arrays = execution
    changed = deepcopy(arrays)
    a, b = changed['terminal'][0]
    norm = math.sqrt(float(np.trace(h._gram(changed['A'][a], changed['B'][b])).real))
    assert 0 < norm < 1e-8
    changed['B'][b] /= norm
    with pytest.raises(ValueError, match='input/reference'):
        h.metrics(changed)


def test_factor_grammar_preserves_native_padding_and_little_endian():
    a = np.zeros((512, 4), complex)
    b = np.zeros((512, 4, 4), complex)
    a[511, 0] = 1j
    b[17, 0, 2] = .25
    value = h._dense(a, b)
    assert value[511+512*17, 2] == .25j
    assert np.count_nonzero(value) == 1


def test_compression_keeps_unnormalised_tiny_failure_and_refuses_large_discard():
    value = np.zeros((1 << 18, 4), complex)
    value[511+512*17, 2] = 1e-18j
    a, b, discarded = h._compress(value)
    np.testing.assert_allclose(h._dense(a, b), value, atol=1e-32, rtol=0)
    assert 0 < np.linalg.norm(b) < 1e-17
    assert discarded < 1e-32
    for i in range(5):
        value[i+512*(i+1), 0] = 1.
    with pytest.raises(ArithmeticError, match='compression'):
        h._compress(value)


def test_compression_diagnostics_recomputed_from_own_binary_arrays(execution):
    _, arrays = execution
    h.verify_compressions(arrays)
    changed = deepcopy(arrays)
    changed['compression_residuals'][:] = 1e-9
    with pytest.raises(ValueError, match='own binary'):
        h.verify_compressions(changed)


def test_local_macro_gate_is_all_input_not_one_preparation(execution):
    _, arrays = execution
    changed = deepcopy(arrays)
    changed['fredkin_macro'][7, 7] += 1e-5
    with pytest.raises(ValueError, match='local source macro'):
        h.verify_compressions(changed)


def test_gauge_of_compressed_tensor_factors_does_not_define_replay(execution):
    report, arrays = execution
    changed = deepcopy(arrays)
    unitary = np.diag([1j, -1., -1j, 1.])[:, [2, 0, 3, 1]]
    changed['A'] = changed['A']@unitary
    changed['B'] = np.einsum('rs,nbsd->nbrd', unitary.conj().T, changed['B'])
    own_report = deepcopy(report)
    own_report['execution'] = h.metrics(changed)
    h.compare_execution(own_report, changed, report, arrays)


@pytest.mark.parametrize('mutation', ['missing_record', 'invalid_label', 'masked', 'nonfinite', 'extra_array'])
def test_complete_record_and_original_array_domain(execution, mutation):
    _, arrays = execution
    changed = deepcopy(arrays)
    if mutation == 'missing_record':
        changed['terminal'] = changed['terminal'][:-1]
        changed['labels'] = changed['labels'][:-1]
    elif mutation == 'invalid_label':
        changed['labels'][0, 0] = 0
    elif mutation == 'masked':
        changed['A'] = np.ma.array(changed['A'])
    elif mutation == 'nonfinite':
        changed['B'][0, 0, 0, 0] = np.nan
    else:
        changed['unused'] = np.zeros(1)
    with pytest.raises(ValueError):
        h.metrics(changed)


def test_recovery_is_actual_native_continuation_not_an_encoded_reset():
    encoded = h.code.code_matrix()[:, :2]
    branches = h.code.repair_branches('C', h.SEAMS['C'], encoded)
    bare = h.code.repair_branches('C', h.SEAMS['C'], encoded, recover=False)
    difference = sum(np.linalg.norm(a['amplitude']-b['amplitude'], 'fro')**2
                     for a, b in zip(branches, bare, strict=True))
    assert difference > .1
    for row in branches:
        if row['abort']:
            assert not row['recovery_applied']
        else:
            target = h.code.ideal_code_matrix()[:, :2]*math.sqrt(float(row['code_branch_weight']))
            assert np.linalg.norm(row['amplitude']-target, 2) < h.SOURCE_ADMISSION


def test_fixed_admission_and_all_input_output_phase_tampering(execution):
    report, arrays = execution
    changed = deepcopy(report)
    changed['admissions']['all_record_source_input_reference_error'] *= 2
    with pytest.raises(ValueError, match='numerical value'):
        h.compare_execution(changed, arrays, report, arrays)
    changed = deepcopy(arrays)
    row = max((i for i, label in enumerate(changed['labels']) if label[0] == 0),
              key=lambda i: float(h._gram(changed['A'][changed['terminal'][i, 0]],
                                         changed['B'][changed['terminal'][i, 1]])[1, 1].real))
    changed['A'][changed['terminal'][row, 0], 511, 0] += 1e-3
    with pytest.raises(ValueError):
        h.metrics(changed)
    changed = deepcopy(report)
    changed['execution']['all_record_dilation_Frobenius_error_upper_diagnostic'] = 0.
    with pytest.raises(ValueError, match='own binary'):
        h.compare_execution(changed, arrays, report, arrays)


def test_worst_path_resource_accounting_and_supplied_service_boundary(execution):
    report, _ = execution
    r = report['resources']
    total = (r['preparation_native_pulses']+r['onsite_native_pulses']+r['exchange_native_pulses']
        +r['detector_native_pulses']+sum(r['repair_and_record_selected_recovery_native_pulses_upper'].values()))
    assert r['one_history_native_pulses_upper'] == total
    assert r['required_proper_code_load_export_service_events_upper'] == 2*total
    assert r['maximum_retained_classical_record_bits'] == 11
    assert r['fixed_zero_padding_bits_in_logical_input_registers']+r['reported_seed_logical_input_blanks']+r['fresh_blank_code_qubits'] == 18
    assert not r['proper_load_export_services_executed_in_this_producer']


def test_complete_cli_roundtrip_uses_retained_arrays_and_refuses_duplicate_npz(execution, tmp_path):
    report, arrays = execution
    path = tmp_path/'result.json'
    artifact = path.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    retained = deepcopy(report)
    retained['source_bindings'] = h.source_bindings()
    retained['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    path.write_text(json.dumps(retained))
    with patch.object(h, 'execute', return_value=execution):
        h.check(path)
    with zipfile.ZipFile(artifact, 'a') as archive:
        # A distinct filename spelling can still decode to a duplicate npz key.
        archive.writestr('labels', b'unused')
    retained['artifact']['sha256'] = hashlib.sha256(artifact.read_bytes()).hexdigest()
    path.write_text(json.dumps(retained))
    with pytest.raises(ValueError, match='duplicate'):
        h.check(path)

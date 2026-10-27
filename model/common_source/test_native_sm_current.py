"""Chosen native current law, same-object records and fail-closed evidence."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from unittest.mock import patch

import numpy as np
import pytest
import sympy as s

from . import native_sm_current as n


@pytest.fixture(scope='module')
def exact():
    return n._exact_certificate()


@pytest.fixture(scope='module')
def execution():
    return n.execute()


def test_actual_native_current_matches_fixture_only_after_checked_port_chart():
    data = n.exact_architecture()
    assert data['native_to_fixture_port'] == (6, 4, 7, 5, 2, 0, 3, 1, 10, 8, 11, 9)
    assert data['fixture_intertwiner'] == s.eye(6)
    assert n._zero(data['geometric_frame'].T*data['geometric_frame']-s.eye(6))
    assert n._zero(sum(data['generators'], s.zeros(6))-data['central'])
    assert n.admit_generator(data['central']) == (1,)*12


def test_fixed_grammar_has_complete_visible_twelve_current_response(exact):
    assert exact['real_current_dimension'] == 12
    assert len(exact['commutator_closure']) == 66
    assert exact['complete_EndM6_tomography_rank'] == 12
    assert exact['block_diagonal_tomography_rank'] == 11
    assert exact['native_port_covariance_identities'] == 720
    assert exact['complex_triplet_is_higher_near_class_eigenvalue']
    assert exact['all_60_proper_actions_use_at_most_two_cayley_factors']
    assert len(exact['actual_order_three_cayley_factors']) == 20
    assert exact['grammar_is_selected_architecture_not_repair_selection_of_native17']
    assert not exact['interface_restriction_proves_physical_inaccessibility_or_gauge_redundancy']


@pytest.mark.parametrize('kind', ['second_phase', 'global_phase', 'mixing', 'nonunitary', 'wrong_dimension', 'float'])
def test_public_tangent_rejects_controls_outside_chosen_fixed_law(kind):
    value = s.zeros(6)
    if kind == 'second_phase':
        value[3, 3] = s.I
    elif kind == 'global_phase':
        value = s.I*s.eye(6)
    elif kind == 'mixing':
        value[0, 3] = 1; value[3, 0] = -1
    elif kind == 'nonunitary':
        value[0, 0] = 1
    elif kind == 'wrong_dimension':
        value = s.zeros(8)
    else:
        value[0, 0] = s.Float(.5)*s.I
    with pytest.raises(ValueError):
        n.admit_generator(value)


@pytest.mark.parametrize('operation', [
    {'kind': 'port', 'port': True, 'parameter': '1'},
    {'kind': 'port', 'port': 12, 'parameter': '1'},
    {'kind': 'port', 'port': 0, 'parameter': '1.0'},
    {'kind': 'port', 'port': 0, 'parameter': 1},
    {'kind': 'uniform', 'port': 0, 'parameter': '1'},
    {'kind': 'exchange17', 'port': 0, 'parameter': '1'},
    {'kind': 'uniform', 'port': None, 'parameter': '1', 'apparatus_tape': []},
])
def test_fixed_repair_api_rejects_unregistered_or_coerced_operations(operation):
    with pytest.raises(ValueError):
        n.request_generator(operation)


def test_every_admitted_primitive_has_exact_commuting_factorization():
    data = n.exact_architecture()
    for generator in (*data['generators'], data['central']):
        assert n._zero(s.re(generator)*s.im(generator)-s.im(generator)*s.re(generator))


def test_source_record_history_retains_all_domains_and_observes_center(execution):
    report, arrays = execution
    result = report['execution']
    assert result['unknown_geometric_input_columns'] == 6
    assert result['source_current_output_rows'] == 64
    assert result['local_native_rows'] == 512
    assert result['local_decoded_input_columns'] == 8
    assert result['retained_high_decoder_syndromes_per_carrier'] == 64
    assert result['joint_raw_native_records'] == 128
    assert result['complete_record_pointer_outcomes'] == 256
    assert arrays['controlled'].shape == (4, 64, 6)
    assert arrays['primitive_source'].shape == (13, 8, 8)
    assert result['source_executed_port_primitives'] == 12
    assert max(result['complete_primitive_menu_dilation_errors_diagnostic']) < n.SOURCE_ADMISSION
    assert arrays['record_effects'].shape == (4, 8, 16, 2, 6, 6)
    assert result['current_dilation_error_diagnostic'] < n.SOURCE_ADMISSION
    assert result['cross_triplet_marker_probabilities_diagnostic'][0] > 1-12*n.SOURCE_ADMISSION
    assert abs(result['cross_triplet_marker_probabilities_diagnostic'][1]-np.cos(.25)**2) < 12*n.SOURCE_ADMISSION
    assert result['ordered_response_difference_diagnostic'] > 1e-5
    assert not result['any_native_raw_record_or_decoder_syndrome_postselected']
    assert not result['full64_preparation_input_basis_source_executed']


def test_record_effect_tensor_contraction_matches_literal_selected_syndromes(execution):
    _, arrays = execution
    A = arrays['decoded_C'][0].reshape(64, 8, 8)
    B = arrays['decoded_E'][0].reshape(64, 8, 8)
    detector = arrays['detector'][0]
    # Independent small subdomain identity; the production calculation sums
    # all64 syndromes, whose complete factors are retained in the artifact.
    a, b = A[:2], B[:3]
    direct = sum((np.kron(x, y).conj().T@detector.conj().T@detector@np.kron(x, y)
                  for x in a for y in b), np.zeros((64, 64), complex))
    ca = np.einsum('hia,hjc->ijac', a.conj(), a)
    cb = np.einsum('hkb,hld->klbd', b.conj(), b)
    effect = (detector.conj().T@detector).reshape(8, 8, 8, 8)
    pulled = np.einsum('ijac,klbd,ikjl->abcd', ca, cb, effect, optimize=True).reshape(64, 64)
    np.testing.assert_allclose(pulled, direct, atol=2e-14, rtol=0)


def test_record_measurement_retains_unknown_input_reference_coherences(execution):
    _, arrays = execution
    for row in arrays['record_effects']:
        np.testing.assert_allclose(row.sum(axis=(0, 1, 2)), np.eye(6), atol=12*n.SOURCE_ADMISSION, rtol=0)
        for effect in row.reshape(-1, 6, 6):
            np.testing.assert_allclose(effect, effect.conj().T, atol=1e-12, rtol=0)
            assert np.linalg.eigvalsh(effect).min() > -1e-12


@pytest.mark.parametrize('kind', ['missing', 'extra', 'nonfinite', 'dtype', 'record', 'phase', 'syndrome', 'effect'])
def test_complete_source_evidence_rejects_tampering(execution, kind):
    report, arrays = execution
    changed = deepcopy(arrays)
    if kind == 'missing':
        changed.pop('detector')
    elif kind == 'extra':
        changed['hidden'] = np.zeros(1)
    elif kind == 'nonfinite':
        changed['decoded_C'][0, 511, 7] = np.nan
    elif kind == 'dtype':
        changed['controlled'] = changed['controlled'].astype(np.complex64)
    elif kind == 'record':
        changed['raw_labels'][0, 0] = 9
    elif kind == 'phase':
        changed['controlled'][:, :, 0] *= 1j
    elif kind == 'syndrome':
        changed['decoded_C'][0, 511, 7] = .1
    else:
        changed['record_effects'][0, 0, 0, 0, 0, 0] += .1
    with pytest.raises(ValueError):
        n.compare_execution(report, changed, report, arrays)


def test_universal_apparatus_is_not_misrepresented_as_passive_repair_derivation(execution):
    report = execution[0]
    assert not report['architecture']['old_seventeen_current_exchange_is_part_of_this_revised_law']
    assert not report['architecture']['passive_CE_derives_the_apparatus_current_primitives']
    assert not report['architecture']['all_physically_possible_controller_interventions_declared_inaccessible']
    assert not report['scope']['full_model_A1_complete_accessible_response_discharged']
    assert report['scope']['complete_current_grammar_is_a_declared_operational_architecture']
    assert not report['scope']['physical_gauge_redundancy_or_particle_masses_derived']


def test_admission_and_exact_boolean_types_are_not_tolerant_claims(execution):
    report, arrays = execution
    changed = deepcopy(report); changed['admissions']['source'] += 1e-12
    with pytest.raises(ValueError):
        n.compare_execution(changed, arrays, report, arrays)
    changed = deepcopy(report); changed['scope']['measured_particle_parameters_used'] = 0
    with pytest.raises(ValueError, match='type'):
        n.compare_execution(changed, arrays, report, arrays)


def test_rehashed_receipt_replay_and_source_drift(execution, tmp_path):
    report, arrays = execution
    path = tmp_path/'current.json'; artifact = path.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    pins = n.source_bindings()
    saved = {**report, 'source_bindings': pins,
             'artifact': {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}}
    path.write_text(json.dumps(saved))
    with patch.object(n, 'execute', return_value=execution):
        n.check(path)
    with patch.object(n, 'execute', return_value=execution), patch.object(n, 'source_bindings',
            side_effect=[pins, {**pins, 'modified': 'changed'}]):
        with pytest.raises(ValueError, match='changed during replay'):
            n.check(path)

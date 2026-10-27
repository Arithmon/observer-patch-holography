"""Joint source charge/gap readout, complements and blind inverse contracts."""
from copy import deepcopy
import math
from pathlib import PureWindowsPath
from unittest.mock import patch

import numpy as np
import pytest

from . import native_particle_records as h


@pytest.fixture(scope='module')
def execution():
    return h.execute()


def test_fixed_pair_frame_comes_from_internal_higgs_geometry():
    a = h.field_architecture()
    assert len(a['pairs']) == 14
    assert sorted(a['images']) == list(range(64))
    np.testing.assert_array_equal(a['signed_pair_frame'].T@a['signed_pair_frame'], np.eye(64))
    for left, right, _ in a['pairs']:
        assert a['charge'][left] == a['charge'][right]
        assert a['color'][left] == a['color'][right]
        assert (left < 32) != (right < 32)
    assert h.calibration._zero(a['P15'].T*a['P15']-h.s.eye(15))


def test_source_register_keeps_coherent_families_and_complements(execution):
    report, arrays = execution
    assert report['execution']['complete_register_input_columns'] == 256
    assert report['execution']['terminal_pointer_outcomes'] == 256
    assert arrays['phases'].shape == (4, 256, 256)
    assert arrays['probabilities'].shape == (42, 4, 2, 256)
    for i in range(4):
        unitary = h._evolution(arrays, i)
        np.testing.assert_allclose(unitary.conj().T@unitary, np.eye(256), atol=h.SOURCE_ADMISSION, rtol=0)
        np.testing.assert_allclose(unitary[192:, 192:], np.eye(64), atol=h.SOURCE_ADMISSION, rtol=0)
    assert not report['calibration']['masses_predicted_without_empirical_parameters']
    assert not report['calibration']['experimental_or_interacting_pole_masses_recovered']


def test_blind_inverse_recovers_charge_and_running_gap_without_names(execution):
    _, arrays = execution
    dictionary = h.record_dictionary()
    with patch.object(h, 'field_architecture', side_effect=AssertionError('inverse accessed action')):
        rows = [h.decode_records(p, dictionary=dictionary) for p in arrays['probabilities']]
    assert all(row['resolved'] and row['identified_particle_name_or_family_label'] is None for row in rows)
    assert {row['electric_Q'] for row in rows} == {'-1', '1', '-2/3', '2/3', '-1/3', '1/3'}
    assert {row['color_C2'] for row in rows} == {'0', '4/3'}


def test_full_pointer_factor_matches_literal_complex_kraus_contraction(execution):
    _, arrays = execution
    rng = np.random.default_rng(71)
    state = rng.normal(size=(256, 2))+1j*rng.normal(size=(256, 2))
    state /= np.linalg.norm(state, axis=0)
    result = h.pointer_probabilities(state, arrays['pointer'])
    for record in (0, 71, 255):
        K = np.ones((1, 1), complex)
        for bit in reversed(range(8)):
            K = np.kron(K, arrays['pointer'][(record >> bit) & 1])
        np.testing.assert_allclose(result[record], np.sum(abs(K@state)**2, axis=0), atol=2e-15, rtol=0)
    np.testing.assert_allclose(result.sum(axis=0), np.ones(2), atol=1e-12)


def test_inverse_rejects_wrong_heldout_gap(execution):
    _, arrays = execution
    values = arrays['probabilities'][0].copy()
    values[3] = arrays['probabilities'][14, 3]
    result = h.decode_records(values, dictionary=h.record_dictionary())
    assert not result['resolved']
    assert 'held-out' in result['reason']


def test_inverse_has_no_hidden_positive_gap_or_charge_postselection(execution):
    _, arrays = execution
    mixed = (arrays['probabilities'][0]+arrays['probabilities'][13])/2
    assert not h.decode_records(mixed, dictionary=h.record_dictionary())['resolved']
    with pytest.raises(ValueError):
        h.decode_records(arrays['probabilities'][0]*.9, dictionary=h.record_dictionary())
    with pytest.raises(ValueError):
        h.decode_records(arrays['probabilities'][0], dictionary=h.record_dictionary(), mass_window=(0, 10000))


@pytest.mark.parametrize('change', ['missing', 'mask', 'nan', 'probe', 'extra'])
def test_complete_source_array_schema_fails_closed(execution, change):
    _, arrays = execution
    changed = deepcopy(arrays)
    if change == 'missing': changed.pop('macro')
    elif change == 'mask': changed['phases'] = np.ma.array(changed['phases'])
    elif change == 'nan': changed['phases'][0, 255, 255] = np.nan
    elif change == 'probe': changed['probe_indices'][0] = 1
    else: changed['untracked'] = np.zeros(1)
    with pytest.raises(ValueError): h._schema(changed)


def test_contract_bounds_are_exact_and_include_all_services(execution):
    report, _ = execution
    r = report['resources']
    assert r['required_proper_load_export_events_upper'] == 2*r['one_protocol_native_pulses_upper']
    assert r['frame_native_pulses'] > 0 and r['pointer_native_pulses'] == 24
    assert report['execution']['fixed_channel_allowance'] < h.RECORD_ERROR
    changed = deepcopy(report)
    changed['execution']['fixed_channel_allowance'] += 1e-12
    with pytest.raises(ValueError): h._compare(report, changed)


@pytest.mark.parametrize('change', ['padding_phase', 'reference_phase', 'macro'])
def test_unprepared_columns_and_source_helpers_are_checked(execution, change):
    _, arrays = execution
    changed = deepcopy(arrays)
    if change == 'padding_phase': changed['phases'][0, 255, 255] *= 1j
    elif change == 'reference_phase': changed['frame'][:, 63] *= 1j
    else: changed['macro'][:, 7] *= -1
    with pytest.raises(ValueError): h.metrics(changed)


def test_resource_accounting_uses_the_same_actual_current_backend(execution):
    with patch.object(h.compiler, '_native', side_effect=AssertionError('old fixture backend')):
        resources = h._source_resources()
    assert resources['one_protocol_native_pulses_upper'] == execution[0]['resources']['one_protocol_native_pulses_upper']


def test_repo_relative_calibration_path_is_platform_independent():
    packet = PureWindowsPath('C:/checkout/model/common_source/calibrated_matter_parameters.json')
    root = PureWindowsPath('C:/checkout')
    with patch.object(h.calibration, 'PACKET', packet), patch.object(h, 'ROOT', root):
        assert h._calibration_packet_name() == 'model/common_source/calibrated_matter_parameters.json'


def test_fresh_source_diagnostic_is_not_the_retained_matrix_diagnostic(execution):
    report, arrays = execution
    fresh = deepcopy(report)
    fresh['execution']['complete_protocol_operator_error_diagnostic'] += 1e-9
    h.compare_execution(report, arrays, fresh, arrays)


def test_retained_diagnostic_cannot_be_normalized_away(execution):
    report, arrays = execution
    changed = deepcopy(report)
    changed['execution']['complete_protocol_operator_error_diagnostic'] = 0.
    with pytest.raises(ValueError, match='own binary'):
        h.compare_execution(changed, arrays, report, arrays)


@pytest.mark.parametrize('kind', ['replay_distance', 'nonfinite', 'extra'])
def test_fresh_outputs_still_obey_complete_replay_admission(execution, kind):
    report, arrays = execution
    changed = deepcopy(arrays)
    if kind == 'replay_distance': changed['frame'][0, 0] += 2*h.REPLAY_ADMISSION
    elif kind == 'nonfinite': changed['frame'][0, 0] = np.nan
    else: changed['extra'] = np.zeros(1)
    with pytest.raises(ValueError): h.compare_execution(report, arrays, report, changed)

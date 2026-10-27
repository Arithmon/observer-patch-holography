"""Native invariance, local resolvents, current algebra and receipt custody."""
from copy import deepcopy
import hashlib
import io
import json
import zipfile

import numpy as np
import pytest
import sympy as s

from . import native_persistent_response as m


@pytest.fixture(scope='module')
def built():
    return m.build()


def test_native_encoder_and_local_corrected_code_projectors():
    u, columns = m.native_port_columns(0, 1)
    pa, pb = m.local_code_projectors(0, 1)
    assert u.shape == (10626, 36)
    assert np.array_equal((u.T@u).toarray(), 36*np.eye(36, dtype=int))
    assert len(columns) == 12
    assert all((u.T@column).nnz == 0 for column in columns)
    assert (pa@pa-6*pa).nnz == (pb@pb-6*pb).nnz == 0
    assert (pa@pb-pb@pa).nnz == 0
    assert (pa@u-6*u).nnz == (pb@u-6*u).nnz == 0


def test_weighted_local_virtual_inverse_on_all_native_outputs():
    _, columns = m.native_port_columns(0, 1)
    pa, pb = m.local_code_projectors(0, 1)
    for column in columns:
        weighted = 6*column+pa@column+pb@column
        # H0=2I-(Pa_integer+Pb_integer)/6 and R0A=weighted/12.
        assert (12*weighted-pa@weighted-pb@weighted-72*column).nnz == 0
        assert (pa@column).nnz and (pb@column).nnz
    total = sum(columns)
    assert (pa@total).nnz == (pb@total).nnz == 0
    flat = m.response_matrices(0, 1)[0]
    local = m.response_matrices(0, 1, 'additive_local_penalty')[0]
    assert any(a != 2*b for a, b in zip(flat, local, strict=True))
    weights = [1 if p == q else 2 for p, q in m.PORT_PAIRS]
    assert sum((w*h for w, h in zip(weights, local, strict=True)), s.zeros(36)) == m.native_interaction(0, 1)/2


@pytest.mark.parametrize('law', ['flat_joint_phase', 'additive_local_penalty'])
def test_preserving_coefficients_act_on_the_complete_36_color_space(law):
    frame = m.geometric_coordinates()['frame']
    full, _, constraints = m.response_matrices(0, 1, law)
    preserved = m.preserving_responses(0, 1, law)
    assert preserved['constraint_rank'] == 58
    assert preserved['coefficient_kernel_dimension'] == 20
    assert preserved['response_dimension'] == 6
    for coefficients, reduced in zip(preserved['source_coefficients'], preserved['basis'], strict=True):
        assert constraints*coefficients.T == s.zeros(216, 1)
        operator = sum((a*b for a, b in zip(coefficients, full, strict=True)), s.zeros(36))
        assert operator*frame == frame*reduced


def test_additive_and_joint_phase_routes_have_the_same_preserving_response_space():
    flat = m.preserving_responses(0, 1)['basis']
    local = m.preserving_responses(0, 1, 'additive_local_penalty')['basis']
    assert m._rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in flat+local))) == 6


@pytest.mark.parametrize('key', ['same_color_response', 'additive_local_same_color_response'])
def test_equal_coloring_response_is_only_scalar_after_native_leakage_constraints(built, key):
    report = built[0][key]
    assert report['preserving_Hermitian_response_dimension'] == 1
    assert report['Lie_real_dimension'] == 1
    assert report['commutator_observable_dimension'] == 0


@pytest.mark.parametrize('key', ['opposite_color_response', 'additive_local_opposite_color_response'])
def test_actual_opposite_coloring_has_17_currents_and_no_relative_triplet_center(built, key):
    report = built[0][key]
    assert (report['Lie_real_skew_dimension'], report['Lie_imaginary_selfadjoint_dimension']) == (6, 11)
    assert report['Lie_real_dimension'] == 17
    assert report['commutator_observable_dimension'] == 16
    assert report['geometric_triplet_phase_traceless_ranks'] == [5, 5]
    assert not report['preserving_phase_response_couples_distinct_triplets']
    coordinates = m.geometric_coordinates()
    three = coordinates['left_inverse']*m.objects.geometric_projectors()['three']*coordinates['frame']
    _, phase = m.lie_closure(m.preserving_responses(0, 1)['basis'], coordinates['metric'])
    assert all(s.expand(s.trace(three*h)-s.trace((s.eye(6)-three)*h)) == 0 for h in phase)
    assert m._rank(s.Matrix.hstack(*(h.reshape(36, 1) for h in phase), three.reshape(36, 1))) == 12


def test_enlarged_interface_visibility_is_not_confused_with_end_m6(built):
    report = built[0]['full_interface_scope']
    assert report['full_preserving_Hermitian_response_rank'] == 9
    assert report['restricted_geometric_M6_response_rank'] == 6
    assert report['full_current_Lie_dimension_lower_bound'] == 20
    assert report['full_End36_commutator_dimension_lower_bound'] == 19
    assert report['global_identity36_belongs_to_preserving_response_span']
    assert report['uniform_B_noncommuting_selected_preserving_lifts'] == 5
    assert not report['uniform_B_is_an_enlarged_interface_central_generator']
    assert report['complement_coherence_can_observe_uniform_B_relative_phase']


def test_actual_finite_rotations_are_endogenous_only_in_the_stated_ideal_current_group(built):
    report = built[0]['native_rotation_endogeneity']
    assert report['actual_native_rotation_count'] == 60
    assert report['order_three_rotation_count'] == 20
    assert report['every_rotation_is_at_most_two_order_three_rotations']
    assert report['order_three_Cayley_generators_belong_to_actual_real_skew_response_span']
    assert not report['relative_class_sum_phase_in_original_current_span']
    assert report['current_dimension_after_promoting_class_sum_to_phase_control'] == 18
    assert not report['finite_source_pulse_and_physical_error_certificate_supplied_by_this_identity']


def test_actual_primary_four_carrier_interaction_retains_native_leakage(built):
    report = built[0]['four_carrier_composite_interaction']
    assert report['matching_color_pattern'] == [0, 1, 0, 1]
    assert report['full_four_carrier_color_dimension'] == 1296
    assert report['leakage_norm_squared_exact'] == '5626523/1620000'
    assert not report['full_action_preserves_geometric_object_pair']
    bell = built[0]['actual_pair_reference_test']
    assert bell['geometric_Bell_block_exact'] == [['92/5', '-2/45'], ['-2/45', '92/5']]
    assert bell['geometric_Bell_eigenvalues_exact'] == ['826/45', '166/9']
    assert bell['upper_Bell_is_unique_projected_maximum']
    assert bell['projected_Gram_gap_exact'] == '1/45'
    assert bell['next_projected_Gram_eigenvalue_exact'] == '829/45'
    mixed = bell['actual_mixed_reference_ray_stabilizer']
    assert mixed['raw_current_dimension'] == 7
    assert mixed['End_M6_conjugation_observable_dimension'] == 6
    supplied = bell['supplied_pure_triplet_reference_control_stabilizer']
    assert supplied['raw_current_dimension'] == 12
    assert supplied['End_M6_conjugation_observable_dimension'] == 11
    assert not bell['pure_triplet_Bell_lines_are_eigenrays']
    assert not bell['full_native_four_carrier_ground_state_selected_by_this_compression']


@pytest.mark.parametrize('key', ['same_color_global_collective_stabilizer', 'projected_composite_global_collective_stabilizer'])
def test_actual_interaction_has_only_scalar_connected_collective_symmetry(built, key):
    report = built[0][key]
    assert report['real_skew']['dimension'] == 0
    assert report['imaginary_symmetric']['dimension'] == 1
    assert s.Matrix(report['imaginary_symmetric']['basis_exact'][0]) == s.eye(6)
    assert not report['local_gauge_redundancy_derived']


def test_declared_operational_and_physical_boundaries(built):
    report = built[0]
    assert report['empirical_inputs_to_prediction'] == []
    assert not report['scope']['desired_gauge_current_matrices_imported']
    assert not report['scope']['Standard_Model_gauge_current_law_selected']
    assert not report['scope']['measured_particle_spectrum_derived']
    control = report['operational_control_scope']
    assert control['negative_Gram_term_requires_inverse_commutator_word']
    assert not control['all_selected_signed_schedules_executed_as_complete_native_histories_here']
    assert not control['exact_autonomous_finite_current_law_selected_by_CE_repair']


@pytest.mark.parametrize('colors', [(2, 0), (0, -1), (True, 0)])
def test_invalid_native_color_choices_rejected(colors):
    with pytest.raises(ValueError):
        m._colors(*colors)


def test_unknown_protection_law_rejected():
    with pytest.raises(ValueError):
        m.response_matrices(0, 1, 'uncorrected_local_Gram')


def _receipt(tmp_path, report, arrays):
    path = tmp_path/'receipt.json'
    artifact = path.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    value = deepcopy(report)
    value['source_bindings'] = m.source_bindings()
    value['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    path.write_text(json.dumps(value, allow_nan=False))
    return path


@pytest.mark.parametrize('tamper', ['missing', 'extra', 'changed', 'nonfinite', 'dtype', 'shape', 'duplicate'])
def test_artifact_tampering_rejected_even_with_rehashed_binary(tmp_path, monkeypatch, tamper):
    report, arrays = {'count': 17}, {'data': np.eye(2)}
    monkeypatch.setattr(m, 'build', lambda: (report, arrays))
    corrupted = {key: value.copy() for key, value in arrays.items()}
    if tamper == 'missing':
        corrupted = {}
    elif tamper == 'extra':
        corrupted['extra'] = np.eye(2)
    elif tamper == 'changed':
        corrupted['data'][0, 1] = 1e-16
    elif tamper == 'nonfinite':
        corrupted['data'][0, 0] = np.nan
    elif tamper == 'dtype':
        corrupted['data'] = corrupted['data'].astype(np.float32)
    elif tamper == 'shape':
        corrupted['data'] = corrupted['data'].reshape(4)
    path = _receipt(tmp_path, report, corrupted)
    if tamper == 'duplicate':
        buffer = io.BytesIO()
        np.save(buffer, arrays['data'])
        with pytest.warns(UserWarning, match='Duplicate name'):
            with zipfile.ZipFile(path.with_suffix('.npz'), 'a') as archive:
                archive.writestr('data.npy', buffer.getvalue())
        value = json.loads(path.read_text())
        value['artifact']['sha256'] = hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest()
        path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        m.check(path)


@pytest.mark.parametrize('corrupted', [{'count': True}, {'count': 17, 'extra': False}, {'count': 18}])
def test_json_claim_tampering_and_boolean_integer_alias_rejected(tmp_path, monkeypatch, corrupted):
    monkeypatch.setattr(m, 'build', lambda: ({'count': 17}, {'data': np.eye(2)}))
    path = _receipt(tmp_path, corrupted, {'data': np.eye(2)})
    with pytest.raises(ValueError):
        m.check(path)


def test_source_change_during_replay_rejected(tmp_path, monkeypatch):
    report, arrays = {'count': 17}, {'data': np.eye(2)}
    path = _receipt(tmp_path, report, arrays)
    pins = m.source_bindings()
    calls = iter([pins, {**pins, 'injected': 'different'}])
    monkeypatch.setattr(m, 'source_bindings', lambda: next(calls))
    monkeypatch.setattr(m, 'build', lambda: (report, arrays))
    with pytest.raises(ValueError, match='sources changed during replay'):
        m.check(path)


def test_exact_receipt_replay_and_stale_source_rejection(tmp_path, monkeypatch):
    report, arrays = {'count': 17}, {'data': np.eye(2)}
    monkeypatch.setattr(m, 'build', lambda: (report, arrays))
    path = _receipt(tmp_path, report, arrays)
    m.check(path)
    value = json.loads(path.read_text())
    value['source_bindings']['model/common_source/native_persistent_response.py'] = '0'*64
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='source bindings'):
        m.check(path)

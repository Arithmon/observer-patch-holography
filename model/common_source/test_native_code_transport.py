"""Same native object, actual port action, retained outcomes and references."""
from copy import deepcopy
import hashlib
import json
from unittest.mock import patch

import numpy as np
import pytest

from . import native_code_transport as h


@pytest.fixture(scope='module')
def execution():
    return h.execute()


def test_native_code_is_derived_from_axis_colorings_and_is_one_fixed_object():
    assert h.code.MATCHINGS == h.object_geometry.equivariant_matchings(0)
    assert h.code.MATCHINGS != h.base.MATCHINGS
    target = h.code.ideal_code_matrix()
    np.testing.assert_allclose(target.conj().T@target, np.eye(6), atol=1e-15)
    for p in h.rotations():
        np.testing.assert_allclose(h._ideal_rotated_code(p), target@h.logical_rotation(p), atol=1e-15)
        np.testing.assert_allclose(h.logical_rotation(p)@np.ones(6), np.ones(6), atol=0)


def test_exact_full_carrier_native_record_covariance():
    proof = h.exact_covariance()
    assert proof['basis_record_equalities_checked'] == 561600
    assert proof['native_integer_carrier_states'] == 78
    assert proof['same_fixed_code_preserved_by_all_rotations']
    assert proof['fixed_code_stabilizer_order'] == 60
    assert proof['distinct_code_projectors'] == 1
    assert not proof['continuous_M6_or_SM_action_follows_from_port_covariance']


def test_orientation_transport_complements_coin_and_keeps_raw_abort_meaning():
    p = h.selected_rotation()
    seam = h.SEAMS['E']
    assert p[seam[0]] > p[seam[1]]
    for law in ('C', 'E'):
        for coin in (0, 1):
            for raw in range(1 << h.circuit.RECORD_BITS[law]):
                _, moved_coin, moved_raw = h.transported_record(law, seam, coin, raw, p)
                assert moved_coin == 1-coin
                assert (moved_raw & 4) == (raw & 4)
                assert (raw in h.circuit.record_values(law, coin)) == (moved_raw in h.circuit.record_values(law, moved_coin))


@pytest.mark.parametrize("local_code", list(h.PAIR_CODES.values()))
def test_all_seam_branch_images_have_scalar_gram_and_disjoint_columns(local_code):
    for law in ('C', 'E'):
        for seam in local_code.SEAMS:
            total = 0
            for coin in (0, 1):
                for raw in range(1 << h.circuit.RECORD_BITS[law]):
                    rows = local_code.branch_images(law, seam, coin, raw)
                    assert len({len(row) for row in rows}) == 1
                    flat = [i for row in rows for i in row]
                    assert len(flat) == len(set(flat))
                    total += local_code.branch_weight(law, seam, coin, raw)
            assert total == 1


def test_actual_source_encoder_and_unknown_input_history(execution):
    report, arrays = execution
    np.testing.assert_allclose(h.code.code_matrix(), h.code.ideal_code_matrix(), atol=2e-10, rtol=0)
    r = report['execution']
    assert r['complete_unknown_input_columns'] == 6
    assert r['terminal_records'] == 158
    assert arrays['outputs'].shape == (158, 512, 6)
    assert set(arrays['labels'][:, 0]) == {0, 1, 2}
    assert 0 < r['raw_halt_Frobenius_probability_diagnostic'] < 1e-20
    assert r['complete_history_dilation_error_diagnostic'] < h.SOURCE_ADMISSION
    assert report['scope']['same_fixed_code_invariant_under_all_proper_rotations']
    assert not report['scope']['arbitrary_code_isometry_used_to_identify_M6_generators']


def test_actual_marker_history_preserves_the_persistent_quintet(execution):
    _, arrays = execution
    q = np.eye(6)-np.ones((6, 6))/6
    code = h.code.ideal_code_matrix()
    incoming = arrays['outputs']@q
    projected = np.einsum('ai,ij,bj->ab', code, q, code.conj())
    residual = incoming-np.einsum('ab,nbc->nac', projected, incoming)
    assert np.linalg.norm(residual) < h.SOURCE_ADMISSION
    assert abs(execution[0]['execution']['logical_marker_wrong_fixed_chart_probability']-.25) < 1e-12


def test_complete_factored_pair_contains_geometric6_and_its_complement(execution):
    report, arrays = execution
    pair = report['two_carrier_geometric_object']
    assert pair['source_matching_choices'] == [0, 1]
    assert h.PAIR_CODES['E'].MATCHINGS == h.object_geometry.equivariant_matchings(1)
    assert not pair['pair_persistence_schedule_includes_candidate_interaction_dynamics']
    assert pair['full_product_logical_input_columns'] == 36
    assert pair['selected_source_geometric_input_columns'] == 6
    assert pair['retained_complementary_product_levels'] == 30
    assert pair['complete_joint_raw_record_outcomes'] == 128
    assert pair['complete_product_Gram_defect_diagnostic'] < 4*h.SOURCE_ADMISSION
    frame = arrays['geometric_frame']
    complement = arrays['product_complement']
    np.testing.assert_allclose(frame@frame.conj().T+complement, np.eye(36), atol=1e-12)
    assert abs(np.trace(complement)-30) < 1e-12
    assert not pair['other_product_sectors_discarded_or_declared_unphysical']
    assert pair['geometric_sector_state_preparation_source_word_executed']
    assert pair['full_product_register_input_columns'] == 64
    assert pair['local_complete_decoded_input_columns'] == 8
    assert pair['retained_padded_product_register_levels'] == 28


def test_actual_simultaneous_port_rotations_preserve_pair_object_but_local_ones_need_not(execution):
    _, arrays = execution
    frame, complement = arrays['geometric_frame'], arrays['product_complement']
    violations = []
    for p in h.rotations():
        r = h.logical_rotation(p)
        assert np.linalg.norm(complement@np.kron(r, r)@frame) < 1e-12
        violations.append(np.linalg.norm(complement@np.kron(r, np.eye(6))@frame))
    assert max(violations) > .1


def test_single_retained_pair_outcome_matches_explicit_small_native_tensor(execution):
    _, arrays = execution
    a, b, frame = arrays['pair_C'][0], arrays['pair_E'][0], arrays['prepared_geometric_inputs']
    # Selected physical output rows suffice to check tensor indexing;
    # the stored factors themselves retain all512 rows on both carriers.
    rows = (0, 1, 17, 511)
    explicit = np.kron(a[list(rows)], b[list(rows)])@frame
    factored = np.einsum('ai,bj,ijc->abc', a[list(rows)], b[list(rows)], frame.reshape(8, 8, 6)).reshape(16, 6)
    np.testing.assert_allclose(explicit, factored, atol=1e-25, rtol=1e-12)


@pytest.mark.parametrize('bad', ['missing', 'mask', 'nonfinite', 'phase', 'halt'])
def test_full_unknown_input_custody_rejects_bad_history(execution, bad):
    _, arrays = execution
    changed = deepcopy(arrays)
    if bad == 'missing':
        changed['outputs'] = changed['outputs'][:-1]
    elif bad == 'mask':
        changed['outputs'] = np.ma.array(changed['outputs'])
    elif bad == 'nonfinite':
        changed['outputs'][0, 0, 0] = np.inf
    elif bad == 'phase':
        changed['outputs'][:, :, 0] *= 1j
    else:
        candidates = np.where(changed['labels'][:, 0] != 0)[0]
        index = max(candidates, key=lambda i: np.linalg.norm(changed['outputs'][i]))
        changed['outputs'][index] /= np.linalg.norm(changed['outputs'][index])
    with pytest.raises(ValueError):
        h.metrics(changed)


@pytest.mark.parametrize('bad', ['input_phase', 'complement', 'record', 'extra', 'preparation_phase', 'preparation_padding', 'extension', 'local_padding'])
def test_pair_custody_keeps_references_complement_and_all_records(execution, bad):
    _, arrays = execution
    changed = deepcopy(arrays)
    if bad == 'input_phase':
        changed['pair_C'][:, :, 0] *= 1j
    elif bad == 'complement':
        changed['product_complement'] *= 0
    elif bad == 'record':
        changed['pair_E_labels'][0, 1] = 7
    elif bad == 'preparation_phase':
        changed['prepared_geometric_inputs'][:, 0] *= 1j
    elif bad == 'preparation_padding':
        changed['prepared_geometric_inputs'][63, 0] = .1
    elif bad == 'extension':
        changed['preparation_orthogonal_extension'][:, 63] *= 1j
    elif bad == 'local_padding':
        changed['pair_C'][:, :, 7] *= 0
    else:
        changed['pair_C'] = changed['pair_C'][:-1]
    with pytest.raises(ValueError):
        h.pair_metrics(changed)


def test_reported_diagnostics_are_bound_to_their_own_binary_arrays(execution):
    report, arrays = execution
    bad = deepcopy(report)
    bad['two_carrier_geometric_object']['tensor_dilation_error_upper_diagnostic'] = 0.
    with pytest.raises(ValueError, match='own binary'):
        h.compare_execution(bad, arrays, report, arrays)


def test_source_service_counts_include_actual_rotations_and_recoveries(execution):
    r = execution[0]['resources']
    for row in r['per_reached_frame'].values():
        assert row['native_pulses_upper'] == (row['C_repair_recovery_upper']+row['E_repair_recovery_upper']
            +row['actual_port_rotation_invocation']+row['encoding_decoding']+row['marker_read'])
    assert r['required_proper_load_export_service_events_upper'] == 2*r['one_history_native_pulses_upper']
    assert r['two_carrier_preparation_persistent_history_proper_services_upper'] == 2*r['two_carrier_preparation_persistent_history_native_pulses_upper']
    prep = r['two_carrier_geometric_preparation']
    assert prep['native_pulses'] > 0
    assert prep['proper_code_service_events'] == 2*prep['native_pulses']
    assert prep['identity_phase_pulses'] > 0
    assert r['geometric_source_preparation_interaction_and_phase_services_counted']
    assert not r['physical_load_export_protection_and_memory_services_executed']


def test_receipt_roundtrip_replays_source_outputs(execution, tmp_path):
    report, arrays = execution
    path = tmp_path/'receipt.json'
    binary = path.with_suffix('.npz')
    np.savez_compressed(binary, **arrays)
    saved = deepcopy(report)
    saved['source_bindings'] = h.source_bindings()
    saved['artifact'] = {'filename': binary.name, 'sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}
    path.write_text(json.dumps(saved))
    with patch.object(h, 'execute', return_value=execution):
        h.check(path)


def test_actual_geometric_preparation_preserves_every_unknown_reference_and_padding(execution):
    report, arrays = execution
    prepared = arrays['prepared_geometric_inputs']
    target = np.zeros((64, 6), complex)
    target[h.product_addresses()] = arrays['geometric_frame']
    assert prepared.shape == (64, 6)
    assert np.linalg.norm(prepared-target) < h.SOURCE_ADMISSION
    assert np.linalg.norm(prepared.conj().T@prepared-np.eye(6)) < h.SOURCE_ADMISSION
    extension = arrays['preparation_orthogonal_extension']
    np.testing.assert_allclose(extension[:, 8*np.arange(6)], target, atol=1e-12)
    np.testing.assert_allclose(extension.conj().T@extension, np.eye(64), atol=1e-12)
    pair = report['two_carrier_geometric_object']
    assert 0 < pair['source_preparation_padded_row_norm_diagnostic'] < h.SOURCE_ADMISSION
    assert pair['preparation_and_tensor_dilation_error_upper_diagnostic'] >= pair['tensor_dilation_error_upper_diagnostic']
    assert pair['preparation_and_tensor_dilation_error_upper_diagnostic'] < pair['fixed_product_dilation_allowance']
    assert pair['complete_orthogonal_extension_columns_source_executed'] == 6
    assert not pair['native_repair_selects_the_geometric_six_subspace']
    assert not pair['preparation_uses_measured_particle_parameters']


def test_complete_local_factors_propagate_padding_instead_of_code_projection(execution):
    _, arrays = execution
    for law in ('C', 'E'):
        maps = arrays['pair_'+law]
        assert maps.shape[-1] == 8
        gram = np.einsum('nac,nad->cd', maps.conj(), maps)
        np.testing.assert_allclose(gram, np.eye(8), atol=h.SOURCE_ADMISSION, rtol=0)
    a, b = (maps[max(range(len(maps)), key=lambda k: np.linalg.norm(maps[k, :, 7]))]
        for maps in (arrays['pair_C'], arrays['pair_E']))
    prepared = arrays['prepared_geometric_inputs'].copy()
    prepared[63, 0] += .1
    rows = np.argsort(np.linalg.norm(a[:, 7:8], axis=1))[-4:]
    columns = np.argsort(np.linalg.norm(b[:, 7:8], axis=1))[-4:]
    actual = np.kron(a[rows], b[columns])@prepared
    projected = prepared.copy()
    projected[np.setdiff1d(np.arange(64), h.product_addresses())] = 0
    assert np.linalg.norm(actual-np.kron(a[rows], b[columns])@projected) > 1e-5


@pytest.mark.parametrize('field', ['fixed_product_dilation_allowance', 'fixed_product_channel_error_allowance'])
def test_publication_bounds_are_exact_contract_values(execution, field):
    report, arrays = execution
    changed = deepcopy(report)
    changed['two_carrier_geometric_object'][field] += 1e-12
    with pytest.raises(ValueError, match='numeric'):
        h.compare_execution(changed, arrays, report, arrays)


@pytest.mark.parametrize('kind', ['missing', 'extra', 'nonfinite', 'altered'])
def test_publication_helper_arrays_cannot_bypass_replay(execution, kind):
    report, arrays = execution
    fresh = deepcopy(arrays)
    if kind == 'missing':
        fresh.pop('transported_prefixes')
    elif kind == 'extra':
        fresh['undeclared_helper'] = np.zeros(1)
    elif kind == 'nonfinite':
        fresh['transported_prefixes'][0, 0, 0] = np.nan
    else:
        fresh['transported_prefixes'][0, 0, 0] += .1
    with pytest.raises(ValueError):
        h.compare_execution(report, arrays, report, fresh)


def test_publication_replay_requires_unchanged_source_pins(execution, tmp_path):
    report, arrays = execution
    binary = tmp_path/'source-pins.npz'
    path = binary.with_suffix('.json')
    np.savez_compressed(binary, **arrays)
    pins = h.source_bindings()
    changed = {**pins, 'altered_source': '0'*64}
    saved = {**report, 'source_bindings': pins,
        'artifact': {'filename': binary.name, 'sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}}
    path.write_text(json.dumps(saved))
    with patch.object(h, 'execute', return_value=execution), patch.object(h, 'source_bindings', side_effect=[pins, changed]):
        with pytest.raises(ValueError, match='changed during replay'):
            h.check(path)


def test_publication_source_pins_cover_exact_object_dependencies():
    assert set(h.object_geometry.source_bindings()) <= set(h.source_bindings())

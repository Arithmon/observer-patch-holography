"""Complete-record phase access, surviving work and attachment boundaries."""
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coherent_native_access_probe as producer
import verify_coherent_native_access_independent as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.m.load(producer.REFERENCE)


def test_source_replay_and_different_arithmetic_route(packet):
    assert producer.build() == packet
    assert verifier.verify(packet)['raw_blocks_replayed'] == 7


def test_positive_phase_witness_and_real_effect_negative_controls(packet):
    assert packet['menu']['settings'] == 4291
    assert packet['menu']['ideal_record_distinguishing_settings'] == {'basis': 0, 'mixer': 0, 'y_analyzer': 2126}
    selected = packet['selected']
    assert selected['pair'] == [0, 6]
    assert selected['ideal_label_0_difference'] == ['0', '-1/100']
    assert selected['ideal_record_TV'] == ['1/8', '13/200']
    assert selected['ideal_plus'][67] == selected['ideal_minus'][67] == ['9/80', '-1/800']
    for law in ('ideal_plus', 'ideal_minus'):
        values = [tuple(map(Q, x)) for x in selected[law]]
        assert len(values) == 68 and tuple(map(sum, zip(*values))) == (1, 0)


def test_source_auxiliary_coupling_retains_phase_in_quantum_work(packet):
    pointer = packet['existing_native_work_pointer']
    assert pointer['source_word_pulses'] == 6
    assert pointer['logical_superoperator_entries_checked'] == 64
    assert len(pointer['bit_reports']) == 9
    assert all(x['ideal_record_marginals_equal'] and x['ideal_surviving_quantum_work_distinguishes']
               and not x['GV_code_invariant_under_bit_projection'] for x in pointer['bit_reports'])
    assert pointer['bit_reports'][7]['retained_entry'] == [0, 6]
    assert pointer['reset_timeout_upper'] == '35/2147483648'
    assert pointer['reset_routing_diameter'] == 4


def test_all_sixty_marked_actions_are_input_actions_not_a_decoder(packet):
    mark = packet['marked_action']
    assert (mark['current_index'], mark['historical_index']) == (4, 1)
    assert mark['input_actions_checked'] == 60
    assert mark['fixed_score_orbit_size']*mark['fixed_score_stabilizer_count'] == 60
    assert not mark['output_triplet_action_supplied']
    assert not packet['selected']['word_preserves_GV_code']
    assert packet['selected']['padding_leakage_control']['probability_lower'] == ['1/80', '-1/200']


@pytest.mark.parametrize('path,value', [
    (('source_pins', 'model/native_interface/pointer.py', 'commit'), '0'*40),
    (('source_pins', 'model/compiled_native_interface/construction.py', 'bytes'), 0),
    (('prior_probe_pins', 'code/a5_closure/receipts/c12_current_attachment_probe.json', 'sha256'), '0'*64),
    (('menu', 'settings'), 4290),
    (('menu', 'ideal_record_distinguishing_settings', 'mixer'), 1),
    (('menu', 'minimum_first_distinguishing_label_gap'), ['0', '0']),
    (('menu', 'complete_ideal_record_differences_sha256'), '0'*64),
    (('menu', 'routes_sha256'), '0'*64),
    (('selected', 'selection_rule'), 'choose the expected golden character'),
    (('selected', 'pair'), [0, 1]),
    (('selected', 'kind'), 'mixer'),
    (('selected', 'ideal_label_0_difference'), ['0', '1/100']),
    (('selected', 'word_preserves_GV_code'), True),
    (('selected', 'padding_leakage_control', 'probability_lower'), ['0', '0']),
    (('marked_action', 'current_index'), 1),
    (('marked_action', 'historical_index'), 4),
    (('marked_action', 'input_actions_checked'), 1),
    (('marked_action', 'output_triplet_action_supplied'), True),
    (('marked_action', 'action_twist_selected'), True),
    (('service_bounds', 'physical_service_error_included'), True),
    (('service_bounds', 'complete_record_half_l1_upper'), '0'),
    (('service_bounds', 'arbitrary_untouched_reference_bound'), False),
    (('service_bounds', 'unknown_input_code_preparation_is_supplied_not_derived'), False),
    (('service_bounds', 'exact_raw_equality_claimed_for_ideal_null_settings'), True),
    (('output_contract', 'record_service_discards_final_native_work'), False),
    (('output_contract', 'finer_controller_timeout_quantum_outputs_claimed_equal'), True),
    (('output_contract', 'tomography_quantum_payload_exported'), True),
    (('output_contract', 'tomography_uses_coherent_native_work_coupling'), True),
    (('output_contract', 'tomography_uses_uncommitted_native_records'), True),
    (('existing_native_work_pointer', 'input_auxiliary_independent_of_unknown_work_and_reference'), False),
    (('existing_native_work_pointer', 'reset_timeout_included_in_complete_reduced_work_comparison'), False),
    (('existing_native_work_pointer', 'finer_quantum_auxiliary_outputs_reconstructed_by_probe'), True),
    (('existing_native_work_pointer', 'finer_quantum_auxiliary_outputs_claimed_equal'), True),
    (('existing_native_work_pointer', 'complete_reduced_work_error_upper'), '0'),
    (('existing_native_work_pointer', 'coherent_GV_decoder_or_internal_triplet_export_supplied'), True),
    (('physical_golden_embedding_selected',), True),
    (('complete_operational_attachment_certified',), True),
    (('universal_coherent_grammar_obstruction_claimed',), True),
    (('phase_0_stop_preserved',), False),
    (('conjugation',), 'Galois-conjugate physical realizations'),
])
def test_hostile_scope_and_source_mutations(packet, path, value):
    bad = deepcopy(packet)
    target = bad
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        verifier.verify(bad)


@pytest.mark.parametrize('attack', ['abort', 'timeout', 'quantum_work', 'normalize_success'])
def test_omitted_or_renormalized_outputs_fail(packet, attack):
    bad = deepcopy(packet)
    if attack in ('abort', 'timeout'):
        at = 67 if attack == 'abort' else 66
        bad['selected']['ideal_plus'].pop(at)
        bad['selected']['ideal_minus'].pop(at)
        bad['output_contract']['labels'].pop(at)
    elif attack == 'quantum_work':
        bad['existing_native_work_pointer']['bit_reports'][7]['ideal_surviving_quantum_work_distinguishes'] = False
    else:
        bad['selected']['ideal_plus'][67] = ['0', '0']
    with pytest.raises(ValueError):
        verifier.verify(bad)


def test_changed_native_route_fails_before_receipt_comparison():
    producer.exact_probe.cache_clear()
    with patch.object(producer.tomography, 'route', return_value=()):
        with pytest.raises(ValueError, match='pair isolation'):
            producer.exact_probe()
    producer.exact_probe.cache_clear()


def test_raw_gate_rejects_missing_phase_or_incorrect_inverse():
    import numpy as np
    # Integer targets give an exact positive norm-control fixture; neither
    # this identity nor the reference is an added operational service.
    identity = np.eye(2, dtype=complex)
    target = {(i, i): (Q(1), Q(), Q(), Q()) for i in range(2)}
    assert verifier.interval_block_gate(identity, target) == 0
    bad = identity.copy(); bad[0, 0] = 1j
    with pytest.raises(ValueError, match='original-binary'):
        verifier.interval_block_gate(bad, target)

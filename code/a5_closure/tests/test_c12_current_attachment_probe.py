"""Source-interface falsifiers; the quantum code states are unit controls."""
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import c12_current_attachment_probe as producer
import verify_c12_current_attachment_independent as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.m.load(producer.REFERENCE)


def test_source_replay_and_independent_proof(packet):
    assert producer.build() == packet
    assert verifier.verify(packet)['verified']


def test_shared_mean_law_does_not_claim_a_state_loader(packet):
    join = packet['mean_observable_join']
    assert join['mean_zero_rank'] == 11
    assert join['native_events_checked'] == 3960
    assert not join['state_encoding_operation_implemented']
    assert packet['slow_code_control']['native_code_Gram_scale'] == ['3600', '720']


def test_actual_preparation_fails_the_natural_centroid_identification(packet):
    row = packet['C12_preparation']['centroid_dictionary_control']
    assert row['native_initial_centroid'] == [['-1/2', '-1/2'], ['0', '0'], ['0', '0']]
    assert row['ideal_success_work_centroid'] == [['0', '0']]*3
    assert not row['automatic_identification_preserved']


def test_full_rank_preparation_retains_success_rejection_and_every_first_failure(packet):
    p = packet['C12_preparation']
    mass = Q(p['request_rejection'])+sum(map(Q, p['success_by_request']))+sum(map(Q, p['first_failure_by_attempt']))
    assert mass == 1
    assert Q(p['minimum_unnormalized_eigenvalue_lower']) > 0
    assert p['rank_of_ideal_success_work_state'] == 12
    assert p['normalized_slow_block'].startswith('P_s/3')
    assert len(p['success_by_request']) == 5 and len(p['first_failure_by_attempt']) == 4
    assert p['selected_slow_Hermitian_observable_rank'] == 3 < p['full_slow_Hermitian_dimension']


def test_source_native_control_census_is_complete(packet):
    rows = packet['C12_preparation']['controller_census']
    assert [r['count'] for r in rows] == [2500, 500, 300, 60, 200, 40]
    assert sum(r['count'] for r in rows) == 3600
    assert packet['native_record_grammar']['matrix_units_checked'] == 182520


def test_terminal_native_coherence_is_a_load_bearing_exclusion(packet):
    # The zero-step identity on the EXISTING native space retains the two
    # different imaginary entries. Dephasing hides them; physical erasure of
    # native coherence is not inferred from inaccessible work-only outputs.
    pair = packet['slow_code_control']
    value = pair['nonzero_imaginary_entry']['plus_coefficient']
    assert any(Q(x) != 0 for x in value)
    assert pair['mutually_orthogonal'] and pair['trace_each'] == pair['purity_each'] == 1
    assert not pair['coherent_native_output_retained_by_receiver']
    assert pair['complete_deterministic_transfer_worst_case_halftrace_floor'] == '1/2'
    assert 'not Galois' in pair['witness_conjugation']


@pytest.mark.parametrize('attack', [
    'kernel_hash', 'G_hash', 'mean_rank', 'state_loader', 'mark', 'witness_scale',
    'witness_diagonal', 'imaginary_entry', 'reference_alias', 'native_output',
    'coherent_coupling', 'preexisting_payload_in_work', 'lost_success_assumption',
    'drop_rejection', 'drop_first_failure', 'renormalize_each_component',
    'pure_slow_state', 'full_slow_tomography', 'feedback_encoding', 'centroid_match',
    'raw_error_certification', 'source_commit', 'source_bytes', 'prior_pin',
    'physical_selection', 'universal_no_go', 'Galois_rebranding', 'scope_extension',
    'native_coherent_control', 'uncommitted_record_access',
])
def test_hostile_mutations_fail_closed(packet, attack):
    bad = deepcopy(packet)
    mean, pair, prep = bad['mean_observable_join'], bad['slow_code_control'], bad['C12_preparation']
    if attack in ('kernel_hash', 'G_hash'):
        mean['kernel_count_sha256' if attack == 'kernel_hash' else 'G_sha256'] = '0'*64
    elif attack == 'mean_rank':
        mean['mean_zero_rank'] = 12
    elif attack == 'state_loader':
        mean['state_encoding_operation_implemented'] = True
    elif attack == 'mark':
        pair['current_marking_index'] = 1
    elif attack == 'witness_scale':
        pair['native_code_Gram_scale'] = ['1', '0']
    elif attack == 'witness_diagonal':
        pair['shared_diagonal_sha256'] = '0'*64
    elif attack == 'imaginary_entry':
        pair['nonzero_imaginary_entry']['plus_coefficient'] = ['0', '0']
    elif attack == 'reference_alias':
        prep['supplied_work_reference_is_native_reference'] = True
    elif attack == 'native_output':
        pair['coherent_native_output_retained_by_receiver'] = True
    elif attack == 'coherent_coupling':
        pair['additional_coherent_native_work_access_admitted'] = True
    elif attack == 'preexisting_payload_in_work':
        pair['work_input_independent_of_unknown_native_payload'] = False
    elif attack == 'lost_success_assumption':
        pair['uniform_herald_success_p_floor'] = 'positive conditional fidelity regardless of success probability'
    elif attack == 'drop_rejection':
        prep['request_rejection'] = '0'
    elif attack == 'drop_first_failure':
        prep['first_failure_by_attempt'].pop()
    elif attack == 'renormalize_each_component':
        prep['success_by_request'] = prep['weights']
    elif attack == 'pure_slow_state':
        prep['rank_of_ideal_success_work_state'] = 3
        prep['outside_slow_population'] = ['0', '0']
    elif attack == 'full_slow_tomography':
        prep['selected_slow_Hermitian_observable_rank'] = 9
    elif attack == 'feedback_encoding':
        prep['ideal_slow_block_encodes_native_feedback'] = True
    elif attack == 'centroid_match':
        prep['centroid_dictionary_control']['native_initial_centroid'] = [['0', '0']]*3
    elif attack == 'raw_error_certification':
        prep['raw_numerical_guards_promoted_to_rigorous_error_bounds'] = True
    elif attack in ('source_commit', 'source_bytes'):
        row = bad['source_pins']['model/observer/native.py']
        row['commit' if attack == 'source_commit' else 'bytes'] = '0'*40 if attack == 'source_commit' else 0
    elif attack == 'prior_pin':
        bad['prior_probe_pins'][producer.PRIOR_FILES[0]]['sha256'] = '0'*64
    elif attack == 'physical_selection':
        bad['physical_golden_embedding_selected'] = True
    elif attack == 'universal_no_go':
        bad['universal_OPH_or_compiler_impossibility'] = True
    elif attack == 'Galois_rebranding':
        pair['witness_conjugation'] = 'physical Galois embedding discrimination'
    elif attack in ('native_coherent_control', 'uncommitted_record_access'):
        key = ('additional_coherent_native_controls_admitted' if attack == 'native_coherent_control'
               else 'coherent_uncommitted_native_record_access_admitted')
        bad['native_record_grammar'][key] = True
    else:
        bad['native_record_grammar']['work_controls'] = 'arbitrary coherent native-work operations'
    with pytest.raises(ValueError):
        verifier.verify(bad)


def test_independent_native_grammar_rejects_coherent_merging_and_lost_coin_branches():
    good = {(0, 0): [(0, 0), (1, 1)], (1, 0): [(0, 0), (1, 1)]}
    verifier.validate_partial_injections(good, 2)
    for bad in ({(0, 0): [(0, 0), (1, 0)], (1, 0): [(0, 0), (1, 1)]},
                {(0, 0): [(0, 0), (1, 1)]}):
        with pytest.raises(ValueError):
            verifier.validate_partial_injections(bad, 2)


def test_a_changed_native_primitive_fails_the_mean_gate_before_any_receipt_comparison():
    with patch.object(producer.native, 'native', side_effect=lambda x, seam, coin: (x, (coin, 0))):
        with pytest.raises(ValueError, match='mean transfer'):
            producer.native_join()

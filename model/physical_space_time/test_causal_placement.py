"""Controls for actual-source geometry, dependency records and cone scope."""
from itertools import product
from copy import deepcopy

import pytest
import sympy as sp

from . import causal_placement as c


def test_gram_is_the_existing_rank_three_response_not_an_added_mesh():
    g = c.gram()
    assert g.trace() == 12
    assert (g*g-4*g).applyfunc(sp.expand) == sp.zeros(12)
    assert set(g.diagonal()) == {1}


def test_complete_tape_parents_include_zero_current_events():
    rows = c.dependencies(((0, 1), (2, 3), (0, 1), (0, 1)))
    assert [r['parents'] for r in rows] == [(), (), (0,), (2,)]
    assert [r['height'] for r in rows] == [0, 0, 1, 2]
    assert rows[-1]['ancestors'] == (0, 2)
    assert c.interval_squared(rows[1], rows[3]) == 0
    assert 1 not in rows[3]['ancestors']


def test_nontrivial_two_endpoint_parent_merge_and_forward_cones():
    seams, _ = c.passage.recovered_carrier()
    first, last = seams[0], seams[-1]
    tape = (*seams, first, last, *reversed(seams))
    rows = c.dependencies(tape)
    assert any(len(r['parents']) == 2 for r in rows)
    for row in rows:
        for ancestor in row['ancestors']:
            assert row['height'] > rows[ancestor]['height']
            assert c.interval_squared(rows[ancestor], row).is_nonnegative


def test_all_native_converse_branches_are_retained_and_counts_are_exact():
    result = c.run()
    control = result['converse_control']
    assert [r['coins'] for r in control['all_coin_branches']] == list(map(list, product((0, 1), repeat=4)))
    assert any(record[1] == 0 for r in control['all_coin_branches'] for record in r['records'])
    assert result['finite_current_readout'] == {
        'total_two_states': 78, 'distinct_rank_three_images': 73,
        'binary_states': 66, 'distinct_binary_images': 61,
        'native_branch_telescoping_checks': 4680,
        'unit_centroid_normalization': 'raw_load_readout / (2*sqrt((5+sqrt(5))/2))',
        'normalized_load_centroid_in_convex_hull_of_unit_ports': True,
        'waiting_on_fixed_total_carrier_realizes_full_signed_control_completion': False}
    pair = result['forward_cone']['sharp_pair']
    assert c.distance_squared(*pair) == c.SPEED_SQUARED
    assert not result['acceptance']['cone_implies_ancestry']
    assert result['Lean_native_branch_binding']['branches_compared'] == 16
    corrupted = deepcopy(control['all_coin_branches'])
    corrupted[0]['records'][0][1] = 0
    with pytest.raises(ArithmeticError, match='Lean native execution table'):
        c.lean_branch_binding(corrupted)


@pytest.mark.parametrize('seam', [(True, 1), (1, 0), (0, 0), (0, 2), (0, 12), '01'])
def test_reject_foreign_seams(seam):
    with pytest.raises(ValueError):
        c.dependencies([seam])


@pytest.mark.parametrize('height', [True, -1, 0.5, 1.0, '1'])
def test_reject_inexact_or_foreign_height(height):
    with pytest.raises(ValueError):
        c.interval_squared({'seam': (0, 1), 'height': height},
                           {'seam': (0, 1), 'height': 2})


INITIAL = tuple(int(i in (0, 2)) for i in range(12))
TAPE = ((0, 1), (2, 3), (0, 1), (0, 1))


def test_actual_native_logs_authenticate_every_coin_branch_and_preserve_zero_currents():
    zero_commits = 0
    for coins in product((0, 1), repeat=4):
        log = c.semantic_execution(INITIAL, TAPE, coins)
        assert c.verify_semantic_execution(INITIAL, log) == c.dependencies(TAPE)
        assert c.authenticated_parent_ports(log[0], log[2]) == (0, 1)
        assert c.authenticated_parent_ports(log[1], log[3]) == ()
        assert c.authenticated_parent_ports(log[2], log[3]) == (0, 1)
        for entry in log:
            if entry['record'][1] == 0:
                zero_commits += 1
                assert entry['before_values'] == entry['after_values']
                assert all(entry['after_writers'][p] == entry['event_id'] for p in entry['write_set'])
    assert zero_commits > 0


def test_authenticated_log_for_merging_endpoint_histories_and_empty_history():
    seams, _ = c.passage.recovered_carrier()
    tape = (*seams, *reversed(seams))
    log = c.semantic_execution(INITIAL, tape, [i % 2 for i in range(len(tape))])
    rows = c.verify_semantic_execution(INITIAL, log)
    assert rows == c.dependencies(tape)
    assert any(len(row['parents']) == 2 for row in rows)
    assert c.semantic_execution(INITIAL, [], []) == []


def test_native_local_support_is_checked_on_the_complete_carrier():
    assert c.semantic_support() == {
        'actual_native_branch_checks': 4680, 'necessary_endpoint_witnesses': 120,
        'support_sound_and_irredundant_within_two_endpoint_read_interface': True,
        'unique_minimal_support_among_global_read_interfaces_claimed': False}


@pytest.mark.parametrize('mutation', [
    'stale_writer', 'forged_value', 'missing_write', 'missing_support',
    'outside_stamp', 'duplicate_id', 'false_current', 'false_coin', 'omit_zero',
])
def test_reject_unthreaded_forged_or_incomplete_semantic_logs(mutation):
    log = deepcopy(c.semantic_execution(INITIAL, TAPE, (1, 1, 1, 1)))
    if mutation == 'stale_writer':
        log[3]['before_writers'][0] = 0
    elif mutation == 'forged_value':
        log[2]['before_values'][0], log[2]['before_values'][1] = 0, 1
    elif mutation == 'missing_write':
        log[0]['write_set'] = [1]
    elif mutation == 'missing_support':
        log[2]['causal_support'] = [1]
    elif mutation == 'outside_stamp':
        log[0]['after_writers'][2] = 0
    elif mutation == 'duplicate_id':
        log[2]['event_id'] = 0
    elif mutation == 'false_current':
        log[0]['record'][1] = 1
    elif mutation == 'false_coin':
        log[0]['coin'] = True
    else:
        del log[2]
    with pytest.raises(ValueError):
        c.verify_semantic_execution(INITIAL, log)


def test_writer_citation_alone_does_not_authenticate_an_edge():
    log = c.semantic_execution(INITIAL, TAPE, (1, 1, 1, 1))
    parent, child = deepcopy(log[0]), deepcopy(log[2])
    parent['write_set'] = []
    assert c.authenticated_parent_ports(parent, child) == ()
    parent = deepcopy(log[0])
    child['before_values'][0], child['before_values'][1] = 0, 1
    assert c.authenticated_parent_ports(parent, child) == ()


def test_external_experiment_coverage_rejects_a_truncated_valid_execution():
    coins = (1, 1, 1, 1)
    log = c.semantic_execution(INITIAL, TAPE, coins)
    assert c.verify_semantic_execution(INITIAL, log[:-1]) == c.dependencies(TAPE[:-1])
    with pytest.raises(ValueError, match='complete declared experiment'):
        c.verify_semantic_execution(INITIAL, log[:-1], tape=TAPE, coins=coins)

"""Matched apparatus schedules, actual CP branches and identification limits."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest

from model.identifiability_checks import consensus_observability as c


@pytest.fixture(scope='module')
def program():
    # Small cutoff tests the solver mechanics, not the retained empirical score.
    return c.selected_program(dimension=4, subdivisions=2)


def test_actual_apparatus_all_four_outcomes_agree_across_completed_schedules(program):
    audit, _ = c.audit_program(program)
    assert audit['complete_terminal_records'] == 4
    assert audit['causal_depth'] == 11
    assert audit['node_count_including_boundary'] == 15
    assert [r['updates'] for r in audit['schedules']] == [14, 154, 154]
    assert all(r['maximum_terminal_matrix_difference'] < 2e-12 for r in audit['schedules'])
    assert not audit['unfair_control']['closed_within_tolerance']


def test_arbitrary_inconsistent_assignments_close_without_terminal_oracle(program, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('terminal execution cannot supply the repair solution')
    monkeypatch.setattr(c.runtime, 'execute', forbidden)
    monkeypatch.setattr(c.data, 'load_validation', forbidden)
    system = c.graph(program)
    rng = np.random.default_rng(27)
    initial = {index: rng.normal(size=program.initial_density.shape)
               +1j*rng.normal(size=program.initial_density.shape)
               for index in range(1, len(system.nodes)+1)}
    forward = c.repair(system, c.schedule(system, 'forward'))
    reverse = c.repair(system, c.schedule(system, 'reverse'), initial)
    assert reverse['closed_within_tolerance']
    for index in forward['assignment']:
        assert np.array_equal(forward['assignment'][index], reverse['assignment'][index])
    assert np.array_equal(reverse['assignment'][0], program.initial_density)


def test_unfinished_schedule_and_changed_boundary_are_not_certified(program):
    system = c.graph(program)
    result = c.repair(system, c.schedule(system, 'reverse', covers=system.depth-1))
    assert not result['closed_within_tolerance']
    with pytest.raises(ValueError):
        c.repair(system, (0,))
    with pytest.raises(ValueError):
        c.repair(system, (), {0: program.initial_density})


def test_different_supplied_actions_both_close_and_change_actual_detector_law(program):
    harmonic = c.selected_program(dimension=4, subdivisions=2, harmonic=True)
    first, second = c.graph(program), c.graph(harmonic)
    p = c.terminal_law(first, c.repair(first, c.schedule(first, 'forward')))
    q = c.terminal_law(second, c.repair(second, c.schedule(second, 'forward')))
    assert p.keys() == q.keys()
    assert abs(sum(p.values())-1) < 2e-11
    assert abs(sum(q.values())-1) < 2e-11
    assert sum(abs(p[k]-q[k]) for k in p)/2 > 1e-4


def test_separate_endpoints_share_only_unread_prefixes():
    apparatus, alpha, thermal = c.apparatus_inputs(dimension=4, subdivisions=2)
    system, groups = c.endpoint_graph(apparatus, alpha, thermal, (0, 1, 2))
    result = c.repair(system, c.schedule(system, 'reverse'))
    assert result['closed_within_tolerance']
    assert len(groups) == 3
    for step, indices in groups.items():
        direct = c.runtime.execute(c.anharmonic.program(apparatus, alpha, thermal, step))
        expected = {b.record[-1]: b.density for b in direct.branches}
        for label, index in enumerate(indices):
            np.testing.assert_allclose(result['assignment'][index], expected[label], atol=2e-12, rtol=0)
        assert abs(sum(np.trace(result['assignment'][i]).real for i in indices)-1) < 2e-11


def test_feedback_selects_from_parent_record_and_retains_all_zero_branches():
    identity = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    projectors = tuple((np.diag(row),) for row in ([1., 0.], [0., 1.]))
    read = c.runtime.Instrument(kraus=projectors, visible=True)
    stay = c.runtime.Instrument(kraus=((identity,),))
    flip = c.runtime.Instrument(kraus=((x,),))
    events = (c.runtime.Event('read', instruments=(read,)),
              c.runtime.Event('feedback', instruments=(stay, flip), control=0),
              c.runtime.Event('read_again', instruments=(read,)))
    program = c.runtime.Program('controlled', 2, identity/2, events)
    system = c.graph(program)
    result = c.repair(system, c.schedule(system, 'reverse'))
    assert c.terminal_law(system, result) == {(0, 0): .5, (0, 1): 0., (1, 0): .5, (1, 1): 0.}
    direct = c.runtime.execute(program)
    for index in system.terminals:
        record = system.nodes[index-1].record
        target = next(branch.density for branch in direct.branches if branch.record == record)
        assert np.array_equal(result['assignment'][index], target)


def retained():
    return json.loads(Path(__file__).with_name('consensus_observability_results.json').read_text())


def test_retained_full_primary_binds_inputs_all_endpoints_and_prediction_seal():
    receipt = retained()
    c.compare_receipt(receipt, receipt)
    for name, digest in receipt['source_sha256'].items():
        assert hashlib.sha256((c.ROOT/name).read_bytes()).hexdigest() == digest
    full = receipt['full_primary_empirical_comparison']
    assert full['endpoint_count'] == len(full['predictions']) == 130
    assert {(row['case'], row['time_ms']) for row in full['predictions']} == {
        (case, step/5) for case in 'abcde' for step in range(0, 51, 2)}
    assert full['scores']['metrics']['combined']['observations'] == 125
    assert full['scores']['metrics']['combined']['position_estimator_rmse'] == pytest.approx(
        .28764278683646793, abs=2e-9)
    assert all(check['maximum_complete_endpoint_matrix_difference'] < 2e-9
               for check in full['case_checks'])
    for row in full['predictions']:
        assert len(row['probabilities']) == 4
        assert sum(row['probabilities']) == pytest.approx(1., abs=2e-9)
    assert not full['new_independent_data']
    assert not receipt['scope']['physical_action_selected']


@pytest.mark.parametrize('corruption', ['prediction', 'seal', 'source_hash', 'missing_branch', 'nonfinite'])
def test_retained_replay_rejects_changed_predictions_counts_and_bindings(corruption):
    expected = retained()
    actual = deepcopy(expected)
    full = actual['full_primary_empirical_comparison']
    if corruption == 'prediction':
        full['predictions'][0]['predicted_position_estimator'] += .01
    elif corruption == 'seal':
        full['prediction_payload_sha256_before_scoring'] = '0'*64
    elif corruption == 'source_hash':
        actual['source_sha256'][next(iter(actual['source_sha256']))] = '0'*64
    elif corruption == 'missing_branch':
        full['predictions'][0]['probabilities'].pop()
    else:
        full['predictions'][0]['probabilities'][0] = float('nan')
    with pytest.raises(ValueError):
        c.compare_receipt(expected, actual)

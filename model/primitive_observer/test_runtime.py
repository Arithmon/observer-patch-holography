"""Independent small channels exercise the complete runtime contract."""
from dataclasses import replace
import json

import numpy as np
import pytest

from . import runtime as r


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
P = (np.diag([1., 0.]), np.diag([0., 1.]))
PLUS = np.ones((2, 2))/2


def unitary(matrix=I, cost=1):
    return r.Instrument(kraus=((matrix,),), resources=(('pulses', cost),))


def read():
    return r.Instrument(kraus=tuple((p,) for p in P), visible=True,
                        resources=(('pulses', 3), ('record_bits', 1)))


def program(events=(), initial=PLUS):
    return r.Program('independent_small_channel', 2, initial, tuple(events),
                     preparation_resources=(('pulses', 2),))


def test_all_zero_and_nonzero_outcomes_and_causal_feedback_remain_complete():
    events = (r.Event('read', instruments=(read(),)),
              r.Event('feedback', instruments=(unitary(I, 5), unitary(X, 5)), control=0),
              r.Event('later_read', instruments=(read(),)))
    history = r.execute(program(events))
    assert history.records == ((0, 0), (0, 1), (1, 0), (1, 1))
    np.testing.assert_allclose([b.probability for b in history.branches], [.5, 0, .5, 0], atol=1e-15)
    for branch in history.branches:
        np.testing.assert_allclose(branch.density, branch.probability*P[0], atol=1e-15)
        assert dict(branch.resources) == {'pulses': 13, 'record_bits': 2}
    assert dict(history.expanded_tree_resources) == {'pulses': 21, 'record_bits': 3}


def test_hidden_alternatives_are_cp_summed_not_added_as_amplitudes_or_records():
    hidden = r.Instrument(kraus=((I/np.sqrt(2), Z/np.sqrt(2)),))
    history = r.execute(program((r.Event('unobserved_phase', instruments=(hidden,)),)))
    assert history.records == ((),)
    np.testing.assert_allclose(history.branches[0].density, I/2, atol=2e-16)


def test_reference_minor_schur_action_agrees_with_independent_lifted_kraus_sum():
    schur = np.array([[1., .3-.2j], [.3+.2j, 1.]])
    hidden = r.Instrument(schur=schur)
    rng = np.random.default_rng(671)
    v = rng.normal(size=6)+1j*rng.normal(size=6)
    v /= np.linalg.norm(v)
    initial = np.outer(v, v.conj())
    history = r.execute(program((r.Event('noise', instruments=(hidden,)),)), initial=initial)
    values, vectors = np.linalg.eigh(schur)
    expected = sum((np.kron(np.diag(np.sqrt(lam)*col), np.eye(3)) @ initial
                    @ np.kron(np.diag(np.sqrt(lam)*col), np.eye(3)).conj().T
                    for lam, col in zip(values, vectors.T, strict=True)), np.zeros_like(initial))
    np.testing.assert_allclose(history.branches[0].density, expected, atol=3e-16)
    before_reference = np.trace(initial.reshape(2, 3, 2, 3), axis1=0, axis2=2)
    after_reference = np.trace(history.branches[0].density.reshape(2, 3, 2, 3), axis1=0, axis2=2)
    np.testing.assert_allclose(before_reference, after_reference, atol=1e-16)
    assert history.reference_dimension == 3
    assert history.preparation_mode == 'supplied'
    assert dict(history.resources) == {}


def test_unitary_reference_and_repeated_hidden_sequence_match_literal_matrix_execution():
    u = np.array([[1., 1j], [1j, 1.]])/np.sqrt(2)
    phase = r.Instrument(schur=np.array([[1., .8], [.8, 1.]]), resources=(('draws', 1),))
    event = r.Event('repeated', sequence=(phase, unitary(u, 7), phase), repeat=4)
    v = np.array([1., 0., 0., 1.])/np.sqrt(2)
    initial = np.outer(v, v)
    history = r.execute(program((event,)), initial=initial)
    expected = initial.astype(complex)
    lifted = np.kron(u, I)
    numbers = np.repeat(np.arange(2), 2)
    s = phase.schur[numbers[:, None], numbers[None, :]]
    for _ in range(4):
        expected *= s
        expected = lifted @ expected @ lifted.conj().T
        expected *= s
    np.testing.assert_allclose(history.branches[0].density, expected, atol=1e-15)
    assert len(history.stages) == 1
    assert dict(history.branches[0].resources) == {'draws': 8, 'pulses': 28}


def test_full_density_checkpoint_continuation_and_json_replay_are_identical():
    p = program((r.Event('read', instruments=(read(),)),
                 r.Event('feedback', instruments=(unitary(), unitary(X)), control=0),
                 r.Event('later', instruments=(read(),))))
    prefix = r.execute(p, stop_after=1)
    saved = r.checkpoint_to_json(prefix)
    restored = r.checkpoint_from_json(p, saved)
    resumed = r.execute(p, checkpoint=restored)
    full = r.execute(p)
    assert r.complete_distance(resumed, full) == 0
    assert resumed.expanded_tree_resources == full.expanded_tree_resources
    assert r.checkpoint_to_json(resumed) == r.checkpoint_to_json(full)


def test_resealed_forged_quantum_state_does_not_replace_replay():
    p = program((r.Event('identity', instruments=(unitary(),)),))
    history = r.execute(p)
    corrupt = replace(history, branches=(replace(history.branches[0], density=I/2),), seal='')
    corrupt = replace(corrupt, seal=r._seal(corrupt))
    saved = r.checkpoint_to_json(corrupt)
    with pytest.raises(ValueError, match='replay'):
        r.checkpoint_from_json(p, saved)


def test_resealed_false_resource_count_does_not_replace_replay():
    p = program((r.Event('identity', instruments=(unitary(),)),))
    history = r.execute(p)
    corrupt = replace(history, expanded_tree_resources=(('pulses', 99),), seal='')
    corrupt = replace(corrupt, seal=r._seal(corrupt))
    with pytest.raises(ValueError, match='replay'):
        r.checkpoint_from_json(p, r.checkpoint_to_json(corrupt))


def test_checkpoint_identity_prevents_source_switch_input_replacement_and_rewind():
    p = program((r.Event('first', instruments=(unitary(),)),))
    q = program((r.Event('first', instruments=(unitary(X),)),))
    history = r.execute(p)
    with pytest.raises(ValueError):
        r.checkpoint_from_json(q, r.checkpoint_to_json(history))
    for kwargs in ({'initial': PLUS}, {'stop_after': 0}):
        with pytest.raises(ValueError):
            r.execute(p, checkpoint=history, **kwargs)


@pytest.mark.parametrize('bad', [True, [[True, 0], [0, 1]], [[float('nan'), 0], [0, 1]],
                                np.ma.array(I, mask=False), [[1, 0, 0], [0, 1, 0]], I*2])
def test_invalid_quantum_maps_fail_closed(bad):
    with pytest.raises(ValueError):
        r.Instrument(kraus=((bad,),))


@pytest.mark.parametrize('bad', [I*2, [[1., 2.], [2., 1.]], [[1., 1j], [1j, 1.]]])
def test_non_cp_or_non_tp_schur_matrices_rejected(bad):
    with pytest.raises(ValueError):
        r.Instrument(schur=bad)


@pytest.mark.parametrize('bad', [(('pulses', True),), (('pulses', -1),),
                                (('pulses', 2), ('pulses', 2)), (('pulses', 10**19),),
                                [('pulses', 2)]])
def test_invalid_resource_ledgers_rejected(bad):
    with pytest.raises(ValueError):
        r.Instrument(kraus=((I,),), resources=bad)


def test_causal_and_complete_outcome_domains_fail_before_execution():
    with pytest.raises(ValueError):
        program((r.Event('future', instruments=(unitary(), unitary()), control=0),))
    with pytest.raises(ValueError):
        r.Event('repeated_measurement', sequence=(read(),), repeat=2)
    with pytest.raises(ValueError):
        r.Event('different_menus', instruments=(unitary(), read()), control=0)
    with pytest.raises(ValueError):
        r.Instrument(kraus=((P[0],),), visible=True)
    p = program(tuple(r.Event(f'read_{i}', instruments=(read(),)) for i in range(8)))
    with pytest.raises(ValueError):
        r.execute(p, max_branches=2)
    with pytest.raises(ValueError):
        program((*p.events, r.Event('too_many', instruments=(read(),))))


def test_immutable_state_and_revalidation_detach_even_noop_resume():
    source = I.copy()
    instrument = unitary(source)
    source[0, 0] = 4
    with pytest.raises(ValueError):
        instrument.kraus[0][0].flags.writeable = True
    p = program((r.Event('first', instruments=(instrument,)),))
    prefix = r.execute(p, stop_after=0)
    resumed = r.execute(p, checkpoint=prefix, stop_after=0)
    assert resumed.program is not prefix.program
    assert resumed.branches[0].density is not prefix.branches[0].density
    p.events[0].instruments[0].__dict__['kraus'] = ((2*I,),)
    with pytest.raises(ValueError):
        r.execute(p)


@pytest.mark.parametrize('payload', ['{"cursor": 1, "cursor": 2}', '{"x": NaN}', '{"x": 1e999}'])
def test_nonfinite_duplicate_or_incomplete_json_fails(payload):
    with pytest.raises(ValueError):
        r.checkpoint_from_json(program(), payload)


def test_boolean_serialized_records_fail_even_after_resealing():
    p = program((r.Event('read', instruments=(read(),)),))
    h = r.execute(p)
    row = json.loads(r.checkpoint_to_json(h))
    row['branches'][0]['record'] = [False]
    with pytest.raises(ValueError):
        r.checkpoint_from_json(p, json.dumps(row))


def test_alternative_experiment_preserves_prefix_and_matches_fresh_complete_execution():
    p = program((r.Event('prepared_phase', instruments=(unitary(Z),)),
                 r.Event('unused_trunk_future', instruments=(unitary(X),))))
    prefix = r.execute(p, stop_after=1)
    original_payload = r.checkpoint_to_json(prefix)
    suffix = (r.Event('terminal', instruments=(read(),)),)
    alternative = r.branch_experiment(prefix, suffix, name='terminal_experiment', binding={'setting': 7})
    fresh = r.execute(alternative.program)
    assert r.complete_distance(alternative, fresh) == 0
    assert alternative.expanded_tree_resources == fresh.expanded_tree_resources
    assert r.checkpoint_to_json(prefix) == original_payload
    assert [e.name for e in alternative.program.events] == ['prepared_phase', 'terminal']
    restored = r.checkpoint_from_json(alternative.program, r.checkpoint_to_json(alternative))
    assert r.complete_distance(restored, alternative) == 0


def test_alternative_after_read_retains_backreaction_and_causal_feedback_records():
    prefix = r.execute(program((r.Event('committed', instruments=(read(),)),)))
    hadamard = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
    suffix = (r.Event('phase_read_rotation', instruments=(unitary(hadamard),)),
              r.Event('later', instruments=(read(),)))
    continued = r.branch_experiment(prefix, suffix, name='measured_continuation')
    assert continued.records == ((0, 0), (0, 1), (1, 0), (1, 1))
    np.testing.assert_allclose([b.probability for b in continued.branches], [.25]*4, atol=1e-15)
    unread = r.execute(program(suffix))
    np.testing.assert_allclose([b.probability for b in unread.branches], [1., 0.], atol=1e-15)
    assert r.complete_distance(continued, r.execute(continued.program)) == 0


def test_alternative_input_reference_and_tamper_guard():
    p = program((r.Event('trunk', instruments=(unitary(Z),)),))
    v = np.array([1., 0., 0., 1.])/np.sqrt(2)
    supplied = np.outer(v, v)
    prefix = r.execute(p, initial=supplied)
    alternative = r.branch_experiment(prefix, (r.Event('read', instruments=(read(),)),), name='ref_read')
    assert alternative.reference_dimension == 2
    assert alternative.preparation_mode == 'supplied'
    assert r.complete_distance(alternative, r.execute(alternative.program, initial=supplied)) == 0
    prefix.program.__dict__['name'] = 'tampered'
    with pytest.raises(ValueError, match='integrity'):
        r.branch_experiment(prefix, (), name='invalid')


def test_reused_compiled_instruments_are_bound_once_without_omitting_event_addresses():
    operation = unitary(X)
    p = program(tuple(r.Event(str(i), sequence=(operation, operation), repeat=3) for i in range(200)))
    specification = p.specification()
    assert len(specification['instruments']) == 1
    assert len(specification['events']) == 200
    assert all(event['sequence'] == [0, 0] and event['repeat'] == 3 for event in specification['events'])


def test_total_resource_overflow_fails_before_any_execution():
    expensive = unitary(cost=r.MAX_RESOURCE)
    with pytest.raises(ValueError):
        program((r.Event('too_expensive', sequence=(expensive,), repeat=2),))

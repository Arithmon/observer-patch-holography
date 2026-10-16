"""Full source seam menus, recorded controllers and source-backed continuation."""
from collections import Counter
import json
import random

import pytest

from ..observer import native
from . import clock_histories as c


@pytest.fixture(scope='module')
def source():
    return c.source()


@pytest.fixture(scope='module')
def example():
    return c.execute(1, ticks=2)


def test_union_census_recovers_each_entire_source_scheduler(source):
    assert source.ports == (15360, 245760)
    assert (len(source.left_edges), len(source.shared), len(source.right_only), source.slots) == (
        40320, 561, 644559, 684879)
    right = tuple(row[1] for row in source.shared)+source.right_only
    for level, menu in zip(c.LEVELS, (source.left_edges, right), strict=True):
        raw = c.screen.screen(level)
        actual = Counter(tuple(sorted((int(a), int(b))))
                         for a, b in zip(raw.seam_a, raw.seam_b, strict=True))
        assert Counter(menu) == actual
        assert set(actual.values()) == {1}
    # Thus each marginal non-holding seam has exactly one uniform slot;
    # there is no replacement of the full scheduler by a selected seam tape.
    assert len(source.incident_edge_indices) == 33
    assert len(source.observed_labels) == 15
    assert all(index in {row[0] for row in source.shared} for index in source.incident_edge_indices)


def test_preparation_is_a_shared_local_bernoulli_product_not_global_two(source):
    initial = c.execute(1, ticks=2, stop_after_slots=0)
    tables = tuple(dict(row) for row in source.matched_addresses)
    assert all(initial.states[0][tables[0][p]] == initial.states[1][tables[1][p]] for p in tables[0])
    assert all(set(state) == {0, 1} for state in initial.states)
    assert all(total > 2 for total in initial.initial_totals)
    assert initial.initial_totals == tuple(map(sum, initial.states))
    assert initial.initial_field_labels == tuple(initial.states[i][tables[i][0]] for i in range(2))


def _independent_history(seed, ticks, cap, feedback):
    """Apply original native rules in canonical order, not stored flip flags."""
    data = c.source()
    initial = c.execute(seed, ticks=ticks, cap=cap, feedback=feedback, stop_after_slots=0)
    states = [bytearray(row) for row in initial.states]
    rng = random.Random(); rng.setstate(initial.rng_state)
    tables = tuple(dict(row) for row in data.matched_addresses)
    shared = {row[0]: row[2] for row in data.shared}
    attempts, repairs, records = [0, 0], [0, 0], [[], []]
    slots = hits = 0

    def apply(side, ordered, coin):
        a, b = ordered
        after, (_, current) = native.native((states[side][a], states[side][b]), (0, 1), coin)
        states[side][a], states[side][b] = after
        return current

    def record(side, kind, edge, coin, current):
        records[side].append((kind, slots, hits, len(records[side])+1, *edge,
                              coin, current, states[side][tables[side][0]]))

    while slots < cap and hits < ticks:
        index, coin = rng.randrange(data.slots), rng.randrange(2)
        slots += 1
        if index >= len(data.left_edges):
            apply(1, data.right_only[index-len(data.left_edges)], coin)
            attempts[1] += 1
        elif index not in shared:
            apply(0, data.left_edges[index], coin)
            attempts[0] += 1
        else:
            edge = shared[index]
            currents = [apply(side, tuple(table[p] for p in edge), coin)
                        for side, table in enumerate(tables)]
            attempts = [value+1 for value in attempts]
            if min(edge) >= 12:
                continue
            hits += 1
            for side, current in enumerate(currents):
                record(side, 'scheduled', edge, coin, current)
            if hits == 1 and feedback:
                chosen_coin = rng.randrange(2)
                for side, table in enumerate(tables):
                    if states[side][table[0]]:
                        current = apply(side, (table[0], table[1]), chosen_coin)
                        repairs[side] += 1
                        record(side, 'feedback', (0, 1), chosen_coin, current)
    return tuple(map(bytes, states)), tuple(tuple(row) for row in records), tuple(attempts), tuple(repairs), rng.getstate()


@pytest.mark.parametrize('feedback', (False, True))
def test_full_states_and_all_recorded_source_branches_match_canonical_oracle(feedback):
    actual = c.execute(1, ticks=2, feedback=feedback)
    expected = _independent_history(1, 2, c.MAX_CAP, feedback)
    assert (actual.states, actual.records, actual.uniform_attempts, actual.feedback_attempts,
            actual.rng_state) == expected


def test_remote_repairs_change_real_states_and_feedback_is_a_clock_tick(source, example):
    initial = c.execute(1, ticks=2, stop_after_slots=0)
    tables = tuple(dict(row) for row in source.matched_addresses)
    for side in (0, 1):
        retained = set(tables[side].values())
        assert any(before != after and port not in retained for port, (before, after)
                   in enumerate(zip(initial.states[side], example.states[side], strict=True)))
        assert example.uniform_attempts[side] > 2
        assert example.feedback_attempts[side] == 1
        assert [row[0] for row in example.records[side]] == ['scheduled', 'feedback', 'scheduled']
        assert [row[3] for row in example.records[side]] == [1, 2, 3]
    report = c.summarize(example)
    assert report['status'] == 'complete'
    for row in report['screens']:
        assert row['total_native_attempts'] == row['uniform_native_attempts']+1
        assert row['total_incident_clock_ticks'] == 3
        assert 'comparison_slot' not in row['record_columns']
        assert all(len(record) == len(row['record_columns']) == 8 for record in row['incident_records'])


def test_timeout_is_retained_and_not_a_successful_local_clock_read():
    result = c.execute(1, ticks=2, cap=1)
    report = c.summarize(result)
    assert report['status'] == 'timeout'
    assert report['scheduled_hits_completed'] < 2
    assert report['comparison_slots_executed'] == 1
    assert all(row['initial_field_label'] in (0, 1) for row in report['screens'])
    assert all(len(row['terminal_observed_bits']) == 15 for row in report['screens'])
    assert c.execute(1, ticks=2, cap=1, checkpoint=result) == result


def test_full_serialized_prefix_continues_without_preparing_or_resetting(example):
    prefix = c.execute(1, ticks=2, stop_after_slots=example.slots//2)
    encoded = c.checkpoint_to_json(prefix)
    recovered = c.checkpoint_from_json(encoded)
    assert recovered.states == prefix.states
    assert c.execute(1, ticks=2, checkpoint=recovered) == example
    # JSON object-key ordering is irrelevant; actual field types are not.
    reordered = dict(reversed(list(json.loads(encoded).items())))
    assert c.checkpoint_from_json(json.dumps(reordered)) == prefix
    with pytest.raises(ValueError):
        c.execute(2, ticks=2, checkpoint=recovered)
    with pytest.raises(ValueError):
        c.execute(1, ticks=2, stop_after_slots=prefix.slots-1, checkpoint=recovered)


@pytest.mark.parametrize('warm', (False, True))
def test_nested_boolean_alias_is_not_an_authentic_source_record(example, warm):
    rows = list(example.records[0])
    first = list(rows[0])
    assert first[-1] == 1
    first[-1] = True
    rows[0] = tuple(first)
    malformed = example._replace(records=(tuple(rows), example.records[1]))
    if not warm:
        c._compute.cache_clear()
    with pytest.raises(ValueError, match='source replay'):
        c.validate_checkpoint(malformed)
    with pytest.raises(ValueError):
        c.validate_checkpoint(example._replace(uniform_attempts=(float(example.uniform_attempts[0]),
                                                               example.uniform_attempts[1])))


def test_forged_states_rng_and_source_binding_are_rejected(example):
    initial = c.execute(1, ticks=2, stop_after_slots=0)
    for changed in (example._replace(states=initial.states), example._replace(rng_state=initial.rng_state),
                    example._replace(source_digest='0'*64), example._replace(scheduled_hits=True)):
        with pytest.raises(ValueError):
            c.validate_checkpoint(changed)
    payload = json.loads(c.checkpoint_to_json(example))
    payload['records'][0][0][-1] = True
    with pytest.raises(ValueError):
        c.checkpoint_from_json(json.dumps(payload))
    payload['records'][0][0][-1] = 1
    payload['states'][0] = payload['states'][1]
    with pytest.raises(ValueError):
        c.checkpoint_from_json(json.dumps(payload))


@pytest.mark.parametrize('override', ({'seed': True}, {'seed': -1}, {'ticks': True}, {'ticks': 0},
    {'cap': True}, {'cap': c.MAX_CAP+1}, {'feedback': 1}, {'stop_after_slots': True}))
def test_parameter_aliases_rejected_before_cached_execution(override):
    parameters = {'seed': 1, 'ticks': 2, 'cap': 1}
    parameters.update(override)
    with pytest.raises(ValueError):
        c.execute(**parameters)

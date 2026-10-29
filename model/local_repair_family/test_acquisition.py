"""Independent event-by-event tape and complete-effect acquisition controls."""
from fractions import Fraction as F
import hashlib
import json

import numpy as np
import pytest

from . import acquisition as a


def parameters():
    return a.local.Parameters(port_weights=(25, 25)+(8,)*10,
                              c_weight=F(1, 2), repair_probability=F(1))


def scalar_oracle(p, n, seed, initial):
    """Every original floor/ceiling micro-event; no producer table or reduction."""
    rng = np.random.Generator(np.random.PCG64(seed))
    state = list(initial); tape = bytearray(a.STREAM_PREFIX); rows = []; counts = {}
    ordinal = 0
    for ei, (left, right) in enumerate(a.records.seams()):
        draws = rng.integers(0, p.port_weights[left]+p.port_weights[right], n, dtype=np.uint64)
        ones = 0; raw_counts = {'C': [[0]*4 for _ in range(2)], 'E': [[0]*8 for _ in range(2)]}
        start = tuple(state)
        for draw in draws:
            coin = int(draw < p.port_weights[left]); ones += coin
            total = state[left]+state[right]
            next_left = (total+coin)//2
            current = state[left]-next_left
            raw = {0: 0, 1: 1, -1: 2}[current]
            law = ordinal % 2
            if law:
                raw += 4*(total % 2)
            raw_counts[('C', 'E')[law]][coin][raw] += 1
            slot = 480*law+16*ei+8*coin+raw
            tape.extend(slot.to_bytes(2, 'little'))
            state[left], state[right] = next_left, (total+1-coin)//2
            ordinal += 1
        counts[left, right] = (n, ones)
        rows.append((start, tuple(state), raw_counts))
    return hashlib.sha256(tape).hexdigest(), tuple(state), counts, rows


@pytest.fixture(scope='module')
def tables():
    return a.transition_tables()


def test_every_native_table_column_matches_independent_endpoint_and_full_raw_menu(tables):
    assert tables['next'].shape == tables['raw'].shape == (2, 30, 2, 78)
    assert tables['certificate']['checked_basis_columns'] == 9360
    assert tables['certificate']['complete_raw_slots'] == 720
    states = tables['states']
    for li, law in enumerate(('C', 'E')):
        for ei, (left, right) in enumerate(tables['seams']):
            for coin in (0, 1):
                bins = [[] for _ in range(4 if law == 'C' else 8)]
                for i, state in enumerate(states):
                    final = list(state); total = state[left]+state[right]
                    final[left], final[right] = (total+coin)//2, (total+1-coin)//2
                    current = state[left]-final[left]
                    raw = (0 if current == 0 else 1 if current == 1 else 2)
                    if law == 'E': raw += 4*(total % 2)
                    assert states[int(tables['next'][li, ei, coin, i])] == tuple(final)
                    assert tables['raw'][li, ei, coin, i] == raw
                    bins[raw].append(tuple(final))
                assert [len(b) for b in bins] == tables['certificate']['raw_domain_sizes_H78'][law][ei][coin]
                assert sum(map(len, bins)) == 78
                assert all(len(set(b)) == len(b) for b in bins)
                assert bins[3] == []
                if law == 'E': assert bins[7] == []


def test_coin_effect_on_complex_entangled_input_is_input_independent(tables):
    rng = np.random.Generator(np.random.PCG64(7741))
    state = rng.normal(size=(78, 3))+1j*rng.normal(size=(78, 3))
    state /= np.linalg.norm(state)
    marginal = state.conj().T@state
    p = parameters()
    for li, law in enumerate(('C', 'E')):
        for ei, seam in enumerate(tables['seams']):
            for coin in (0, 1):
                probability = float(a.local.coin_probability(p, seam, coin))
                reference = np.zeros((3, 3), complex)
                for raw in range(4 if law == 'C' else 8):
                    chosen = np.flatnonzero(tables['raw'][li, ei, coin] == raw)
                    output = np.zeros_like(state)
                    output[tables['next'][li, ei, coin, chosen]] = np.sqrt(probability)*state[chosen]
                    reference += output.conj().T@output
                assert np.allclose(reference, probability*marginal, atol=2e-16, rtol=0)


def test_tables_and_certificate_have_no_public_cache_alias(tables):
    with pytest.raises(ValueError): tables['next'].flags.writeable = True
    with pytest.raises(ValueError): tables['raw'].flags.writeable = True
    first = a.transition_tables()
    first['certificate']['raw_domain_sizes_H78']['C'][0][0][0] = -1
    assert a.transition_tables()['certificate']['raw_domain_sizes_H78']['C'][0][0][0] >= 0


def test_native_and_source_maps_revalidated_after_previous_success(monkeypatch, tables):
    original = a.native.native
    def wrong(state, seam, coin):
        output, record = original(state, seam, coin)
        if seam == (0, 1) and coin == 0:
            return output, (coin, record[1]+1)
        return output, record
    monkeypatch.setattr(a.native, 'native', wrong)
    with pytest.raises(ArithmeticError): a.transition_tables()


def test_wrong_source_dilation_rejected_even_with_correct_endpoint_function(monkeypatch):
    swaps = dict(a.circuit.REPAIR_SWAPS)
    swaps['C', 0] = ((2, 10), (3, 20), (5, 12))
    monkeypatch.setattr(a.circuit, 'REPAIR_SWAPS', swaps)
    with pytest.raises(ArithmeticError, match='source record'): a.transition_tables()


def test_every_small_training_microevent_matches_independent_serial_execution():
    p = parameters(); initial = (0, 0, 1, 1)+(0,)*8
    expected_hash, terminal, expected_counts, expected_blocks = scalar_oracle(p, 37, 1843, initial)
    report, counts = a.acquire(p, samples_per_seam=37, seed=1843, initial_state=initial, chunk_size=11)
    assert report['stream_sha256'] == expected_hash
    assert tuple(report['private_replay']['terminal_state']) == terminal
    assert counts == expected_counts
    assert report['private_replay']['terminal_state_index'] == a.native.states(binary=True).index(terminal)
    for actual, (start, end, raw_counts) in zip(report['blocks'], expected_blocks, strict=True):
        assert a.native.states()[actual['start_state_index_H78']] == start
        assert a.native.states()[actual['terminal_state_index_H78']] == end
        assert actual['raw_counts'] == raw_counts
    json.dumps(report, allow_nan=False)


def test_chunk_boundaries_do_not_change_any_tape_event_or_coin_count():
    p = parameters()
    first, counts = a.acquire(p, samples_per_seam=33, seed=99, chunk_size=1)
    second, other = a.acquire(p, samples_per_seam=33, seed=99, chunk_size=32)
    assert counts == other
    assert first['blocks'] == second['blocks']
    assert first['stream_sha256'] == second['stream_sha256']
    assert first['private_replay'] == second['private_replay']
    assert first['resources']['host_vectorized_rng_calls'] != second['resources']['host_vectorized_rng_calls']


def test_coin_counts_do_not_use_hidden_source_state_but_raw_tape_does():
    p = parameters()
    left, c = a.acquire(p, samples_per_seam=11, seed=9, initial_state=(1, 1)+(0,)*10)
    right, d = a.acquire(p, samples_per_seam=11, seed=9, initial_state=(0, 0, 1, 1)+(0,)*8)
    assert c == d
    assert left['stream_sha256'] != right['stream_sha256']
    assert all(type(edge) is tuple and len(edge) == 2 and all(type(n) is int for n in row)
               for edge, row in c.items())


def test_same_carrier_chain_and_unrounded_realized_counts():
    report, counts = a.acquire(parameters(), samples_per_seam=101, seed=3)
    blocks = report['blocks']
    assert all(x['terminal_state_index_H78'] == y['start_state_index_H78']
               for x, y in zip(blocks, blocks[1:]))
    assert len({b['terminal_state_index_H78'] for b in blocks}) > 2
    assert any(count != round(n*float(a.local.coin_probability(parameters(), edge, 1)))
               for edge, (n, count) in counts.items())
    for block in blocks:
        assert sum(sum(row) for law in block['raw_counts'].values() for row in law) == 101
        assert all(row[3] == 0 for row in block['raw_counts']['C'])
        assert all(row[3] == row[7] == 0 for row in block['raw_counts']['E'])


def test_resources_cover_complete_archive_and_literal_source_counts():
    n = 13; row = a.resource_counts(n, chunk_size=5)
    assert row['native_endpoint_invocations'] == row['ideal_biased_coin_invocations'] == 390
    assert row['source_resets_during_training'] == 0
    assert row['FIFO_bits']+row['exported_record_bits'] == row['total_record_bits'] == 3900
    assert row['event_counter_bits'] == (390).bit_length()
    assert row['all30_coin_counter_bits'] == 30*(13).bit_length()
    assert row['host_vectorized_rng_calls'] == 90
    costs = row['compiled_word_costs']
    assert costs['native_pulses'] == sum(
        r['C_events']*r['compiled_C_pulses_per_event']+r['E_events']*r['compiled_E_pulses_per_event']
        for r in costs['per_seam'])
    assert costs['proper_code_load_calls'] == costs['proper_code_export_calls'] == costs['native_pulses']
    assert all(r['compiled_E_pulses_per_event'] == r['compiled_C_pulses_per_event']+3 for r in costs['per_seam'])


def test_requested_native_word_costs_match_independently_emitted_routes():
    requested = ((0, 1), (0, 2), (1, 7), (5, 9))
    actual = a.native_word_costs(requested)
    assert tuple(actual) == requested
    for edge in requested:
        # Independent permutation by simultaneous port assignment, rather
        # than calling the producer's route function or swapping load lists.
        ports = list(range(12))
        ports[0], ports[edge[0]] = ports[edge[0]], ports[0]
        ports[1], ports[edge[1]] = ports[edge[1]], ports[1]
        images = list(range(512))
        for state in a.native.states():
            images[a.circuit.encode(state)] = a.circuit.encode(tuple(state[p] for p in ports))
        gray = a.routing.permutation_edges(tuple(images))
        macro = a.routing._canonical_word()
        route = sum(len(macro)+2*pattern.bit_count() for _, pattern in gray)
        for law in ('C', 'E'):
            assert actual[edge][law] == 2*route+1+len(a.circuit.repair_word(law, 0))
    assert actual[0, 1]['C'] == 1+len(a.circuit.repair_word('C', 0))
    actual[0, 1]['C'] = -1
    assert a.native_word_costs(requested)[0, 1]['C'] > 0


@pytest.mark.parametrize('edges', [[], (), [(0, 1)], ((0, 1), (0, 1)),
                                  ((False, 1),), ((0, 0),), ((1, 0),),
                                  ((0, 12),), ((0., 1),)])
def test_requested_source_cost_pair_domain(edges):
    with pytest.raises(ValueError): a.native_word_costs(edges)


def test_source_cost_rejects_wrong_route_even_if_it_is_a_permutation(monkeypatch):
    original = a.record_code._route_permutation
    def wrong(edge):
        images = list(original(edge))
        images[510], images[511] = images[511], images[510]
        return tuple(images)
    monkeypatch.setattr(a.record_code, '_route_permutation', wrong)
    with pytest.raises(ArithmeticError, match='padding'):
        a.native_word_costs(((0, 1),))


@pytest.mark.parametrize('key,value', [('samples_per_seam', True), ('samples_per_seam', 0),
    ('samples_per_seam', (1 << 20)+1), ('samples_per_seam', 3.0), ('seed', -1), ('seed', True),
    ('seed', 1 << 64), ('chunk_size', 0), ('chunk_size', False), ('chunk_size', float('nan'))])
def test_original_scalar_domains_rejected(key, value):
    with pytest.raises(ValueError): a.acquire(parameters(), **{key: value})


@pytest.mark.parametrize('value', [(True, 1)+(0,)*10, [1, 1]+[0]*10,
                                 (2,)+(0,)*11, (1,)+(0,)*11])
def test_invalid_initial_source_rejected(value):
    with pytest.raises(ValueError): a.acquire(parameters(), samples_per_seam=1, initial_state=value)


def test_mutated_fraction_and_parameter_slots_rejected_before_sampling():
    p = parameters(); p.c_weight._denominator = -2
    with pytest.raises(ValueError): a.acquire(p, samples_per_seam=1)
    p = parameters(); object.__setattr__(p, 'port_weights', (True,)+(8,)*11)
    with pytest.raises(ValueError): a.acquire(p, samples_per_seam=1)


def test_returned_report_and_counts_are_detached():
    first, counts = a.acquire(parameters(), samples_per_seam=1)
    counts[next(iter(counts))] = (-1, -1)
    first['blocks'][0]['raw_counts']['C'][0][0] = -4
    second, later = a.acquire(parameters(), samples_per_seam=1)
    assert all(row[0] == 1 for row in later.values())
    assert second['blocks'][0]['raw_counts']['C'][0][0] >= 0


def test_full_production_horizon_has_all_actual_tape_events():
    report, counts = a.acquire(parameters())
    total = 30*(1 << 20)
    assert report['protocol']['total_events'] == total
    assert all(row[0] == 1 << 20 for row in counts.values())
    assert sum(b['trials'] for b in report['blocks']) == total
    assert sum(sum(sum(row) for row in raw) for b in report['blocks']
               for raw in b['raw_counts'].values()) == total
    assert report['resources']['total_record_bits'] == 10*total
    assert report['resources']['exported_records'] == total-1
    assert report['protocol']['seed_supplies_statistical_guarantee'] is False

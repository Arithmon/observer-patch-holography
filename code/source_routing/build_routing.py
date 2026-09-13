"""Canonical mean-only path calibration with destination-local readback."""
from __future__ import annotations

import hashlib
import json
from collections import deque
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def graph(support):
    adj = [set() for _ in range(support['carriers'] * 12)]
    for c in range(support['carriers']):
        for a, b in support['intra_carrier_seams']:
            adj[12*c+a].add(12*c+b)
            adj[12*c+b].add(12*c+a)
    for c, a, e, b in support['glued_pairs']:
        adj[12*c+a].add(12*e+b)
        adj[12*e+b].add(12*c+a)
    return adj


def paths(support, distances):
    adj = graph(support)
    c, a, _, _ = support['glued_pairs'][0]
    start = 12*c+a
    previous = {start: None}
    distance = {start: 0}
    queue = deque([start])
    selected = {}
    while queue:
        u = queue.popleft()
        if distance[u] in distances and distance[u] not in selected:
            route = [u]
            while previous[route[-1]] is not None:
                route.append(previous[route[-1]])
            selected[distance[u]] = list(reversed(route))
        for v in sorted(adj[u]):
            if v not in previous:
                previous[v] = u
                distance[v] = distance[u] + 1
                queue.append(v)
    if set(selected) != set(distances):
        raise ValueError('the supplied support does not realize every requested path')
    return [selected[d] for d in distances]


def decode_local(samples):
    """Only destination observations, chronological; no trace/pre-state input."""
    original = []
    state = []
    for final in samples:
        source = final
        for baseline in reversed(state):
            source = 2*source - baseline
        original.insert(0, source)
        next_state = []
        carry = source
        for baseline in state:
            carry = (carry + baseline)/2
            next_state.append(carry)
        next_state.append(carry)
        state = next_state
    return original, state


def rational_bits(value):
    return abs(value.numerator).bit_length() + value.denominator.bit_length() + 1


def episode(route, variant):
    n = len(route)
    d = n - 1
    initial = [Fraction(1000 + ((7919*r + 13) % 100003)) for r in route]
    if variant == 'source_plus_one':
        initial[0] += 1
    if variant == 'nonsource_minus_one':
        initial[max(1, d//2)] -= 1
    events, state, observed = [], {}, []
    previous_hash = '0'*64

    def event(op, reads, writes, **extra):
        nonlocal previous_hash
        row = {'id': len(events), 'previous_hash': previous_hash, 'op': op,
               'reads': [dict(register=r, **state[r]) for r in reads],
               'writes': [], **extra}
        for r, value in writes:
            version = state[r]['version'] + 1 if r in state else 0
            cell = {'version': version, 'writer': row['id'], 'value': value}
            state[r] = cell
            row['writes'].append(dict(register=r, **cell))
        row['event_hash'] = digest(row)
        previous_hash = row['event_hash']
        events.append(row)
        return row

    for r, value in zip(route, initial):
        event('initialize', [], [(f'p{r}', str(value))])

    def observe():
        k = len(observed)
        value = state[f'p{route[-1]}']['value']
        event('local_readback', [f'p{route[-1]}'], [(f's{k}', value)])
        observed.append(Fraction(value))

    observe()
    for start in range(d-1, -1, -1):
        for k in range(start, d):
            a, b = (f'p{route[k]}', f'p{route[k+1]}')
            x, y = Fraction(state[a]['value']), Fraction(state[b]['value'])
            mean = (x+y)/2
            event('pair_mean', [a, b], [(a, str(mean)), (b, str(mean))],
                  active=x != y, quadratic_decrement=str((x-y)**2/2))
        observe()
    recovered, final = decode_local(observed)
    assert recovered == initial
    assert final == [Fraction(state[f'p{r}']['value']) for r in route]
    event('decode', [f's{k}' for k in range(n)], [('decoded', [str(v) for v in recovered])])
    inverse_rows = [[] for _ in range(n)]
    for j in range(n):
        result, _ = decode_local([Fraction(int(i == j)) for i in range(n)])
        for i, value in enumerate(result):
            inverse_rows[i].append(value)
    source_gain = sum(abs(v) for v in inverse_rows[0])
    rational_values = []
    for row in events:
        for write in row['writes']:
            values = write['value'] if isinstance(write['value'], list) else [write['value']]
            rational_values.extend(Fraction(v) for v in values)
    means = [e for e in events if e['op'] == 'pair_mean']
    return {
        'path': route, 'variant': variant, 'events': events,
        'final_event_hash': previous_hash,
        'local_samples': [str(v) for v in observed],
        'source_inverse_coefficients': [str(v) for v in inverse_rows[0]],
        'cost': {
            'path_edges': d, 'initial_port_writes': n,
            'seam_mean_operations': len(means),
            'active_seam_means': sum(e['active'] for e in means),
            'calibration_seam_operations': d*(d-1)//2,
            'final_transfer_seam_operations': d,
            'seam_reads': 2*len(means), 'seam_writes': 2*len(means),
            'local_readbacks': n, 'protected_sample_words': n,
            'decoded_output_words': n,
            'all_events': len(events),
            'all_register_reads': sum(len(e['reads']) for e in events),
            'all_register_write_records': sum(len(e['writes']) for e in events),
            'max_rational_storage_bits': max(map(rational_bits, rational_values)),
            'written_rational_payload_bits': sum(map(rational_bits, rational_values)),
            'protected_sample_payload_bits': sum(map(rational_bits, observed)),
        },
        'conditioning': {
            'final_sample_source_response_to_plus_one': str(Fraction(1, 2**d)),
            'source_error_per_uniform_sample_error': str(source_gain),
            'meaning': 'For exact mean dynamics and independently perturbed local samples of absolute error at most epsilon, the decoded initial source error is at most this l1 coefficient norm times epsilon. No arithmetic or dynamical error is included.',
        },
    }


def make():
    specification = json.loads((HERE/'specification.json').read_text())
    support = json.loads((HERE/'support_w12_l3.json').read_text())
    selected = paths(support, specification['path_lengths'])
    episodes = [episode(route, variant) for route in selected for variant in specification['variants']]
    refs = ['specification.json', 'support_w12_l3.json', 'build_routing.py']
    return {
        'schema': 'oph.source_routing.path_tomography.v1',
        'pins': {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in refs},
        'support': {'level': 3, 'carriers': support['carriers'], 'glued_pairs': len(support['glued_pairs'])},
        'episodes': episodes,
        'scope': {
            'canonical_scalar_pair_means_only': True,
            'decoder_uses_destination_local_samples_only': True,
            'path_and_schedule_declared': True,
            'initial_loads_declared': True,
            'protected_history_declared': True,
            'all_intermediate_events_retained': True,
            'hash_chain_is_custody_not_causal_edges': True,
            'full_metric_neighbor_refinement': False,
            'q13_q21_complete_routing': False,
            'source_population_produced': False,
            'physical_capacity_or_clock_identified': False,
            'physical_premise_discharged': False,
        },
    }


if __name__ == '__main__':
    packet = make()
    (HERE/'runtime/path_tomography_receipt.json').write_bytes(canonical(packet))
    for item in packet['episodes']:
        if item['variant'] == 'baseline':
            print(item['cost'], item['conditioning'])

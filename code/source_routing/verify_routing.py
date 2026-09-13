"""Independent exact matrix and authenticated-event verifier; no producer import.

The decoder is checked by constructing the local observation matrix directly
from the specified canonical seam operations and inverting it by elimination.
All graph, value, custody, count and conditioning assertions are recomputed.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def raw(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def reject(value):
        raise ValueError('non-exact JSON number: '+value)
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_float=reject, parse_constant=reject)


def invert(matrix):
    n = len(matrix)
    a = [list(row) + [F(int(i == j)) for j in range(n)] for i, row in enumerate(matrix)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        require(pivot is not None, 'singular local observation map')
        a[k], a[pivot] = a[pivot], a[k]
        scale = a[k][k]
        a[k] = [x/scale for x in a[k]]
        for i in range(n):
            if i != k:
                scale = a[i][k]
                a[i] = [x-scale*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


def observation_matrix(d):
    n = d+1
    state = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    observations = [state[-1][:]]
    for width in range(1, n):
        for k in range(n-1-width, n-1):
            mean = [(a+b)/2 for a, b in zip(state[k], state[k+1])]
            state[k], state[k+1] = mean[:], mean[:]
        observations.append(state[-1][:])
    return observations, state


def verify_support(support):
    require(support['schema'] == 'oph.source_routing.w12_support.v1', 'support schema')
    require(support['level'] == 3 and support['carriers'] == 1280, 'support size')
    faces = support['faces']
    require(len(faces) == 1280 and all(len(set(f)) == 3 for f in faces), 'mesh faces')
    require(len({tuple(sorted(f)) for f in faces}) == 1280, 'duplicate face')
    by_vertex = defaultdict(list)
    for cell, face in enumerate(faces):
        for v in face:
            by_vertex[v].append(cell)
    require(len(by_vertex) == 642, 'level-three vertex count')
    require(Counter(map(len, by_vertex.values())) == Counter({6:630, 5:12}), 'geodesic vertex valences')
    expected_pairs = set()
    for cells in by_vertex.values():
        for i, a in enumerate(cells):
            for b in cells[i+1:]:
                expected_pairs.add(tuple(sorted((a, b))))
    actual_pairs, slots, edges = set(), set(), set()
    for c, a, e, b in support['glued_pairs']:
        require(0 <= c < e < 1280 and 0 <= a < 12 and 0 <= b < 12, 'glued endpoint domain')
        require((c,e) not in actual_pairs, 'duplicate glued pair')
        require((c,a) not in slots and (e,b) not in slots, 'glued slot reused')
        actual_pairs.add((c,e)); slots.update(((c,a),(e,b)))
        edges.add(tuple(sorted((12*c+a,12*e+b))))
    require(actual_pairs == expected_pairs and len(actual_pairs) == 7650, 'full shared-vertex support')
    seams = {tuple(x) for x in support['intra_carrier_seams']}
    # The committed twelve-port carrier's fixed combinatorial icosahedron.
    carrier_faces = [(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
                     (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
                     (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
                     (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    expected_seams = {tuple(sorted((face[i],face[(i+1)%3]))) for face in carrier_faces for i in range(3)}
    require(seams == expected_seams, 'canonical intra-carrier seams')
    for c in range(1280):
        edges.update((12*c+a,12*c+b) for a,b in seams)
    adj = defaultdict(set)
    for a,b in edges:
        adj[a].add(b); adj[b].add(a)
    visited, todo = set(), [0]
    while todo:
        u = todo.pop()
        if u not in visited:
            visited.add(u); todo.extend(adj[u]-visited)
    require(len(visited) == 15360, 'connected port support')
    return edges


def bits(x):
    return x.numerator.bit_length() + x.denominator.bit_length() + 1


def verify_episode(ep, edges):
    route = ep['path']; n = len(route); d = n-1
    require(len(set(route)) == n and n >= 2, 'simple nontrivial path')
    require(all(tuple(sorted(pair)) in edges for pair in zip(route, route[1:])), 'path uses an absent seam')
    variant = ep['variant']
    initial = [F(1000+((7919*r+13)%100003)) for r in route]
    require(variant in ('baseline','source_plus_one','nonsource_minus_one'), 'variant')
    if variant == 'source_plus_one': initial[0] += 1
    if variant == 'nonsource_minus_one': initial[max(1, d//2)] -= 1
    observation, final_matrix = observation_matrix(d)
    inverse = invert(observation)
    # Independently verify both matrix products, not just the elimination trace.
    for a,b in ((inverse,observation),(observation,inverse)):
        require(all(sum(a[i][k]*b[k][j] for k in range(n)) == int(i==j)
                    for i in range(n) for j in range(n)), 'inverse identity')
    schedule = [('initialize', [f'p{r}']) for r in route]
    schedule.append(('local_readback', ['s0']))
    for phase in range(1,n):
        schedule.extend(('pair_mean', [f'p{route[k]}',f'p{route[k+1]}']) for k in range(d-phase,d))
        schedule.append(('local_readback',[f's{phase}']))
    schedule.append(('decode',['decoded']))
    require(len(ep['events']) == len(schedule), 'event loss or insertion')
    state, sample_values, rational_values = {}, [], []
    previous_hash = '0'*64
    ancestry = []
    for eid,(event, (op, registers)) in enumerate(zip(ep['events'], schedule)):
        require(type(event['id']) is int and event['id'] == eid and event['previous_hash'] == previous_hash, 'event chronology')
        unhashed = {k:v for k,v in event.items() if k != 'event_hash'}
        require(event['event_hash'] == hashlib.sha256(raw(unhashed)).hexdigest(), 'event hash')
        previous_hash = event['event_hash']
        require(event['op'] == op, 'declared calibration schedule')
        mask = 0
        for read in event['reads']:
            require(read['register'] in state and {k:v for k,v in read.items() if k!='register'} == state[read['register']], 'version, writer or value mismatch')
            require(read['writer'] < eid, 'future read')
            mask |= ancestry[read['writer']] | (1 << read['writer'])
        ancestry.append(mask)
        if op == 'initialize':
            require(event['reads'] == [], 'preparation dependency')
            expected = [str(initial[eid])]
        elif op == 'pair_mean':
            require([r['register'] for r in event['reads']] == registers, 'seam read footprint')
            x,y = (F(r['value']) for r in event['reads'])
            expected = [str((x+y)/2)]*2
            require(event['active'] == (x!=y), 'active repair count')
            require(F(event['quadratic_decrement']) == (x-y)**2/2, 'exact quadratic decrement')
            require(x*x+y*y-2*((x+y)/2)**2 == F(event['quadratic_decrement']), 'energy identity')
        elif op == 'local_readback':
            require([r['register'] for r in event['reads']] == [f'p{route[-1]}'], 'nonlocal decoder side-channel')
            expected = [state[f'p{route[-1]}']['value']]
            sample_values.append(F(expected[0]))
        else:
            require([r['register'] for r in event['reads']] == [f's{k}' for k in range(n)], 'decoder access firewall')
            recovered = [sum(row[j]*sample_values[j] for j in range(n)) for row in inverse]
            require(recovered == initial, 'initial records not recovered')
            expected = [[str(v) for v in recovered]]
        require([r['register'] for r in event['writes']] == registers, 'write footprint')
        require([r['value'] for r in event['writes']] == expected, 'wrong law or readback')
        for write in event['writes']:
            r=write['register']; version=state[r]['version']+1 if r in state else 0
            require(type(write['writer']) is int and type(write['version']) is int and write['writer'] == eid and write['version'] == version, 'write custody')
            state[r]={k:v for k,v in write.items() if k!='register'}
            values=write['value'] if isinstance(write['value'],list) else [write['value']]
            rational_values.extend(F(v) for v in values)
    require(ep['final_event_hash'] == previous_hash, 'final hash')
    require([F(x) for x in ep['local_samples']] == sample_values, 'sample mirror')
    require(sample_values == [sum(row[j]*initial[j] for j in range(n)) for row in observation], 'local observation matrix')
    require([F(state[f'p{r}']['value']) for r in route] == [sum(row[j]*initial[j] for j in range(n)) for row in final_matrix], 'changed relay state')
    # Every original path preparation is an ancestor of the decoder; every
    # transport/calibration event remains an authenticated event in this order.
    require(ancestry[-1] == (1 << (len(ep['events'])-1))-1,
            'an initialization, repair or readback is absent from final ancestry')
    coefficients = inverse[0]
    require([F(v) for v in ep['source_inverse_coefficients']] == coefficients, 'inverse coefficients')
    gain=sum(abs(v) for v in coefficients)
    require(F(ep['conditioning']['source_error_per_uniform_sample_error']) == gain, 'sample error amplification')
    require(F(ep['conditioning']['final_sample_source_response_to_plus_one']) == observation[-1][0] == F(1,2**d), 'source attenuation')
    # The stated uniform-error norm is sharp: the signed corner attains it.
    require(sum(c*(1 if c>0 else -1 if c<0 else 0) for c in coefficients) == gain, 'noise corner')
    means=[e for e in ep['events'] if e['op']=='pair_mean']; t=d*(d+1)//2
    expected_cost={
        'path_edges':d,'initial_port_writes':n,'seam_mean_operations':t,
        'active_seam_means':sum(e['active'] for e in means),
        'calibration_seam_operations':d*(d-1)//2,'final_transfer_seam_operations':d,
        'seam_reads':2*t,'seam_writes':2*t,'local_readbacks':n,
        'protected_sample_words':n,'decoded_output_words':n,'all_events':t+2*n+1,
        'all_register_reads':2*t+2*n,'all_register_write_records':2*t+2*n+1,
        'max_rational_storage_bits':max(map(bits,rational_values)),
        'written_rational_payload_bits':sum(map(bits,rational_values)),
        'protected_sample_payload_bits':sum(map(bits,sample_values)),
    }
    require(raw(ep['cost'])==raw(expected_cost), 'cost or capacity accounting')
    return {'path_edges':d,'variant':variant,'seam_means':t,'events':len(ep['events']),
            'source_noise_gain':str(gain),'source_samples':sample_values,
            'original':initial,'provenance_edges':sum(len({r['writer'] for r in e['reads']}) for e in ep['events'])}


def verify(packet=None, support=None, check_pins=True):
    if packet is None: packet=load(HERE/'runtime/path_tomography_receipt.json')
    if support is None: support=load(HERE/'support_w12_l3.json')
    spec=load(HERE/'specification.json')
    require(packet['schema']=='oph.source_routing.path_tomography.v1','receipt schema')
    if check_pins:
        require(set(packet['pins'])=={'specification.json','support_w12_l3.json','build_routing.py'},'pin coverage')
        for p,h in packet['pins'].items(): require(hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h,'pin mismatch')
    require(raw(packet['scope'])==raw({
        'canonical_scalar_pair_means_only':True,'decoder_uses_destination_local_samples_only':True,
        'path_and_schedule_declared':True,'initial_loads_declared':True,'protected_history_declared':True,
        'all_intermediate_events_retained':True,'hash_chain_is_custody_not_causal_edges':True,
        'full_metric_neighbor_refinement':False,'q13_q21_complete_routing':False,
        'source_population_produced':False,'physical_capacity_or_clock_identified':False,
        'physical_premise_discharged':False}),'scope promotion')
    require(packet['support']=={'level':3,'carriers':1280,'glued_pairs':7650},'support summary')
    expected=[(d,v) for d in spec['path_lengths'] for v in spec['variants']]
    require([(len(e['path'])-1,e['variant']) for e in packet['episodes']]==expected,'missing scheduled outcome')
    edges=verify_support(support)
    results=[verify_episode(e,edges) for e in packet['episodes']]
    for k in range(0,len(results),3):
        baseline,source,relay=results[k:k+3]
        require(all(e['path']==packet['episodes'][k]['path'] for e in packet['episodes'][k:k+3]),'intervention changes path')
        require(source['source_samples'][:-1]==baseline['source_samples'][:-1],'source enters before its seam read')
        require(source['source_samples'][-1]-baseline['source_samples'][-1]==F(1,2**baseline['path_edges']),'intervention response')
        require(source['original'][0]-baseline['original'][0]==1 and source['original'][1:]==baseline['original'][1:],'source intervention decoding')
        require(sum(a!=b for a,b in zip(relay['original'],baseline['original']))==1,'relay intervention decoding')
    return [{k:v for k,v in row.items() if k not in ('source_samples','original')} for row in results]


if __name__=='__main__':
    print(json.dumps({'verified':True,'episodes':verify()},indent=2))

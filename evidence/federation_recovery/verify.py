"""Independent finite-word and complete-gluing ambiguity checks.

The event verifier imports neither the recovery decoder nor the producer.
Geometry and native execution enter only through the outer certificate
boundary. The accepted counterexample concerns the declared finite gluing
class, not two canonical geometric subdivisions or all possible logs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import struct


NONE = (1 << 64) - 1
NAMES = ('ZERO', 'SCRATCH', 'PUBLISH', 'EXPORT', 'RESET', 'MEAN',
         'CAPTURE', 'START', 'ADD', 'COMMIT', 'ADDRESS')
ASSUMPTIONS = [
    'authenticated_complete_retained_event_rows',
    'declared_unsigned64_rows_signed_half_unit_values',
    'declared_M1_opcode_and_register_ownership_semantics',
    'candidate_class_connected_finite_twelve_port_gluings',
    'level_three_provenance_supplied_outside_recovery',
]
NONCLAIMS = [
    'no_canonical_geometric_alternative_asserted',
    'no_source_selection_of_feedback_law',
    'no_unseen_population_or_routing_law_recovery',
    'no_complete_gluing_fixed_point',
    'no_physical_clock_or_quantum_instrument',
    'no_universe_level_existence_or_uniqueness',
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True) + '\n').encode('ascii')


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def same(actual, declared, message):
    # JSON bytes distinguish False from 0, unlike Python container equality.
    require(canonical(actual) == canonical(declared), message)


def check_rows(rows):
    """Reconstruct the observed finite program from actual consumed writers.

    This checks the baseline register grammar. No support, radius menu,
    source coordinate generator or externally supplied population is read.
    The fixed receipt is checked separately against its authenticated bytes.
    """
    rows = list(rows)
    require(0 < len(rows) <= 100000, 'bounded event word')
    memory, ancestor_masks = {}, []
    ports, zero, accumulators, addresses, initial = {}, {}, {}, {}, {}
    ports_finished = False
    counts = Counter()
    reads = writes = maximum = 0
    versions, pending_reads, logical_reads, hops = [], {}, [], []
    site_by_owner, generation = {}, {}
    starts, commits, phase = set(), set(), 0
    active_hop = None
    gluing, program = set(), []
    partners, carrier_seams = {}, {}

    for eid, raw in enumerate(rows):
        require(isinstance(raw, (tuple, list)) and len(raw) == 8, 'event width')
        require(all(type(x) is int and 0 <= x <= NONE for x in raw[:7]),
                'unsigned event interface field')
        require(type(raw[7]) is int and -(1 << 63) <= raw[7] <= NONE,
                'signed scalar representation')
        op, a, b, out, owner, wa, wb, encoded = raw
        require(op < len(NAMES), 'unknown primitive')
        value = encoded - (1 << 64) if encoded >= (1 << 63) else encoded
        counts[NAMES[op]] += 1
        maximum = max(maximum, abs(value))
        operands, mask = [], 0
        for register, writer in ((a, wa), (b, wb)):
            if register == NONE:
                require(writer == NONE, 'writer on absent operand')
                operands.append(None)
            else:
                require(register in memory, 'absent input register')
                cell = memory[register]
                require(cell['writer'] == writer < eid, 'stale consumed writer')
                require(op == 5 or cell['owner'] == owner, 'nonlocal unary read')
                operands.append(cell)
                mask |= ancestor_masks[writer]
                reads += 1
        x, y = operands
        origin, path = None, None
        role, immutable = None, False
        preparation = op in (0, 1, 2, 10)
        if preparation:
            require(phase == 0 and a == b == NONE and out == len(memory),
                    'preparation interface')
            require(owner != NONE, 'preparation owner')
            if op == 1 and not ports_finished:
                ports.setdefault(owner, []).append(out)
                role = 'port'
            elif op == 0:
                ports_finished = True
                require(owner in ports and owner not in zero and value == 0,
                        'protected zero preparation')
                zero[owner] = out
                role, immutable = 'zero', True
            elif op == 1:
                require(owner in zero and owner not in accumulators and value == 0,
                        'accumulator preparation')
                accumulators[owner] = out
                role = 'accumulator'
            elif op == 10:
                require(owner in accumulators and value % 2 == 0,
                        'integer address preparation')
                addresses.setdefault(owner, []).append((out, value // 2))
                require(len(addresses[owner]) <= 6, 'address width')
                role, immutable = 'address', True
            else:
                require(owner in addresses and len(addresses[owner]) == 6
                        and owner not in initial, 'initial logical record')
                site = len(initial)
                require(value == 2 * (site + 1), 'baseline initial payload')
                site_by_owner[owner] = site
                generation[owner] = 0
                initial[owner] = out
                origin, path = len(versions), [owner]
                mask |= 1 << origin
                versions.append({'index': origin, 'site': site, 'layer': 0,
                                 'event': eid, 'register': out, 'scaled_value': value})
                role, immutable = 'archive', True
        else:
            require(bool(initial), 'dynamics without population')
            require(owner != NONE or op == 5, 'non-mean absent owner')
            require(op == 5 or owner in ports, 'unknown dynamic owner')
            if active_hop is not None:
                expected_op = (3, 4, 5, 6, 4, 4)[len(active_hop['events'])]
                require(op == expected_op, 'interrupted transport word')
            if op == 3:
                require(active_hop is None and x is not None and b == NONE
                        and x['role'] == 'archive' and x['origin'] is not None,
                        'export interface')
                require(out in memory and memory[out]['role'] == 'port', 'export destination')
                require(value == x['value'], 'incorrect export law')
                origin, path = x['origin'], list(x['path'])
                active_hop = {'events': [], 'source_version': origin,
                              'source_register': a, 'source_port': out}
            elif op == 4:
                require(active_hop is not None and x is not None and b == NONE
                        and x['role'] == 'zero' and x['value'] == value == 0,
                        'incorrect reset law')
                position = len(active_hop['events'])
                if position == 1:
                    active_hop['target_port'] = out
                else:
                    require(out == active_hop['source_port' if position == 4 else 'target_port'],
                            'transport cleanup destination')
                require(out in memory and memory[out]['role'] == 'port', 'reset destination')
            elif op == 5:
                require(active_hop is not None and x is not None and y is not None
                        and a != b and owner == NONE and out == a,
                        'mean interface')
                require(a == active_hop['source_port'] and b == active_hop['target_port'],
                        'transport mean endpoints')
                require(x['role'] == y['role'] == 'port' and x['owner'] != y['owner'],
                        'mean carrier ports')
                require(y['value'] == 0 and y['origin'] is None,
                        'mean receiver not reset')
                require(x['value'] + y['value'] == 2 * value, 'incorrect mean law')
                origin, path = x['origin'], list(x['path'])
                seam = tuple(sorted((a, b)))
                carriers = tuple(sorted((x['owner'], y['owner'])))
                require(partners.get(a, b) == b and partners.get(b, a) == a,
                        'observed port has two gluing partners')
                require(carrier_seams.get(carriers, seam) == seam,
                        'two seams between the same carrier pair')
                partners[a], partners[b] = b, a
                carrier_seams[carriers] = seam
                gluing.add(seam)
            elif op == 6:
                require(active_hop is not None and a == active_hop['target_port']
                        and x is not None and b == NONE and out == len(memory),
                        'capture interface')
                require(value == 2 * x['value'], 'incorrect capture law')
                origin, path = x['origin'], list(x['path']) + [owner]
                role, immutable = 'archive', True
                active_hop['capture_register'] = out
            elif op == 7:
                require(active_hop is None and x is not None and b == NONE
                        and x['role'] == 'zero' and x['value'] == 0 and value == 2
                        and out == accumulators.get(owner), 'incorrect start law')
                if not starts or commits == set(initial):
                    starts, commits = set(), set()
                    phase += 1
                require(owner not in starts and not commits, 'duplicate or late start')
                starts.add(owner)
                pending_reads[owner] = []
            elif op == 8:
                require(starts == set(initial) and not commits and x is not None
                        and y is not None and a == out == accumulators.get(owner)
                        and x['role'] == 'accumulator' and y['role'] == 'archive',
                        'logical addition interface')
                require(value == x['value'] + y['value'], 'incorrect addition law')
                source = y['origin']
                require(source is not None and versions[source]['layer'] == phase - 1,
                        'wrong logical source generation')
                require(y['path'][-1] == owner, 'read path endpoint')
                require(source not in [r['source'] for r in pending_reads[owner]],
                        'duplicate logical source read')
                pending_reads[owner].append({'event': eid, 'source': source,
                                             'via_register': b, 'path': list(y['path'])})
            elif op == 9:
                require(starts == set(initial) and owner not in commits
                        and x is not None and b == NONE and a == accumulators.get(owner)
                        and out == len(memory), 'logical commit interface')
                require(value == x['value'], 'baseline commit law')
                require(generation[owner] == phase - 1, 'commit generation')
                generation[owner] = phase
                origin, path = len(versions), [owner]
                mask |= 1 << origin
                versions.append({'index': origin, 'site': site_by_owner[owner],
                                 'layer': phase, 'event': eid, 'register': out,
                                 'scaled_value': value})
                for row in pending_reads[owner]:
                    logical_reads.append(dict(row, target=origin))
                commits.add(owner)
                role, immutable = 'archive', True
            else:
                raise ValueError('unsupported dynamic operation')

        ancestor_masks.append(mask)
        require(len(versions) <= 512 and len(memory) <= 30000, 'finite recovery resource bound')
        if op == 5:
            targets = (a, b)
        else:
            targets = (out,)
        for register in targets:
            if register in memory:
                old = memory[register]
                require(not old['immutable'], 'protected register overwrite')
                require(op == 5 or old['owner'] == owner, 'destination owner changed')
                next_role, next_owner = old['role'], old['owner']
            else:
                require(register == len(memory) and role is not None, 'noncontiguous allocation')
                next_role, next_owner = role, owner
            memory[register] = {'writer': eid, 'value': value, 'owner': next_owner,
                                'role': next_role, 'immutable': immutable,
                                'origin': origin, 'path': path}
            writes += 1
        if active_hop is not None:
            active_hop['events'].append(eid)
            if len(active_hop['events']) == 6:
                hops.append(active_hop)
                active_hop = None
        program.append([op, a, b, out, owner, value if preparation else None])

    require(active_hop is None and phase > 0 and commits == set(initial),
            'incomplete finite program')
    require(sorted(ports) == list(range(len(ports))) and set(zero) == set(ports),
            'carrier inventory')
    require(all(len(v) == 12 for v in ports.values()), 'twelve-port inventory')
    require(set(initial) == set(addresses) == set(accumulators), 'population inventory')
    population = []
    for owner, site in sorted(site_by_owner.items(), key=lambda pair: pair[1]):
        pairs = addresses[owner]
        require(len(pairs) == 6, 'six address coefficients')
        population.append({'site': site, 'owner': owner, 'accumulator': accumulators[owner],
                           'address_registers': [r for r, _ in pairs],
                           'address_coefficients': [v for _, v in pairs],
                           'initial_payload_register': initial[owner]})
    # Check the declared positive/negative prepared-port interpretation by
    # looking at the actual initial writer rows, not final mutable contents.
    for owner, registers in ports.items():
        coefficients = [v for _, v in addresses.get(owner, [])]
        expected = [0] * 12
        for j, (plus, minus) in enumerate(zip((0, 1, 4, 5, 8, 9), (3, 2, 7, 6, 11, 10))):
            if coefficients:
                expected[plus] = 2 * max(coefficients[j], 0)
                expected[minus] = 2 * max(-coefficients[j], 0)
        for register, value in zip(registers, expected):
            encoded = rows[register][7]
            actual = encoded - (1 << 64) if encoded >= (1 << 63) else encoded
            require(rows[register][3] == register and actual == value,
                    'baseline address port split')
    for version in versions:
        mask = ancestor_masks[version['event']]
        version['ancestors'] = [i for i in range(len(versions)) if mask & (1 << i)]
    logical_reads.sort(key=lambda row: row['event'])
    edges = sorted({(r['source'], r['target']) for r in logical_reads})
    menus = {}
    for source, target in edges:
        key = (versions[target]['layer'], versions[target]['site'])
        menus.setdefault(key, set()).add(versions[source]['site'])
    require(all(menus.get((layer, site)) == menus.get((1, site))
                for layer in range(1, phase + 1) for site in range(len(initial))),
            'changing observed layer menu')
    require(all(site in menus[(1, site)] for site in range(len(initial))),
            'missing logical self read')
    ancestors = [set(row['ancestors']) for row in versions]
    strict = [(a, b) for b, chain in enumerate(ancestors) for a in sorted(chain) if a != b]
    histogram = Counter(sum(a in ancestors[x] for x in ancestors[b] - {a, b}) for a, b in strict)
    heights = []
    for i, chain in enumerate(ancestors):
        heights.append(1 + max((heights[a] for a in chain if a != i), default=0))
    layer_sizes = [sum(v['layer'] == layer for v in versions) for layer in range(phase + 1)]
    return {
        'schema': 'oph.federation_recovery.observed.v1',
        'grammar': {'cohort': 'baseline', 'value_encoding': 'signed_int64_half_units',
                    'ports_per_owner': 12, 'address_coefficients_per_site': 6,
                    'primitive_semantics': 'supplied finite register grammar'},
        'scope': {'event_only': True, 'observed_routing_recovered': True,
                  'population_records_recovered': True, 'logical_poset_recovered': True,
                  'full_port_gluing_recovered': False, 'general_population_law_recovered': False,
                  'general_routing_algorithm_recovered': False, 'physical_identification': False},
        'census': {'events': len(rows), 'registers': len(memory), 'scalar_reads': reads,
                   'scalar_writes': writes, 'operation_counts': dict(counts),
                   'max_abs_scaled_scalar': maximum},
        'carrier_ports': [ports[owner] for owner in range(len(ports))],
        'population': population, 'observed_gluing': [list(edge) for edge in sorted(gluing)],
        'transport_hops': hops, 'logical_versions': versions, 'logical_reads': logical_reads,
        'logical_edges': [list(edge) for edge in edges],
        'poset': {'nodes': len(versions), 'strict_relations': len(strict),
                  'reflexive_relations': len(strict) + len(versions), 'layer_sizes': layer_sizes,
                  'max_chain_length': max(heights),
                  'open_interval_cardinality_histogram': {str(k): v for k, v in sorted(histogram.items())}},
        'invariants': {'owners': len(ports), 'ports': 12 * len(ports),
                       'population_sites': len(initial), 'rounds': phase,
                       'observed_glued_seams': len(gluing), 'logical_nodes': len(versions),
                       'logical_reads': len(logical_reads), 'logical_strict_relations': len(strict),
                       'logical_height': max(heights)},
        'program': program,
    }


def alternative_edges(edges):
    """Apply the fixed witness switch independently of the topology builder."""
    replacement = {(120, 5, 123, 4): (120, 5, 135, 11),
                   (132, 9, 135, 11): (123, 4, 132, 9)}
    rows = [tuple(row) for row in edges]
    require(all(rows.count(edge) == 1 for edge in replacement), 'missing switch source')
    return [list(replacement.get(row, row)) for row in rows]


def _graph(edges):
    neighbours, occupied, pairs = {}, set(), set()
    rows = []
    for row in edges:
        require(isinstance(row, (tuple, list)) and len(row) == 4
                and all(type(x) is int for x in row), 'gluing row type')
        a, p, b, q = row
        require(0 <= a < b and 0 <= p < 12 and 0 <= q < 12, 'gluing domain')
        require((a, b) not in pairs and (a, p) not in occupied and (b, q) not in occupied,
                'gluing repeats an edge or port')
        pairs.add((a, b)); occupied.update(((a, p), (b, q)))
        neighbours.setdefault(a, []).append(b); neighbours.setdefault(b, []).append(a)
        rows.append(tuple(row))
    require(bool(neighbours) and sorted(neighbours) == list(range(len(neighbours))),
            'gluing owner inventory')
    adjacency = [sorted(neighbours[i]) for i in range(len(neighbours))]
    # Count a triangle once by its unique ascending sequence a < b < c.
    triangles = sum(1 for a, row in enumerate(adjacency) for b in row if a < b
                    for c in adjacency[b] if b < c and c in row)
    degrees = [len(row) for row in adjacency]
    summary = {'carriers': len(adjacency), 'glued_edges': len(rows), 'triangles': triangles,
               'degree_histogram': [list(x) for x in sorted(Counter(degrees).items())],
               'degree_vector_sha256': digest(degrees),
               'occupied_ports_sha256': digest(sorted(occupied)),
               'canonical_edges_sha256': digest(sorted(rows))}
    return summary, adjacency, occupied


def _search(graph, root):
    parents = [-1] * len(graph)
    parents[root] = root
    queue = [root]
    for vertex in queue:
        for neighbour in graph[vertex]:
            if parents[neighbour] < 0:
                parents[neighbour] = vertex
                queue.append(neighbour)
    return parents, queue


def check_completions(header, hosts, edges, alternate=None):
    require(tuple(header) == (3, 27, 1280, 2, 13, 7650, 311, 1)
            and all(type(x) is int for x in header), 'fixed input header')
    require(len(hosts) == 27 and all(type(x) is int for x in hosts)
            and sorted(hosts) == list(range(27)), 'fixed host inventory')
    expected_alternate = alternative_edges(edges)
    alternate = expected_alternate if alternate is None else alternate
    same(expected_alternate, [list(e) for e in alternate], 'undeclared topology switch')
    old, old_graph, old_ports = _graph(edges)
    new, new_graph, new_ports = _graph(alternate)
    require(old['carriers'] == new['carriers'] == 1280
            and old['glued_edges'] == new['glued_edges'] == 7650,
            'fixed complete graph census')
    require(old['triangles'] == 14000 and new['triangles'] == 13990,
            'triangle nonisomorphism witness')
    roots = sorted(hosts)
    first, second = hashlib.sha256(), hashlib.sha256()
    for root in roots:
        p, q = _search(old_graph, root)
        pp, qq = _search(new_graph, root)
        require(len(q) == len(qq) == 1280, 'disconnected completion')
        require(p == pp and q == qq, 'changed source-host BFS')
        first.update(canonical({'root': root, 'parents': p, 'queue': q}))
        second.update(canonical({'root': root, 'parents': pp, 'queue': qq}))
    require([len(x) for x in old_graph] == [len(x) for x in new_graph]
            and old_ports == new_ports, 'changed degree or port roster')
    return {'schema': 'oph.federation_recovery.topology_comparison.v1',
            'original': old, 'alternate': new, 'hosts': roots,
            'comparison': {key: True for key in (
                'carrier_counts_equal', 'glued_edge_counts_equal', 'per_carrier_degrees_equal',
                'occupied_carrier_ports_equal', 'original_connected', 'alternate_connected',
                'all_host_BFS_parent_maps_equal', 'all_host_BFS_queue_orders_equal',
                'nonisomorphic_by_triangle_count')},
            'original_BFS_sha256': first.hexdigest(), 'alternate_BFS_sha256': second.hexdigest(),
            'first_BFS_difference': None, 'witness_valid': True,
            'scope': 'Declared finite twelve-port gluing graphs; no event-log recovery, '
                     'native execution, or canonical geometric realization is certified here.'}


def execute_program(program):
    """Independently compute an observed program's values and writer IDs.

    Only preparation instructions carry literals.  The returned rows still
    require the full grammar, ownership and gluing checks in ``check_rows``.
    """
    require(isinstance(program, list) and 0 < len(program) <= 100000,
            'bounded instruction program')
    memory, rows = {}, []
    for eid, instruction in enumerate(program):
        require(isinstance(instruction, list) and len(instruction) == 6,
                'instruction shape')
        require(all(type(x) is int and 0 <= x <= NONE for x in instruction[:5]),
                'instruction integer interface')
        op, a, b, out, owner, literal = instruction
        require(op < 11, 'unknown instruction opcode')
        if op in (0, 1, 2, 10):
            require(a == b == NONE and type(literal) is int, 'preparation literal interface')
            value = literal
        else:
            require(literal is None and a in memory, 'computed instruction literal or input')
            if op in (5, 8):
                require(b in memory, 'binary instruction input')
            else:
                require(b == NONE, 'unary instruction input')
            x = memory[a][0]
            if op in (3, 9):
                value = x
            elif op == 4:
                require(x == 0, 'reset instruction source')
                value = 0
            elif op == 5:
                total = x + memory[b][0]
                require(total % 2 == 0, 'non-half-unit mean')
                value = total // 2
            elif op == 6:
                value = 2*x
            elif op == 7:
                require(x == 0, 'start instruction source')
                value = 2
            elif op == 8:
                value = x + memory[b][0]
        require(-(1 << 63) <= value < (1 << 63), 'instruction scalar overflow')
        wa = NONE if a == NONE else memory[a][1]
        wb = NONE if b == NONE else memory[b][1]
        rows.append((op, a, b, out, owner, wa, wb, value & NONE))
        for destination in ((a, b) if op == 5 else (out,)):
            require(destination in memory or destination == len(memory),
                    'instruction allocation order')
            memory[destination] = (value, eid)
        require(len(memory) <= 30000, 'instruction register bound')
    return rows


def check_roundtrip(recovered, original_rows):
    raw = b''.join(struct.pack('<8Q', *row[:7], row[7] & NONE) for row in original_rows)
    current, iterations = recovered, []
    for iteration in (1, 2):
        rows = execute_program(current['program'])
        encoded = b''.join(struct.pack('<8Q', *row) for row in rows)
        require(raw == encoded, 'independent finite-program replay bytes')
        current = check_rows(rows)
        same(recovered, current, 'independent finite-program recovery changes')
        iterations.append({'iteration': iteration, 'decoded_sha256': hashlib.sha256(encoded).hexdigest(),
                           'recovered_specification_equal': True})
    return {'scope': 'observed_instruction_program_only', 'iterations': iterations}


def check_alternative_input(data, header, hosts, edges):
    require(isinstance(data, bytes), 'raw input bytes required')
    prefix = struct.pack('<8Q', *header) + struct.pack('<27I', *hosts)
    edge_bytes = b''.join(struct.pack('<4I', *row) for row in edges)
    require(data.startswith(prefix + edge_bytes), 'parsed input does not match raw bytes')
    replacement = b''.join(struct.pack('<4I', *row) for row in alternative_edges(edges))
    changed = prefix + replacement + data[len(prefix) + len(edge_bytes):]
    return hashlib.sha256(changed).hexdigest()


def verify_packet(packet, *, rows, header, hosts, edges, input_bytes, native_runs, source_pins):
    require(type(packet) is dict and set(packet) == {
        'schema', 'verdict', 'first_failed_invariant', 'fixed_point_attained',
        'assumptions', 'nonclaims', 'recovered', 'completions', 'native_replays',
        'finite_program_roundtrip', 'alternative_input_sha256', 'source_hashes'}, 'receipt field inventory')
    require(packet.get('schema') == 'oph.federation_recovery.v1', 'receipt schema')
    require(packet.get('verdict') == 'complete_gluing_not_identified'
            and packet.get('first_failed_invariant') == 'complete_gluing_triangle_count'
            and packet.get('fixed_point_attained') is False, 'false full-recovery assertion')
    same(ASSUMPTIONS, packet.get('assumptions'), 'assumption contract')
    same(NONCLAIMS, packet.get('nonclaims'), 'nonclaim contract')
    same(source_pins, packet.get('source_hashes'), 'source custody')
    require(hashlib.sha256(input_bytes).hexdigest() == source_pins.get('code/source_read_routing/controls/q3.input'),
            'input bytes custody')
    rows = list(rows)
    recovered = check_rows(rows)
    same(recovered, packet.get('recovered'), 'independent event recovery differs')
    inv = recovered['invariants']
    require(recovered['census']['events'] == 22818
            and inv['owners'] == 1280 and inv['population_sites'] == 27
            and inv['logical_nodes'] == 81 and inv['logical_strict_relations'] == 1331
            and inv['observed_glued_seams'] == 86, 'baseline salient census')
    completions = check_completions(header, hosts, edges)
    same(completions, packet.get('completions'), 'independent complete-graph comparison differs')
    same(check_alternative_input(input_bytes, header, hosts, edges),
         packet['alternative_input_sha256'], 'alternative input custody')
    roundtrip = check_roundtrip(recovered, rows)
    same(roundtrip, packet['finite_program_roundtrip'],
         'independent finite-program roundtrip differs')
    original_ports = {tuple(sorted((12*a+p, 12*b+q))) for a, p, b, q in edges}
    alternate_ports = {tuple(sorted((12*a+p, 12*b+q))) for a, p, b, q in alternative_edges(edges)}
    require(all(tuple(edge) in original_ports & alternate_ports
                for edge in recovered['observed_gluing']), 'observed seam absent from completion')
    same(native_runs, packet.get('native_replays'), 'fresh native replay differs')
    expected = {
        'baseline': 'b4fc16b1c168d1805aadd1219426777fcabebedf71b6fd26a68cb2aa357d6a57',
        'source': '47f281345b0a03a6369ad3b74852d632f8c4fc58f81c1e596e420b8e39101bc5',
        'branch': 'fd719adc1208ea8c4c2288fece5053df89a1607d6de105c0dd0ba66792bff921',
        'scratch': '047d7c261ea7ad0a96f7f9b6af097f21ca1fec5ffc87ae7b8facc9ede85d541f',
    }
    require(set(native_runs) == set(expected), 'complete native control inventory')
    require(all(row['decoded_sha256'] == expected['baseline'] for row in roundtrip['iterations']),
            'recovered word differs from authenticated native baseline')
    for variant, pin in expected.items():
        same({'events': 22818, 'decoded_sha256': pin, 'original_equals_alternative': True,
              'original_matches_retained': True}, native_runs[variant], 'native event identity')
    return {'verified': True, 'verdict': 'complete_gluing_not_identified',
            'fixed_point_attained': False, 'events': 22818,
            'triangle_counts': [14000, 13990], 'physical_identification': False}


def load_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def reject(value):
        raise ValueError('nonintegral JSON number: ' + value)
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs,
                      parse_float=reject, parse_constant=reject)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path,
                        default=Path(__file__).with_name('receipt.json'))
    args = parser.parse_args()
    import custody
    header, hosts, edges = custody.read_input()
    result = verify_packet(load_json(args.receipt), rows=custody.retained_rows(),
                           header=header, hosts=hosts, edges=edges,
                           input_bytes=(custody.CONTROL / 'q3.input').read_bytes(),
                           native_runs=custody.native_replays(), source_pins=custody.source_hashes())
    print(json.dumps(result, sort_keys=True))

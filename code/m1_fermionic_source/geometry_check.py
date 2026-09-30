"""Independent coordinate graph, sparse Pauli and high-precision ledger."""

from collections import deque
import hashlib
import itertools
import json
import math

import mpmath as mp
import numpy as np

from .check import close, complex_array, exact, keys, need


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode('ascii')).hexdigest()


def geometry(q):
    # Direct integer coordinates, not the producer's tuple-to-index builder.
    cells = list(itertools.product(range(q), repeat=3))
    def address(x, bank, mode):
        return 32*((x[0]*q+x[1])*q+x[2])+16*bank+mode
    edges, flights = set(), []
    directions = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    for x in cells:
        base = address(x, 0, 0)
        for mode in range(16):
            for other in range(mode+1, 16):
                if (mode//4 == other//4 and other-mode == 1) or mode ^ other in (4, 8):
                    edges.add((base+mode, base+other))
            edges.add((base+mode, base+16+mode))
            direction = np.array(directions[mode % 4])*(1 if mode % 8 < 4 else -1)
            y = tuple(int(a+b) for a, b in zip(x, direction))
            if all(0 <= t < q for t in y):
                target = address(y, 1, mode)
            else:
                target = base+16+(mode ^ 4)
            edges.add((base+mode, target) if base+mode < target else (target, base+mode))
            flights.append((base+mode, target))
        for axis in range(3):
            if x[axis]+1 < q:
                y = tuple(v+(k == axis) for k, v in enumerate(x))
                edges.add((base, address(y, 0, 0)))
    edges = sorted(edges)
    # A permutation between complete active and blank banks; reflecting
    # endpoints are local, unlike a long uncharged periodic wrap in E.
    need(len({v for _, v in flights}) == 16*q**3, 'bijective flight destinations')
    for i, j in edges:
        a, b = cells[i//32], cells[j//32]
        need(max(abs(s-t) for s, t in zip(a, b)) <= 1, 'finite geometric range')
    return edges, flights


TABLE = {('X', 'Y'): ('Z', 1), ('Y', 'Z'): ('X', 1), ('Z', 'X'): ('Y', 1),
         ('Y', 'X'): ('Z', 3), ('Z', 'Y'): ('X', 3), ('X', 'Z'): ('Y', 3)}


def loop_word(path, edges, incident):
    # A sparse tensor of I/X/Y/Z and its scalar, not symplectic multiplication.
    factors, phase = {}, (len(path)-1) % 4
    lookup = {tuple(e): k for k, e in enumerate(edges)}
    for left, right in zip(path, path[1:]):
        e = lookup[tuple(sorted((left, right)))]
        phase += 0 if left < right else 2
        local = [(e, 'X')]+[(j, 'Z') for v in (left, right) for j in incident[v] if j < e]
        for bit, letter in local:
            if bit not in factors:
                factors[bit] = letter
            elif factors[bit] == letter:
                del factors[bit]
            else:
                factors[bit], extra = TABLE[factors[bit], letter]
                phase += extra
    x = sum(1 << k for k, p in factors.items() if p in ('X', 'Y'))
    z = sum(1 << k for k, p in factors.items() if p in ('Z', 'Y'))
    phase += sum(p == 'Y' for p in factors.values())
    return [str(x), str(z), phase % 4], len(factors)


def graph_summary(q):
    edges, flights = geometry(q)
    n = 32*q**3
    neighbors, incident = [[] for _ in range(n)], [[] for _ in range(n)]
    for e, (u, v) in enumerate(edges):
        neighbors[u].append(v)
        neighbors[v].append(u)
        incident[u].append(e)
        incident[v].append(e)
    need(max(map(len, incident)) <= 12, 'degree twelve bound')
    parent, depth, queue = {0: None}, {0: 0}, deque([0])
    while queue:
        v = queue.popleft()
        for u in sorted(neighbors[v]):
            if u not in parent:
                parent[u], depth[u] = v, depth[v]+1
                queue.append(u)
    need(len(parent) == n, 'connected encoding graph')
    tree = {tuple(sorted((v, p))) for v, p in parent.items() if p is not None}
    words, weights = [], []
    for a, b in edges:
        if (a, b) in tree:
            continue
        x, y, left, right = a, b, [a], [b]
        while x != y:
            if depth[x] >= depth[y]:
                x = parent[x]
                left.append(x)
            else:
                y = parent[y]
                right.append(y)
        cycle = left+list(reversed(right[:-1]))+[a]
        raw, weight = loop_word(cycle, edges, incident)
        words.append(raw)
        weights.append(weight)
    return dict(side=q, vertices=n, edges=len(edges), graph_sha256=digest(edges),
                flight_sha256=digest(flights), maximum_degree=max(map(len, incident)),
                cycles=len(edges)-n+1, maximum_tree_depth=max(depth.values()),
                loop_weight_sum=sum(weights), maximum_loop_weight=max(weights, default=0),
                loop_words_sha256=digest(words))


def verify_graphs(rows):
    need(type(rows) is list and len(rows) == 4, 'all geometry cutoffs')
    for row, q in zip(rows, (1, 2, 3, 4)):
        exact(row, graph_summary(q))


def expected_budget(q, spacing, speed=3.):
    with mp.workdps(70):
        a, c = mp.mpf(str(spacing)), mp.mpf(str(speed))
        d, w, colors = 12, 23, 64
        cells = q**3
        vertices, edges = 32*cells, 192*cells
        height = 3*(q-1)+16
        length = 16
        cycles = edges-vertices+1
        weight = 368
        omega, eta = mp.pi*c/a, a/c
        ell = mp.sqrt(3)*a
        pulse_count, event_count = 231, 648
        gadget = 184*ell/c+pulse_count*mp.pi/omega+event_count*eta
        rotations = 16+32*6+32*8
        tick = 64*rotations*gadget
        prep = 11*eta+3072*((2944*mp.sqrt(3)+13989)*eta)+(2*q+2)*(mp.sqrt(3)+1024)*eta
        interfaces = 3*edges+2*(cells+1)
        out = dict(side=q, spacing=a, speed=c, cells=cells, vertices=vertices,
                   edge_qubits_upper=edges, helpers=cells+1, cycles_upper=cycles,
                   tree_depth_upper=height, cycle_length_upper=length, cycle_weight_upper=weight,
                   preparation_slots=3072, decoder_rounds=2*q+2,
                   native_drive_bound=omega, event_time=eta, colors=colors, pauli_weight_upper=w,
                   pulses_per_rotation_upper=pulse_count, code_events_per_rotation_upper=event_count,
                   rotations_per_cell_tick=rotations,
                   native_pulses_per_tick_upper=rotations*cells*pulse_count,
                   code_events_per_tick_upper=rotations*cells*event_count,
                   flights_per_tick_upper=rotations*cells*4*w, rotation_time_upper=gadget,
                   wall_tick=tick, wall_factor=tick/(mp.sqrt(3)*a/c),
                   preparation_time_upper=prep, noisy_interfaces_upper=interfaces,
                   phase_rate_for_point_zero_zero_one=mp.mpf('.001')/(interfaces*(prep+tick)))
        return {k: float(v) if isinstance(v, mp.mpf) else v for k, v in out.items()}


def verify_budgets(rows):
    need(type(rows) is list and len(rows) == 3, 'budget catalogue')
    for row, (q, a) in zip(rows, ((2, .1), (4, .03), (8, .01))):
        exact(row, expected_budget(q, a))


def closed_walk(q, spacing):
    # Direct precoin, reflecting permutation, postcoin, mass coupling.
    # This contains no QR decomposition and no blank-bank circuit.
    c = np.array([[0, 1, 1, 1], [1, 0, 1j, -1j], [1, -1j, 0, 1j], [1, 1j, -1j, 0]])/np.sqrt(3)
    n = 16*q**3
    before, after, coupling = np.eye(n, dtype=complex), np.eye(n, dtype=complex), np.eye(n, dtype=complex)
    for cell in range(q**3):
        base = 16*cell
        for mass_index, mass in enumerate((.7, 1.1)):
            offset = base+8*mass_index
            before[offset:offset+4, offset:offset+4] = c
            after[offset+4:offset+8, offset+4:offset+8] = c
            mu = mass*spacing/math.sqrt(3)
            for mode in range(4):
                pair = [offset+mode, offset+mode+4]
                coupling[np.ix_(pair, pair)] = [[math.cos(mu), -1j*math.sin(mu)], [-1j*math.sin(mu), math.cos(mu)]]
    _, flights = geometry(q)
    perm = np.zeros((n, n))
    for source, target in flights:
        i = (source//32)*16+source % 16
        j = (target//32)*16+target % 16
        perm[j, i] = 1
    return coupling@after@perm@before


def verify_walks(rows):
    need(type(rows) is list and len(rows) == 3, 'finite clock catalogue')
    for row, (q, a) in zip(rows, ((1, .07), (2, .13), (3, .03))):
        keys(row, 'side spacing layers pauli_rotations program output banks')
        exact([row['side'], row['spacing'], row['layers'], row['pauli_rotations']],
              [q, a, dict(phase=16, mix=32, ferry=16, restore=16), 464])
        n = 16*q**3
        expected = closed_walk(q, a)
        edges, flights = geometry(q)
        active = [32*x+j for x in range(q**3) for j in range(16)]
        blank = [v+16 for v in active]
        actual = np.zeros((2*n, n), complex)
        actual[active] = np.eye(n)
        program = row['program']
        need(type(program) is list and len(program) == 80, 'complete clock program')
        counts = dict(phase=0, mix=0, ferry=0, restore=0)
        for instruction in program:
            keys(instruction, 'kind modes parameter')
            kind, modes, parameter = (instruction[k] for k in ('kind', 'modes', 'parameter'))
            need(type(kind) is str and kind in counts, 'clock instruction kind')
            counts[kind] += 1
            need(type(modes) is list and len(modes) == (2 if kind == 'mix' else 1)
                 and all(type(v) is int and 0 <= v < 16 for v in modes), 'clock mode domain')
            if kind == 'mix':
                need(modes[0] != modes[1], 'distinct mixed modes')
                g = complex_array(parameter, (2, 2))
                close(g.conj().T@g, np.eye(2), 'clock SU2')
                need(abs(np.linalg.det(g)-1) < 1e-12, 'clock mixer determinant')
            elif kind == 'phase':
                need(type(parameter) is float and math.isfinite(parameter), 'finite clock phase')
            else:
                need(parameter is None, 'parameter-free flight')
            pairs = [(u, v) for u, v in flights if u % 16 == modes[0]] if kind == 'ferry' else [
                tuple(32*cell+v for v in (modes if kind == 'mix' else [modes[0], modes[0]+16]))
                for cell in range(q**3)]
            for u, v in pairs:
                if kind == 'phase':
                    actual[u] *= np.exp(1j*parameter)
                else:
                    need(tuple(sorted((u, v))) in edges, 'available spatial encoding edge')
                    actual[[u, v]] = (g@actual[[u, v]] if kind == 'mix' else actual[[v, u]])
        exact(counts, row['layers'])
        close(actual[active], expected, 'all entries of actual clock program')
        close(actual[blank], np.zeros((n, n)), 'all blank-bank columns restored')
        probe = np.array([complex(math.cos(j*j/17), math.sin(j*j/17)) for j in range(n)])/math.sqrt(n)
        close(complex_array(row['output'], (n,)), expected@probe, 'reflecting massive walk')
        close(complex_array(row['banks'], (n,)), np.zeros(n), 'all blank banks restored')



def verify_clock_witness(row):
    from m1_operational_clocks.check import clock_interval_certificate
    clock_interval_certificate()
    with mp.workdps(70):
        q, a, c = 2**41, mp.mpf('1e-9'), mp.mpf(3)
        cells, height = q**3, 3*(q-1)+16
        ell, tau = mp.sqrt(3)*a, a/mp.sqrt(3)
        eta = a/c
        gadget = (184*mp.sqrt(3)+879)*a/c
        tick = 29696*gadget
        prep = 11*eta+3072*(2944*mp.sqrt(3)+13989)*eta+(2*q+2)*(mp.sqrt(3)+1024)*eta
        packet = 147456*height*gadget
        read = 4608*gadget
        report = ell*(q-1)/c+2624*eta
        ticks = int(mp.ceil(25*mp.pi/tau))
        total = prep+packet+ticks*tick+read+report
        interfaces = 578*cells+2
        energy = 32*cells*mp.pi/tau
        trace_error = min(mp.mpf('1e-5'), mp.mpf('.01')/energy)/2
        rate = trace_error/(interfaces*total)
        setup_ops = 1920*cells+3072*15461*cells+(2*q+2)*1024*cells+2
        records = setup_ops+(ticks*464+2304*height+72)*65536*cells+64*cells
        bits = (records+32*cells+1).bit_length()
        expected = dict(workspace_side=q, spacing=a, speed=c, native_tick=tau, ticks=ticks,
            wall_factor=tick/tau, wall_tick=tick, gauge_preparation_time_upper=prep,
            packet_preparation_time_upper=packet, read_time_upper=read, report_time_upper=report,
            total_exposure_upper=total, noisy_interfaces_upper=interfaces, accounting_range_upper=energy,
            phase_rate_cap=rate, trace_error_upper=trace_error, accounting_error_upper=energy*trace_error,
            swing_lower=mp.mpf('.6753')-2*trace_error, wall_beat_frequency=mp.mpf('.08')*tau/tick,
            primitive_events_upper=records, record_identifier_bits=bits,
            central_record_slots_upper=(1+8*bits)*records)
        exact(row, {k: float(v) if isinstance(v, mp.mpf) else v for k, v in expected.items()})
        need(rate > 0 and energy*trace_error < .01 and expected['swing_lower'] > .6752, 'informative finite source witness')

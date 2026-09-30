"""Finite reflecting clock graph and a conflict-free source schedule."""

import hashlib
import itertools
import json
import math

from .model import Encoding


SIGNS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
DEGREE_BOUND = 12
WEIGHT_BOUND = 2*DEGREE_BOUND-1


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode('ascii')).hexdigest()


def graph(side):
    if type(side) is not int or side < 1:
        raise ValueError('positive integer side required')
    cells = list(itertools.product(range(side), repeat=3))
    vertices = [(x, bank, mode) for x in cells for bank in (0, 1) for mode in range(16)]
    ids = {key: i for i, key in enumerate(vertices)}
    pairs, moves = set(), []
    def edge(left, right):
        pairs.add(tuple(sorted((ids[left], ids[right]))))
    for x in cells:
        for group in range(4):
            for j in range(3):
                edge((x, 0, 4*group+j), (x, 0, 4*group+j+1))
        for mass in (0, 1):
            for j in range(4):
                edge((x, 0, 8*mass+j), (x, 0, 8*mass+j+4))
        for j in range(8):
            edge((x, 0, j), (x, 0, j+8))
        for mode in range(16):
            edge((x, 0, mode), (x, 1, mode))
            chirality = (mode % 8)//4
            delta = SIGNS[mode % 4]
            target = tuple(x[k]+(1 if chirality == 0 else -1)*delta[k] for k in range(3))
            if all(0 <= t < side for t in target):
                destination = (target, 1, mode)
            else:
                destination = (x, 1, mode ^ 4)
            moves.append((ids[(x, 0, mode)], ids[destination]))
            edge((x, 0, mode), destination)
        # Encoding-only links connect the full workspace; they are never
        # inserted as physical hopping terms into the clock Hamiltonian.
        for axis in range(3):
            y = list(x)
            y[axis] += 1
            if y[axis] < side:
                edge((x, 0, 0), (tuple(y), 0, 0))
    return vertices, sorted(pairs), moves


def graphs():
    rows = []
    for side in (1, 2, 3, 4):
        vertices, edges, moves = graph(side)
        enc = Encoding(len(vertices), edges)
        depth = [0]*len(vertices)
        for v in enc.parent:
            if v:
                depth[v] = depth[enc.parent[v]]+1
        weights = [s.weight() for s in enc.loops]
        rows.append(dict(side=side, vertices=len(vertices), edges=len(edges),
                         graph_sha256=digest(edges), flight_sha256=digest(moves),
                         maximum_degree=max(map(len, enc.incident)), cycles=len(enc.loops),
                         maximum_tree_depth=max(depth), loop_weight_sum=sum(weights),
                         maximum_loop_weight=max(weights, default=0),
                         loop_words_sha256=digest([s.row() for s in enc.loops])))
    return rows


def budget(side, spacing, speed=3.):
    """Conservative global bounds, with every setup and clock cost positive."""
    if type(side) is not int or side < 1 or any(type(t) not in (int, float) or
            isinstance(t, bool) or not math.isfinite(t) or t <= 0 for t in (spacing, speed)):
        raise ValueError('finite positive geometry required')
    w = WEIGHT_BOUND
    cells = side**3
    vertices = 32*cells
    # In-cell diameter <= 16; cell-zero modes are linked along the axes.
    tree_depth = 3*(side-1)+16
    loop_length = 16
    edge_upper = DEGREE_BOUND*vertices//2
    cycles_upper = edge_upper-vertices+1
    loop_weight = loop_length*(2*DEGREE_BOUND-1)
    omega = math.pi*speed/spacing
    eta = spacing/speed
    flight_length = math.sqrt(3)*spacing
    # Four tours / event-pair charges are padded even for trivial phases.
    pulse_bound = 10*w+1
    event_bound = 28*w+4
    gadget = 8*w*flight_length/speed+pulse_bound*math.pi/omega+event_bound*eta
    colors = 64
    rotations_per_cell = 4*(4+6*6)+8*6+32*8
    tick = colors*rotations_per_cell*gadget
    factor = tick/(math.sqrt(3)*spacing/speed)
    # Local loop instruments run in 64 colors and at most 48 anchor slots.
    # The exact prefix decoder uses 2(q-1) rounds plus four local rounds.
    # Classical messages and conditional corrections are charged as well.
    loop_gadget = (8*loop_weight*flight_length/speed
                   +(10*loop_weight+1)*math.pi/omega+(28*loop_weight+4)*eta)
    decoder_rounds = 2*(side-1)+4
    prepare = (10*eta+math.pi/omega+64*48*loop_gadget
               +decoder_rounds*(flight_length/speed+1024*eta))
    # One helper per cell plus one setup helper; all physical qubits count.
    # Data and visitor buffers at edge owners; helper buffers; every active
    # M6 processor. Each is one dephasing channel, regardless of dimension.
    noise_interfaces = 3*edge_upper+2*(cells+1)
    return dict(side=side, spacing=float(spacing), speed=float(speed), cells=cells,
                vertices=vertices, edge_qubits_upper=edge_upper, helpers=cells+1,
                cycles_upper=cycles_upper, tree_depth_upper=tree_depth,
                cycle_length_upper=loop_length, cycle_weight_upper=loop_weight,
                preparation_slots=3072, decoder_rounds=decoder_rounds,
                native_drive_bound=omega, event_time=eta, colors=colors,
                pauli_weight_upper=w, pulses_per_rotation_upper=pulse_bound,
                code_events_per_rotation_upper=event_bound,
                rotations_per_cell_tick=rotations_per_cell,
                native_pulses_per_tick_upper=rotations_per_cell*cells*pulse_bound,
                code_events_per_tick_upper=rotations_per_cell*cells*event_bound,
                flights_per_tick_upper=rotations_per_cell*cells*4*w,
                rotation_time_upper=gadget, wall_tick=tick, wall_factor=factor,
                preparation_time_upper=prepare, noisy_interfaces_upper=noise_interfaces,
                phase_rate_for_point_zero_zero_one=0.001/(noise_interfaces*(prepare+tick)))


def budgets():
    return [budget(q, a) for q, a in ((2, .1), (4, .03), (8, .01))]


def clock_witness():
    # Conservative binary reports; the independent high-precision checker
    # verifies the direction and every composed inequality, not just closeness.
    def up(x):
        for _ in range(32):
            x = math.nextafter(x, math.inf)
        return x
    def down(x):
        for _ in range(32):
            x = math.nextafter(x, 0.)
        return x
    q, a, c = 2**41, 1e-9, 3.
    r = budget(q, a, c)
    tau = math.sqrt(3)*a/c
    ticks = math.ceil(25*math.pi/tau)
    prep, tick = up(r['preparation_time_upper']), up(r['wall_tick'])
    packet = up(6*12*32*64*r['tree_depth_upper']*r['rotation_time_upper'])
    read = up(72*64*r['rotation_time_upper'])
    report = up(math.sqrt(3)*a*(q-1)/c+64*41*r['event_time'])
    total = up(prep+packet+ticks*tick+read+report)
    # Extend the fixed positive Fock accounting observable with its upper
    # spectral endpoint on invalid-code outcomes. No discarded bad syndrome.
    energy_range = up(32*q**3*math.pi/tau)
    error_budget = min(1e-5, .01/energy_range)
    rate = down(error_budget/(2*r['noisy_interfaces_upper']*total))
    failure = up(rate*r['noisy_interfaces_upper']*total)
    # All elementary quantum events, idle slots and classical OR reports
    # receive retained records. This is a conservative finite allocation.
    setup_ops = (10*r['edge_qubits_upper']+64*48*(38*368+5+4*368)*q**3
                 +(2*q+2)*1024*q**3+2)
    quantum_ops = (ticks*464+6*12*32*r['tree_depth_upper']+72)*64*q**3*1024
    records = setup_ops+quantum_ops+64*q**3
    bits = (records+32*q**3+1).bit_length()
    return dict(workspace_side=q, spacing=a, speed=c, native_tick=tau, ticks=ticks,
                wall_factor=r['wall_factor'], wall_tick=tick,
                gauge_preparation_time_upper=prep,
                packet_preparation_time_upper=packet, read_time_upper=read,
                report_time_upper=report, total_exposure_upper=total,
                noisy_interfaces_upper=r['noisy_interfaces_upper'],
                accounting_range_upper=energy_range, phase_rate_cap=rate,
                trace_error_upper=failure, accounting_error_upper=up(energy_range*failure),
                swing_lower=down(.6753-2*failure), wall_beat_frequency=.08/r['wall_factor'],
                primitive_events_upper=records, record_identifier_bits=bits,
                central_record_slots_upper=(1+8*bits)*records)

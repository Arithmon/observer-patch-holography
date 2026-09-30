"""Executed proper-code gates, spatial flights, preparation and resource ledger."""

import itertools
import math

import numpy as np
from scipy.linalg import expm

from .model import decompose


SIGNS = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)])
COIN = np.array([[0, 1, 1, 1], [1, 0, 1j, -1j], [1, -1j, 0, 1j],
                 [1, 1j, -1j, 0]], complex)/np.sqrt(3)


def encode(z):
    return np.stack((np.asarray(z).real, np.asarray(z).imag), axis=-1).tolist()


def native_pair(g):
    """Occupation order 00,01,10,11; one-particle order first, second."""
    native = np.eye(6, dtype=complex)
    native[1:3, 1:3] = g[::-1, ::-1]
    return native


def coin_instructions(offset):
    diagonal, moves = decompose(COIN)
    tape = [('phase', [offset+j], np.array([[z]])) for j, z in enumerate(diagonal)]
    tape += [('pair', [offset+i, offset+j], native_pair(g)) for i, j, g in moves]
    return tape


def instructions(a, mass):
    tau = math.sqrt(3)*a/3
    mix = expm(-1j*mass*tau*np.array([[0, 1], [1, 0]]))
    return coin_instructions(0)+[('flight', [], None)]+coin_instructions(4)+[
        ('pair', [j, j+4], native_pair(mix)) for j in range(4)]


def execute(q, a, mass):
    """Every basis input of the entire periodic one-particle sector."""
    size = 8*q**3
    state = np.eye(size, dtype=complex).reshape(q, q, q, 8, size)
    for kind, ports, native in instructions(a, mass):
        if kind == 'flight':
            old = state.copy()
            for j in range(8):
                shift = SIGNS[j % 4]*(1 if j < 4 else -1)
                state[..., j, :] = np.roll(old[..., j, :], tuple(shift), axis=(0, 1, 2))
        elif kind == 'phase':
            state[..., ports[0], :] *= native[0, 0]
        else:
            # Encode into |10>,|01>, pulse on the full six-level processor,
            # decode. The proper-code complement is tested separately.
            lifted = np.zeros(state.shape[:3]+(6, size), complex)
            lifted[..., 2, :] = state[..., ports[0], :]
            lifted[..., 1, :] = state[..., ports[1], :]
            lifted = np.einsum('ij,...jk->...ik', native, lifted)
            state[..., ports[0], :], state[..., ports[1], :] = lifted[..., 2, :], lifted[..., 1, :]
    return state.reshape(size, size)


def spatial():
    rows = []
    for q, a, mass in itertools.product((2, 3), (.07, .013), (.7, 1.1)):
        u = execute(q, a, mass)
        tape = instructions(a, mass)
        pulses = sum(kind != 'flight' for kind, _, _ in tape)
        flights = 8*sum(kind == 'flight' for kind, _, _ in tape)
        size = len(u)
        probe = np.exp(1j*np.arange(size)**2/17)/np.sqrt(size)
        out = probe.copy()
        for _ in range(3):
            out = u@out
        rows.append(dict(q=q, a=a, mass=mass, modes=size, instructions=len(tape),
                         pulses_per_cell=pulses, code_events_per_cell=2*(pulses+flights), flights_per_cell=flights,
                         column_norms=np.sum(abs(u)**2, axis=0).tolist(),
                         three_step_probe=encode(out)))
    return rows


def tree_prepare(target):
    """Real amplitude splitting on a complete octree, then local phases.

    Returns the actual leaf amplitudes and counts all allocated branches,
    including zero-weight ones. The small executions are not the huge witness.
    """
    side = target.shape[0]
    if target.shape != (side, side, side) or side < 2 or side & (side-1):
        raise ValueError('complete dyadic cube required')
    if not np.isclose(np.linalg.norm(target), 1):
        raise ValueError('normalized target required')
    result = np.zeros_like(target)
    pulses = flights = 0

    def visit(origin, width, amplitude):
        nonlocal pulses, flights
        if width == 1:
            z = target[origin]
            result[origin] = amplitude*(z/abs(z) if z else 1)
            pulses += 1  # Include the zero-angle phase instruction.
            return
        half = width//2
        children = [tuple(origin[j]+half*b[j] for j in range(3))
                    for b in itertools.product((0, 1), repeat=3)]
        norms = np.array([np.linalg.norm(target[tuple(slice(x, x+half) for x in child)])
                          for child in children])
        desired = norms/np.linalg.norm(norms) if np.linalg.norm(norms) else np.eye(8)[0]
        amplitudes = np.zeros(8)
        amplitudes[0] = amplitude
        for j in range(7, 0, -1):
            rest = np.linalg.norm(desired[:j+1])
            sine = desired[j]/rest if rest else 0
            cosine = np.linalg.norm(desired[:j])/rest if rest else 1
            amplitudes[[0, j]] = np.array([[cosine, -sine], [sine, cosine]])@amplitudes[[0, j]]
            pulses += 1
        flights += 8
        for child, value in zip(children, amplitudes):
            visit(child, half, value)

    visit((0, 0, 0), side, 1.)
    return result, pulses, flights


def preparations():
    rows = []
    for side in (2, 4, 8):
        x = np.indices((side, side, side))
        for holes in (False, True):
            target = np.exp(-sum((z-side/2+.5)**2 for z in x)/side+1j*(x[0]+2*x[1]-x[2])/3)
            if holes:
                target[:side//2, :side//2, :side//2] = 0
            target /= np.linalg.norm(target)
            actual, pulses, flights = tree_prepare(target)
            rows.append(dict(side=side, holes=holes, pulses=pulses, flights=flights,
                             amplitudes=encode(actual)))
    return rows


def noise():
    rows = []
    for modes in (4, 8, 16):
        psi = np.exp(1j*np.arange(modes)/3)/np.sqrt(modes)
        rho = np.outer(psi, psi.conj())
        groups = np.arange(modes)//2
        for exposure in (.03, .4):
            s = math.exp(-exposure)
            actual = rho.copy()
            # Apply each pair's local code dephasing, including its vacuum.
            for group in set(groups):
                labels = np.where(groups == group, 1+np.arange(modes) % 2, 0)
                actual *= np.where(labels[:, None] == labels[None, :], 1., s)
            distance = np.linalg.eigvalsh(actual-rho)
            rows.append(dict(modes=modes, exposure=exposure, output=encode(actual),
                             half_trace_distance=float(sum(abs(distance))/2),
                             population_free_upper=1-s*s,
                             vacuum_bitflip_failure=1-(1-.02)**modes))
    return rows


def resource_witness():
    a, c, overhead, rate = 1e-9, 3., .01, 1e-14
    tau = math.sqrt(3)*a/c
    wall_tick = (1+overhead)*tau
    ticks = math.ceil(2*math.pi*1.25/(.1*tau))
    side = 2**math.ceil(math.log2(2*(90+a)/a))
    depth = side.bit_length()-1
    # 48 pulses, 96 pack/unpack and 32 departure/arrival events per tick.
    tapes = instructions(a, 100.)+instructions(a, 100.1)
    pulse_count = sum(kind != 'flight' for kind, _, _ in tapes)
    flight_count = 8*sum(kind == 'flight' for kind, _, _ in tapes)
    code_count = 2*(pulse_count+flight_count)
    omega = 2*pulse_count*math.pi/(overhead*tau)
    service = overhead*tau/(2*code_count)
    # Seven split rotations per tree level, then 15 internal rotations and
    # 16 phases at each leaf, all branches performed in parallel.
    prep_layers = 7*depth+31
    layer_time = math.pi/omega+2*service
    prep_flight = math.sqrt(3)*a*(side-1)/2/c
    prep_time = prep_flight+prep_layers*layer_time+(2*depth+2)*service
    # Eight mass-beam rotations, phases, local read and record timestamp.
    read_time = 16*layer_time+service
    run_time = ticks*wall_tick
    exposure = prep_time+run_time+read_time
    probability_noise = 2*rate*exposure
    workspace_side = 2**math.ceil(math.log2(2*1024/a))
    buffers = 32*workspace_side**3+16*(side**3-1)//7
    prep_pulses, prep_flights = 32*side**3-1, 8*(side**3-1)//7
    identifier_bits = (1+buffers+prep_pulses+prep_flights+ticks*workspace_side**3).bit_length()
    record_bound = (1+8*identifier_bits)*(4*buffers+4*prep_pulses+4*prep_flights
                     +(256*ticks+100)*workspace_side**3)
    return dict(a=a, c=c, overhead_fraction=overhead, phase_rate=rate,
                flight_tick=tau, wall_tick=wall_tick, ticks=ticks,
                preparation_side_cells=side, preparation_leaves=side**3,
                preparation_depth=depth, preparation_pulses=side**3-1+31*side**3,
                preparation_flights=8*(side**3-1)//7,
                workspace_side_cells=workspace_side,
                mode_buffers=32*workspace_side**3+16*(side**3-1)//7,
                active_processors=workspace_side**3+(side**3-1)//7,
                record_identifier_bits=identifier_bits, central_record_slots_upper=record_bound,
                pulses_per_cell_tick=pulse_count, code_events_per_cell_tick=code_count,
                register_flights_per_cell_tick=flight_count, drive_norm_bound=omega,
                code_event_time=service, preparation_time=prep_time,
                run_time=run_time, read_time=read_time,
                matter_speed=1/(1+overhead), clock_velocity=.6/(1+overhead),
                beat_frequency=.1/(1.25*(1+overhead)),
                probability_noise_upper=probability_noise,
                accounting_energy_error_upper=math.pi/tau*probability_noise,
                observable_swing_lower=.6753-2*probability_noise,
                reporting_delay_upper=math.sqrt(3)*600/c)


def candidate():
    return dict(spatial=spatial(), preparations=preparations(), noise=noise(), resources=resource_witness())

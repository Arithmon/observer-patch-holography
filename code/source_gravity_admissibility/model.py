"""Producer: execute the inherited twelve-port controls on proper codes."""
from fractions import Fraction as F
from functools import lru_cache
import itertools
import math

import numpy as np
from scipy.linalg import expm
import sympy as sp

from m1_source_realization.model import controls, source_generators, instruments


STRENGTHS = ('0', '1/20', '1/10', '1/5')
POLICIES = ('flat', 'curved')


def encode(a):
    return np.stack((np.asarray(a).real, np.asarray(a).imag), axis=-1).tolist()


@lru_cache(maxsize=1)
def control_data():
    return controls()


@lru_cache(maxsize=1)
def generators():
    gs = source_generators()
    out = {}
    for row in control_data()['rows']:
        coeff = [sp.Rational(a)+sp.sqrt(5)*sp.Rational(b) for a, b in row['coefficients']]
        out[row['name']] = np.array(sum((x*g for x, g in zip(coeff, gs)), sp.zeros(6)), complex)
    return out


@lru_cache(maxsize=128)
def pulse(kind, angle):
    return expm(angle*generators()['phase_1' if kind == 'phase' else 'real_12'])


def basis_word(axis):
    # Chronological product diag(I,sigma_axis) (Hadamard tensor I).
    out = [('phase', [2], F(1)), ('real', [0, 2], F(1, 4)),
           ('phase', [3], F(1)), ('real', [1, 3], F(1, 4))]
    if axis == 0:
        out += [('phase', [3], F(1)), ('real', [2, 3], F(1, 2))]
    elif axis == 1:
        out += [('real', [2, 3], F(1, 2)), ('phase', [2], F(1, 2)),
                ('phase', [3], F(1, 2))]
    else:
        out += [('phase', [3], F(1))]
    return out


def word_packet():
    return [[dict(kind=k, modes=m, pi=str(x)) for k, m, x in basis_word(axis)] for axis in range(3)]


def profile(shape, strength):
    """Exact grounded source solve, not a prescribed detector profile."""
    cells = list(itertools.product(*(range(n) for n in shape)))
    interior = [x for x in cells if all(0 < x[j] < shape[j]-1 for j in range(len(shape)))]
    middle = tuple(n//2 for n in shape)
    mat = sp.zeros(len(interior))
    for i, x in enumerate(interior):
        mat[i, i] = 2*len(shape)
        for j, y in enumerate(interior):
            if sum(abs(a-b) for a, b in zip(x, y)) == 1:
                mat[i, j] = -1
    rhs = sp.Matrix([sp.Rational(str(strength)) if x == middle else 0 for x in interior])
    solution = mat.inv()*rhs
    values = dict(zip(interior, map(F, solution)))
    return [values.get(x, F(0)) for x in cells]


def local(state, site, kind, modes, angle):
    """Load the proper code, execute a full M6 pulse, and unpack."""
    native = pulse(kind, float(angle))
    lifted = np.zeros((6, state.shape[-1]), complex)
    if kind == 'phase':
        lifted[1] = state[site, modes[0]]
        lifted = native@lifted
        state[site, modes[0]] = lifted[1]
    else:
        lifted[2], lifted[1] = state[site, modes[0]], state[site, modes[1]]
        lifted = native@lifted
        state[site, modes[0]], state[site, modes[1]] = lifted[2], lifted[1]


def flight(state, shape, axis):
    """All registers, including empty ones; actual edges or local reflection."""
    cells = list(itertools.product(*(range(n) for n in shape)))
    addresses = {x: i for i, x in enumerate(cells)}
    moved = np.zeros_like(state)
    # In a line the selected Clifford axis uses its single spatial direction.
    direction = 0 if len(shape) == 1 else axis
    for i, x in enumerate(cells):
        for mode in range(4):
            delta = 1 if mode < 2 else -1
            y = list(x)
            y[direction] += delta
            target = mode
            if not 0 <= y[direction] < shape[direction]:
                y = list(x)
                target = (mode+2) % 4
            moved[addresses[tuple(y)], target] = state[i, mode]
    return moved


def program(axes):
    """A fixed finite-state controller word; no outcome-dependent scheduling."""
    out = []
    for axis in axes:
        word = basis_word(axis)
        inverse = [(k, modes, -a) for k, modes, a in reversed(word)]
        for sign in (-1, 1):
            out.extend(dict(kind=k, modes=modes, angle='pi', coefficient=str(a)) for k, modes, a in inverse)
            out.append(dict(kind='flight', axis=axis))
            out.extend(dict(kind=k, modes=modes, angle='pi', coefficient=str(a)) for k, modes, a in word)
            out.extend(dict(kind='phase', modes=[mode], angle='source_acos',
                            coefficient=str(sign*(1 if mode < 2 else -1))) for mode in range(4))
    out.extend(dict(kind='phase', modes=[mode], angle='mass_lapse',
                    coefficient=str(-1 if mode < 2 else 1)) for mode in range(4))
    return out


def timeline(axes):
    word = program(axes)
    pulse_count = sum(op['kind'] != 'flight' for op in word)
    flight_count = len(word)-pulse_count
    delta = F(2, 3*pulse_count+2*flight_count)
    time = F(0)
    events = []
    for index, op in enumerate(word):
        # Local mode owners load flight slots in parallel. Gates on each
        # site's workspace are serial; all sites use the same stage counter.
        duration = 1+2*delta if op['kind'] == 'flight' else 3*delta
        events.append(dict(stage=index, instruction=op, start=str(time),
                           load_done=str(time+delta),
                           action_done=str(time+duration-delta), end=str(time+duration)))
        time += duration
    return dict(axes=list(axes), service_slot=str(delta), stages=events,
                interrogation_end=str(time), receiver_record_ready=str(time+1),
                returned_record_ready=str(time+3))


def execute(shape, values, policy, axes, mass_pi=F(0)):
    """Every column of the finite reflecting one-particle evolution."""
    sites = math.prod(shape)
    state = np.eye(4*sites, dtype=complex).reshape(sites, 4, 4*sites)
    r = [F(1, 2)/(1+u)**(2 if policy == 'curved' else 0) for u in values]
    angles = [math.acos(float(x)) for x in r]
    for op in program(axes):
        if op['kind'] == 'flight':
            state = flight(state, shape, op['axis'])
            continue
        for site, u in enumerate(values):
            base = (math.pi if op['angle'] == 'pi' else angles[site]
                    if op['angle'] == 'source_acos' else math.pi*float(mass_pi/(1+u)))
            local(state, site, op['kind'], op['modes'], base*float(F(op['coefficient'])))
    return state.reshape(4*sites, 4*sites)


def initial(shape, axis=None):
    sites = math.prod(shape)
    if axis is None:
        psi = np.exp(1j*np.arange(4*sites)**2/17)/np.sqrt(4*sites)
    else:
        # Prepare |site=1> tensor the first positive alpha_axis eigenvector
        # by the same explicit native basis word used by the flight compiler.
        state = np.zeros((sites, 4, 1), complex)
        state[1, 0, 0] = 1
        for kind, modes, angle in basis_word(axis):
            local(state, 1, kind, modes, math.pi*float(angle))
        psi = state.reshape(4*sites)
    return psi


def ledger(shape, axes):
    sites = math.prod(shape)
    pulses = sites*(sum(4*len(basis_word(axis))+8 for axis in axes)+4)
    flights = 8*sites*len(axes)
    # All boundary slots consume one duration even though their displacement is zero.
    waits = sum(8*sites//shape[0 if len(shape) == 1 else axis] for axis in axes)
    return dict(sites=sites, mode_slots=4*sites, native_pulses=pulses,
                proper_code_events=2*(pulses+flights),
                flight_slots=flights, boundary_wait_slots=waits,
                moving_slots=flights-waits, flight_layers=2*len(axes),
                service_time='2', propagation_time=str(2*len(axes)+2),
                receiver_read_time='1', report_distance='2', report_time='2',
                total_probe_and_report_time=str(2*len(axes)+5))


def native_basis(axis):
    state = np.eye(4, dtype=complex).reshape(1, 4, 4)
    for kind, modes, angle in basis_word(axis):
        local(state, 0, kind, modes, math.pi*float(angle))
    return state[0]


def phase_matrix(angles):
    state = np.eye(4, dtype=complex).reshape(1, 4, 4)
    for mode, angle in enumerate(angles):
        local(state, 0, 'phase', [mode], angle)
    return state[0]


def clock_evidence():
    rows = []
    for strength, elapsed, policy in itertools.product(STRENGTHS, ('1/7', '2/7'), POLICIES):
        a = 1/(1+F(strength))
        matrix = phase_matrix([-math.pi*float(F(elapsed)*a)*b for b in (1, 1, -1, -1)])
        plus = np.array([1, 0, 1, 0])/math.sqrt(2)
        rows.append(dict(strength=strength, mass_time_pi=elapsed, policy=policy,
                         matrix=encode(matrix.reshape(-1)), fringe=float(abs(np.vdot(plus, matrix@plus))**2)))
    return rows


def symbols():
    rows = []
    beta = np.diag([1., 1., -1., -1.]).astype(complex)
    for strength, policy, spacing in itertools.product(('0', '1/10'), POLICIES, ('1/20', '1/40')):
        lapse = 1/(1+F(strength))
        a = F(spacing)
        tau = 8*a  # Six paid flights and a common positive 2a service budget.
        r = F(1, 2)*(lapse**2 if policy == 'curved' else 1)
        theta = math.acos(float(r))
        rotation = phase_matrix([theta*b for b in (1, 1, -1, -1)])
        half = phase_matrix([theta*b/2 for b in (1, 1, -1, -1)])
        matrix = np.eye(4, dtype=complex)
        linear = float(F(1, 5)*tau*lapse)*beta
        momenta = (F(1, 3), F(1, 5), F(-1, 7))
        for axis, k in enumerate(momenta):
            b = native_basis(axis)
            alpha = b@beta@b.conj().T
            flight = b@phase_matrix([-float(a*k)*x for x in (1, 1, -1, -1)])@b.conj().T
            matrix = rotation@flight@rotation.conj().T@flight@matrix
            linear += float(2*a*r*k)*(half@alpha@half.conj().T)
        matrix = phase_matrix([-float(F(1, 5)*tau*lapse)*b for b in (1, 1, -1, -1)])@matrix
        length = F(1, 5)*tau*lapse+2*a*sum(map(abs, momenta))
        rows.append(dict(strength=strength, policy=policy, spacing=spacing, period=str(tau),
                         matrix=encode(matrix.reshape(-1)),
                         speed_relative=str(2*r), spatial_scale=str(lapse/(2*r)),
                         remainder_bound=str(length*length/2),
                         linear_error=float(np.linalg.norm(matrix-np.eye(4)+1j*linear, ord=2))))
    return rows


def bands():
    rows = []
    for strength, policy, momentum_pi, mass_pi in itertools.product(
            ('0', '1/10', '1/5'), POLICIES, ('0', '1/7', '1/3'), ('0', '1/11')):
        lapse = 1/(1+F(strength))
        r = F(1, 2)*(lapse**2 if policy == 'curved' else 1)
        theta = math.acos(float(r))
        rotation = phase_matrix([theta*b for b in (1, 1, -1, -1)])
        basis = native_basis(0)
        p = math.pi*float(F(momentum_pi))
        flight = basis@phase_matrix([-p*b for b in (1, 1, -1, -1)])@basis.conj().T
        mu = math.pi*float(F(mass_pi)*lapse)
        matrix = phase_matrix([-mu*b for b in (1, 1, -1, -1)])@rotation@flight@rotation.conj().T@flight
        rows.append(dict(strength=strength, policy=policy, momentum_pi=momentum_pi, mass_pi=mass_pi,
                         matrix=encode(matrix.reshape(-1)),
                         band_cosine=float(np.trace(matrix).real/4),
                         eigenphases=sorted(np.angle(np.linalg.eigvals(matrix)).tolist())))
    return rows


def candidate():
    rows = []
    for strength, axis, policy in itertools.product(STRENGTHS, range(3), POLICIES):
        s = F(strength)
        values = profile((5,), s)
        unitary = execute((5,), values, policy, [axis])
        result = unitary@initial((5,), axis)
        rows.append(dict(strength=strength, axis=axis, policy=policy,
                         source=list(map(str, values)),
                         cosine=[str(F(1, 2)/(1+u)**(2 if policy == 'curved' else 0)) for u in values],
                         clock_ratio=str(1/(1+s)),
                         result=encode(result), click=float(sum(abs(result[12:16])**2)),
                         receiver=[3], report_to=[1],
                         exact_click=str(F(1, 4)/(1+s)**(4 if policy == 'curved' else 0)),
                         ledger=ledger((5,), [axis])))
    cubes = []
    for strength, policy, mass in itertools.product(('1/10', '1/5'), POLICIES, ('0', '1/7')):
        values = profile((3, 3, 3), F(strength))
        unitary = execute((3, 3, 3), values, policy, [0, 1, 2], F(mass))
        result = unitary@initial((3, 3, 3))
        cubes.append(dict(strength=strength, policy=policy, mass_pi=mass,
                          source=list(map(str, values)),
                          result=encode(result), receiver=[2, 0, 0], report_to=[0, 0, 0],
                          click=float(sum(abs(result[72:76])**2)),
                          ledger=ledger((3, 3, 3), [0, 1, 2])))
    # Exact finite budgets, chosen before comparing approximate executions.
    robustness = []
    for strength in STRENGTHS[1:]:
        gap = F(1, 4)*(1-1/(1+F(strength))**4)
        budget = gap/8
        robustness.append(dict(strength=strength, exact_gap=str(gap),
                               per_branch_total_error=str(budget),
                               remaining_separation=str(gap-2*budget)))
    return dict(controls=control_data(), instruments=instruments(),
                basis_words=word_packet(), line=rows, cube=cubes, robustness=robustness,
                clocks=clock_evidence(), symbols=symbols(),
                timelines=[timeline(axes) for axes in ([0], [1], [2], [0, 1, 2])], bands=bands())

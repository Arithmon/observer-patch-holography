"""Independent source equations, projector flights and exact detector proof.

No import of the producer, its compiler, source solver or gate constructor.
The inherited independent checker supplies exact port-response arithmetic.
"""
from fractions import Fraction as F
import itertools
import math
import re

import numpy as np

from m1_source_realization.check import verify_controls, expected_instruments


def need(ok, message):
    if not ok:
        raise ValueError(message)


def keys(row, names):
    need(type(row) is dict and row.keys() == set(names.split()), 'exact object keys')


def fraction(value):
    need(type(value) is str and len(value) < 100 and
         re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', value) is not None,
         'bounded rational string')
    x = F(value)
    need(str(x) == value and max(x.numerator.bit_length(), x.denominator.bit_length()) < 256,
         'canonical rational')
    return x


def same(value, expected, label='evidence'):
    need(type(value) is type(expected), label+': type')
    if type(expected) is dict:
        need(value.keys() == expected.keys(), label+': keys')
        for k in expected:
            same(value[k], expected[k], label+'.'+k)
    elif type(expected) is list:
        need(len(value) == len(expected), label+': length')
        for i, (x, y) in enumerate(zip(value, expected)):
            same(x, y, label+f'[{i}]')
    elif type(expected) is float:
        need(math.isfinite(value) and abs(value-expected) <= 4e-12, label+': numerical replay')
    else:
        need(value == expected, label+': value')


def decode(value, size):
    need(type(value) is list and len(value) == size, 'full vector census')
    for cell in value:
        need(type(cell) is list and len(cell) == 2 and
             all(type(x) is float and math.isfinite(x) for x in cell), 'finite complex entry')
    array = np.asarray(value)
    return array[:, 0]+1j*array[:, 1]


I = np.eye(4, dtype=complex)
BETA = np.diag([1., 1., -1., -1.]).astype(complex)
PAULI = (np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.]))
ALPHA = tuple(np.block([[np.zeros((2, 2)), p], [p, np.zeros((2, 2))]]) for p in PAULI)


def basis(axis):
    return np.block([[np.eye(2), np.eye(2)], [PAULI[axis], -PAULI[axis]]])/math.sqrt(2)


def verify_words(words):
    need(type(words) is list and len(words) == 3, 'three complete Clifford-axis words')
    out = []
    for axis, word in enumerate(words):
        need(type(word) is list and len(word) == (6, 7, 5)[axis], 'complete native basis word')
        matrix = I.copy()
        for step in word:
            keys(step, 'kind modes pi')
            need(step['kind'] in ('phase', 'real'), 'native instruction kind')
            modes = step['modes']
            need(type(modes) is list and len(modes) == (1 if step['kind'] == 'phase' else 2)
                 and len(set(modes)) == len(modes)
                 and all(type(m) is int and 0 <= m < 4 for m in modes), 'proper-code addresses')
            angle = fraction(step['pi'])
            need(abs(angle) <= 2, 'bounded basis angle')
            t = math.pi*float(angle)
            gate = I.copy()
            if step['kind'] == 'phase':
                gate[modes[0], modes[0]] = complex(math.cos(t), math.sin(t))
            else:
                c, s = math.cos(t), math.sin(t)
                gate[np.ix_(modes, modes)] = [[c, -s], [s, c]]
            matrix = gate@matrix
        need(np.linalg.norm(matrix-basis(axis)) < 2e-14, 'native word changes the complete flight basis')
        out.append(matrix)
    return out


def source_check(strings, shape, strength):
    sites = list(itertools.product(*(range(n) for n in shape)))
    need(type(strings) is list and len(strings) == len(sites), 'complete source field')
    values = {x: fraction(s) for x, s in zip(sites, strings)}
    middle = tuple(n//2 for n in shape)
    for x, value in values.items():
        if any(x[j] in (0, shape[j]-1) for j in range(len(shape))):
            need(value == 0, 'actual grounded boundary')
        else:
            residual = 2*len(shape)*value
            for j in range(len(shape)):
                for direction in (-1, 1):
                    y = list(x)
                    y[j] += direction
                    residual -= values[tuple(y)]
            need(residual == (strength if x == middle else 0), 'exact local source equation')
        need(value >= 0, 'positive source maximum principle')
    return list(values.values())


def projector_flight(shape, axis):
    """Closed coordinate/projector formula, with an explicit reflected edge."""
    cells = list(itertools.product(*(range(n) for n in shape)))
    dim = len(cells)*4
    out = np.zeros((dim, dim), complex)
    plus, minus = (I+ALPHA[axis])/2, (I-ALPHA[axis])/2
    coordinate = 0 if len(shape) == 1 else axis
    for col, x in enumerate(cells):
        for sign, projector in ((1, plus), (-1, minus)):
            y = list(x)
            y[coordinate] += sign
            block = projector
            if y[coordinate] < 0 or y[coordinate] >= shape[coordinate]:
                y = list(x)
                block = BETA@projector
            row = cells.index(tuple(y))
            out[4*row:4*row+4, 4*col:4*col+4] += block
    need(np.linalg.norm(out.conj().T@out-np.eye(dim)) < 2e-13, 'complete reflected flight unitarity')
    return out


def word_flight(shape, axis, word_basis):
    """Replay the retained code-basis word and primitive register permutation."""
    cells = list(itertools.product(*(range(n) for n in shape)))
    p = np.zeros((4*len(cells), 4*len(cells)), complex)
    coordinate = 0 if len(shape) == 1 else axis
    for i, x in enumerate(cells):
        for channel in range(4):
            target = list(x)
            target[coordinate] += (-1 if channel >= 2 else 1)
            mode = channel
            if target[coordinate] not in range(shape[coordinate]):
                target = list(x)
                mode = channel ^ 2
            p[4*cells.index(tuple(target))+mode, 4*i+channel] = 1
    need(np.array_equal(p.sum(axis=0), np.ones(len(p))) and
         np.array_equal(p.sum(axis=1), np.ones(len(p))), 'all register branches retained')
    b = np.kron(np.eye(len(cells)), word_basis)
    return b@p@b.conj().T


def closed(shape, values, policy, axes, mass_pi=F(0), words=None):
    rs = [F(1, 2)/(1+u)**(2 if policy == 'curved' else 0) for u in values]
    phases = []
    for r in rs:
        z = complex(float(r), math.sqrt(float(1-r*r)))
        phases.extend((z, z, z.conjugate(), z.conjugate()))
    r = np.diag(phases)
    matrix = np.eye(4*len(values), dtype=complex)
    for axis in axes:
        s = projector_flight(shape, axis)
        if words is not None:
            need(np.linalg.norm(word_flight(shape, axis, words[axis])-s) < 3e-13,
                 'complete native flight differs from projector law')
        matrix = r@s@r.conj().T@s@matrix
    mass = [np.exp(-1j*math.pi*float(mass_pi/(1+u))*b) for u in values for b in (1, 1, -1, -1)]
    return np.asarray(mass)[:, None]*matrix


def expected_ledger(shape, axes):
    n = math.prod(shape)
    local_layers = sum(4*(6, 7, 5)[j]+8 for j in axes)+4
    flight_slots = n*4*2*len(axes)
    waits = sum(8*n//shape[0 if len(shape) == 1 else j] for j in axes)
    return dict(sites=n, mode_slots=4*n, native_pulses=n*local_layers,
                proper_code_events=2*(n*local_layers+flight_slots),
                flight_slots=flight_slots, boundary_wait_slots=waits,
                moving_slots=flight_slots-waits, flight_layers=2*len(axes),
                service_time='2', propagation_time=str(2*len(axes)+2),
                receiver_read_time='1', report_distance='2', report_time='2',
                total_probe_and_report_time=str(2*len(axes)+5))


def line_case(row, expected, words):
    keys(row, 'strength axis policy source cosine clock_ratio result click exact_click ledger receiver report_to')
    strength, axis, policy = expected
    same([row['strength'], row['axis'], row['policy']], [strength, axis, policy], 'line case census')
    s = F(strength)
    values = source_check(row['source'], (5,), s)
    same(row['cosine'], [str(F(1, 2)/(1+u)**(2 if policy == 'curved' else 0)) for u in values])
    same(row['clock_ratio'], str(1/(1+s)), 'common clock normalization')
    same(row['receiver'], [3], 'fixed local detector address')
    same(row['report_to'], [1], 'fixed returned record address')
    matrix = closed((5,), values, policy, [axis], words=words)
    state = np.zeros(20, complex)
    state[4:8] = basis(axis)[:, 0]
    result = matrix@state
    actual = decode(row['result'], 20)
    need(np.linalg.norm(actual-result) < 2e-12, 'all detector amplitudes replay')
    need(abs(np.vdot(actual, actual)-1) < 2e-12, 'no discarded output branch')
    p = F(1, 4)/(1+s)**(4 if policy == 'curved' else 0)
    same(row['exact_click'], str(p), 'exact source-conditioned detector theorem')
    same(row['click'], float(p), 'finite detector probability')
    need(abs(sum(abs(actual[12:16])**2)-float(p)) < 2e-12, 'actual local detector read')
    same(row['ledger'], expected_ledger((5,), [axis]), 'charged finite protocol')


def cube_case(row, expected, words):
    keys(row, 'strength policy mass_pi source result ledger receiver report_to click')
    strength, policy, mass = expected
    same([row['strength'], row['policy'], row['mass_pi']], [strength, policy, mass], 'cube case census')
    values = source_check(row['source'], (3, 3, 3), F(strength))
    matrix = closed((3, 3, 3), values, policy, [0, 1, 2], F(mass), words)
    need(np.linalg.norm(matrix.conj().T@matrix-np.eye(108)) < 3e-12, 'full 3D evolution, every mode')
    psi = np.array([complex(math.cos(j*j/17), math.sin(j*j/17)) for j in range(108)])/math.sqrt(108)
    actual = decode(row['result'], 108)
    need(np.linalg.norm(actual-matrix@psi) < 3e-12, 'complete 3D finite evolution')
    same(row['receiver'], [2, 0, 0], 'fixed local cube detector')
    same(row['report_to'], [0, 0, 0], 'fixed cube report geometry')
    same(row['click'], float(sum(abs((matrix@psi)[72:76])**2)), 'local cube probability')
    same(row['ledger'], expected_ledger((3, 3, 3), [0, 1, 2]), 'three-dimensional timing')


def clock_cases(rows):
    cases = list(itertools.product(('0', '1/20', '1/10', '1/5'), ('1/7', '2/7'), ('flat', 'curved')))
    need(type(rows) is list and len(rows) == len(cases), 'complete clock census')
    for row, (strength, elapsed, policy) in zip(rows, cases):
        keys(row, 'strength mass_time_pi policy matrix fringe')
        same([row['strength'], row['mass_time_pi'], row['policy']], [strength, elapsed, policy])
        phase = math.pi*float(F(elapsed)/(1+F(strength)))
        matrix = decode(row['matrix'], 16).reshape(4, 4)
        expected = np.diag([complex(math.cos(phase), -b*math.sin(phase)) for b in (1, 1, -1, -1)])
        need(np.linalg.norm(matrix-expected) < 1e-12, 'complete native clock channel')
        same(row['fringe'], math.cos(phase)**2, 'actual Ramsey read')


def symbol_cases(rows):
    cases = list(itertools.product(('0', '1/10'), ('flat', 'curved'), ('1/20', '1/40')))
    need(type(rows) is list and len(rows) == len(cases), 'complete symbol census')
    for row, (strength, policy, spacing) in zip(rows, cases):
        keys(row, 'strength policy spacing period matrix speed_relative spatial_scale remainder_bound linear_error')
        same([row['strength'], row['policy'], row['spacing']], [strength, policy, spacing])
        a = F(spacing)
        tau = 8*a
        lapse = 1/(1+F(strength))
        r = F(1, 2)*(lapse**2 if policy == 'curved' else 1)
        sine = math.sqrt(float(1-r*r))
        rotation = float(r)*I+1j*sine*BETA
        matrix = I.copy()
        linear = float(F(1, 5)*tau*lapse)*BETA
        ks = (F(1, 3), F(1, 5), F(-1, 7))
        for k, alpha in zip(ks, ALPHA):
            t = float(a*k)
            s = math.cos(t)*I-1j*math.sin(t)*alpha
            matrix = rotation@s@rotation.conj().T@s@matrix
            prime = float(r)*alpha+1j*sine*BETA@alpha
            linear += float(2*a*r*k)*prime
        mass_phase = float(F(1, 5)*tau*lapse)
        matrix = (math.cos(mass_phase)*I-1j*math.sin(mass_phase)*BETA)@matrix
        actual = decode(row['matrix'], 16).reshape(4, 4)
        need(np.linalg.norm(actual-matrix) < 2e-12, 'complete native Fourier symbol')
        same(row['period'], str(tau), 'charged continuum period')
        same(row['speed_relative'], str(2*r), 'actual transport derivative')
        same(row['spatial_scale'], str(lapse/(2*r)), 'reconstructed spatial scale')
        length = F(1, 5)*tau*lapse+2*a*sum(map(abs, ks))
        same(row['remainder_bound'], str(length*length/2), 'analytic second-derivative bound')
        err = float(np.linalg.norm(actual-I+1j*linear, ord=2))
        same(row['linear_error'], err)
        need(err <= float(length*length/2), 'finite word expansion bound')
        # Test every Clifford anticommutator, not only one dispersion direction.
        primes = [float(r)*p+1j*sine*BETA@p for p in ALPHA]
        for i, p in enumerate(primes):
            need(np.linalg.norm(p@BETA+BETA@p) < 2e-14, 'mass anticommutes with every kinetic direction')
            for j, q in enumerate(primes):
                need(np.linalg.norm(p@q+q@p-(2*I if i == j else 0)) < 2e-14, 'complete principal Clifford algebra')


def check_timelines(rows, basis_words):
    need(type(rows) is list and len(rows) == 4, 'all physical controller timelines')
    for row, axes in zip(rows, ([0], [1], [2], [0, 1, 2])):
        keys(row, 'axes service_slot stages interrogation_end receiver_record_ready returned_record_ready')
        same(row['axes'], axes)
        # Derive the chronological right-to-left factors of R S R* S
        # and B permutation B*, separately from the producer's controller.
        instructions = []
        for axis in axes:
            basis_word = basis_words[axis]
            for sign in (-1, 1):
                for gate in reversed(basis_word):
                    instructions.append(dict(kind=gate['kind'], modes=gate['modes'],
                                             angle='pi', coefficient=str(-fraction(gate['pi']))))
                instructions.append(dict(kind='flight', axis=axis))
                for gate in basis_word:
                    instructions.append(dict(kind=gate['kind'], modes=gate['modes'],
                                             angle='pi', coefficient=gate['pi']))
                for component in range(4):
                    instructions.append(dict(kind='phase', modes=[component], angle='source_acos',
                                             coefficient=str(sign*(1 if component < 2 else -1))))
        for component in range(4):
            instructions.append(dict(kind='phase', modes=[component], angle='mass_lapse',
                                     coefficient=str(-1 if component < 2 else 1)))
        count_flights = 2*len(axes)
        service = F(2, 3*(len(instructions)-count_flights)+2*count_flights)
        same(row['service_slot'], str(service), 'strictly positive load/pulse/unload slots')
        need(type(row['stages']) is list and len(row['stages']) == len(instructions), 'no hidden controller stage')
        time = F(0)
        for i, (stage, op) in enumerate(zip(row['stages'], instructions)):
            duration = 1+2*service if op['kind'] == 'flight' else 3*service
            expected = dict(stage=i, instruction=op, start=str(time), load_done=str(time+service),
                            action_done=str(time+duration-service), end=str(time+duration))
            same(stage, expected, 'ordered charged controller')
            time += duration
        need(time == count_flights+2, 'elapsed time includes every service and flight')
        same(row['interrogation_end'], str(time))
        same(row['receiver_record_ready'], str(time+1))
        same(row['returned_record_ready'], str(time+3))


def check_bands(rows):
    cases = list(itertools.product(('0', '1/10', '1/5'), ('flat', 'curved'),
                                   ('0', '1/7', '1/3'), ('0', '1/11')))
    need(type(rows) is list and len(rows) == len(cases), 'all finite bands and masses')
    for row, (strength, policy, momentum, mass) in zip(rows, cases):
        keys(row, 'strength policy momentum_pi mass_pi matrix band_cosine eigenphases')
        same([row['strength'], row['policy'], row['momentum_pi'], row['mass_pi']],
             [strength, policy, momentum, mass], 'unfiltered finite spectrum census')
        lapse = 1/(1+F(strength))
        r = float(F(1, 2)*(lapse**2 if policy == 'curved' else 1))
        q = math.sqrt(1-r*r)
        p, mu = math.pi*float(F(momentum)), math.pi*float(F(mass)*lapse)
        prime = r*ALPHA[0]+1j*q*BETA@ALPHA[0]
        w = (1-2*r*r*math.sin(p)**2)*I-1j*(
            2*r*math.cos(p)*math.sin(p)*prime+2*r*q*math.sin(p)**2*BETA)
        expected = (math.cos(mu)*I-1j*math.sin(mu)*BETA)@w
        actual = decode(row['matrix'], 16).reshape(4, 4)
        need(np.linalg.norm(actual-expected) < 2e-12, 'exact finite Clifford word')
        f = math.cos(mu)*(1-2*r*r*math.sin(p)**2)-2*r*q*math.sin(mu)*math.sin(p)**2
        need(-1 <= f <= 1, 'finite unitary band domain')
        same(row['band_cosine'], f, 'untruncated massive dispersion')
        omega = math.acos(f)
        same(row['eigenphases'], [-omega, -omega, omega, omega], 'every band multiplicity')
        need(np.linalg.norm(actual@actual-2*f*actual+I) < 3e-12, 'full matrix minimal polynomial')


def verify_evidence(packet):
    keys(packet, 'controls instruments basis_words line cube robustness clocks symbols timelines bands')
    verify_controls(packet['controls'])
    same(packet['instruments'], expected_instruments(), 'complete inherited physical instruments')
    words = verify_words(packet['basis_words'])
    check_timelines(packet['timelines'], packet['basis_words'])
    expected = list(itertools.product(('0', '1/20', '1/10', '1/5'), range(3), ('flat', 'curved')))
    need(type(packet['line']) is list and len(packet['line']) == len(expected), 'complete line experiment coverage')
    for row, case in zip(packet['line'], expected):
        line_case(row, case, words)
    cubes = list(itertools.product(('1/10', '1/5'), ('flat', 'curved'), ('0', '1/7')))
    need(type(packet['cube']) is list and len(packet['cube']) == len(cubes), 'complete cube experiment coverage')
    for row, case in zip(packet['cube'], cubes):
        cube_case(row, case, words)
    robust = []
    for strength in ('1/20', '1/10', '1/5'):
        gap = F(1, 4)*(1-1/(1+F(strength))**4)
        robust.append(dict(strength=strength, exact_gap=str(gap), per_branch_total_error=str(gap/8),
                           remaining_separation=str(3*gap/4)))
    same(packet['robustness'], robust, 'exact separated error intervals')
    clock_cases(packet['clocks'])
    symbol_cases(packet['symbols'])
    check_bands(packet['bands'])
    return True

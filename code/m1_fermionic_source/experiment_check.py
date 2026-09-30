"""Independent full-instrument, interaction, conversion and noise controls."""

import itertools
import math

import numpy as np

from .check import X, Z, annihilator, at, close, complex_array, exact, exterior, keys, need, word


def replay_native(tape, n):
    need(type(tape) is list and len(tape) > 0, 'nonempty native instrument')
    h = (X+Z)/math.sqrt(2)
    result = np.eye(1 << n, dtype=complex)
    for gate in tape:
        need(type(gate) is list and len(gate) in (2, 3) and type(gate[1]) is int and 0 <= gate[1] < n, 'native instruction')
        kind, q, *rest = gate
        if kind == 'cz':
            need(len(rest) == 1 and type(rest[0]) is int and 0 <= rest[0] < n and q != rest[0], 'CZ instruction')
            u = np.diag([1-2*((s >> q) & 1)*((s >> rest[0]) & 1) for s in range(1 << n)])
        else:
            need(not rest and kind in ('x', 'h', 's', 'sdg'), 'Clifford instruction')
            u = at({'x': X, 'h': h, 's': np.diag([1., 1j]), 'sdg': np.diag([1., -1j])}[kind], q, n)
        result = u@result
    return result


def verify_measurements(rows):
    challenges = [(3, ['3', '6', 1]), (3, ['0', '7', 0]), (1, ['0', '1', 2])]
    need(type(rows) is list and len(rows) == len(challenges), 'all QND challenges')
    for row, (n, raw) in zip(rows, challenges):
        keys(row, 'qubits word tape unitary pulses pulse_code_events')
        exact([row['qubits'], row['word']], [n, raw])
        actual = replay_native(row['tape'], n+1)
        p = word(raw, n)
        identity = np.eye(1 << n)
        plus, minus = (identity+p)/2, (identity-p)/2
        expected = np.kron(np.eye(2), plus)+np.kron(X, minus)
        close(actual, expected, 'complete QND dilation')
        close(complex_array(row['unitary'], actual.shape), expected, 'actual native QND execution')
        # These are the complete, trace-preserving retained-outcome Kraus maps.
        close(actual[:1 << n, :1 << n], plus, 'positive read branch')
        close(actual[1 << n:, :1 << n], minus, 'negative read branch')
        close(plus.conj().T@plus+minus.conj().T@minus, identity, 'both outcomes retained')
        need(len(row['tape']) <= 7*n+1, 'bounded QND pulse count')
        exact(row['pulses'], len(row['tape']))
        exact(row['pulse_code_events'], 2*len(row['tape'])+2)


def graph_operators(n, edges, parity):
    m = len(edges)
    eye = np.eye(1 << m)
    b = []
    for v in range(n):
        value = eye*(-1 if parity and v == 0 else 1)
        for e, ends in enumerate(edges):
            if v in ends:
                value = value@at(Z, e, m)
        b.append(value)
    a = []
    for k, (i, j) in enumerate(edges):
        value = at(X, k, m)
        for v in (i, j):
            for e, ends in enumerate(edges[:k]):
                if v in ends:
                    value = value@at(Z, e, m)
        a.append(value)
    return b, a


def verify_conversions(rows):
    expected = [('add_vertex', 2, [[0, 1]], 3, [[0, 1], [1, 2]]),
                ('add_cycle', 3, [[0, 1], [1, 2]], 3, [[0, 1], [1, 2], [0, 2]])]
    need(type(rows) is list and len(rows) == 2, 'conversion catalogue')
    for row, (kind, n0, e0, n1, e1) in zip(rows, expected):
        keys(row, 'kind old_vertices old_edges new_vertices new_edges corrected_branches')
        exact([row[k] for k in ('kind', 'old_vertices', 'old_edges', 'new_vertices', 'new_edges')],
              [kind, n0, e0, n1, e1])
        m0, m1 = len(e0), len(e1)
        append = np.kron(np.array([[1.], [0.]]), np.eye(1 << m0))
        if kind == 'add_vertex':
            expected_maps = [append]
        else:
            _, a = graph_operators(n1, e1, 0)
            s = 1j**3*a[0]@a[1]@(-a[2])
            expected_maps = [(np.eye(1 << m1)+s)@append/2,
                             at(Z, m0, m1)@(np.eye(1 << m1)-s)@append/2]
        need(type(row['corrected_branches']) is list and len(row['corrected_branches']) == len(expected_maps), 'all conversion outcomes')
        maps = [complex_array(raw, (1 << m1, 1 << m0)) for raw in row['corrected_branches']]
        for actual, expected_map in zip(maps, expected_maps):
            close(actual, expected_map, 'physical conversion instrument')
            close(actual.conj().T@actual, np.eye(1 << m0)/len(maps), 'input-independent branch probability')
            close(actual, maps[0], 'every corrected branch gives the same isometry')
        for parity in (0, 1):
            b0, a0 = graph_operators(n0, e0, parity)
            b1, a1 = graph_operators(n1, e1, parity)
            for left, right in zip(b0+a0, b1[:n0]+a1[:len(a0)]):
                close(right@maps[0], maps[0]@left, 'unknown-state conversion intertwining')
            if kind == 'add_vertex':
                close(b1[-1]@maps[0], maps[0], 'added fermion mode is empty')


def replay_rotations(circuit, m, loops, counts=(1, 6), allowed=None):
    keys(circuit, 'scalar rotations')
    scalar = complex_array(circuit['scalar'], ())
    need(abs(abs(scalar)-1) < 1e-12, 'scalar phase')
    need(type(circuit['rotations']) is list and len(circuit['rotations']) in counts, 'detector rotation count')
    out = np.eye(1 << m, dtype=complex)
    for row in circuit['rotations']:
        keys(row, 'word theta')
        p = word(row['word'], m)
        if allowed is not None:
            need((int(row['word'][0]) | int(row['word'][1])) & ~allowed == 0, 'local circuit support')
        close(p, p.conj().T, 'Hermitian readout generator')
        for loop in loops:
            close(p@loop, loop@p, 'detector preserves code')
        need(type(row['theta']) is float and math.isfinite(row['theta']), 'finite readout angle')
        t = row['theta']
        out = (math.cos(t)*np.eye(len(out))+1j*math.sin(t)*p)@out
    return scalar*out


def verify_detector(rows, cache):
    g = cache['complete_four']
    need(type(rows) is list and len(rows) == 2, 'both detector parity blocks')
    one_particle = np.zeros((4, 4), complex)
    theta, weights = .47, [.37, .81]
    for i, j, w in ((0, 2, weights[0]), (1, 3, weights[1])):
        f = np.zeros(4, complex)
        f[i], f[j] = 1/math.sqrt(2), np.exp(1j*theta)/math.sqrt(2)
        one_particle += w*np.outer(f, f.conj())
    expected = np.eye(16)-exterior(np.eye(4)-one_particle)
    for parity, row in enumerate(rows):
        keys(row, 'parity theta weights circuits code_effect')
        exact([row['parity'], row['theta'], row['weights']], [parity, theta, weights])
        need(type(row['circuits']) is list and len(row['circuits']) == 4, 'every detector phase and mixer')
        unitary = np.eye(1 << g['m'], dtype=complex)
        for k, circuit in enumerate(row['circuits']):
            pair = ((0, 2), (1, 3))[k//2]
            modes = [pair[1]] if k % 2 == 0 else pair
            allowed = sum(1 << e for e, ends in enumerate(g['edges']) if any(v in ends for v in modes))
            unitary = replay_rotations(circuit, g['m'], g['loops'], ((1,) if k % 2 == 0 else (6,)), allowed)@unitary
        # Independent no-click effect is the product of commuting local
        # occupation factors, including all occupied configurations.
        no_click = np.eye(len(unitary))
        for mode, weight in enumerate(weights):
            b = (-1 if parity and mode == 0 else 1)*g['b'][mode]
            occupation = (np.eye(len(unitary))-b)/2
            no_click = no_click@(np.eye(len(unitary))-weight*occupation)
        effect = np.eye(len(unitary))-unitary.conj().T@no_click@unitary
        iso, basis = g['isometries'][parity], g['bases'][parity]
        target = expected[np.ix_(basis, basis)]
        close(iso.conj().T@effect@iso, target, 'full fermionic click effect')
        close(complex_array(row['code_effect'], target.shape), target, 'executed serial weighted detector')
    need(np.linalg.eigvalsh(expected).min() > -1e-12 and np.linalg.eigvalsh(expected).max() <= 1+1e-12, 'bounded all-sector effect')


def verify_interaction(row, cache):
    keys(row, 'parity angle modes initial encoded_output circuit')
    exact([row['parity'], row['angle'], row['modes']], [0, float(math.pi), [1, 3]])
    before = np.zeros(16, complex)
    before[[5, 6, 9, 10]] = .5
    close(complex_array(row['initial'], (16,)), before, 'two occupied dual rails')
    after = before.copy()
    after[10] *= -1
    g = cache['complete_four']
    iso, basis = g['isometries'][0], g['bases'][0]
    emitted = complex_array(row['encoded_output'], (1 << g['m'],))
    allowed = sum(1 << e for e, ends in enumerate(g['edges']) if 1 in ends or 3 in ends)
    actual = replay_rotations(row['circuit'], g['m'], g['loops'], (3,), allowed)
    target_gate = np.diag([(-1.)**(((s >> 1) & 1)*((s >> 3) & 1)) for s in basis])
    close(actual@iso, iso@target_gate, 'entire interacting Fock operation')
    close(actual@iso@before[basis], emitted, 'executed density circuit')
    close(emitted, iso@after[basis], 'actual non-Gaussian density operation')
    a = [annihilator(4, i) for i in range(4)]
    number = [v.conj().T@v for v in a]
    left_z, right_z = number[0]-number[1], number[2]-number[3]
    left_x = a[0].conj().T@a[1]+a[1].conj().T@a[0]
    right_x = a[2].conj().T@a[3]+a[3].conj().T@a[2]
    chsh = math.sqrt(2)*(left_z@right_x+left_x@right_z)
    need(abs(np.vdot(after, chsh@after).real-2*math.sqrt(2)) < 1e-12, 'interacting dual-rail CHSH')
    gamma = np.array([[np.vdot(after, a[i].conj().T@a[j]@after) for j in range(4)] for i in range(4)])
    close(gamma, np.eye(4)/2, 'non-Slater one-body density')


def verify_noise(rows, cache):
    names = ['triangle', 'two_cycles', 'complete_four']
    need(type(rows) is list and len(rows) == len(names), 'noise catalogue')
    for row, name in zip(rows, names):
        keys(row, 'name exposure phase_flip_probability cycle_rank chord_syndromes vacuum_fidelity failure_lower')
        g = cache[name]
        rank, m = len(g['loops']), g['m']
        exact([row['name'], row['exposure'], row['cycle_rank']], [name, .2, rank])
        exact(row['chord_syndromes'], np.eye(rank, dtype=int).tolist())
        p = -math.expm1(-.2)/2
        exact(row['phase_flip_probability'], p)
        # The diagonal channel's action on the vacuum density is evaluated
        # entrywise from Hamming distance, not by the producer's error sum.
        vacuum = g['vacuum']
        rho = np.outer(vacuum, vacuum.conj())
        for i, j in itertools.product(range(1 << m), repeat=2):
            rho[i, j] *= math.exp(-.2*(i ^ j).bit_count())
        fidelity = float(np.vdot(vacuum, rho@vacuum).real)
        exact(row['vacuum_fidelity'], fidelity)
        lower = 1-(1-p)**rank
        exact(row['failure_lower'], lower)
        need(1-fidelity >= lower-1e-12 and fidelity < .9, 'physical vacuum decoherence')


def verify_resolved_reads(rows, cache):
    g = cache['cycle_with_tail']
    n, m = g['n'], g['m']
    targets = [('complex', np.array([1., .3j, 0., -.4+.2j, .1])),
               ('zero_subtree', np.array([1., .4j, 0., 0., 0.]))]
    need(type(rows) is list and len(rows) == 2, 'resolved read catalogue')
    for row, (name, f) in zip(rows, targets):
        keys(row, 'name kernel blocks forward_rotations qnd_rotations_with_return')
        exact([row['name'], row['forward_rotations'], row['qnd_rotations_with_return']], [name, 24, 48])
        f /= np.linalg.norm(f)
        close(complex_array(row['kernel'], (n,)), f, 'complex spatial read kernel')
        annih = sum(f[j].conjugate()*annihilator(n, j) for j in range(n))
        target = annih.conj().T@annih
        need(type(row['blocks']) is list and len(row['blocks']) == 2, 'all read parity blocks')
        for parity, block in enumerate(row['blocks']):
            keys(block, 'parity circuits code_effect')
            exact(block['parity'], parity)
            need(type(block['circuits']) is list and len(block['circuits']) == n-1, 'zero subtrees also charged')
            u = np.eye(1 << m, dtype=complex)
            single = np.eye(n, dtype=complex)
            for step in block['circuits']:
                keys(step, 'modes matrix scalar rotations')
                modes = step['modes']
                need(type(modes) is list and len(modes) == 2 and all(type(i) is int for i in modes)
                     and sorted(modes) in g['edges'], 'local tree mixing')
                two = complex_array(step['matrix'], (2, 2))
                close(two.conj().T@two, np.eye(2), 'tree SU2')
                need(abs(np.linalg.det(two)-1) < 1e-12, 'tree determinant')
                single[modes] = two@single[modes]
                allowed = sum(1 << e for e, ends in enumerate(g['edges']) if any(v in ends for v in modes))
                u = replay_rotations({k: step[k] for k in ('scalar', 'rotations')}, m, g['loops'], (6,), allowed)@u
            close(single@f, np.eye(n)[0], 'tree really gathers the requested mode')
            iso, basis = g['isometries'][parity], g['bases'][parity]
            close(u@iso, iso@exterior(single)[np.ix_(basis, basis)], 'complete tree Fock implementation')
            occupation = (np.eye(len(u))-(-1 if parity else 1)*g['b'][0])/2
            effect = iso.conj().T@u.conj().T@occupation@u@iso
            close(effect, target[np.ix_(basis, basis)], 'actual CAR number read')
            close(complex_array(block['code_effect'], effect.shape), effect, 'emitted CAR read')
            close(effect@effect, effect, 'repeatable number instrument after returning the frame')

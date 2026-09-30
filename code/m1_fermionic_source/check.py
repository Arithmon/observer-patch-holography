"""Independent occupation-basis, tensor-Pauli and circuit verification.

No producer module is imported. Fermion signs come directly from the
occupation rule and many-particle unitaries from determinants of minors.
"""

import itertools
import math

import numpy as np


def need(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, fields):
    need(type(value) is dict and set(value) == set(fields.split()), 'schema: '+fields)


def exact(actual, expected):
    need(type(actual) is type(expected), 'JSON type')
    if isinstance(expected, dict):
        need(actual.keys() == expected.keys(), 'keys')
        for k in expected:
            exact(actual[k], expected[k])
    elif isinstance(expected, list):
        need(len(actual) == len(expected), 'list length')
        for a, e in zip(actual, expected):
            exact(a, e)
    elif isinstance(expected, float):
        need(math.isfinite(actual) and abs(actual-expected) <= (2e-10*abs(expected) if expected else 1e-13), 'scalar')
    else:
        need(actual == expected, 'value')


def complex_array(raw, shape):
    def floats(x):
        return all(floats(t) for t in x) if type(x) is list else type(x) is float and math.isfinite(x)
    need(type(raw) is list and floats(raw), 'finite JSON float pairs')
    a = np.asarray(raw)
    need(a.shape == shape+(2,), 'complex shape')
    return a[..., 0]+1j*a[..., 1]


def close(a, b, label):
    need(np.shape(a) == np.shape(b) and np.max(np.abs(a-b)) < 2e-10, label)


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def at(matrix, bit, size):
    out = np.array([[1.]], complex)
    for j in reversed(range(size)):
        out = np.kron(out, matrix if j == bit else I2)
    return out


def word(raw, qubits):
    need(type(raw) is list and len(raw) == 3, 'Pauli word')
    x, z, phase = raw
    need(all(type(s) is str and s.isascii() and s.isdigit() and str(int(s)) == s for s in (x, z)), 'canonical masks')
    need(type(phase) is int and 0 <= phase < 4, 'Pauli phase')
    x, z = int(x), int(z)
    need(x < 1 << qubits and z < 1 << qubits, 'Pauli support')
    out = np.eye(1 << qubits, dtype=complex)*1j**phase
    for j in range(qubits):
        if x >> j & 1:
            out = out@at(X, j, qubits)
    for j in range(qubits):
        if z >> j & 1:
            out = out@at(Z, j, qubits)
    return out


def annihilator(n, mode):
    out = np.zeros((1 << n, 1 << n), complex)
    for occupied in range(1 << n):
        if occupied >> mode & 1:
            lower = sum((occupied >> k) & 1 for k in range(mode))
            out[occupied ^ (1 << mode), occupied] = (-1)**lower
    return out


def exterior(u):
    n = len(u)
    occupied = [[j for j in range(n) if bits >> j & 1] for bits in range(1 << n)]
    result = np.zeros((1 << n, 1 << n), complex)
    for r, a in enumerate(occupied):
        for s, b in enumerate(occupied):
            if len(a) == len(b):
                result[r, s] = np.linalg.det(u[np.ix_(a, b)]) if a else 1.
    return result


def graph_catalog():
    return [('triangle', 3, [[0, 1], [0, 2], [1, 2]]),
            ('square', 4, [[0, 1], [0, 3], [1, 2], [2, 3]]),
            ('two_cycles', 4, [[0, 1], [0, 2], [0, 3], [1, 2], [2, 3]]),
            ('complete_four', 4, [[i, j] for i in range(4) for j in range(i+1, 4)]),
            ('cycle_with_tail', 5, [[0, 1], [0, 2], [1, 2], [2, 3], [3, 4]])]


def verify_small(rows):
    need(type(rows) is list and len(rows) == len(graph_catalog()), 'graph coverage')
    cache = {}
    for row, (name, n, edges) in zip(rows, graph_catalog()):
        keys(row, 'name vertices edges parent chords cycles loops edge_operators branches blocks')
        exact([row['name'], row['vertices'], row['edges']], [name, n, edges])
        m = len(edges)
        identity = np.eye(1 << m)
        parent = row['parent']
        need(type(parent) is list and len(parent) == n and parent[0] is None, 'rooted tree')
        for i in range(1, n):
            seen, v = set(), i
            while v != 0:
                need(v not in seen and type(parent[v]) is int and 0 <= parent[v] < n, 'tree cycle')
                seen.add(v)
                need(sorted([v, parent[v]]) in edges, 'tree edge')
                v = parent[v]
        tree = [sorted([i, parent[i]]) for i in range(1, n)]
        chords = [k for k, e in enumerate(edges) if e not in tree]
        exact(row['chords'], chords)
        need(type(row['cycles']) is list and type(row['loops']) is list and
             len(row['cycles']) == len(row['loops']) == m-n+1, 'all fundamental cycles')
        bq = []
        for v in range(n):
            b = identity.copy()
            for e, pair in enumerate(edges):
                if v in pair:
                    b = b@at(Z, e, m)
            bq.append(b)
        need(type(row['edge_operators']) is list and len(row['edge_operators']) == m, 'all edge operators')
        aq = []
        for e, (i, j) in enumerate(edges):
            a = at(X, e, m)
            for t in range(e):
                if i in edges[t] or j in edges[t]:
                    a = a@at(Z, t, m)
            close(word(row['edge_operators'][e], m), a, 'edge ordering signs')
            aq.append(a)
        loops = []
        for chord, cycle, raw in zip(chords, row['cycles'], row['loops']):
            need(type(cycle) is list and len(cycle) >= 4 and cycle[0] == cycle[-1], 'closed simple cycle')
            need(all(type(v) is int and 0 <= v < n for v in cycle), 'cycle vertices')
            need(len(set(cycle[:-1])) == len(cycle)-1, 'simple cycle')
            used = [sorted([a, b]) for a, b in zip(cycle, cycle[1:])]
            need(edges[chord] in used and all(e in tree or e == edges[chord] for e in used), 'fundamental cycle')
            s = identity.astype(complex)*1j**len(used)
            for a, b in zip(cycle, cycle[1:]):
                s = s@((1 if a < b else -1)*aq[edges.index(sorted([a, b]))])
            close(word(raw, m), s, 'signed loop relation')
            close(s@s, identity, 'loop involution')
            for operator in bq+aq+loops:
                close(s@operator, operator@s, 'loop commutation')
            loops.append(s)
        branches = row['branches']
        outcomes = list(itertools.product((0, 1), repeat=len(loops)))
        need(type(branches) is list and len(branches) == len(outcomes), 'all preparation branches')
        for branch, bits in zip(branches, outcomes):
            keys(branch, 'outcomes probability state')
            exact(branch['outcomes'], list(bits))
            state = np.eye(1 << m, dtype=complex)[:, 0]
            for bit, loop in zip(bits, loops):
                state = (identity+(-1)**bit*loop)@state/2
            probability = float(np.vdot(state, state).real)
            exact(branch['probability'], float(2.**(-len(loops))))
            need(abs(probability-branch['probability']) < 1e-12, 'actual outcome probability')
            state /= np.sqrt(probability)
            for bit, e in zip(bits, chords):
                if bit:
                    state = at(Z, e, m)@state
            emitted = complex_array(branch['state'], (1 << m,))
            close(emitted, state, 'complete corrected branch')
            for stabilizer in loops+bq:
                close(stabilizer@emitted, emitted, 'prepared vacuum constraints')
        need(type(row['blocks']) is list and len(row['blocks']) == 2, 'both parity blocks')
        annih = [annihilator(n, i) for i in range(n)]
        majorana = [a+a.conj().T for a in annih]
        bf = [np.eye(1 << n)-2*a.conj().T@a for a in annih]
        af = [-1j*majorana[i]@majorana[j] for i, j in edges]
        cache[name] = dict(n=n, m=m, edges=edges, b=bq, a=aq, loops=loops,
                           vacuum=complex_array(branches[0]['state'], (1 << m,)),
                           isometries=[], bases=[])
        for parity, block in enumerate(row['blocks']):
            keys(block, 'parity basis isometry gates')
            exact(block['parity'], parity)
            basis = [x for x in range(1 << n) if sum((x >> k) & 1 for k in range(n)) % 2 == parity]
            exact(block['basis'], basis)
            iso = complex_array(block['isometry'], (1 << m, len(basis)))
            cache[name]['isometries'].append(iso)
            cache[name]['bases'].append(basis)
            close(iso.conj().T@iso, np.eye(len(basis)), 'isometry on entire parity sector')
            for v in range(n):
                close((-1 if v == 0 and parity else 1)*bq[v]@iso,
                      iso@bf[v][np.ix_(basis, basis)], 'all occupation operators')
            for a, f in zip(aq, af):
                close(a@iso, iso@f[np.ix_(basis, basis)], 'all Majorana edges')
            for s in loops:
                close(s@iso, iso, 'full code image')
            verify_gates(block['gates'], edges[-1], n, m, iso, basis, loops, edges)
    return cache


def verify_gates(rows, pair, n, m, iso, basis, loops, edges):
    need(type(rows) is list and len(rows) == 4, 'gate coverage')
    i, j = pair
    g = np.array([[np.exp(.23j)*np.cos(.37), np.exp(-.51j)*np.sin(.37)],
                  [-np.exp(.51j)*np.sin(.37), np.exp(-.23j)*np.cos(.37)]])
    for row, kind, count in zip(rows, ('phase', 'mix', 'swap', 'density'), (1, 6, 8, 3)):
        keys(row, 'kind modes parameter scalar rotations')
        exact(row['kind'], kind)
        exact(row['modes'], [i] if kind == 'phase' else [i, j])
        single = np.eye(n, dtype=complex)
        if kind == 'phase':
            exact(row['parameter'], .41)
            single[i, i] = np.exp(.41j)
        elif kind == 'mix':
            close(complex_array(row['parameter'], (2, 2)), g, 'fixed complex mixing challenge')
            single[np.ix_([i, j], [i, j])] = g
        elif kind == 'swap':
            exact(row['parameter'], None)
            single[[i, j]] = single[[j, i]]
        else:
            exact(row['parameter'], .73)
        target = (np.diag([np.exp(-.73j*((s >> i) & 1)*((s >> j) & 1)) for s in range(1 << n)])
                  if kind == 'density' else exterior(single))
        need(type(row['rotations']) is list and len(row['rotations']) == count, 'rotation count')
        u = np.eye(1 << m, dtype=complex)
        allowed = sum(1 << e for e, ends in enumerate(edges) if any(v in ends for v in row['modes']))
        for step in row['rotations']:
            keys(step, 'word theta')
            p = word(step['word'], m)
            need((int(step['word'][0]) | int(step['word'][1])) & ~allowed == 0, 'local gate support')
            need(type(step['theta']) is float and math.isfinite(step['theta']), 'finite rotation angle')
            close(p, p.conj().T, 'Hermitian generator')
            close(p@p, np.eye(1 << m), 'involution')
            for loop in loops:
                close(loop@p, p@loop, 'code-preserving generator')
            theta = step['theta']
            u = (math.cos(theta)*np.eye(1 << m)+1j*math.sin(theta)*p)@u
        scalar = complex_array(row['scalar'], ())
        need(abs(abs(scalar)-1) < 1e-12, 'global phase only')
        close(scalar*u@iso, iso@target[np.ix_(basis, basis)], 'full Fock gate and no leakage')


def verify_gadgets(rows):
    cases = [(2, ['1', '2', 0], .37), (3, ['3', '6', 1], -.51),
             (3, ['5', '7', 0], .29), (1, ['0', '1', 2], .63)]
    need(type(rows) is list and len(rows) == len(cases), 'native gadget coverage')
    h = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    for row, (n, raw, theta) in zip(rows, cases):
        keys(row, 'qubits word theta tape unitary pulses pulse_code_events')
        exact([row['qubits'], row['word'], row['theta']], [n, raw, theta])
        need(type(row['tape']) is list and 0 < len(row['tape']) <= 10*n+1, 'bounded pulse tape')
        actual = np.eye(1 << (n+1), dtype=complex)
        for gate in row['tape']:
            need(type(gate) is list and len(gate) >= 2 and type(gate[1]) is int and 0 <= gate[1] <= n, 'native pulse')
            kind, q, *rest = gate
            if kind == 'cz':
                need(len(rest) == 1 and type(rest[0]) is int and 0 <= rest[0] <= n and q != rest[0], 'CZ domain')
                event = np.diag([(-1.)**(((s >> q) & 1)*((s >> rest[0]) & 1)) for s in range(len(actual))])
            elif kind == 'rz':
                need(len(rest) == 1 and type(rest[0]) is float and math.isfinite(rest[0]), 'rotation domain')
                event = at(np.diag([np.exp(1j*rest[0]), np.exp(-1j*rest[0])]), q, n+1)
            else:
                need(not rest and kind in ('h', 's', 'sdg'), 'Clifford pulse domain')
                event = at({'h': h, 's': np.diag([1., 1j]), 'sdg': np.diag([1., -1j])}[kind], q, n+1)
            actual = event@actual
        expected = math.cos(theta)*np.eye(len(actual))+1j*math.sin(theta)*np.kron(Z, word(raw, n))
        close(actual, expected, 'full ancilla circuit, including uncomputation')
        close(complex_array(row['unitary'], actual.shape), expected, 'native M6 execution')
        exact(row['pulses'], len(row['tape']))
        exact(row['pulse_code_events'], 2*len(row['tape']))

"""Edge-code producer, deterministic preparation and fermionic gate tapes.

The encoding is Bravyi--Kitaev's superfast encoding, not an OPH invention.
The independent checker constructs fermions directly from occupation signs.
"""

from collections import deque
import itertools

import numpy as np

from .pauli import Word, rotation


CATALOG = (
    ('triangle', 3, ((0, 1), (0, 2), (1, 2))),
    ('square', 4, ((0, 1), (0, 3), (1, 2), (2, 3))),
    ('two_cycles', 4, ((0, 1), (0, 2), (0, 3), (1, 2), (2, 3))),
    ('complete_four', 4, tuple(itertools.combinations(range(4), 2))),
    ('cycle_with_tail', 5, ((0, 1), (0, 2), (1, 2), (2, 3), (3, 4))),
)


class Encoding:
    def __init__(self, vertices, edges, orders=None):
        if type(vertices) is not int or vertices < 2:
            raise ValueError('at least two vertices required')
        if any(len(e) != 2 or any(type(v) is not int or not 0 <= v < vertices for v in e)
               or e[0] >= e[1] for e in edges):
            raise ValueError('simple increasing graph edges required')
        self.n = vertices
        # Edge order is part of the encoding convention. Appending an edge
        # leaves every existing local order unchanged during code extension.
        self.edges = tuple(tuple(e) for e in edges)
        if len(set(self.edges)) != len(self.edges):
            raise ValueError('duplicate edge')
        self.m = len(self.edges)
        self.index = {e: k for k, e in enumerate(self.edges)}
        self.neighbors = [[] for _ in range(vertices)]
        self.incident = [[] for _ in range(vertices)]
        for k, (i, j) in enumerate(self.edges):
            self.neighbors[i].append(j)
            self.neighbors[j].append(i)
            self.incident[i].append(k)
            self.incident[j].append(k)
        if orders is not None:
            if len(orders) != vertices or any(any(type(k) is not int for k in order) or
                    sorted(order) != expected for order, expected in zip(orders, self.incident)):
                raise ValueError('each local order must permute its incident edges')
            self.incident = [list(order) for order in orders]
        self.parent = {0: None}
        queue = deque([0])
        while queue:
            i = queue.popleft()
            for j in self.neighbors[i]:
                if j not in self.parent:
                    self.parent[j] = i
                    queue.append(j)
        if len(self.parent) != vertices:
            raise ValueError('connected graph required')
        tree = {tuple(sorted((i, j))) for i, j in self.parent.items() if j is not None}
        self.chords = [k for k, e in enumerate(self.edges) if e not in tree]
        self.cycles = [self.path(*self.edges[k])+[self.edges[k][0]] for k in self.chords]
        self.loops = [self.loop(path) for path in self.cycles]

    def path(self, i, j):
        left, right = [i], [j]
        while self.parent[left[-1]] is not None:
            left.append(self.parent[left[-1]])
        while self.parent[right[-1]] is not None:
            right.append(self.parent[right[-1]])
        while len(left) > 1 and len(right) > 1 and left[-2] == right[-2]:
            left.pop()
            right.pop()
        return left + list(reversed(right[:-1]))

    def b(self, i, parity=0):
        return Word(0, sum(1 << k for k in self.incident[i]), 2*parity if i == 0 else 0)

    def a(self, i, j):
        k = self.index[tuple(sorted((i, j)))]
        z = sum(1 << t for v in (i, j) for t in self.incident[v][:self.incident[v].index(k)])
        return Word(1 << k, z, 0 if i < j else 2)

    def loop(self, cycle):
        result = Word(phase=(len(cycle)-1) % 4)
        for i, j in zip(cycle, cycle[1:]):
            result = result*self.a(i, j)
        return result

    def prepare(self, outcomes=None):
        """Retain every branch; negative outcomes get a single chord-Z repair."""
        if outcomes is None:
            outcomes = [0]*len(self.loops)
        if len(outcomes) != len(self.loops) or any(type(x) is not int or x not in (0, 1) for x in outcomes):
            raise ValueError('one binary outcome per fundamental cycle required')
        state = np.eye(1 << self.m, dtype=complex)[:, 0]
        probability = 1.
        for bit, loop in zip(outcomes, self.loops):
            state = (state+(-1)**bit*loop.matrix(self.m)@state)/2
            p = float(np.vdot(state, state).real)
            probability *= p
            state /= np.sqrt(p)
        for bit, chord in zip(outcomes, self.chords):
            if bit:
                state = Word(z=1 << chord).matrix(self.m)@state
        return state, probability

    def isometry(self, parity):
        if type(parity) is not int or parity not in (0, 1):
            raise ValueError('parity must be zero or one')
        seed, _ = self.prepare()
        basis = [j for j in range(1 << self.n) if j.bit_count() % 2 == parity]
        columns = []
        for bits in basis:
            physical, fermionic = Word(), Word()
            for v in range(1, self.n):
                if bits >> v & 1:
                    path = self.path(0, v)
                    for i, j in zip(path, path[1:]):
                        physical = self.a(i, j)*physical
                        gi, gj = Word(1 << i, (1 << i)-1), Word(1 << j, (1 << j)-1)
                        fermionic = (gi*gj).times_i(3)*fermionic
            input_bits = parity
            if input_bits ^ fermionic.x != bits:
                raise ValueError('occupation path failed')
            phase = 1j**fermionic.phase*(-1)**((fermionic.z & input_bits).bit_count())
            columns.append(physical.matrix(self.m)@seed/phase)
        return basis, np.array(columns).T

    def phase(self, i, angle, parity=0):
        return np.exp(.5j*angle), [(self.b(i, parity), -angle/2)]

    def su2(self, i, j, g, parity=0):
        g = np.asarray(g, complex)
        if g.shape != (2, 2) or not np.all(np.isfinite(g)) or np.linalg.norm(g.conj().T@g-np.eye(2)) > 1e-12 or abs(np.linalg.det(g)-1) > 1e-12:
            raise ValueError('SU2 matrix required')
        a, b = g[0]
        theta = np.arctan2(abs(b), abs(a))
        xi, zeta = (np.angle(a)+np.angle(b))/2, (np.angle(a)-np.angle(b))/2
        bi, bj, edge = self.b(i, parity), self.b(j, parity), self.a(i, j)
        # D(xi) exp(i theta sigma_y) D(zeta), rightmost gate first.
        return 1.+0j, [(bi, -zeta/2), (bj, zeta/2), (edge, theta/2),
                      (edge*bi*bj, -theta/2), (bi, -xi/2), (bj, xi/2)]

    def fswap(self, i, j, parity=0):
        scalar, tape = self.su2(i, j, -1j*np.array([[0, 1], [1, 0]]), parity)
        for v in (i, j):
            phase, part = self.phase(v, np.pi/2, parity)
            scalar *= phase
            tape += part
        return scalar, tape

    def density(self, i, j, angle, parity=0):
        bi, bj = self.b(i, parity), self.b(j, parity)
        return np.exp(-.25j*angle), [(bi, angle/4), (bj, angle/4), (bi*bj, -angle/4)]

    def execute(self, circuit):
        scalar, tape = circuit
        result = np.eye(1 << self.m, dtype=complex)
        for word, theta in tape:
            result = rotation(word, theta, self.m)@result
        return scalar*result


def encode(value):
    z = np.asarray(value, complex)
    return np.stack((z.real, z.imag), axis=-1).tolist()


def gate_cases(enc, parity):
    i, j = enc.edges[-1]
    g = np.array([[np.exp(.23j)*np.cos(.37), np.exp(-.51j)*np.sin(.37)],
                  [-np.exp(.51j)*np.sin(.37), np.exp(-.23j)*np.cos(.37)]])
    return [('phase', [i], .41, enc.phase(i, .41, parity)),
            ('mix', [i, j], encode(g), enc.su2(i, j, g, parity)),
            ('swap', [i, j], None, enc.fswap(i, j, parity)),
            ('density', [i, j], .73, enc.density(i, j, .73, parity))]


def small_graphs():
    rows = []
    for name, n, edges in CATALOG:
        enc = Encoding(n, edges)
        branches = []
        for outcomes in itertools.product((0, 1), repeat=len(enc.loops)):
            state, probability = enc.prepare(list(outcomes))
            branches.append(dict(outcomes=list(outcomes), probability=probability, state=encode(state)))
        blocks = []
        for parity in (0, 1):
            basis, iso = enc.isometry(parity)
            gates = []
            for kind, modes, parameter, (scalar, tape) in gate_cases(enc, parity):
                gates.append(dict(kind=kind, modes=modes, parameter=parameter, scalar=encode(scalar),
                                  rotations=[dict(word=w.row(), theta=float(t)) for w, t in tape]))
            blocks.append(dict(parity=parity, basis=basis, isometry=encode(iso), gates=gates))
        rows.append(dict(name=name, vertices=n, edges=[list(e) for e in enc.edges],
                         parent=[enc.parent[i] for i in range(n)], chords=enc.chords,
                         cycles=enc.cycles, loops=[w.row() for w in enc.loops],
                         edge_operators=[enc.a(i, j).row() for i, j in enc.edges],
                         branches=branches, blocks=blocks))
    return rows

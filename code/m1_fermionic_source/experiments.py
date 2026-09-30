"""Code conversion, actual detector, quartic interaction and vacuum faults."""

import itertools

import numpy as np

from .model import Encoding, encode
from .pauli import Word
from .circuits import execute_weighted, weighted_program


def conversions():
    cases = [('add_vertex', 2, [(0, 1)], 3, [(0, 1), (1, 2)]),
             ('add_cycle', 3, [(0, 1), (1, 2)], 3, [(0, 1), (1, 2), (0, 2)])]
    rows = []
    for kind, n0, e0, n1, e1 in cases:
        old, new = Encoding(n0, e0), Encoding(n1, e1)
        append_zero = np.vstack((np.eye(1 << old.m), np.zeros((1 << old.m, 1 << old.m))))
        if kind == 'add_vertex':
            maps = [append_zero.astype(complex)]
        else:
            cycle = new.loops[-1].matrix(new.m)
            identity = np.eye(1 << new.m)
            correction = Word(z=1 << old.m).matrix(new.m)
            maps = [(identity+cycle)@append_zero/2,
                    correction@(identity-cycle)@append_zero/2]
        rows.append(dict(kind=kind, old_vertices=n0, old_edges=[list(e) for e in e0],
                         new_vertices=n1, new_edges=[list(e) for e in e1],
                         corrected_branches=[encode(k) for k in maps]))
    return rows


def detector():
    enc = Encoding(4, list(itertools.combinations(range(4), 2)))
    theta, weights = .47, (.37, .81)
    rows = []
    for parity in (0, 1):
        _, iso = enc.isometry(parity)
        u = np.eye(1 << enc.m, dtype=complex)
        tapes = []
        for i, j in ((0, 2), (1, 3)):
            for circuit in (enc.phase(j, -theta, parity),
                            enc.su2(i, j, np.array([[1., 1.], [-1., 1.]])/np.sqrt(2), parity)):
                scalar, tape = circuit
                u = enc.execute(circuit)@u
                tapes.append(dict(scalar=encode(scalar), rotations=[dict(word=w.row(), theta=float(t)) for w, t in tape]))
        effect = np.zeros_like(u)
        instruments, programs = [], []
        for i, weight in enumerate(weights):
            program = weighted_program(enc.b(i, parity), weight, enc.m)
            programs.append(program)
            maps = execute_weighted(program)
            instruments.append((maps[3], maps[:3]))
        for click, no_click in reversed(instruments):
            effect = click.conj().T@click+sum(k.conj().T@effect@k for k in no_click)
        effect = u.conj().T@effect@u
        rows.append(dict(parity=parity, theta=theta, weights=list(weights), circuits=tapes,
                         instruments=programs, code_effect=encode(iso.conj().T@effect@iso)))
    return rows


def interaction():
    enc = Encoding(4, list(itertools.combinations(range(4), 2)))
    basis, iso = enc.isometry(0)
    psi = np.zeros(16, complex)
    psi[[5, 6, 9, 10]] = .5
    circuit = enc.density(1, 3, np.pi, 0)
    after = enc.execute(circuit)@iso@psi[basis]
    scalar, tape = circuit
    return dict(parity=0, angle=float(np.pi), modes=[1, 3], initial=encode(psi),
                encoded_output=encode(after),
                circuit=dict(scalar=encode(scalar), rotations=[dict(word=w.row(), theta=float(t)) for w, t in tape]))


def noise():
    rows = []
    for name, n, edges in (('triangle', 3, [(0, 1), (0, 2), (1, 2)]),
                           ('two_cycles', 4, [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]),
                           ('complete_four', 4, list(itertools.combinations(range(4), 2)))):
        enc = Encoding(n, edges)
        vacuum, _ = enc.prepare()
        exposure = .2
        p = -np.expm1(-exposure)/2
        fidelity = 0.
        syndromes = []
        for errors in itertools.product((0, 1), repeat=len(edges)):
            z = sum(bit << i for i, bit in enumerate(errors))
            error = Word(z=z)
            state = error.matrix(len(edges))@vacuum
            probability = float(p**sum(errors)*(1-p)**(len(edges)-sum(errors)))
            fidelity += probability*abs(np.vdot(vacuum, state))**2
        for chord in enc.chords:
            e = Word(z=1 << chord)
            syndromes.append([int(not e.commutes(s)) for s in enc.loops])
        rows.append(dict(name=name, exposure=exposure, phase_flip_probability=float(p),
                         cycle_rank=len(enc.loops), chord_syndromes=syndromes,
                         vacuum_fidelity=float(fidelity),
                         failure_lower=float(1-(1-p)**len(enc.loops))))
    return rows


def mode_rotation(enc, target):
    """Tree Givens elimination, including every zero subtree, sends f to e0."""
    vector = np.asarray(target, complex).copy()
    if vector.shape != (enc.n,) or not np.all(np.isfinite(vector)) or abs(np.linalg.norm(vector)-1) > 1e-12:
        raise ValueError('normalized finite mode kernel required')
    tape = []
    for child in reversed(list(enc.parent)):
        if child == 0:
            continue
        parent = enc.parent[child]
        x, y = vector[parent], vector[child]
        norm = float(np.hypot(abs(x), abs(y)))
        g = (np.array([[x.conjugate(), y.conjugate()], [-y, x]])/norm
             if norm else np.eye(2, dtype=complex))
        vector[[parent, child]] = g@vector[[parent, child]]
        tape.append((parent, child, g))
    return tape


def resolved_reads():
    enc = Encoding(5, [(0, 1), (0, 2), (1, 2), (2, 3), (3, 4)])
    rows = []
    for name, target in [('complex', np.array([1., .3j, 0., -.4+.2j, .1])),
                         ('zero_subtree', np.array([1., .4j, 0., 0., 0.]))]:
        target = target/np.linalg.norm(target)
        instructions = mode_rotation(enc, target)
        blocks = []
        for parity in (0, 1):
            _, iso = enc.isometry(parity)
            u = np.eye(1 << enc.m, dtype=complex)
            tapes = []
            for i, j, g in instructions:
                circuit = enc.su2(i, j, g, parity)
                u = enc.execute(circuit)@u
                scalar, tape = circuit
                tapes.append(dict(modes=[i, j], matrix=encode(g), scalar=encode(scalar),
                                  rotations=[dict(word=w.row(), theta=float(t)) for w, t in tape]))
            effect = u.conj().T@(np.eye(len(u))-enc.b(0, parity).matrix(enc.m))@u/2
            blocks.append(dict(parity=parity, circuits=tapes,
                               code_effect=encode(iso.conj().T@effect@iso)))
        rows.append(dict(name=name, kernel=encode(target), blocks=blocks,
                         forward_rotations=6*(enc.n-1), qnd_rotations_with_return=12*(enc.n-1)))
    return rows

"""Deterministic Slater preparation and a returned source QND read."""

import itertools
import numpy as np
from m1_source_realization.model import decompose
from m1_fermionic_source.model import Encoding
from m1_fermionic_source.experiments import mode_rotation
from m1_fermionic_source.circuits import measurement_tape, execute
from .format import pack


def circuit_row(enc, modes, matrix, circuit):
    scalar, tape = circuit
    return dict(modes=list(modes), matrix=pack(matrix), scalar=pack(scalar),
                rotations=[dict(word=w.row(), theta=float(t)) for w, t in tape])


def compile_unitary(enc, u, parity):
    diagonal, moves = decompose(u)
    rows = []
    for i, value in enumerate(diagonal):
        g = np.array([[value]])
        rows.append(circuit_row(enc, [i], g, enc.phase(i, float(np.angle(value)), parity)))
    swap = np.array([[0., 1.], [1., 0.]], complex)
    for i, j, g in moves:
        path = enc.path(i, j)
        swaps = list(zip(path[:-2], path[1:-1]))
        for left, right in swaps:
            rows.append(circuit_row(enc, [left, right], swap, enc.fswap(left, right, parity)))
        left, right = path[-2:]
        rows.append(circuit_row(enc, [left, right], g, enc.su2(left, right, g, parity)))
        for left, right in reversed(swaps):
            rows.append(circuit_row(enc, [left, right], swap, enc.fswap(left, right, parity)))
    return rows


def execute_rows(enc, rows):
    # Producer contracts the emitted Pauli words with its own Word implementation.
    from m1_fermionic_source.pauli import Word
    result = np.eye(1 << enc.m, dtype=complex)
    for row in rows:
        scalar = complex(*row['scalar'])
        tape = [(Word(int(r['word'][0]), int(r['word'][1]), r['word'][2]), r['theta']) for r in row['rotations']]
        result = enc.execute((scalar, tape))@result
    return result


def case(rank, gauge=None):
    enc = Encoding(4, [(0, 1), (0, 2), (0, 3), (1, 2)])
    z = np.exp(2j*np.pi*np.outer(np.arange(4), np.arange(4))/4)/2
    z = np.exp(.17j*np.arange(4)**2)[:, None]*z
    target = z[:, :rank]@z[:, :rank].conj().T
    h = (z*np.array([-2., -1., .9, 1.4]))@z.conj().T
    if gauge is not None:
        z = z@gauge
    parity = rank % 2
    basis, iso = enc.isometry(parity)
    gates = compile_unitary(enc, z, parity)
    prepared = execute_rows(enc, gates)@iso[:, basis.index((1 << rank)-1)]
    f = np.array([1., .3j, -.4+.2j, .1], complex)
    f /= np.linalg.norm(f)
    reflected = np.eye(4)-2*np.outer(f, f.conj())
    mean_read_energy = float(np.trace(h@(target+reflected@target@reflected)/2).real+3)
    read = [circuit_row(enc, [i, j], g, enc.su2(i, j, g, parity)) for i, j, g in mode_rotation(enc, f)]
    gathering = execute_rows(enc, read)
    qnd = measurement_tape(enc.b(0, parity), enc.m)
    instrument = np.kron(np.eye(2), gathering.conj().T)@execute(qnd, enc.m+1)@np.kron(np.eye(2), gathering)
    maps = [instrument[b*(1 << enc.m):(b+1)*(1 << enc.m), :1 << enc.m] for b in (0, 1)]
    branches = []
    for outcomes in itertools.product((0, 1), repeat=len(enc.loops)):
        state, probability = enc.prepare(list(outcomes))
        # Pair toggles implement occupation preparation, not an asserted encoder.
        for v in range(1, rank):
            for i, j in zip(enc.path(0, v), enc.path(0, v)[1:]):
                state = enc.a(i, j).matrix(enc.m)@state
        state = execute_rows(enc, gates)@state
        branches.append(dict(outcomes=list(outcomes), probability=probability, density=pack(np.outer(state, state.conj()))))
    return dict(rank=rank, vertices=4, edges=[list(e) for e in enc.edges], parity=parity,
                orbitals=pack(z), covariance=pack(target), basis=basis, isometry=pack(iso),
                gates=gates, branches=branches, kernel=pack(f), read=read, qnd=qnd,
                kraus=[pack(k) for k in maps], probabilities=[float(np.vdot(k@prepared, k@prepared).real) for k in maps],
                mean_read_energy=mean_read_energy,
                pulse_rotations=sum(len(g['rotations']) for g in gates)+2*sum(len(g['rotations']) for g in read))


def candidate():
    return [case(2), case(3)]

"""Native proper-code Clifford/parity circuits; no free encoded rotations."""

import numpy as np
from scipy.linalg import expm

from .model import encode
from .pauli import Word


H = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2)


def pauli_tape(word, theta, qubits):
    if word*word != Word() or (word.x | word.z) >= 1 << qubits:
        raise ValueError('Hermitian supported Pauli required')
    sign = 1j**((word.phase-(word.x & word.z).bit_count()) % 4)
    if abs(sign.imag) > 0:
        raise ValueError('non-Hermitian Pauli')
    support = [j for j in range(qubits) if (word.x | word.z) >> j & 1]
    if not support:
        raise ValueError('nonempty Pauli support required')
    basis, inverse = [], []
    for j in support:
        if word.x >> j & 1:
            if word.z >> j & 1:
                basis.append(['sdg', j])
                inverse.insert(0, ['s', j])
            basis.append(['h', j])
            inverse.insert(0, ['h', j])
    def cx(j):
        return [['h', qubits], ['cz', j, qubits], ['h', qubits]]
    compute = [gate for j in support for gate in cx(j)]
    uncompute = [gate for j in reversed(support) for gate in cx(j)]
    return basis+compute+[['rz', qubits, float(theta*sign.real)]]+uncompute+inverse


def native_one(kind, angle=None):
    u = dict(h=H, s=np.diag([1., 1j]), sdg=np.diag([1., -1j]), x=np.array([[0., 1.], [1., 0.]])).get(kind)
    if kind == 'rz':
        u = np.diag(np.exp(1j*angle*np.array([1., -1.])))
    elif kind == 'ry':
        u = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    if u is None:
        raise ValueError('unknown native one-code pulse')
    processor = np.eye(6, dtype=complex)
    processor[:2, :2] = u
    return processor


def native_cz():
    # The audited source sum D_p is diag(i I3,0). Its pi pulse, up to one
    # scalar on the pair code, is precisely CZ. Use the actual M6 pulse.
    return -expm(-np.pi*np.diag([1j, 1j, 1j, 0., 0., 0.]))[:4, :4]


def execute(tape, qubits):
    out = np.eye(1 << qubits, dtype=complex)
    for gate in tape:
        kind, j, *tail = gate
        if kind == 'cz':
            k = tail[0]
            diagonal = np.diag(native_cz())
            out *= np.array([diagonal[2*((r >> j) & 1)+((r >> k) & 1)]
                             for r in range(len(out))])[:, None]
        else:
            u = native_one(kind, tail[0] if tail else None)[:2, :2]
            for r in range(len(out)):
                if not r >> j & 1:
                    pair = [r, r | (1 << j)]
                    out[pair, :] = u@out[pair, :]
    return out


def gadgets():
    cases = [(2, Word(1, 2), .37), (3, Word(3, 6, 1), -.51),
             (3, Word(5, 7, 0), .29), (1, Word(0, 1, 2), .63)]
    rows = []
    for n, p, theta in cases:
        tape = pauli_tape(p, theta, n)
        rows.append(dict(qubits=n, word=p.row(), theta=theta, tape=tape,
                         unitary=encode(execute(tape, n+1)), pulses=len(tape),
                         pulse_code_events=2*len(tape)))
    return rows


def measurement_tape(word, qubits):
    full = pauli_tape(word, .2, qubits)
    middle = next(k for k, g in enumerate(full) if g[0] == 'rz')
    # After parity accumulation, read the ancilla and undo only the data
    # basis changes. Its classical bit is retained rather than uncomputed.
    support = [j for j in range(qubits) if (word.x | word.z) >> j & 1]
    basis_count = sum(1+((word.z >> j) & 1) for j in support if word.x >> j & 1)
    tape = full[:middle]+full[-basis_count:] if basis_count else full[:middle]
    if (word.phase-(word.x & word.z).bit_count()) % 4 == 2:
        tape += [['x', qubits]]
    return tape


def measurements():
    rows = []
    for n, p in ((3, Word(3, 6, 1)), (3, Word(0, 7)), (1, Word(0, 1, 2))):
        tape = measurement_tape(p, n)
        rows.append(dict(qubits=n, word=p.row(), tape=tape,
                         unitary=encode(execute(tape, n+1)), pulses=len(tape),
                         pulse_code_events=2*len(tape)+2))
    return rows


def weighted_program(word, weight, qubits):
    """QND parity read, actual helper reset, conditional native pulse, read."""
    if type(weight) not in (int, float) or not np.isfinite(weight) or not 0 <= weight <= 1:
        raise ValueError('finite acceptance probability in [0,1] required')
    tape = measurement_tape(word, qubits)
    return dict(qubits=qubits, word=word.row(), weight=float(weight), qnd_tape=tape,
                reset_before_acceptance=True,
                acceptance_angles=[0., float(np.arctan2(np.sqrt(weight), np.sqrt(1-weight)))],
                outcomes=[[n, b] for n in (0, 1) for b in (0, 1)],
                pulses_per_path=len(tape)+1, code_events_per_path=2*(len(tape)+1)+8)


def execute_weighted(program):
    n = program['qubits']
    size = 1 << n
    qnd = execute(program['qnd_tape'], n+1)
    maps = []
    for occupied, accepted in program['outcomes']:
        data = qnd[occupied*size:(occupied+1)*size, :size]
        helper = np.eye(2)[:, occupied]
        if program['reset_before_acceptance']:
            # The measured helper is |occupied>. This is its explicit
            # proper-code reset branch |0><occupied|, not a free new ancilla.
            helper = np.outer(np.eye(2)[:, 0], np.eye(2)[occupied])@helper
        pulse = native_one('ry', program['acceptance_angles'][occupied])[:2, :2]
        maps.append((pulse@helper)[accepted]*data)
    return maps


def weighted_measurements():
    rows = []
    for p, weight in ((Word(3, 6, 1), .37), (Word(0, 7, 2), .81),
                      (Word(0, 7), 0.), (Word(0, 7), 1.)):
        program = weighted_program(p, weight, 3)
        rows.append(dict(program=program, kraus=[encode(k) for k in execute_weighted(program)]))
    return rows

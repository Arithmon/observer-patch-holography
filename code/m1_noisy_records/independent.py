"""Projector-instrument oracle and independently embedded elementary matrices.

No producer imports. All externally supplied tapes are validated before use.
Channel tests use the full Choi matrix, not a favorable input state or marginal.
"""

import itertools
import math
import numpy as np


def need(ok, message):
    if not ok:
        raise ValueError(message)


def keys(obj, expected):
    need(type(obj) is dict and set(obj) == set(expected.split()), 'exact object keys')


def integer(value, lo, hi):
    need(type(value) is int and lo <= value <= hi, 'integer range')


def close(a, b, message):
    a, b = np.asarray(a), np.asarray(b)
    need(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
         and np.max(np.abs(a-b), initial=0) < 3e-12, message)


def ry(t):
    return np.array([[np.cos(t/2), -np.sin(t/2)], [np.sin(t/2), np.cos(t/2)]])


def single(gate, target, n):
    result = np.array([[1.]])
    for bit in reversed(range(n)):
        result = np.kron(result, gate if bit == target else np.eye(2))
    return result


def elementary(op, n):
    need(type(op) is list and len(op) in (2, 3) and type(op[0]) is str, 'gate syntax')
    name = op[0]
    need(name in ('h', 't', 'tdg', 'ry', 'cx', 'cz'), 'supported gate')
    need(len(op) == (3 if name in ('ry', 'cx', 'cz') else 2), 'gate arity')
    integer(op[1], 0, n-1)
    if name in ('cx', 'cz'):
        integer(op[2], 0, n-1)
        need(op[1] != op[2], 'distinct operands')
        p1 = single(np.diag([0, 1]), op[1], n)
        action = np.array([[0, 1], [1, 0]]) if name == 'cx' else np.diag([1, -1])
        return np.eye(2**n)-p1+p1 @ single(action, op[2], n)
    if name == 'ry':
        need(type(op[2]) in (int, float) and math.isfinite(op[2]), 'finite angle')
        gate = ry(op[2])
    elif name == 'h':
        gate = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    else:
        gate = np.diag([1, np.exp((1 if name == 't' else -1)*1j*np.pi/4)])
    return single(gate, op[1], n)


def isometry(tape, n=7, inputs=2):
    need(type(tape) is list and 0 < len(tape) <= 128, 'bounded nonempty tape')
    state = np.eye(2**n, dtype=complex)[:, :2**inputs]
    for op in tape:
        state = elementary(op, n) @ state
    return state


def choi_from_isometry(v):
    # Output includes data, both history bits AND the abort bit; only e0,e1
    # are traced. An omitted environmental copy creates detectable coherences.
    need(v.shape == (128, 4), 'isometry dimensions')
    out = np.zeros((128, 128), complex)
    for environment in range(4):
        k = v[32*environment:32*(environment+1)]
        w = k.reshape(-1)
        out += np.outer(w, w.conj())
    return out


def reference(angles):
    theta, phi, alpha = angles
    h = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    x, z = np.array([[0, 1], [1, 0]]), np.diag([1, -1])
    cx = sum((np.outer(np.eye(4)[:, j ^ (2 if j & 1 else 0)], np.eye(4)[:, j])
              for j in range(4)), np.zeros((4, 4)))
    start = cx @ np.kron(np.eye(2), ry(theta))
    out = np.zeros((128, 128), complex)
    kraus = []
    for r0, r1 in itertools.product(range(2), repeat=2):
        p0 = np.kron(np.eye(2), np.diag([1-r0, r0]))
        p1 = np.kron(np.diag([1-r1, r1]), np.eye(2))
        middle = np.kron(h if r0 else np.eye(2), ry(phi))
        feedback = np.linalg.matrix_power(z, r1) @ np.linalg.matrix_power(x, r0)
        abort = r0*r1
        final = np.kron(np.eye(2), (ry(alpha) if not abort else np.eye(2)) @ feedback)
        k = final @ p1 @ middle @ p0 @ start
        full = np.zeros((32, 4), complex)
        offset = 4*r0+8*r1+16*abort
        full[offset:offset+4] = k
        v = full.reshape(-1)
        out += np.outer(v, v.conj())
        kraus.append(k)
    close(sum(k.conj().T @ k for k in kraus), np.eye(4), 'reference TP')
    need(all(np.linalg.norm(k) > .01 for k in kraus), 'every history, including abort, present')
    return out


def verify_case(row, expected_angles):
    keys(row, 'angles tape')
    need(type(row['angles']) is list and len(row['angles']) == 3
         and all(type(x) in (int, float) and math.isfinite(x) for x in row['angles']), 'angles')
    close(row['angles'], expected_angles, 'frozen experiment')
    v = isometry(row['tape'])
    close(v.conj().T @ v, np.eye(4), 'compiled TP')
    close(choi_from_isometry(v), reference(expected_angles), 'complete history instrument')


def verify_toffoli(tape):
    actual = isometry(tape, n=3, inputs=3)
    ideal = np.zeros((8, 8))
    for j in range(8):
        ideal[j ^ (4 if j & 3 == 3 else 0), j] = 1
    close(actual, ideal, 'Toffoli, including coherent phases')


def verify_budget(row):
    fields = ('record_bits_q_degree record_bits_log_degree depth_q_degree depth_log_degree '
              'volume_q_degree volume_log_degree diagnostic_bits_q_degree diagnostic_bits_log_degree '
              'physical_volume_q_degree physical_volume_log_degree levels fault_power fault_error_q_degree '
              'synthesis_per_gate_q_degree synthesis_error_q_degree accounting_q_degree '
              'total_error_q_degree accounting_error_q_degree')
    keys(row, fields)
    need(all(type(x) is int for x in row.values()), 'integer exponent ledger')
    # Freeze the inputs proved in the note, then independently derive all totals.
    for name, value in [('record_bits_q_degree', 4), ('record_bits_log_degree', 8),
                        ('depth_q_degree', 1), ('depth_log_degree', 8), ('levels', 5),
                        ('synthesis_per_gate_q_degree', -24), ('accounting_q_degree', 4)]:
        need(row[name] == value, 'proved budget input '+name)
    volume = row['record_bits_q_degree'] + row['depth_q_degree']
    power = 1
    for _ in range(row['levels']):
        power *= 2
    synthesis = row['record_bits_q_degree'] + row['synthesis_per_gate_q_degree']
    noise = volume-power
    for name, expected in [('volume_q_degree', volume), ('volume_log_degree', 16),
                           ('diagnostic_bits_q_degree', volume), ('diagnostic_bits_log_degree', 17),
                           ('physical_volume_q_degree', volume+row['depth_q_degree']),
                           ('physical_volume_log_degree', 17+row['depth_log_degree']),
                           ('fault_power', power), ('fault_error_q_degree', noise),
                           ('synthesis_error_q_degree', synthesis),
                           ('total_error_q_degree', max(noise, synthesis)),
                           ('accounting_error_q_degree', max(noise, synthesis)+4)]:
        need(row[name] == expected, 'derived budget '+name)


def verify_terminal(row):
    keys(row, 'repetition horizon majority_error_coefficients raw_fanout')
    integer(row['repetition'], 5, 5)
    integer(row['raw_fanout'], 7, 7)
    need(type(row['horizon']) is list and row['horizon'] == [1, 2]
         and all(type(x) is int for x in row['horizon']), 'terminal sample horizon')
    coeff = [0]*6
    # Enumerate every five-bit error pattern and expand p^w (1-p)^(5-w).
    for mask in range(32):
        w = mask.bit_count()
        if w >= 3:
            for j in range(6-w):
                coeff[w+j] += (-1)**j*math.comb(5-w, j)
    need(type(row['majority_error_coefficients']) is list
         and all(type(x) is int for x in row['majority_error_coefficients'])
         and row['majority_error_coefficients'] == coeff, 'majority polynomial')
    # A corrupt bare broadcast bit changes every target, not one position.
    control, targets = 1, [0]*row['raw_fanout']
    for i in range(len(targets)):
        targets[i] ^= control
    need(sum(targets) == 7, 'raw feedback counterexample')


def verify_export(row):
    keys(row, 'copies tape')
    integer(row['copies'], 5, 5)
    v = isometry(row['tape'])
    # Input is data+classical record, possibly reference-entangled on data.
    # Check every allowed matrix unit and every pattern of bad output islands.
    # Physical reads dephase probe bits; these columns are already basis-coded.
    for label in (0, 1):
        for mask in range(32):
            indices = np.arange(128) ^ (mask << 2)
            noisy = v[indices]
            expected = np.zeros_like(noisy)
            for data in (0, 1):
                out = data+2*label+4*((31*label) ^ mask)
                expected[out, data+2*label] = 1
            for i, j in itertools.product(range(2), repeat=2):
                a, b = i+2*label, j+2*label
                close(np.outer(noisy[:, a], noisy[:, b].conj()),
                      np.outer(expected[:, a], expected[:, b].conj()),
                      'public export with spectators and island faults')


def verify_evidence(evidence):
    from .archive_check import verify as verify_archive
    from .computation_code_check import verify as verify_code
    from .instrument_check import verify as verify_instruments
    keys(evidence, 'cases native_instruments computation_code toffoli archive export budget terminal')
    verify_instruments(evidence['native_instruments'])
    verify_code(evidence['computation_code'])
    cases = ((.37, -.61, .43), (math.pi/2, math.pi/3, -math.pi/4), (-.83, .29, 1.17))
    need(type(evidence['cases']) is list and len(evidence['cases']) == len(cases), 'complete cases')
    for row, angles in zip(evidence['cases'], cases):
        verify_case(row, angles)
    verify_toffoli(evidence['toffoli'])
    verify_archive(evidence['archive'])
    verify_export(evidence['export'])
    verify_budget(evidence['budget'])
    verify_terminal(evidence['terminal'])

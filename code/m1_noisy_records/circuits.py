"""Produce elementary tapes; the verifier reconstructs their channels independently."""

import math
import numpy as np


CASES = ((0.37, -0.61, 0.43), (math.pi/2, math.pi/3, -math.pi/4),
         (-0.83, 0.29, 1.17))


def toffoli(a, b, target):
    """Exact Clifford+T decomposition, including its relative phases."""
    return [['h', target], ['cx', b, target], ['tdg', target],
            ['cx', a, target], ['t', target], ['cx', b, target],
            ['tdg', target], ['cx', a, target], ['t', b], ['t', target],
            ['h', target], ['cx', a, b], ['t', a], ['tdg', b], ['cx', a, b]]


def tape(theta, phi, alpha):
    # d0,d1,r0,r1,abort,e0,e1. All ancillas start in zero.
    # The environmental copies are discarded, never coherently erased.
    ops = [['ry', 0, theta], ['cx', 0, 1], ['cx', 0, 2], ['cx', 2, 5],
           ['ry', 1, -math.pi/4], ['cz', 2, 1], ['ry', 1, math.pi/4],
           ['ry', 0, phi], ['cx', 1, 3], ['cx', 3, 6]]
    ops += toffoli(2, 3, 4)
    ops += [['cx', 2, 0], ['cz', 3, 0],
            # controlled Ry(-alpha), followed by unconditional Ry(alpha):
            # accepted branches rotate, the abort branch retains its output.
            ['ry', 0, -alpha/2], ['cx', 4, 0], ['ry', 0, alpha/2],
            ['cx', 4, 0], ['ry', 0, alpha]]
    return ops


def execute(ops, qubits=7, inputs=2):
    """Forward state-vector/isometry execution (little-endian qubits)."""
    if (type(qubits) is not int or type(inputs) is not int or
            not 0 <= inputs <= qubits or not 1 <= qubits or qubits+inputs > 20):
        raise ValueError('finite executor dimensions')
    if type(ops) is not list or not 0 < len(ops) <= 1024:
        raise ValueError('bounded nonempty tape required')
    for op in ops:
        if type(op) is not list or len(op) < 2 or op[0] not in ('h', 't', 'tdg', 'ry', 'cx', 'cz'):
            raise ValueError('unknown or malformed gate')
        if len(op) != (3 if op[0] in ('cx', 'cz', 'ry') else 2):
            raise ValueError('gate arity')
        operands = op[1:] if op[0] in ('cx', 'cz') else op[1:2]
        if (any(type(q) is not int or not 0 <= q < qubits for q in operands)
                or len(set(operands)) != len(operands)):
            raise ValueError('distinct valid operands required')
        if op[0] == 'ry' and (type(op[2]) not in (int, float) or not math.isfinite(op[2])):
            raise ValueError('finite angle required')
    state = np.zeros((2**qubits, 2**inputs), complex)
    state[:2**inputs] = np.eye(2**inputs)
    for op in ops:
        name, target = op[0], op[-1] if op[0] in ('cx', 'cz') else op[1]
        if name in ('cx', 'cz'):
            control = op[1]
            for row in range(2**qubits):
                if row >> control & 1 and not row >> target & 1:
                    other = row | 1 << target
                    if name == 'cx':
                        state[[row, other]] = state[[other, row]]
                    else:
                        state[other] *= -1
        else:
            if name == 'ry':
                c, s = math.cos(op[2]/2), math.sin(op[2]/2)
                gate = np.array([[c, -s], [s, c]])
            elif name == 'h':
                gate = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
            else:
                gate = np.diag([1, np.exp((1 if name == 't' else -1)*1j*math.pi/4)])
            for row in range(2**qubits):
                if not row >> target & 1:
                    pair = [row, row | 1 << target]
                    state[pair] = gate @ state[pair]
    return state


def candidate():
    from .archive import evidence
    from .computation_code import candidate as code_candidate
    return {'cases': [{'angles': list(angles), 'tape': tape(*angles)} for angles in CASES],
            'computation_code': code_candidate(),
            'toffoli': toffoli(0, 1, 2),
            'archive': evidence(),
            'export': {'copies': 5, 'tape': [['cx', 1, j] for j in range(2, 7)]},
            'budget': {'record_bits_q_degree': 4, 'record_bits_log_degree': 8,
                       'depth_q_degree': 1, 'depth_log_degree': 8,
                       'volume_q_degree': 5, 'volume_log_degree': 16,
                       'diagnostic_bits_q_degree': 5, 'diagnostic_bits_log_degree': 17,
                       'physical_volume_q_degree': 6, 'physical_volume_log_degree': 25,
                       'levels': 5, 'fault_power': 32,
                       'fault_error_q_degree': -27, 'synthesis_per_gate_q_degree': -24,
                       'synthesis_error_q_degree': -20, 'accounting_q_degree': 4,
                       'total_error_q_degree': -20, 'accounting_error_q_degree': -16},
            'terminal': {'repetition': 5, 'horizon': [1, 2],
                         'majority_error_coefficients': [0, 0, 0, 10, -15, 6],
                         'raw_fanout': 7}}

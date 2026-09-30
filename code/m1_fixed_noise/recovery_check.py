"""Independent reverse (Heisenberg) propagation for the recovery certificate.

This module never imports the producer. It builds every stabilizer support
from the nonzero three-bit columns and propagates output sensitivities BACK
through the circuit; the producer propagates faults forward as polynomials.
"""

import hashlib
import itertools
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(obj):
    return hashlib.sha256(json.dumps(obj, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def canonical(config):
    gates, flags = [], []
    next_bit = 14+6*config['rounds']
    for repetition in range(config['rounds']):
        for basis, bit in itertools.product(('X', 'Z'), range(3)):
            block_flags = []
            for bank in ((7, 8, 9, 10), (11, 12, 13, 14)):
                gates += [['r', [q], -1] for q in bank]
                gates += [['h', [bank[0]], -1]]
                gates += [['cx', [bank[0], q], -1] for q in bank[1:]]
                recorded = []
                for i, j in config['verification']:
                    gates += [['r', [15], -1], ['cx', [bank[i], 15], -1],
                              ['cx', [bank[j], 15], -1], ['mz', [15], next_bit]]
                    recorded.append(next_bit)
                    next_bit += 1
                block_flags.append(recorded)
            flags.append(block_flags)
            support = [column-1 for column in range(1, 8) if column & (1 << bit)]
            gates += [['cx' if basis == 'X' else 'cz', [7+j, q], -1] for j, q in enumerate(support)]
            read = 14+repetition*6+bit+(3 if basis == 'Z' else 0)
            for q in range(7, 11):
                gates += [['h', [q], -1], ['mz', [q], read]]
    for op in ('correct_x', 'correct_z'):
        gates += [[op, [q], -1] for q in range(7)]
    expanded = []
    for gate in gates:
        expanded.append(gate)
        expanded += [['i', [q], -1] for q in sorted(set(range(16))-set(gate[1]))]
    return expanded, flags, next_bit


def reverse_columns(tape):
    xx = [1 << i if i < 7 else 0 for i in range(16)]
    zz = [1 << (i+7) if i < 7 else 0 for i in range(16)]
    reversed_blocks = []
    for op, operands, record in reversed(tape):
        if op == 'mz':
            xx[operands[0]], zz[operands[0]] = 1 << record, 0
        reversed_blocks.append([v for q in operands for v in (xx[q], zz[q])])
        if op == 'h':
            a = operands[0]
            xx[a], zz[a] = zz[a], xx[a]
        elif op == 'cx':
            a, b = operands
            xx[a] ^= xx[b]
            zz[b] ^= zz[a]
        elif op == 'cz':
            a, b = operands
            xx[a] ^= zz[b]
            xx[b] ^= zz[a]
        elif op == 'r':
            xx[operands[0]] = zz[operands[0]] = 0
        else:
            need(op in ('mz', 'i', 'correct_x', 'correct_z'), 'known circuit location')
    return [x for block in reversed(reversed_blocks) for x in block]


def stabilizer_set():
    words = [(0, 0)]
    for axis in (0, 1):
        for digit in range(3):
            support = sum(1 << (i-1) for i in range(1, 8) if i >> digit & 1)
            new = []
            for a, b in words:
                new.append((a ^ support, b) if axis == 0 else (a, b ^ support))
            words += new
    return set(words)


def residual(signature, flags, rounds):
    invalid_first = False
    for first, second in flags:
        a = sum((signature >> i & 1) for i in first)
        b = sum((signature >> i & 1) for i in second)
        need(not (a and b), 'one fault cannot reject both independent candidates')
        invalid_first |= bool(a)
    if invalid_first:
        return (0, 0)
    history = [(signature >> (14+6*r)) % 64 for r in range(rounds)]
    value = next((history[i] for i in range(1, rounds) if history[i] == history[i-1]), None)
    if rounds == 1:
        value = history[0]
    if value is None:
        return None
    x, z = signature % 128, (signature >> 7) % 128
    if value // 8:
        x ^= 2**(value//8-1)
    if value % 8:
        z ^= 2**(value % 8-1)
    return x, z


def replay(config):
    tape, flags, _ = canonical(config)
    columns = reverse_columns(tape)
    stabilizers = stabilizer_set()
    allowed = set(stabilizers)
    for a, b in stabilizers:
        for q in range(7):
            for dx, dz in ((1, 0), (0, 1), (1, 1)):
                allowed.add((a ^ (dx << q), b ^ (dz << q)))
    offset, count, histogram, bad = 0, 0, {}, []
    for index, (_, wires, _) in enumerate(tape):
        for errors in itertools.product((0, 1, 2, 3), repeat=len(wires)):
            if all(e == 0 for e in errors):
                continue
            effect = 0
            for j, error in enumerate(errors):
                for bit in range(2):
                    if error >> bit & 1:
                        effect ^= columns[offset+2*j+bit]
            output = residual(effect, flags, config['rounds'])
            key = 'abort' if output is None else ':'.join(map(str, output))
            histogram[key] = histogram.get(key, 0)+1
            count += 1
            if output not in allowed:
                bad.append([index, list(errors), key])
        offset += 2*len(wires)
    return dict(config=config, locations=len(tape), one_location_faults=count,
                primitive_slots=sum(g[0] != 'i' for g in tape),
                idle_locations=sum(g[0] == 'i' for g in tape),
                tape_sha256=digest(tape), response_sha256=digest(columns),
                residual_histogram=histogram, failures=bad)


def verify_recovery(row):
    from .instrument_check import check_program
    need(type(row) is dict, 'recovery object')
    required = dict(rounds=4, candidates=2, verification=[[0, 1], [1, 2], [2, 3]])
    need(row.get('config') == required, 'complete bounded recovery schedule')
    expected = replay(required)
    check_program(canonical(required)[0])
    need(not expected['failures'], 'independent single-fault recovery proof')
    need(row == expected, 'complete recovery evidence replay')
    # Reject booleans-as-integers and noncanonical numeric spellings too.
    need(json.dumps(row, sort_keys=True) == json.dumps(expected, sort_keys=True), 'strict recovery types')

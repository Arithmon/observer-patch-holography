"""Bounded Shor recovery; exact binary propagation of every one-location fault.

The two ancilla candidates are BOTH prepared and verified. The first accepted
one is used. Each stabilizer is read in four rounds, with the first consecutive
agreement selected. All slots, unused candidates and later rounds are paid.
"""

import hashlib
import itertools
import json

import numpy as np

CHECKS = (85, 102, 120)
CONFIG = dict(rounds=4, candidates=2, verification=[[0, 1], [1, 2], [2, 3]])


def digest(obj):
    return hashlib.sha256(json.dumps(obj, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def parity(x):
    return x.bit_count() & 1


def syndrome(x, z):
    return sum(parity(z & h) << k for k, h in enumerate(CHECKS)) | sum(
        parity(x & h) << (k+3) for k, h in enumerate(CHECKS))


def correction(s):
    phase, bit = s & 7, s >> 3
    return (1 << (bit-1) if bit else 0, 1 << (phase-1) if phase else 0)


def stabilizers():
    span = [0]
    for h in CHECKS:
        span += [x ^ h for x in span]
    return {(x, z) for x in span for z in span}


def correctable():
    errors = {(0, 0)} | {(x << i, z << i) for i in range(7) for x, z in ((1, 0), (0, 1), (1, 1))}
    return {(x ^ a, z ^ b) for x, z in errors for a, b in stabilizers()}


def program(config=CONFIG):
    """Entries: [opcode, operands, recorded parity bit or -1]. Idle is explicit."""
    if (type(config) is not dict or set(config) != {'rounds', 'candidates', 'verification'}
            or type(config['rounds']) is not int or not 1 <= config['rounds'] <= 4
            or type(config['candidates']) is not int or config['candidates'] != 2
            or type(config['verification']) is not list
            or any(type(pair) is not list or len(pair) != 2
                   or any(type(q) is not int or not 0 <= q < 4 for q in pair)
                   or pair[0] == pair[1] for pair in config['verification'])):
        raise ValueError('bounded recovery configuration required')
    tape = []
    output = 14 + 6*config['rounds']
    flag_groups = []
    def step(op, operands, record=-1):
        tape.append([op, operands, record])
        tape.extend(['i', [q], -1] for q in range(16) if q not in operands)
    for r in range(config['rounds']):
        for g in range(6):
            local_flags = []
            for candidate in range(2):
                cat = list(range(7+4*candidate, 11+4*candidate))
                for q in cat:
                    step('r', [q])
                step('h', [cat[0]])
                for q in cat[1:]:
                    step('cx', [cat[0], q])
                bits = []
                for a, b in config['verification']:
                    step('r', [15])
                    step('cx', [cat[a], 15])
                    step('cx', [cat[b], 15])
                    step('mz', [15], output)
                    bits.append(output)
                    output += 1
                local_flags.append(bits)
            flag_groups.append(local_flags)
            targets = [i for i in range(7) if CHECKS[g % 3] >> i & 1]
            # Linear reference route uses candidate A. A flagged A fault is
            # handled by the exact one-fault branch reduction in interpret().
            for ancilla, data in zip(range(7, 11), targets):
                step('cx' if g < 3 else 'cz', [ancilla, data])
            for q in range(7, 11):
                step('h', [q])
                step('mz', [q], 14+6*r+g)
    for op in ('correct_x', 'correct_z'):
        for q in range(7):
            step(op, [q])
    return tape, flag_groups, output


def responses(tape, outputs):
    """Forward symbolic Pauli propagation; integers are GF(2) polynomials."""
    x, z = [0]*16, [0]*16
    records = [0]*outputs
    variables, spans = 0, []
    for op, qs, rec in tape:
        spans.append((variables, len(qs)))
        if op == 'h':
            q = qs[0]
            x[q], z[q] = z[q], x[q]
        elif op == 'cx':
            a, b = qs
            x[b] ^= x[a]
            z[a] ^= z[b]
        elif op == 'cz':
            a, b = qs
            z[b] ^= x[a]
            z[a] ^= x[b]
        elif op == 'r':
            x[qs[0]] = z[qs[0]] = 0
        elif op not in ('mz', 'i', 'correct_x', 'correct_z'):
            raise ValueError('unknown recovery operation')
        # Faults follow gates/preparations, and precede destructive reads.
        for q in qs:
            x[q] ^= 1 << variables
            z[q] ^= 1 << (variables+1)
            variables += 2
        if op == 'mz':
            q = qs[0]
            records[rec] ^= x[q]
            x[q] = z[q] = 0
    records[:14] = x[:7]+z[:7]
    width = (variables+7)//8
    bits = np.unpackbits(np.array([list(r.to_bytes(width, 'little')) for r in records], dtype=np.uint8),
                         axis=1, bitorder='little')[:, :variables]
    columns = np.packbits(bits.T, axis=1, bitorder='little')
    return [int.from_bytes(row.tobytes(), 'little') for row in columns], spans


def interpret(signature, flags, rounds):
    observed = [(signature >> (14+6*r)) & 63 for r in range(rounds)]
    rejected = [tuple(any(signature >> bit & 1 for bit in group) for group in candidates)
                for candidates in flags]
    # A single fault that rejects candidate A occurred strictly before any
    # data coupling, solely in A's preparation/verification. B is fresh and
    # fault-free. The actual selected path then leaves the data untouched.
    if any(a for a, _ in rejected):
        if any(b for a, b in rejected if a):
            return None
        observed, x, z = [0]*rounds, 0, 0
    else:
        x, z = signature & 127, (signature >> 7) & 127
    selected = next((a for a, b in zip(observed, observed[1:]) if a == b), None)
    if rounds == 1:
        selected = observed[0]
    if selected is None:
        return None
    dx, dz = correction(selected)
    return x ^ dx, z ^ dz


def single_faults(columns, spans):
    for index, (start, width) in enumerate(spans):
        for labels in itertools.product(range(4), repeat=width):
            if not any(labels):
                continue
            response = 0
            for j, label in enumerate(labels):
                if label & 1:
                    response ^= columns[start+2*j]
                if label & 2:
                    response ^= columns[start+2*j+1]
            yield index, list(labels), response


def audit(config=CONFIG):
    tape, flags, count = program(config)
    columns, spans = responses(tape, count)
    allowed = correctable()
    histogram, failures, cases = {}, [], 0
    for index, labels, response in single_faults(columns, spans):
        result = interpret(response, flags, config['rounds'])
        cases += 1
        key = 'abort' if result is None else f'{result[0]}:{result[1]}'
        histogram[key] = histogram.get(key, 0)+1
        if result not in allowed:
            failures.append([index, labels, key])
    return dict(config=config, locations=len(tape), one_location_faults=cases,
                primitive_slots=sum(op != 'i' for op, _, _ in tape),
                idle_locations=sum(op == 'i' for op, _, _ in tape),
                tape_sha256=digest(tape), response_sha256=digest(columns),
                residual_histogram=histogram, failures=failures)


if __name__ == '__main__':
    result = audit()
    print({k: v for k, v in result.items() if k not in ('residual_histogram', 'failures')})
    print('failures:', result['failures'][:5], 'total', len(result['failures']))

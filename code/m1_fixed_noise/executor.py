"""Direct adaptive Pauli execution, used to cross-check the fast fault census.

Unlike the symbolic census this actually chooses the second cat when the
first is rejected. Measurement parity is independent of the random even
cat-read string, so a representative even string suffices for Pauli frames.
"""

from .recovery import CHECKS, CONFIG, correction, program


def execute(fault=None, incoming=(0, 0), config=CONFIG, faults=()):
    tape, _, _ = program(config)
    if (type(incoming) not in (tuple, list) or len(incoming) != 2
            or any(type(v) is not int or not 0 <= v < 128 for v in incoming)
            or type(faults) not in (tuple, list)):
        raise ValueError('seven-qubit Pauli input and explicit fault list required')
    injected = {}
    for entry in list(faults)+([] if fault is None else [fault]):
        if type(entry) not in (tuple, list) or len(entry) != 2:
            raise ValueError('fault index and exact Pauli labels required')
        index, labels = entry
        if (type(index) is not int or not 0 <= index < len(tape) or index in injected
                or type(labels) not in (tuple, list) or len(labels) != len(tape[index][1])
                or any(type(label) is not int or not 0 <= label < 4 for label in labels)
                or not any(labels)):
            raise ValueError('unique in-range location and nonidentity supported Pauli required')
        injected[index] = labels
    x, z = incoming
    index, records, rejected = 0, [], []
    def location(op, qs):
        nonlocal x, z, index
        if op == 'h':
            q = qs[0]
            if (x >> q ^ z >> q) & 1:
                x ^= 1 << q
                z ^= 1 << q
        elif op == 'cx':
            a, b = qs
            x ^= ((x >> a) & 1) << b
            z ^= ((z >> b) & 1) << a
        elif op == 'cz':
            a, b = qs
            z ^= ((x >> a) & 1) << b
            z ^= ((x >> b) & 1) << a
        elif op == 'r':
            x &= ~(1 << qs[0])
            z &= ~(1 << qs[0])
        if index in injected:
            for q, label in zip(qs, injected[index]):
                x ^= (label & 1) << q
                z ^= (label >> 1) << q
        index += 1
        result = 0
        if op == 'mz':
            result = (x >> qs[0]) & 1
            x &= ~(1 << qs[0])
            z &= ~(1 << qs[0])
        return result
    def step(op, qs):
        result = location(op, qs)
        for q in range(16):
            if q not in qs:
                location('i', [q])
        return result
    abort = False
    for _ in range(config['rounds']):
        measured = 0
        for g in range(6):
            flags = []
            for start in (7, 11):
                for q in range(start, start+4):
                    step('r', [q])
                step('h', [start])
                for q in range(start+1, start+4):
                    step('cx', [start, q])
                outcomes = []
                for a, b in config['verification']:
                    step('r', [15])
                    step('cx', [start+a, 15])
                    step('cx', [start+b, 15])
                    outcomes.append(step('mz', [15]))
                flags.append(any(outcomes))
            rejected.append(flags)
            abort |= all(flags)
            start = 11 if flags[0] else 7
            support = [i for i in range(7) if CHECKS[g % 3] >> i & 1]
            for a, q in enumerate(support):
                step('cx' if g < 3 else 'cz', [start+a, q])
            bit = 0
            for q in range(start, start+4):
                step('h', [q])
                bit ^= step('mz', [q])
            measured |= bit << g
        records.append(measured)
    selected = next((a for a, b in zip(records, records[1:]) if a == b), None)
    if config['rounds'] == 1:
        selected = records[0]
    abort |= selected is None
    dx, dz = correction(selected or 0)
    x ^= dx
    z ^= dz
    for op in ('correct_x', 'correct_z'):
        for q in range(7):
            step(op, [q])
    # Even exhaustion leaves a complete local quantum output. At a lower
    # concatenation level its flag is diagnostic, not an application abort.
    return dict(residual=None if abort else (x & 127, z & 127),
                continued_frame=(x & 127, z & 127), local_failure=bool(abort),
                syndrome_history=records, candidate_rejections=rejected, locations=index)

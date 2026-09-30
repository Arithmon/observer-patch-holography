"""Full ideal cat-instrument contraction, independent of Pauli-frame replay.

Four active data qubits, two four-qubit cats and one verifier suffice. The
remaining three data qubits are spectators. Every cat-read/discard branch
is retained as a Kraus map on the full 16-dimensional active input space.
"""

import numpy as np


def need(condition, message):
    if not condition:
        raise ValueError(message)


def check_gadget(gates, generator):
    support = [j-1 for j in range(1, 8) if j & (1 << (generator % 3))]
    wires = support+list(range(7, 16))
    remap = {q: i for i, q in enumerate(wires)}
    indices = np.arange(1 << 13)
    state = np.zeros((1 << 13, 16))
    state[:16] = np.eye(16)
    read = set()
    verification_reads = 0
    for op, physical, _ in gates:
        need(all(q in remap for q in physical), 'correct stabilizer support')
        qs = [remap[q] for q in physical]
        need(not read.intersection(qs), 'no use after deferred destructive read')
        q = qs[0]
        if op == 'h':
            low = indices[(indices >> q & 1) == 0]
            high = low | (1 << q)
            a, b = state[low].copy(), state[high].copy()
            state[low], state[high] = (a+b)/np.sqrt(2), (a-b)/np.sqrt(2)
        elif op == 'cx':
            state = state[indices ^ (((indices >> q) & 1) << qs[1])]
        elif op == 'cz':
            state *= (1-2*((indices >> q) & (indices >> qs[1]) & 1))[:, None]
        elif op == 'r' or (op == 'mz' and physical == [15]):
            # Only initialized cats and the just-read verifier are reset.
            # The rejected ideal verification branch is the ZERO operator
            # on every data input, not an outcome selected and renormalized.
            need(np.max(np.abs(state[(indices >> q & 1) != 0])) < 2e-12,
                 'ideal reset/verification has no discarded nonzero branch')
            if op == 'mz':
                verification_reads += 1
        elif op == 'mz':
            need(physical[0] in range(7, 11), 'first accepted cat is read')
            read.add(q)
        else:
            raise ValueError('known cat-instrument operation')
    need(read == set(range(4, 8)) and verification_reads == 6,
         'all cat and verification outcomes present')
    tensor = state.reshape(2, 16, 16, 16, 16)
    need(np.max(np.abs(tensor[1])) < 2e-12, 'verifier ends in its reset state')
    identity = np.eye(16)
    stabilizer = identity[np.arange(16) ^ 15] if generator < 3 else np.diag(
        [(-1.)**j.bit_count() for j in range(16)])
    complete = np.zeros((16, 16))
    for first in range(16):
        for unused in range(16):
            actual = tensor[0, unused, first]
            expected = ((identity+(-1.)**first.bit_count()*stabilizer)/8
                        if unused in (0, 15) else np.zeros((16, 16)))
            need(np.max(np.abs(actual-expected)) < 2e-12, 'full cat Kraus output')
            complete += actual.T@actual
    need(np.max(np.abs(complete-identity)) < 2e-12, 'complete cat instrument')


def check_program(tape):
    # Split the first full syndrome round from the concrete schedule. Idle
    # locations have identity ideal action and remain in the fault census.
    pieces = []
    for gate in tape:
        op, qs, _ = gate
        if op == 'i':
            continue
        if op == 'r' and qs == [7]:
            if len(pieces) == 6:
                break
            pieces.append([])
        if pieces:
            pieces[-1].append(gate)
    need(len(pieces) == 6, 'all six stabilizer instruments')
    for generator, gates in enumerate(pieces):
        check_gadget(gates, generator)

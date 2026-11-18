"""Production adapter for the separate written-contract reconstruction.

No receipt is an oracle input. Run independent.py first for import isolation.
"""

from fractions import Fraction as F
import hashlib
import json
from math import comb, sqrt
from pathlib import Path

from . import independent as audit


def custody():
    from tools.pin_scope import bound_pins

    contract, protocol = audit.admit()
    root = Path(__file__).resolve().parents[2]
    # Proof digests document the audited snapshot. Editorial changes are
    # exempt from live custody checks under the shared repository policy.
    for name, digest in bound_pins(contract["audited_sha256"]).items():
        audit.require(hashlib.sha256((root / name).read_bytes()).hexdigest() == digest,
                      "audited source custody: " + name)
    return protocol


def own_data(data):
    return audit.Registers(data.buffers, data.addresses, data.couriers,
                           data.pointer_a, data.pointer_b, data.old_bits, data.second_address)


def producer_data(data, program):
    return program.Data(*data)


def config(background, variant, program):
    return program.Program(background, "uphill" if variant == "uphill" else "native",
                           variant != "no_communication", variant != "no_feedback",
                           variant != "no_publish")


def same_columns(actual, expected):
    audit.require(actual == expected, "complete-register isometry column mismatch")
    audit.require(audit.pointer(actual) == audit.pointer(expected), "public outcome partition mismatch")


def complete_weights(reported, original):
    audit.require(len(reported) == len(original) == 5, "all unfinished prefixes retained")
    audit.require(abs(sum(reported) - 1) < 2e-10, "complete-law normalization")
    audit.require(max(abs(x - y) for x, y in zip(reported, original)) < 2e-10,
                  "no silent completion postselection")


def rejected(predicate, action):
    try:
        action()
    except ValueError as error:
        audit.require(predicate in str(error), "mutation reached wrong predicate: " + str(error))
        return str(error)
    raise ValueError("mutation survived: " + predicate)


def mutations():
    bg = audit.BACKGROUNDS[0]
    good = audit.prefixes(1, 0, bg)[1]
    bad_word = audit.event(1, 1, bg)[0]
    bad = good._replace(sources=(0, bad_word, 0, 0, 0))
    swapped = rejected("complete-register isometry", lambda: same_columns(bad, good))
    source, address, memory = 0, 0, 1
    forward = audit.reversible(source, address, memory, bg)
    # Wrong inverse ignores the output memory and declares a blank input.
    wrong_inverse = (forward[0], 0)
    inverse = rejected("dirty inverse identity", lambda: audit.require(
        wrong_inverse == (source, memory), "dirty inverse identity"))
    # Overwrite the successor instead of XORing dirty input memory.
    overwritten = (audit.event(source, address, bg)[0], (source >> address) & 1)
    recovered = audit.reversible(overwritten[0], address, overwritten[1], bg, inverse=True)
    overwrite = rejected("dirty inverse identity", lambda: audit.require(
        recovered == (source, memory), "dirty inverse identity"))
    weights = [F(1, 2), F(1, 8), F(1, 16), F(1, 16), F(1, 4)]
    omitted = rejected("complete-law normalization", lambda: complete_weights([0, 0, 0, 0, weights[-1]], weights))
    postselected = rejected("no silent completion postselection", lambda: complete_weights([0, 0, 0, 0, 1], weights))
    coherence = audit.coherence_witness()
    dephased = rejected("retained coherent block", lambda: audit.require(
        coherence["block_density_difference_frobenius"] == 0, "retained coherent block differs after dephasing"))
    return {
        "swapped_address_bit": {"input": [1, 0, list(bg)], "native_target": good.sources[1],
                                "mutated_target": bad_word, "predicate": swapped},
        "wrong_inverse": {"dirty_input": [source, address, memory], "forward": forward,
                          "wrong_recovery": wrong_inverse, "predicate": inverse},
        "dirty_overwrite": {"dirty_input": [source, address, memory], "wrong_recovery": recovered,
                            "predicate": overwrite},
        "missing_unfinished": {"retained_trace": "1/4", "predicate": omitted},
        "silent_postselection": {"completion_weight": "1/4", "predicate": postselected},
        "private_dephasing": {"public_probabilities_equal": True, "predicate": dephased,
                              "density_difference_frobenius": coherence["block_density_difference_frobenius"]}}


def small_rail(protocol, program, packet):
    """Literal own gates, forward and adjoint, on a closed prepared span.

    The reference is a two-column spectator, never included in the generator.
    Dense diagonalization is confined to the 48 by 48 scalar clock.
    """
    import numpy as np
    from scipy.sparse import coo_matrix, diags
    from scipy.sparse.linalg import expm_multiply

    n, offset, m, time = (protocol["literal_check"][k] for k in ("sites", "offset", "order", "time"))
    start = offset + m
    stage_at = lambda position: min(4, max(0, position - start))
    phi = np.zeros(n, complex)
    for j in range(m + 1):
        phi[offset + j] = (-1j) ** j * comb(m, j) / sqrt(comb(2 * m, m))
    # Declared complex entangled witness, with source0/source3 and reference2.
    coefficients = {0: np.array([.5, .5j]), 3: np.array([.5, .5])}
    histories = {(x, a): audit.prefixes(x, a, audit.BACKGROUNDS[0])
                 for x in coefficients for a in range(12)}
    states = [(position, history[stage_at(position)])
              for history in histories.values() for position in range(n)]
    lookup = {state: i for i, state in enumerate(states)}
    audit.require(len(lookup) == len(states) == n * 24, "small reducing span injectivity")
    row, col, values = [], [], []
    initial = np.zeros((len(states), 2), complex)
    producer_config = program.Program()
    rail = packet.Rail(n, offset, m)
    for column, (position, data) in enumerate(states):
        own_column = {}
        for step in (-1, 1):
            target = position + step
            if 0 <= target < n:
                edge = min(position, target) - start
                changed = audit.gate(data, edge, audit.BACKGROUNDS[0], inverse=step < 0) if 0 <= edge < 4 else data
                audit.require((target, changed) in lookup, "literal column closure")
                own_column[target, changed] = .5
                row.append(lookup[target, changed])
                col.append(column)
                values.append(.5)
        actual = packet.column(rail, producer_config, (position, producer_data(data, program)))
        converted = {(pos, own_data(out)): value for (pos, out), value in actual.items()}
        audit.require(converted == own_column, "literal full column comparison")
    for (source, address), history in histories.items():
        for position in range(offset, start + 1):
            initial[lookup[position, history[0]]] = phi[position] * coefficients[source] / sqrt(12)
    h = coo_matrix((values, (row, col)), shape=(len(states), len(states))).tocsr()
    audit.require((h - h.conjugate().T).nnz == 0, "literal adjoints")
    literal = expm_multiply(-1j * time * h, initial)
    clock = diags([np.full(n - 1, .5)] * 2, [-1, 1]).toarray()
    eigenvalues, vectors = np.linalg.eigh(clock)
    sine_eigenvalues = np.cos(np.pi * np.arange(1, n + 1) / (n + 1))
    spectrum_error = float(np.max(abs(eigenvalues - sine_eigenvalues[::-1])))
    alpha = vectors @ (np.exp(-1j * time * eigenvalues) * (vectors.T @ phi))
    analytic = np.zeros_like(initial)
    for (source, address), history in histories.items():
        for position in range(n):
            analytic[lookup[position, history[stage_at(position)]]] = alpha[position] * coefficients[source] / sqrt(12)
    residual = float(np.linalg.norm(literal - analytic))
    producer_clock_error = float(np.linalg.norm(alpha - packet.evolve_clock(rail, time)))
    weights = [float(sum(abs(alpha[pos]) ** 2 for pos in range(n) if stage_at(pos) == j)) for j in range(5)]
    complete_weights(weights, weights)
    omission = rejected("complete-law normalization", lambda: complete_weights([0, 0, 0, 0, weights[-1]], weights))
    postselection = rejected("postselection", lambda: complete_weights([0, 0, 0, 0, 1], weights))
    audit.require(max(residual, spectrum_error, producer_clock_error) < 2e-10, "small literal propagation residual")
    # Retain unnormalized pointer-only blocks; never measure private histories.
    blocks, reference = {}, np.zeros((2, 2), complex)
    mixture_blocks = {}
    for position in range(n):
        selected = [(data, literal[index]) for (pos, data), index in lookup.items() if pos == position]
        groups = {}
        for data, amplitudes in selected:
            groups.setdefault(audit.pointer(data), []).append((data, amplitudes))
            reference += np.outer(amplitudes, amplitudes.conjugate())
        for pointer, entries in groups.items():
            block = blocks.setdefault(pointer, {})
            expected_block = mixture_blocks.setdefault(pointer, {})
            for left, x in entries:
                for right, y in entries:
                    for r in range(2):
                        for s in range(2):
                            key = (audit.identity(left), r, audit.identity(right), s)
                            block[key] = block.get(key, 0j) + x[r] * y[s].conjugate()
    for j, weight in enumerate(weights):
        groups = {}
        for (source, address), history in histories.items():
            data = history[j]
            groups.setdefault(audit.pointer(data), []).append((data, coefficients[source] / sqrt(12)))
        for pointer, entries in groups.items():
            block = mixture_blocks.setdefault(pointer, {})
            for left, x in entries:
                for right, y in entries:
                    for r in range(2):
                        for s in range(2):
                            key = (audit.identity(left), r, audit.identity(right), s)
                            block[key] = block.get(key, 0j) + weight * x[r] * y[s].conjugate()
    block_error = max(abs(blocks.get(o, {}).get(k, 0j) - mixture_blocks.get(o, {}).get(k, 0j))
                      for o in blocks.keys() | mixture_blocks.keys()
                      for k in blocks.get(o, {}).keys() | mixture_blocks.get(o, {}).keys())
    audit.require(block_error < 2e-10, "unnormalized complete pointer mixture")
    initial_reference = sum((np.outer(v, v.conjugate()) for v in coefficients.values()), np.zeros((2, 2), complex))
    audit.require(np.max(abs(reference - initial_reference)) < 2e-10, "small rail untouched reference")
    return {"sites": n, "joint_dimension": len(states), "reference_dimension": 2,
            "time": time, "literal_vector_l2": residual, "spectrum_error": spectrum_error,
            "producer_clock_vector_l2": producer_clock_error,
            "unnormalized_block_max_error": block_error, "prefix_weights": weights,
            "completion_weight": weights[-1], "half_diamond_from_flag_argument": 1 - weights[-1],
            "missing_unfinished_predicate": omission, "silent_postselection_predicate": postselection,
            "unconditional_S4_source85_TV": {"feedback": weights[-1] / 2, "uphill": weights[-1]},
            "scope": "N48 composition diagnostic only; no numerical transfer to N4096"}


def compare():
    protocol = custody()
    # These imports are restricted to this adapter. The independent oracle
    # has completed its standalone run before any adapter execution.
    from model.native_memory_clock import native
    from . import packet, program

    entries = 0
    for background in audit.BACKGROUNDS:
        for variant in audit.VARIANTS:
            producer_config = config(background, variant, program)
            for source in range(4096):
                for address in range(12):
                    expected = audit.prefixes(source, address, background, variant)
                    actual = program.prefixes(producer_config, program.preparation(source, address))
                    for wanted, observed in zip(expected, actual):
                        same_columns(own_data(observed), wanted)
                        entries += 1
            for data in audit.dirty_states():
                for j in range(4):
                    for reverse in (False, True):
                        observed = program.edge(producer_config, producer_data(data, program), j, -1 if reverse else 1)
                        same_columns(own_data(observed), audit.gate(data, j, background, variant, reverse))
        for law in ("native", "uphill"):
            for source in range(4096):
                for address in range(12):
                    expected, old, current = audit.event(source, address, background, law)
                    observed = native.native_event(source, address, background, law)
                    audit.require((observed.target, observed.old_bit, observed.current) == (expected, old, current),
                                  "native event mismatch")
                    audit.require(observed.before_loads == audit.loads(source, background)
                                  and observed.after_loads == audit.loads(expected, background), "native load mismatch")
                    for memory in range(2):
                        for reverse in (False, True):
                            actual = (native.inverse if reverse else native.forward)(source, address, memory, background, law)
                            audit.require(actual == audit.reversible(source, address, memory, background, law, reverse),
                                          "dirty event producer comparison")
    bound = audit.clock_bound(protocol)
    for direction in (1, -1):
        vector = packet.packet(packet.Rail(), direction)
        expected = [0j] * 4096
        for j in range(65):
            expected[2048 + j] = (-1j * direction) ** j * comb(64, j) / sqrt(comb(128, 64))
        audit.require(all(x == y for x, y in zip(vector, expected)), "actual signed binary packet")
    certificate = packet.interval_certificate(packet.Rail())
    frac = lambda r: F(int(r["numerator"]), int(r["denominator"]))
    audit.require(frac(certificate["cantelli_failure"]) == F(bound["cantelli"]), "producer rational readiness")
    return {"status": "VERIFIED_WITHIN_SCOPE", "prepared_register_columns_compared": entries,
            "arbitrary_reference": "all columns and public partitions equal with common real amplitude 1/sqrt(12)",
            "small_rail": small_rail(protocol, program, packet), "mutations": mutations()}


if __name__ == "__main__":
    print(json.dumps(compare(), indent=2))

"""Independent affine-matrix replay; never imports the tape producer."""
from fractions import Fraction as F
import hashlib
if __package__:
    from . import codec
else:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_encoded_memory import codec


def need(condition, label):
    if not condition:
        raise ValueError(label)


def expected_cases():
    # Independent experimental contract: not data selected by the producer.
    return [
        ("empty_route", "move", 0, "3/2", "captured_square"),
        ("one_copy", "copy", 1, "3/2", "captured_square"),
        ("reread_positive", "read", 16, "3/2", "captured_square"),
        ("reread_negative", "read", 16, "-3/2", "captured_square"),
        ("reuse_positive", "move", 64, "3/2", "captured_square"),
        ("reuse_negative", "move", 64, "-3/2", "captured_square"),
        ("line_positive", "move", 24, "3/2", "declared_ladder"),
        ("line_negative", "move", 24, "-3/2", "declared_ladder"),
    ]


def edge_for(mode, topology, k):
    cycle, stage = divmod(k, 3)
    s = cycle if topology == "declared_ladder" else cycle % 2
    if mode in ("read", "copy"):
        s = 0
    t = s+1 if topology == "declared_ladder" else 1-s
    if stage < 2:
        return 2*s+stage, 2*t+stage
    cleared = t if mode == "read" else s
    return 2*cleared, 2*cleared+1


def values(matrix, initial):
    return [sum((weight*x for weight, x in zip(row, initial) if weight), F(0))
            for row in matrix]


def check_case(case, spec, captured_edges):
    name, mode, depth, signed, topology = spec
    cells = 25 if topology == "declared_ladder" else 2
    n = 2*cells
    a, b = F(signed), F(2)
    global_ports = [0, 1, 5, 9] if topology == "captured_square" else None
    initial = [b+a, b-a]+[b]*(n-2)
    for value in case["initial"]:
        codec.rational(value)
    codec.equal(case["initial"], list(map(str, initial)), "prepared input")
    count = 2 if mode == "copy" else 3*depth
    need(len(case["events"]) == count, "missing or extra native operations")
    matrix = [[F(i == j) for j in range(n)] for i in range(n)]
    writers = [f"init:{i}" for i in range(n)]
    loss, peak_bits = F(0), 0
    output_events, checkpoints = [], []

    def check_state(x):
        nonlocal peak_bits
        need(all(v >= 0 for v in x), "negative raw load")
        need(sum(x) == n*b, "total load changed")
        q = sum((v-b)**2 for v in x)
        need(q+loss == 2*a*a, "quadratic ledger changed")
        peak_bits = max(peak_bits, sum(abs(v.numerator).bit_length()
                                       + v.denominator.bit_length() for v in x))
    check_state(initial)
    for k, event in enumerate(case["events"]):
        u, v = edge_for(mode, topology, k)
        if global_ports is not None:
            need(tuple(sorted((global_ports[u], global_ports[v]))) in captured_edges,
                 "operation is not a captured seam")
        else:
            cu, ru = divmod(u, 2)
            cv, rv = divmod(v, 2)
            need((cu == cv and ru != rv) or (abs(cu-cv) == 1 and ru == rv),
                 "operation is outside declared ladder")
        previous = values(matrix, initial)
        row = [(matrix[u][j]+matrix[v][j])/2 for j in range(n)]
        matrix[u], matrix[v] = row.copy(), row.copy()
        current = values(matrix, initial)
        defect = sum((z-b)**2 for z in previous)-sum((z-b)**2 for z in current)
        need(defect == (previous[u]-previous[v])**2/2, "wrong native defect")
        output_events.append({"index": k, "op": "mean", "edge": [u, v],
                              "inputs": [str(previous[u]), str(previous[v])],
                              "writers": [writers[u], writers[v]],
                              "output": str(current[u]), "quadratic_loss": str(defect)})
        codec.equal(event, output_events[-1], "native event and writer custody")
        writers[u] = writers[v] = f"mean:{k}"
        loss += defect
        check_state(current)
        if k % 3 == 1:
            cycle = k//3
            s, t = u//2, v//2
            amp = a/F(2)**(cycle+1)
            expected = [b]*n
            for cell in (s, t):
                expected[2*cell], expected[2*cell+1] = b+amp, b-amp
            need(current == expected, "logical copy or blank ancilla failed")
            checkpoints.append({"after_mean": k, "rails": list(map(str, current))})
    final = values(matrix, initial)
    amp = a/F(2)**depth
    active = [0, 1] if mode == "copy" else [
        0 if mode == "read" else depth if topology == "declared_ladder" else depth % 2]
    expected = [b]*n
    for cell in active:
        expected[2*cell], expected[2*cell+1] = b+amp, b-amp
    need(final == expected, "cleanup or attenuation formula failed")
    reconstructed = {
        "name": name, "mode": mode, "copies": depth, "topology": topology,
        "baseline": "2", "amplitude": signed, "cells": cells,
        "global_ports": global_ports, "initial": list(map(str, initial)),
        "events": output_events, "copy_checkpoints": checkpoints,
        "final": list(map(str, final)), "final_writers": writers,
        "resources": {"scalar_registers": n, "preparation_writes": n,
                      "means": count, "scalar_reads": 2*count,
                      "scalar_writes": n+2*count, "peak_scalar_value_bits": peak_bits,
                      "initial_quadratic": str(2*a*a),
                      "final_quadratic": str(sum((z-b)**2 for z in final)),
                      "loss": str(loss), "total_load": str(n*b)}}
    codec.equal(case, reconstructed, "complete experiment and resource account")
    return count, len(checkpoints)


def check_margins(rows):
    expected = []
    for h in [0, 1, 4, 12, 14, 15, 24, 64]:
        amp = F(3, 2**(h+1))
        error = F(17+3*h, 2**20)
        # [2,2] is within exactly amp of BOTH opposite final codewords.
        for sign in (-1, 1):
            need(abs(F(2)-(F(2)+sign*amp)) == amp, "ambiguity control")
        expected.append({"hops": h, "amplitude": str(amp),
                         "preparation_bound": "1/1048576", "per_mean_bound": "1/1048576",
                         "readout_bound": "1/65536", "total_bound": str(error),
                         "sufficient_sign_margin": amp > error,
                         "ambiguous_observation": ["2", "2"],
                         "opposite_bits_terminal_ambiguity_radius": str(amp)})
    codec.equal(rows, expected, "analytic margins and ambiguity")


def verify(packet):
    codec.equal(sorted(packet), sorted(["schema", "source_sha256", "executions",
                                      "analytic_margin_controls_not_noisy_executions"]),
                "top-level schema")
    codec.equal(packet["schema"], "oph.source_encoded_memory.controls.v1", "schema version")
    codec.equal(packet["source_sha256"], codec.pins(), "source pins")
    support = codec.load(codec.ROOT/codec.SUPPORT)
    edges = {tuple(sorted(e)) for e in support["intra_carrier_seams"]}
    need(len(packet["executions"]) == 8, "experiment count")
    counts = [check_case(c, s, edges) for c, s in zip(packet["executions"], expected_cases())]
    check_margins(packet["analytic_margin_controls_not_noisy_executions"])
    return {"schema": "oph.source_encoded_memory.receipt.v1",
            "controls_sha256": hashlib.sha256(codec.canonical(packet)).hexdigest(),
            "source_sha256": codec.pins(), "executions": len(counts),
            "executed_means": sum(c[0] for c in counts),
            "checked_copy_states": sum(c[1] for c in counts),
            "captured_scope": "four ports in one carrier; no inter-carrier route",
            "declared_ladder_scope": "25 prepared cells; no W12 embedding claimed",
            "noise_rows": "analytic sufficient bounds and terminal ambiguity; not noisy executions",
            "physical_M1_derived": False, "amplitude_refresh_derived": False,
            "physical_comparator_derived": False, "scientific_status_promoted": False}


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    data = (codec.HERE/"controls.json").read_bytes()
    packet = codec.load(codec.HERE/"controls.json")
    need(data == codec.canonical(packet), "noncanonical controls bytes")
    receipt = codec.canonical(verify(packet))
    path = codec.HERE/"receipt.json"
    if args.write_receipt:
        path.write_bytes(receipt)
    else:
        need(path.read_bytes() == receipt, "retained receipt mismatch")
    print("Verified 8 executions, 626 native means and 209 copy states.")


if __name__ == "__main__":
    main()

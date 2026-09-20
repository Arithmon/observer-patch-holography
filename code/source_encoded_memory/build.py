"""Produce complete rational tapes. The independent verifier certifies them."""
from fractions import Fraction as F
if __package__:
    from . import codec
else:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_encoded_memory import codec


def specifications():
    return [
        ("empty_route", "move", 0, F(3, 2), "captured_square"),
        ("one_copy", "copy", 1, F(3, 2), "captured_square"),
        ("reread_positive", "read", 16, F(3, 2), "captured_square"),
        ("reread_negative", "read", 16, F(-3, 2), "captured_square"),
        ("reuse_positive", "move", 64, F(3, 2), "captured_square"),
        ("reuse_negative", "move", 64, F(-3, 2), "captured_square"),
        ("line_positive", "move", 24, F(3, 2), "declared_ladder"),
        ("line_negative", "move", 24, F(-3, 2), "declared_ladder"),
    ]


def execute(name, mode, depth, amplitude, support):
    cells = depth+1 if support == "declared_ladder" else 2
    baseline = F(2)
    x = [baseline]*(2*cells)
    x[0] += amplitude
    x[1] -= amplitude
    initial = x.copy()
    writers = [f"init:{i}" for i in range(len(x))]
    events, checkpoints = [], []
    loss, peak_bits = F(0), 0
    word = []
    source = 0
    for cycle in range(depth):
        target = source+1 if support == "declared_ladder" else 1-source
        word += [(2*source, 2*target), (2*source+1, 2*target+1)]
        if mode != "copy":
            cleared = target if mode == "read" else source
            word.append((2*cleared, 2*cleared+1))
        if mode == "move":
            source = target

    def charge():
        nonlocal peak_bits
        # A scalar-value encoding count only, not total machine storage.
        peak_bits = max(peak_bits, sum(abs(v.numerator).bit_length()
                                       + v.denominator.bit_length() for v in x))
    charge()
    for k, (u, v) in enumerate(word):
        before = [x[u], x[v]]
        output = sum(before)/2
        defect = (before[0]-before[1])**2/2
        events.append({"index": k, "op": "mean", "edge": [u, v],
                       "inputs": list(map(str, before)),
                       "writers": [writers[u], writers[v]],
                       "output": str(output), "quadratic_loss": str(defect)})
        x[u] = x[v] = output
        writers[u] = writers[v] = f"mean:{k}"
        loss += defect
        charge()
        if k % 3 == 1:
            checkpoints.append({"after_mean": k, "rails": list(map(str, x))})
    return {"name": name, "mode": mode, "copies": depth,
            "topology": support, "baseline": str(baseline),
            "amplitude": str(amplitude), "cells": cells,
            "global_ports": [0, 1, 5, 9] if support == "captured_square" else None,
            "initial": list(map(str, initial)), "events": events,
            "copy_checkpoints": checkpoints,
            "final": list(map(str, x)), "final_writers": writers,
            "resources": {"scalar_registers": len(x), "preparation_writes": len(x),
                          "means": len(word), "scalar_reads": 2*len(word),
                          "scalar_writes": len(x)+2*len(word),
                          "peak_scalar_value_bits": peak_bits,
                          "initial_quadratic": str(sum((v-baseline)**2 for v in initial)),
                          "final_quadratic": str(sum((v-baseline)**2 for v in x)),
                          "loss": str(loss), "total_load": str(sum(x))}}


def margins():
    rows = []
    for h in (0, 1, 4, 12, 14, 15, 24, 64):
        amplitude = F(3, 2)/2**h
        preparation, disturbance, readout = F(1, 2**20), F(1, 2**20), F(1, 2**16)
        bound = preparation+3*h*disturbance+readout
        rows.append({"hops": h, "amplitude": str(amplitude),
                     "preparation_bound": str(preparation),
                     "per_mean_bound": str(disturbance), "readout_bound": str(readout),
                     "total_bound": str(bound), "sufficient_sign_margin": amplitude > bound,
                     "ambiguous_observation": ["2", "2"],
                     "opposite_bits_terminal_ambiguity_radius": str(amplitude)})
    return rows


def build():
    return {"schema": "oph.source_encoded_memory.controls.v1",
            "source_sha256": codec.pins(),
            "executions": [execute(*s) for s in specifications()],
            "analytic_margin_controls_not_noisy_executions": margins()}


if __name__ == "__main__":
    (codec.HERE/"controls.json").write_bytes(codec.canonical(build()))
    print("Wrote controls.json; independent verification is required.")

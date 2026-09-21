#!/usr/bin/env python3
"""Reproduce the finite three-observer operator-join evidence packet.

No experiment, quantum outcome history or physical clock is produced.
Universal operator statements are proved in TripleCarrierOperatorJoin.lean.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARENTS = (
    "Lean/QFT/SourceOperatorGeneration.lean",
    "Lean/QFT/TripleCarrierJoin.lean",
    "Lean/QFT/TowerAnchoredCorrelation.lean",
)


def source_rows(root: Path = ROOT) -> list[list[int]]:
    text = (root / PARENTS[0]).read_text()
    rows = []
    for observer in (86, 88, 247):
        block = text.split(f"def obs{observer} :", 1)[1].split("  cnt :=", 1)[0]
        raw = block.split("  path := !", 1)[1].strip()
        rows.append(ast.literal_eval(raw))
    return rows


def swap(n: int, a: int, b: int) -> list[int]:
    p = list(range(n))
    p[a], p[b] = p[b], p[a]
    return p


def counts(rows: list[list[int]]) -> list[list[int]]:
    return [list(key) + [value] for key, value in sorted(Counter(zip(*rows)).items())]


def produce(root: Path = ROOT) -> dict:
    rows = source_rows(root)
    return {
        "schema": "oph.source_operator_join.v1",
        "parents": {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in PARENTS},
        "carrier": {
            "observers": [86, 88, 247], "slot_dimensions": [13, 14, 13],
            "pair_dimensions": [182, 169], "triple_dimension": 2366,
            "common_full_matrix_dimension": 13,
            "common_complex_vector_dimension": 169,
            "joint_complex_vector_dimension": 5597956,
        },
        "embeddings": {"pair88_axes": [0, 1], "pair247_axes": [0, 2],
                       "unused_factor": "identity", "common_axes": [0]},
        "state": {
            "denominator": 32,
            "triple_counts": counts(rows),
            "pair88_counts": counts(rows[:2]),
            "pair247_counts": counts([rows[0], rows[2]]),
        },
        "transitions": [
            {"row": t, "next_row": t + 1,
             "permutations": [swap(n, r[t], r[t + 1])
                              for n, r in zip((13, 14, 13), rows)]}
            for t in range(31)
        ],
        "interpretation": {
            "tensor_assembly": "declared",
            "source_labels": "post_hoc_committed",
            "empirical_state": "source_counted_diagonal",
            "tower_state": "independent_uniform_not_identified",
            "admitted_generation": "full_source_operator_generators",
            "overlapping_pair_algebras_commute": False,
            "physical_regions": False, "physical_clock": False,
            "observed_quantum_outcomes": False,
            "continuum_refinement": False,
            "scalar64_instrument_same_carrier": False,
            "universal_algebraic_pushout": False,
        },
    }


def canonical(packet: dict) -> bytes:
    return (json.dumps(packet, sort_keys=True, indent=2) + "\n").encode()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_bytes(canonical(produce()))

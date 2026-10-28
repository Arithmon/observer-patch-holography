"""Exact stationarity test for the transported clock reference and native C/E.

The complete native repair law is tested without conditioning on a particular
history. This does not test the predecessor's branchwise transported-reference
identity, which uses an updated reference and retains each exported record.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path

from . import continuation, defect_history
from .joins_contract import binding, check_receipt, source_bindings
from ..observer import record_recovery as records

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "reference_balance_results.json"
CORES = ((0, 5), (0, 8))
PROFILE = (1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0)


def kernel(law: str) -> tuple[tuple[F, ...], ...]:
    """Every source state, seam and coin, including all zero-current slots."""
    if law not in ("C", "E"):
        raise ValueError("a complete native C or E law is required")
    basis = continuation.BASIS
    index = {state: i for i, state in enumerate(basis)}
    counts = [[0 for _ in basis] for _ in basis]
    for i, state in enumerate(basis):
        for seam in records.seams():
            for coin in (0, 1):
                output, _ = records.recorded_step(law, state, seam, coin)
                counts[i][index[tuple(output)]] += 1
    result = tuple(tuple(F(n, 60) for n in row) for row in counts)
    if any(sum(row) != 1 for row in result):
        raise ArithmeticError("native kernel lost a seam/coin outcome")
    return result


def row(core: tuple[int, int], transition: tuple[tuple[F, ...], ...]) -> dict:
    """The fixed clock reference under one unconditioned native repair."""
    if core not in CORES:
        raise ValueError("only the two preregistered core controls are admitted")
    basis, edges, *_ = defect_history.geometry()
    sigma = continuation.reference(core)
    if len(sigma) != len(basis) or sum(sigma) != 1:
        raise ArithmeticError("faithful clock reference is not normalized")
    if len(transition) != len(basis) or any(len(x) != len(basis) for x in transition):
        raise ValueError("a complete binary-shell native kernel is required")
    after = tuple(sum((sigma[i] * transition[i][j] for i in range(66)), F())
                  for j in range(66))
    drift = tuple(b-a for a, b in zip(sigma, after, strict=True))
    if sum(drift) != 0 or any(x <= 0 for x in sigma):
        raise ArithmeticError("native source or faithful reference changed")
    mode = tuple(sum(PROFILE[p] for p in pair) for pair in basis)
    mean = sum((sigma[i]*mode[i] for i in range(66)), F())
    variance = sum((sigma[i]*(mode[i]-mean)**2 for i in range(66)), F())
    form = sum((sigma[i]*transition[i][j]*(mode[j]-mode[i])**2
                for i in range(66) for j in range(66)), F()) / 2
    if variance <= 0 or form <= 0:
        raise ArithmeticError("the fixed-reference observable lost its response")
    violated = sum(sigma[i]*transition[i][j] != sigma[j]*transition[j][i]
                   for i, j in edges)
    return {
        "core": list(core),
        "clock_reference_weight_multiplicities": {
            str(weight): count for weight, count in sorted(Counter(sigma).items())
        },
        "one_attempt_total_variation_from_fixed_reference": str(sum(map(abs, drift))/2),
        "maximum_absolute_probability_drift": str(max(map(abs, drift))),
        "native_edges_violating_fixed_reference_detailed_balance": violated,
        "native_edges_total": len(edges),
        "mode_variance_at_fixed_reference": str(variance),
        "mode_quadratic_change_at_fixed_reference": str(form),
        "mode_quadratic_change_over_variance": str(form/variance),
        "fixed_reference_stationary": all(value == 0 for value in drift),
    }


def report() -> dict:
    pairs, edges, *_ = defect_history.geometry()
    c, e = kernel("C"), kernel("E")
    if c != e or any(c[i][j] != c[j][i] for i in range(66) for j in range(66)):
        raise ArithmeticError("the two complete population kernels are not symmetric")
    native_edges = {(i, j) for i, j in edges}
    if ({(i, j) for i in range(66) for j in range(i+1, 66) if c[i][j]}
            != native_edges or any(c[i][j] != F(1, 60) for i, j in edges)):
        raise ArithmeticError("source transition support differs from actual native incidence")
    uniform = (F(1, 66),)*66
    if any(sum((uniform[i]*c[i][j] for i in range(66)), F()) != uniform[j]
           for j in range(66)):
        raise ArithmeticError("uniform native control is not stationary")
    rows = [row(core, c) for core in CORES]
    if any(item["fixed_reference_stationary"] for item in rows):
        raise ArithmeticError("transported clock reference unexpectedly stationary")
    paths = (
        "model/clock_defect_energy/reference_balance.py",
        "model/clock_defect_energy/test_reference_balance.py",
        "model/clock_defect_energy/REFERENCE_BALANCE.md",
        "model/clock_defect_energy/continuation.py",
        "model/clock_defect_energy/defect_history.py",
        "model/observer/record_recovery.py",
        "model/native_modular_reference/history_transport.py",
        "Lean/Geometry/SourceActionMeasure.lean",
    )
    return {
        "schema": "oph.finite_theorem_joins.reference_balance.v1",
        "campaign": binding(),
        "domain": "unconditioned complete native C/E population repair on binary H66; fixed faithful diagonal reference for two declared cores",
        "native_states": len(pairs),
        "native_movement_edges": len(edges),
        "native_transition_probability_per_directed_edge": "1/60",
        "native_uniform_reference_is_stationary": True,
        "native_uniform_reference_unique_by_connected_reciprocity": True,
        "results": rows,
        "theorem_application": {
            "source": "Lean/Geometry/SourceActionMeasure.lean",
            "result": "normalized_symmetrizer_unique",
            "actual_hypotheses": "connected 66-state native movement graph; symmetric nonzero probability 1/60 on every edge; normalized positive uniform stationary weight",
            "conclusion": "the fixed nonuniform transported reference cannot be another reversible symmetrizer of the same bare population kernel",
        },
        "scope": {
            "branchwise_transported_reference_covariance_preserved": True,
            "native_law_or_clock_service_changed": False,
            "clock_or_energy_physical_identification_claimed": False,
            "integrated_gate_accepted": False,
        },
        "source_sha256": source_bindings(paths),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    else:
        check_receipt(OUTPUT, result)
        print("fixed-reference balance receipt OK")


if __name__ == "__main__":
    main()

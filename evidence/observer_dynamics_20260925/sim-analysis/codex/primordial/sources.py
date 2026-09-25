"""Conditional 3D history-mixture curvature spectra, ready for CAMB.

CAMB's scalar source is dimensionless Delta_R^2, not dimensional P_R.
This module performs no cosmological fitting and imports no CAMB package.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.special import gammainc, gammaincc

HERE = Path(__file__).resolve().parent


def _positive(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def gamma_window(k, beta: float, k_ir: float, k_uv: float):
    """Regularized integral between (k/k_uv)^2 and (k/k_ir)^2.

    Choose complementary Gamma tails in the UV to avoid subtracting two
    numbers rounded to one. Physical exponential tails may underflow to zero.
    """
    beta = _positive("beta", beta)
    k_ir = _positive("k_ir", k_ir)
    k_uv = _positive("k_uv", k_uv)
    if k_ir >= k_uv:
        raise ValueError("k_ir must be smaller than k_uv")
    k = np.asarray(k, dtype=float)
    if np.any(~np.isfinite(k)) or np.any(k <= 0):
        raise ValueError("all wave numbers must be positive and finite")
    with np.errstate(over="ignore"):
        lower, upper = (k / k_uv) ** 2, (k / k_ir) ** 2
    window = np.where(
        lower <= beta,
        gammainc(beta, upper) - gammainc(beta, lower),
        gammaincc(beta, lower) - gammaincc(beta, upper),
    )
    return np.maximum(window, 0.0)


def history_delta2(k, *, As: float, beta: float, k_ir: float, k_uv: float,
                   k_pivot: float = 0.05):
    """Dimensionless curvature power, normalized to As at the pivot.

    P_R(k) is proportional to integral t^(beta-1) exp(-2 kappa k^2 t) dt.
    k_ir=(2 kappa t_max)^(-1/2), k_uv=(2 kappa t_min)^(-1/2).
    Wave numbers must all use the same units (Mpc^-1 for CAMB).
    """
    As, k_pivot = _positive("As", As), _positive("k_pivot", k_pivot)
    window = gamma_window(k, beta, k_ir, k_uv)
    pivot_window = float(gamma_window(k_pivot, beta, k_ir, k_uv))
    if pivot_window <= 0:
        raise ValueError("pivot window underflows; choose a resolved pivot or cutoffs")
    return As * (np.asarray(k, dtype=float) / k_pivot) ** (3 - 2 * beta) * window / pivot_window


def powerlaw_delta2(k, *, As: float, ns: float, k_pivot: float = 0.05):
    As, k_pivot = _positive("As", As), _positive("k_pivot", k_pivot)
    k = np.asarray(k, dtype=float)
    if np.any(~np.isfinite(k)) or np.any(k <= 0) or not math.isfinite(ns):
        raise ValueError("wave numbers must be positive and finite; ns must be finite")
    return As * (k / k_pivot) ** (ns - 1)


def dimensional_curvature_power(k, delta2):
    """P_R in Mpc^3 if k is Mpc^-1; never pass this value to CAMB."""
    k = np.asarray(k, dtype=float)
    if np.any(~np.isfinite(k)) or np.any(k <= 0):
        raise ValueError("wave numbers must be positive and finite")
    return 2 * np.pi**2 * np.asarray(delta2) / k**3


def load_spec(path: Path | str = HERE / "spec.json") -> dict:
    return json.loads(Path(path).read_text())


def source_callable(candidate_id: str, *, spec: dict | None = None,
                    numerical_floor: float | None = None) -> Callable:
    """Return vectorized Delta_R^2(k) for CAMB set_initial_power_function.

    The optional floor is explicit and is for log interpolation only. The
    exact analytical spectrum is obtained with numerical_floor=None.
    """
    spec = load_spec() if spec is None else spec
    candidates = {row["id"]: row for row in spec["candidates"]}
    if candidate_id not in candidates:
        raise KeyError(f"unknown candidate {candidate_id!r}")
    candidate, baseline = candidates[candidate_id], spec["baseline"]
    kwargs = {"As": baseline["As"],
              "k_pivot": baseline["k_pivot_Mpc_inverse"]}
    if numerical_floor is not None:
        numerical_floor = _positive("numerical_floor", numerical_floor)

    def evaluate(k):
        if candidate["family"] == "powerlaw":
            result = powerlaw_delta2(k, ns=baseline["ns"], **kwargs)
        elif candidate["family"] == "finite_history":
            result = history_delta2(k, beta=candidate["beta"],
                                    k_ir=candidate["k_ir_Mpc_inverse"],
                                    k_uv=candidate["k_uv_Mpc_inverse"], **kwargs)
        else:
            raise ValueError(f"unknown family {candidate['family']!r}")
        return result if numerical_floor is None else np.maximum(result, numerical_floor)

    return evaluate


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_receipts() -> tuple[dict, dict]:
    spec = load_spec()
    params = (HERE / spec["baseline"]["parameter_source"]).resolve()
    if sha256(params) != spec["baseline"]["parameter_source_sha256"]:
        raise ValueError("official parameter source hash mismatch")
    grid = spec["source_grid"]
    k = np.geomspace(grid["k_min_Mpc_inverse"], grid["k_max_Mpc_inverse"], grid["log_samples"])
    pivot = spec["baseline"]["k_pivot_Mpc_inverse"]
    if grid["insert_pivot"]:
        k = np.unique(np.append(k, pivot))
    floor = grid["numerical_floor"]
    values, summaries = {}, {}
    probe_k = np.array([1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 0.002, 0.05, 0.15, 0.2, 1.0])
    for row in spec["candidates"]:
        key = row["id"]
        f = source_callable(key, spec=spec)
        exact = f(k)
        values[key] = np.maximum(exact, floor).tolist()
        summaries[key] = {
            "delta2_at_pivot": float(f(pivot)),
            "floored_grid_points": int(np.count_nonzero(exact < floor)),
            "minimum_exact_grid_delta2": float(np.min(exact)),
            "maximum_exact_grid_delta2": float(np.max(exact)),
            "probe_exact_delta2": f(probe_k).tolist(),
        }
        if row["family"] == "finite_history":
            summaries[key].update({
                "intermediate_ns": 4 - 2 * row["beta"],
                "age_variance_measure_exponent": row["beta"] - 1,
                "t_max_over_t_min": (row["k_uv_Mpc_inverse"] / row["k_ir_Mpc_inverse"]) ** 2,
                "pivot_gamma_window": float(gamma_window(pivot, row["beta"], row["k_ir_Mpc_inverse"], row["k_uv_Mpc_inverse"])),
            })
    receipt = {
        "schema": "oph.conditional-primordial-source.receipt.v1",
        "source_sha256": sha256(Path(__file__)), "spec_sha256": sha256(HERE / "spec.json"),
        "parameter_source_sha256": sha256(params),
        "As": spec["baseline"]["As"],
        "grid_points": len(k), "numerical_floor": floor,
        "probe_k_Mpc_inverse": probe_k.tolist(), "candidates": summaries,
        "interpretation": "Analytical extra-assumption source construction only; no Boltzmann output or data fit is contained in this receipt."
    }
    spectra = {"schema": "oph.dimensionless-curvature-grid.v1", "spec_sha256": receipt["spec_sha256"],
               "k_Mpc_inverse": k.tolist(), "numerical_floor": floor, "delta_R_squared": values}
    return receipt, spectra


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    receipt, spectra = build_receipts()
    for name, value in [("source_receipt.json", receipt), ("source_grid.json", spectra)]:
        target = HERE / name
        if args.verify:
            if json.loads(target.read_text()) != value:
                raise SystemExit(f"receipt differs: {target}")
        else:
            target.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    print("Verified source receipt and grid" if args.verify else "Wrote source receipt and grid")


if __name__ == "__main__":
    main()

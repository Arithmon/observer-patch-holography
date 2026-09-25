"""Small, reproducible audit of actual CAMB primordial interpolation.

No background, transfer, Boltzmann or nonlinear calculation is performed.
The grid combines dense logarithmic probes with samples inside every
initial interpolation interval. Positivity here is a finite sampling check,
not a mathematical guarantee between every pair of probes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "primordial"))
from sources import source_callable  # noqa: E402


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def probe_case(knots, floor, kind="max", *, dense_samples=180001):
    import camb

    exact = source_callable("history_tilt_calibrated_uv")
    if kind == "max":
        source = lambda k: np.maximum(exact(k), floor)
    elif kind == "add":
        source = lambda k: exact(k) + floor
    else:
        raise ValueError("Unknown regularization")
    ks = np.geomspace(1e-7, 100, knots)
    interior = (ks[:-1, None] + (ks[1:] - ks[:-1])[:, None]
                * np.linspace(0, 1, 17)[None, :]).ravel()
    probe = np.unique(np.concatenate([np.geomspace(1e-7, 100, dense_samples), interior]))
    params = camb.CAMBparams()
    params.set_initial_power_function(source, kmin=1e-7, kmax=100,
                                     N_min=knots, rtol=5e-6,
                                     effective_ns_for_nonlinear=0.9660499)
    interpolated = params.scalar_power(probe)
    target = source(probe)
    if np.any(~np.isfinite(interpolated)):
        raise ValueError("Nonfinite CAMB interpolant")
    negative = interpolated < 0
    core = (probe >= 1e-5) & (probe <= 0.5)
    return {
        "N_min": knots, "floor": floor, "regularization": kind,
        "samples": len(probe), "min_spline": float(interpolated.min()),
        "negative_count": int(negative.sum()),
        "first_negative_k_Mpc_inverse": float(probe[negative].min()) if negative.any() else None,
        "core_relative_interpolation_error_max": float(np.max(np.abs(interpolated[core] / target[core] - 1))),
        "max_absolute_interpolation_error": float(np.max(np.abs(interpolated - target))),
    }


def build_audit():
    import camb

    expected = json.loads((HERE / "spec.json").read_text())["camb_version"]
    if camb.__version__ != expected:
        raise ValueError(f"Expected CAMB {expected}, found {camb.__version__}")
    cases = [(1800, 1e-300, "max"), (1800, 1e-40, "max"),
             (1800, 1e-50, "max"), (1800, 1e-40, "add"),
             (1800, 1e-50, "add"), (7200, 1e-300, "max"),
             (7200, 1e-40, "max"), (7200, 1e-50, "max")]
    return {
        "schema": "oph.camb-interpolation-audit.v1",
        "camb": camb.__version__, "numpy": np.__version__,
        "inputs_sha256": {name: sha(HERE / name) for name in [
            "interpolation_audit.py", "spec.json", "../primordial/sources.py",
            "../primordial/spec.json"]},
        "probe": {"dense_log_samples": 180001, "samples_per_initial_interval": 17,
                  "k_min_Mpc_inverse": 1e-7, "k_max_Mpc_inverse": 100,
                  "core_k_range_Mpc_inverse": [1e-5, 0.5]},
        "cases": [probe_case(*case) for case in cases],
        "scope": "Actual CAMB scalar_power interpolation only. No cosmological evolution. Finite-grid positivity is not an interval proof; floors are numerical regularizations, not physical source terms.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    result = build_audit()
    target = HERE / "interpolation_audit.json"
    if args.verify:
        if json.loads(target.read_text()) != result:
            raise SystemExit("Interpolation receipt differs; inspect before replacing it")
    else:
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print("Verified interpolation audit" if args.verify else "Wrote interpolation audit")


if __name__ == "__main__":
    main()

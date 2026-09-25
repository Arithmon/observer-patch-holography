"""Reproducible conditional curvature-source -> lensed/unlensed TT, TE, EE.

Run each case in a fresh CAMB calculation, keeping its own parameters and hashes.
Does not treat agreement with the imported Planck source as an OPH prediction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "primordial"))
from sources import load_spec, source_callable  # noqa: E402


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def comparison(predicted, reference, ell_min=30, ell_max=2500):
    """Input columns ell, TT, TE, EE; TE normalization is safe at its zeros."""
    ell = reference[:, 0].astype(int)
    use = (ell >= ell_min) & (ell <= ell_max)
    chosen = ell[use]
    indices = np.searchsorted(predicted[:, 0], chosen)
    if np.any(indices >= len(predicted)) or not np.array_equal(predicted[indices, 0], chosen):
        raise ValueError("Prediction does not contain all reference multipoles")
    ref, model = reference[use, 1:4], predicted[indices, 1:4]
    if np.any(ref[:, [0, 2]] <= 0):
        raise ValueError("Reference auto-power must be positive")
    denominator = np.column_stack([ref[:, 0], np.sqrt(ref[:, 0] * ref[:, 2]), ref[:, 2]])
    residual = (model - ref) / denominator
    return {
        "ell_min": int(chosen[0]), "ell_max": int(chosen[-1]), "multipoles": len(chosen),
        "TT_fractional_rms": float(np.sqrt(np.mean(residual[:, 0]**2))),
        "TE_auto_normalized_rms": float(np.sqrt(np.mean(residual[:, 1]**2))),
        "EE_fractional_rms": float(np.sqrt(np.mean(residual[:, 2]**2))),
        "TT_fractional_max_abs": float(np.max(np.abs(residual[:, 0]))),
        "TE_auto_normalized_max_abs": float(np.max(np.abs(residual[:, 1]))),
        "EE_fractional_max_abs": float(np.max(np.abs(residual[:, 2]))),
    }


def run_case(case, output):
    import camb
    import scipy

    started = time.monotonic()
    spec = json.loads((HERE / "spec.json").read_text())
    if camb.__version__ != spec["camb_version"]:
        raise RuntimeError(f"Expected CAMB {spec['camb_version']}, found {camb.__version__}")
    ini, theory = HERE / spec["ini"], HERE / spec["theory"]
    for path, expected in [(ini, spec["ini_sha256"]), (theory, spec["theory_sha256"])]:
        if digest(path) != expected:
            raise ValueError(f"Input hash changed: {path.name}")
    if case not in spec["cases"]:
        raise ValueError(case)
    params = camb.read_ini(str(ini))
    if case != "native_planck":
        source_spec = load_spec()
        source_row = next(row for row in source_spec["candidates"] if row["id"] == case)
        effective_ns = (source_spec["baseline"]["ns"] if source_row["family"] == "powerlaw"
                        else 4 - 2 * source_row["beta"])
        spline = dict(spec["spline"])
        source = source_callable(case, numerical_floor=spline.pop("numerical_floor"))
        params.set_initial_power_function(source, effective_ns_for_nonlinear=effective_ns, **spline)
    camb.set_feedback_level(0)
    print(f"Starting {case}: CAMB {camb.__version__}, fresh transfer calculation", flush=True)
    result = camb.get_results(params)
    lmax = spec["output_lmax"]
    lensed = result.get_total_cls(lmax=lmax, CMB_unit="muK")
    unlensed = result.get_unlensed_scalar_cls(lmax=lmax, CMB_unit="muK")
    # CAMB is TT, EE, BB, TE. Export is ell, TT, TE, EE for both variants.
    table = np.column_stack([np.arange(2, lmax + 1), lensed[2:, 0], lensed[2:, 3], lensed[2:, 1],
                             unlensed[2:, 0], unlensed[2:, 3], unlensed[2:, 1]])
    if np.any(~np.isfinite(table)) or np.any(table[:, [1, 3, 4, 6]] <= 0):
        raise ValueError("Invalid auto-power or nonfinite CAMB output")
    if np.any(table[:, 2]**2 > table[:, 1] * table[:, 3] * (1 + 1e-8)):
        raise ValueError("TT/TE/EE covariance positivity failed")
    output.mkdir(parents=True, exist_ok=True)
    data_path = output / f"{case}.txt"
    np.savetxt(data_path, table, fmt=["%d"] + ["%.12e"] * 6,
               header="ell lensed_TT lensed_TE lensed_EE unlensed_TT unlensed_TE unlensed_EE; D_ell=ell(ell+1)C_ell/(2pi) in microK^2")
    param_path = output / f"{case}_params.txt"
    param_path.write_text(str(result.Params))
    reference = np.loadtxt(theory)[:, :4]
    validation = comparison(table, reference)
    gates = spec["baseline_validation"]
    baseline_pass = all(validation[key] <= gates[limit] for key, limit in [
        ("TT_fractional_rms", "TT_EE_fractional_rms_max"),
        ("EE_fractional_rms", "TT_EE_fractional_rms_max"),
        ("TE_auto_normalized_rms", "TE_rms_normalized_by_sqrt_TT_EE_max")])
    receipt = {
        "schema": "oph.conditional-boltzmann.receipt.v1", "case": case,
        "python": platform.python_version(), "camb": camb.__version__, "numpy": np.__version__,
        "scipy": scipy.__version__, "platform": platform.platform(),
        "thread_environment": {k: os.environ.get(k) for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
        "runtime_seconds": time.monotonic() - started,
        "inputs_sha256": {"run.py": digest(__file__), "spec.json": digest(HERE / "spec.json"),
            "primordial/sources.py": digest(HERE.parent / "primordial/sources.py"),
            "primordial/spec.json": digest(HERE.parent / "primordial/spec.json"),
            "planck_2018.ini": digest(ini), "official_theory": digest(theory)},
        "spectrum_sha256": digest(data_path), "parameters_sha256": digest(param_path),
        "raw_official_theory_comparison": validation,
        "baseline_engineering_gate_pass": baseline_pass if case == "native_planck" else None,
        "calPlanck": spec["calPlanck"],
        "derived_background": {key: float(value) for key, value in result.get_derived_params().items()},
        "interpretation": spec["interpretation"],
    }
    (output / f"{case}_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"case": case, "seconds": receipt["runtime_seconds"], "comparison": validation,
                      "baseline_gate": receipt["baseline_engineering_gate_pass"]}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default="native_planck")
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    cases = json.loads((HERE / "spec.json").read_text())["cases"] if args.case == "all" else [args.case]
    for case in cases:
        run_case(case, args.out)


if __name__ == "__main__":
    main()

"""Declared numerical controls after the retained UV HMcode failure."""
from __future__ import annotations

import json
import os
from pathlib import Path
import platform
import time

import camb
import numpy as np
import scipy

from run import HERE, comparison, digest, load_spec, source_callable


def calculate(case, nonlinear, floor, suffix):
    started = time.monotonic()
    spec = json.loads((HERE / "spec.json").read_text())
    control = json.loads((HERE / "control_spec.json").read_text())
    if camb.__version__ != spec["camb_version"]:
        raise ValueError("CAMB version differs from specification")
    ini, theory = HERE / spec["ini"], HERE / spec["theory"]
    if digest(ini) != spec["ini_sha256"] or digest(theory) != spec["theory_sha256"]:
        raise ValueError("Reference input hash differs")
    output = HERE / "control_results"
    output.mkdir(exist_ok=True)
    name = f"{case}_{suffix}"
    params = camb.read_ini(str(ini))
    if not nonlinear:
        params.NonLinear = camb.model.NonLinear_none
    if case != "native_planck":
        source_spec = load_spec()
        row = next(row for row in source_spec["candidates"] if row["id"] == case)
        ns = source_spec["baseline"]["ns"] if row["family"] == "powerlaw" else 4 - 2 * row["beta"]
        spline = {k: v for k, v in spec["spline"].items() if k != "numerical_floor"}
        spline["N_min"] = control["spline_N_min"]
        params.set_initial_power_function(source_callable(case, numerical_floor=floor),
                                          effective_ns_for_nonlinear=ns, **spline)
    grid = control["positivity_check"]
    ks = np.geomspace(grid["k_min"], grid["k_max"], grid["samples"])
    interpolated = params.scalar_power(ks)
    if np.any(~np.isfinite(interpolated)) or np.any(interpolated <= 0):
        raise ValueError("Actual CAMB initial-power interpolant is not positive on validation grid")
    print(f"Starting {name}", flush=True)
    receipt = {
        "schema": "oph.boltzmann-control.receipt.v1", "case": case, "name": name,
        "nonlinear_matter_lensing": nonlinear, "numerical_floor": floor,
        "camb": camb.__version__, "numpy": np.__version__, "scipy": scipy.__version__,
        "python": platform.python_version(),
        "thread_environment": {k: os.environ.get(k) for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]},
        "interpolant_positivity_probe": {"samples": len(ks), "k_min": float(ks[0]), "k_max": float(ks[-1]),
                                         "min_power": float(interpolated.min()), "negative_points": 0},
        "inputs_sha256": {path: digest(HERE / path) for path in [
            "controls.py", "control_spec.json", "run.py", "spec.json", "../primordial/sources.py", "../primordial/spec.json",
            spec["ini"], spec["theory"]]},
        "interpretation": control["scope"],
    }
    try:
        camb.set_feedback_level(0)
        results = camb.get_results(params)
        lmax = spec["output_lmax"]
        lensed = results.get_total_cls(lmax=lmax, CMB_unit="muK")
        unlensed = results.get_unlensed_scalar_cls(lmax=lmax, CMB_unit="muK")
        table = np.column_stack([np.arange(2, lmax + 1), lensed[2:, 0], lensed[2:, 3], lensed[2:, 1],
                                 unlensed[2:, 0], unlensed[2:, 3], unlensed[2:, 1]])
        if np.any(~np.isfinite(table)) or np.any(table[:, [1, 3, 4, 6]] <= 0):
            raise ValueError("Nonpositive auto-power or nonfinite output")
        target = output / f"{name}.txt"
        np.savetxt(target, table, fmt=["%d"] + ["%.12e"] * 6,
                   header="ell lensed_TT lensed_TE lensed_EE unlensed_TT unlensed_TE unlensed_EE; D_ell in microK^2")
        parameters = output / f"{name}_params.txt"
        parameters.write_text(str(results.Params))
        receipt.update({"status": "computed", "spectrum_sha256": digest(target),
                        "parameters_sha256": digest(parameters),
                        "raw_official_theory_comparison": comparison(table, np.loadtxt(theory)[:, :4])})
    except camb.CAMBError as exc:
        receipt.update({"status": "CAMB_failure", "error": str(exc)})
    receipt["runtime_seconds"] = time.monotonic() - started
    (output / f"{name}_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k in ["name", "status", "error", "runtime_seconds", "raw_official_theory_comparison"]}), flush=True)


if __name__ == "__main__":
    spec = json.loads((HERE / "spec.json").read_text())
    control = json.loads((HERE / "control_spec.json").read_text())
    for case in spec["cases"]:
        calculate(case, nonlinear=False, floor=control["primary_floor"], suffix="linear")
    calculate("history_tilt_calibrated_uv", nonlinear=False, floor=control["sensitivity_floor"], suffix="linear_floor40")
    calculate("history_tilt_calibrated_uv", nonlinear=True, floor=control["primary_floor"], suffix="hmcode_floor50")

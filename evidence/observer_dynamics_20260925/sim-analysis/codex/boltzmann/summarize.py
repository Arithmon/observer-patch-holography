"""Audit forward outputs; descriptive residuals are not a Planck likelihood."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from run import HERE, comparison, digest


def observed_diagnostics(table, calibration):
    output = {}
    for col, name in enumerate(["TT", "TE", "EE"], start=1):
        folder = "observation" if name == "TT" else "boltzmann_data"
        data = np.loadtxt(HERE.parent / folder / f"COM_PowerSpect_CMB-{name}-full_R3.01.txt")
        data = data[(data[:, 0] >= 30) & (data[:, 0] <= table[-1, 0])]
        idx = np.searchsorted(table[:, 0], data[:, 0])
        if not np.array_equal(table[idx, 0], data[:, 0]):
            raise ValueError("Data/prediction multipoles differ")
        sigma = (data[:, 2] + data[:, 3]) / 2
        residual = (data[:, 1] - table[idx, col] / calibration**2) / sigma
        output[name] = {"points": len(data), "diagonal_squared_residual_sum": float(residual @ residual),
                        "rms_in_marginal_error_units": float(np.sqrt(np.mean(residual**2)))}
    return output


def mean_cl(table, col, lower, upper):
    use = (table[:, 0] >= lower) & (table[:, 0] <= upper)
    ell = table[use, 0]
    return float(np.average(table[use, col] * 2 * np.pi / (ell * (ell + 1)), weights=2 * ell + 1))


def checked_table(directory, name):
    receipt = json.loads((directory / f"{name}_receipt.json").read_text())
    path = directory / f"{name}.txt"
    if receipt.get("status") == "CAMB_failure":
        return None, receipt
    if digest(path) != receipt["spectrum_sha256"]:
        raise ValueError(f"Changed result: {name}")
    params = directory / f"{name}_params.txt"
    if digest(params) != receipt["parameters_sha256"]:
        raise ValueError(f"Changed parameters: {name}")
    return np.loadtxt(path), receipt


def build():
    spec = json.loads((HERE / "spec.json").read_text())
    output = {"schema": "oph.boltzmann-summary.v1", "producer_sha256": digest(__file__),
              "spec_sha256": digest(HERE / "spec.json"), "cases": {}, "linear_controls": {},
              "scope": "Reference-spectrum reproduction and descriptive coadded-data residuals only. No covariance, foreground, nuisance refit, official likelihood, p value or model preference is computed."}
    native, _ = checked_table(HERE / "results", "native_planck")
    for name in spec["cases"]:
        if name == "history_tilt_calibrated_uv":
            continue  # Its first nonlinear calculation failed; retained below.
        table, receipt = checked_table(HERE / "results", name)
        output["cases"][name] = {
            "official_raw_theory": receipt["raw_official_theory_comparison"],
            "versus_native": comparison(table, native),
            "observed_full_spectrum_diagonal_diagnostic": observed_diagnostics(table, spec["calPlanck"]),
            "TT_peak_locations_in_fixed_windows": [int(table[(table[:, 0] >= lo) & (table[:, 0] <= hi)][
                np.argmax(table[(table[:, 0] >= lo) & (table[:, 0] <= hi), 1]), 0])
                for lo, hi in [(100, 350), (350, 650), (650, 950)]],
        }
    native_linear, _ = checked_table(HERE / "control_results", "native_planck_linear")
    for name in spec["cases"]:
        table, receipt = checked_table(HERE / "control_results", name + "_linear")
        output["linear_controls"][name] = {
            "versus_native_linear": comparison(table, native_linear),
            "TT_quadrupole_ratio": float(table[0, 1] / native_linear[0, 1]),
            "TT_low_ell_2_29_mean_Cl_ratio": mean_cl(table, 1, 2, 29) / mean_cl(native_linear, 1, 2, 29),
            "TT_high_ell_1800_2500_mean_Cl_ratio": mean_cl(table, 1, 1800, 2500) / mean_cl(native_linear, 1, 1800, 2500),
            "EE_high_ell_1800_2500_mean_Cl_ratio": mean_cl(table, 3, 1800, 2500) / mean_cl(native_linear, 3, 1800, 2500),
        }
    uv50, receipt50 = checked_table(HERE / "control_results", "history_tilt_calibrated_uv_linear")
    uv40, receipt40 = checked_table(HERE / "control_results", "history_tilt_calibrated_uv_linear_floor40")
    output["UV_floor_sensitivity"] = {"comparison": comparison(uv40, uv50, ell_min=2, ell_max=2508),
                                      "identical_saved_spectra": receipt40["spectrum_sha256"] == receipt50["spectrum_sha256"]}
    _, failure = checked_table(HERE / "control_results", "history_tilt_calibrated_uv_hmcode_floor50")
    output["UV_nonlinear_regularized_attempt"] = failure
    output["initial_UV_nonlinear_failure"] = {"log": "initial_attempt/run.log", "sha256": digest(HERE / "initial_attempt/run.log"),
        "error": "HMCode INTEGRATE, Integration timed out", "interpretation": "Original interpolant had tiny negative ultraviolet power; not a physical model exclusion."}
    output["observational_input_hashes"] = {str(path.relative_to(HERE.parent)): digest(path) for path in [
        HERE.parent / "observation/COM_PowerSpect_CMB-TT-full_R3.01.txt",
        HERE.parent / "boltzmann_data/COM_PowerSpect_CMB-TE-full_R3.01.txt",
        HERE.parent / "boltzmann_data/COM_PowerSpect_CMB-EE-full_R3.01.txt"]}
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    value = build()
    destination = HERE / "summary.json"
    if args.verify:
        if json.loads(destination.read_text()) != value:
            raise SystemExit("Summary differs from outputs")
        print("Verified forward summary and result hashes")
    else:
        destination.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
        print(destination)

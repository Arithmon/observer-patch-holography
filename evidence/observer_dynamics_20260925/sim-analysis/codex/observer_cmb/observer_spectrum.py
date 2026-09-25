#!/usr/bin/env python3
"""An ideal present-day observer's conditional TT/TE/EE sampling distribution.

Consumes retained, hash-pinned Boltzmann output; it never runs CAMB. The
Wishart model is exact for Gaussian T/E harmonics, but only an approximation
to a lensed sky because connected lensing covariance is omitted.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import sys

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
CODEX = HERE.parent
SPECTRUM = CODEX / "boltzmann/results/history_tilt_calibrated_wide.txt"
RECEIPT = SPECTRUM.with_name(SPECTRUM.stem + "_receipt.json")
INI = CODEX / "boltzmann_data/planck_2018.ini"
EXPECTED_SPECTRUM_SHA256 = "ade55adaa9056b6bcee2162ae2c45b37cd5db8b3a9afd0e7fba263bbda3c96ac"
EXPECTED_RECEIPT_SHA256 = "e60a982ef67a7e1a7a6a7a6fbbf4376474a985bc4570a2647ba51942c9cdbe8a"
EXPECTED_INI_SHA256 = "dc314edd1a79c2c41156159411969759ab509ddd16bda3d48c26a5fb2372fcc9"
DEFAULT_SEED = 20260925
ELL_MIN, ELL_MAX = 2, 2500
PEAK_WINDOWS = ((100, 350), (400, 650), (700, 950))
CHANNELS = ("TT", "TE", "EE")
COV_PAIRS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_spectra(ell, spectra):
    """Validate [TT, TE, EE], with consecutive integer multipoles >= 2.

    Zero auto-power and singular positive-semidefinite covariance are valid.
    A tiny round-off allowance applies to TE's correlation coefficient.
    """
    ell = np.asarray(ell, dtype=float)
    spectra = np.asarray(spectra, dtype=float)
    if ell.ndim != 1 or not len(ell) or spectra.shape != (len(ell), 3):
        raise ValueError("Expected nonempty ell and an (N, 3) TT/TE/EE array")
    if not np.isfinite(ell).all() or not np.isfinite(spectra).all():
        raise ValueError("Multipoles and spectra must be finite")
    if (ell < 2).any() or (ell != np.floor(ell)).any() or (np.diff(ell) != 1).any():
        raise ValueError("Multipoles must be consecutive increasing integers >= 2")
    tt, te, ee = spectra.T
    if (tt < 0).any() or (ee < 0).any():
        raise ValueError("Negative auto-power is not positive semidefinite")
    scale = np.sqrt(tt) * np.sqrt(ee)
    if (np.abs(te) > scale * (1 + 1e-12)).any():
        raise ValueError("TT/TE/EE covariance is not positive semidefinite")
    return ell.astype(np.int64), spectra


def d_to_c(ell, spectra):
    """D_l = l(l+1) C_l/(2 pi); input and output remain in microK^2."""
    ell, spectra = validate_spectra(ell, spectra)
    return spectra * (2 * np.pi / (ell * (ell + 1)))[:, None]


def gaussian_estimator_covariance(ell, spectra):
    """Covariance of estimated [TT, TE, EE], in square input units.

    Cov(C^XY_l,C^WZ_l)=(C^XW_l C^YZ_l+C^XZ_l C^YW_l)/(2l+1).
    The formula also applies directly to D_l when its units are used for
    every input. Different multipoles are independent in this approximation.
    """
    ell, spectra = validate_spectra(ell, spectra)
    tt, te, ee = spectra.T
    cov = np.empty((len(ell), 3, 3))
    cov[:, 0, 0] = 2 * tt**2
    cov[:, 0, 1] = cov[:, 1, 0] = 2 * tt * te
    cov[:, 0, 2] = cov[:, 2, 0] = 2 * te**2
    cov[:, 1, 1] = tt * ee + te**2
    cov[:, 1, 2] = cov[:, 2, 1] = 2 * ee * te
    cov[:, 2, 2] = 2 * ee**2
    return cov / (2 * ell + 1)[:, None, None]


def auto_sampling_interval(ell, auto_power, coverage):
    """Central marginal sampling interval for an estimator, not a posterior.

    (2l+1) * estimated_C_l / true_C_l ~ chi-square(2l+1).
    These pointwise intervals are not a simultaneous whole-spectrum band.
    """
    ell = np.asarray(ell)
    auto_power = np.asarray(auto_power, dtype=float)
    validate_spectra(ell, np.column_stack((auto_power, np.zeros_like(auto_power), auto_power)))
    if not 0 < coverage < 1:
        raise ValueError("Coverage must lie strictly between zero and one")
    alpha = (1 - coverage) / 2
    df = 2 * ell + 1
    quantiles = chi2.ppf(np.array([alpha, 1 - alpha])[:, None], df)
    return (quantiles * auto_power / df).T


def draw_wishart_spectra(ell, spectra, seed=DEFAULT_SEED):
    """One correlated ideal-observer spectrum using 2x2 Bartlett factors.

    In a real orthonormal harmonic basis there are 2l+1 independent T/E
    pairs. Their scatter matrix follows W_2(Sigma,2l+1); divide by 2l+1
    to obtain the power estimator. This shortcut avoids retaining a sky.
    """
    ell, spectra = validate_spectra(ell, spectra)
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    rng = np.random.Generator(np.random.PCG64(seed))
    df = 2 * ell + 1
    # Draw order is part of the reproducibility contract.
    a11 = np.sqrt(rng.chisquare(df))
    a21 = rng.standard_normal(len(ell))
    a22 = np.sqrt(rng.chisquare(df - 1))
    tt, te, ee = spectra.T
    l11 = np.sqrt(tt)
    l21 = np.divide(te, l11, out=np.zeros_like(te), where=l11 > 0)
    l22 = np.sqrt(np.maximum(ee - l21**2, 0))
    b11 = l11 * a11
    b21 = l21 * a11 + l22 * a21
    b22 = l22 * a22
    return np.column_stack((b11**2, b11 * b21, b21**2 + b22**2)) / df[:, None]


def load_inputs(spectrum=SPECTRUM, receipt_path=RECEIPT, ini_path=INI):
    """Reject changes to the frozen input even if its receipt is also edited."""
    hashes = {}
    for name, path, expected in (
        ("spectrum", spectrum, EXPECTED_SPECTRUM_SHA256),
        ("source_receipt", receipt_path, EXPECTED_RECEIPT_SHA256),
        ("planck_2018.ini", ini_path, EXPECTED_INI_SHA256),
    ):
        hashes[name] = sha256(path)
        if hashes[name] != expected:
            raise ValueError(f"Pinned {name} SHA-256 mismatch")
    receipt = json.loads(receipt_path.read_text())
    if receipt["spectrum_sha256"] != hashes["spectrum"]:
        raise ValueError("Source receipt spectrum SHA-256 mismatch")
    if receipt["inputs_sha256"]["planck_2018.ini"] != hashes["planck_2018.ini"]:
        raise ValueError("Source receipt parameter-file SHA-256 mismatch")
    table = np.loadtxt(spectrum)
    if table.ndim != 2 or table.shape[1] != 7:
        raise ValueError("Retained spectrum must have seven columns")
    ell, spectra = validate_spectra(table[:, 0], table[:, 1:4])
    selected = (ell >= ELL_MIN) & (ell <= ELL_MAX)
    ell, spectra = ell[selected], spectra[selected]
    if not np.array_equal(ell, np.arange(ELL_MIN, ELL_MAX + 1)):
        raise ValueError("Retained spectrum must cover ell=2..2500")
    temperature = re.findall(r"^\s*temp_cmb\s*=\s*([0-9.eE+-]+)\s*(?:#.*)?$", ini_path.read_text(), re.M)
    if len(temperature) != 1:
        raise ValueError("Expected exactly one imported temp_cmb parameter")
    return ell, spectra, receipt, float(temperature[0]), hashes


def blackbody_summary(temperature):
    """Frequency-radiance Planck-law peak, with exact SI h and k_B.

    This is an imported blackbody assumption and temperature, not a result
    derived from patch histories. It says nothing about spectral distortions.
    """
    h, k_b = 6.62607015e-34, 1.380649e-23
    x_peak = brentq(lambda x: 3 * (-np.expm1(-x)) - x, 2.0, 4.0, xtol=1e-14)
    return {
        "temperature_K": temperature,
        "temperature_status": "Imported Planck 2018 background parameter; not inferred from simulation data.",
        "assumption": "Ideal Planck blackbody; no mu/y or other spectral distortions modeled.",
        "frequency_radiance_peak_GHz": float(x_peak * k_b * temperature / h / 1e9),
        "peak_definition": "Maximum of B_nu = 2 h nu^3 / (c^2 expm1(h nu/(k_B T))); not B_lambda.",
        "h_J_s_exact_SI": h,
        "k_B_J_per_K_exact_SI": k_b,
    }


def make_output(seed=DEFAULT_SEED):
    ell, mean_d, source_receipt, temperature, hashes = load_inputs()
    mean_c = d_to_c(ell, mean_d)
    draw_d = draw_wishart_spectra(ell, mean_d, seed)
    draw_c = d_to_c(ell, draw_d)
    cov_d = gaussian_estimator_covariance(ell, mean_d)
    sigma_d = np.sqrt(np.diagonal(cov_d, axis1=1, axis2=2))
    columns = {"ell": ell.tolist()}
    for kind, array in (("mean_D", mean_d), ("mean_C", mean_c), ("draw_D", draw_d), ("draw_C", draw_c)):
        for i, channel in enumerate(CHANNELS):
            columns[f"{kind}_{channel}_uK2"] = array[:, i].tolist()
    for index, channel in ((0, "TT"), (2, "EE")):
        for coverage in (0.68, 0.95):
            bounds = auto_sampling_interval(ell, mean_d[:, index], coverage)
            for j, side in enumerate(("lower", "upper")):
                columns[f"{channel}_D_{round(100 * coverage)}_{side}_uK2"] = bounds[:, j].tolist()
    for i, channel in enumerate(CHANNELS):
        columns[f"sigma_D_{channel}_uK2"] = sigma_d[:, i].tolist()
    for i, j in COV_PAIRS:
        columns[f"cov_D_{CHANNELS[i]}_{CHANNELS[j]}_uK4"] = cov_d[:, i, j].tolist()
    peaks = []
    for lower, upper in PEAK_WINDOWS:
        indices = np.flatnonzero((ell >= lower) & (ell <= upper))
        index = indices[np.argmax(mean_d[indices, 0])]
        peaks.append({"window_ell_inclusive": [lower, upper], "ell": int(ell[index]), "D_TT_uK2": float(mean_d[index, 0])})
    expected_variance = float(np.sum((2 * ell + 1) * mean_c[:, 0]) / (4 * np.pi))
    drawn_variance = float(np.sum((2 * ell + 1) * draw_c[:, 0]) / (4 * np.pi))
    return {
        "schema": "oph.conditional-observer-cmb.v1",
        "case": "history_tilt_calibrated_wide",
        "inputs_sha256": hashes,
        "inputs_relative_to_codex": {
            "spectrum": str(SPECTRUM.relative_to(CODEX)),
            "source_receipt": str(RECEIPT.relative_to(CODEX)),
            "planck_2018.ini": str(INI.relative_to(CODEX)),
        },
        "producer": {
            "script": "observer_spectrum.py",
            "sha256": sha256(Path(__file__)),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "python": sys.version.split()[0],
            "random_generator": "numpy.random.Generator(PCG64)",
            "seed": int(seed),
            "draw_order": "Vector chi2(df); vector standard_normal; vector chi2(df-1); ell ascending.",
        },
        "observer": {
            "redshift": 0,
            "frame": "CMB rest frame",
            "sky_fraction": 1.0,
            "ell_min": ELL_MIN, "ell_max": ELL_MAX,
            "monopole_and_dipole": "Removed from angular spectra; monopole temperature reported separately as imported.",
            "instrument": "Ideal full-sky noiseless observer with no beam smoothing.",
        },
        "statistical_model": {
            "label": "Gaussian harmonic approximation using retained lensed TT/TE/EE two-point functions",
            "covariance_order": list(CHANNELS),
            "covariance_formula": "Cov(C_XY,C_WZ)=(C_XW*C_YZ+C_XZ*C_YW)/(2*ell+1)",
            "normalization": "D_ell=ell*(ell+1)*C_ell/(2*pi); both are in microK^2.",
            "sampling": "One normalized 2x2 Bartlett Wishart draw at each ell, df=2*ell+1; TT/TE/EE drawn jointly.",
            "intervals": "Pointwise central 68% and 95% chi-square sampling intervals for TT/EE estimators at fixed theoretical spectrum; not parameter posteriors or simultaneous bands.",
            "TE_sigma": "Standard deviation from Gaussian-harmonic estimator covariance; TE estimator itself is not generally Gaussian and is not chi-square.",
            "approximation": "Lensing connected covariance and correlations between multipoles are omitted. This is not an exact lensed-sky realization.",
            "excludes": ["mask", "instrument noise", "beam", "foregrounds", "peculiar velocity and kinematic dipole", "lensing connected covariance", "parameter uncertainty", "B modes", "data-derived sky phases"],
            "interpretation": "Conditional calibrated postdiction, not a new first-principles CMB derivation or a prediction of the actual phases of our sky.",
            "source_interpretation": source_receipt["interpretation"],
        },
        "summary": {
            "temperature_variance_uK2": expected_variance,
            "temperature_rms_uK": float(np.sqrt(expected_variance)),
            "draw_temperature_rms_uK": float(np.sqrt(drawn_variance)),
            "rms_definition": "sqrt(sum_{ell=2}^{2500} (2ell+1)*C_TT/(4pi)); no beam and truncated at ell=2500.",
            "fractional_temperature_rms": float(np.sqrt(expected_variance) * 1e-6 / temperature),
            "TT_peaks": peaks,
            "TT_EE_cosmic_variance_fraction": [
                {"ell": l, "sigma_over_mean": float(np.sqrt(2 / (2 * l + 1)))}
                for l in (2, 10, 30, 100, 220, 1000, 2500)
            ],
            "monopole": blackbody_summary(temperature),
        },
        "columns": columns,
    }


def table_bytes(output):
    columns = output["columns"]
    stream = io.StringIO()
    np.savetxt(stream, np.column_stack(list(columns.values())),
               fmt=["%d"] + ["%.14e"] * (len(columns) - 1),
               header=" ".join(columns))
    return stream.getvalue().encode()


def write_output(output_dir=HERE, seed=DEFAULT_SEED):
    output = make_output(seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "observer_spectrum.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    (output_dir / "observer_spectrum.txt").write_bytes(table_bytes(output))
    return output


def _compare(actual, expected, path="root"):
    """Strict structure/metadata comparison, tolerant only of float round-off."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f"{path}: object keys mismatch")
        for key in expected:
            _compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(f"{path}: list length mismatch")
        for i, (a, b) in enumerate(zip(actual, expected)):
            _compare(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (float, int)) or not np.isclose(actual, expected, rtol=2e-11, atol=1e-14):
            raise ValueError(f"{path}: numeric mismatch")
    elif actual != expected or type(actual) is not type(expected):
        raise ValueError(f"{path}: value mismatch")


def verify_output(output_dir=HERE, seed=DEFAULT_SEED):
    expected = make_output(seed)
    retained = json.loads((output_dir / "observer_spectrum.json").read_text())
    # Environment versions are provenance, not numerical acceptance criteria.
    # The producer hash, seed, input pins, assumptions and all data remain strict.
    for key in ("numpy", "scipy", "python"):
        if not isinstance(retained.get("producer", {}).get(key), str):
            raise ValueError(f"Missing producer {key} provenance")
        expected["producer"][key] = retained["producer"][key]
    _compare(retained, expected)
    text_path = output_dir / "observer_spectrum.txt"
    header = text_path.open().readline().strip()
    if header != "# " + " ".join(expected["columns"]):
        raise ValueError("Text table column header mismatch")
    actual_table = np.loadtxt(text_path)
    expected_table = np.column_stack(list(expected["columns"].values()))
    if actual_table.shape != expected_table.shape or not np.allclose(actual_table, expected_table, rtol=2e-11, atol=1e-14):
        raise ValueError("Text table numeric mismatch")
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Recompute and compare retained outputs; write nothing")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    output = verify_output(args.output_dir, args.seed) if args.verify else write_output(args.output_dir, args.seed)
    print(json.dumps({"verified" if args.verify else "written": str(args.output_dir), "summary": output["summary"]}, indent=2))


if __name__ == "__main__":
    main()

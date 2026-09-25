#!/usr/bin/env python3
"""Audit downloaded Planck PR3 plotting products without fitting a model.

Uses only Python's standard library. No network access, CAMB execution, or
modification of the observation lane. Diagonal scores are not likelihoods.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OBSERVATION = HERE.parent / "observation"
MODEL_STEM = "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum"


def digest(path):
    data = Path(path).read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def load_table(path, columns):
    lines = Path(path).read_text().splitlines()
    rows = [[float(x) for x in line.split()] for line in lines
            if line.strip() and not line.lstrip().startswith("#")]
    if not rows or any(len(row) != columns for row in rows):
        raise ValueError(f"Wrong column count in {path}")
    if not all(math.isfinite(x) for row in rows for x in row):
        raise ValueError(f"Non-finite value in {path}")
    if any(a[0] >= b[0] for a, b in zip(rows, rows[1:])):
        raise ValueError(f"Multipoles are not increasing in {path}")
    return rows


def parse_bestfit(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0].isdigit():
            values[fields[2]] = float(fields[1])
    return values


def parse_ini(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def coadded_theory(raw_theory, cal_planck):
    """Only raw CMB TT/TE/EE theory, never already-coadded data or BestFit."""
    if cal_planck <= 0:
        raise ValueError("Calibration must be positive")
    return raw_theory / cal_planck**2


def diagonal_summary(data, theory):
    """Signed data-model residuals, using the error toward the model.

    This convenient diagnostic discards covariance, frequency spectra,
    nuisance parameters, and non-Gaussian low-ell likelihood information.
    """
    if len(data) != len(theory) or not data:
        raise ValueError("Data/theory length mismatch or empty input")
    residuals = []
    for row, model in zip(data, theory):
        if min(row[2:4]) <= 0:
            raise ValueError("Uncertainty must be strictly positive")
        sigma = row[3] if model > row[1] else row[2]
        residuals.append((row[1] - model) / sigma)
    score = sum(r*r for r in residuals)
    return {
        "n_points": len(data), "ell_min": data[0][0], "ell_max": data[-1][0],
        "sum_squared_standardized_residuals": score,
        "score_per_point": score / len(data),
        "mean_signed_standardized_residual": sum(residuals) / len(data),
        "max_abs_standardized_residual": max(map(abs, residuals)),
        "is_official_likelihood": False,
        "p_value_computed": False,
    }


def verify_downloads():
    inputs = []
    for directory, manifest_name in [(HERE, "download_manifest.json"),
                                     (OBSERVATION, "download_manifest.json")]:
        manifest = json.loads((directory / manifest_name).read_text())
        for record in manifest["files"]:
            if not record["filename"].startswith("COM_PowerSpect_CMB"):
                continue
            path = directory / record["filename"]
            actual = digest(path)
            if any(actual[key] != record[key] for key in actual):
                raise ValueError(f"Download does not match receipt: {path}")
            inputs.append({"path_relative_to_codex": str(path.relative_to(HERE.parent)),
                           "url": record["url"], **actual})
    record = json.loads((HERE / "configuration_download.json").read_text())
    actual = digest(HERE / record["filename"])
    if any(actual[key] != record[key] for key in actual):
        raise ValueError("Upstream CAMB ini does not match download receipt")
    inputs.append({"path_relative_to_codex": "boltzmann_data/planck_2018.ini",
                   "url": record["url"], **actual})
    return inputs


def make_receipt():
    inputs = verify_downloads()
    best = parse_bestfit(OBSERVATION / f"{MODEL_STEM}_R3.01.txt")
    ini = parse_ini(HERE / "planck_2018.ini")
    # Compare semantically equivalent values, not unused placeholder fields.
    equivalences = {"hubble": "H0", "ombh2": "omegabh2", "omch2": "omegach2",
                    "omnuh2": "omeganuh2", "helium_fraction": "yheused",
                    "re_optical_depth": "tau", "scalar_spectral_index": "ns"}
    for key, field in equivalences.items():
        if float(ini[key]) != best[field]:
            raise ValueError(f"Best fit versus upstream ini mismatch: {key}")
    theory_path = OBSERVATION / f"{MODEL_STEM}-theory_R3.01.txt"
    if theory_path.read_text().splitlines()[0].split() != ["#", "L", "TT", "TE", "EE", "BB", "PP"]:
        raise ValueError("Unexpected raw theory column header")
    raw_theory = load_table(theory_path, 6)
    theory_by_ell = {int(row[0]): row for row in raw_theory}
    schemas, diagnostics = {}, {}
    for spectrum, column in [("TT", 1), ("TE", 2), ("EE", 3)]:
        directory = OBSERVATION if spectrum == "TT" else HERE
        full_name = f"COM_PowerSpect_CMB-{spectrum}-full_R3.01.txt"
        bin_version = "R3.01" if spectrum == "TT" else "R3.02"
        bin_name = f"COM_PowerSpect_CMB-{spectrum}-binned_{bin_version}.txt"
        full = load_table(directory / full_name, 4)
        binned = load_table(directory / bin_name, 5)
        if any(row[0] != i+2 for i, row in enumerate(full)):
            raise ValueError("Full spectrum multipoles must run consecutively from 2")
        high = [row for row in full if row[0] >= 30]
        predicted = [coadded_theory(theory_by_ell[int(row[0])][column], best["calPlanck"])
                     for row in high]
        interpolation_errors = []
        for ell, _, minus, plus, bin_model in binned:
            lo = math.floor(ell)
            frac = ell-lo
            point_model = coadded_theory(
                (1-frac)*theory_by_ell[lo][column] + frac*theory_by_ell[lo+1][column],
                best["calPlanck"])
            interpolation_errors.append((point_model-bin_model)/((minus+plus)/2))
        schemas[spectrum] = {
            "full": {"filename": full_name, "columns": ["ell", "Dl", "error_minus", "error_plus"],
                     "rows": len(full), "ell_min": full[0][0], "ell_max": full[-1][0]},
            "binned": {"filename": bin_name,
                       "columns": ["ell_effective", "Dl", "error_minus", "error_plus", "BestFit"],
                       "rows": len(binned), "ell_min": binned[0][0], "ell_max": binned[-1][0],
                       "contains_bin_edges": False, "contains_window_matrix": False},
            "units": "Dl = ell*(ell+1)*Cl/(2*pi), microkelvin squared",
            "full_high_ell_error_bars_symmetric": all(r[2] == r[3] for r in high),
            "binned_error_bars_symmetric": all(r[2] == r[3] for r in binned),
        }
        diagnostics[spectrum] = {
            "full_high_ell_against_calibrated_official_theory": diagonal_summary(high, predicted),
            "binned_against_own_official_BestFit_column": diagonal_summary(binned, [r[4] for r in binned]),
            "bin_center_interpolation_is_not_actual_binning": {
                "rms_difference_from_official_BestFit_in_sigma": math.sqrt(sum(x*x for x in interpolation_errors)/len(binned)),
                "max_abs_difference_from_official_BestFit_in_sigma": max(map(abs, interpolation_errors)),
            },
        }
    mapping = {
        "preferred_configuration": "camb.read_ini('planck_2018.ini') with runtime threads capped separately",
        "set_cosmology_equivalent_values": {
            "H0": best["H0"], "ombh2": best["omegabh2"], "omch2": best["omegach2"],
            "omk": best["omegak"], "mnu": best["mnu"], "nnu": best["nnu"],
            "num_massive_neutrinos": 1, "neutrino_hierarchy": "degenerate",
            "tau": best["tau"], "YHe": best["yheused"], "TCMB": float(ini["temp_cmb"]),
        },
        "primordial": {"As_upstream_ini": float(ini["scalar_amp"]),
                       "As_from_printed_logA": math.exp(best["logA"])*1e-10,
                       "As_from_printed_derived_A": best["A"]*1e-9,
                       "ns": best["ns"], "nrun": 0, "nrunrun": 0, "r": 0,
                       "pivot_scalar_Mpc_inverse": float(ini["pivot_scalar"])},
        "neutrino_density_target": best["omeganuh2"],
        "neutrino_ini": {key: ini[key] for key in ["massless_neutrinos", "nu_mass_eigenstates", "massive_neutrinos", "share_delta_neff", "nu_mass_fractions"]},
        "helium": {"mass_fraction_used": best["yheused"], "unused_yhe_placeholder": best["yhe"],
                   "BBN_nucleon_fraction_not_CAMB_mass_fraction": best["YpBBN"]},
        "reionization_width": float(ini["re_delta_redshift"]),
        "Alens": best["Alens"], "w": best["w"], "wa": best["wa"],
        "nonlinear_model": "HMcode2016; halofit_version=5 or 'mead2016'/'mead'",
        "nonlinear_ini_mode": int(ini["do_nonlinear"]),
        "recombination": {key: ini[key] for key in ["recombination_model", "RECFAST_H_fudge", "RECFAST_fudge_He", "RECFAST_Heswitch", "RECFAST_Hswitch", "RECFAST_He_rate_correction"]},
        "calPlanck": best["calPlanck"],
        "raw_CMB_theory_to_coadded_factor": best["calPlanck"]**-2,
        "caveats": [
            "The upstream ini sets initial_ratio=1 but get_tensor_cls=F; it does not imply r=1.",
            "The physical neutrino density target and effective species are more reproducible than trusting version-dependent mass conversion defaults.",
            "Raw-theory validation uses no calPlanck rescaling; coadded TT/TE/EE comparison divides raw theory once by calPlanck squared.",
            "Do not rescale observation data, the binned BestFit column, or the lensing PP column using this CMB calibration conversion.",
            "The ini pins physical inputs, not the CAMB executable. Current CAMB need not reproduce the old release bit for bit.",
        ],
    }
    return {"schema": "oph.codex.planck-pr3-input-audit.v1", "inputs": inputs,
            "parameter_mapping": mapping, "spectra_schemas": schemas,
            "theory_schema": {"columns": ["L", "TT", "TE", "EE", "BB", "PP"],
                              "rows": len(raw_theory), "ell_min": raw_theory[0][0], "ell_max": raw_theory[-1][0],
                              "CAMB_cmb_array_columns": ["TT", "EE", "BB", "TE"],
                              "PP_used_in_this_audit": False},
            "diagnostics": diagnostics, "candidate_model_fitted": False,
            "interpretation": "Descriptive diagonal residuals only; missing covariance and nuisance likelihood prohibit official chi-square significance or model evidence."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Recompute and compare existing receipt without writing")
    args = parser.parse_args()
    receipt = make_receipt()
    output = HERE / "input_audit.json"
    if args.check:
        if json.loads(output.read_text()) != receipt:
            raise SystemExit("Input audit receipt mismatch")
        print("Input hashes, schemas, parameter mapping, and diagnostic receipt verified")
    else:
        output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
        print(output)


if __name__ == "__main__":
    main()

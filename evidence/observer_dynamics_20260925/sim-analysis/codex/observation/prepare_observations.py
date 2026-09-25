"""Verify small official source files and expose explicit observation units.

Run with no arguments to regenerate JSON in this directory, or --verify to
compare existing JSON to the pinned raw files. Does not fit any model.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TT = "COM_PowerSpect_CMB-TT-full_R3.01.txt"
THEORY = "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt"
PARAMS = "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum_R3.01.txt"
FIRAS = "firas_monopole_spec_v1.txt"


def table(path):
    return [[float(v) for v in line.split()] for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]


def verify_inputs(root=ROOT):
    manifest = json.loads((root / "download_manifest.json").read_text())
    for row in manifest["files"]:
        data = (root / row["filename"]).read_bytes()
        if len(data) != row["bytes"] or hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"source mismatch: {row['filename']}")
    return manifest


def prepare(root=ROOT):
    manifest = verify_inputs(root)
    pins = {row["filename"]: row for row in manifest["files"]}
    theory = {int(row[0]): row[1] for row in table(root / THEORY)}
    calibration = None
    for line in (root / PARAMS).read_text().splitlines():
        cols = line.split()
        if len(cols) > 2 and cols[2] == "calPlanck":
            calibration = float(cols[1])
    if calibration is None:
        raise ValueError("missing best-fit calibration")
    full = table(root / TT)
    if [int(r[0]) for r in full] != list(range(2, 2509)):
        raise ValueError("unexpected Planck multipole coverage")
    rows = []
    for ell, dl, minus, plus in full:
        ell = int(ell)
        if ell > 40:
            continue
        factor = 2 * math.pi / (ell * (ell + 1))
        rows.append({
            "ell": ell,
            "Dl_uK2": dl, "Dl_error_minus_uK2": minus, "Dl_error_plus_uK2": plus,
            "Cl_uK2": dl * factor, "Cl_error_minus_uK2": minus * factor,
            "Cl_error_plus_uK2": plus * factor,
            "component": "Commander_lowell" if ell <= 29 else "Plik_highell",
            "primary_lowell": ell <= 29,
            "lcdm_bestfit_Dl_uK2_raw": theory[ell],
            "lcdm_bestfit_Dl_uK2_calibrated": theory[ell] / calibration**2,
        })
    planck = {
        "schema": "oph.observation.planck-pr3-tt-small.v1",
        "rows": rows,
        "definition": "Dl = ell*(ell+1)*Cl/(2*pi), in thermodynamic microkelvin squared",
        "bestfit_calPlanck": calibration,
        "source_files": [pins[name] for name in (TT, THEORY, PARAMS)],
        "source_description": "https://esdcdoi.esac.esa.int/doi/html/data/astronomy/planck/Cosmology.html",
        "release_doi": "https://doi.org/10.5270/esa-gb3sw1a",
        "likelihood_reference": "https://www.aanda.org/articles/aa/full_html/2020/09/aa36386-19/aa36386-19.html",
        "interpretation": {
            "ell_2_29": "Commander component separation, 86% of sky; asymmetric marginal 68% confidence limits include foreground-subtraction uncertainty",
            "ell_30_40": "Plik cross-half-mission; foreground and nuisance parameters fixed to base-LCDM best fit, symmetric 1-sigma errors include beam uncertainty",
            "theory": "baseline TTTEEE+lowl+lowE+lensing best-fit theory; calibrated column divides original TT by calPlanck squared, as ESA documents",
        },
        "limitations": [
            "These are published spectrum estimates and marginal intervals, not a likelihood or its full covariance.",
            "Ideal full-sky chi-square sampling or f_sky=0.86 rescaling is a diagnostic approximation, not the Planck likelihood.",
            "Do not combine posterior error bars and an additional cosmic-variance term as if independent errors.",
            "Do not report a primordial scalar index from a directly fitted two-dimensional temperature spectrum.",
            "Fitting only an amplitude can reject a chosen field-to-temperature bridge; it cannot reject all OPH models.",
        ],
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    firas = {
        "schema": "oph.observation.cobe-firas-monopole-small.v1",
        "rows": [{"wavenumber_per_cm": k, "frequency_GHz": k * 29.9792458,
                  "monopole_MJy_per_sr": intensity, "residual_kJy_per_sr": residual,
                  "sigma_kJy_per_sr": sigma, "galaxy_model_kJy_per_sr": galaxy}
                 for k, intensity, residual, sigma, galaxy in table(root / FIRAS)],
        "source_file": pins[FIRAS],
        "source_description": "https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html",
        "original_paper": "https://arxiv.org/abs/astro-ph/9605054",
        "published_95_percent_limits": {"absolute_mu": 9e-5, "absolute_y": 1.5e-5},
        "interpretation": "NASA file adds Table4 residuals to a 2.725K blackbody; residual/sigma/galaxy columns are kJy/sr, while total intensity is MJy/sr",
        "limitations": ["Table contains marginal uncertainties, not the complete correlated spectral covariance.",
                        "Published mu/y limits are cited from the primary analysis, not refitted from diagonal table errors.",
                        "Repair energy is not photon heat without an extra physical conversion and redshift assignment."],
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    return {"planck_tt_l2_40.json": planck, "firas_monopole.json": firas}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    for filename, data in prepare().items():
        path = ROOT / filename
        if args.verify:
            if json.loads(path.read_text()) != data:
                raise SystemExit(f"derived product mismatch: {filename}")
        else:
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        print(("verified " if args.verify else "wrote ") + str(path))


if __name__ == "__main__":
    main()

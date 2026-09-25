"""Compare the retained conditional observer spectrum with Planck PR3 data.

Binned points are visual overlays only. Residuals use exact full multipoles,
the official coadded calibration, and the quoted marginal error toward the
model. No parameters are fitted and no likelihood significance is computed.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import observer_spectrum as observer

HERE = Path(__file__).resolve().parent
CODEX = HERE.parent
MANIFEST_PINS = {
    "observation": "45e5b586ed61745b66036af37fd2701c0b09786d2a76271c013653ea936b3939",
    "boltzmann_data": "796bc95515aabc41dab6c0273349d06d1349de7ca8743a3a6049d4fcea2285ba",
}
BACKGROUND, INK, MUTED = "#131922", "#edf2f7", "#b3bece"
MODEL, DATA = "#68dbd1", "#efa381"
STYLE = {"figure.facecolor": BACKGROUND, "axes.facecolor": BACKGROUND,
         "savefig.facecolor": BACKGROUND, "axes.edgecolor": "#667384",
         "text.color": INK, "axes.labelcolor": MUTED, "xtick.color": MUTED,
         "ytick.color": MUTED, "font.size": 11, "axes.titlesize": 14}


def load_measurements():
    manifests = {}
    for folder, digest in MANIFEST_PINS.items():
        path = CODEX / folder / "download_manifest.json"
        if observer.sha256(path) != digest:
            raise ValueError(f"Changed Planck download manifest: {folder}")
        manifests[folder] = json.loads(path.read_text())
    tables, inputs = {}, []
    for channel in observer.CHANNELS:
        folder = "observation" if channel == "TT" else "boltzmann_data"
        for kind in ("full", "binned"):
            version = "3.02" if kind == "binned" and channel != "TT" else "3.01"
            name = f"COM_PowerSpect_CMB-{channel}-{kind}_R{version}.txt"
            row = next(r for r in manifests[folder]["files"] if r["filename"] == name)
            path = CODEX / folder / name
            if observer.sha256(path) != row["sha256"] or path.stat().st_size != row["bytes"]:
                raise ValueError(f"Changed measured spectrum: {name}")
            data = np.loadtxt(path)
            expected_cols = 4 if kind == "full" else 5
            if data.ndim != 2 or data.shape[1] != expected_cols or not np.isfinite(data).all():
                raise ValueError(f"Invalid measured table: {name}")
            if (np.diff(data[:, 0]) <= 0).any() or (data[:, 2:4] <= 0).any():
                raise ValueError(f"Invalid multipoles or error magnitudes: {name}")
            tables[channel, kind] = data
            inputs.append({"path": str(path.relative_to(CODEX)), "url": row["url"],
                           "sha256": row["sha256"], "bytes": row["bytes"]})
    return tables, inputs


def exact_residuals(data, ell, raw_model, calibration, ell_min=30):
    """Return ell, observation, coadded model, error minus/plus, residual/sigma."""
    if not np.isfinite(calibration) or calibration <= 0:
        raise ValueError("Calibration must be finite and positive")
    use = (data[:, 0] >= ell_min) & (data[:, 0] <= ell[-1])
    data = data[use]
    indices = np.searchsorted(ell, data[:, 0])
    if not len(data) or np.any(indices >= len(ell)) or not np.array_equal(ell[indices], data[:, 0]):
        raise ValueError("Residuals require exact matching multipoles; no bin-center interpolation")
    model = raw_model[indices] / calibration**2
    sigma = np.where(model > data[:, 1], data[:, 3], data[:, 2])
    if not np.isfinite(sigma).all() or np.any(sigma <= 0):
        raise ValueError("Quoted uncertainties must be finite and positive")
    return np.column_stack((data[:, 0], data[:, 1], model, data[:, 2], data[:, 3],
                            (data[:, 1] - model) / sigma))


def prepare():
    source = observer.verify_output()
    columns = {k: np.asarray(v) for k, v in source["columns"].items()}
    tables, inputs = load_measurements()
    cal = json.loads(observer.RECEIPT.read_text())["calPlanck"]
    residuals = {
        ch: exact_residuals(tables[ch, "full"], columns["ell"], columns[f"mean_D_{ch}_uK2"], cal)
        for ch in observer.CHANNELS
    }
    diagnostics = {}
    for ch, rows in residuals.items():
        z = rows[:, -1]
        diagnostics[ch] = {"multipoles": len(rows), "ell_min": int(rows[0, 0]),
                           "ell_max": int(rows[-1, 0]), "mean_standardized_residual": float(z.mean()),
                           "rms_in_quoted_error_units": float(np.sqrt(np.mean(z*z))),
                           "diagonal_squared_residual_sum": float(z @ z)}
    receipt = {
        "schema": "oph.conditional-observer-planck-comparison.v1",
        "producer_sha256": observer.sha256(Path(__file__)),
        "observer_receipt_sha256": observer.sha256(HERE / "observer_spectrum.json"),
        "source_spectrum_sha256": observer.EXPECTED_SPECTRUM_SHA256,
        "measurement_inputs": inputs, "download_manifest_pins": MANIFEST_PINS,
        "calPlanck": cal, "raw_theory_to_coadded_factor": 1/cal**2,
        "calibration_rule": "Divide raw model and ideal-sky bands by calPlanck^2 once; never rescale measurements.",
        "residual_definition": "(measurement - coadded conditional model) / quoted error toward model, exact full multipoles only.",
        "diagnostics": diagnostics,
        "scope": [
            "Mean from calibrated history source plus standard transfer; amplitude, tilt and background imported.",
            "Official binned plotting data and asymmetric marginal error bars overlaid with a continuous model.",
            "Binned BestFit columns are not used as measurements. No bin-center residuals are calculated.",
            "Full-multipole residuals at ell>=30 use quoted marginal errors; they omit covariance and nuisance uncertainty.",
            "Diagnostic sums are not an official likelihood, p value, independent degrees of freedom or model evidence.",
            "Large-angle TT bands are ideal full-sky Gaussian sampling bands, not Planck error bars; lensed connected covariance omitted.",
            "No fitting, new cosmological evolution or new simulated sky realization is performed.",
        ],
    }
    return columns, tables, cal, residuals, receipt


def table_bytes(rows):
    stream = io.StringIO()
    np.savetxt(stream, rows, fmt=["%d"] + ["%.14e"]*5,
               header="ell observed_D_uK2 conditional_model_D_coadded_uK2 error_minus_uK2 error_plus_uK2 residual_over_quoted_error")
    return stream.getvalue().encode()


def finish_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=.12)


def save(fig, name):
    for suffix in ("png", "svg"):
        fig.savefig(HERE / f"{name}.{suffix}", dpi=175)
    plt.close(fig)


def plot_comparison(columns, tables, cal, residuals):
    ell = columns["ell"]
    titles = {"TT": "Temperature · TT", "TE": "Temperature × polarization · TE", "EE": "E-mode polarization · EE"}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(3, 2, figsize=(15, 11.5), gridspec_kw={"width_ratios": [1.65, 1]})
        for (ch, pair) in zip(observer.CHANNELS, axes):
            left, right = pair
            binned = tables[ch, "binned"]
            left.plot(ell, columns[f"mean_D_{ch}_uK2"]/cal**2, color=MODEL, lw=1.8,
                      label="Conditional observer spectrum")
            left.errorbar(binned[:, 0], binned[:, 1], yerr=binned[:, 2:4].T,
                          fmt="o", ms=3, color=DATA, lw=.75, capsize=0, label="Planck PR3 measured bandpowers")
            left.set(title=titles[ch], xlim=(30, 2500), ylabel=rf"$D_\ell^{{{ch}}}$ [$\mu$K$^2$]")
            if ch == "TT":
                left.legend(loc="upper right", frameon=False, fontsize=10)
            rows = residuals[ch]
            right.axhspan(-1, 1, color=MUTED, alpha=.09)
            right.axhline(0, color=MODEL, lw=1)
            right.scatter(rows[:, 0], rows[:, -1], s=4, color=DATA, alpha=.43, linewidths=0, rasterized=True)
            limit = max(4., float(np.max(np.abs(rows[:, -1]))) + .25)
            right.set(xlim=(30, rows[-1, 0]), ylim=(-limit, limit), ylabel="(Measured − model) / quoted error")
            right.set_title(f"Individual multipoles · {len(rows):,} measurements", loc="left", fontsize=12)
            for ax in pair:
                ax.set_xlabel(r"Multipole $\ell$")
                finish_axis(ax)
                ax.title.set_ha("left")
                ax.title.set_position((0, 1.0))
        fig.suptitle("Conditional observer spectrum vs measured CMB", x=.075, ha="left", y=.975, fontsize=24)
        fig.text(.075, .931, "Temperature and polarization from Planck PR3  |  Curves use a calibrated history source and standard cosmological transfer", color=MUTED, fontsize=12)
        fig.text(.075, .038, "Left: published binned measurements with error bars; continuous model. Right: residuals at exact full multipoles, not bin centers.\n"
                 "Amplitude, tilt and background are supplied. Residuals use marginal errors only; this comparison is not a Planck likelihood or independent OPH test.",
                 color=MUTED, fontsize=10, linespacing=1.7)
        fig.subplots_adjust(left=.075, right=.97, top=.87, bottom=.13, hspace=.56, wspace=.26)
        save(fig, "measured_comparison")


def plot_large_angles(columns, tables, cal):
    select = columns["ell"] <= 29
    ell = columns["ell"][select]
    data = tables["TT", "full"]
    data = data[data[:, 0] <= 29]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(12, 6.8))
        for pct, alpha in ((95, .12), (68, .25)):
            ax.fill_between(ell, columns[f"TT_D_{pct}_lower_uK2"][select]/cal**2,
                            columns[f"TT_D_{pct}_upper_uK2"][select]/cal**2,
                            color=MODEL, alpha=alpha, label=f"{pct}% ideal-sky pointwise sampling band")
        ax.plot(ell, columns["mean_D_TT_uK2"][select]/cal**2, lw=2, color=MODEL,
                label="Conditional ensemble mean")
        ax.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2:4].T, color=DATA,
                    fmt="o", ms=4, lw=1, label="Planck measured TT + quoted errors")
        ax.set(xlim=(1.5, 29.5), ylim=(0, None), xlabel=r"Multipole $\ell$", ylabel=r"$D_\ell^{TT}$ [$\mu$K$^2$]")
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, frameon=False, loc="upper left", ncol=2,
                   fontsize=10, bbox_to_anchor=(.08, .847))
        finish_axis(ax)
        fig.suptitle("The largest angular scales", x=.09, ha="left", y=.96, fontsize=23)
        fig.text(.09, .893, "Actual Planck measurements against the conditional mean and ideal full-sky sampling range", color=MUTED, fontsize=11)
        fig.text(.09, .052, "Shaded bands describe hypothetical Gaussian full-sky draws; error bars belong to the measured Planck estimates.\n"
                 "They are different uncertainty objects. Lensing covariance is omitted; no low-multipole likelihood is evaluated.",
                 color=MUTED, fontsize=10, linespacing=1.6)
        fig.subplots_adjust(left=.09, right=.97, top=.74, bottom=.20)
        save(fig, "measured_large_angles")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Recompute numerical comparisons and verify all retained output hashes")
    args = parser.parse_args()
    columns, tables, cal, residuals, receipt = prepare()
    target = HERE / "measured_comparison.json"
    if args.verify:
        retained = json.loads(target.read_text())
        outputs = retained.pop("outputs_sha256")
        observer._compare(retained, receipt)
        for ch, rows in residuals.items():
            if (HERE / f"measured_residuals_{ch}.txt").read_bytes() != table_bytes(rows):
                raise ValueError(f"Changed residual table: {ch}")
        expected_outputs = {f"measured_residuals_{ch}.txt" for ch in observer.CHANNELS}
        expected_outputs |= {f"{name}.{ext}" for name in ("measured_comparison", "measured_large_angles") for ext in ("png", "svg")}
        if set(outputs) != expected_outputs:
            raise ValueError("Output inventory mismatch")
        for name, digest in outputs.items():
            if observer.sha256(HERE / name) != digest:
                raise ValueError(f"Output hash mismatch: {name}")
    else:
        for ch, rows in residuals.items():
            (HERE / f"measured_residuals_{ch}.txt").write_bytes(table_bytes(rows))
        plot_comparison(columns, tables, cal, residuals)
        plot_large_angles(columns, tables, cal)
        names = [f"measured_residuals_{ch}.txt" for ch in observer.CHANNELS]
        names += [f"{name}.{ext}" for name in ("measured_comparison", "measured_large_angles") for ext in ("png", "svg")]
        receipt["outputs_sha256"] = {name: observer.sha256(HERE / name) for name in names}
        target.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"verified" if args.verify else "written": str(target), "diagnostics": receipt["diagnostics"]}, indent=2))


if __name__ == "__main__":
    main()

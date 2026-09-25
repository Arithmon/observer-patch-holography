"""Plot retained conditional observer spectra; never fit or evolve a source."""
from pathlib import Path
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from observer_spectrum import verify_output


def main():
    here = Path(__file__).resolve().parent
    result = verify_output()
    c = {key: np.asarray(value) for key, value in result["columns"].items()}
    ell = c["ell"]
    background, ink, muted = "#131922", "#edf2f7", "#b3bece"
    cyan, coral = "#68dbd1", "#efa381"
    with plt.rc_context({
        "figure.facecolor": background, "axes.facecolor": background,
        "axes.edgecolor": "#667384", "text.color": ink,
        "axes.labelcolor": muted, "xtick.color": muted, "ytick.color": muted,
        "font.size": 11, "axes.titlesize": 13, "legend.fontsize": 9,
        "savefig.facecolor": background,
    }):
        fig = plt.figure(figsize=(13, 9))
        grid = fig.add_gridspec(2, 2, width_ratios=[1.65, 1], hspace=.42, wspace=.27)
        tt, low, te, ee = [fig.add_subplot(grid[i, j]) for i, j in ((0, 0), (0, 1), (1, 0), (1, 1))]
        for ax, channel, small in ((tt, "TT", False), (low, "TT", True), (ee, "EE", False)):
            sel = (ell <= 30) if small else (ell >= 30)
            x = ell[sel]
            ax.fill_between(x, c[f"{channel}_D_95_lower_uK2"][sel], c[f"{channel}_D_95_upper_uK2"][sel],
                            color=cyan, alpha=.10, label="95% pointwise sampling band")
            ax.fill_between(x, c[f"{channel}_D_68_lower_uK2"][sel], c[f"{channel}_D_68_upper_uK2"][sel],
                            color=cyan, alpha=.22, label="68% pointwise sampling band")
            ax.plot(x, c[f"draw_D_{channel}_uK2"][sel], color=coral,
                    lw=.6, alpha=.6, marker="." if small else None, ms=4,
                    label="One correlated spectrum draw")
            ax.plot(x, c[f"mean_D_{channel}_uK2"][sel], color=cyan, lw=1.7,
                    label="Conditional ensemble mean")
            ax.set_xlim(x[0], x[-1])
            ax.set_ylim(bottom=0)
        te_sel = ell >= 30
        x = ell[te_sel]
        mean = c["mean_D_TE_uK2"][te_sel]
        sigma = c["sigma_D_TE_uK2"][te_sel]
        te.fill_between(x, mean-sigma, mean+sigma, color=cyan, alpha=.22)
        te.plot(x, c["draw_D_TE_uK2"][te_sel], color=coral, lw=.6, alpha=.5)
        te.plot(x, mean, color=cyan, lw=1.7)
        te.axhline(0, lw=.6, color=muted, alpha=.5)
        te.set_xlim(30, 2500)
        for ax, title in ((tt, "Temperature · TT"), (low, "Largest angular scales · TT"),
                          (te, "Temperature × polarization · TE"), (ee, "E-mode polarization · EE")):
            ax.set_title(title, loc="left", pad=12)
            ax.set_xlabel(r"Multipole $\ell$")
            ax.set_ylabel(r"$D_\ell$ [$\mu$K$^2$]")
            ax.spines[["right", "top"]].set_visible(False)
            ax.grid(alpha=.10)
        low.set_xticks([2, 10, 20, 30])
        low.legend(frameon=False, loc="upper right", fontsize=8)
        for peak in result["summary"]["TT_peaks"]:
            tt.annotate(str(peak["ell"]), (peak["ell"], peak["D_TT_uK2"]),
                        xytext=(0, 18), textcoords="offset points", ha="center", fontsize=9,
                        arrowprops={"arrowstyle": "-", "color": muted, "lw": .6})
        fig.suptitle("The CMB for an ideal observer today", x=.08, ha="left", y=.965, fontsize=23)
        fig.text(.08, .915, "Calibrated history source + standard cosmological transfer  |  CMB rest frame, redshift zero", color=muted)
        fig.text(.08, .06, "Lensed mean spectra; Gaussian full-sky sampling approximation. Lensing covariance, noise and foregrounds omitted.\n"
                 "Amplitude, tilt and background are imported. The draw is illustrative; it does not reconstruct our sky. TE shading is ±1 SD.",
                 color=muted, fontsize=10, linespacing=1.6)
        fig.subplots_adjust(left=.08, right=.97, bottom=.16, top=.84)
        output = here / "observer_spectrum.png"
        fig.savefig(output, dpi=170)
        plt.close(fig)
        print(json.dumps({"figure": str(output), "temperature_rms_uK": result["summary"]["temperature_rms_uK"]}))


if __name__ == "__main__":
    main()

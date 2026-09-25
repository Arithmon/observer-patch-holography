"""Export standalone CMB comparison figures from verified small results."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "primordial"))
from sources import source_callable

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.facecolor": "white"})


def save(fig, name):
    fig.savefig(HERE / f"{name}.png", dpi=180)
    fig.savefig(HERE / f"{name}.pdf")
    plt.close(fig)


def main():
    native = np.loadtxt(HERE / "results/native_planck.txt")
    history = np.loadtxt(HERE / "results/history_tilt_calibrated_wide.txt")
    cal = 1.000442**2
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for col, (name, ax) in enumerate(zip(["TT", "TE", "EE"], axes), 1):
        folder = "observation" if name == "TT" else "boltzmann_data"
        version = "3.01" if name == "TT" else "3.02"
        data = np.loadtxt(HERE.parent / folder / f"COM_PowerSpect_CMB-{name}-binned_R{version}.txt")
        ax.plot(history[:, 0], history[:, col] / cal, color="#bd3654", lw=2.2, label="History source + standard transfer (calibrated)")
        ax.plot(native[:, 0], native[:, col] / cal, color="#216a89", lw=1.3, ls=(0, (5, 4)), label="Standard Planck-source calculation")
        ax.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2:4].T, fmt=".", color="#182c3d", ms=4,
                    lw=.6, capsize=0, alpha=.75, label="Planck PR3 measurements")
        ax.set_ylabel(rf"$D_\ell^{{{name}}}$ [$\mu$K$^2$]")
        ax.grid(alpha=.15)
        ax.set_xlim(2, 2508)
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    axes[2].set_xlabel(r"Angular multipole $\ell$")
    fig.suptitle("A primordial history source can reproduce the acoustic spectra", fontsize=16, x=.08, ha="left", y=.975)
    fig.text(.08, .941, "Conditional construction: Planck-calibrated amplitude, tilt and background", fontsize=11, color="#506375")
    fig.text(.08, .025, "TT: temperature • TE: temperature–polarization correlation • EE: polarization\n"
             "Curves are continuous predictions; points are published plotting products. No Planck likelihood fit is performed.", fontsize=9, color="#506375")
    fig.subplots_adjust(left=.10, right=.97, top=.90, bottom=.105, hspace=.15)
    save(fig, "cmb_forward")

    rows = [
        ("history_scale_invariant_wide", "Untilted histories", "#b18520"),
        ("history_tilt_calibrated_wide", "Calibrated histories", "#bd3654"),
        ("history_tilt_calibrated_ir", "Finite oldest history (IR)", "#238879"),
        ("history_tilt_calibrated_uv", "Finite youngest history (UV)", "#7750a5"),
    ]
    base = np.loadtxt(HERE / "control_results/native_planck_linear.txt")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    k = np.geomspace(1e-6, .5, 1000)
    for key, label, color in rows:
        axes[0, 0].plot(k, source_callable(key)(k) / 2.100549e-9, color=color, label=label)
        table = np.loadtxt(HERE / f"control_results/{key}_linear.txt")
        use = table[:, 0] <= 60
        axes[0, 1].plot(table[use, 0], table[use, 1] / base[use, 1], color=color)
        for ax, col in [(axes[1, 0], 1), (axes[1, 1], 3)]:
            ax.plot(table[:, 0], 100 * (table[:, col] / base[:, col] - 1), color=color)
    axes[0, 0].plot(k, source_callable("planck_powerlaw")(k) / 2.100549e-9, color="black", ls="--", lw=1, label="Standard power law")
    axes[0, 0].set(xscale="log", yscale="log", ylim=(.01, 2), xlabel=r"Comoving $k$ [Mpc$^{-1}$]", ylabel=r"Primordial $\Delta_{\mathcal{R}}^2 / A_s$")
    axes[0, 1].set(xlabel=r"Multipole $\ell$", ylabel="Low-multipole TT / baseline", xlim=(2, 60), ylim=(.5, 1.1))
    for ax, label in zip(axes[1], ["TT", "EE"]):
        ax.set(xlabel=r"Multipole $\ell$", ylabel=f"{label} change from baseline [%]", xlim=(30, 2500))
        ax.axhline(0, color="black", ls="--", lw=.7, alpha=.6)
    for ax in axes.ravel():
        ax.grid(alpha=.15)
    fig.suptitle("History assumptions have distinct, testable spectral effects", x=.08, ha="left", y=.98, fontsize=16)
    fig.text(.08, .941, "Fixed source menu; common standard background and linear-matter lensing in every curve", fontsize=11, color="#506375")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(.5, .022))
    fig.subplots_adjust(left=.08, right=.98, bottom=.17, top=.88, hspace=.3, wspace=.28)
    save(fig, "source_sensitivity")


if __name__ == "__main__":
    main()

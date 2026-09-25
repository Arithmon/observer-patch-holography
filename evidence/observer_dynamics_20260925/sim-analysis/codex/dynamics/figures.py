"""Standalone figures for the frozen finite dynamics experiment."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
r = json.loads((HERE / "receipt.json").read_text())
levels = r["levels"]
fig, ax = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
ax[0, 0].plot([v["carriers"] for v in levels], [v["gap_times_carriers"] for v in levels], "o-")
ax[0, 0].set(xscale="log", xlabel="Carriers N", ylabel="N × first positive eigenvalue", title="Resolved finite diffusion gap")
for v in levels[2:]:
    bands = v["sphere_band_diagnostics"]
    ax[0, 1].plot([b["sphere_index_band"] for b in bands], [b["mean_times_carriers_over_ell_ellplus1"] for b in bands], "o-", label=f"L{v['level']}")
ax[0, 1].set(xlabel="Sphere-index band ℓ (diagnostic)", ylabel="N × mean λ / [ℓ(ℓ+1)]", title="Low bands admit a common diffusion scale")
ax[0, 1].legend()
for control in r["controls"]:
    rows = control["heat"]
    ax[1, 0].plot([p["sweeps"] for p in rows], [p["dimension_lower"] for p in rows], "o--", alpha=.6, label=f"Torus d={control['known_dimension']}")
rows = [p for p in levels[2]["heat"] if p["truncation_accepted"] and .03 <= p["gap_scaled_time"] <= .3]
second = ax[1, 0].twiny()
second.plot([p["sweeps"] for p in rows], [(p["dimension_lower"] + p["dimension_upper"]) / 2 for p in rows], "o-", color="black", label="OPH screen L2")
second.set_xlabel("OPH L2 sweep time (upper axis)")
second.legend(loc="center right")
ax[1, 0].set(xlabel="Torus sweep time (lower axis)", ylabel="Heat spectral dimension", title="Screen transport versus controls", ylim=(0, 3.5))
ax[1, 0].legend(loc="lower left", fontsize=8)
h = levels[-1]["histories"]
ax[1, 1].plot(h["dimensionless_times"], h["diffusion"], label="Declared mean repair")
ax[1, 1].plot(h["dimensionless_times"], h["monotone_log_clock"], label="Monotone clock change")
ax[1, 1].plot(h["dimensionless_times"], h["wave_extension_zero_initial_velocity"], "--", label="Added inertial wave control")
ax[1, 1].set(xlabel="Dimensionless mode time", ylabel="Mode amplitude", title="A clock change does not supply inertia")
ax[1, 1].legend(fontsize=8)
fig.suptitle("Finite seam-graph dynamics: diffusion, hidden-port correction and clock boundary", fontsize=13)
for suffix in ("png", "pdf"):
    fig.savefig(HERE / f"dynamics.{suffix}", dpi=160)

geometry = json.loads((HERE / "geometry_receipt.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
xs = [v["level"] for v in levels]
axes[0].semilogy(xs, [v["sphere_band_diagnostics"][2]["splitting_fraction"] for v in levels], "o-", label="Actual twelve-port seam graph")
axes[0].semilogy(xs, [v["cell_dual"]["bands"][2]["splitting_fraction"] for v in geometry["levels"]], "s-", label="Unweighted cell-dual graph")
fem = geometry["levels"][1:]
axes[0].semilogy(xs[1:], [v["round_embedding_fem"]["bands"][2]["splitting_fraction"] for v in fem], "^-", label="Added round-metric FEM operator")
axes[0].set(xlabel="Tower level", ylabel="ℓ=3 index-band range / mean", title="The metric operator changes the symmetry trend")
axes[0].legend(fontsize=8)
last = geometry["levels"][-1]
for entry, label, color in ((levels[-1]["low_eigenvalues"][9:16], "Actual ports", "C0"),
                            (last["cell_dual"]["values"][9:16], "Cell dual", "C1"),
                            (last["round_embedding_fem"]["values"][9:16], "Round FEM", "C2")):
    normalized = np.asarray(entry) / np.mean(entry)
    axes[1].plot(range(1, 8), normalized, "o-", label=label, color=color)
axes[1].axhline(1., color="gray", linewidth=.7)
axes[1].set(xlabel="Index within the seven-mode band", ylabel="Eigenvalue / band mean", title="L4: 3+4 cluster versus nearly sevenfold equality")
axes[1].legend(fontsize=8)
for suffix in ("png", "pdf"):
    fig.savefig(HERE / f"geometry_control.{suffix}", dpi=160)

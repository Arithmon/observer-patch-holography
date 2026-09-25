"""Render the frozen clock experiment, retaining the negative control failures."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / "receipt.json").read_text())
rows = [r for r in data["rows"] if r["q"] == 34 and r["law"] == "euclidean" and r["status"] == "timelike_and_graph_related"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), gridspec_kw={"width_ratios": [1.2, 1]})
ax = axes[0]
beta = np.linspace(0, .94, 200)
ax.plot(beta, np.sqrt(1-beta**2), color="#162f41", label="Imported Minkowski proper-time ratio", lw=2)
ax.axhline(1, color="#a44a3f", ls="--", lw=1.5, label="Uniform history / rank clock")
ax.scatter([r["actual_beta"] for r in rows], [r["volume_clock_ratio"] for r in rows], marker="o", facecolors="none", edgecolors="#8a9ea8", s=42, label="Raw volume clock")
ax.scatter([r["actual_beta"] for r in rows], [r["corrected_volume_clock_ratio"] for r in rows], marker="^", color="#177e89", s=42, label="Layer-corrected volume clock")
ax.set(xlabel="Actual endpoint displacement / coordinate duration", ylabel="Moving / rest clock ratio", title="Same five layers, different clock readings", xlim=(-.025,.95), ylim=(.3,1.08))
ax.legend(frameon=False, fontsize=8, loc="lower left")
ax.text(.02,.96,"q = 34 · three directions",transform=ax.transAxes,va="top",fontsize=9)

ax = axes[1]
summary = data["summary"]["q34_primary_comparison"]
values = [summary["euclidean"][k]*100 for k in ["median_abs_rank_relative_error", "median_abs_volume_relative_error", "median_abs_corrected_volume_relative_error"]]
values.append(summary["axis_only_control"]["median_abs_corrected_volume_relative_error"]*100)
labels = ["History / rank", "Raw volume", "Corrected\nvolume", "Axis-only\ncontrol"]
bars = ax.bar(range(4), values, color=["#a44a3f", "#8a9ea8", "#177e89", "#c89952"], width=.67)
for bar, value in zip(bars, values):
    ax.text(bar.get_x()+bar.get_width()/2, value+.4, f"{value:.2f}%", ha="center", fontsize=10)
ax.set_xticks(range(4), labels)
ax.set(ylabel="Median absolute relative error (%)", title="Calibration separates the read laws", ylim=(0,20))
ax.text(.02,.98,"Primary targets: β = 0.25, 0.5, 0.75\nEuclidean: 9/9 related · control: 6/9 related",transform=ax.transAxes,va="top",fontsize=8)
fig.suptitle("Observer clock experiment: a conditional failure and a calibrated alternative",fontsize=13,fontweight="bold",y=1.01)
fig.text(.01,-.04,"Exploratory, fixed specification. Supplied spatial geometry and read law; no physical clock identification. Near-cone failures retained in the receipt.",fontsize=8,color="#465761")
fig.tight_layout()
fig.savefig(ROOT/"clock_comparison.png",dpi=180,bbox_inches="tight")
fig.savefig(ROOT/"clock_comparison.svg",bbox_inches="tight")

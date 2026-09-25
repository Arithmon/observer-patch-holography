"""Plot the receipt and clearly labeled analytic dimension controls."""
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("native_history_plot_formula",HERE/"run.py")
native=importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)
receipt=json.loads((HERE/"receipt.json").read_text())
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
fig,axes=plt.subplots(2,2,figsize=(11,8),layout="constrained")
ax=axes[0,0]
graph=receipt["actual_graphs"][-1]
eigen=np.array(graph["positive_eigenvalues"])
ax.loglog(eigen,graph["mode_gff_limit_variances"],"--",color="#333333",label="Infinite-window limit")
for row,color in zip(graph["windows"][1:],["#7d98ad","#118a89","#df813c"]):
    ax.loglog(eigen,row["mode_variances"],".-",color=color,label=f"T = {row['window_over_slowest_relaxation']:g} relaxation times")
ax.axhline(graph["stationary_mode_variance"],color="#a7a7a7",linestyle=":",label="Instantaneous load")
ax.set(title="Actual L2 seam graph: a new readback field",xlabel="Graph eigenvalue λ",ylabel="Mode variance")
ax.legend(fontsize=8)
ax=axes[0,1]
u=np.geomspace(.001,1000,400)
fraction=np.array([native.occupation_variance(1,x)/2 for x in u])
ax.semilogx(u,fraction,color="#118a89",linewidth=2)
ax.axhline(.99,color="#888888",linestyle="--",linewidth=1)
ax.axvline(100,color="#888888",linestyle=":",linewidth=1)
ax.set(title="Long memory is required",xlabel="Mode age  aT = λT/2",ylabel="Fraction of inverse-Laplacian variance",ylim=(0,1.035))
ax.annotate("99% at aT ≈ 100",xy=(100,.99),xytext=(.09,.75),arrowprops={"arrowstyle":"->","color":"#555555"})
ax=axes[1,0]
for chain,color in zip(receipt["chains"],["#118a89","#df813c"]):
    ax.semilogx([x["aT"] for x in chain["windows"]],[x["kurtosis"] for x in chain["windows"]],"o-",label=chain["name"],color=color)
ax.axhline(3,color="#333333",linestyle="--",label="Gaussian")
ax.set(title="Exact finite-chain fourth moments",xlabel="Mode age  aT",ylabel="Kurtosis",ylim=(2.15,3.06))
ax.legend(fontsize=9)
ax=axes[1,1]
k=np.geomspace(.001,3,400)
power=np.array([native.occupation_variance(q*q/2,10000) for q in k])
pivot=native.occupation_variance(.5,10000)
for dimension,color in [(2,"#118a89"),(3,"#df813c")]:
    ax.loglog(k,k**dimension*power/pivot,label=f"Assumed {dimension}D continuum",color=color,linewidth=2)
ax.set(title="Same temporal mechanism, different dimension",xlabel="Wave number k (arbitrary units)",ylabel="Power per log k, normalized at k=1")
ax.legend(fontsize=9)
fig.suptitle("Conservative repair + temporal readback",fontsize=16)
fig.savefig(HERE/"native_history.png",dpi=180)
fig.savefig(HERE/"native_history.pdf")
print("Wrote native-history figures")

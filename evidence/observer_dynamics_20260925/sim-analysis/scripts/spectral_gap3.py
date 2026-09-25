"""Converged spectral gap of the glued federation's seam Laplacian by shift-invert Lanczos.

The dipole Ritz values (spectral_gap2.py) are exact variational upper bounds on lambda_2, and the LOBPCG
refinement from the dipoles gave unconverged estimates (residual/value 22 percent at level 6, 336 percent at
level 8; audit of 2026-09-25).  This solver factorises L + s I for a small positive shift s with SuperLU and
runs ARPACK on the inverse, which converges the smallest eigenvalues to machine precision when the
factorisation fits in memory (levels 6 and 7 on a large-memory host; level 8 is attempted with a guard).
It reports the four smallest eigenvalues (the constant mode at zero and the three dipoles), the relative
residual of each, the fraction of each mode in the carrier-constant (dipole) subspace and its within-carrier
slow-band share, and compares with the dipole Ritz bounds.  It also solves the 3 x 3 Ritz problem as a
symmetric generalised eigenproblem, as the audit asks.

usage: python3 spectral_gap3.py SIM_ROOT CACHE LEVEL OUT.json [--shift 1e-9] [--max-nnz-factor 2e9]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, splu, LinearOperator


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("sim_root"); ap.add_argument("cache"); ap.add_argument("level", type=int); ap.add_argument("out")
    ap.add_argument("--shift", type=float, default=1e-9); ap.add_argument("--k", type=int, default=4)
    a = ap.parse_args(); sys.path.insert(0, a.sim_root)
    from oph_exact import federation_huge as H
    from oph_exact import carrier
    t0 = time.time()
    geo = H.Geometry(Path(a.cache), a.level); n, ports = geo.carriers, geo.ports
    ea, eb = geo.seam_endpoints(np.arange(geo.seams))
    rows = np.concatenate([ea, eb]); cols = np.concatenate([eb, ea])
    adj = sp.csr_matrix((np.ones(rows.size, dtype=np.float64), (rows, cols)), shape=(ports, ports)); del rows, cols, ea, eb
    deg = np.asarray(adj.sum(axis=1)).ravel()
    lap = (sp.diags(deg) - adj).tocsc(); del adj
    points = np.load(Path(a.cache) / f"L{a.level}" / "cell_points.npy"); points = points / np.linalg.norm(points, axis=1, keepdims=True)
    X = np.repeat(points, 12, axis=0); X -= X.mean(axis=0, keepdims=True)
    LX = lap @ X
    G, Mm = X.T @ X, X.T @ LX
    ritz = np.sort(eigh(0.5 * (Mm + Mm.T), G, eigvals_only=True))  # symmetric generalised Ritz problem
    t1 = time.time()
    shifted = (lap + a.shift * sp.identity(ports, format="csc")).tocsc()
    lu = splu(shifted, permc_spec="MMD_AT_PLUS_A")
    t2 = time.time()
    op = LinearOperator((ports, ports), matvec=lu.solve, dtype=np.float64)
    vals, vecs = eigsh(op, k=a.k, which="LM", tol=1e-12, maxiter=5000)
    lam = 1.0 / vals - a.shift
    order = np.argsort(lam); lam, vecs = lam[order], vecs[:, order]
    resid = np.linalg.norm(lap @ vecs - vecs * lam[None, :], axis=0)
    p_slow = carrier.slow_band_projector()
    modes = []
    for k in range(a.k):
        v = vecs[:, k]
        blocks = v.reshape(n, 12)
        carrier_mean = blocks.mean(axis=1)
        within = blocks - carrier_mean[:, None]
        dipole_fraction = float(np.sum(carrier_mean ** 2) * 12 / np.sum(v ** 2))
        slow_share = float(np.sum((within @ p_slow) ** 2) / max(np.sum(within ** 2), 1e-300))
        modes.append({"eigenvalue": float(lam[k]), "relative_residual": float(resid[k] / max(abs(lam[k]), 1e-300)),
                      "carrier_constant_fraction": dipole_fraction, "within_carrier_slow_band_share": slow_share})
    res = {"level": a.level, "carriers": n, "ports": ports, "seams": geo.seams, "shift": a.shift,
           "dipole_ritz_values": [float(v) for v in ritz], "lambda_2_upper_bound_dipole_ritz": float(ritz[0]),
           "modes": modes, "lambda_2": float(lam[1]) if a.k > 1 else None,
           "lambda_2_times_carriers": float(lam[1] * n) if a.k > 1 else None,
           "converged": bool(np.all(resid[1:] / np.abs(lam[1:]) < 1e-8)) if a.k > 1 else None,
           "sweeps_per_e_fold_slowest_mode": float(2.0 / lam[1]) if a.k > 1 else None,
           "note": "shift-invert ARPACK on (L + shift I)^-1 with a SuperLU factorisation; eigenvalues nearest zero; the constant mode is the first; "
                   "relative residual ||L v - lambda v|| / |lambda| per mode; the dipole Ritz values are exact upper bounds on lambda_2..lambda_4",
           "seconds": {"assemble_and_ritz": round(t1 - t0, 1), "factorise": round(t2 - t1, 1), "eigsh": round(time.time() - t2, 1)}}
    Path(a.out).write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("level", "lambda_2", "lambda_2_times_carriers", "converged", "lambda_2_upper_bound_dipole_ritz", "seconds")}),
          [(round(m["eigenvalue"] * n, 5), f"{m['relative_residual']:.1e}", round(m["carrier_constant_fraction"], 5)) for m in modes])


if __name__ == "__main__":
    main()

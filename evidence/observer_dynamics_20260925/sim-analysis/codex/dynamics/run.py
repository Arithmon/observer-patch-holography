"""Bounded finite seam-graph dynamics experiment; writes small JSON only.

The spec is fixed before execution. Production geometry is imported only after
checking dependency bytes against the pinned Git commit. The federation module
itself is evaluated from that commit rather than the mutable working tree.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent
SOURCES = ("oph_exact/federation.py", "oph_exact/carrier.py", "oph_fpe/core/icosahedral.py",
           "oph_fpe/core/screen_ports.py", "oph_fpe/gauge/covariant_overlap.py",
           "oph_fpe/finite_groups.py", "oph_fpe/dynamics/self_readback_repair_closure.py")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def pinned_source(repo, commit):
    pins, blobs = {}, {}
    for path in SOURCES:
        blob = subprocess.check_output(["git", "-C", str(repo), "show", f"{commit}:{path}"])
        pins[path] = sha(blob)
        blobs[path] = blob
        if path != SOURCES[0] and (repo / path).read_bytes() != blob:
            raise RuntimeError(f"Working dependency differs from frozen source: {path}")
    sys.path.insert(0, str(repo))
    name = "codex_pinned_federation"
    module = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, loader=None))
    module.__file__ = str(repo / SOURCES[0])
    sys.modules[name] = module
    exec(compile(blobs[SOURCES[0]], f"{commit}:{SOURCES[0]}", "exec"), module.__dict__)
    return module, pins


def heat_rows(eigenvalues, size, max_eigenvalue, times, *, tail_tolerance=.001):
    """Return trace/dimension enclosures for ordered lowest eigenvalues.

    Missing eigenvalues are >= the last retained value. Residual errors of the
    supplied eigensolver are recorded separately: these are analytic truncation
    bounds on numerical eigenvalues, not interval-arithmetic certificates.
    """
    eigenvalues = np.maximum(np.asarray(eigenvalues), 0.)
    missing = size - len(eigenvalues)
    rows = []
    for t in times:
        s = t / 2
        ex = np.exp(-s * eigenvalues)
        known, numerator = float(ex.sum()), float(np.dot(eigenvalues, ex))
        tail = missing * math.exp(-s * eigenvalues[-1])
        maximizer = min(max(eigenvalues[-1], 1 / s), max_eigenvalue)
        numerator_tail = missing * maximizer * math.exp(-s * maximizer)
        rows.append({"sweeps": float(t), "gap_scaled_time": float(s * eigenvalues[1]),
                     "return_probability_lower": known / size,
                     "return_probability_upper": (known + tail) / size,
                     "tail_over_known_trace_bound": tail / known,
                     "dimension_lower": 2 * s * numerator / (known + tail),
                     "dimension_upper": 2 * s * (numerator + numerator_tail) / known,
                     "truncation_accepted": bool(tail / known <= tail_tolerance)})
    return rows


def torus_eigenvalues(side, dimension):
    one = 2 - 2 * np.cos(2 * np.pi * np.arange(side) / side)
    values = one
    for _ in range(dimension - 1):
        values = (values[:, None] + one).ravel()
    return np.sort(values)


def controls():
    result = []
    for side, dimension in ((256, 1), (32, 2), (12, 3)):
        eigen = torus_eigenvalues(side, dimension)
        rows = heat_rows(eigen, len(eigen), 4 * dimension, [4., 8., 16., 32.])
        result.append({"graph": f"flat_torus_{side}^{dimension}", "known_dimension": dimension, "heat": rows})
    return result


def static_schur(lap, carriers):
    matrix = lap.toarray()
    u = np.kron(np.eye(carriers), np.ones((12, 1)) / np.sqrt(12))
    local_fast = la.null_space(np.ones((1, 12)))
    r = np.kron(np.eye(carriers), local_fast)
    a, b, c = u.T @ matrix @ u, u.T @ matrix @ r, r.T @ matrix @ r
    solve = la.solve(c, b.T, assume_a="pos")
    effective = a - b @ solve
    lifted = u - r @ solve
    residual = np.linalg.norm(r.T @ matrix @ lifted, ord=2)
    exact = la.eigvalsh(matrix)
    small = la.eigvalsh(effective)
    naive = la.eigvalsh(a)
    return {"carriers": carriers, "fast_block_min_eigenvalue": float(la.eigvalsh(c)[0]),
            "fast_stationarity_residual": float(residual),
            "constant_mode_residual": float(np.linalg.norm(effective @ np.ones(carriers))),
            "effective_symmetry_residual": float(np.linalg.norm(effective - effective.T)),
            "naive_coarse_first_positive": naive[1:10].tolist(),
            "static_schur_first_positive": small[1:10].tolist(),
            "full_first_positive": exact[1:10].tolist(),
            "first_triplet_naive_relative_error": (naive[1:4] / exact[1:4] - 1).tolist(),
            "first_triplet_schur_relative_error": (small[1:4] / exact[1:4] - 1).tolist(),
            "scope": "Static elimination at level0; frequency-dependent memory and refinement limit not inferred."}


def experiment_level(federation, level, spec):
    start = time.monotonic()
    f = federation.build_federation(level, "port_pair")
    lap = f.laplacian.astype(float)
    k = spec["low_eigenpairs_by_level"][str(level)]
    if level == 0:
        vals, vecs = la.eigh(lap.toarray(), subset_by_index=(0, k - 1))
    else:
        vals, vecs = sla.eigsh(lap, k=k, sigma=-1e-6, which="LM", tol=1e-10,
                             maxiter=4000, v0=np.random.default_rng(725 + level).normal(size=f.ports))
        order = np.argsort(vals)
        vals, vecs = vals[order], vecs[:, order]
    residuals = np.linalg.norm(lap @ vecs - vecs * vals[None, :], axis=0)
    relative = residuals[1:] / vals[1:]
    passed = bool(np.all(relative <= spec["eigensolver_relative_residual_max"]))
    if vals[1] <= 0:
        raise RuntimeError("Connected-graph gap not resolved")
    points = federation.geodesic_icosahedral_patch_arrays(level, patch_basis="cells")[0]
    trial = np.repeat(points, 12, axis=0)
    trial -= trial.mean(axis=0)
    gram, stiffness = trial.T @ trial, trial.T @ (lap @ trial)
    ritz = la.eigvalsh(stiffness, gram)
    # Euclidean projection on all carrier means; no 'exactly constant' assumption.
    modes = vecs[:, 1:4].reshape(f.carriers, 12, 3)
    fast = modes - modes.mean(axis=1, keepdims=True)
    fast_fractions = np.sum(fast ** 2, axis=(0, 1)) / np.sum(modes ** 2, axis=(0, 1))
    time_values = np.geomspace(.01 / vals[1], 10 / vals[1], 121)
    heat = heat_rows(vals, f.ports, 2 * float(lap.diagonal().max()), time_values)
    window = [r for r in heat if r["truncation_accepted"] and .03 <= r["gap_scaled_time"] <= .3]
    bands = []
    for ell in (1, 2, 3):
        lo, hi = ell ** 2, (ell + 1) ** 2
        cluster = vals[lo:hi]
        bands.append({"sphere_index_band": ell, "indices_including_zero_at_index0": [lo, hi - 1],
                      "count": len(cluster), "mean": float(cluster.mean()),
                      "splitting_fraction": float(np.ptp(cluster) / cluster.mean()),
                      "mean_times_carriers_over_ell_ellplus1": float(cluster.mean() * f.carriers / (ell * (ell + 1)))})
    schur = static_schur(lap, f.carriers) if level == 0 else None
    isolated = None
    if level == 0:
        control = federation.build_federation(0, "isolated")
        spectrum = la.eigvalsh(control.laplacian.toarray())
        isolated = {"components": control.components, "zero_eigenvalues": int(np.sum(abs(spectrum) < 1e-10)),
                    "positive_min": float(spectrum[spectrum > 1e-10][0]),
                    "expected_carrier_gap": 5 - math.sqrt(5)}
    eigenrate = float(vals[1] / 2)
    tau = np.linspace(0, 10, 41)
    # Dimensionless histories use the same first positive spatial eigenmode.
    # The wave control adds an independent velocity coordinate; not a repair-law output.
    histories = {"dimensionless_times": tau.tolist(), "diffusion": np.exp(-tau).tolist(),
                 "monotone_quadratic_clock": np.exp(-tau ** 2).tolist(),
                 "monotone_log_clock": np.exp(-np.log1p(tau)).tolist(),
                 "wave_extension_zero_initial_velocity": np.cos(tau).tolist(),
                 "diffusion_gap_rate_per_sweep": eigenrate,
                 "wave_control_scope": "Mode-normalized illustrative added dynamics d2u/dt2=-c2 L u; c and velocity are additional inputs."}
    return {"level": level, "carriers": f.carriers, "ports": f.ports, "seams": f.seams,
            "seam_endpoints_sha256": sha(np.stack([f.seam_a, f.seam_b], axis=1).astype("<i8").tobytes()),
            "sparse_matrix_bytes": lap.data.nbytes + lap.indices.nbytes + lap.indptr.nbytes,
            "low_eigenvalues": vals.tolist(), "residual_norms": residuals.tolist(),
            "maximum_positive_relative_residual": float(relative.max()), "numerically_resolved": passed,
            "orthogonality_max_error": float(np.max(abs(vecs.T @ vecs - np.eye(k)))),
            "gap_times_carriers": float(vals[1] * f.carriers),
            "dipole_ritz_values": ritz.tolist(), "dipole_ritz_to_true_low_mode_ratios": (ritz / vals[1:4]).tolist(),
            "true_low_mode_within_carrier_energy_fractions": fast_fractions.tolist(),
            "sphere_band_diagnostics": bands, "heat": heat,
            "spectral_dimension_window": {"points": len(window),
                "lower_min": min((r["dimension_lower"] for r in window), default=None),
                "upper_max": max((r["dimension_upper"] for r in window), default=None),
                "midpoint_mean": float(np.mean([(r["dimension_lower"] + r["dimension_upper"]) / 2 for r in window])) if window else None,
                "status": "resolved finite window" if window and passed else "insufficient window or unresolved eigensolve"},
            "static_schur": schur, "isolated_control": isolated, "histories": histories,
            "seconds": time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=HERE.parents[2] / "oph-physics-sim")
    parser.add_argument("--spec", type=Path, default=HERE / "spec.json")
    parser.add_argument("--out", type=Path, default=HERE / "receipt.json")
    args = parser.parse_args()
    spec = json.loads(args.spec.read_text())
    producer, pins = pinned_source(args.repo.resolve(), spec["source_commit"])
    receipt = {"schema": "oph.codex-observables.dynamics.v1", "spec": spec,
               "spec_sha256": sha(args.spec.read_bytes()), "script_sha256": sha(Path(__file__).read_bytes()),
               "source_pins": pins, "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
               "controls": controls(), "levels": []}
    for level in spec["levels"]:
        try:
            result = experiment_level(producer, level, spec)
        except Exception as error:
            result = {"level": level, "status": "numerical failure", "exception": repr(error)}
        receipt["levels"].append(result)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        print(json.dumps({key: result.get(key) for key in ("level", "numerically_resolved", "gap_times_carriers", "maximum_positive_relative_residual", "spectral_dimension_window", "seconds", "exception")}), flush=True)


if __name__ == "__main__":
    main()

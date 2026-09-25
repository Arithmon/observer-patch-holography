"""Recover native second moments from archived refinement regression receipts.

This replays only the declared fresh RNG, never the repair simulation. Spatial
OLS residuals are not nonlinear conditional innovations or geometric curvature.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIM = HERE.parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def replay_fresh(carriers, seed, mode, sizes, chunk_cells=65536, deadline=None):
    """Exact integer accumulation; no load array larger than one chunk.

    Returned variance uses 1/G, as does the original ratio estimator. Group
    variance is computed around the realized mean, not the ensemble mean 2.
    """
    if mode not in {"chain", "shuffled"}:
        raise ValueError("unsupported archived mode")
    if any(carriers % size or chunk_cells % size for size in sizes):
        raise ValueError("every group size must divide carriers and chunk size")
    if chunk_cells * 12 * 8 > 25_000_000 or carriers // 4 * 8 > 25_000_000:
        raise ValueError("working-array resource bound exceeded")
    rng = np.random.default_rng(seed)
    if mode == "shuffled":
        # The producer permutes entire parent carriers before drawing f.
        permutation = rng.permutation(carriers // 4)
        del permutation
    moments = {size: [0, 0, 0] for size in sizes}
    stream_hash = hashlib.sha256()
    for start in range(0, carriers, chunk_cells):
        if deadline is not None and time.monotonic() > deadline:
            raise TimeoutError("frozen local extraction time limit exceeded")
        count = min(chunk_cells, carriers - start)
        loads = rng.integers(0, 5, size=count * 12)
        stream_hash.update(memoryview(loads).cast("B"))
        loads -= 2
        for size in sizes:
            sums = loads.reshape(-1, size * 12).sum(axis=1, dtype=np.int64)
            groups, total, squares = moments[size]
            moments[size] = [groups + len(sums), total + int(sums.sum()),
                             squares + int(sums @ sums)]
    results = {}
    for size, (groups, total, squares) in moments.items():
        ports = 12 * size
        mean_delta = total / (groups * ports)
        variance = squares / (groups * ports * ports) - mean_delta**2
        results[size] = {
            "groups": groups, "fresh_mean": 2 + mean_delta,
            "fresh_variance": variance,
            "iid_expected_centered_variance": (1 - 1 / groups) * 2 / ports,
            "integer_group_sum_total": total,
            "integer_group_sum_squares": squares,
        }
    return results, stream_hash.hexdigest()


def recover(row, fresh_variance):
    """Invert the two normal equations and retain absolute covariance units."""
    b, B, C = (float(row[k]) for k in
               ("response_b", "joint_response_b", "joint_fresh_c"))
    rff = float(row["fresh_power_over_inherited"])
    if not all(math.isfinite(v) for v in [b, B, C, rff, fresh_variance]):
        raise ValueError("nonfinite covariance input")
    if rff <= 0 or fresh_variance <= 0 or abs(C) < 1e-12:
        raise ValueError("degenerate joint regression")
    vx = fresh_variance / rff
    xf = vx * (b - B) / C
    xy = b * vx
    fy = B * xf + C * fresh_variance
    vy = float(row["total_power_ratio"]) * vx
    innovation = float(row["innovation_power_ratio"]) * vx
    joint = float(row["joint_innovation_power_ratio"]) * vx
    covariance = np.array([[vx, xf, xy], [xf, fresh_variance, fy], [xy, fy, vy]])
    scale = float(np.max(np.abs(covariance)))
    tolerance = 2e-9 * scale
    if np.linalg.eigvalsh(covariance).min() < -tolerance:
        raise ValueError("covariance is not positive semidefinite")
    if innovation < -tolerance or joint < -tolerance:
        raise ValueError("negative residual variance")
    inherited_piece = b * b * vx
    fresh_orthogonal_piece = C * C * (fresh_variance - xf * xf / vx)
    errors = {
        "simple_decomposition": abs(vy - inherited_piece - innovation),
        "joint_decomposition": abs(vy - B * B * vx - C * C * fresh_variance
                                   - 2 * B * C * xf - joint),
        "innovation_fresh_plus_residual": abs(innovation - fresh_orthogonal_piece - joint),
        "coherence_identity": abs(float(row["squared_coherence"]) * vx * vy - xy * xy) / scale,
    }
    if max(errors.values()) > tolerance:
        raise ValueError("archived moment/decomposition inconsistency")
    return {
        "cells_per_group": int(row["cells_per_group"]), "groups": int(row["groups"]),
        "covariance_order": ["inherited_x", "fresh_f", "settled_y"],
        "covariance_native_load_squared": covariance.tolist(),
        "inherited_variance": vx, "fresh_variance": fresh_variance,
        "settled_variance": vy, "response_b": b,
        "joint_response_b": B, "joint_fresh_c": C,
        "inherited_linear_residual_variance": innovation,
        "inherited_linear_residual_rms": math.sqrt(max(innovation, 0)),
        "joint_linear_residual_variance": joint,
        "fresh_orthogonal_transferred_variance": fresh_orthogonal_piece,
        "inherited_response_power": inherited_piece,
        "linear_residual_fraction_of_settled": innovation / vy,
        "max_relative_algebra_error": max(errors.values()) / scale,
    }


def descriptive(values):
    values = np.asarray(values, dtype=float)
    return {"mean": float(values.mean()), "minimum": float(values.min()),
            "maximum": float(values.max()),
            "sample_sd": float(values.std(ddof=1)) if len(values) > 1 else None}


def aggregate(chains):
    """Descriptive parent-schedule clusters, never iid cosmological errors."""
    metrics = ["inherited_variance", "settled_variance", "inherited_linear_residual_variance",
               "joint_linear_residual_variance", "fresh_orthogonal_transferred_variance",
               "joint_response_b", "joint_fresh_c"]
    result = []
    for level, mode in sorted({(c["parent_level"], c["mode"]) for c in chains}):
        selected = [c for c in chains if (c["parent_level"], c["mode"]) == (level, mode)]
        block = {"parent_level": level, "child_level": level + 1, "mode": mode,
                 "chains": len(selected), "parent_schedule_clusters": len({c["parent_seed"] for c in selected}),
                 "independent_initial_load_realizations": 1,
                 "uncertainty": "Descriptive spread across parent schedules; common parent initial load, no population SE.",
                 "by_scale": []}
        for size in sorted({r["cells_per_group"] for c in selected for r in c["by_scale"]}):
            parent_means = {}
            for c in selected:
                row = next(r for r in c["by_scale"] if r["cells_per_group"] == size)
                parent_means.setdefault(c["parent_seed"], []).append(row)
            summary = {"cells_per_group": size}
            for key in metrics:
                means = [np.mean([r[key] for r in rows]) for rows in parent_means.values()]
                summary[key] = descriptive(means)
            block["by_scale"].append(summary)
        rows = block["by_scale"]
        for lo, hi in zip(rows[:-1], rows[1:]):
            # This is a block-size slope, not a harmonic/primordial tilt.
            hi["adjacent_log_variance_slope_vs_group_cells"] = math.log(
                hi["inherited_linear_residual_variance"]["mean"] /
                lo["inherited_linear_residual_variance"]["mean"]
            ) / math.log(hi["cells_per_group"] / lo["cells_per_group"])
        result.append(block)
    return result


def build(spec_path=HERE / "spec.json"):
    spec = json.loads(Path(spec_path).read_text())
    for entry in spec["inputs"]:
        if sha(SIM / entry["path"]) != entry["sha256"]:
            raise ValueError(f"input digest mismatch: {entry['path']}")
    actual = sorted(p.relative_to(SIM).as_posix() for p in (SIM / "data/refine_ensemble").glob("chain_*.json"))
    wanted = [e["path"] for e in spec["inputs"] if e["role"] == "refinement_chain"]
    if actual != wanted:
        raise ValueError("complete archive inventory changed")
    anchor = json.loads((SIM / "data/codex_audit_20260925/scale_ensemble.json").read_text())
    parent_variances = {(lev["level"], row["seed"], int(scale["name"].split("_")[1])): row["terminal_variance"]
                        for lev in anchor["levels"] for scale in lev["scales"]
                        if scale["name"].startswith("ancestor_") for row in scale["per_schedule"]}
    chains, started = [], time.monotonic()
    max_anchor_error = 0.
    checked_anchor_rows = 0
    for entry in (e for e in spec["inputs"] if e["role"] == "refinement_chain"):
        data = json.loads((SIM / entry["path"]).read_text())
        if not data["terminated"]:
            raise ValueError("archived chain did not terminate")
        sizes = [r["cells_per_group"] for r in data["ancestral_groups"]]
        if sizes != spec["group_cells"]:
            raise ValueError("frozen complete scale inventory differs")
        moments, stream_hash = replay_fresh(data["carriers"], data["fresh_seed"], data["mode"], sizes,
                                           spec["chunk_cells"], started + spec["wall_limit_seconds"])
        chain = {k: data[k] for k in ["parent_level", "child_level", "parent_seed", "fresh_seed", "mode",
                                      "child_schedule_seed", "carriers", "child_terminal_sha256"]}
        chain.update({"source_path": entry["path"], "fresh_replay_int64_stream_sha256": stream_hash, "by_scale": []})
        for row in data["ancestral_groups"]:
            size = row["cells_per_group"]
            if row["groups"] != moments[size]["groups"]:
                raise ValueError("inconsistent group count")
            recovered = recover(row, moments[size]["fresh_variance"])
            recovered["fresh_replay"] = moments[size]
            key = (data["parent_level"], data["parent_seed"], size // 4)
            if data["mode"] == "chain" and key in parent_variances:
                parent_variance = parent_variances[key]
                relative_error = abs(recovered["inherited_variance"] / parent_variance - 1)
                if relative_error > spec["parent_anchor_relative_tolerance"]:
                    raise ValueError("fresh RNG reconstruction disagrees with independent parent variance")
                recovered["independent_parent_anchor_relative_error"] = relative_error
                checked_anchor_rows += 1
                max_anchor_error = max(max_anchor_error, relative_error)
            chain["by_scale"].append(recovered)
        chains.append(chain)
        print(f"replayed {len(chains)}/{len(wanted)} receipts, {time.monotonic()-started:.1f}s", flush=True)
    return {
        "schema": "oph.native-source-innovation-extraction.v1",
        "spec_sha256": sha(spec_path), "producer_sha256": sha(__file__),
        "versions": {"numpy": np.__version__, "python": platform.python_version()},
        "status": "NATIVE_LINEAR_RESIDUAL_MOMENTS_RECOVERED_NOT_PRIMORDIAL_CURVATURE",
        "normalization": "Exact finite-seed replay conditional on archived producer RNG semantics; population centering 1/G.",
        "units": "squared integer port-load units after cell and ancestral averaging",
        "boundaries": spec["boundaries"],
        "independent_parent_anchor_rows": checked_anchor_rows,
        "independent_parent_anchor_max_relative_error": max_anchor_error,
        "chains": chains, "aggregate": aggregate(chains),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Recompute all RNG moments and compare receipt byte for byte")
    args = parser.parse_args()
    result = build()
    destination = HERE / "receipt.json"
    if args.check:
        if json.loads(destination.read_text()) != result:
            raise SystemExit("receipt differs from deterministic recomputation")
        print("receipt verified")
    else:
        write_json(destination, result)


if __name__ == "__main__":
    main()

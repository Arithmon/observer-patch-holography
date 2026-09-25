"""Read-only, source-pinned scale/ensemble diagnostics of first-hitting states.

Run on the existing AWS host, with BLAS thread counts set to one. Only the
small JSON receipt is written. Inputs remain in their original directories.
Ancestral groups have equal cell counts; cube groups provide a weighting
control. Equal-solid-angle z/phi pixels provide an independent angular
estimator at two resolutions. No field-to-CMB identification is assumed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.special import sph_harm_y


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def array_sha(array):
    h = hashlib.sha256()
    flat = array.reshape(-1)
    for start in range(0, flat.size, 8 * 1024 * 1024):
        h.update(np.ascontiguousarray(flat[start:start + 8 * 1024 * 1024]).tobytes())
    return h.hexdigest()


def centered(x, weights):
    return x - np.dot(weights, x) / weights.sum()


def decomposition(x, y, weights):
    """An intercept and the SAME weights/centering for every statistic."""
    x, y = centered(x, weights), centered(y, weights)
    w = weights / weights.sum()
    xx, xy, yy = np.dot(w, x * x), np.dot(w, x * y), np.dot(w, y * y)
    if xx <= 0 or yy <= 0:
        raise ValueError("nonpositive empirical variance")
    slope = xy / xx
    residual = np.dot(w, (y - slope * x) ** 2)
    assert np.isclose(yy, slope * slope * xx + residual, rtol=1e-9, atol=1e-14)
    return {
        "response": float(slope), "power_ratio": float(yy / xx),
        "coherence_squared": float(xy * xy / (xx * yy)),
        "residual_power_ratio": float(residual / xx),
        "initial_variance": float(xx), "terminal_variance": float(yy),
    }


def exchangeable_null_variance(group_ports, weights, total_ports, fraction):
    """E[weighted centered variance] under uniform fixed-total binary states.

    Cov(group_i,group_j)=c*(delta_ij/group_ports_i - 1/total_ports),
    c=p(1-p)*total_ports/(total_ports-1). Centering cancels the constant.
    Groups must be disjoint. They need not cover all ports.
    """
    a = weights / weights.sum()
    c = fraction * (1 - fraction) * total_ports / (total_ports - 1)
    return float(c * np.sum((a - a * a) / group_ports))


def iid_null_variance(group_ports, weights, variance=35 / 12):
    a = weights / weights.sum()
    return float(variance * np.sum((a - a * a) / group_ports))


class Scale:
    def __init__(self, name, groups, initial, weights, ports, fraction):
        self.name, self.groups, self.weights = name, groups, weights
        self.initial = initial
        self.sum_y = np.zeros_like(initial)
        self.sum_var = 0.0
        self.rows = []
        self.equilibrium = exchangeable_null_variance(groups * 12, weights, ports, fraction)
        self.initial_null = iid_null_variance(groups * 12, weights)

    def add(self, y, seed, global_x, global_y):
        row = decomposition(self.initial, y, self.weights)
        # This separate number reproduces the former global-centered,
        # unweighted definition; it is NOT combined with centered correlation.
        xg, yg = self.initial - global_x, y - global_y
        row["legacy_unweighted_global_centered_response"] = float(np.dot(xg, yg) / np.dot(xg, xg))
        row["seed"] = int(seed)
        row["variance_over_exchangeable_equilibrium"] = row["terminal_variance"] / self.equilibrium
        self.rows.append(row)
        self.sum_y += centered(y, self.weights)
        self.sum_var += row["terminal_variance"]

    def result(self):
        s = len(self.rows)
        assert s >= 2
        a = self.weights / self.weights.sum()
        sum_power = float(np.dot(a, self.sum_y ** 2))
        noise = (self.sum_var - sum_power / s) / (s - 1)
        common = (sum_power - self.sum_var) / (s * (s - 1))
        total = self.sum_var / s
        assert np.isclose(total, common + noise, rtol=1e-9, atol=1e-14)
        responses = np.array([r["response"] for r in self.rows])
        return {
            "name": self.name, "groups": int(self.groups.size),
            "cells_min": int(self.groups.min()), "cells_max": int(self.groups.max()),
            "cells_arithmetic_mean": float(self.groups.mean()),
            "cells_harmonic_mean": float(1 / np.mean(1 / self.groups)),
            "response_mean": float(responses.mean()),
            "response_schedule_sd": float(responses.std(ddof=1)),
            "response_schedule_se": float(responses.std(ddof=1) / np.sqrt(s)),
            "initial_variance_over_iid_null": self.rows[0]["initial_variance"] / self.initial_null,
            "equilibrium_variance": self.equilibrium,
            "total_terminal_variance": total,
            "conditional_schedule_variance_unbiased": noise,
            "cross_schedule_common_power_unbiased": common,
            "schedule_noise_fraction": noise / total,
            "per_schedule": self.rows,
        }


def cube_partition(directions, resolution):
    key = np.minimum(((directions + 1) * resolution / 2).astype(np.int32), resolution - 1)
    labels = key[:, 0] * resolution ** 2 + key[:, 1] * resolution + key[:, 2]
    counts = np.bincount(labels, minlength=resolution ** 3)
    keep = counts >= 4
    return labels, counts, keep


def angular_partition(directions, nz, nphi):
    iz = np.minimum(((directions[:, 2] + 1) * nz / 2).astype(np.int32), nz - 1)
    phi = np.mod(np.arctan2(directions[:, 1], directions[:, 0]), 2 * np.pi)
    ip = np.minimum((phi * nphi / (2 * np.pi)).astype(np.int32), nphi - 1)
    labels = iz * nphi + ip
    counts = np.bincount(labels, minlength=nz * nphi)
    if np.any(counts == 0):
        raise ValueError("angular grid has empty pixels; no interpolation permitted")
    return labels, counts


def angular_coefficients(pixel_means, nz, nphi, lmax):
    """Midpoint solid-angle quadrature, all m>=0, with exact longitude FFT."""
    grid = pixel_means.reshape(nz, nphi)
    grid = grid - grid.mean()
    ft = np.fft.rfft(grid, axis=1)
    theta = np.arccos(-1 + (np.arange(nz) + 0.5) * 2 / nz)
    out = []
    for ell in range(1, lmax + 1):
        a = []
        for m in range(ell + 1):
            y = sph_harm_y(ell, m, theta, np.zeros(nz)).real
            value = 4 * np.pi / (nz * nphi) * np.dot(y, ft[:, m]) * np.exp(-1j * m * np.pi / nphi)
            a.append(value)
        out.append(np.asarray(a))
    return out


def angular_summary(coefficients):
    x, terminals = coefficients[0], coefficients[1:]
    rows = []
    for i, initial in enumerate(x):
        ell = i + 1
        w = np.ones(ell + 1); w[1:] = 2
        xx = float(np.dot(w, abs(initial) ** 2) / (2 * ell + 1))
        ys = np.stack([a[i] for a in terminals])
        yy = np.sum(w * abs(ys) ** 2, axis=1) / (2 * ell + 1)
        xy = np.sum(w * (ys * initial.conj()).real, axis=1) / (2 * ell + 1)
        response = xy / xx
        residual = yy - xy ** 2 / xx
        assert residual.min() >= -1e-12
        noise = float(np.sum(w * np.var(ys, axis=0, ddof=1)) / (2 * ell + 1))
        rows.append({
            "ell": ell, "C_initial": xx, "C_terminal_mean": float(yy.mean()),
            "cross_response_mean": float(response.mean()),
            "cross_response_schedule_sd": float(response.std(ddof=1)),
            "total_power_ratio_mean": float(yy.mean() / xx),
            "coherence_squared_mean": float(np.mean(xy ** 2 / (xx * yy))),
            "orthogonal_residual_power_ratio_mean": float(residual.mean() / xx),
            "conditional_schedule_noise_power": noise,
            "cross_schedule_common_power": float(yy.mean() - noise),
            "per_schedule_response": response.tolist(),
            "per_schedule_power_ratio": (yy / xx).tolist(),
        })
    return rows


def build_level(run_root, geo_root, level, spec):
    started = time.monotonic()
    run = Path(run_root) / f"L{level}"
    rp = run / "receipt.json"
    receipt_bytes = rp.read_bytes()
    receipt = json.loads(receipt_bytes)
    n, ports = receipt["carriers"], receipt["ports"]
    entries = receipt["integer_law"]["entries"]
    assert len(entries) == 16 and ports == 12 * n
    expected_sum = receipt["expected"]["component_total"][0]
    q, remainder = divmod(expected_sum, ports)
    fraction = remainder / ports
    points_path = Path(geo_root) / f"L{level}/cell_points.npy"
    assert file_sha(points_path) == receipt["geometry"]["arrays"]["cell_points"]["sha256"]
    points = np.load(points_path, mmap_mode="r")
    directions = points / np.linalg.norm(points, axis=1, keepdims=True)
    loads = np.random.default_rng(20260909 + level).integers(0, 6, size=ports)
    assert int(loads.sum()) == expected_sum
    assert int(np.dot(loads, loads)) == entries[0]["V_initial"]
    initial_hash = array_sha(loads)
    initial = loads.reshape(n, 12).mean(axis=1)
    del loads
    global_mean = expected_sum / ports
    scales, ancestor, cube, angular, alms = {}, {}, {}, {}, {}
    for size in spec["primary_ancestral_group_cells"]:
        if n % size or n // size < 20:
            continue
        means = initial.reshape(-1, size).mean(axis=1)
        counts = np.full(means.size, size, dtype=np.int64)
        key = f"ancestor_{size}"
        scales[key] = Scale(key, counts, means, np.ones(means.size), ports, fraction)
        ancestor[size] = key
    for resolution in spec["cube_control_resolutions"]:
        labels, counts, keep = cube_partition(directions, resolution)
        means = np.bincount(labels, weights=initial, minlength=counts.size)[keep] / counts[keep]
        cube[resolution] = (labels, counts, keep)
        for mode, weights in (("weighted", counts[keep].astype(float)), ("unweighted", np.ones(keep.sum()))):
            key = f"cube_{resolution}_{mode}"
            scales[key] = Scale(key, counts[keep], means, weights, ports, fraction)
    for nz, nphi in spec["angular_grids"]:
        labels, counts = angular_partition(directions, nz, nphi)
        key = f"z{nz}_phi{nphi}"
        angular[key] = (labels, counts, nz, nphi)
        means = np.bincount(labels, weights=initial, minlength=counts.size) / counts
        alms[key] = [angular_coefficients(means, nz, nphi, spec["angular_lmax"])]
    del directions, points
    input_states = []
    for e in entries:
        path = run / f"integer_{e['seed']}" / e["terminal_state"]["path"]
        x = np.load(path, mmap_mode="r")
        digest = array_sha(x)
        assert digest == e["terminal_state"]["sha256"], path
        assert x.dtype == np.int8 and x.size == ports
        assert int(x.min()) == q and int(x.max()) == q + 1
        assert int(x.sum(dtype=np.int64)) == expected_sum
        field = x.reshape(n, 12).mean(axis=1)
        for size, key in ancestor.items():
            scales[key].add(field.reshape(-1, size).mean(axis=1), e["seed"], global_mean, global_mean)
        for resolution, (labels, counts, keep) in cube.items():
            means = np.bincount(labels, weights=field, minlength=counts.size)[keep] / counts[keep]
            for mode in ("weighted", "unweighted"):
                scales[f"cube_{resolution}_{mode}"].add(means, e["seed"], global_mean, global_mean)
        for key, (labels, counts, nz, nphi) in angular.items():
            means = np.bincount(labels, weights=field, minlength=counts.size) / counts
            alms[key].append(angular_coefficients(means, nz, nphi, spec["angular_lmax"]))
        input_states.append({"seed": e["seed"], "array_sha256": digest, "relative_path": str(path.relative_to(Path(run_root)))})
        del field, x
        print(f"L{level} verified and measured schedule {e['seed']} in {time.monotonic()-started:.1f}s", flush=True)
    assert rp.read_bytes() == receipt_bytes, "source receipt changed during analysis"
    return {
        "level": level, "carriers": n, "ports": ports, "schedules": len(entries),
        "initial_load_seed": 20260909 + level, "initial_load_dtype": str(np.dtype(np.int64)),
        "initial_load_array_sha256": initial_hash,
        "source_receipt_sha256": hashlib.sha256(receipt_bytes).hexdigest(),
        "geometry_file_sha256": file_sha(points_path), "inputs": input_states,
        "balanced_fraction": fraction, "scales": [v.result() for v in scales.values()],
        "angular": {key: {"counts_min": int(angular[key][1].min()), "counts_max": int(angular[key][1].max()), "rows": angular_summary(value)} for key, value in alms.items()},
        "seconds": time.monotonic() - started,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-root", type=Path, required=True)
    ap.add_argument("--geo-root", type=Path, required=True)
    ap.add_argument("--spec", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--levels", type=int, nargs="+")
    args = ap.parse_args()
    spec = json.loads(args.spec.read_text())
    result = {"schema": "oph.sim.codex-scale-ensemble.v1", "spec": spec,
              "spec_sha256": file_sha(args.spec), "script_sha256": file_sha(__file__),
              "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}, "levels": []}
    for level in args.levels or spec["levels"]:
        result["levels"].append(build_level(args.run_root, args.geo_root, level, spec))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.out.with_suffix(".tmp")
        tmp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        tmp.replace(args.out)


if __name__ == "__main__":
    main()

"""Metric-sensitive spectral controls on the same pinned triangular tower."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import time

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("codex_dynamics", HERE / "run.py")
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


def cotangent_fem(vertices, faces):
    """Linear FEM stiffness on planar chord triangles, with lumped area mass."""
    rows, columns, values = [], [], []
    mass = np.zeros(len(vertices))
    tri = vertices[faces]
    twice_area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    for corner in range(3):
        i, j, k = corner, (corner + 1) % 3, (corner + 2) % 3
        cotangent = np.sum((tri[:, j] - tri[:, i]) * (tri[:, k] - tri[:, i]), axis=1) / twice_area
        weight = .5 * cotangent
        a, b = faces[:, j], faces[:, k]
        rows.extend((a, b, a, b)); columns.extend((a, b, b, a)); values.extend((weight, weight, -weight, -weight))
        np.add.at(mass, faces[:, corner], twice_area / 6)
    stiffness = sp.coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))), shape=(len(vertices), len(vertices))).tocsr()
    return stiffness, sp.diags(mass).tocsr()


def bands(values):
    return [{"ell_index_band": ell, "values": values[ell ** 2:(ell + 1) ** 2].tolist(),
             "mean": float(values[ell ** 2:(ell + 1) ** 2].mean()),
             "splitting_fraction": float(np.ptp(values[ell ** 2:(ell + 1) ** 2]) / values[ell ** 2:(ell + 1) ** 2].mean())}
            for ell in (1, 2, 3)]


def solve(stiffness, mass=None):
    # Oversampling and a larger Krylov space protect complete symmetry multiplets.
    # The initial k=16/default-ncv diagnostic is retained separately; small
    # eigenpair residuals alone failed to establish the requested eigenvalue index.
    size = stiffness.shape[0]
    requested = min(32, size - 2)
    values, vectors = sla.eigsh(stiffness, M=mass, k=requested, ncv=min(size, 128), sigma=-1e-6, which="LM", tol=1e-11,
                               v0=np.random.default_rng(20260925).normal(size=stiffness.shape[0]))
    order = np.argsort(values); values, vectors = values[order], vectors[:, order]
    values, vectors = values[:16], vectors[:, :16]
    mv = vectors if mass is None else mass @ vectors
    residual = np.linalg.norm(stiffness @ vectors - mv * values, axis=0)
    relative = residual[1:] / np.linalg.norm(mv[:, 1:] * values[1:], axis=0)
    if size <= 1280:
        reference = la.eigvalsh(stiffness.toarray(), None if mass is None else mass.toarray(), subset_by_index=(0, 15))
        reference_method = "dense symmetric generalized eigh, first16"
    else:
        reference = np.sort(sla.eigsh(stiffness, M=mass, k=48, ncv=min(size, 192), sigma=-1e-6, which="LM", tol=1e-12,
                                     v0=np.random.default_rng(20260926).normal(size=size), return_eigenvectors=False))[:16]
        reference_method = "independent48-mode192-vector Krylov solve and starting vector"
    agreement = float(np.max(abs(values[1:] / reference[1:] - 1)))
    if agreement > 1e-7:
        raise RuntimeError(f"Low-spectrum completeness cross-check failed: relative disagreement {agreement}")
    return {"values": values.tolist(), "max_relative_residual": float(relative.max()), "bands": bands(values),
            "index_crosscheck": {"method": reference_method, "maximum_positive_relative_disagreement": agreement}}


def main():
    spec_path = HERE / "geometry_spec.json"
    spec = json.loads(spec_path.read_text())
    source, pins = RUN.pinned_source(HERE.parents[2] / "oph-physics-sim", spec["source_commit"])
    output = {"schema": "oph.codex-observables.geometry-control.v1", "spec": spec,
              "spec_sha256": RUN.sha(spec_path.read_bytes()), "script_sha256": RUN.sha(Path(__file__).read_bytes()),
              "source_pins": pins, "levels": [],
              "numerical_repair": "Initial default Krylov space missed exact-multiplet members despite tiny residuals; retained geometry_initial_solver_diagnostic.json. Final first16 spectrum uses32 modes/ncv128 and dense or independent48-mode completeness cross-checks."}
    for level in spec["levels"]:
        start = time.monotonic()
        mesh = source.build_geodesic_icosahedral_tower(level).levels[level]
        valence = np.bincount(mesh.edges.ravel(), minlength=len(mesh.vertices))
        values, counts = np.unique(valence, return_counts=True)
        deficits = (6 - valence) * np.pi / 3
        points, left, right = source.geodesic_icosahedral_patch_arrays(level, patch_basis="cells")
        adjacency = sp.coo_matrix((np.ones(2 * len(left)), (np.r_[left, right], np.r_[right, left])), shape=(len(points), len(points))).tocsr()
        lap = sp.diags(np.asarray(adjacency.sum(axis=1)).ravel()) - adjacency
        dual = solve(lap)
        fem_lap, mass = cotangent_fem(mesh.vertices, mesh.faces)
        # Level0 has only12 vertices, insufficient to isolate ell3; retain this
        # declared resolution failure rather than silently changing band indices.
        fem = solve(fem_lap, mass) if len(mesh.vertices) > 16 else {"status": "too few vertices for16 eigenpairs"}
        record = {"level": level, "vertices": len(mesh.vertices), "faces": len(mesh.faces),
                  "valence_counts": {str(v): int(c) for v, c in zip(values, counts)},
                  "nonzero_defect_vertices": np.flatnonzero(abs(deficits) > 1e-14).tolist(),
                  "equilateral_deficit_per_original_vertex": deficits[:12].tolist(),
                  "equilateral_total_deficit": float(deficits.sum()),
                  "lumped_chord_area": float(mass.diagonal().sum()),
                  "cell_dual": dual, "round_embedding_fem": fem, "seconds": time.monotonic() - start}
        output["levels"].append(record)
        (HERE / "geometry_receipt.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        print(json.dumps({"level": level, "defects": len(record["nonzero_defect_vertices"]),
                          "dual_ell3_split": dual["bands"][2]["splitting_fraction"],
                          "fem_ell3_split": fem.get("bands", [{}, {}, {}])[2].get("splitting_fraction"),
                          "seconds": record["seconds"]}), flush=True)


if __name__ == "__main__":
    main()

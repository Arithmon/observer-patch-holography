"""Conditional finite Maxwell stress/energy execution on the existing OPH cone.

This is an explicitly supplied Maxwell-action experiment, not a native repair
equation of state. The spatial geometry, unit constitutive coefficients,
canonical phase space and model time are inputs. Pressure means canonical
virtual work under metric deformation, or equivalently the volume-averaged
Maxwell momentum-flux tensor. No temperature, entropy, equilibration, material
wall energy or observer-memory thermodynamics is inferred.

For x -> L*x in three spatial dimensions, barycentric gradients scale L^-1
and volume scales L^3. Consequently the *assembled* Whitney one-form mass
M(L)=L*M(1) and two-form mass N(L)=L^-1*N(1). With edge potential a,
canonical momentum pi=-M*E, and curvature b=C*a,

    H = (pi^T M^-1 pi + b^T N b)/2.

At fixed canonical (a,pi), H(L)=H(1)/L, V(L)=L^3*V(1), hence
-dH/dV=H/(3V) for H>0. The code does not use this identity to produce
pressure: it rebuilds the geometry and differentiates its energy, then
compares with independently integrated stress. The result is a mean-stress
relation valid for anisotropic fields as well as isotropic preparations.
Directional pressure need not be positive. H=0 leaves p/rho undefined.

The finite closed Hamiltonian uses the natural variational boundary encoded
by the complete edge space, with no external source covectors. The constant
magnetic initial field is not assumed a stationary boundary solution.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import lu_factor, lu_solve, null_space, eigh

HERE = Path(__file__).resolve().parent
RER = HERE.parents[2]
EM = RER / "code/electromagnetism"
sys.path.insert(0, str(EM))
import cone_whitney_bridge as source_geometry
import verify_whitney_maxwell_dynamics as source_replay

SCHEMA = "oph.exploratory_whitney_maxwell_eos.v1"
SCOPE = "SUPPLIED_MAXWELL_ACTION__CANONICAL_METRIC_WORK__NO_NATIVE_REPAIR_EOS"
SOURCE_PATHS = [
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "code/electromagnetism/cone_whitney_bridge.py",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "code/electromagnetism/verify_whitney_maxwell_dynamics.py",
    "code/electromagnetism/runtime/whitney_maxwell_dynamics_receipt.json",
]
STEPS = 80
DT = 0.025
EPSILONS = [0.01, 0.001, 0.0001]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mesh_and_forms():
    edges, faces, d, c, exact_vertices = source_geometry.load_geometry()
    maps = source_geometry.cone_maps(d, c)
    mesh = {
        "vertices": [[0.0] * 3] + [[float(x) for x in v] for v in exact_vertices],
        "edges": [[0, i + 1] for i in range(12)] + [[u + 1, v + 1] for u, v in edges],
        "faces": [[i + 1 for i in f] for f in faces] + [[0, u + 1, v + 1] for u, v in edges],
        "tetrahedra": [[0, *(i + 1 for i in f)] for f in faces],
    }
    mesh["D"] = np.array(maps["D"].tolist(), dtype=int).tolist()
    mesh["C"] = np.array(maps["C"].tolist(), dtype=int).tolist()
    return mesh, forms(mesh, np.eye(3))


def forms(mesh, deformation):
    vertices = np.array(mesh["vertices"]) @ np.asarray(deformation).T
    m, n, local = source_geometry.whitney_matrices(
        vertices, mesh["edges"], mesh["faces"], mesh["tetrahedra"])
    return m, n, local, sum(row["volume"] for row in local)


def energy(a, pi, m, n, c):
    b = c @ a
    electric = float(pi @ np.linalg.solve(m, pi) / 2)
    magnetic = float(b @ n @ b / 2)
    return electric + magnetic, electric, magnetic


def stress_integral(e, b, local):
    """Analytic barycentric moments; verifier instead uses simplex quadrature."""
    tensor = np.zeros((3, 3))
    total = 0.0
    for row in local:
        ec = np.einsum("e,eic->ic", e[row["edge_indices"]], row["one_coefficients"])
        bc = np.einsum("f,fic->ic", b[row["face_indices"]], row["two_coefficients"])
        moments = row["volume"] * (np.ones((4, 4)) + np.eye(4)) / 20
        quadratic = ec.T @ moments @ ec + bc.T @ moments @ bc
        density_integral = float(np.trace(quadratic) / 2)
        tensor += density_integral * np.eye(3) - quadratic
        total += density_integral
    return tensor, total


def readout(a, pi, m, n, c, local, volume):
    e, b = -np.linalg.solve(m, pi), c @ a
    total, ue, ub = energy(a, pi, m, n, c)
    tensor, field_total = stress_integral(e, b, local)
    pressure = tensor / volume
    mean = float(np.trace(pressure) / 3)
    return {"energy": total, "electric_energy": ue, "magnetic_energy": ub,
            "field_integral_energy": field_total, "rho": total / volume,
            "stress_tensor": pressure.tolist(), "mean_pressure": mean,
            "w_mean_stress": mean * volume / total if total > 0 else None,
            "stress_eigenvalues": np.linalg.eigvalsh(pressure).tolist(),
            "anisotropic_stress_frobenius": float(np.linalg.norm(pressure - mean * np.eye(3))),
            "E": e.tolist(), "B": b.tolist()}


def metric_probe(a, pi, mesh, m, n, volume):
    """Volume work from freshly assembled deformed forms, without an EoS formula."""
    c = np.array(mesh["C"])
    total, ue, ub = energy(a, pi, m, n, c)
    isotropic = []
    for eps in EPSILONS:
        values = []
        for sign in (-1, 1):
            mm, nn, _, vv = forms(mesh, (1 + sign * eps) * np.eye(3))
            values.append({"length_factor": 1 + sign * eps, "volume": vv,
                           "energy": energy(a, pi, mm, nn, c)[0]})
        pressure = -(values[1]["energy"] - values[0]["energy"]) / (values[1]["volume"] - values[0]["volume"])
        isotropic.append({"epsilon": eps, "minus": values[0], "plus": values[1],
                          "work_pressure": pressure,
                          "work_w": pressure * volume / total if total > 0 else None})
    # Matrix Frechet derivative by symmetric differences; shear pairs use 1/2
    # on each off-diagonal, so -dH/ds is the corresponding stress component.
    directional = []
    eps = EPSILONS[-1]
    for i, j in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
        g = np.zeros((3, 3))
        g[i, j] = g[j, i] = 1.0 if i == j else 0.5
        vals = []
        for sign in (-1, 1):
            mm, nn, _, _ = forms(mesh, np.eye(3) + sign * eps * g)
            vals.append(energy(a, pi, mm, nn, c)[0])
        directional.append({"component": [i, j], "minus_energy": vals[0], "plus_energy": vals[1],
                            "pressure": -(vals[1] - vals[0]) / (2 * eps * volume)})
    return {"isotropic": isotropic, "directional": directional}


def initial_data(mesh):
    x, c = np.array(mesh["vertices"]), np.array(mesh["C"])
    # Exact line integrals of A=(B cross x)/2, B=(0,0,1); not fitted to w.
    b0 = np.array([0.0, 0.0, 1.0])
    a = np.array([np.cross(b0, (x[u] + x[v]) / 2) @ (x[v] - x[u]) / 2
                  for u, v in mesh["edges"]])
    # Incidence boundary-of-boundary gives D^T pi=0, without projection.
    face_seed = np.array([((7 * f + 3) % 17 - 8) / 100 for f in range(50)])
    pi = c.T @ face_seed
    return a, pi, face_seed


def source_free_episode(mesh, base):
    m, n, local, volume = base
    c, d = np.array(mesh["C"]), np.array(mesh["D"])
    k, invm = c.T @ n @ c, np.linalg.inv(m)
    a, pi, seed = initial_data(mesh)
    hessian_flow = np.block([[np.zeros_like(m), invm], [-k, np.zeros_like(m)]])
    left, right = np.eye(84) - DT * hessian_flow / 2, np.eye(84) + DT * hessian_flow / 2
    factor = lu_factor(left)
    state = np.concatenate((a, pi))
    frames, checks = [], []
    for step in range(STEPS + 1):
        a, pi = state[:42], state[42:]
        frame = {"step": step, "model_time": step * DT, "A": a.tolist(), "pi": pi.tolist(),
                 **readout(a, pi, m, n, c, local, volume)}
        frames.append(frame)
        if step == STEPS:
            break
        next_state = lu_solve(factor, right @ state)
        aa, pp = next_state[:42], next_state[42:]
        amid, pmid = (a + aa) / 2, (pi + pp) / 2
        emid, bmid = -invm @ pmid, c @ amid
        disp = -(pp - pi) / DT
        curl = c.T @ n @ bmid
        checks.append({"step": step,
            "hamilton_position_residual": float(np.max(np.abs((aa - a) / DT + emid))),
            "ampere_residual": float(np.max(np.abs(disp - curl))),
            "faraday_residual": float(np.max(np.abs((c @ aa - c @ a) / DT + c @ emid))),
            "gauss_load_max_abs": float(np.max(np.abs(d.T @ pp))),
            "displacement_current_l2": float(np.linalg.norm(disp)),
            "conduction_current_l2": 0.0,
            "energy_change": energy(aa, pp, m, n, c)[0] - frame["energy"],
            "external_work": 0.0})
        state = next_state
    return {"initial_data_rule": "A is line integral of (ez cross x)/2; pi=C^T*s; s[f]=((7*f+3)%17-8)/100",
            "face_seed": seed.tolist(), "step_size": DT, "steps": STEPS,
            "integrator": "implicit midpoint in canonical (A,pi); temporal gauge; J=rho_load=0",
            "readout_custody": "new numerical episode; complete states retained; not a new self-reading observer instrument",
            "frames": frames, "step_checks": checks,
            "metric_probes": [{"step": step, **metric_probe(np.array(frames[step]["A"]),
                                np.array(frames[step]["pi"]), mesh, m, n, volume)} for step in (0, 40, 80)]}


def inherited_observer_readout(mesh, base):
    """Use authenticated decoded A,phi; do not import copied field tables."""
    parent = source_replay.load()
    verification = source_replay.verify(parent)
    m, n, local, volume = base
    c, d = np.array(mesh["C"]), np.array(mesh["D"])
    rows = []
    for ex in parent["executions"]:
        decoded = np.array([[float(Fraction(x)) for x in frame] for frame in ex["decoded"]])
        for slab in (0, 1):
            phi = decoded[slab, :13]
            a = (decoded[slab, 13:] + decoded[slab + 1, 13:]) / 2
            e = -(decoded[slab + 1, 13:] - decoded[slab, 13:]) / 0.5 - d @ phi
            pi = -m @ e
            rows.append({"gauge": ex["gauge"], "slab": slab,
                "decode_event_ids": ex["decode_event_ids"][slab:slab + 2],
                "time_location": "slab midpoint of the declared piecewise linear A, constant phi history",
                "A": a.tolist(), "pi": pi.tolist(), **readout(a, pi, m, n, c, local, volume)})
    return {"parent_verification": verification,
            "scope": "field-only stress from authenticated original sourced readout; excludes source, clock, memory and wall energies",
            "energy_boundary": "instantaneous continuum-time Whitney field energy; not the parent's modified conserved integrator energy",
            "rows": rows}


def controls(mesh, base, episode):
    m, n, local, volume = base
    c, d = np.array(mesh["C"]), np.array(mesh["D"])
    first = episode["frames"][0]
    a, pi = np.array(first["A"]), np.array(first["pi"])
    chi = np.array([(u * u + 3 * u) % 11 - 5 for u in range(13)]) / 7
    shifted = a + d @ chi
    pure_magnetic = readout(a, np.zeros(42), m, n, c, local, volume)
    e, b = -np.linalg.solve(m, pi), c @ a
    eps = EPSILONS[-1]
    wrong, vv = [], []
    for sign in (-1, 1):
        mm, nn, _, v = forms(mesh, (1 + sign * eps) * np.eye(3))
        wrong.append(float((e @ mm @ e + b @ nn @ b) / 2))
        vv.append(v)
    wrong_pressure = -(wrong[1] - wrong[0]) / (vv[1] - vv[0])
    field_energy = first["energy"]
    # Explicit comparison Hamiltonians expose missing energy-sector closure.
    offset = field_energy
    vacuum_density = field_energy / volume
    return {"gauge": {"chi": chi.tolist(), "A_shifted": shifted.tolist(),
                      "energy_difference": energy(shifted, pi, m, n, c)[0] - field_energy},
            "zero_field": readout(np.zeros(42), np.zeros(42), m, n, c, local, volume),
            "anisotropic_pure_magnetic": pure_magnetic,
            "wrong_fixed_electric_cochain": {"description": "invalid canonical virtual work: E held fixed instead of pi",
                "work_pressure": wrong_pressure, "work_w": wrong_pressure * volume / field_energy,
                "minus_energy": wrong[0], "plus_energy": wrong[1]},
            "frozen_metric": {"description": "if the abstract update ignores geometry, dH/dV is zero by convention",
                "work_pressure": 0.0, "field_stress_mean": first["mean_pressure"]},
            "added_constant_energy": {"description": "counterexample Hamiltonian H+E0; not a claimed OPH sector",
                "E0": offset, "rho_total": (field_energy + offset) / volume,
                "pressure_total": first["mean_pressure"],
                "w_total": first["mean_pressure"] * volume / (field_energy + offset)},
            "added_volume_energy": {"description": "counterexample Hamiltonian H+rho_v*V; not a claimed OPH sector",
                "rho_v": vacuum_density, "rho_total": field_energy / volume + vacuum_density,
                "pressure_total": first["mean_pressure"] - vacuum_density,
                "w_total": (first["mean_pressure"] - vacuum_density) / (field_energy / volume + vacuum_density)}}


def classical_gibbs(mesh):
    """Finite classical Gibbs ensemble on a fixed canonical gauge quotient.

    Q is a Euclidean orthonormal complement of gradients, chosen from the
    topology once, independent of L. Set A=Q*q and pi=Q*p: pi.dA=p.dq.
    The Hamiltonian forms are G=Q^T M^-1 Q and Kq=Q^T C^T N C Q.
    The Gaussian integral gives log Z = n log T - 1/2 log(det G det Kq),
    up to a volume-independent phase-space cell normalization. No
    L-dependent eigenvector Jacobian is discarded. Temperature is supplied;
    Hamiltonian execution alone neither prepares Gibbs nor thermalizes.
    """
    d, c = np.array(mesh["D"]), np.array(mesh["C"])
    qbasis = null_space(d.T)
    require(qbasis.shape == (42, 30), "canonical quotient dimension")
    cases = []
    for length in (0.5, 1.0, 2.0):
        m, n, local, volume = forms(mesh, length * np.eye(3))
        invm = np.linalg.inv(m)
        g = qbasis.T @ invm @ qbasis
        k = qbasis.T @ c.T @ n @ c @ qbasis
        freq = np.sqrt(eigh(k, np.linalg.inv(g), eigvals_only=True))
        covariance_a_unit = qbasis @ np.linalg.solve(k, qbasis.T)
        covariance_pi_unit = qbasis @ np.linalg.solve(g, qbasis.T)
        require(np.min(freq) > 0.1, "positive physical Gibbs frequencies")
        deformation_data = []
        for sign in (-1, 1):
            factor = length * (1 + sign * 0.0001)
            mm, nn, _, vv = forms(mesh, factor * np.eye(3))
            gg = qbasis.T @ np.linalg.solve(mm, qbasis)
            kk = qbasis.T @ c.T @ nn @ c @ qbasis
            logdet = np.linalg.slogdet(gg)[1] + np.linalg.slogdet(kk)[1]
            deformation_data.append({"length": factor, "volume": vv,
                                     "log_det_product": float(logdet)})
        rows = []
        for temperature in (0.25, 1.0, 4.0):
            covariance_a = temperature * covariance_a_unit
            covariance_pi = temperature * covariance_pi_unit
            ecov = invm @ covariance_pi @ invm
            bcov = c @ covariance_a @ c.T
            stress, total_integral = np.zeros((3, 3)), 0.0
            for row in local:
                ec = np.array(row["one_coefficients"])
                bc = np.array(row["two_coefficients"])
                moment = row["volume"] * (np.ones((4, 4)) + np.eye(4)) / 20
                local_e = ecov[np.ix_(row["edge_indices"], row["edge_indices"])]
                local_b = bcov[np.ix_(row["face_indices"], row["face_indices"])]
                qq = np.einsum("ef,eia,ij,fjb->ab", local_e, ec, moment, ec)
                qq += np.einsum("ef,eia,ij,fjb->ab", local_b, bc, moment, bc)
                uu = float(np.trace(qq) / 2)
                stress += uu * np.eye(3) - qq
                total_integral += uu
            stress /= volume
            ue = float(np.trace(invm @ covariance_pi) / 2)
            ub = float(np.trace(c.T @ n @ c @ covariance_a) / 2)
            partition = [30 * np.log(temperature) - row["log_det_product"] / 2 for row in deformation_data]
            work_pressure = temperature * (partition[1] - partition[0]) / (deformation_data[1]["volume"] - deformation_data[0]["volume"])
            mean_pressure = float(np.trace(stress) / 3)
            rows.append({"temperature": temperature, "covariance_A": covariance_a.tolist(),
                "covariance_pi": covariance_pi.tolist(), "electric_energy": ue, "magnetic_energy": ub,
                "energy": ue + ub, "field_integral_energy": total_integral,
                "rho": (ue + ub) / volume, "stress_tensor": stress.tolist(),
                "mean_pressure": mean_pressure, "w": mean_pressure * volume / (ue + ub),
                "anisotropic_stress_frobenius": float(np.linalg.norm(stress - mean_pressure * np.eye(3))),
                "log_Z_without_volume_independent_cell_constant": float(30 * np.log(temperature) - np.log(freq).sum()),
                "partition_minus": float(partition[0]), "partition_plus": float(partition[1]),
                "partition_work_pressure": float(work_pressure)})
        cases.append({"length": length, "volume": volume, "frequencies": freq.tolist(),
                      "partition_metric_stencil": deformation_data, "rows": rows})
    return {"scope": "SUPPLIED_CLASSICAL_GIBBS_ENSEMBLE_OF_FINITE_MAXWELL_GAUGE_QUOTIENT",
        "ensemble_prepared_by_native_repair": False,
        "thermalization_demonstrated": False,
        "quantum_or_continuum_blackbody_law": False,
        "canonical_measure": "fixed topology-only gauge quotient; Liouville dq dp; cell normalization independent of geometry",
        "temperature_units": "supplied classical temperature with k_B=1; no SI calibration",
        "physical_modes": 30, "zero_gauge_modes_excluded": 12,
        "scale_cases": cases}


def build():
    mesh, base = mesh_and_forms()
    m, n, local, volume = base
    episode = source_free_episode(mesh, base)
    scale_checks = []
    for length in (0.5, 1.5, 2.0):
        mm, nn, _, v = forms(mesh, length * np.eye(3))
        scale_checks.append({"length": length, "volume": v,
            "mass_one_scaling_error": float(np.max(np.abs(mm - length * m))),
            "mass_two_scaling_error": float(np.max(np.abs(nn - n / length)))})
    return {"schema": SCHEMA, "scope": SCOPE,
        "contract": {"native_repair_eos": False,
            "equilibrium_thermodynamics": "separate supplied classical Gibbs branch only",
            "geometry": "supplied three-dimensional Euclidean tetrahedral cone over twelve ports",
            "action": "supplied unit-coefficient Whitney Maxwell action",
            "pressure": "canonical metric virtual work and averaged field momentum flux",
            "phase_space": "edge-integrated A and conjugate pi fixed under metric deformation",
            "energy_scope": "Maxwell field only; excludes matter, source, clock, wall and memory energies",
            "boundary": "closed finite natural variational boundary; no explicit exterior flux channel",
            "clock": "supplied model time; no physical-unit calibration",
            "precision": "float64 with independent quadrature and numeric tolerances; not interval arithmetic"},
        "source_pins": {p: sha(RER / p) for p in SOURCE_PATHS},
        "implementation_pins": {p: sha(HERE / p) for p in ("maxwell_eos.py", "verify_maxwell_eos.py", "test_maxwell_eos.py")},
        "mesh": mesh, "forms": {"M1": m.tolist(), "M2": n.tolist(), "volume": volume,
            "rank_D": int(np.linalg.matrix_rank(mesh["D"])), "rank_C": int(np.linalg.matrix_rank(mesh["C"]))},
        "scaling_reassembly": scale_checks,
        "source_free": episode,
        "authenticated_observer_fields": inherited_observer_readout(mesh, base),
        "controls": controls(mesh, base, episode),
        "classical_gibbs": classical_gibbs(mesh)}


def write_exports(packet, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(canonical(packet))
    columns = ["step", "model_time", "energy", "electric_energy", "magnetic_energy", "rho", "mean_pressure", "w_mean_stress", "anisotropic_stress_frobenius"]
    with output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: frame[key] for key in columns} for frame in packet["source_free"]["frames"])
    thermal_columns = ["length", "volume", "temperature", "energy", "electric_energy",
        "magnetic_energy", "rho", "mean_pressure", "w", "partition_work_pressure", "anisotropic_stress_frobenius"]
    with output.with_name("maxwell_classical_gibbs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=thermal_columns)
        writer.writeheader()
        for case in packet["classical_gibbs"]["scale_cases"]:
            for row in case["rows"]:
                writer.writerow({key: case[key] if key in ("length", "volume") else row[key] for key in thermal_columns})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "runs/maxwell_eos_receipt.json")
    args = parser.parse_args()
    receipt = build()
    write_exports(receipt, args.output)
    first = receipt["source_free"]["frames"][0]
    print(json.dumps({"receipt": str(args.output), "energy": first["energy"], "rho": first["rho"],
        "mean_pressure": first["mean_pressure"], "w": first["w_mean_stress"], "scope": SCOPE}, sort_keys=True))

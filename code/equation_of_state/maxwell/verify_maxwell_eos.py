"""Independent quadrature, variational-work and trajectory audit.

Does not import the EoS producer. Reconstructs the committed geometry using
the existing independent cone verifier, evaluates Whitney fields at positive
degree-two quadrature nodes, and replays Hamilton's equations directly.
Source observer states must come from the authenticated decode operations.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
RER = HERE.parents[2]
sys.path.insert(0, str(RER / "code/electromagnetism"))
import verify_cone_whitney_bridge as geometry
import verify_whitney_maxwell_dynamics as original

SCOPE = "SUPPLIED_MAXWELL_ACTION__CANONICAL_METRIC_WORK__NO_NATIVE_REPAIR_EOS"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def close(actual, expected, name, atol=1e-9, rtol=1e-9):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
            and np.allclose(a, b, atol=atol, rtol=rtol), name)


def load(path):
    def reject_duplicates(pairs):
        result = {}
        for k, v in pairs:
            require(k not in result, "duplicate JSON key")
            result[k] = v
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=reject_duplicates,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))


def build_mesh():
    v, e, f = geometry.source_mesh()
    mesh = {"vertices": [[0.0] * 3] + [[float(x) for x in row] for row in v],
            "edges": [[0, u + 1] for u in range(12)] + [[u + 1, w + 1] for u, w in e],
            "faces": [[u + 1 for u in face] for face in f] + [[0, u + 1, w + 1] for u, w in e],
            "tetrahedra": [[0, *(u + 1 for u in face)] for face in f]}
    mesh["D"] = geometry.coboundary([[u] for u in range(13)], mesh["edges"]).tolist()
    mesh["C"] = geometry.coboundary(mesh["edges"], mesh["faces"]).tolist()
    return mesh


def quadrature(mesh, deformation=None):
    deformation = np.eye(3) if deformation is None else deformation
    vertices = np.array(mesh["vertices"]) @ deformation.T
    return geometry.quadrature(vertices, mesh["edges"], mesh["faces"], mesh["tetrahedra"])


def hamiltonian(a, pi, c, quad):
    _, _, _, _, m, n = quad
    b = c @ a
    return float((pi @ np.linalg.solve(m, pi) + b @ n @ b) / 2)


def check_readout(row, a, pi, c, quad):
    one, two, weights, volumes, m, n = quad
    e = -np.linalg.solve(m, pi)
    b = c @ a
    ef = np.einsum("e,qec->qc", e, one)
    bf = np.einsum("f,qfc->qc", b, two)
    uu_e, uu_b = np.sum(ef * ef, axis=1) / 2, np.sum(bf * bf, axis=1) / 2
    volume = float(volumes.sum())
    ue, ub = float(weights @ uu_e), float(weights @ uu_b)
    total = ue + ub
    tensor = total * np.eye(3) - np.einsum("q,qi,qj->ij", weights, ef, ef) - np.einsum("q,qi,qj->ij", weights, bf, bf)
    tensor /= volume
    mean = float(np.trace(tensor) / 3)
    expected = {"E": e, "B": b, "energy": total, "electric_energy": ue,
                "magnetic_energy": ub, "field_integral_energy": total,
                "rho": total / volume, "stress_tensor": tensor, "mean_pressure": mean,
                "stress_eigenvalues": np.linalg.eigvalsh(tensor),
                "anisotropic_stress_frobenius": np.linalg.norm(tensor - mean * np.eye(3))}
    for k, v in expected.items():
        close(row[k], v, "field quadrature " + k)
    if total == 0:
        require(row["w_mean_stress"] is None, "zero energy has undefined ratio")
    else:
        close(row["w_mean_stress"], mean * volume / total, "ratio from integrated stress")
    return expected


def check_probe(probe, a, pi, mesh, base):
    c = np.array(mesh["C"])
    volume = float(base[3].sum())
    h0 = hamiltonian(a, pi, c, base)
    require(len(probe["isotropic"]) == 3, "isotropic stencil census")
    errors = []
    for eps, row in zip((0.01, 0.001, 0.0001), probe["isotropic"], strict=True):
        require(row["epsilon"] == eps, "fixed metric perturbations")
        data = []
        for sign, side in ((-1, "minus"), (1, "plus")):
            length = 1 + sign * eps
            qq = quadrature(mesh, length * np.eye(3))
            vv, hh = float(qq[3].sum()), hamiltonian(a, pi, c, qq)
            close(row[side]["length_factor"], length, "dilation factor")
            close(row[side]["volume"], vv, "reassembled dilation volume")
            close(row[side]["energy"], hh, "reassembled canonical energy")
            data.append((vv, hh))
        p = -(data[1][1] - data[0][1]) / (data[1][0] - data[0][0])
        close(row["work_pressure"], p, "metric-work pressure", atol=3e-9)
        close(row["work_w"], p * volume / h0, "metric-work ratio", atol=3e-9)
        errors.append(abs(p - h0 / (3 * volume)))
    require(errors[-1] < 1e-7 * h0 / volume, "metric-work / stress-trace agreement")
    require(errors[2] < errors[1] / 20 and errors[1] < errors[0] / 20, "dilation stencil convergence")
    e, b = -np.linalg.solve(base[4], pi), c @ a
    ef, bf = np.einsum("e,qec->qc", e, base[0]), np.einsum("f,qfc->qc", b, base[1])
    stress = (h0 * np.eye(3) - np.einsum("q,qi,qj->ij", base[2], ef, ef)
              - np.einsum("q,qi,qj->ij", base[2], bf, bf)) / volume
    require(len(probe["directional"]) == 6, "directional census")
    for (i, j), row in zip(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)), probe["directional"], strict=True):
        require(row["component"] == [i, j], "directional component")
        gg = np.zeros((3, 3))
        gg[i, j] = gg[j, i] = 1 if i == j else 0.5
        hm = hamiltonian(a, pi, c, quadrature(mesh, np.eye(3) - 0.0001 * gg))
        hp = hamiltonian(a, pi, c, quadrature(mesh, np.eye(3) + 0.0001 * gg))
        close(row["minus_energy"], hm, "directional minus energy")
        close(row["plus_energy"], hp, "directional plus energy")
        pressure = -(hp - hm) / (0.0002 * volume)
        close(row["pressure"], pressure, "directional work", atol=3e-9)
        close(row["pressure"], stress[i, j], "directional Maxwell stress", atol=2e-7 * h0 / volume)


def check_observer(rows, original_packet, c, d, base):
    require(len(rows) == 4, "observer readout census")
    m = base[4]
    for gi, ex in enumerate(original_packet["executions"]):
        # Read actual writes, not public copied field/decoded summaries.
        decoded = np.zeros((3, 55))
        ids = []
        for event in ex["events"]:
            if event["op"] == "decode":
                ids.append(event["id"])
                for key, value in event["writes"].items():
                    if key.startswith("d/"):
                        _, nn, ss = key.split("/")
                        decoded[int(nn), int(ss)] = float(Fraction(value))
        require(len(ids) == 3, "original decode census")
        for slab in range(2):
            row = rows[gi * 2 + slab]
            require(row["gauge"] == bool(gi) and row["slab"] == slab, "observer state identity")
            require(row["decode_event_ids"] == ids[slab:slab + 2], "actual decode provenance")
            a = (decoded[slab, 13:] + decoded[slab + 1, 13:]) / 2
            e = -2 * (decoded[slab + 1, 13:] - decoded[slab, 13:]) - d @ decoded[slab, :13]
            pi = -m @ e
            close(row["A"], a, "decoded midpoint potential")
            close(row["pi"], pi, "decoded canonical field momentum")
            check_readout(row, a, pi, c, base)
    for slab in range(2):
        for key in ("E", "B", "energy", "stress_tensor", "w_mean_stress"):
            close(rows[slab][key], rows[slab + 2][key], "observer gauge invariance " + key)


def check_controls(ctrl, mesh, base, frame):
    c, d, volume = np.array(mesh["C"]), np.array(mesh["D"]), float(base[3].sum())
    a, pi = np.array(frame["A"]), np.array(frame["pi"])
    chi = np.array([(u * u + 3 * u) % 11 - 5 for u in range(13)]) / 7
    close(ctrl["gauge"]["chi"], chi, "gauge control parameter")
    close(ctrl["gauge"]["A_shifted"], a + d @ chi, "gauge action")
    require(np.linalg.norm(d @ chi) > 1, "nontrivial gauge perturbation")
    difference = hamiltonian(a + d @ chi, pi, c, base) - frame["energy"]
    close(ctrl["gauge"]["energy_difference"], difference, "gauge energy")
    close(difference, 0, "gauge invariant Hamiltonian")
    check_readout(ctrl["zero_field"], np.zeros(42), np.zeros(42), c, base)
    pure = check_readout(ctrl["anisotropic_pure_magnetic"], a, np.zeros(42), c, base)
    require(pure["stress_eigenvalues"][0] < -0.9 * pure["rho"], "directional tension control")
    require(pure["stress_eigenvalues"][-1] > 0.9 * pure["rho"], "directional positive pressure control")
    wrong = ctrl["wrong_fixed_electric_cochain"]
    require(wrong["description"].startswith("invalid canonical virtual work"), "wrong-variable control interpretation")
    e, b = np.array(frame["E"]), np.array(frame["B"])
    vals, volumes = [], []
    for sign in (-1, 1):
        qq = quadrature(mesh, (1 + sign * 0.0001) * np.eye(3))
        vals.append(float((e @ qq[4] @ e + b @ qq[5] @ b) / 2))
        volumes.append(float(qq[3].sum()))
    pp = -(vals[1] - vals[0]) / (volumes[1] - volumes[0])
    close(wrong["minus_energy"], vals[0], "wrong ensemble minus")
    close(wrong["plus_energy"], vals[1], "wrong ensemble plus")
    close(wrong["work_pressure"], pp, "wrong ensemble work", atol=3e-9)
    close(wrong["work_w"], pp * volume / frame["energy"], "wrong ensemble ratio", atol=3e-9)
    require(abs(wrong["work_w"] - frame["w_mean_stress"]) > 0.01, "wrong ensemble must separate")
    close(ctrl["frozen_metric"]["work_pressure"], 0, "frozen abstract metric work")
    close(ctrl["frozen_metric"]["field_stress_mean"], frame["mean_pressure"], "frozen-metric mismatch")
    extra = ctrl["added_constant_energy"]
    require(extra["description"].startswith("counterexample Hamiltonian"), "constant-energy scope")
    close(extra["E0"], frame["energy"], "added constant input")
    close(extra["rho_total"], 2 * frame["rho"], "constant-energy total density")
    close(extra["pressure_total"], frame["mean_pressure"], "constant-energy total stress")
    close(extra["w_total"], extra["pressure_total"] / extra["rho_total"], "constant-energy ratio")
    extra = ctrl["added_volume_energy"]
    require(extra["description"].startswith("counterexample Hamiltonian"), "volume-energy scope")
    close(extra["rho_v"], frame["rho"], "added volume-energy input")
    close(extra["rho_total"], 2 * frame["rho"], "volume-energy total density")
    close(extra["pressure_total"], frame["mean_pressure"] - frame["rho"], "volume-energy work")
    close(extra["w_total"], extra["pressure_total"] / extra["rho_total"], "volume-energy ratio")


def check_gibbs(gibbs, mesh):
    """Independent tree/cotree gauge; no shared orthonormal quotient routine.

    Radial edges form a spanning tree. Fix their potentials to zero: A=R*q.
    The Gauss constraint gives pi=S*p, with S=(-D_boundary^T,I)^T and R^T S=I.
    This reduction is canonical and independent of geometry. It therefore
    uses the same Liouville measure up to a volume-independent constant.
    """
    require(gibbs["scope"] == "SUPPLIED_CLASSICAL_GIBBS_ENSEMBLE_OF_FINITE_MAXWELL_GAUGE_QUOTIENT", "Gibbs scope")
    for flag in ("ensemble_prepared_by_native_repair", "thermalization_demonstrated", "quantum_or_continuum_blackbody_law"):
        require(gibbs[flag] is False, "Gibbs nonclaim " + flag)
    require(gibbs["canonical_measure"] == "fixed topology-only gauge quotient; Liouville dq dp; cell normalization independent of geometry", "fixed Liouville measure")
    require(gibbs["temperature_units"] == "supplied classical temperature with k_B=1; no SI calibration", "temperature assumption")
    require(gibbs["physical_modes"] == 30 and gibbs["zero_gauge_modes_excluded"] == 12, "physical Gibbs mode count")
    d, c = np.array(mesh["D"]), np.array(mesh["C"])
    r = np.vstack((np.zeros((12, 30)), np.eye(30)))
    s = np.vstack((-d[12:, 1:].T, np.eye(30)))
    close(r.T @ s, np.eye(30), "canonical tree-cotree pairing")
    close(d.T @ s, np.zeros((13, 30)), "tree-cotree Gauss constraint")
    grad = d[:, 1:]
    euclidean_projector = np.eye(42) - grad @ np.linalg.solve(grad.T @ grad, grad.T)
    require(len(gibbs["scale_cases"]) == 3, "Gibbs scale census")
    maximum_anisotropy = 0.0
    for length, case in zip((0.5, 1.0, 2.0), gibbs["scale_cases"], strict=True):
        require(case["length"] == length, "Gibbs length")
        quad = quadrature(mesh, length * np.eye(3))
        one, two, weights, volumes, m, n = quad
        volume = float(volumes.sum())
        close(case["volume"], volume, "Gibbs volume")
        g = s.T @ np.linalg.solve(m, s)
        k = r.T @ c.T @ n @ c @ r
        require(np.linalg.eigvalsh(g).min() > 0 and np.linalg.eigvalsh(k).min() > 0, "normalizable reduced Gibbs form")
        eig = np.linalg.eigvals(g @ k)
        require(np.max(np.abs(eig.imag)) < 1e-9 and np.min(eig.real) > 0, "physical Hamiltonian frequencies")
        freq = np.sort(np.sqrt(eig.real))
        close(case["frequencies"], freq, "independent gauge-quotient spectrum")
        ga = euclidean_projector @ r
        covariance_a_unit = ga @ np.linalg.solve(k, ga.T)
        covariance_pi_unit = s @ np.linalg.solve(g, s.T)
        logdet = float(np.linalg.slogdet(g)[1] + np.linalg.slogdet(k)[1])
        stencil = []
        require(len(case["partition_metric_stencil"]) == 2, "partition stencil census")
        for sign, row in zip((-1, 1), case["partition_metric_stencil"], strict=True):
            ll = length * (1 + sign * 0.0001)
            qq = quadrature(mesh, ll * np.eye(3))
            gg = s.T @ np.linalg.solve(qq[4], s)
            kk = r.T @ c.T @ qq[5] @ c @ r
            determinant = float(np.linalg.slogdet(gg)[1] + np.linalg.slogdet(kk)[1])
            close(row["length"], ll, "partition dilation")
            close(row["volume"], qq[3].sum(), "partition deformed volume")
            close(row["log_det_product"], determinant, "canonical Gaussian determinant")
            stencil.append((float(qq[3].sum()), determinant))
        require(len(case["rows"]) == 3, "temperature census")
        for temperature, row in zip((0.25, 1.0, 4.0), case["rows"], strict=True):
            require(row["temperature"] == temperature, "prepared temperature")
            ca, cp = temperature * covariance_a_unit, temperature * covariance_pi_unit
            close(row["covariance_A"], ca, "gauge-fixed thermal potential covariance")
            close(row["covariance_pi"], cp, "thermal momentum covariance")
            ec = np.linalg.solve(m, np.linalg.solve(m, cp).T).T
            bc = c @ ca @ c.T
            ee = np.einsum("ef,qei,qfj->qij", ec, one, one)
            bb = np.einsum("ef,qei,qfj->qij", bc, two, two)
            ue, ub = float(np.einsum("q,qii->", weights, ee) / 2), float(np.einsum("q,qii->", weights, bb) / 2)
            total = ue + ub
            stress = (total * np.eye(3) - np.einsum("q,qij->ij", weights, ee + bb)) / volume
            pressure = float(np.trace(stress) / 3)
            anisotropy = float(np.linalg.norm(stress - pressure * np.eye(3)))
            partitions = [30 * np.log(temperature) - v[1] / 2 for v in stencil]
            work_pressure = temperature * (partitions[1] - partitions[0]) / (stencil[1][0] - stencil[0][0])
            expected = {"electric_energy": ue, "magnetic_energy": ub, "energy": total,
                "field_integral_energy": total, "rho": total / volume, "stress_tensor": stress,
                "mean_pressure": pressure, "w": pressure * volume / total,
                "anisotropic_stress_frobenius": anisotropy,
                "log_Z_without_volume_independent_cell_constant": 30 * np.log(temperature) - logdet / 2,
                "partition_minus": partitions[0], "partition_plus": partitions[1],
                "partition_work_pressure": work_pressure}
            for name, value in expected.items():
                close(row[name], value, "Gibbs " + name, atol=1e-7, rtol=1e-8)
            close(total, 30 * temperature, "finite classical equipartition")
            close(work_pressure, pressure, "partition work / thermal stress", atol=1e-7, rtol=1e-8)
            require(anisotropy < 1e-8 * pressure, "thermal stress isotropic on full icosahedral cone")
            maximum_anisotropy = max(maximum_anisotropy, anisotropy)
    return maximum_anisotropy


def verify(packet, replay_original=True):
    require(packet["schema"] == "oph.exploratory_whitney_maxwell_eos.v1", "schema")
    require(packet["scope"] == SCOPE, "interpretation scope")
    contract = packet["contract"]
    require(contract == {"native_repair_eos": False,
        "equilibrium_thermodynamics": "separate supplied classical Gibbs branch only",
        "geometry": "supplied three-dimensional Euclidean tetrahedral cone over twelve ports",
        "action": "supplied unit-coefficient Whitney Maxwell action",
        "pressure": "canonical metric virtual work and averaged field momentum flux",
        "phase_space": "edge-integrated A and conjugate pi fixed under metric deformation",
        "energy_scope": "Maxwell field only; excludes matter, source, clock, wall and memory energies",
        "boundary": "closed finite natural variational boundary; no explicit exterior flux channel",
        "clock": "supplied model time; no physical-unit calibration",
        "precision": "float64 with independent quadrature and numeric tolerances; not interval arithmetic"}, "scope and consumed assumptions")
    required_pins = {"Lean/Screen/SeamCurrentCarrierQuotient.lean", "Lean/ObserverPatchHolography/CoreAxioms.lean",
        "Lean/Screen/SeamCurrentEdge30Moment.lean", "code/electromagnetism/cone_whitney_bridge.py",
        "code/electromagnetism/verify_cone_whitney_bridge.py", "code/electromagnetism/verify_whitney_maxwell_dynamics.py",
        "code/electromagnetism/runtime/whitney_maxwell_dynamics_receipt.json"}
    require(set(packet["source_pins"]) == required_pins, "source pin census")
    require(set(packet["implementation_pins"]) == {"maxwell_eos.py", "verify_maxwell_eos.py", "test_maxwell_eos.py"}, "implementation pin census")
    for root, pins in ((RER, packet["source_pins"]), (HERE, packet["implementation_pins"])):
        for path, digest in pins.items():
            require(hashlib.sha256((root / path).read_bytes()).hexdigest() == digest, "source digest " + path)
    mesh = build_mesh()
    require(packet["mesh"] == mesh, "committed oriented mesh")
    base = quadrature(mesh)
    c, d = np.array(mesh["C"]), np.array(mesh["D"])
    m, n, volume = base[4], base[5], float(base[3].sum())
    require(np.array_equal(c @ d, np.zeros((50, 13))), "boundary-of-boundary identity")
    close(packet["forms"]["M1"], m, "independent electric form")
    close(packet["forms"]["M2"], n, "independent magnetic form")
    close(packet["forms"]["volume"], volume, "independent geometric volume")
    require(packet["forms"]["rank_D"] == np.linalg.matrix_rank(d) == 12, "gauge rank")
    require(packet["forms"]["rank_C"] == np.linalg.matrix_rank(c) == 30, "full-volume curvature rank")
    for length, row in zip((0.5, 1.5, 2), packet["scaling_reassembly"], strict=True):
        require(row["length"] == length, "scale contract")
        qq = quadrature(mesh, length * np.eye(3))
        close(row["volume"], qq[3].sum(), "scaled volume")
        close(row["mass_one_scaling_error"], np.max(np.abs(qq[4] - length * m)), "one-form dilation")
        close(row["mass_two_scaling_error"], np.max(np.abs(qq[5] - n / length)), "two-form dilation")
        require(max(row["mass_one_scaling_error"], row["mass_two_scaling_error"]) < 1e-10, "form homogeneity")
    episode = packet["source_free"]
    require(episode["steps"] == 80 and episode["step_size"] == 0.025, "run budget")
    require(episode["integrator"] == "implicit midpoint in canonical (A,pi); temporal gauge; J=rho_load=0", "integrator contract")
    require(episode["readout_custody"] == "new numerical episode; complete states retained; not a new self-reading observer instrument", "new episode custody")
    frames, checks = episode["frames"], episode["step_checks"]
    require(len(frames) == 81 and len(checks) == 80, "full trajectory retention")
    x = np.array(mesh["vertices"])
    # Different explicit equivalent line-integral formula for A=(-y/2,x/2,0).
    a0 = np.array([(x[u, 0] * x[v, 1] - x[u, 1] * x[v, 0]) / 2 for u, v in mesh["edges"]])
    seed = np.array([((7 * f + 3) % 17 - 8) / 100 for f in range(50)])
    close(episode["face_seed"], seed, "unfitted input seed")
    close(frames[0]["A"], a0, "magnetic initial data")
    close(frames[0]["pi"], c.T @ seed, "source-free initial momentum")
    require(np.linalg.norm(c.T @ seed) > 0.1, "nonzero electric preparation")
    h0 = hamiltonian(a0, c.T @ seed, c, base)
    max_displacement = 0.0
    for i, row in enumerate(frames):
        require(row["step"] == i and row["model_time"] == i * 0.025, "ordered model times")
        a, pi = np.array(row["A"]), np.array(row["pi"])
        require(a.shape == pi.shape == (42,), "full state dimension")
        check_readout(row, a, pi, c, base)
        close(row["energy"], h0, "closed-Hamiltonian energy balance")
        close(d.T @ pi, np.zeros(13), "zero charge constraint")
        if i == 80:
            continue
        aa, pp = np.array(frames[i + 1]["A"]), np.array(frames[i + 1]["pi"])
        emid = -np.linalg.solve(m, (pi + pp) / 2)
        bmid = c @ ((a + aa) / 2)
        displacement = -(pp - pi) / 0.025
        faraday = (c @ aa - c @ a) / 0.025 + c @ emid
        ampere = displacement - c.T @ n @ bmid
        position = (aa - a) / 0.025 + emid
        require(max(np.max(np.abs(position)), np.max(np.abs(ampere)), np.max(np.abs(faraday))) < 1e-9, "source-free Maxwell Hamilton equations")
        check = checks[i]
        require(check["step"] == i, "step check identity")
        expected = {"hamilton_position_residual": np.max(np.abs(position)), "ampere_residual": np.max(np.abs(ampere)),
            "faraday_residual": np.max(np.abs(faraday)), "gauss_load_max_abs": np.max(np.abs(d.T @ pp)),
            "displacement_current_l2": np.linalg.norm(displacement), "conduction_current_l2": 0.0,
            "energy_change": hamiltonian(aa, pp, c, base) - hamiltonian(a, pi, c, base), "external_work": 0.0}
        for key, value in expected.items():
            close(check[key], value, "step equation/work " + key)
        max_displacement = max(max_displacement, float(expected["displacement_current_l2"]))
    require(max_displacement > 0.1, "nonzero displacement with zero conduction current")
    require(len(episode["metric_probes"]) == 3, "metric probe census")
    for step, probe in zip((0, 40, 80), episode["metric_probes"], strict=True):
        require(probe["step"] == step, "metric state identity")
        check_probe(probe, np.array(frames[step]["A"]), np.array(frames[step]["pi"]), mesh, base)
    parent = original.load()
    auth = packet["authenticated_observer_fields"]
    if replay_original:
        fresh = original.verify(parent)
        recorded = auth["parent_verification"]
        require(set(recorded) == set(fresh), "original source replay summary census")
        for key, value in fresh.items():
            if key in ("gauss_max_abs", "ampere_max_abs"):
                # Round-off residuals of the parent replay differ by platform;
                # both the recorded and the fresh value must sit at the
                # float64 noise floor.
                require(0.0 <= float(recorded[key]) <= 1e-12 and 0.0 <= float(value) <= 1e-12,
                        "original source replay residual " + key)
            else:
                require(recorded[key] == value, "original source replay summary " + key)
    require(auth["scope"] == "field-only stress from authenticated original sourced readout; excludes source, clock, memory and wall energies", "sourced readout energy boundary")
    require(auth["energy_boundary"] == "instantaneous continuum-time Whitney field energy; not the parent's modified conserved integrator energy", "finite-step energy distinction")
    check_observer(auth["rows"], parent, c, d, base)
    check_controls(packet["controls"], mesh, base, frames[0])
    thermal_anisotropy = check_gibbs(packet["classical_gibbs"], mesh)
    return {"scope": SCOPE, "frames": 81, "source_free_steps": 80, "authenticated_sourced_field_rows": 4,
            "source_original_fully_replayed": replay_original,
            "volume": volume, "initial_energy": h0, "rho": frames[0]["rho"],
            "mean_pressure": frames[0]["mean_pressure"], "w_mean_stress": frames[0]["w_mean_stress"],
            "max_energy_drift": max(abs(row["energy"] - h0) for row in frames),
            "max_displacement_current_l2": max_displacement,
            "largest_metric_work_ratio_error": max(abs(p["isotropic"][-1]["work_w"] - 1 / 3) for p in episode["metric_probes"]),
            "classical_gibbs_cases": 9, "gibbs_max_anisotropic_stress": thermal_anisotropy,
            "independent_quadrature": True, "all_checks_passed": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", nargs="?", type=Path, default=HERE / "runs/maxwell_eos_receipt.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(verify(load(args.receipt)), sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result + "\n")
    print(result)

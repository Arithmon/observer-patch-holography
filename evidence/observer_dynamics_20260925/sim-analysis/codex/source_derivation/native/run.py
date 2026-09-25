"""Absolute finite-window source power from the pinned native repair rule.

General state-space readouts are supported without a density-only assumption.
The source calculation never imports a measured spectrum or fitted amplitude.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from scipy.integrate import quad_vec
from scipy.special import sph_harm_y

HERE = Path(__file__).resolve().parent
CODEX = HERE.parents[1]
WORKSPACE = HERE.parents[3]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def average_response(u):
    """h(u)=2*(u-1+exp(-u))/u², continuous at u=0; no CLT rescaling."""
    values = np.asarray(u, dtype=float)
    if np.any(values < -1e-12) or np.any(~np.isfinite(values)):
        raise ValueError("Nonnegative finite relaxation exponent required")
    values = np.maximum(values, 0.0)
    small = values < 1e-3
    out = np.empty_like(values)
    x = values[small]
    out[small] = 1-x/3+x*x/12-x**3/60+x**4/360-x**5/2520
    x = values[~small]
    out[~small] = 2*(x+np.expm1(-x))/(x*x)
    return out


def chain(ports, raised, edges, rule=None):
    if not 0 < raised < ports:
        raise ValueError("Require 0<raised<ports")
    states = np.array([tuple(int(i in s) for i in range(ports))
                       for s in itertools.combinations(range(ports), raised)])
    where = {tuple(row): i for i, row in enumerate(states)}
    q = np.zeros((len(states), len(states)))
    lap = np.zeros((ports, ports))
    for a, b, weight in edges:
        if a == b or min(a, b) < 0 or max(a, b) >= ports or weight <= 0 or not np.isfinite(weight):
            raise ValueError("Valid endpoints and positive finite edge weight required")
        lap[a, a] += weight; lap[b, b] += weight
        lap[a, b] -= weight; lap[b, a] -= weight
        for i, state in enumerate(states):
            for coin in (False, True):
                after = state.copy()
                if rule is None:
                    low, high = divmod(int(state[a]+state[b]), 2)
                    after[a], after[b] = (low+high, low) if coin else (low, low+high)
                else:
                    after[a], after[b] = rule(int(state[a]), int(state[b]), ceiling_to_first=coin)
                q[i, where[tuple(after)]] += weight/2
                q[i, i] -= weight/2
    if la.eigvalsh(lap)[1] <= 1e-12:
        raise ValueError("Connected graph required")
    return states, q, lap


def spectral_covariance(eigenvalues, eigenvectors, readout, window):
    """Stationary covariance of a literal finite average of any state readout."""
    if window <= 0 or not np.isfinite(window):
        raise ValueError("Positive finite averaging window required")
    centered = np.asarray(readout, dtype=float)
    if centered.ndim == 1:
        centered = centered[:, None]
    if np.max(np.abs(centered.mean(axis=0))) > 1e-10:
        raise ValueError("Readout must be centered under the uniform stationary law")
    coefficients = eigenvectors.T @ centered
    return (coefficients.T * average_response(eigenvalues*window)) @ coefficients / len(centered)


def graph_covariance(lap, kappa, window, rate=1.0):
    if window <= 0 or rate <= 0 or kappa < 0 or not np.all(np.isfinite([kappa, window, rate])):
        raise ValueError("Positive window/rate and nonnegative kappa required")
    lam, vec = la.eigh(lap)
    return (vec[:, 1:] * (kappa*average_response(rate*lam[1:]*window/2))) @ vec[:, 1:].T


def moment_variance(q, observable, window):
    """Independent Feynman-Kac moment ODE, not graph covariance insertion."""
    size = len(q)
    generator, zero = sp.csr_matrix(q), sp.csr_matrix((size, size))
    mark = sp.diags(observable)
    block = sp.bmat([[generator, zero, zero], [mark, generator, zero],
                     [zero, 2*mark, generator]], format="csr")
    initial = np.zeros(3*size); initial[:size] = 1.0
    result = expm_multiply(window*block, initial)
    mean = result[size:2*size].mean()/window
    return float(result[2*size:].mean()/window**2-mean**2)


def density_split(readout, centered_load, kappa):
    f = np.asarray(readout, dtype=float)
    if f.ndim == 1:
        f = f[:, None]
    b = f.T @ centered_load/(len(f)*kappa)
    residual = f-centered_load@b.T
    return b, residual


def harmonic_coefficient_bounds(directions, multipoles, kappa, gain_l1=1.0):
    """Addition-theorem bound for equal-weight empirical harmonic power.

    Directions may repeat and need not form an exact quadrature rule.
    The temporal kernel must be common, deterministic and trajectory-independent.
    """
    directions = np.asarray(directions, dtype=float)
    if directions.ndim != 2 or directions.shape[1] != 3 or not len(directions):
        raise ValueError("Nonempty array of three-dimensional directions required")
    norms = np.linalg.norm(directions, axis=1)
    if not np.all(np.isfinite(directions)) or np.any(norms <= 0) or kappa < 0 or gain_l1 < 0:
        raise ValueError("Finite nonzero directions and nonnegative variance/gain required")
    directions = directions/norms[:, None]
    theta = np.arccos(np.clip(directions[:,2], -1, 1))
    phi = np.arctan2(directions[:,1], directions[:,0])
    ports = len(directions)
    result = []
    for ell in multipoles:
        if ell < 0 or int(ell) != ell:
            raise ValueError("Nonnegative integer multipoles required")
        coefficients = np.array([(4*np.pi/ports)*sph_harm_y(ell, m, theta, phi).conj()
                                 for m in range(-ell, ell+1)])
        projected = coefficients-coefficients.mean(axis=1, keepdims=True)
        raw_norm = float(np.sum(np.abs(coefficients)**2)/(2*ell+1))
        projected_norm = float(np.sum(np.abs(projected)**2)/(2*ell+1))
        result.append({"ell": ell, "mean_coefficient_squared_norm": raw_norm,
                       "mean_projected_coefficient_squared_norm": projected_norm,
                       "addition_theorem_absolute_error": abs(raw_norm-4*np.pi/ports),
                       "constant_projection_norm_increase": max(0.0, projected_norm-raw_norm),
                       "projected_pseudo_Cl_upper_bound": kappa*gain_l1**2*projected_norm,
                       "universal_pseudo_Cl_upper_bound": kappa*gain_l1**2*4*np.pi/ports})
    return result


def experiment(name, ports, raised, edges, windows, rule=None, quadrature=False):
    states, q, lap = chain(ports, raised, edges, rule)
    z = states-raised/ports
    kappa = raised*(ports-raised)/(ports*(ports-1))
    pi = np.eye(ports)-np.ones((ports, ports))/ports
    operator_eigen, operator_vec = la.eigh(-q)
    lam, vec = la.eigh(lap)
    first_a, first_b = edges[0][:2]
    product = states[:, first_a]*states[:, first_b]
    product = product-product.mean()
    # The extra readout is a directly accessible local two-port observable.
    # It is an independent control, never a claimed geometric volume.
    f = np.column_stack([z[:, 0], z@lap[:, 0], product])
    b, residual = density_split(f, z, kappa)
    rows = []
    for window in windows:
        expected = graph_covariance(lap, kappa, window)
        observed = spectral_covariance(operator_eigen, operator_vec, z, window)
        total = spectral_covariance(operator_eigen, operator_vec, f, window)
        rest = spectral_covariance(operator_eigen, operator_vec, residual, window)
        component = b@expected@b.T
        checks = {"graph_vs_configuration_absolute_error": float(np.max(np.abs(expected-observed))),
                  "density_residual_split_absolute_error": float(np.max(np.abs(total-component-rest))),
                  "residual_min_eigenvalue": float(la.eigvalsh(rest)[0]),
                  "moment_ode_load_variance": moment_variance(q, f[:, 0], window),
                  "moment_ode_product_variance": moment_variance(q, f[:, 2], window)}
        if quadrature:
            direct, error = quad_vec(lambda t:2*(window-t)/window**2*(z.T@la.expm(q*t)@z)/len(z),
                                     0, window, epsabs=1e-11, epsrel=1e-11)
            checks["quadrature_covariance_absolute_error"] = float(np.max(np.abs(expected-direct)))
            checks["quadrature_error_estimate"] = float(error)
        mode_variance = kappa*average_response(lam[1:]*window/2)
        energy = float(np.dot(lam[1:], mode_variance)/2)
        rows.append({"window_native_time": window,
                     "mode_variances_native": mode_variance.tolist(),
                     "mode_lambda_times_variance": (lam[1:]*mode_variance).tolist(),
                     "time_average_port_rms": float(np.sqrt(np.trace(expected)/ports)),
                     "expected_load_Dirichlet_energy": energy,
                     "energy_per_positive_mode_times_two": 2*energy/(ports-1),
                     "long_window_Green_amplitude": 4*kappa/window,
                     "readback_covariance": total.tolist(),
                     "nonlinear_residual_covariance": rest.tolist(),
                     "checks": checks})
    # Equal-area quadrature is bookkeeping; it is not proof of geometric attachment.
    mass = 4*np.pi/ports
    c = vec[:, 1:].T@z[0]
    weighted_c = np.sqrt(mass)*c
    normalization = {"declared_equal_area_weight": mass,
                     "coefficient_rule": "a_j=sqrt(w)*c_j; V(a_j)=w*V(c_j)",
                     "operator_rule": "If K=L is the screen stiffness and M=w I, mu_j=lambda_j/w",
                     "weighted_Parseval_absolute_error": float(abs(np.dot(weighted_c, weighted_c)-mass*np.dot(z[0], z[0]))),
                     "weighted_energy_absolute_error": float(abs(np.dot(lam[1:]/mass, weighted_c**2)-z[0]@lap@z[0])),
                     "qualification": "These algebraic identities do not identify native K with the unit-sphere metric stiffness, or load with geometric q."}
    return {"name": name, "ports": ports, "raised": raised, "configurations": len(states),
            "seams": len(edges), "edge_rates": [list(x) for x in edges],
            "kappa": kappa, "positive_graph_eigenvalues": lam[1:].tolist(),
            "generator_source": "pinned production rule" if rule else "independent direct nearest-integer arithmetic",
            "generator_sha256": hashlib.sha256(q.astype("<f8").tobytes()).hexdigest(),
            "stationarity_absolute_error": float(np.max(np.abs(q.sum(axis=0)))),
            "detailed_balance_absolute_error": float(np.max(np.abs(q-q.T))),
            "drift_absolute_error": float(np.max(np.abs(q@z+z@lap/2))),
            "static_covariance_absolute_error": float(np.max(np.abs(z.T@z/len(z)-kappa*pi))),
            "local_readback_names": ["centered_load_at_port0", "laplacian_drive_at_port0", "centered_first_seam_occupancy_product"],
            "nonlinear_readback_density_coefficients": b[2].tolist(),
            "density_residual_orthogonality_error": float(np.max(np.abs(residual.T@z/len(z)))),
            "windows": rows, "continuum_normalization": normalization}


def build():
    spec = json.loads((HERE/"spec.json").read_text())
    helper_path = CODEX/"dynamics/run.py"
    module_spec = importlib.util.spec_from_file_location("source_native_pins", helper_path)
    helper = importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(helper)
    federation, pins = helper.pinned_source(WORKSPACE/"oph-physics-sim", spec["source_commit"])
    from oph_exact import carrier
    carrier_spec = spec["carrier"]
    source = experiment("native_twelve_port_carrier", carrier_spec["ports"], carrier_spec["raised"],
                        [(a,b,1.0) for a,b in carrier.seams()], carrier_spec["windows"],
                        carrier.integer_nearest_agreement)
    controls = [experiment(row["name"], row["ports"], row["raised"], row["edges"],
                           spec["control_windows"], quadrature=True) for row in spec["controls"]]
    harmonics = []
    for level in spec["harmonic_bound_control"]["federation_levels"]:
        centers = federation.geodesic_icosahedral_patch_arrays(level, patch_basis="cells")[0]
        directions = np.repeat(centers, 12, axis=0)
        p = len(directions)
        harmonics.append({"level": level, "carriers": len(centers), "ports": p,
                          "directions_sha256": hashlib.sha256(directions.astype("<f8").tobytes()).hexdigest(),
                          "bounds": harmonic_coefficient_bounds(directions, spec["harmonic_bound_control"]["multipoles"], p/(4*(p-1)))})
    result = {"schema": "oph.native-finite-source.receipt.v1",
              "inputs_sha256": {"spec.json": sha(HERE/"spec.json"), "run.py": sha(HERE/"run.py"),
                                "../../dynamics/run.py": sha(helper_path)},
              "production_source_commit": spec["source_commit"], "production_source_sha256": pins,
              "native_carrier": source, "controls": controls,
              "empirical_harmonic_bounds": harmonics,
              "fixed_area_refinement_bound": [
                  {"ports": p, "raised": p//2, "area_weight": 4*np.pi/p,
                   "kappa": p/(4*(p-1)),
                   "all_windows_all_rates_mode_variance_upper_bound": np.pi/(p-1)}
                  for p in spec["fixed_area_refinement_control"]["even_port_counts"]],
              "observational_data_read": False,
              "curvature_identification": "UNIDENTIFIED: no attached J_X/Jbar_X observed",
              "nonclaims": spec["nonclaims"]}
    validate(result, spec["absolute_tolerance"])
    result["numerical_identity_status"] = "PASS"
    return result


def validate(result, tolerance=1e-9):
    if result["observational_data_read"] is not False:
        raise AssertionError("Source-only input boundary violated")
    for graph in result["empirical_harmonic_bounds"]:
        for row in graph["bounds"]:
            if row["addition_theorem_absolute_error"] > tolerance or row["constant_projection_norm_increase"] > tolerance:
                raise AssertionError("Addition theorem or projection bound failed")
    for graph in [result["native_carrier"]]+result["controls"]:
        for key in ["stationarity_absolute_error", "detailed_balance_absolute_error", "drift_absolute_error",
                    "static_covariance_absolute_error", "density_residual_orthogonality_error"]:
            if not abs(graph[key]) <= tolerance:
                raise AssertionError(key)
        for row in graph["windows"]:
            checks = row["checks"]
            for key, value in checks.items():
                if "absolute_error" in key and not abs(value) <= tolerance:
                    raise AssertionError(key)
            if checks["residual_min_eigenvalue"] < -tolerance:
                raise AssertionError("Nonpositive residual covariance")
            covariance = np.asarray(row["readback_covariance"])
            for index, key in [(0, "moment_ode_load_variance"), (2, "moment_ode_product_variance")]:
                if not abs(covariance[index,index]-checks[key]) <= tolerance:
                    raise AssertionError("Independent finite-time moment mismatch")


def compare_reproduction(stored, fresh, path="receipt"):
    """Check every receipt field, allowing roundoff in numerical diagnostics."""
    if isinstance(fresh, dict):
        if set(stored) != set(fresh):
            raise AssertionError(f"Different fields: {path}")
        for key in fresh:
            compare_reproduction(stored[key], fresh[key], f"{path}.{key}")
    elif isinstance(fresh, list):
        if len(stored) != len(fresh):
            raise AssertionError(f"Different lengths: {path}")
        for index, (old, new) in enumerate(zip(stored, fresh)):
            compare_reproduction(old, new, f"{path}[{index}]")
    elif isinstance(fresh, float):
        if not np.isclose(stored, fresh, rtol=1e-9, atol=1e-11):
            raise AssertionError(f"Changed numerical result: {path}")
    elif stored != fresh:
        raise AssertionError(f"Changed result: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        stored = json.loads((HERE/"receipt.json").read_text())
        for path, value in stored["inputs_sha256"].items():
            if sha(HERE/path) != value:
                raise SystemExit(f"Input hash mismatch: {path}")
        fresh = build()
        validate(stored)
        compare_reproduction(stored, fresh)
        print("Verified input hashes, source-only covariance identities and fresh absolute spectra")
    else:
        result = build()
        (HERE/"receipt.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        print("Wrote source-only finite-window native spectrum receipt")


if __name__ == "__main__":
    main()

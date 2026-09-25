"""Native balanced-repair occupation readback; bounded exact-chain checks.

Run with BLAS/OMP threads capped to one. No Monte Carlo or CMB fitting is used.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.linalg as la
from scipy.integrate import quad
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent
CODEX = HERE.parents[1]
WORKSPACE = HERE.parents[3]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def laplacian(n, edges):
    result = np.zeros((n, n), dtype=int)
    for i, j in edges:
        result[i, i] += 1
        result[j, j] += 1
        result[i, j] -= 1
        result[j, i] -= 1
    return result


def configuration_generator(n, raised, edges):
    """Exact generator: each seam attempts at rate1 and swaps with chance1/2."""
    states = [tuple(int(i in s) for i in range(n))
              for s in itertools.combinations(range(n), raised)]
    where = {state: i for i, state in enumerate(states)}
    q = [[F(0) for _ in states] for _ in states]
    for i, state in enumerate(states):
        for a, b in edges:
            changed = list(state)
            changed[a], changed[b] = changed[b], changed[a]
            j = where[tuple(changed)]
            q[i][j] += F(1, 2)
            q[i][i] -= F(1, 2)
    return states, q


def exact_checks(states, q, edges):
    n, size = len(states[0]), len(states)
    raised = sum(states[0])
    lap = laplacian(n, edges)
    kappa = F(raised * (n-raised), n*(n-1))
    assert all(sum(row) == 0 for row in q)
    assert all(q[i][j] == q[j][i] for i in range(size) for j in range(size))
    assert all(q[i][j] >= 0 for i in range(size) for j in range(size) if i != j)
    for i, state in enumerate(states):
        for a in range(n):
            drift = sum(q[i][j] * states[j][a] for j in range(size))
            assert drift == -sum(F(int(lap[a, b]), 2) * state[b] for b in range(n))
    for a in range(n):
        for b in range(n):
            covariance = sum((F(x[a])-F(raised,n)) * (F(x[b])-F(raised,n)) for x in states) / size
            assert covariance == kappa * (int(a == b)-F(1,n))
    return kappa


def occupation_variance(rate, window, variance=1.0):
    """Variance of T^-1/2 integral f, for C_f(t)=variance*exp(-rate*t)."""
    if rate < 0 or window <= 0 or variance < 0:
        raise ValueError("Nonnegative rate/variance and positive window required")
    u = rate * window
    if u < 1e-3:
        factor = 1-u/3+u*u/12-u**3/60+u**4/360-u**5/2520
        return variance * window * factor
    return 2 * variance / rate * (1 + np.expm1(-u) / u)


def leaky_variance(rate, leak, variance=1.0):
    """Variance of sqrt(2*leak) integral exp(-leak*s) f(t-s) ds."""
    if rate < 0 or leak <= 0 or variance < 0:
        raise ValueError("Nonnegative rate/variance and positive leak required")
    return 2 * variance / (rate + leak)


def occupation_moments(q, observable, window, order=4):
    """Independent Feynman-Kac moment ODE, solved by one block exponential."""
    size = len(q)
    operator = np.zeros(((order+1)*size, (order+1)*size))
    d = np.diag(observable)
    for degree in range(order+1):
        rows = slice(degree*size, (degree+1)*size)
        operator[rows, rows] = q
        if degree:
            operator[rows, slice((degree-1)*size, degree*size)] = degree*d
    initial = np.zeros((order+1)*size)
    initial[:size] = 1
    values = la.expm(operator*window) @ initial
    return np.array([values[j*size:(j+1)*size].mean() / window**(j/2)
                     for j in range(order+1)])


def discrete_occupation_variance(rho, attempts, seams, variance=1.0):
    """Exact rectangle-sum readout at equally spaced attempted moves."""
    lags = np.arange(1, attempts)
    return variance / seams * (1 + 2*np.sum((1-lags/attempts)*rho**lags))


def chain_experiment(description, spec):
    n, raised, edges = description["ports"], description["raised"], description["edges"]
    states, rational_q = configuration_generator(n, raised, edges)
    kappa = exact_checks(states, rational_q, edges)
    q = np.array(rational_q, dtype=float)
    centered = np.array(states, dtype=float) - raised/n
    lap = laplacian(n, edges).astype(float)
    eigen, vectors = la.eigh(lap)
    modes = centered @ vectors[:, 1:]
    covariance_errors = []
    for mode, eigenvalue in enumerate(eigen[1:]):
        rate = eigenvalue/2
        f = modes[:, mode]
        for u in spec["dimensionless_lags"]:
            observed = f @ (la.expm(q*(u/rate)) @ f) / len(states)
            covariance_errors.append(abs(observed-float(kappa)*np.exp(-u)))
    clt_covariance = 2 * centered.T @ la.pinvh(-q) @ centered / len(states)
    graph_covariance = 4 * float(kappa) * la.pinvh(lap)
    first, rate = modes[:, 0], eigen[1]/2
    rows = []
    for u in spec["dimensionless_windows"]:
        window = u/rate
        analytic = occupation_variance(rate, window, float(kappa))
        integrand = lambda lag: 2*(window-lag)*(first @ la.expm(q*lag) @ first)/len(states)/window
        integrated, quad_error = quad(integrand, 0, window, epsabs=1e-10, epsrel=1e-10, limit=200)
        moments = occupation_moments(q, first, window)
        rows.append({"aT": u, "window_sweeps": float(window),
                     "variance_analytic": float(analytic), "variance_semigroup_quadrature": float(integrated),
                     "quadrature_error_estimate": float(quad_error), "variance_moment_ode": float(moments[2]),
                     "fourth_moment_normalized_readout": float(moments[4]),
                     "kurtosis": float(moments[4]/moments[2]**2),
                     "variance_over_gff_limit": float(analytic / (2*float(kappa)/rate))})
    leaks = []
    for ratio in spec["leak_over_mode_rate"]:
        leak = ratio*rate
        observed = 2 * first @ la.solve(leak*np.eye(len(q))-q, first, assume_a="pos") / len(states)
        analytic = leaky_variance(rate, leak, float(kappa))
        leaks.append({"gamma_over_a": ratio, "variance_analytic": float(analytic),
                      "variance_configuration_resolvent": float(observed)})
    discrete = []
    transition = np.eye(len(q)) + q/len(edges)
    rho = 1-rate/len(edges)
    for attempts in (4,16,64):
        current = first.copy()
        variance = float(kappa)
        for lag in range(1, attempts):
            current = transition @ current
            variance += 2*(1-lag/attempts)*(first@current)/len(states)
        observed = variance/len(edges)
        analytic = discrete_occupation_variance(rho, attempts, len(edges), float(kappa))
        discrete.append({"attempts": attempts, "variance_chain": float(observed),
                         "variance_formula": float(analytic),
                         "independent_reshuffle_white_variance": float(kappa)/len(edges)})
    return {"name": description["name"], "ports": n, "raised": raised,
            "states": states, "generator_rationals": [[str(x) for x in row] for row in rational_q],
            "exact_stationarity_drift_and_covariance_checks": True, "kappa_exact": str(kappa),
            "graph_eigenvalues": eigen.tolist(), "first_mode_vector": vectors[:,1].tolist(),
            "max_semigroup_covariance_absolute_error": float(max(covariance_errors)),
            "clt_covariance_max_absolute_error": float(np.max(np.abs(clt_covariance-graph_covariance))),
            "windows": rows, "leak_controls": leaks, "discrete_controls": discrete}


def graph_experiment(federation, level, spec):
    fed = federation.build_federation(level, "port_pair")
    lap = fed.laplacian.astype(float)
    count = spec["actual_graph_modes"]
    if level == 0:
        eigen, vectors = la.eigh(lap.toarray(), subset_by_index=(0,count-1))
    else:
        eigen, vectors = sla.eigsh(lap, k=count, sigma=-1e-6, which="LM", tol=1e-11,
                                   ncv=96, v0=np.random.default_rng(20260925+level).normal(size=fed.ports))
        order = np.argsort(eigen)
        eigen, vectors = eigen[order], vectors[:,order]
    residual = np.linalg.norm(lap @ vectors - vectors*eigen[None,:], axis=0)
    raised = fed.ports//2
    kappa = raised*(fed.ports-raised)/(fed.ports*(fed.ports-1))
    rates = eigen[1:]/2
    windows = []
    for factor in spec["actual_graph_window_over_slowest_relaxation"]:
        window = factor/rates[0]
        values = [occupation_variance(a, window, kappa) for a in rates]
        windows.append({"window_over_slowest_relaxation": factor, "window_sweeps": float(window),
                        "mode_variances": list(map(float,values)),
                        "mode_variance_over_gff_limit": (np.array(values)*rates/(2*kappa)).tolist()})
    leak_rows = []
    for ratio in spec["actual_graph_leak_over_slowest_rate"]:
        leak = ratio*rates[0]
        leak_rows.append({"leak_over_slowest_rate": ratio, "gamma_per_sweep": float(leak),
                          "mode_variances": [float(leaky_variance(a,leak,kappa)) for a in rates],
                          "readout_absolute_support_bound": float(np.sqrt(2/leak)/2)})
    degrees, counts = np.unique(lap.diagonal(), return_counts=True)
    return {"level": level, "carriers": fed.carriers, "ports": fed.ports, "seams": fed.seams,
            "seam_endpoints_sha256": hashlib.sha256(np.stack([fed.seam_a,fed.seam_b],axis=1).astype("<i8").tobytes()).hexdigest(),
            "raised": raised, "stationary_mode_variance": kappa,
            "positive_eigenvalues": eigen[1:].tolist(),
            "max_eigenpair_relative_residual": float(np.max(residual[1:]/eigen[1:])),
            "slowest_mode_relaxation_sweeps": float(1/rates[0]),
            "mode_gff_limit_variances": (2*kappa/rates).tolist(),
            "windows": windows, "leaks": leak_rows,
            "port_attempt_rates_per_sweep": {str(int(d)):int(c) for d,c in zip(degrees,counts)},
            "slowest_mode_discrete_vs_poisson_relative_limit_correction": float(eigen[1]/(4*fed.seams)),
            "scope": "Analytical covariance prediction evaluated on actual graph eigenmodes; no large configuration chain or trajectory ensemble simulated."}


def build():
    spec = json.loads((HERE/"spec.json").read_text())
    helper_path = CODEX/"dynamics/run.py"
    module_spec = importlib.util.spec_from_file_location("codex_native_history_geometry",helper_path)
    helper = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(helper)
    federation, pins = helper.pinned_source(WORKSPACE/"oph-physics-sim",spec["source_commit"])
    return {"schema":"oph.native-occupation-readback.receipt.v1",
            "inputs_sha256":{"spec.json":sha(HERE/"spec.json"),"run.py":sha(Path(__file__)),
                               "../../dynamics/run.py":sha(helper_path)},
            "production_source_commit":spec["source_commit"], "production_source_sha256":pins,
            "chains":[chain_experiment(row,spec) for row in spec["exact_chains"]],
            "actual_graphs":[graph_experiment(federation,level,spec) for level in spec["actual_graph_levels"]],
            "formulae":{"occupation_mode":"2*kappa/a * [1-(1-exp(-a*T))/(a*T)], a=lambda/2",
                        "fixed_graph_clt_covariance":"4*kappa*L^+",
                        "finite_leak_mode":"4*kappa/(lambda+2*gamma)",
                        "discrete_attempt_limit":"4*kappa/lambda-kappa/m"},
            "nonclaims":spec["nonclaims"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    result=build()
    (HERE/"receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("Wrote native occupation-readback receipt")


if __name__=="__main__":
    main()

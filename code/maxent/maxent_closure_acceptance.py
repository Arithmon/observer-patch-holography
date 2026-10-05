#!/usr/bin/env python3
"""Acceptance test for GitHub issue #539: MaxEnt multiplier counting and closure defect.

Issue #539 requires, on two successive regulator lattices:

1. The number of independent constraints and of Lagrange multipliers must match the
   parameter-space dimension displayed by the papers (N_con homogeneous global sums plus
   N_glob optional global charges), independent of the number of regulator cells. The
   rejected per-cell reading is counted alongside to show its dimension grows with cells.
2. The coarse-grained MaxEnt state must lie in the claimed family with a proved residual
   bound. This is the I-projection residual bound (spacetime and Einstein paper, Lemma
   ``lem:closure-residual``; synthesis paper, Lemma 2.6b): the moment-matching
   I-projection R(lambda) onto the coarse homogeneous family exists and is unique, the
   closure defect is eps = D(sigma || omega_L(R(lambda))), and
   || sigma - omega_L(R(lambda)) ||_1 <= sqrt(2 * eps).

The model is the transverse-field Ising family on a spin-1/2 ring: density labels
O_1(x) = Z_x Z_{x+1} and O_2(x) = X_x, constrained through their homogeneous global sums,
so N_con = 2 and N_glob = 0. The refinement channel is decimation (partial trace over odd
sites), a completely positive trace-preserving coarse-graining.

Both regimes relevant to the separate optimizer-pushforward interface are exhibited:

- generic multipliers: the closure defect is strictly positive, so closure under one fixed
  finite exponential family is a substantive renormalization condition, and the residual
  bound quantifies exactly how far the coarse-grained state sits from the family;
- the transverse-field product subfamily (lambda_ZZ = 0): decimation lands exactly in the
  coarse family, the defect vanishes, and R acts as the identity on the multipliers.

Run directly for a full receipt (written to runs/maxent_closure_acceptance_receipt.json):

    python3 maxent_closure_acceptance.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from quantum_information import (
    faithful_density_matrix, finite_real_scalar, partial_trace,
    relative_entropy as _relative_entropy,
)

HERE = Path(__file__).resolve().parent

PAULI_I = np.eye(2, dtype=complex)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def site_operator(op: np.ndarray, site: int, n_sites: int) -> np.ndarray:
    """op acting on `site` of an n_sites spin-1/2 chain, identity elsewhere."""
    out = np.array([[1.0 + 0.0j]])
    for k in range(n_sites):
        out = np.kron(out, op if k == site else PAULI_I)
    return out


def cell_densities(n_sites: int) -> list[list[np.ndarray]]:
    """Per-cell density operators [O_1(x)]_x, [O_2(x)]_x on a ring of n_sites cells."""
    zz = [
        site_operator(PAULI_Z, x, n_sites) @ site_operator(PAULI_Z, (x + 1) % n_sites, n_sites)
        for x in range(n_sites)
    ]
    xs = [site_operator(PAULI_X, x, n_sites) for x in range(n_sites)]
    return [zz, xs]


def global_sum_constraints(n_sites: int) -> list[np.ndarray]:
    """The homogeneous global-sum constrained operators S_a = sum_x O_a(x)."""
    return [sum(ops) for ops in cell_densities(n_sites)]


def _validated_constraints(constraints):
    if not constraints:
        raise ValueError("nonempty constraint family required")
    operators = [np.asarray(s, dtype=complex) for s in constraints]
    shape = operators[0].shape
    if len(shape) != 2 or shape[0] == 0 or shape[0] != shape[1]:
        raise ValueError("square constraint operators required")
    for s in operators:
        if (s.shape != shape or not np.all(np.isfinite(s))
                or np.linalg.norm(s-s.conj().T, ord="fro") > 1e-12):
            raise ValueError("constraints must be finite Hermitian operators on one algebra")
    return [s/2+s.conj().T/2 for s in operators]


def _validated_multipliers(lam, count):
    multipliers = np.asarray(lam)
    if (multipliers.shape != (count,) or multipliers.dtype.kind not in "iuf"
            or not np.all(np.isfinite(multipliers))):
        raise ValueError("one finite real multiplier is required per constraint")
    return multipliers


def _centered_constraints(operators):
    """Remove each identity term before sums or spectral calculations lose it.

    Subtract a diagonal anchor before taking the mean. This preserves a
    representable small diagonal difference beside a large identity shift.
    """
    centered, offsets = [], []
    for operator in operators:
        a = operator.copy()
        diagonal = np.diag_indices_from(a)
        anchor = float(a[0, 0].real)
        with np.errstate(over="ignore", invalid="ignore"):
            a[diagonal] -= anchor
            mean = float(np.sum(a.diagonal().real/len(a)))
            a[diagonal] -= mean
        if not np.all(np.isfinite(a)) or not np.isfinite(mean):
            raise ValueError("centered constraints exceed finite numerical range")
        centered.append(a)
        # Keep the small mean separate: other identity anchors may cancel.
        offsets.append((anchor, mean))
    return centered, np.asarray(offsets)


def _constraint_coordinates(operators):
    """Real HS singular vectors of unit-scaled, traceless observables.

    Unlike a Gram matrix, this does not square the condition number. Exact
    dependence and numerically unresolved independence share the rank guard.
    """
    centered, _ = _centered_constraints(operators)
    scales = np.array([np.max(np.abs(s)) for s in centered])
    scaled = np.array([s.real/scale + 1j*(s.imag/scale) if scale else s
                       for s, scale in zip(centered, scales)])
    count, dim, _ = scaled.shape
    rows = np.concatenate((scaled.real.reshape(count, -1),
                           scaled.imag.reshape(count, -1)), axis=1)
    u, singular, vh = np.linalg.svd(rows, full_matrices=False)
    threshold = np.finfo(float).eps * max(rows.shape) * singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    basis = (vh[:rank, :dim*dim] + 1j*vh[:rank, dim*dim:]).reshape(rank, dim, dim)
    return list(basis), scales, u, singular


def independent_operator_count(operators: list[np.ndarray]) -> int:
    """Numerically resolved real dimension modulo identity, independent of units."""
    basis, _, _, _ = _constraint_coordinates(_validated_constraints(operators))
    return len(basis)


def constrained_hamiltonian(constraints, lam):
    """Reject missing multipliers and invalid observables before summing."""
    operators = _validated_constraints(constraints)
    multipliers = _validated_multipliers(lam, len(operators))
    with np.errstate(over="ignore", invalid="ignore"):
        ham = sum(l*s for l, s in zip(multipliers, operators))
    if not np.all(np.isfinite(ham)):
        raise ValueError("constraint combination exceeds finite numerical range")
    return ham


def _gibbs_data(constraints, lam):
    operators = _validated_constraints(constraints)
    multipliers = _validated_multipliers(lam, len(operators))
    centered, offsets = _centered_constraints(operators)
    with np.errstate(over="ignore", invalid="ignore"):
        ham = sum(l*s for l, s in zip(multipliers, centered))
    try:
        offset = math.fsum(float(l)*float(part) for l, pair in zip(multipliers, offsets)
                           for part in pair)
    except (OverflowError, ValueError) as error:
        raise ValueError("constraint combination exceeds finite numerical range") from error
    if not np.all(np.isfinite(ham)) or not np.isfinite(offset):
        raise ValueError("constraint combination exceeds finite numerical range")
    energies, vectors = np.linalg.eigh(ham)
    with np.errstate(over="ignore"):
        weights = np.exp(-(energies - energies.min()))
    if np.any(weights == 0):
        raise ValueError("Gibbs spectrum underflow; faithful-state precision is insufficient")
    log_z = math.log(weights.sum()) - energies.min() - offset
    if not np.isfinite(log_z):
        raise ValueError("Gibbs log partition exceeds finite numerical range")
    probs = weights/weights.sum()
    if np.any(probs == 0):
        raise ValueError("normalized Gibbs spectrum underflow; precision is insufficient")
    rho = faithful_density_matrix((vectors * probs) @ vectors.conj().T)
    return rho, log_z, probs, vectors, centered


def gibbs_state(constraints: list[np.ndarray], lam: np.ndarray) -> tuple[np.ndarray, float]:
    """omega(lambda) and log Z, with identity energies removed before diagonalizing."""
    rho, log_z, _, _, _ = _gibbs_data(constraints, lam)
    return rho, log_z


def decimate(rho: np.ndarray, n_sites: int) -> np.ndarray:
    """Partial trace over the odd sites of a ring, keeping sites 0, 2, 4, ..."""
    if type(n_sites) is not int or n_sites <= 0 or n_sites % 2:
        raise ValueError("decimation channel needs an even number of fine sites")
    return partial_trace(rho, [2]*n_sites, list(range(0,n_sites,2)))


def relative_entropy(sigma: np.ndarray, rho: np.ndarray) -> float:
    """D(sigma || rho); preserve the historical keyword-argument convention."""
    return _relative_entropy(sigma,rho)


def trace_norm(delta: np.ndarray) -> float:
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta + delta.conj().T) / 2))))


def duhamel_covariance(constraints: list[np.ndarray], lam: np.ndarray) -> np.ndarray:
    """Kubo-Mori covariance matrix of the constraints: the Hessian of log Z(lambda).

    K_ab = d^2 log Z / (d lambda_a d lambda_b); positive definite iff the constrained
    operators together with the identity are linearly independent, which is the strict
    convexity input of the I-projection lemma.
    """
    rho, _, probs, vectors, operators = _gibbs_data(constraints, lam)
    centered = [
        vectors.conj().T @ (s - np.real(np.trace(rho @ s)) * np.eye(s.shape[0])) @ vectors
        for s in operators
    ]
    logp = np.log(probs)
    # Stable logarithmic mean, including nearly equal probabilities.
    delta = np.abs(logp[:, None] - logp[None, :])
    ratio = np.divide(-np.expm1(-delta), delta, out=np.ones_like(delta), where=delta != 0)
    kernel = np.maximum(probs[:, None], probs[None, :]) * ratio
    n_con = len(constraints)
    cov = np.empty((n_con, n_con))
    with np.errstate(over="ignore", invalid="ignore"):
        for a in range(n_con):
            for b in range(n_con):
                cov[a, b] = float(np.real(np.sum(kernel * centered[a] * centered[b].T)))
    if not np.all(np.isfinite(cov)):
        raise ValueError("Duhamel covariance exceeds finite numerical range")
    return cov/2 + cov.T/2


def i_projection(
    sigma: np.ndarray, constraints: list[np.ndarray], tol: float = 1e-11, max_iter: int = 200
) -> tuple[np.ndarray, float]:
    """Unique Gibbs I-projection, solved in orthonormal traceless HS coordinates.

    Requires a faithful target and numerically independent constraints modulo
    identity. Centering, scaling and an SVD change only multiplier coordinates,
    not the exponential family. tol and the returned moment residual use this
    dimensionless basis, so observable units cannot manufacture convergence.
    Multipliers are returned in the caller's original coordinates; their Gibbs
    state must pass the same residual check after conversion. Unresolved rank,
    conversion precision, or failed convergence raises instead of certifying.
    """
    sigma = faithful_density_matrix(sigma)
    constraints = _validated_constraints(constraints)
    if constraints[0].shape != sigma.shape:
        raise ValueError("target and constraints must use the same algebra")
    basis, scales, u, singular = _constraint_coordinates(constraints)
    if len(basis) != len(constraints):
        raise ValueError("constraints must be independent modulo identity at numerical precision")
    tol = finite_real_scalar(tol, "convergence tolerance")
    if tol <= 0:
        raise ValueError("finite positive convergence tolerance required")
    if type(max_iter) is not int or max_iter <= 0:
        raise ValueError("positive integer iteration budget required")

    def moments(state):
        return np.array([float(np.trace(state @ s).real) for s in basis])

    targets = moments(sigma)
    theta = np.zeros(len(basis))
    rho, log_z = gibbs_state(basis, theta)
    grad = moments(sigma-rho)
    for _ in range(max_iter):
        residual = float(np.linalg.norm(grad))
        if residual < tol:
            break
        hess = duhamel_covariance(basis, theta)
        step = np.linalg.solve(hess, -grad)
        slope = float(grad @ step)
        if not np.all(np.isfinite(step)) or slope >= 0:
            raise RuntimeError("information projection has no resolved descent direction")
        base = log_z + float(theta @ targets)
        scale = 1.0
        for _ in range(60):
            candidate = theta + scale*step
            try:
                trial, trial_log_z = gibbs_state(basis, candidate)
            except ValueError:
                # A full Newton step may leave representable faithful states.
                # Backtrack rather than accepting it or aborting a valid solve.
                scale /= 2
                continue
            trial_grad = moments(sigma-trial)
            value = trial_log_z + float(candidate @ targets)
            roundoff = 16*np.finfo(float).eps*max(1., abs(base))
            if (value <= base + 1e-4*scale*slope
                    or (abs(value-base) <= roundoff
                        and np.linalg.norm(trial_grad) < residual)):
                theta, rho, log_z, grad = candidate, trial, trial_log_z, trial_grad
                break
            scale /= 2
        else:
            raise RuntimeError("information projection line search did not converge")
    if np.linalg.norm(grad) >= tol:
        raise RuntimeError("information projection did not converge within its budget")

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        lam = (u @ (theta/singular))/scales
    if not np.all(np.isfinite(lam)):
        raise ValueError("projection multiplier coordinates are numerically unresolved")
    recovered, _ = gibbs_state(constraints, lam)
    residual = float(np.linalg.norm(moments(sigma-recovered)))
    if not np.isfinite(residual) or residual >= tol:
        raise RuntimeError("projection multiplier conversion did not preserve convergence")
    return lam, residual


def run_lattice_pair(n_fine: int, lam_fine: np.ndarray) -> dict:
    """Run the full issue-#539 acceptance test on the lattice pair (n_fine, n_fine/2)."""
    n_coarse = n_fine // 2
    fine = global_sum_constraints(n_fine)
    coarse = global_sum_constraints(n_coarse)
    n_con = len(fine)

    fine_count = independent_operator_count(fine)
    coarse_count = independent_operator_count(coarse)
    per_cell_fine = sum(len(ops) for ops in cell_densities(n_fine))
    per_cell_coarse = sum(len(ops) for ops in cell_densities(n_coarse))

    omega_fine, _ = gibbs_state(fine, lam_fine)
    sigma = decimate(omega_fine, n_fine)
    lam_star, _ = i_projection(sigma, coarse)
    omega_coarse, _ = gibbs_state(coarse, lam_star)
    # Keep the historical receipt in original observable units; the solver
    # separately certifies its dimensionless convergence residual.
    moment_residual = float(np.linalg.norm([
        np.trace((sigma-omega_coarse) @ s).real for s in coarse
    ]))

    defect = max(relative_entropy(sigma, omega_coarse), 0.0)
    residual = trace_norm(sigma - omega_coarse)
    pinsker_bound = math.sqrt(2 * defect)
    hess_floor = float(np.linalg.eigvalsh(duhamel_covariance(coarse, lam_star)).min())

    return {
        "lattice_pair": [n_fine, n_coarse],
        "fine_multipliers": list(map(float, lam_fine)),
        "displayed_dimension_N_con_plus_N_glob": n_con,
        "independent_global_sum_constraints": {"fine": fine_count, "coarse": coarse_count},
        "rejected_per_cell_constraint_count": {"fine": per_cell_fine, "coarse": per_cell_coarse},
        "induced_map_R_multipliers": list(map(float, lam_star)),
        "moment_matching_residual": moment_residual,
        "duhamel_hessian_min_eigenvalue": hess_floor,
        "closure_defect_nats": defect,
        "trace_norm_residual": residual,
        "pinsker_residual_bound": pinsker_bound,
        "counts_match_displayed_dimension": fine_count == coarse_count == n_con,
        "projection_unique": hess_floor > 1e-9,
        "residual_bound_holds": residual <= pinsker_bound + 1e-9,
    }


def run_acceptance() -> dict:
    generic = np.array([0.7, 0.4])
    closed = np.array([0.0, 0.4])
    results = {
        "issue": 539,
        "model": "transverse-field Ising densities O_1(x)=Z_x Z_{x+1}, O_2(x)=X_x on a ring;"
        " homogeneous global-sum constraints; decimation refinement channel",
        "generic_branch": [run_lattice_pair(6, generic), run_lattice_pair(8, generic)],
        "closed_subfamily": run_lattice_pair(6, closed),
    }
    closed_run = results["closed_subfamily"]
    checks = {
        "multiplier_and_constraint_counts_cutoff_independent": all(
            r["counts_match_displayed_dimension"]
            for r in results["generic_branch"] + [closed_run]
        ),
        "i_projection_unique_everywhere": all(
            r["projection_unique"] for r in results["generic_branch"] + [closed_run]
        ),
        "pinsker_residual_bound_holds_everywhere": all(
            r["residual_bound_holds"] for r in results["generic_branch"] + [closed_run]
        ),
        "generic_closure_defect_strictly_positive": all(
            r["closure_defect_nats"] > 1e-6 for r in results["generic_branch"]
        ),
        "transverse_field_subfamily_exactly_closed": closed_run["closure_defect_nats"] < 1e-10
        and abs(closed_run["induced_map_R_multipliers"][0]) < 1e-8
        and abs(closed_run["induced_map_R_multipliers"][1] - 0.4) < 1e-8,
    }
    results["checks"] = checks
    results["all_checks_pass"] = all(checks.values())
    return results


def main() -> int:
    results = run_acceptance()
    receipt = HERE / "runs" / "maxent_closure_acceptance_receipt.json"
    receipt.parent.mkdir(exist_ok=True)
    receipt.write_text(json.dumps(results, indent=2) + "\n")
    for name, passed in results["checks"].items():
        print(f"{'PASS' if passed else 'FAIL'}  {name}")
    for run in results["generic_branch"]:
        print(
            f"lattices {run['lattice_pair']}: constraints/multipliers = "
            f"{run['independent_global_sum_constraints']} (displayed "
            f"{run['displayed_dimension_N_con_plus_N_glob']}; per-cell reading would be "
            f"{run['rejected_per_cell_constraint_count']}), closure defect = "
            f"{run['closure_defect_nats']:.6f} nats, trace-norm residual "
            f"{run['trace_norm_residual']:.6f} <= bound {run['pinsker_residual_bound']:.6f}"
        )
    print(f"receipt: {receipt}")
    return 0 if results["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

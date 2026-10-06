"""Finite Gibbs inference in the real observable space modulo the identity.

Unit changes, invertible real changes of basis and scalar energy offsets do
not change that space. Numerical rank and support must still be resolved.
See PROJECTION_COORDINATES.md for the proofs and precision boundary.
"""

from dataclasses import dataclass
import math

import numpy as np

from quantum_information import faithful_density_matrix, finite_real_scalar
from quantum_information.gibbs import _finite, _numeric, _observables, _thermal
from .hamiltonian_assembly import _assemble_hamiltonian


def _multipliers(lam, count):
    result = _numeric(lam, "multipliers", real=True)
    if result.shape != (count,):
        raise ValueError("one finite real multiplier is required per constraint")
    return result


def _combination(operators, coefficients):
    return _assemble_hamiltonian(operators, coefficients, centered=False).matrix


def _hamiltonian_parts(constraints, lam):
    centered, scales, _ = _observables(constraints)
    lam = _multipliers(lam, len(centered))
    # Normalized observables supply derivative coordinates, not the original
    # Hamiltonian: normalization and product rounding can lose cancellation.
    supplied = [_numeric(a, "constraints") for a in constraints]
    assembly = _assemble_hamiltonian(supplied, lam, centered=True)
    return assembly.matrix, assembly.scalar, centered, scales


def hamiltonian_assembly_diagnostics(constraints, lam):
    """Exact-input assembly error; excludes diagonalization and state rounding.

    Returns a centered matrix and scalar, with rigorous rational-to-binary64
    error bounds. The ideal Gibbs bounds concern these Hamiltonians, not a
    certified floating-point Gibbs state. See HAMILTONIAN_ASSEMBLY.md.
    """
    _observables(constraints)
    lam = _multipliers(lam, len(constraints))
    supplied = [_numeric(a, "constraints") for a in constraints]
    if any(not np.array_equal(a, a.conj().T) for a in supplied):
        raise ValueError("exact Hermitian inputs required for assembly certificate")
    return _assemble_hamiltonian(supplied, lam, centered=True)


def constrained_hamiltonian(constraints, lam):
    """Full Hamiltonian; Gibbs evaluation keeps its scalar part separate."""
    _observables(constraints)
    lam = _multipliers(lam, len(constraints))
    return _combination([_numeric(a, "constraints") for a in constraints], lam).copy()


def gibbs_state(constraints, lam):
    """Return exp(-H)/Tr exp(-H) and log Z, without scalar-offset cancellation."""
    ham, shift, _, _ = _hamiltonian_parts(constraints, lam)
    rho, log_z, _, _ = _thermal(ham)
    return rho, float(_finite(log_z-shift, "log partition function"))


def _covariance(operators, probs, vectors):
    rotated = [vectors.conj().T @ a @ vectors for a in operators]
    for a in rotated:
        a[np.diag_indices(len(a))] -= np.dot(probs, np.diag(a)).real
    # Logarithmic mean: expm1 avoids cancellation; use max(p_i,p_j) to
    # avoid overflow even when the smaller positive eigenvalue is subnormal.
    logs = np.log(probs)
    distance = np.abs(logs[:, None]-logs[None, :])
    ratio = np.ones_like(distance)
    np.divide(-np.expm1(-distance), distance, out=ratio, where=distance != 0)
    kernel = np.maximum(probs[:, None], probs[None, :])*ratio
    result = np.array([[np.sum(kernel*a*b.T).real for b in rotated] for a in rotated])
    return _finite((result+result.T)/2, "Duhamel covariance")


def duhamel_covariance(constraints, lam):
    """Kubo-Mori covariance in the caller's units, with scalar parts removed."""
    ham, _, centered, scales = _hamiltonian_parts(constraints, lam)
    _, _, probs, vectors = _thermal(ham)
    covariance = _covariance(centered, probs, vectors)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = (covariance*scales[:, None])*scales[None, :]
    _finite(result, "Duhamel covariance")
    if np.any((covariance != 0) & (result == 0)):
        raise ValueError("Duhamel covariance underflow in the requested units")
    return result


def _coordinates(constraints):
    centered, scales, _ = _observables(constraints)
    d = len(centered[0])
    # Real coordinates in a traceless Hermitian orthonormal basis. Working
    # inside this space prevents SVD roundoff in a nearly dependent family
    # from creating a spurious trace or anti-Hermitian direction.
    real = np.column_stack([_hermitian_coordinates(a) for a in centered])
    norms = np.linalg.norm(real, axis=0)
    nonzero = norms != 0
    normalized = real[:, nonzero]/norms[nonzero]
    if not normalized.size:
        return [], None, scales, 0
    u, singular, vh = np.linalg.svd(normalized, full_matrices=False)
    floor = 64*np.finfo(float).eps*max(normalized.shape)*singular[0]
    rank = int(np.count_nonzero(singular > floor))
    basis = [_hermitian_matrix(u[:, i], d) for i in range(rank)]
    transform = (singular[:, None]*vh)*norms[nonzero]
    return basis, transform, scales, rank


def _hermitian_coordinates(a):
    d = len(a)
    diagonal = a.diagonal().real
    k = np.arange(1, d)
    helmert = (np.cumsum(diagonal)[:-1]-k*diagonal[1:])/np.sqrt(k*(k+1))
    upper = a[np.triu_indices(d, 1)]
    return np.concatenate((helmert, np.sqrt(2)*upper.real, np.sqrt(2)*upper.imag))


def _hermitian_matrix(coordinates, d):
    result = np.zeros((d, d), complex)
    diagonal = np.zeros(d)
    for k in range(1, d):
        value = coordinates[k-1]/math.sqrt(k*(k+1))
        diagonal[:k] += value
        diagonal[k] -= k*value
    result[np.diag_indices(d)] = diagonal
    count = d*(d-1)//2
    upper = (coordinates[d-1:d-1+count]+1j*coordinates[d-1+count:])/math.sqrt(2)
    indices = np.triu_indices(d, 1)
    result[indices] = upper
    result[(indices[1], indices[0])] = upper.conj()
    return result


def independent_operator_count(operators):
    """Scale-normalized numerical rank modulo identity (not an exact-rank proof)."""
    return _coordinates(operators)[3]


def _moments(state, operators):
    return np.array([np.trace(state @ a).real for a in operators])


@dataclass(frozen=True)
class ProjectionResult:
    multipliers: np.ndarray
    state: np.ndarray
    normalized_residual: float
    raw_moment_residual: float
    iterations: int
    optimality_gap_bound: float
    trace_distance_bound: float


def projection_diagnostics(sigma, constraints, lam):
    """Recompute state, invariant residual and global error bound for any candidate.

    This checks the candidate independently of the optimizer's status. The
    bound is a floating-point evaluation of the proved inequality, not an
    outward-rounded certificate; see the note for numerical limitations.
    """
    sigma = faithful_density_matrix(_numeric(sigma, "target"))
    basis, transform, scales, rank = _coordinates(constraints)
    if rank != len(constraints):
        raise ValueError("constraints must be independent modulo identity; rank is unresolved")
    if basis[0].shape != sigma.shape:
        raise ValueError("target and constraints must use the same algebra")
    lam = _multipliers(lam, len(constraints))
    rho, _ = gibbs_state(constraints, lam)
    gradient = _moments(sigma-rho, basis)
    residual = math.hypot(*gradient)
    theta = transform @ (lam*scales)
    d = len(sigma)
    mu = float(np.linalg.eigvalsh(sigma)[0])
    # Coercivity gives ||theta_*|| <= R. Convexity then bounds f-f_*.
    radius = math.log(d)*math.sqrt((d-1)/d)/mu
    bound = float(_finite(residual*(math.hypot(*theta)+radius), "optimality bound"))
    centered, scales, _ = _observables(constraints)
    with np.errstate(over="ignore", invalid="ignore"):
        raw = _finite(_moments(sigma-rho, centered)*scales, "raw moment residual")
    # math.hypot, unlike sqrt(dot(raw,raw)), does not square the units first.
    raw_residual = float(_finite(math.hypot(*raw), "raw moment residual"))
    return ProjectionResult(lam.copy(), rho, residual, raw_residual, 0,
                            bound, min(2., math.sqrt(2*bound)))


def project_information(sigma, constraints, tol=1e-11, max_iter=200):
    """Unique faithful-target Gibbs projection, solved in orthonormal coordinates.

    tol bounds the HS norm of the moment error projected onto the traceless
    constraint span. It is invariant under equivalent observable coordinates.
    No Hessian ridge, artificial support floor or discarded constraint is used.
    The returned original multipliers are replayed before acceptance.
    """
    sigma = faithful_density_matrix(_numeric(sigma, "target"))
    tol = finite_real_scalar(tol, "convergence tolerance")
    if not 0 < tol <= 1e-6:
        raise ValueError("convergence tolerance must lie in (0, 1e-6]")
    if type(max_iter) is not int or max_iter <= 0:
        raise ValueError("positive integer iteration budget required")
    basis, transform, scales, rank = _coordinates(constraints)
    if rank != len(constraints):
        raise ValueError("constraints must be independent modulo identity; rank is unresolved")
    if basis[0].shape != sigma.shape:
        raise ValueError("target and constraints must use the same algebra")
    targets = _moments(sigma, basis)

    def evaluate(theta):
        rho, log_z, probs, vectors = _thermal(_combination(basis, theta))
        gradient = _moments(sigma-rho, basis)
        return log_z + float(theta @ targets), gradient, rho, probs, vectors

    theta = np.zeros(rank)
    for iteration in range(max_iter+1):
        base, gradient, _, probs, vectors = evaluate(theta)
        norm = math.hypot(*gradient)
        if norm < tol:
            with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                lam = np.linalg.solve(transform, theta)/scales
            result = projection_diagnostics(sigma, constraints, lam)
            if result.normalized_residual >= tol:
                raise RuntimeError("original multiplier replay does not converge at this precision")
            return ProjectionResult(result.multipliers, result.state,
                                    result.normalized_residual, result.raw_moment_residual,
                                    iteration, result.optimality_gap_bound,
                                    result.trace_distance_bound)
        if iteration == max_iter:
            break
        hessian = _covariance(basis, probs, vectors)
        spectrum = np.linalg.eigvalsh(hessian)
        if spectrum[0] <= 64*np.finfo(float).eps*rank*spectrum[-1]:
            raise ValueError("positive Duhamel Hessian is numerically unresolved")
        step = np.linalg.solve(hessian, -gradient)
        descent = float(gradient @ step)
        accepted = False
        for backtrack in range(60):
            factor = 2.**(-backtrack)
            try:
                candidate = evaluate(theta + factor*step)
            except ValueError:
                # A trial can leave the resolved faithful region; reduce it.
                continue
            value, new_gradient = candidate[:2]
            roundoff = 16*np.finfo(float).eps*max(1., abs(base))
            armijo = base+1e-4*factor*descent
            # Equality with base after rounding is not objective decrease.
            if ((value < base and value <= armijo)
                    or (abs(value-base) <= roundoff
                        and math.hypot(*new_gradient) < norm/2)):
                theta = theta + factor*step
                accepted = True
                break
        if not accepted:
            raise RuntimeError("information projection line search did not converge")
    raise RuntimeError("information projection did not converge within its budget")


def i_projection(sigma, constraints, tol=1e-11, max_iter=200):
    """Compatibility pair (multipliers, raw moment residual).

    Stopping uses the coordinate-invariant residual, so the returned raw norm
    need not be smaller than tol in arbitrary units. project_information also
    exposes the invariant residual and an optimum-error diagnostic.
    """
    result = project_information(sigma, constraints, tol, max_iter)
    return result.multipliers, result.raw_moment_residual

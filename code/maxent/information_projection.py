"""Finite Gibbs inference in the real observable space modulo the identity.

Unit changes, invertible real changes of basis and scalar energy offsets do
not change that space. Numerical rank and support must still be resolved.
See PROJECTION_COORDINATES.md for the proofs and precision boundary.
"""

from dataclasses import dataclass
import math

import numpy as np

from quantum_information import faithful_density_matrix, finite_real_scalar


def _finite(value, name):
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} exceeds finite numerical range")
    return value


def _numeric(value, name, *, real=False):
    # Check mixed Python containers before NumPy coerces True to 1.
    if isinstance(value, (list, tuple)):
        for part in value:
            _numeric(part, name, real=real)
    raw = np.asarray(value)
    if raw.dtype.kind not in ("iuf" if real else "iufc"):
        raise ValueError(f"{name} must contain finite numeric values, not Booleans")
    _finite(raw, name)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = raw.astype(float if real else complex)
    _finite(result, name)
    if (np.any((raw.real != 0) & (result.real == 0))
            or np.any((raw.imag != 0) & (result.imag == 0))):
        raise ValueError(f"{name} loses nonzero components in float64 conversion")
    return result


def _max_component(a):
    return max(float(np.max(np.abs(a.real))), float(np.max(np.abs(a.imag))))


def _scaled(a, scale):
    # Real division avoids the complex reciprocal overflowing for tiny units.
    result = a.real/scale + 1j*(a.imag/scale)
    if (np.any((a.real != 0) & (result.real == 0))
            or np.any((a.imag != 0) & (result.imag == 0))):
        raise ValueError("observable normalization loses nonzero components")
    return _finite(result, "observable normalization")


def _observables(constraints):
    """Return S_a = offset_a I + scale_a C_a with bounded traceless C_a.

Subtract a diagonal anchor before scaling: a huge scalar must neither erase
off-diagonal entries nor dominate the Hermiticity check on the variation.
"""
    if not isinstance(constraints, (list, tuple)) or not constraints:
        raise ValueError("nonempty constraint family required")
    centered, scales, offsets = [], [], []
    shape = None
    for value in constraints:
        a = _numeric(value, "constraints")
        if a.ndim != 2 or not a.shape[0] or a.shape[0] != a.shape[1]:
            raise ValueError("square constraint operators required")
        if shape is not None and a.shape != shape:
            raise ValueError("constraints must use the same algebra")
        shape = a.shape
        d = len(a)
        anchor = float(a[0, 0].real)
        b = a.copy()
        with np.errstate(over="ignore", invalid="ignore"):
            b[np.diag_indices(d)] -= anchor
        _finite(b, "observable differences")
        scale = _max_component(b)
        if scale == 0:
            centered.append(b)
            scales.append(1.)
            offsets.append(anchor)
            continue
        b = _scaled(b, scale)
        if np.linalg.norm(b-b.conj().T) > 1e-12:
            raise ValueError("constraints must be Hermitian")
        b = (b+b.conj().T)/2
        mean = float(np.trace(b).real/d)
        b[np.diag_indices(d)] -= mean
        centered.append(b)
        scales.append(scale)
        offsets.append(_finite(anchor + scale*mean, "observable offset"))
    return centered, np.array(scales), np.array(offsets)


def _multipliers(lam, count):
    result = _numeric(lam, "multipliers", real=True)
    if result.shape != (count,):
        raise ValueError("one finite real multiplier is required per constraint")
    return result


def _combination(operators, coefficients):
    with np.errstate(over="ignore", invalid="ignore"):
        terms = [c*a for c, a in zip(coefficients, operators)]
    for c, a, term in zip(coefficients, operators, terms):
        _finite(term, "constraint term")
        if c != 0 and (np.any((a.real != 0) & (term.real == 0))
                       or np.any((a.imag != 0) & (term.imag == 0))):
            raise ValueError("constraint term underflow; precision is insufficient")
    with np.errstate(over="ignore", invalid="ignore"):
        result = sum(terms)
    return _finite(result, "constraint combination")


def _hamiltonian_parts(constraints, lam):
    centered, scales, offsets = _observables(constraints)
    lam = _multipliers(lam, len(centered))
    with np.errstate(over="ignore", invalid="ignore"):
        coefficients = _finite(lam*scales, "scaled multipliers")
        terms = _finite(lam*offsets, "scalar energy")
    if (np.any((lam != 0) & (coefficients == 0))
            or np.any((lam != 0) & (offsets != 0) & (terms == 0))):
        raise ValueError("Hamiltonian coefficient underflow; precision is insufficient")
    try:
        shift = math.fsum(terms)
    except OverflowError as exc:
        raise ValueError("scalar energy exceeds finite numerical range") from exc
    return _combination(centered, coefficients), shift, centered, scales


def constrained_hamiltonian(constraints, lam):
    """Full Hamiltonian; Gibbs evaluation keeps its scalar part separate."""
    ham, shift, _, _ = _hamiltonian_parts(constraints, lam)
    with np.errstate(over="ignore", invalid="ignore"):
        return _finite(ham + shift*np.eye(len(ham)), "constraint combination")


def _thermal(ham):
    energies, vectors = np.linalg.eigh(ham)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        shifted = energies-energies[0]
        weights = np.exp(-shifted)
        probs = weights/weights.sum()
    if np.any(probs == 0) or not np.all(np.isfinite(probs)):
        raise ValueError("Gibbs spectrum underflow; faithful-state precision is insufficient")
    rho = faithful_density_matrix((vectors*probs) @ vectors.conj().T)
    return rho, math.log(weights.sum())-energies[0], probs, vectors


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

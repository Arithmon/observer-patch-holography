"""Construct and independently recognize finite state-preserving expectations.

A map is supplied as a complete linear superoperator, not a callable whose
linearity could be inferred from samples. vec(X) is row-major. The observable
map E is unital; its Hilbert-Schmidt adjoint E_* is the state channel.
Constructing E from a chosen algebra does not prove that a simulator uses E.
"""

import numpy as np

from .algebras import FiniteAlgebra, algebra_intersection, operator
from .states import faithful_density_matrix, finite_real_scalar


def _tolerance(value):
    value = finite_real_scalar(value, "numerical tolerance")
    if not 0 < value <= 1e-6:
        raise ValueError("numerical tolerance must lie in (0, 1e-6]")
    return value


def _data(algebra, reference):
    if not isinstance(algebra, FiniteAlgebra):
        raise ValueError("validated finite algebra required")
    rho = faithful_density_matrix(reference)
    if len(rho) != algebra.size:
        raise ValueError("reference and algebra must use the same operator space")
    return rho, np.kron(np.eye(algebra.size), rho.T)


def apply_map(superoperator, value, *, dual=False):
    a = operator(value)
    s = operator(superoperator, len(a)**2)
    if type(dual) is not bool:
        raise ValueError("dual must be Boolean")
    return ((s.conj().T if dual else s) @ a.reshape(-1)).reshape(a.shape)


def choi_matrix(superoperator, size):
    """Choi matrix in input tensor output order, without normalization by d."""
    from .states import dimensions
    size, = dimensions((size,))
    s = operator(superoperator, size*size)
    return s.reshape(size, size, size, size).transpose(2, 0, 3, 1).reshape(s.shape)


def gns_projection(algebra, reference):
    """Weighted orthogonal projection; it need not be positive or *-preserving.

    This deliberately exposes the algebraic candidate for negative controls.
    Use state_preserving_expectation to obtain a validated observable map.
    """
    _, weight = _data(algebra, reference)
    q = algebra._q
    gram = q.conj().T @ weight @ q
    spectrum = np.linalg.eigvalsh((gram + gram.conj().T)/2)
    if spectrum[0] <= 64*np.finfo(float).eps*len(gram)*spectrum[-1]:
        raise ValueError("weighted Gram matrix is numerically unresolved")
    result = q @ np.linalg.solve(gram, q.conj().T @ weight)
    return operator(result, algebra.size**2)


def expectation_diagnostics(superoperator, algebra, reference, *, tol=1e-9):
    """Recognize a supplied map by complete finite-matrix conditions.

    No comparison with the constructor's formula is used. CP, unitality,
    range in B, fixation of B, and preservation of rho are recomputed.
    Idempotence and weighted self-adjointness are also reported and checked.
    The ordinary trace-preservation defect is informational: E need not be
    trace preserving, whereas its adjoint must be a channel on states.
    """
    tol = _tolerance(tol)
    rho, weight = _data(algebra, reference)
    d, q = algebra.size, algebra._q
    s = operator(superoperator, d*d)
    choi = choi_matrix(s, d)
    norm = np.linalg.norm
    defects = {
        "choi_hermiticity": float(norm(choi - choi.conj().T)),
        "choi_min_eigenvalue": float(np.linalg.eigvalsh((choi+choi.conj().T)/2)[0]),
        "unitality": float(norm(apply_map(s, np.eye(d))-np.eye(d))),
        "range": float(norm(s - q @ (q.conj().T @ s))),
        "fixes_algebra": float(norm(s @ q - q)),
        "reference_preservation": float(norm(apply_map(s, rho, dual=True)-rho)),
        "idempotence": float(norm(s @ s - s)),
        "gns_self_adjointness": float(norm(weight @ s - s.conj().T @ weight)),
        "ordinary_trace_preservation": float(norm(apply_map(s, np.eye(d), dual=True)-np.eye(d))),
    }
    if not all(np.isfinite(value) for value in defects.values()):
        raise ValueError("map diagnostics exceed finite numerical range")
    required = ("choi_hermiticity", "unitality", "range", "fixes_algebra",
                "reference_preservation", "idempotence", "gns_self_adjointness")
    return {**defects, "tolerance": tol,
            "passed": (defects["choi_min_eigenvalue"] >= -tol
                       and all(defects[key] <= tol for key in required))}


def state_preserving_expectation(algebra, reference, *, tol=1e-9):
    """Unique reference-preserving expectation, rejecting a nonexistent one."""
    tol = _tolerance(tol)
    _data(algebra, reference)
    if algebra.modular_invariance_defect(reference) > tol:
        raise ValueError("reference does not preserve the algebra under modular flow")
    result = gns_projection(algebra, reference)
    if not expectation_diagnostics(result, algebra, reference, tol=tol)["passed"]:
        raise ValueError("conditional-expectation conditions are numerically unresolved")
    return result


def repair_generator(superoperators, algebras, reference, rates, *, tol=1e-9):
    """Recognize primitive maps and form L=sum rate_m (E_m-I), observables.

    Returns the generator and a numerical intersection/gap diagnostic. The
    finite theorem needs positive rates but does not require commuting E_m.
    Rates are supplied data; this does not select their values or a clock.
    """
    tol = _tolerance(tol)
    if (not isinstance(superoperators, (list, tuple)) or not superoperators
            or not isinstance(algebras, (list, tuple))
            or not isinstance(rates, (list, tuple, np.ndarray))
            or (isinstance(rates, np.ndarray) and rates.ndim != 1)
            or len(superoperators) != len(algebras) or len(rates) != len(algebras)):
        raise ValueError("nonempty matched maps, algebras and rates required")
    rates = [finite_real_scalar(rate, "repair rate") for rate in rates]
    if any(rate <= 0 for rate in rates):
        raise ValueError("every primitive repair rate must be positive")
    common = algebra_intersection(algebras)
    rho, weight = _data(common, reference)
    d = common.size
    maps = [operator(s, d*d) for s in superoperators]
    for s, algebra in zip(maps, algebras):
        if not expectation_diagnostics(s, algebra, rho, tol=tol)["passed"]:
            raise ValueError("supplied primitive is not a reference-preserving expectation")
    identity = np.eye(d*d)
    generator = operator(sum(rate*(s-identity) for rate, s in zip(rates, maps)), d*d)
    target = state_preserving_expectation(common, rho, tol=tol)
    # R^* R = W turns the weighted Dirichlet form into an ordinary Hermitian
    # matrix. Its kernel dimension is supplied by the actual intersection,
    # rather than inferred by deleting arbitrary small positive eigenvalues.
    r = np.linalg.cholesky(weight).conj().T
    symmetric = operator(np.linalg.solve(r.T, (-r @ generator).T).T, d*d)
    hermitian_defect = float(np.linalg.norm(symmetric-symmetric.conj().T))
    spectrum = np.linalg.eigvalsh((symmetric+symmetric.conj().T)/2)
    count = common.dimension
    kernel_defect = float(np.max(np.abs(spectrum[:count])))
    stationarity = float(np.linalg.norm(apply_map(generator, rho, dual=True)))
    target_residuals = [float(np.linalg.norm(generator @ target)),
                        float(np.linalg.norm(target @ generator))]
    residuals = [hermitian_defect, kernel_defect, stationarity, *target_residuals]
    if not all(np.isfinite(value) for value in residuals):
        raise ValueError("repair diagnostics exceed finite numerical range")
    target_defect = max(target_residuals)
    if max(residuals) > tol:
        raise ValueError("repair generator identities are numerically unresolved")
    gap = None if count == d*d else float(spectrum[count])
    if gap is not None and gap <= 64*np.finfo(float).eps*d*d*max(1.,float(spectrum[-1])):
        raise ValueError("positive repair gap is numerically unresolved")
    return generator, {
        "intersection_dimension": count,
        "complement_dimension": d*d-count,
        "gap": gap,
        "weighted_hermiticity_defect": hermitian_defect,
        "kernel_defect": kernel_defect,
        "stationarity_defect": stationarity,
        "intersection_fixation_defect": target_defect,
        "tolerance": tol,
    }

"""Validated finite null tomography; see NULL_TOMOGRAPHY.md for the contract.

Charges are quadratic in the supplied ray. Fits and error budgets use the
representative k/k[0], in a fixed observer frame. A fit is a finite linear-
algebra diagnostic, not a derivation of local stress or unsampled charges.
"""

from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np

ETA = np.diag([-1., 1., 1., 1.])
_EPS = np.finfo(float).eps


def _real_array(value, name):
    try:
        original = np.asarray(value)
        if original.dtype.kind not in "iuf" or not np.all(np.isfinite(original)):
            raise ValueError
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            result = original.astype(float)
        if (not np.all(np.isfinite(result))
                or np.any((original != 0) & (result == 0))):
            raise ValueError
        # Comparisons after NumPy's int-to-float promotion can hide lost bits.
        if original.dtype.kind in "iu":
            if any(int(x) != int(y) for x, y in zip(original.flat, result.flat)):
                raise ValueError
        elif np.any(result.astype(original.dtype) != original):
            raise ValueError
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain exactly representable finite real binary64 data") from exc
    return result


def _finite(value, name):
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} is outside binary64 range")
    return value


def _norm(value):
    # Unlike sqrt(dot(x,x)), hypot does not square tiny residuals to zero.
    return float(_finite(math.hypot(*np.ravel(value)), "norm"))


def _budget(value):
    value = _real_array(value, "error budget")
    if value.ndim != 0 or value < 0:
        raise ValueError("error budget must be a finite nonnegative scalar")
    return float(value)


def _rounded_fraction(value):
    """Round once; do not erase a nonzero result or return infinity."""
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError("quadratic charge is outside binary64 range") from exc
    if not math.isfinite(result) or (value != 0 and result == 0):
        raise ValueError("quadratic charge is outside binary64 range")
    return result


def null_vector(direction):
    """Return (1,n) for a finite, nonzero spatial direction of any scale."""
    d = _real_array(direction, "direction")
    if d.shape != (3,) or not np.any(d):
        raise ValueError("direction must be a nonzero three-vector")
    scaled = d / np.max(np.abs(d))
    normalized = scaled / _norm(scaled)
    if np.any((d != 0) & (normalized == 0)):
        raise ValueError("direction dynamic range is unresolved in binary64")
    return np.r_[1., normalized]


def _rays(null_dirs):
    rays = _real_array(null_dirs, "null rays")
    if rays.ndim != 2 or rays.shape[1] != 4 or len(rays) == 0:
        raise ValueError("null rays must be a nonempty (n,4) array")
    times = rays[:, 0]
    if np.any(times == 0):
        raise ValueError("each null ray must have a nonzero time component")
    with np.errstate(over="ignore", invalid="ignore"):
        normalized = rays / times[:, None]
    _finite(normalized, "normalized null rays")
    if np.any((rays != 0) & (normalized == 0)):
        raise ValueError("null ray dynamic range is unresolved in binary64")
    for k in normalized:
        if abs(_norm(k[1:]) - 1.) > 64 * _EPS:
            raise ValueError("each ray must be Minkowski-null")
    return normalized, times


def _tensor(t_matrix):
    t = _real_array(t_matrix, "tensor")
    if t.shape != (4, 4) or not np.array_equal(t, t.T):
        raise ValueError("tensor must be a symmetric real 4x4 matrix")
    return t


def sym_basis():
    """Legacy ten-coordinate basis of symmetric 4x4 matrices."""
    basis = []
    for i in range(4):
        for j in range(i, 4):
            b = np.zeros((4, 4))
            b[i, j] = b[j, i] = 1.
            basis.append(b)
    return basis


def tracefree_basis():
    """Frobenius-orthonormal basis of eta-trace-free symmetric tensors."""
    basis = [np.diag([3., 1., 1., 1.]) / math.sqrt(12),
             np.diag([0., 1., -1., 0.]) / math.sqrt(2),
             np.diag([0., 1., 1., -2.]) / math.sqrt(6)]
    for i in range(4):
        for j in range(i + 1, 4):
            b = np.zeros((4, 4))
            b[i, j] = b[j, i] = 1 / math.sqrt(2)
            basis.append(b)
    return np.array(basis)


def _design(rays, basis):
    return np.einsum("ni,aij,nj->na", rays, basis, rays)


def design_matrix(null_dirs):
    """Raw quadratic rows in sym_basis coordinates, after null validation.

    For very large/small ray weights these rows may be unrepresentable; the
    normalized solver does not construct them.
    """
    _rays(null_dirs)
    rows = []
    for ray in np.asarray(null_dirs, dtype=float):
        k = list(map(Fraction, ray))
        rows.append([_rounded_fraction(k[i] * k[j] * (1 if i == j else 2))
                     for i in range(4) for j in range(i, 4)])
    return np.array(rows)


def charges_of(t_matrix, null_dirs):
    """Raw T(k,k) charges; reject malformed tensors and unrepresentable data."""
    t = _tensor(t_matrix)
    _rays(null_dirs)
    # Sixteen exact rational products are cheap at this finite audit scale.
    # A common floating scale would erase small terms surviving cancellation,
    # e.g. T00=-T22=1e200, T02=1e-200, k=(1,0,1,0).
    exact_t = [[Fraction(x) for x in row] for row in t]
    charges = []
    for ray in np.asarray(null_dirs, dtype=float):
        k = list(map(Fraction, ray))
        q = sum(k[i] * exact_t[i][j] * k[j] for i in range(4) for j in range(4))
        charges.append(_rounded_fraction(q))
    return np.array(charges)


def eta_project_out(t_matrix):
    """Remove the metric component in the Frobenius inner product."""
    t = _tensor(t_matrix)
    # Divide before summing to avoid overflow in the trace.
    coefficient = math.fsum(np.diag(t) * np.diag(ETA) / 4)
    with np.errstate(over="ignore"):
        return _finite(t - coefficient * ETA, "trace-free tensor")


def tomography_directions():
    """The six axes and three body diagonals already proved in Tensor.lean."""
    spatial = [np.eye(3)[i] * sign for i in range(3) for sign in (1, -1)]
    spatial += [np.array(s) for s in ((1, 1, 1), (1, 1, -1), (1, -1, 1))]
    return np.array([null_vector(n) for n in spatial])


@dataclass(frozen=True)
class NullTomographyFit:
    """Finite normalized least-squares diagnostic, never a physical source.

    Floating-point singular values/bounds are estimates, not interval proofs.
    witness is the unit residual direction (zero if the residual is zero).
    Its design defect reports how nearly it satisfies A.T @ witness = 0.
    """

    tensor: np.ndarray
    residual_norm: float
    singular_values: np.ndarray
    witness: np.ndarray
    witness_design_defect: float
    witness_charge: float

    @property
    def noise_amplification(self):
        """Sharp Frobenius tensor / normalized-charge Euclidean norm gain."""
        return 1. / float(self.singular_values[-1])

    def tensor_error_bound(self, charge_error_bound):
        """||T_hat - T_true||_F <= delta/s_min for true charges within delta.

        Both tensors are in the eta-trace-free gauge; geometry is fixed.
        This data-error bound does not include floating-point solve error.
        """
        return _rounded_fraction(Fraction(_budget(charge_error_bound))
                                  / Fraction(float(self.singular_values[-1])))

    def require_consistent(self, *, error_budget):
        """Require distance to the sampled tensor image <= explicit budget.

        Budget includes any measurement and numerical allowance. There is
        no hidden absolute tolerance and no inference about unsampled rays.
        """
        if self.residual_norm > _budget(error_budget):
            raise ValueError("charges violate a dependent-family relation beyond the error budget")
        return self.tensor.copy()


def fit_null_charges(charges, null_dirs):
    """Fit in nine trace-free coordinates; reject unresolved designs.

    Raw inputs transform as (k,q)->(c*k,c**2*q). All diagnostics instead use
    k/k[0] and q/k[0]**2, making them independent of these ray representatives.
    Rank uses 64*eps*max(n,9)*s_max. No rank-deficient pseudoinverse is returned.
    """
    rays, times = _rays(null_dirs)
    q = _real_array(charges, "charges")
    if q.shape != (len(rays),):
        raise ValueError("one scalar charge is required per null ray")
    if len(rays) < 9:
        raise ValueError("tomography requires at least nine independent null rays")
    q = np.array([_rounded_fraction(Fraction(x) / Fraction(t)**2)
                  for x, t in zip(q, times)])
    basis = tracefree_basis()
    a = _design(rays, basis)
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    if s[-1] <= 64 * _EPS * max(a.shape) * s[0]:
        raise ValueError("null tomography design has unresolved rank below nine")
    scale = float(np.max(np.abs(q)))
    y = q / scale if scale else q
    if np.any((q != 0) & (y == 0)):
        raise ValueError("charge dynamic range is unresolved in binary64")
    coefficients = vt.T @ ((u.T @ y) / s)
    residual = y - a @ coefficients
    residual_size = _norm(residual)
    witness = residual / residual_size if residual_size else np.zeros(len(q))
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        tensor_unit = np.einsum("a,aij->ij", coefficients, basis)
        tensor = tensor_unit * scale
        residual_norm = residual_size * scale
        witness_charge = float(witness @ y) * scale
    _finite(tensor, "reconstructed tensor")
    _finite([residual_norm, witness_charge], "residual diagnostics")
    if (np.any((tensor_unit != 0) & (tensor == 0)) and scale
            or (residual_size and scale and residual_norm == 0)):
        raise ValueError("tomography result is outside binary64 range")
    return NullTomographyFit(tensor, float(residual_norm), s, witness,
                             _norm(a.T @ witness), witness_charge)


def reconstruct_from_charges(charges, null_dirs):
    """Compatibility diagnostic: (tensor, normalized residual), not acceptance.

    Use fit_null_charges(...).require_consistent(error_budget=...) to enforce
    consistency. Nonunit ray representatives now use normalized residuals.
    """
    fit = fit_null_charges(charges, null_dirs)
    return fit.tensor, fit.residual_norm

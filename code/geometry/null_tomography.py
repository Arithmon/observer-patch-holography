"""Validated finite null tomography; see NULL_TOMOGRAPHY.md for the contract.

Charges are quadratic in the supplied ray. Fits and error budgets use the
representative k/k[0], in a fixed observer frame. A fit is a finite linear-
algebra diagnostic, not a derivation of local stress or unsampled charges.
"""

from dataclasses import dataclass, field
from fractions import Fraction
import math

import numpy as np

ETA = np.diag([-1., 1., 1., 1.])
_EPS = np.finfo(float).eps


def _real_array(value, name):
    try:
        # Inspect elements BEFORE NumPy promotes mixed integer/float sequences.
        # Otherwise [2**60, 2**60+1, 0.] has already lost its differing bit.
        original = np.asarray(value, dtype=object)
        if any(isinstance(x, (bool, np.bool_)) or not isinstance(
                x, (int, float, np.integer, np.floating)) for x in original.flat):
            raise ValueError
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            result = original.astype(float)
        if not np.all(np.isfinite(result)):
            raise ValueError
        for x, y in zip(original.flat, result.flat):
            if isinstance(x, (int, np.integer)):
                if int(x) != int(y):
                    raise ValueError
            elif x != y:
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


def _sqrt_fraction(value):
    """Scale an exact squared norm before taking its numerical square root."""
    if not value:
        return 0.
    exponent = (value.numerator.bit_length() - value.denominator.bit_length()) // 2
    scale = Fraction(2)**exponent
    return _rounded_fraction(Fraction(math.sqrt(float(value / scale**2))) * scale)


def _quadratic(t, k):
    return sum(k[i] * t[i][j] * k[j] for i in range(4) for j in range(4))


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
    return normalized, times, rays


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
        q = _quadratic(exact_t, k)
        charges.append(_rounded_fraction(q))
    return np.array(charges)


def eta_project_out(t_matrix):
    """Remove the metric component in the Frobenius inner product."""
    t = _tensor(t_matrix)
    # Dividing the summands first would erase a subnormal metric coefficient.
    signs = (-1, 1, 1, 1)
    coefficient = sum(Fraction(t[i, i]) * signs[i] for i in range(4)) / 4
    return np.array([[_rounded_fraction(Fraction(t[i, j]) - (
        coefficient * signs[i] if i == j else 0)) for j in range(4)] for i in range(4)])


def tomography_directions():
    """The six axes and three body diagonals already proved in Tensor.lean."""
    spatial = [np.eye(3)[i] * sign for i in range(3) for sign in (1, -1)]
    spatial += [np.array(s) for s in ((1, 1, 1), (1, 1, -1), (1, -1, 1))]
    return np.array([null_vector(n) for n in spatial])


@dataclass(frozen=True)
class NullTomographyFit:
    """Finite normalized least-squares diagnostic, never a physical source.

    Singular values/bounds are estimates. Acceptance separately replays the
    returned tensor with exact rational arithmetic on the original inputs.
    witness is the numerical unit residual direction (zero for zero residual).
    Its design defect reports how nearly it satisfies A.T @ witness = 0.
    """

    tensor: np.ndarray
    residual_norm: float
    singular_values: np.ndarray
    witness: np.ndarray
    witness_design_defect: float
    witness_charge: float
    _residual_squared: Fraction = field(repr=False)

    def __post_init__(self):
        # A frozen dataclass alone leaves its arrays writable. Immutable bytes
        # prevent both in-place edits and re-enabling NumPy's WRITEABLE flag.
        for name in ("tensor", "singular_values", "witness"):
            array = getattr(self, name)
            readonly = np.frombuffer(array.tobytes(), dtype=float).reshape(array.shape)
            object.__setattr__(self, name, readonly)

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
        """Replay the returned tensor against original inputs within the budget.

        Compare exact rational squared residual and budget, including solve
        and normalization roundoff. A failure can mean numerical resolution
        is insufficient; it alone does not prove nonexistence of a better fit.
        """
        if self._residual_squared > Fraction(_budget(error_budget))**2:
            raise ValueError("returned tensor exceeds the consistency error budget")
        return self.tensor.copy()


def fit_null_charges(charges, null_dirs):
    """Fit in nine trace-free coordinates; reject unresolved designs.

    Raw inputs transform as (k,q)->(c*k,c**2*q). All diagnostics instead use
    k/k[0] and q/k[0]**2, making them independent of these ray representatives.
    Rank uses 64*eps*max(n,9)*s_max. No rank-deficient pseudoinverse is returned.
    """
    rays, times, raw_rays = _rays(null_dirs)
    raw_q = _real_array(charges, "charges")
    if raw_q.shape != (len(rays),):
        raise ValueError("one scalar charge is required per null ray")
    if len(rays) < 9:
        raise ValueError("tomography requires at least nine independent null rays")
    exact_q = [Fraction(x) / Fraction(t)**2 for x, t in zip(raw_q, times)]
    q = np.array([_rounded_fraction(x) for x in exact_q])
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
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        tensor_unit = np.einsum("a,aij->ij", coefficients, basis)
        tensor = tensor_unit * scale
    _finite(tensor, "reconstructed tensor")
    if scale and np.any((tensor_unit != 0) & (tensor == 0)):
        raise ValueError("tomography result is outside binary64 range")

    # Acceptance concerns the returned binary tensor and ORIGINAL raw data,
    # not an intermediate factorization or rounded normalized inputs.
    exact_t = [[Fraction(x) for x in row] for row in tensor]
    residual = []
    for reading, ray in zip(raw_q, raw_rays):
        k = list(map(Fraction, ray))
        residual.append((Fraction(reading) - _quadratic(exact_t, k)) / k[0]**2)
    squared = sum((r*r for r in residual), Fraction(0))
    residual_norm = _sqrt_fraction(squared)
    witness = (np.array([float(r / Fraction(residual_norm)) for r in residual])
               if residual_norm else np.zeros(len(q)))
    witness_charge = _rounded_fraction(sum(Fraction(w)*x for w, x in zip(witness, exact_q)))
    return NullTomographyFit(tensor, float(residual_norm), s, witness,
                             _norm(a.T @ witness), witness_charge, squared)


def reconstruct_from_charges(charges, null_dirs):
    """Compatibility diagnostic: (tensor, normalized residual), not acceptance.

    Use fit_null_charges(...).require_consistent(error_budget=...) to enforce
    consistency. Nonunit ray representatives now use normalized residuals.
    """
    fit = fit_null_charges(charges, null_dirs)
    return fit.tensor.copy(), fit.residual_norm

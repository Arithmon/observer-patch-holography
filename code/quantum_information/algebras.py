"""Finite operator algebras, explicit tensor slots and separation diagnostics.

Matrices use ordinary trace and row-major vectorization. Numerical rank
decisions are diagnostics; unresolved support or basis rank raises.
"""

from math import prod

import numpy as np

from .gibbs import _numeric
from .states import (
    _indices, _spectrum, _unresolved_positive_spectrum,
    dimensions, faithful_log, partial_trace,
)

ALGEBRA_ATOL = 1e-10


def operator(value, size=None):
    raw = _numeric(value, "operator")
    if (raw.ndim != 2 or raw.shape[0] == 0 or raw.shape[0] != raw.shape[1]
            or raw.dtype.kind not in "iufc"):
        raise ValueError("nonempty square numeric operator required")
    a = raw
    if not np.all(np.isfinite(a)) or (size is not None and a.shape != (size, size)):
        raise ValueError("finite operator on the declared space required")
    return a


def orthonormal_operator_basis(basis):
    """HS-orthonormalize the supplied independent binary64 operator span.

    A well-conditioned floating chart is cheap. Ill-conditioned presentations
    are preconditioned over the original exact entries, before roundoff can discard
    their independent directions. Unresolved normalization or conversion raises.
    """
    if not isinstance(basis, (list, tuple)) or not basis:
        raise ValueError("nonempty operator basis required")
    matrices = [operator(b) for b in basis]
    d = len(matrices[0])
    if any(b.shape != (d, d) for b in matrices):
        raise ValueError("basis operators must act on one space")
    if len(basis) > d*d:
        raise ValueError("operator basis is dependent")
    from .operator_spans import numerical_chart
    columns = np.column_stack([b.reshape(-1) for b in matrices])
    q = numerical_chart(columns)
    if q is None:
        from .operator_spans import exact_columns, numerical_basis
        return numerical_basis(exact_columns(matrices), d)
    return [q[:, i].reshape(d, d) for i in range(len(basis))]


def embed_operator(value, dims, region):
    """Embed in any subset of tensor slots; retained slots use original order."""
    dims = dimensions(dims)
    region = _indices(region, len(dims))
    rest = tuple(i for i in range(len(dims)) if i not in region)
    a = operator(value, prod(dims[i] for i in region))
    grouped = region + rest
    inverse = tuple(int(i) for i in np.argsort(grouped))
    grouped_dims = tuple(dims[i] for i in grouped)
    full = np.kron(a, np.eye(prod(dims[i] for i in rest)))
    tensor = full.reshape(grouped_dims + grouped_dims)
    order = inverse + tuple(i + len(dims) for i in inverse)
    return tensor.transpose(order).reshape(prod(dims), prod(dims))


def matrix_units(size):
    size, = dimensions((size,))
    units = []
    for i in range(size):
        for j in range(size):
            unit = np.zeros((size, size), complex)
            unit[i, j] = 1
            units.append(unit)
    return units


def separation_modulus(basis, omega):
    """inf ||a omega|| over a in span(basis) with ambient HS norm one.

    Orthonormalizing the operator span solves the generalized Gram problem
    H c = lambda G c. Nonorthogonal and unequally scaled bases are allowed;
    dependent bases are rejected. No algebra closure is implied by this API.
    General singular values are numerical, not interval bounds. Exact source
    arithmetic resolves singleton actions and kernels; unresolved small
    positive multi-column actions raise rather than becoming a false zero.
    """
    orthogonal = orthonormal_operator_basis(basis)
    raw = _numeric(omega, "reference vector")
    if (raw.shape != (len(orthogonal[0]),) or raw.dtype.kind not in "iufc"
            or not np.all(np.isfinite(raw))):
        raise ValueError("finite reference vector on the operator space required")
    # Validate in the arithmetic used by the action. Integer dot products can
    # wrap to one, and low-precision dot products can erase a nonunit norm.
    raw = raw.astype(complex)
    squared_norm = np.vdot(raw, raw)
    if not np.isfinite(squared_norm) or abs(squared_norm - 1) > 1e-12:
        raise ValueError("reference vector must be normalized")
    if len(orthogonal) == 1:
        from .operator_spans import small_separation_modulus
        return small_separation_modulus([operator(b) for b in basis], raw)
    action = np.column_stack([b @ raw for b in orthogonal])
    if action.shape[1] > action.shape[0]:
        return 0.0
    singular = np.linalg.svd(action, compute_uv=False)
    # Unit HS columns and a unit vector bound each action column by one.
    # An absolute chart-error scale also detects an entirely rounded kernel.
    floor = 64*np.finfo(float).eps*max(action.shape)
    if singular[-1] <= floor:
        # Chart and action roundoff cannot certify a true zero. Resolve a kernel
        # from the original entries, without turning small positives into zero.
        from .operator_spans import small_separation_modulus
        return small_separation_modulus([operator(b) for b in basis], raw)
    return float(singular[-1])


def resolved_state_spectrum(rho):
    a, eigenvalues, vectors = _spectrum(rho)
    if _unresolved_positive_spectrum(a, eigenvalues):
        raise ValueError("state support is numerically unresolved")
    kernel = vectors[:, eigenvalues == 0]
    if kernel.size and np.any(a @ kernel != 0):
        raise ValueError("state support is numerically unresolved")
    return eigenvalues


def tensor_separation_modulus(rho, dims, region):
    """Separation modulus of L(M_region tensor I_rest) on HS(H).

    Its square is lambda_min(rho_region)/(dim(rest)*dim(H)), because
    ||L(a_region tensor I_rest)||_HS^2 = dim(H)*dim(rest)*||a_region||_HS^2.
    This computes no dim(H)^2-by-dim(H)^2 representation matrices.
    """
    dims = dimensions(dims)
    region = _indices(region, len(dims))
    reduced = partial_trace(rho, dims, region)
    smallest = resolved_state_spectrum(reduced)[0]
    total = prod(dims)
    spectator = total // len(reduced)
    return float(np.sqrt(smallest) / np.sqrt(total * spectator))


class FiniteAlgebra:
    """A validated unital complex *-subalgebra of M_d supplied by a basis.

    Closure is checked on every pair of basis elements, not random probes.
    The tolerance concerns numerical closure, not a theorem of exact rank.
    """

    def __init__(self, basis):
        if not isinstance(basis, (list, tuple)) or not basis:
            raise ValueError("nonempty operator basis required")
        supplied = tuple(operator(b) for b in basis)
        normalized = orthonormal_operator_basis(list(supplied))
        self.size = len(normalized[0])
        self.dimension = len(normalized)
        self._q = np.column_stack([b.reshape(-1) for b in normalized])
        # Closure checks already query the original span. Intersections can
        # later replace its exact columns without retaining a stale projector.
        self._supplied = supplied
        self._exact = None
        self._projector = None
        if not self._project(np.eye(self.size), within=ALGEBRA_ATOL*np.sqrt(self.size)):
            raise ValueError("algebra must contain the identity")
        for a in normalized:
            if not self._project(a.conj().T, within=ALGEBRA_ATOL):
                raise ValueError("algebra must be closed under adjoint")
            for b in normalized:
                if not self._project(a @ b, within=ALGEBRA_ATOL):
                    raise ValueError("algebra must be closed under multiplication")
        self._q.setflags(write=False)
        # A numerical chart is not evidence of exact linear dependence.
        # Keep the supplied span separately for conservation decisions.
        for b in supplied:
            b.setflags(write=False)

    def _exact_columns(self):
        from .operator_spans import exact_columns
        if self._exact is None:
            self._exact = exact_columns(self._supplied)
        return self._exact

    @property
    def basis(self):
        return [self._q[:, i].reshape(self.size, self.size).copy()
                for i in range(self.dimension)]

    def project(self, value):
        """HS projection onto the original span; not generally rho-preserving."""
        a = operator(value, self.size)
        return self._project(a)

    def distance(self, value):
        """HS distance before projection rounding; refuse final range loss."""
        a = operator(value, self.size)
        return self._project(a, distance=True)

    def _project(self, value, *, distance=False, within=None):
        if self.dimension == self.size**2:
            if within is not None:
                return True
            return 0. if distance else value.copy()
        from .linear_actions import SpanProjector
        columns = self._exact_columns()
        if self._projector is None or self._projector.columns is not columns:
            self._projector = SpanProjector(columns)
        return self._projector.evaluate(value, distance=distance, within=within)

    def modular_invariance_defect(self, rho):
        k = faithful_log(operator(rho, self.size))
        return max(self.distance(k @ b - b @ k) for b in self.basis)


def tensor_factor_algebra(dims, region):
    dims = dimensions(dims)
    region = _indices(region, len(dims))
    return FiniteAlgebra([embed_operator(b, dims, region)
                          for b in matrix_units(prod(dims[i] for i in region))])


def algebra_intersection(algebras):
    """Intersect the supplied spans over exact binary64 rational entries.

    The returned numerical basis must still be resolved and algebra-closed.
    Rounded presentations of the same ideal algebra may have a different
    exact intersection; no tolerance enlarges it into conserved observables.
    """
    if (not isinstance(algebras, (list, tuple)) or not algebras
            or any(not isinstance(a, FiniteAlgebra) for a in algebras)):
        raise ValueError("nonempty family of finite algebras required")
    if any(a.size != algebras[0].size for a in algebras):
        raise ValueError("intersection requires a common operator space")
    from .operator_spans import contains_identity, intersect_columns, numerical_basis
    common = algebras[0]._exact_columns()
    for algebra in algebras[1:]:
        common = intersect_columns(common, algebra._exact_columns())
    d = algebras[0].size
    if not contains_identity(common, d):
        raise ValueError("supplied spans have no exact common identity; intersection is unresolved")
    result = FiniteAlgebra(numerical_basis(common, d))
    # Preserve exact output across repeated intersections; converting this
    # rational basis to its floating chart must not change later rank tests.
    result._exact = common
    return result

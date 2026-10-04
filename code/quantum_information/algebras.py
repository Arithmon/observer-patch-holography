"""Finite operator algebras, explicit tensor slots and separation diagnostics.

Matrices use ordinary trace and row-major vectorization. Numerical rank
decisions are diagnostics; unresolved support or basis rank raises.
"""

from math import prod

import numpy as np

from .states import (
    _indices, _spectrum, _unresolved_positive_spectrum,
    dimensions, faithful_log, partial_trace,
)

ALGEBRA_ATOL = 1e-10


def operator(value, size=None):
    raw = np.asarray(value)
    if (raw.ndim != 2 or raw.shape[0] == 0 or raw.shape[0] != raw.shape[1]
            or raw.dtype.kind not in "iufc"):
        raise ValueError("nonempty square numeric operator required")
    a = raw.astype(complex)
    if not np.all(np.isfinite(a)) or (size is not None and a.shape != (size, size)):
        raise ValueError("finite operator on the declared space required")
    return a


def orthonormal_operator_basis(basis):
    """HS-orthonormalize an independent basis, without assuming equal norms."""
    if not isinstance(basis, (list, tuple)) or not basis:
        raise ValueError("nonempty operator basis required")
    matrices = [operator(b) for b in basis]
    d = len(matrices[0])
    if any(b.shape != (d, d) for b in matrices):
        raise ValueError("basis operators must act on one space")
    columns = np.column_stack([b.reshape(-1) for b in matrices])
    norms = np.linalg.norm(columns, axis=0)
    if np.any(norms == 0) or not np.all(np.isfinite(norms)):
        raise ValueError("nonzero finite basis norms required")
    columns = columns / norms
    singular = np.linalg.svd(columns, compute_uv=False)
    floor = 64*np.finfo(float).eps*max(columns.shape)*singular[0]
    if len(basis) > d*d or singular[-1] <= floor:
        raise ValueError("operator basis is dependent or numerically unresolved")
    q, _ = np.linalg.qr(columns, mode="reduced")
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
    """
    orthogonal = orthonormal_operator_basis(basis)
    raw = np.asarray(omega)
    if (raw.shape != (len(orthogonal[0]),) or raw.dtype.kind not in "iufc"
            or not np.all(np.isfinite(raw))):
        raise ValueError("finite reference vector on the operator space required")
    squared_norm = np.vdot(raw, raw)
    if not np.isfinite(squared_norm) or abs(squared_norm - 1) > 1e-12:
        raise ValueError("reference vector must be normalized")
    action = np.column_stack([b @ raw for b in orthogonal])
    if action.shape[1] > action.shape[0]:
        return 0.0
    return float(np.linalg.svd(action, compute_uv=False)[-1])


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
        normalized = orthonormal_operator_basis(basis)
        self.size = len(normalized[0])
        self.dimension = len(normalized)
        self._q = np.column_stack([b.reshape(-1) for b in normalized])
        if self.distance(np.eye(self.size)) > ALGEBRA_ATOL*np.sqrt(self.size):
            raise ValueError("algebra must contain the identity")
        for a in normalized:
            if self.distance(a.conj().T) > ALGEBRA_ATOL:
                raise ValueError("algebra must be closed under adjoint")
            for b in normalized:
                if self.distance(a @ b) > ALGEBRA_ATOL:
                    raise ValueError("algebra must be closed under multiplication")
        self._q.setflags(write=False)

    @property
    def basis(self):
        return [self._q[:, i].reshape(self.size, self.size).copy()
                for i in range(self.dimension)]

    def project(self, value):
        """Ordinary-trace HS projection; not generally rho-preserving."""
        a = operator(value, self.size)
        return (self._q @ (self._q.conj().T @ a.reshape(-1))).reshape(a.shape)

    def distance(self, value):
        a = operator(value, self.size)
        return float(np.linalg.norm(a - self.project(a), ord="fro"))

    def modular_invariance_defect(self, rho):
        k = operator(faithful_log(rho), self.size)
        return max(self.distance(k @ b - b @ k) for b in self.basis)


def tensor_factor_algebra(dims, region):
    dims = dimensions(dims)
    region = _indices(region, len(dims))
    return FiniteAlgebra([embed_operator(b, dims, region)
                          for b in matrix_units(prod(dims[i] for i in region))])


def algebra_intersection(algebras):
    if (not isinstance(algebras, (list, tuple)) or not algebras
            or any(not isinstance(a, FiniteAlgebra) for a in algebras)):
        raise ValueError("nonempty family of finite algebras required")
    if any(a.size != algebras[0].size for a in algebras):
        raise ValueError("intersection requires a common operator space")
    q = algebras[0]._q.copy()
    for algebra in algebras[1:]:
        residual = q - algebra._q @ (algebra._q.conj().T @ q)
        _, singular, vh = np.linalg.svd(residual, full_matrices=False)
        # Closure tolerance is not a rank cutoff. In particular, two distinct
        # almost coincident MASAs must not acquire a spurious large gap by
        # silently enlarging their intersection at the looser closure scale.
        rank_roundoff = 64*np.finfo(float).eps*max(residual.shape)
        q = q @ vh.conj().T[:, singular <= rank_roundoff]
    d = algebras[0].size
    return FiniteAlgebra([q[:, i].reshape(d, d) for i in range(q.shape[1])])

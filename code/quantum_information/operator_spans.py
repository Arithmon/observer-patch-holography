"""Exact supplied operator spans, separate from their numerical QR charts."""

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def unit_columns(columns):
    """Normalize complex columns without squaring their supplied units."""
    scales = np.maximum(np.max(np.abs(columns.real), axis=0),
                        np.max(np.abs(columns.imag), axis=0))
    if np.any(scales == 0):
        raise ValueError("nonzero basis operators required")
    result = np.empty_like(columns)
    # Real component division works even when a complex reciprocal would
    # overflow at subnormal scales. Check both conversion steps for lost data.
    with np.errstate(under="ignore"):
        result.real = columns.real/scales
        result.imag = columns.imag/scales
        norms = np.linalg.norm(result, axis=0)
        result.real /= norms
        result.imag /= norms
    if (np.any((columns.real != 0) & (result.real == 0))
            or np.any((columns.imag != 0) & (result.imag == 0))):
        raise ValueError("operator chart loses nonzero components; span is unresolved")
    return result


def exact_columns(basis):
    entries = sp.Matrix.hstack(*[
        sp.Matrix([sp.Rational(float(x.real))+sp.I*sp.Rational(float(x.imag))
                   for x in b.flat]) for b in basis])
    return DomainMatrix.from_Matrix(entries).convert_to(sp.QQ_I)


def _inner(left, right):
    return sum((sp.QQ_I(x.x, -x.y)*y for x, y in zip(left, right)), sp.QQ_I.zero)


def orthogonal_columns(columns):
    """Hermitian-orthogonal rational directions of an independent source span.

    Canonical row reduction removes supplied units and shears. Gram-Schmidt
    then uses exact QQ(i) products, so conversion to floats happens only after
    cancellation and conditioning in the supplied presentation are removed.
    """
    reduced, pivots = columns.transpose().rref()
    if len(pivots) != columns.shape[1]:
        raise ValueError("operator basis is dependent")
    orthogonal = []
    for row in reduced.to_list():
        value = row[:]
        for prior in orthogonal:
            coefficient = _inner(prior, value)/_inner(prior, prior)
            value = [x - coefficient*y for x, y in zip(value, prior)]
        orthogonal.append(value)
    return DomainMatrix(orthogonal, reduced.shape, sp.QQ_I).transpose()


def _action_columns(columns, omega):
    rows = columns.to_list()
    vector = exact_columns([omega]).to_list()
    size, count = len(omega), columns.shape[1]
    result = [[sum((rows[i*size+j][k]*vector[j][0] for j in range(size)),
                   sp.QQ_I.zero) for k in range(count)] for i in range(size)]
    return DomainMatrix(result, (size, count), sp.QQ_I)


def small_separation_modulus(basis, omega):
    """Resolve exact kernels and directly evaluable tiny positive actions.

    A full-rank action outside the resolved SVD domain is refused unless exact
    orthogonal action directions give its singular values without cancellation.
    """
    supplied = exact_columns(basis)
    action = _action_columns(supplied, omega)
    if action.rank() < supplied.shape[1]:
        return 0.0
    orthogonal = supplied if supplied.shape[1] == 1 else orthogonal_columns(supplied)
    directions = orthogonal.transpose().to_list()
    images = _action_columns(orthogonal, omega).transpose().to_list()
    if any(_inner(a, b) != 0 for i, a in enumerate(images) for b in images[:i]):
        raise ValueError("positive separation modulus is numerically unresolved")
    squared = min((_inner(a, a)/_inner(b, b)).x
                  for a, b in zip(images, directions))
    result = float(sp.sqrt(sp.QQ.to_sympy(squared)))
    if result == 0:
        raise ValueError("positive separation modulus underflows; result is unresolved")
    return result


def intersect_columns(left, right):
    # Uu=Vv characterizes the intersection. Independent columns in U and V
    # make the map (u,v)->Uu injective on this exact nullspace.
    null = left.hstack(-right).nullspace()
    coefficients = null.extract(range(null.shape[0]), range(left.shape[1])).transpose()
    return left.matmul(coefficients)


def contains_identity(columns, size):
    """Test exact membership, allowing any coefficients in the given basis."""
    identity = exact_columns([np.eye(size)])
    return columns.hstack(identity).rank() == columns.shape[1]


def numerical_basis(columns, size):
    """A resolved HS-orthonormal chart, preconditioned before float conversion."""
    entries = orthogonal_columns(columns).to_Matrix()
    result = []
    for j in range(entries.cols):
        pairs = [x.as_real_imag() for x in entries[:, j]]
        scale = max(abs(x) for pair in pairs for x in pair)
        if not scale:
            raise ValueError("exact operator span returned a dependent basis")
        values = []
        for r, i in pairs:
            real, imag = float(r/scale), float(i/scale)
            if (r != 0 and real == 0) or (i != 0 and imag == 0):
                raise ValueError("operator basis underflow; numerical representation is unresolved")
            values.append(complex(real, imag))
        result.append(np.array(values).reshape(size, size))
    normalized = unit_columns(np.column_stack([b.reshape(-1) for b in result]))
    return [normalized[:, j].reshape(size, size) for j in range(len(result))]

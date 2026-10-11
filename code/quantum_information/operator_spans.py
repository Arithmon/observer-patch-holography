"""Exact supplied operator spans, separate from their numerical charts."""

import math

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


def numerical_chart(columns):
    """Return a guarded numerical chart, or request exact preconditioning."""
    try:
        normalized = unit_columns(columns)
    except ValueError:
        return None
    singular = np.linalg.svd(normalized, compute_uv=False)
    # Accuracy gate, not a rank cutoff. Small source directions remain part
    # of the exact span when this presentation needs preconditioning.
    if singular[-1]*32 < singular[0]:
        return None
    # Householder subtraction can erase small leading coordinates, including
    # part of an occupied row. Reorthogonalized modified Gram-Schmidt avoids
    # that subtraction and controls loss of orthogonality on this gated path.
    q = np.empty_like(normalized)
    q[:, 0] = normalized[:, 0]
    try:
        with np.errstate(under="ignore"):
            for j in range(1, columns.shape[1]):
                value = normalized[:, j].copy()
                for _ in range(2):
                    for k in range(j):
                        value -= q[:, k]*np.vdot(q[:, k], value)
                q[:, j] = unit_columns(value[:, None])[:, 0]
    except ValueError:
        return None
    if np.any(np.any(columns != 0, axis=1) & ~np.any(q != 0, axis=1)):
        return None
    return q


def _inner(left, right):
    return sum((sp.QQ_I(x.x, -x.y)*y for x, y in zip(left, right)), sp.QQ_I.zero)


def _independent_columns(columns):
    reduced, pivots = columns.transpose().rref()
    if len(pivots) != columns.shape[1]:
        raise ValueError("operator basis is dependent")
    return reduced.transpose()


def orthogonal_columns(columns):
    """Hermitian-orthogonal rational directions of an independent source span.

    Canonical row reduction removes supplied units and shears. Gram-Schmidt
    then uses exact QQ(i) products, so conversion to floats happens only after
    cancellation and conditioning in the supplied presentation are removed.
    """
    reduced = _independent_columns(columns).transpose()
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


def _rounded_positive_sqrt(squared):
    """Round sqrt of a positive rational directly to binary64, ties to even.

    An intermediate rounded square root can erase which side of a subnormal
    midpoint the exact value occupies. Instead choose the binary64 spacing,
    take an integer square root, and compare the exact squared midpoint.
    """
    numerator, denominator = int(squared.numerator), int(squared.denominator)
    exponent = (numerator.bit_length() - denominator.bit_length()) // 2
    if ((numerator < denominator << (2*exponent)) if exponent >= 0
            else (numerator << (-2*exponent) < denominator)):
        exponent -= 1
    spacing = max(exponent - 52, -1074)
    if spacing < 0:
        numerator <<= -2*spacing
    else:
        denominator <<= 2*spacing
    lower = math.isqrt(numerator // denominator)
    comparison = 4*numerator - denominator*(2*lower + 1)**2
    if comparison > 0 or (comparison == 0 and lower % 2):
        lower += 1
    if lower == 0:
        raise ValueError("positive separation modulus underflows; result is unresolved")
    return math.ldexp(float(lower), spacing)


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
    if any(_inner(a, b) != sp.QQ_I.zero
           for i, a in enumerate(images) for b in images[:i]):
        raise ValueError("positive separation modulus is numerically unresolved")
    squared = min((_inner(a, a)/_inner(b, b)).x
                  for a, b in zip(images, directions))
    return _rounded_positive_sqrt(squared)


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


def _scaled_exact_columns(columns):
    """Convert exact directions only after removing their arbitrary units."""
    entries = columns.to_Matrix()
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
        result.append(values)
    return np.array(result, dtype=complex).T


def numerical_basis(columns, size):
    """A guarded chart of the exact source span, with exact preconditioning."""
    canonical = _independent_columns(columns)
    try:
        q = numerical_chart(_scaled_exact_columns(canonical))
    except ValueError:
        q = None
    if q is None:
        # Orthogonalize exactly only if the canonical chart needs it: exact
        # Gram-Schmidt can itself introduce unrepresentable tiny coordinates
        # even when a well-conditioned numerical chart is available.
        q = unit_columns(_scaled_exact_columns(orthogonal_columns(canonical)))
    return [q[:, j].reshape(size, size) for j in range(columns.shape[1])]

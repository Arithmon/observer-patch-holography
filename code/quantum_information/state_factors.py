"""Exact source factors for finite-state matrix functions and fidelity.

Binary64 conversion and the shared input tolerances are checked first. The
Hermitian average and a pivoted LDL decomposition then use exact rational
complex arithmetic. Positive pivots are never removed by an eigenvalue
cutoff. Multiprecision refinement below is a numerical accuracy guard, not
an interval enclosure or an exact certificate for irrational answers.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt, nextafter

import mpmath
import numpy as np

from .gibbs import _numeric
from .states import ATOL


_ZERO = (Fraction(0), Fraction(0))
_ONE = (Fraction(1), Fraction(0))


class UnresolvedStatePrecision(ValueError):
    """The requested state calculation is unresolved at available precision."""


def _add(a, b):
    return a[0] + b[0], a[1] + b[1]


def _subtract(a, b):
    return a[0] - b[0], a[1] - b[1]


def _conjugate(a):
    return a[0], -a[1]


def _multiply(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def _scale(a, scale):
    return a[0] * scale, a[1] * scale


def _norm_squared(a):
    return a[0] ** 2 + a[1] ** 2


def _divide(a, b):
    return _scale(_multiply(a, _conjugate(b)), 1 / _norm_squared(b))


def factor_entries(entries):
    """Return exact LDL columns and positive pivots of a Hermitian matrix.

    Columns retain the original coordinate order. A largest remaining
    diagonal pivot avoids a needlessly ill-conditioned chart. A PSD matrix
    with zero diagonal has a zero corresponding row; an exact nonzero
    residual with no positive pivot is therefore indefinite.
    """
    residual = [list(row) for row in entries]
    size = len(residual)
    if any(len(row) != size for row in residual) or any(
            residual[i][j] != _conjugate(residual[j][i])
            for i in range(size) for j in range(size)):
        raise ValueError("state factorization requires a Hermitian matrix")
    columns, pivots = [], []
    for _ in range(size):
        diagonal = [residual[i][i][0] for i in range(size)]
        if any(value < 0 for value in diagonal):
            raise ValueError("state must be positive semidefinite")
        pivot_index = max(range(size), key=diagonal.__getitem__)
        pivot = diagonal[pivot_index]
        if not pivot:
            if any(value != _ZERO for row in residual for value in row):
                raise ValueError("state must be positive semidefinite")
            break
        column = tuple(_scale(residual[i][pivot_index], 1 / pivot)
                       for i in range(size))
        columns.append(column)
        pivots.append(pivot)
        for i in range(size):
            for j in range(i, size):
                updated = _subtract(residual[i][j],
                                    _scale(_multiply(column[i], _conjugate(column[j])), pivot))
                residual[i][j] = updated
                residual[j][i] = _conjugate(updated)
    return tuple(columns), tuple(pivots)


@dataclass(frozen=True)
class SourceState:
    """Original converted entries and their exact Hermitian PSD factors."""

    raw: np.ndarray
    entries: tuple
    columns: tuple
    pivots: tuple
    trace: Fraction

    @property
    def rank(self):
        return len(self.pivots)


def source_state(value):
    """Validate without letting rounded averaging determine source support."""
    raw = _numeric(value, "state")
    if raw.ndim != 2 or not raw.shape[0] or raw.shape[0] != raw.shape[1]:
        raise ValueError("nonempty square numeric density matrix required")
    size = len(raw)
    original = tuple(tuple((Fraction(float(z.real)), Fraction(float(z.imag)))
                           for z in row) for row in raw)
    # Apply the shared tolerances to the original entries, without a
    # throwaway binary64 norm, trace reduction, average, or eigensolve.
    # Squared exact comparisons also preserve the caller's NumPy policy.
    tolerance_squared = Fraction(ATOL)**2
    hermitian_defect = sum((_norm_squared(_subtract(original[i][j],
                              _conjugate(original[j][i])))
                            for i in range(size) for j in range(size)), Fraction(0))
    if hermitian_defect > tolerance_squared:
        raise ValueError("density matrix must be Hermitian")
    raw_trace = _ZERO
    for i in range(size):
        raw_trace = _add(raw_trace, original[i][i])
    if _norm_squared(_subtract(raw_trace, _ONE)) > tolerance_squared:
        raise ValueError("density matrix must have trace one")
    entries = tuple(tuple(_scale(_add(original[i][j], _conjugate(original[j][i])),
                                 Fraction(1, 2))
                          for j in range(size)) for i in range(size))
    columns, pivots = factor_entries(entries)
    raw.setflags(write=False)
    return SourceState(raw, entries, columns, pivots,
                       sum((entries[i][i][0] for i in range(size)), Fraction(0)))


def _mp_real(ctx, value):
    return ctx.mpf(value.numerator) / value.denominator


def mp_matrix(ctx, entries):
    return ctx.matrix([[ctx.mpc(_mp_real(ctx, real), _mp_real(ctx, imag))
                        for real, imag in row] for row in entries])


def fraction_float(value):
    """Correctly round an exact rational, refusing lost nonzero output."""
    try:
        result = float(value)
    except OverflowError as error:
        raise UnresolvedStatePrecision("state output exceeds finite precision") from error
    if not np.isfinite(result) or (value and result == 0):
        raise UnresolvedStatePrecision("state output underflow or insufficient precision")
    return result


def _mp_fraction(value):
    # An mpf is a finite dyadic value. Converting that value to Fraction first
    # avoids an intermediate float conversion double-rounding a subnormal.
    sign, mantissa, exponent, _ = value._mpf_
    numerator = -int(mantissa) if sign else int(mantissa)
    if exponent >= 0:
        return Fraction(numerator << exponent)
    return Fraction(numerator, 1 << -exponent)


def refined_float(ctx, value, uncertainty):
    """Round a refined real result only away from unresolved float midpoints.

    ``uncertainty`` is a numerical refinement guard, not a proven enclosure.
    Exact rational results should use ``fraction_float`` instead.
    """
    if not ctx.isfinite(value) or not ctx.isfinite(uncertainty):
        raise UnresolvedStatePrecision("nonfinite state precision estimate")
    exact, error = _mp_fraction(value), abs(_mp_fraction(uncertainty))
    result = fraction_float(exact)
    for neighbor in (nextafter(result, -np.inf), nextafter(result, np.inf)):
        if np.isfinite(neighbor):
            midpoint = (Fraction(result) + Fraction(float(neighbor))) / 2
            if abs(exact - midpoint) <= error:
                raise UnresolvedStatePrecision("state output rounding midpoint is unresolved")
    return result


def rounded_mp(ctx, value, resolution=0):
    """Shared conversion for a real multiprecision matrix-function output."""
    return refined_float(ctx, value, ctx.mpf(resolution))


def _overlap_matrix(left, right):
    # Complex pairs cannot use Python's tuple-concatenating sum.
    rows = []
    for column_a in left.columns:
        row = []
        for column_b in right.columns:
            value = _ZERO
            for a, b in zip(column_a, column_b):
                value = _add(value, _multiply(_conjugate(a), b))
            row.append(value)
        rows.append(tuple(row))
    return tuple(rows)


def _rank(entries):
    if not entries:
        return 0
    remaining = [list(row) for row in entries]
    rank = 0
    for column in range(len(remaining[0])):
        pivot = next((i for i in range(rank, len(remaining))
                      if remaining[i][column] != _ZERO), None)
        if pivot is None:
            continue
        remaining[rank], remaining[pivot] = remaining[pivot], remaining[rank]
        for i in range(rank + 1, len(remaining)):
            if remaining[i][column] == _ZERO:
                continue
            ratio = _divide(remaining[i][column], remaining[rank][column])
            for j in range(column, len(remaining[0])):
                remaining[i][j] = _subtract(remaining[i][j],
                                           _multiply(ratio, remaining[rank][j]))
        rank += 1
        if rank == len(remaining):
            break
    return rank


def _trace_product(left, right, overlap):
    return sum((p * q * _norm_squared(overlap[i][j])
                for i, p in enumerate(left.pivots)
                for j, q in enumerate(right.pivots)), Fraction(0))


def _rank_two_moment(left, right, overlap):
    """Second singular-value moment from exact weighted 2x2 minors.

    Cauchy--Binet gives e2=sum_{i<k,j<l}|minor(R)|^2 for
    R=sqrt(D) C sqrt(E). At rank two it is the product of the two
    nonzero squared singular values, hence F=Tr(A B)+2*sqrt(e2).
    """
    moment = Fraction(0)
    for i, p in enumerate(left.pivots):
        for k in range(i + 1, left.rank):
            for j, q in enumerate(right.pivots):
                for l in range(j + 1, right.rank):
                    minor = _subtract(_multiply(overlap[i][j], overlap[k][l]),
                                      _multiply(overlap[i][l], overlap[k][j]))
                    moment += p * left.pivots[k] * q * right.pivots[l] * _norm_squared(minor)
    return moment


def _rational_sqrt(value):
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    if numerator ** 2 == value.numerator and denominator ** 2 == value.denominator:
        return Fraction(numerator, denominator)
    return None


def _fidelity_at_precision(left, right, overlap, rank, precision, rank_two=None,
                           diagonal=None):
    ctx = mpmath.mp.clone()
    ctx.dps = precision
    if diagonal is not None:
        return ctx, ctx.fsum(ctx.sqrt(_mp_real(ctx, p)) for p in diagonal)**2
    if rank_two is not None:
        trace, moment = rank_two
        return ctx, _mp_real(ctx, trace) + 2 * ctx.sqrt(_mp_real(ctx, moment))
    product = mp_matrix(ctx, overlap)
    for i, p in enumerate(left.pivots):
        for j, q in enumerate(right.pivots):
            product[i, j] *= ctx.sqrt(_mp_real(ctx, p * q))
    try:
        singular = sorted(ctx.svd_c(product, compute_uv=False), reverse=True)
    except RuntimeError as error:
        if not str(error).startswith("svd: no convergence"):
            raise
        raise UnresolvedStatePrecision("fidelity singular values failed to converge") from error
    scale = max(singular)
    noise = ctx.eps * scale * max(product.rows, product.cols) ** 3 * 1000
    # The exact overlap rank determines how many values are positive; the
    # numerical budget only decides whether that known support is resolved.
    if (any(value <= noise * 10 ** 10 for value in singular[:rank])
            or any(abs(value) > noise for value in singular[rank:])):
        raise UnresolvedStatePrecision("fidelity support is unresolved at this precision")
    value = ctx.fsum(singular[:rank]) ** 2
    if not ctx.isfinite(value) or value <= 0:
        raise UnresolvedStatePrecision("fidelity is unresolved at this precision")
    return ctx, value


def fidelity_factors(left, right):
    """Squared homogeneous fidelity from exact source LDL overlaps."""
    if len(left.entries) != len(right.entries):
        raise ValueError("fidelity requires states on the same space")
    if left.entries == right.entries:
        return fraction_float(left.trace ** 2)
    diagonal = None
    rank_two = None
    size = len(left.entries)
    if all(source.entries[i][j] == _ZERO for source in (left, right)
           for i in range(size) for j in range(size) if i != j):
        diagonal = tuple(left.entries[i][i][0] * right.entries[i][i][0]
                         for i in range(size))
        # F=(sum sqrt(p_i))^2. Even irrational individual roots can have
        # rational cross products, giving an exact arbitrary-rank answer
        # at binary64 midpoints. No SVD or positive-mass threshold is needed.
        exact = sum(diagonal, Fraction(0))
        for i, p in enumerate(diagonal):
            for q in diagonal[i+1:]:
                root = _rational_sqrt(p*q)
                if root is None:
                    exact = None
                    break
                exact += 2*root
            if exact is None:
                break
        if exact is not None:
            return fraction_float(exact)
        overlap, rank = None, None
    else:
        overlap = _overlap_matrix(left, right)
        rank = _rank(overlap)
        if not rank:
            return 0.0
        if rank == 1:
            # A single positive singular value has nuclear norm squared equal
            # to the exact Frobenius norm squared: Tr(A B).
            return fraction_float(_trace_product(left, right, overlap))
        if rank == 2:
            trace = _trace_product(left, right, overlap)
            moment = _rank_two_moment(left, right, overlap)
            root = _rational_sqrt(moment)
            if root is not None:
                return fraction_float(trace + 2 * root)
            rank_two = trace, moment
    failure = "fidelity remains unresolved at available precision"
    for precision in (80, 240, 640, 1280):
        try:
            coarse, value = _fidelity_at_precision(left, right, overlap, rank, precision,
                                                  rank_two, diagonal)
            ctx, check = _fidelity_at_precision(left, right, overlap, rank, precision + 40,
                                              rank_two, diagonal)
            difference = abs(ctx.mpf(value) - check)
            if difference > abs(check) * ctx.mpf("1e-25"):
                continue
            # Agreement alone can hide an erased term at both precisions.
            # Include the coarser arithmetic's scale before selecting which
            # side of a binary64 rounding midpoint contains the result.
            uncertainty = max(8 * difference,
                              ctx.mpf(coarse.eps) * abs(check)
                              * len(left.entries) ** 3 * 1000)
            return refined_float(ctx, check, uncertainty)
        except (UnresolvedStatePrecision, ZeroDivisionError) as error:
            failure = str(error)
    raise UnresolvedStatePrecision(failure)

"""Source-preserving finite matrix functions with explicit precision refusal.

Exact arithmetic determines support and structural zeros. Refinement checks
the numerical function on that support; it is not an interval proof.
"""
from fractions import Fraction
from math import isqrt

import mpmath
import numpy as np

from .state_factors import (
    SourceState, UnresolvedStatePrecision, _ZERO, _ONE, _add, _multiply, _scale, _conjugate,
    _subtract, factor_entries, fraction_float, mp_matrix, rounded_mp,
)


def _product(a, b):
    result = []
    for row in a:
        out = []
        for column in zip(*b):
            value = _ZERO
            for x, y in zip(row, column):
                value = _add(value, _multiply(x, y))
            out.append(value)
        result.append(tuple(out))
    return tuple(result)


def _function_zeros(entries):
    """Zeros in A,...,A^n persist for every spectral function with f(0)=0."""
    n = len(entries)
    power = entries
    zeros = np.ones((n, n, 2), bool)
    for _ in range(n):
        for i in range(n):
            for j in range(n):
                for part in range(2):
                    zeros[i, j, part] &= power[i][j][part] == 0
        power = _product(power, entries)
    return zeros


def _entries(array):
    return tuple(tuple((Fraction(float(z.real)), Fraction(float(z.imag)))
                       for z in row) for row in array)


def _rational_spectrum(entries):
    import sympy as sp
    matrix = sp.Matrix([[sp.Rational(r.numerator, r.denominator) +
                         sp.I * sp.Rational(s.numerator, s.denominator)
                         for r, s in row] for row in entries])
    roots = matrix.charpoly().as_poly().ground_roots()
    if sum(roots.values()) != len(entries):
        return None
    return {Fraction(int(value.p), int(value.q)): count for value, count in roots.items()}


def _rational_function_zeros(entries, inverse):
    """Exact spectral cancellations when the source eigenvalues are rational.

    Polynomial spectral projectors are rational complex matrices. Grouping
    square roots by rational square class certifies cancellations without
    numerical eigenvectors, integer factorization or symbolic tree equality.
    """
    roots = _rational_spectrum(entries)
    if roots is None:
        return None
    n = len(entries)
    identity = tuple(tuple(_ONE if i == j else _ZERO for j in range(n)) for i in range(n))
    classes = {}
    for value in roots:
        if not value:
            continue
        projector = identity
        for other in roots:
            if other == value:
                continue
            factor = tuple(tuple(_scale(_subtract(entries[i][j],
                                (other, Fraction(0)) if i == j else _ZERO), 1 / (value - other))
                                 for j in range(n)) for i in range(n))
            projector = _product(projector, factor)
        representative, ratio = value, Fraction(1)
        for candidate in classes:
            quotient = value / candidate
            num, den = isqrt(quotient.numerator), isqrt(quotient.denominator)
            if num * num == quotient.numerator and den * den == quotient.denominator:
                representative, ratio = candidate, Fraction(num, den)
                break
        coefficient = 1 / (ratio * representative) if inverse else ratio
        previous = classes.get(representative, tuple(tuple(_ZERO for _ in range(n)) for _ in range(n)))
        classes[representative] = tuple(tuple(_add(previous[i][j], _scale(projector[i][j], coefficient))
                                             for j in range(n)) for i in range(n))
    return np.array([[[all(block[i][j][part] == 0 for block in classes.values())
                       for part in range(2)] for j in range(n)] for i in range(n)])


def _exact_root(candidate, source, unresolved):
    # Certify only the invariant blocks containing uncertain components.
    # An unrelated irrational root must not prevent an exact block identity.
    selected = set(unresolved)
    for _ in range(len(source.entries)):
        selected.update(j for i in tuple(selected)
                        for j, value in enumerate(source.entries[i]) if value != _ZERO)
    indices = sorted(selected)
    entries = _entries(candidate[np.ix_(indices, indices)])
    target = tuple(tuple(source.entries[i][j] for j in indices) for i in indices)
    if _product(entries, entries) != target:
        return False
    try:
        factor_entries(entries)
    except ValueError:
        return False
    return True


def _eigenvalues(ctx, matrix, *, vectors=False):
    try:
        return ctx.eighe(matrix, eigvals_only=not vectors)
    except RuntimeError as error:
        if not str(error).startswith("tridiag_eigen: no convergence"):
            raise
        raise UnresolvedStatePrecision("state eigensolver failed to converge") from error


def _root_at_precision(source, inverse, precision):
    ctx = mpmath.mp.clone()
    ctx.dps = precision
    values, vectors = _eigenvalues(ctx, mp_matrix(ctx, source.entries), vectors=True)
    n, rank = len(source.entries), source.rank
    noise = ctx.eps * max(abs(v) for v in values) * n ** 3 * 1000
    split = n - rank
    if (any(abs(values[i]) > noise for i in range(split)) or
            any(values[i] <= noise * 10 ** 10 for i in range(split, n))):
        raise UnresolvedStatePrecision("source square-root support is unresolved")
    weights = ctx.diag([0 if i < split else
                        (1 / ctx.sqrt(values[i]) if inverse else ctx.sqrt(values[i]))
                        for i in range(n)])
    return ctx, vectors * weights * vectors.H


def source_root(source, *, inverse=False):
    """PSD root (or support inverse root) of the exact Hermitian source."""
    n = len(source.entries)
    zeros = _function_zeros(source.entries)
    spectral_zeros_checked = False
    failure = "source square root remains unresolved at available precision"
    for precision in (80, 240, 640, 1280):
        try:
            coarse, first = _root_at_precision(source, inverse, precision)
            ctx, second = _root_at_precision(source, inverse, precision + 40)
            scale = max(abs(z) for z in second)
            floor = ctx.mpf(coarse.eps) * scale * n ** 3 * 1000
            result = np.zeros((n, n), complex)
            unresolved = set()
            uncertain_parts = []
            for i in range(n):
                for j in range(n):
                    components = []
                    for part, component in enumerate((ctx.re, ctx.im)):
                        value = component(second[i, j])
                        delta = abs(value - component(first[i, j]))
                        if zeros[i, j, part]:
                            if abs(value) > floor:
                                raise UnresolvedStatePrecision("structural root zero is unresolved")
                            components.append(0.)
                        elif (abs(value) <= floor * 10 ** 10 or
                              delta > abs(value) * ctx.mpf("1e-25")):
                            # Zero is only a candidate until the exact root
                            # identity below certifies the entire matrix.
                            unresolved.update((i, j))
                            uncertain_parts.append((i, j, part))
                            components.append(0.)
                        else:
                            components.append(rounded_mp(ctx, value, max(floor, 8 * delta)))
                    result[i, j] = complex(*components)
            if not unresolved or (not inverse and _exact_root(result, source, unresolved)):
                return result
            if not spectral_zeros_checked:
                extra = _rational_function_zeros(source.entries, inverse)
                spectral_zeros_checked = True
                if extra is not None:
                    zeros |= extra
            if all(zeros[index] for index in uncertain_parts):
                return result
        except (UnresolvedStatePrecision, ZeroDivisionError) as error:
            failure = str(error)
    raise UnresolvedStatePrecision(failure)


def source_reduction(source, dims, keep):
    """Partial trace before binary64 summation or Hermitian averaging."""
    size = len(source.entries)
    if int(np.prod(dims)) != size:
        raise ValueError("state and tensor dimensions do not match")
    traced = tuple(i for i in range(len(dims)) if i not in keep)
    coordinates = list(np.ndindex(*dims))
    kept_dims = tuple(dims[i] for i in keep)
    out_size = int(np.prod(kept_dims))
    entries = [[_ZERO for _ in range(out_size)] for _ in range(out_size)]
    for i, row in enumerate(coordinates):
        ri = np.ravel_multi_index(tuple(row[k] for k in keep), kept_dims)
        for j, col in enumerate(coordinates):
            if any(row[k] != col[k] for k in traced):
                continue
            cj = np.ravel_multi_index(tuple(col[k] for k in keep), kept_dims)
            entries[ri][cj] = _add(entries[ri][cj], source.entries[i][j])
    entries = tuple(tuple(row) for row in entries)
    columns, pivots = factor_entries(entries)
    # This view is only for the conservative binary64 inverse-domain check.
    # Numerical functions still consume the exact entries retained above.
    raw = np.array([[complex(float(r), float(s)) for r, s in row] for row in entries])
    raw.setflags(write=False)
    return SourceState(raw, entries, columns, pivots, source.trace)


def entries_array(entries):
    """Convert a source matrix, refusing to erase a nonzero component."""
    return np.array([[complex(fraction_float(r), fraction_float(s))
                      for r, s in row] for row in entries])


def channel_entries(kraus, entries):
    """Exact action of the returned binary64 Kraus factors on source entries."""
    size = len(kraus[0])
    result = tuple(tuple(_ZERO for _ in range(size)) for _ in range(size))
    for matrix in kraus:
        factor = _entries(matrix)
        adjoint = tuple(tuple(_conjugate(z) for z in column) for column in zip(*factor))
        term = _product(_product(factor, entries), adjoint)
        result = tuple(tuple(_add(a, b) for a, b in zip(ra, rb)) for ra, rb in zip(result, term))
    return result


def source_distance(left, right):
    """Half trace norm, forming the difference in original-entry arithmetic."""
    n = len(left.entries)
    if n != len(right.entries):
        raise ValueError("trace distance requires states on the same space")
    entries = tuple(tuple(_subtract(a, b) for a, b in zip(ra, rb))
                    for ra, rb in zip(left.entries, right.entries))
    if all(entries[i][j] == _ZERO for i in range(n) for j in range(n) if i != j):
        return fraction_float(sum((abs(entries[i][i][0]) for i in range(n)), Fraction(0)) / 2)
    failure = "trace distance remains unresolved at available precision"
    for precision in (80, 240, 640, 1280):
        answers = []
        try:
            for digits in (precision, precision + 40):
                ctx = mpmath.mp.clone()
                ctx.dps = digits
                values = _eigenvalues(ctx, mp_matrix(ctx, entries))
                answers.append((ctx, ctx.fsum(abs(v) for v in values) / 2))
        except UnresolvedStatePrecision as error:
            failure = str(error)
            continue
        coarse, first = answers[0]
        ctx, second = answers[1]
        delta = abs(second - first)
        if second <= 0 or delta > abs(second) * ctx.mpf("1e-25"):
            continue
        try:
            return rounded_mp(ctx, second, max(8 * delta, ctx.mpf(coarse.eps) * second * n ** 3 * 1000))
        except UnresolvedStatePrecision as error:
            failure = str(error)
        # Resolve exact rational ties independently of numerical eigensolvers.
        if precision == 80:
            roots = _rational_spectrum(entries)
            if roots is not None:
                exact = sum(abs(value) * count for value, count in roots.items()) / 2
                return fraction_float(exact)
    raise UnresolvedStatePrecision(failure)

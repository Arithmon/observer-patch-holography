"""Outward rounding of a supplied numerical Gram factor, without a PSD floor.

The returned matrix is exactly PSD as a matrix of binary64 entries. The
bound concerns rounding the supplied factor's exact Gram matrix, not the
factor's accuracy or statistical uncertainty in a tomographic estimate.
"""

from fractions import Fraction
import math

import numpy as np

from .gibbs import _numeric
from .state_factors import factor_entries


def _parts(a):
    return ([[Fraction(float(z.real)) for z in row] for row in a],
            [[Fraction(float(z.imag)) for z in row] for row in a])


def _is_exact_psd(a):
    """Exact Hermitian Schur elimination on validated binary64 entries.

    A zero diagonal in a PSD Schur complement forces its entire row to
    vanish. Positive pivots preserve PSD iff the next Schur complement is
    PSD. No structural symbolic equality or numerical rank threshold enters.
    """
    if not np.array_equal(a, a.conj().T):
        return False
    try:
        factor_entries(tuple(tuple((Fraction(float(z.real)), Fraction(float(z.imag)))
                                   for z in row) for row in a))
    except ValueError:
        return False
    return True


def _upward(value):
    try:
        rounded = float(value)
    except OverflowError as exc:
        raise ValueError("Gram rounding exceeds finite numerical range") from exc
    if not math.isfinite(rounded):
        raise ValueError("Gram rounding exceeds finite numerical range")
    if Fraction(rounded) < value:
        with np.errstate(over="ignore", under="ignore"):
            rounded = float(np.nextafter(rounded, np.inf))
    if not math.isfinite(rounded):
        raise ValueError("Gram rounding exceeds finite numerical range")
    return rounded


def round_positive_gram(factor):
    """Return M >= B B* and an outward bound on Tr(M-B B*), exactly.

    Off-diagonals are rounded once, with exact rational error accounting.
    Upward diagonals dominate each row's off-diagonal rounding errors.
    Thus M-B B* is Hermitian diagonally dominant with nonnegative diagonal,
    hence PSD. Its trace bounds its trace norm and operator norm. No final
    renormalization is permitted: it would require another rounding proof.
    """
    a = _numeric(factor, "Gram factor")
    if a.ndim != 2 or not all(a.shape):
        raise ValueError("nonempty rectangular Gram factor required")
    real, imag = _parts(a)
    size, width = a.shape
    result = np.zeros((size, size), complex)
    budgets = [Fraction(0) for _ in range(size)]
    diagonal = [sum((r*r+s*s for r, s in zip(rr, ss)), Fraction(0))
                for rr, ss in zip(real, imag)]
    for i in range(size):
        for j in range(i):
            r = sum((real[i][k]*real[j][k]+imag[i][k]*imag[j][k]
                     for k in range(width)), Fraction(0))
            s = sum((imag[i][k]*real[j][k]-real[i][k]*imag[j][k]
                     for k in range(width)), Fraction(0))
            try:
                z = complex(float(r), float(s))
            except OverflowError as exc:
                raise ValueError("Gram rounding exceeds finite numerical range") from exc
            if not (math.isfinite(z.real) and math.isfinite(z.imag)):
                raise ValueError("Gram rounding exceeds finite numerical range")
            result[i, j], result[j, i] = z, z.conjugate()
            error = abs(Fraction(z.real)-r)+abs(Fraction(z.imag)-s)
            budgets[i] += error
            budgets[j] += error
    for i in range(size):
        result[i, i] = _upward(diagonal[i]+budgets[i])
    increase = sum((Fraction(float(result[i, i].real))-diagonal[i]
                    for i in range(size)), Fraction(0))
    return result, _upward(increase)

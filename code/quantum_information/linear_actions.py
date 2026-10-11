"""Evaluate finite linear operations without discarding cancelled observables.

Exact fallbacks concern supplied binary64 entries. Algebra projections use
their original span, whereas a supplied superoperator keeps its own entries.
"""

from fractions import Fraction

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def _exact_matrix(value):
    def rational(x):
        return sp.QQ(*float(x).as_integer_ratio())
    return DomainMatrix([
        [sp.QQ_I(rational(x.real), rational(x.imag)) for x in row]
        for row in value], value.shape, sp.QQ_I).to_sparse()


def _adjoint(matrix):
    return DomainMatrix([
        [sp.QQ_I(x.x, -x.y) for x in row]
        for row in matrix.transpose().to_list()], matrix.shape[::-1], sp.QQ_I).to_sparse()


def _rounded_component(value):
    try:
        result = float(Fraction(int(value.numerator), int(value.denominator)))
    except OverflowError as error:
        raise ValueError("linear result exceeds finite numerical range") from error
    if value and result == 0:
        raise ValueError("nonzero linear result underflows; precision is insufficient")
    return result


def _rounded_vector(vector):
    return np.array([complex(_rounded_component(row[0].x),
                             _rounded_component(row[0].y))
                     for row in vector.to_list()])


def _sqrt_ratio(squared):
    """Reuse exact span rounding, with the distance API's range errors."""
    if not squared:
        return 0.
    from .operator_spans import _rounded_positive_sqrt
    try:
        return _rounded_positive_sqrt(squared)
    except ValueError as error:
        raise ValueError("positive distance underflows; precision is insufficient") from error
    except OverflowError as error:
        raise ValueError("distance exceeds finite numerical range") from error


class SpanProjector:
    """Exact factored HS projection; no ambient dense projector is needed."""

    def __init__(self, columns):
        self.columns = columns
        adjoint = _adjoint(columns)
        self.coordinates = adjoint.matmul(columns).inv().matmul(adjoint)

    def evaluate(self, value, *, distance=False, within=None):
        vector = _exact_matrix(value.reshape(-1, 1))
        projected = self.columns.matmul(self.coordinates.matmul(vector))
        if not distance and within is None:
            return _rounded_vector(projected).reshape(value.shape)
        # Subtract and sum squares before rounding. The rounded projection
        # can erase a residual, and individual residual coordinates can be
        # unrepresentable even when their scalar norm is representable.
        residual = vector - projected
        squared = sum((x.x*x.x + x.y*x.y for row in residual.to_list()
                       for x in row), sp.QQ.zero)
        if within is not None:
            # A tolerance decision needs no representable scalar distance.
            threshold = sp.QQ(*float(within).as_integer_ratio())
            return squared <= threshold*threshold
        return _sqrt_ratio(squared)


def apply_linear(matrix, vector):
    """Use a modestly conditioned dot product, else evaluate original entries.

    The guard is componentwise: a large real part cannot hide a lost imaginary
    response. Possible subnormal products take the exact path before any
    reduction, including when several individually underflowing terms add up.
    """
    magnitudes = [np.concatenate((np.abs(a.real).ravel(), np.abs(a.imag).ravel()))
                  for a in (matrix, vector)]
    nonzero = [a[a != 0] for a in magnitudes]
    possible_underflow = (all(a.size for a in nonzero)
                          and sum(int(np.frexp(a.min())[1]) for a in nonzero) < -1020)
    try:
        if possible_underflow:
            raise FloatingPointError
        with np.errstate(over="raise", under="raise", invalid="raise", divide="raise"):
            result = matrix @ vector
            real_sum = np.abs(matrix.real) @ np.abs(vector.real) + np.abs(matrix.imag) @ np.abs(vector.imag)
            imag_sum = np.abs(matrix.real) @ np.abs(vector.imag) + np.abs(matrix.imag) @ np.abs(vector.real)
            if (not np.all(np.isfinite(result))
                    or np.any(np.abs(result.real) < real_sum/32)
                    or np.any(np.abs(result.imag) < imag_sum/32)):
                raise FloatingPointError
        return result
    except FloatingPointError:
        exact = _exact_matrix(matrix).matmul(_exact_matrix(vector.reshape(-1, 1)))
        return _rounded_vector(exact)

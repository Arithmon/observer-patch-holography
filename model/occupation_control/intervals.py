"""Small outward interval kernel; no BLAS operations enter a certificate.

Basic binary64 arithmetic and square root use IEEE-754 round-to-nearest with
gradual underflow. Each elementary result is widened by one adjacent float.
All reductions use explicitly outward pairwise sums. Nonfinite arithmetic
fails closed. Input floats denote exact binary rationals, not decimal values.
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np


def down(x):
    return np.nextafter(x, -np.inf)


def up(x):
    return np.nextafter(x, np.inf)


def exact(x):
    x = np.asarray(x, dtype=float)
    if not np.isfinite(x).all():
        raise ValueError('finite interval input required')
    return x, x


def rational(x):
    x = Fraction(x)
    value = float(x)
    return (value if Fraction(value) <= x else float(down(value)),
            value if Fraction(value) >= x else float(up(value)))


def add(a, b):
    return down(a[0]+b[0]), up(a[1]+b[1])


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    candidates = np.array([a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]])
    return down(np.min(candidates, axis=0)), up(np.max(candidates, axis=0))


def scale(a, x):
    return mul(a, exact(x))


def square_root(a):
    if np.any(np.asarray(a[0]) < 0):
        raise ValueError('nonnegative interval required')
    return np.maximum(0., down(np.sqrt(a[0]))), up(np.sqrt(a[1]))


def sum_up(values):
    values = np.asarray(values, dtype=float).reshape(-1)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ArithmeticError('finite nonnegative summands required')
    if not len(values):
        return 0.
    while len(values) > 1:
        if len(values) % 2:
            values = np.r_[values, 0.]
        values = up(values[::2]+values[1::2])
        if not np.isfinite(values).all():
            raise ArithmeticError('nonfinite outward sum')
    return float(values[0])


def frobenius_upper(real, imag):
    ar = np.maximum(abs(real[0]), abs(real[1]))
    ai = np.maximum(abs(imag[0]), abs(imag[1]))
    total = sum_up(np.r_[up(ar*ar).ravel(), up(ai*ai).ravel()])
    value = float(up(math.sqrt(total)))
    if not math.isfinite(value):
        raise ArithmeticError('nonfinite norm bound')
    return value


def trace_upper(real, imag):
    """Nuclear norm <= sqrt(rank bound) times Frobenius norm."""
    rank = min(np.shape(real[0]))
    return float(up(up(math.sqrt(rank))*frobenius_upper(real, imag)))


def complex_exact(matrix):
    return exact(matrix.real), exact(matrix.imag)


def complex_sub(a, b):
    return sub(a[0], b[0]), sub(a[1], b[1])


def complex_scale(a, factor):
    return mul(a[0], factor), mul(a[1], factor)

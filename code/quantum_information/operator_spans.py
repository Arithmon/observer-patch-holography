"""Exact supplied operator spans, separate from their numerical QR charts."""

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def exact_columns(basis):
    entries = sp.Matrix.hstack(*[
        sp.Matrix([sp.Rational(float(x.real))+sp.I*sp.Rational(float(x.imag))
                   for x in b.flat]) for b in basis])
    return DomainMatrix.from_Matrix(entries).convert_to(sp.QQ_I)


def intersect_columns(left, right):
    # Uu=Vv characterizes the intersection. Independent columns in U and V
    # make the map (u,v)->Uu injective on this exact nullspace.
    null = left.hstack(-right).nullspace()
    coefficients = null.extract(range(null.shape[0]), range(left.shape[1])).transpose()
    return left.matmul(coefficients)


def numerical_basis(columns, size):
    entries = columns.to_Matrix()
    result = []
    for j in range(entries.cols):
        pairs = [x.as_real_imag() for x in entries[:, j]]
        scale = max(abs(x) for pair in pairs for x in pair)
        if not scale:
            raise ValueError("exact intersection returned a dependent basis")
        values = []
        for r, i in pairs:
            real, imag = float(r/scale), float(i/scale)
            if (r != 0 and real == 0) or (i != 0 and imag == 0):
                raise ValueError("intersection basis underflow; numerical representation is unresolved")
            values.append(complex(real, imag))
        result.append(np.array(values).reshape(size, size))
    return result

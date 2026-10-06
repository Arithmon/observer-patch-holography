"""Accumulate supplied finite Hamiltonians before rounding products or sums.

The caller validates real coefficients and Hermitian operators. Binary64
components are rational; only assembly uses exact arithmetic. Diagonalization
and the returned Gibbs states remain numerical. The usual tolerated input
Hermiticity roundoff is removed by taking the Hermitian part exactly.
"""

from fractions import Fraction
from dataclasses import dataclass
import math

import numpy as np


def _sum_components(operators, coefficients, component):
    total = [Fraction(0)]*operators[0].size
    for operator, coefficient in zip(operators, coefficients):
        if coefficient == 0:
            continue
        coefficient = Fraction(float(coefficient))
        values = getattr(operator, component).reshape(-1)
        # Local operators are often sparse and have repeated entries. Cache
        # products, not rounded products, and avoid work on structural zeros.
        nonzero = np.flatnonzero(values)
        products = {value: coefficient*Fraction(float(value))
                    for value in np.unique(values[nonzero])}
        for index in nonzero:
            total[index] += products[values[index]]
    return np.array(total, dtype=object).reshape(operators[0].shape)


def _rounded(value):
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError("Hamiltonian exceeds finite numerical range") from exc
    if not math.isfinite(result):
        raise ValueError("Hamiltonian exceeds finite numerical range")
    if value != 0 and result == 0:
        raise ValueError("Hamiltonian underflow; precision is insufficient")
    return result


def _upper_float(value):
    result = float(value)
    if Fraction(result) < value:
        result = np.nextafter(result, math.inf)
    return result


@dataclass(frozen=True)
class HamiltonianAssembly:
    matrix: np.ndarray
    scalar: float
    operator_rounding_bound: float
    scalar_rounding_bound: float
    ideal_gibbs_trace_distance_bound: float
    ideal_log_partition_error_bound: float


def _assemble_hamiltonian(operators, coefficients, *, centered):
    """Return the assembled matrix, scalar and outward-rounded error bounds.

    With centered=True, remove Tr(H)/d before conversion. Otherwise return
    the full matrix and scalar zero. Round only the final components, reject
    erased nonzero components and overflow, and never drop a summand.
    Bounds compare ideal Gibbs functions of the exact sum and returned
    matrix/scalar; they do NOT bound a subsequent numerical eigensolve.
    This private numerical primitive expects validated nonempty matched data.
    """
    real = _sum_components(operators, coefficients, 'real')
    imag = _sum_components(operators, coefficients, 'imag')
    real = (real+real.T)/2
    imag = (imag-imag.T)/2
    size = len(real)
    scalar = sum(real.diagonal(), Fraction(0))/size if centered else Fraction(0)
    for i in range(size):
        real[i, i] -= scalar
    cache = {}
    for value in {*real.flat, *imag.flat, scalar}:
        rounded = _rounded(value)
        cache[value] = (rounded, abs(Fraction(rounded)-value))
    rounded_real = np.array([cache[x][0] for x in real.flat]).reshape(real.shape)
    rounded_imag = np.array([cache[x][0] for x in imag.flat]).reshape(imag.shape)
    row_errors = [sum((cache[r][1]+cache[i][1] for r, i in zip(row_r, row_i)), Fraction(0))
                  for row_r, row_i in zip(real, imag)]
    # For Hermitian error E, ||E||_op <= max_i sum_j |E_ij|. Bounding a
    # complex modulus by |Re|+|Im| keeps this accumulation entirely rational.
    operator_error = max(row_errors)
    scalar_error = cache[scalar][1]
    matrix = rounded_real+1j*rounded_imag
    matrix.setflags(write=False)
    return HamiltonianAssembly(
        matrix, cache[scalar][0],
        _upper_float(operator_error), _upper_float(scalar_error),
        _upper_float(min(Fraction(2), operator_error)),
        _upper_float(operator_error+scalar_error))

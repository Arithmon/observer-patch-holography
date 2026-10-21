"""Outward spectral and thermal traces from exact factorial-metric forms.

LAPACK supplies candidates only. Dyadic source rounding, Gram products and
eigen-equation residuals are certified with integer arithmetic. Four parity
blocks keep the exact Gram calculation small; checked int64 limb products
are reconstructed as Python integers before any wide accumulation.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np

from ..scalar_observations.bounds import sqrt_upper


def _rational(value, *, positive=False):
    if (type(value) not in (int, F) or (positive and value <= 0)
            or (type(value) is F and (type(value.numerator) is not int
                or type(value.denominator) is not int or value.denominator <= 0
                or math.gcd(value.numerator, value.denominator) != 1))):
        raise ValueError('an original canonical rational is required')
    return F(value)


def _rounded(value, scale, upper=False):
    numerator = value.numerator*scale
    integer = (-((-numerator)//value.denominator) if upper
               else numerator//value.denominator)
    return F(integer, scale)


def exponential_interval(value, bits=96):
    """Enclose exp(value) by rational Taylor series and directed squaring.

    The supported exponent range prevents an untrusted input from requesting
    an unbounded-size integer. No floating exponential or libm assumption is
    used. Rounding after each squaring keeps denominator sizes bounded.
    """
    value = _rational(value)
    if type(bits) is not int or not 32 <= bits <= 256 or abs(value) > 16384:
        raise ValueError('32..256 arithmetic bits and exponent magnitude <= 16384 required')
    if not value:
        return F(1), F(1)
    reduced, squares = abs(value), 0
    while reduced > F(1, 2):
        reduced /= 2
        squares += 1
    term = total = F(1)
    degree = 48
    for order in range(1, degree+1):
        term *= reduced/order
        total += term
    following = term*reduced/(degree+1)
    remainder = following/(1-reduced/(degree+2))
    lower, upper = total, total+remainder
    if value < 0:
        lower, upper = 1/upper, 1/lower
    scale = 1 << bits
    lower, upper = _rounded(lower, scale), _rounded(upper, scale, True)
    for _ in range(squares):
        lower = _rounded(lower*lower, scale)
        upper = _rounded(upper*upper, scale, True)
    return lower, upper


def exact_gram(integer_basis, limb_bits=21):
    """Exact U.T U for int64 candidates, with independently checked limbs."""
    if (type(integer_basis) is not np.ndarray or integer_basis.ndim != 2
            or integer_basis.dtype != np.dtype('int64')
            or min(integer_basis.shape) < 1):
        raise ValueError('a nonempty original int64 matrix is required')
    if type(limb_bits) is not int or not 1 <= limb_bits <= 30:
        raise ValueError('integer limb width from 1 through 30 required')
    n = integer_basis.shape[0]
    limit = np.iinfo(np.int64).max
    if n*((1 << limb_bits)-1)**2 > limit:
        raise ValueError('int64 limb dot products can overflow')
    if np.any(integer_basis == np.iinfo(np.int64).min):
        raise ValueError('the minimum int64 value has no signed absolute representation')
    magnitude, sign = np.abs(integer_basis), np.sign(integer_basis)
    maximum = int(magnitude.max())
    limbs = max(1, (maximum.bit_length()+limb_bits-1)//limb_bits)
    chunks = [((magnitude >> (index*limb_bits)) & ((1 << limb_bits)-1))*sign
              for index in range(limbs)]
    result = np.zeros((integer_basis.shape[1],)*2, dtype=object)
    for i, left in enumerate(chunks):
        for j in range(i, len(chunks)):
            right = chunks[j]
            product = left.T@right  # The absolute sum bound above prevents wrap.
            wide = product.astype(object)*(1 << (limb_bits*(i+j)))
            result += wide
            if i != j:
                result += wide.T
    return result


def _form_data(form):
    d, entries, metric = form.occupation, form.entries, form.metric
    if type(d) is not int or not 2 <= d <= 64 or form.dimension != d*d:
        raise ValueError('a square two-mode occupation from 2 through 64 required')
    n = d*d
    if (type(entries) is not tuple or type(metric) is not tuple or len(metric) != n
            or any(type(weight) is not int or weight <= 0 for weight in metric)):
        raise ValueError('immutable exact sparse entries and a positive integer metric required')
    matrix = {}
    for entry in entries:
        if (type(entry) is not tuple or len(entry) != 3
                or any(type(index) is not int or not 0 <= index < n for index in entry[:2])):
            raise ValueError('original in-range sparse row and column addresses required')
        row, column, value = entry
        value = _rational(value)
        if (row, column) in matrix:
            raise ValueError('duplicate sparse entries are ambiguous')
        matrix[row, column] = value
    for (row, column), value in matrix.items():
        if value*metric[row] != matrix.get((column, row), F())*metric[column]:
            raise ValueError('the exact source form is not self-adjoint in its metric')
        if value and ((row % d) % 2, (row//d) % 2) != ((column % d) % 2, (column//d) % 2):
            raise ValueError('the declared form does not preserve both occupation parities')
    return d, matrix, metric


def _dyadic_block(addresses, entries, metric, scale):
    lookup = {address: index for index, address in enumerate(addresses)}
    rows = [dict() for _ in addresses]
    nonexact = 0
    for row in addresses:
        for column, j in lookup.items():
            if row > column:
                continue
            value = entries.get((row, column), F())
            if not value:
                continue
            square = value*value*F(metric[row], metric[column])
            scaled = square*scale*scale
            integer = math.isqrt(scaled.numerator//scaled.denominator)
            exact = integer*integer*scaled.denominator == scaled.numerator
            if not exact and 4*scaled.numerator > (2*integer+1)**2*scaled.denominator:
                integer += 1
            if value < 0:
                integer = -integer
            i = lookup[row]
            rows[i][j] = rows[j][i] = integer
            if not exact:
                nonexact += 1 if row == column else 2
    source_radius = sqrt_upper(F(nonexact, 4*scale*scale))
    return rows, source_radius


def _candidate(values, vectors, size, scale, norm_bound):
    for value, shape in ((values, (size,)), (vectors, (size, size))):
        if (type(value) is not np.ndarray or value.shape != shape
                or value.dtype.kind not in 'iuf' or not np.isfinite(value).all()):
            raise ValueError('complete finite real candidate eigenvalues and eigenvectors required')
    values, vectors = np.asarray(values, dtype=float), np.asarray(vectors, dtype=float)
    if (not np.isfinite(values).all() or not np.isfinite(vectors).all()
            or np.max(np.abs(vectors)) > 2
            or any(abs(F(float(value))) > 2*norm_bound+1 for value in values)):
        raise ValueError('candidate magnitude exceeds the finite source domain')
    basis = np.rint(vectors*scale).astype(np.int64)
    diagonal = tuple(round(F(float(value))*scale) for value in values)
    return basis, diagonal


def _certify_block(rows, source_radius, candidate, scale, return_candidates=False):
    size = len(rows)
    source_norm = max((sum(abs(value) for value in row.values()) for row in rows), default=0)
    source_norm = F(source_norm, scale)
    if candidate is None:
        approximate = np.zeros((size, size))
        for i, row in enumerate(rows):
            for j, value in row.items():
                approximate[i, j] = float(F(value, scale))
        candidate = np.linalg.eigh(approximate)
    if not isinstance(candidate, (tuple, list)) or len(candidate) != 2:
        raise ValueError('candidate must contain eigenvalues and eigenvectors')
    basis, diagonal = _candidate(*candidate, size, scale, source_norm)
    gram = exact_gram(basis)
    squared = 0
    for i in range(size):
        for j in range(size):
            entry = int(gram[i, j])-(scale*scale if i == j else 0)
            squared += entry*entry
    epsilon = sqrt_upper(F(squared, scale**4))
    if epsilon >= 1:
        raise ArithmeticError('candidate does not certify an invertible near-orthonormal basis')
    residual_squared = 0
    for i, row in enumerate(rows):
        for j, eigenvalue in enumerate(diagonal):
            residual = sum(coefficient*int(basis[k, j]) for k, coefficient in row.items())
            residual -= int(basis[i, j])*eigenvalue
            residual_squared += residual*residual
    residual = sqrt_upper(F(residual_squared, scale**4))
    maximum = F(max(map(abs, diagonal), default=0), scale)
    # Polar V of U has ||U-V||_F <= ||U.T U-I||_F=epsilon<1.
    # Hence ||A0-V Lambda V.T|| <= residual+(||A0||+||Lambda||)*epsilon.
    radius = source_radius+residual+(source_norm+maximum)*epsilon
    values = tuple(sorted(F(value, scale) for value in diagonal))
    result = {'dimension': size, 'candidate_values': values,
            'source_rounding_frobenius_upper': source_radius,
            'gram_frobenius_upper': epsilon,
            'eigen_equation_frobenius_upper': residual,
            'dyadic_source_operator_norm_upper': source_norm,
            'eigenvalue_error_upper': radius}
    if return_candidates:
        result.update(dyadic_basis=tuple(tuple(int(value) for value in row) for row in basis),
                      dyadic_eigenvalues=diagonal, dyadic_scale=scale)
    return result


def certify_spectrum(form, beta=F(1, 4), *, bits=40, candidate=None,
                     return_candidates=False):
    """Certify all eigenvalues and Tr exp(-beta H) of the exact finite form.

    An optional candidate maps each parity pair to (values, vectors), with
    entries in that block's increasing original-address order. Every value
    returned for an error or interval endpoint is an exact Fraction. With
    return_candidates=True, each block additionally retains its immutable
    integer dyadic basis, unsorted paired eigenvalues, scale and addresses.
    """
    beta = _rational(beta, positive=True)
    if type(return_candidates) is not bool:
        raise ValueError('an original Boolean candidate-export flag required')
    if type(bits) is not int or not 16 <= bits <= 52:
        raise ValueError('an original dyadic precision from 16 through 52 required')
    d, entries, metric = _form_data(form)
    groups = {parity: tuple(i for i in range(d*d) if (i % d % 2, i//d % 2) == parity)
              for parity in ((0, 0), (0, 1), (1, 0), (1, 1))}
    if candidate is not None and (type(candidate) is not dict or candidate.keys() != groups.keys()):
        raise ValueError('one candidate for each of the four parity blocks required')
    scale = 1 << bits
    blocks, intervals = [], []
    lower = upper = F()
    for parity, addresses in groups.items():
        rows, source_radius = _dyadic_block(addresses, entries, metric, scale)
        block = _certify_block(rows, source_radius,
                               None if candidate is None else candidate[parity], scale,
                               return_candidates)
        block['parity'] = parity
        if return_candidates:
            block['addresses'] = addresses
        blocks.append(block)
        radius = block['eigenvalue_error_upper']
        intervals.extend((value-radius, value+radius) for value in block['candidate_values'])
        center_lower = center_upper = F()
        for value in block['candidate_values']:
            lo, hi = exponential_interval(-beta*value)
            center_lower += lo
            center_upper += hi
        shrink = exponential_interval(-beta*radius)[0]
        growth = exponential_interval(beta*radius)[1]
        lower += shrink*center_lower
        upper += growth*center_upper
    ordered_intervals = tuple(zip(sorted(lo for lo, _ in intervals),
                                  sorted(hi for _, hi in intervals), strict=True))
    return {'occupation': d, 'dimension': d*d, 'inverse_temperature': beta,
            'dyadic_bits': bits, 'blocks': tuple(blocks),
            'eigenvalue_intervals': ordered_intervals,
            'partition_function_interval': (lower, upper),
            'eigensolver_accuracy_assumed': False,
            'arithmetic': 'exact integers, checked int64 limb dots, rational Taylor bounds'}

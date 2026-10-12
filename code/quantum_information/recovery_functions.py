"""Petz outputs evaluated from original exact marginals before float export.

Exact classical conditional products and sufficient polynomial zero
certificates complement private multiprecision evaluation. Refinement and
conditioning guards are numerical safeguards, not interval bounds.
"""
from fractions import Fraction
from math import prod

import numpy as np

from .state_factors import (
    UnresolvedStatePrecision, _ZERO, _ONE, _add, _conjugate, _multiply,
    _mp_real, mp_matrix, rounded_mp, source_state,
)
from .state_functions import _root_at_precision, entries_array, source_reduction


def _diagonal(entries):
    return all(value == _ZERO for i, row in enumerate(entries)
               for j, value in enumerate(row) if i != j)


def _identity(size):
    return tuple(tuple(_ONE if i == j else _ZERO for j in range(size))
                 for i in range(size))


def _product(left, right):
    result = []
    for row in left:
        active = [(k, value) for k, value in enumerate(row) if value != _ZERO]
        result.append(tuple(
            _sum(_multiply(value, right[k][j]) for k, value in active if right[k][j] != _ZERO)
            for j in range(len(right[0]))))
    return tuple(result)


def _sum(values):
    total = _ZERO
    for value in values:
        total = _add(total, value)
    return total


def _tensor(left, right):
    return tuple(tuple(_multiply(a, b) for a in row_left for b in row_right)
                 for row_left in left for row_right in right)


def _scale(entries, coefficient):
    return tuple(tuple((r*coefficient, s*coefficient) for r, s in row)
                 for row in entries)


def _powers(matrix):
    value = matrix
    for _ in range(len(matrix)):
        yield value
        value = _product(value, matrix)


def _component_zeros(ab, bc, middle, dims):
    """Sufficient exact zeros for the full Petz polynomial function span.

    Both Hermitian spectral functions have real polynomial coefficients
    and f(0)=0. For V_ij=BC^i (B^j tensor I), paired unordered terms
    V_ij D V_kl*+V_kl D V_ij* have real coefficients. A component absent
    from every such term is zero in the actual source Petz output.
    """
    a, b, c = dims
    if len(bc.entries) * len(middle.entries) > 64:
        raise UnresolvedStatePrecision("Petz zero certificate exceeds the bounded polynomial span")
    ident_a, ident_c = _identity(a), _identity(c)
    powers_bc, powers_b = list(_powers(bc.entries)), list(_powers(middle.entries))
    terms = []
    seen = set()
    for x in powers_bc:
        for y in powers_b:
            value = _tensor(ident_a, _product(x, _tensor(y, ident_c)))
            if value not in seen:
                seen.add(value)
                terms.append(value)
    supplied = _tensor(ab.entries, ident_c)
    acted = [_product(value, supplied) for value in terms]
    size = a * b * c
    zeros = np.ones((size, size, 2), dtype=bool)
    for i in range(size):
        for j in range(size):
            for t, left in enumerate(acted):
                for u in range(t, len(terms)):
                    value = _sum(_multiply(x, _conjugate(y))
                                 for x, y in zip(left[i], terms[u][j])
                                 if x != _ZERO and y != _ZERO)
                    if u != t:
                        value = _add(value, _sum(
                            _multiply(x, _conjugate(y))
                            for x, y in zip(acted[u][i], terms[t][j])
                            if x != _ZERO and y != _ZERO))
                    for part in range(2):
                        if value[part]:
                            zeros[i, j, part] = False
    return zeros


def _mp_tensor(ctx, left, right):
    return ctx.matrix([[left[i // right.rows, j // right.cols]
                        * right[i % right.rows, j % right.cols]
                        for j in range(left.cols * right.cols)]
                       for i in range(left.rows * right.rows)])


def _at_precision(ab, bc, middle, dims, precision):
    ctx, root = _root_at_precision(bc, False, precision)
    _, inverse = _root_at_precision(middle, True, precision)
    a, _, c = dims
    conditional = root * _mp_tensor(ctx, inverse, ctx.eye(c))
    lifted = _mp_tensor(ctx, ctx.eye(a), conditional)
    supplied = _mp_tensor(ctx, mp_matrix(ctx, ab.entries), ctx.eye(c))
    return ctx, lifted * supplied * lifted.H


def _conditioning(source):
    # With pivoted unit-pivot LDL columns, det(L*L)>=1 by Cauchy--Binet.
    # Thus the least positive eigenvalue is at least product(D)/Tr(A)^(r-1).
    # This deliberately conservative scale also covers singular supports.
    return source.trace ** source.rank / prod(source.pivots)


def _certified_markov_candidate(result, ab, bc, dims):
    """Prove a whole proposed output using its exact state and marginals.

    A Markov state is its own Petz recovery. Exact equality with the supplied
    AB and BC therefore proves it is the required marginal-determined output.
    The order-16 limit bounds computational work, never physical support.
    """
    if len(result) > 16:
        return False
    try:
        candidate = source_state(result)
    except ValueError:
        # A provisional numerical candidate has no PSD or trace guarantee.
        return False
    if (source_reduction(candidate, dims, [0, 1]).entries != ab.entries
            or source_reduction(candidate, dims, [1, 2]).entries != bc.entries):
        return False
    from .information import is_markov_exact
    return is_markov_exact(result, dims, [0], [1], [2])


def source_petz_output(ab, bc, dims):
    """Evaluate the ideal supported Petz action, then export its components."""
    a, b, c = dims
    middle = source_reduction(bc, (b, c), [0])
    size = a * b * c
    # A product reference BC gives X -> X tensor C/Tr(C) on supported B
    # inputs, including the AB marginal of the same positive source.
    # A product AB input gives A tensor BC/Tr(BC). These identities use
    # the marginals only, so also recover non-Markov ABC inputs correctly.
    end_c = source_reduction(bc, (b, c), [1])
    if _scale(bc.entries, bc.trace) == _tensor(middle.entries, end_c.entries):
        return entries_array(_scale(_tensor(ab.entries, end_c.entries), 1/bc.trace))
    end_a = source_reduction(ab, (a, b), [0])
    if _scale(ab.entries, ab.trace) == _tensor(end_a.entries, middle.entries):
        return entries_array(_scale(_tensor(end_a.entries, bc.entries), 1/ab.trace))
    if _diagonal(ab.entries) and _diagonal(bc.entries):
        values = []
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    mass = middle.entries[j][j][0]
                    values.append(ab.entries[i*b+j][i*b+j][0]
                                  * bc.entries[j*c+k][j*c+k][0] / mass if mass else Fraction(0))
        return entries_array(tuple(tuple((values[i], Fraction(0)) if i == j else _ZERO
                                         for j in range(size)) for i in range(size)))
    conditioning = max(Fraction(1), _conditioning(bc), _conditioning(middle))
    zeros = None
    failure = "Petz output remains unresolved at available source precision"
    for precision in (80, 240, 640, 1280):
        try:
            coarse, first = _at_precision(ab, bc, middle, dims, precision)
            ctx, second = _at_precision(ab, bc, middle, dims, precision + 40)
            # Include source support conditioning, intermediate multiplication
            # growth and the coarser arithmetic even when both runs agree.
            floor = (ctx.mpf(coarse.eps) * max(1, max(abs(z) for z in second))
                     * _mp_real(ctx, conditioning) ** 2 * size ** 6 * 10000)
            result = np.zeros((size, size), complex)
            uncertain = []
            for i in range(size):
                for j in range(size):
                    parts = []
                    for part, component in enumerate((ctx.re, ctx.im)):
                        value = component(second[i, j])
                        delta = abs(value - component(first[i, j]))
                        if zeros is not None and zeros[i, j, part]:
                            if abs(value) > floor:
                                raise UnresolvedStatePrecision("Petz structural zero is numerically unresolved")
                            parts.append(0.)
                        elif (abs(value) <= floor * 10 ** 10
                              or delta > abs(value) * ctx.mpf("1e-25")):
                            uncertain.append((i, j, part))
                            parts.append(0.)
                        else:
                            parts.append(rounded_mp(ctx, value, max(floor, 8 * delta)))
                    result[i, j] = complex(*parts)
            if uncertain and zeros is None:
                zeros = _component_zeros(ab, bc, middle, dims)
            if all(zeros[index] for index in uncertain):
                return result
            # These provisional zeros are accepted only if an independent
            # exact certificate proves the entire output, including its
            # marginals. This also recovers a non-Markov input when its
            # marginals have this exact Markov extension.
            if _certified_markov_candidate(result, ab, bc, dims):
                return result
        except (UnresolvedStatePrecision, ZeroDivisionError) as error:
            failure = str(error)
    raise UnresolvedStatePrecision(failure)

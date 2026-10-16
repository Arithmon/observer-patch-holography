"""Exact encoded native targets and reference-uniform reported-map bounds.

An ideal sparse map stores ``(row, column): amplitude**2`` as Fractions;
the amplitudes are nonnegative square roots, never binary approximations.
The target is derived from the integer endpoint law and the four-configuration
logical code. The original reported ndarray is compared on all 512 output
rows, including padding. This comparison certifies supplied numbers; it does
not establish that a source program produced them or certify physical gates.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math

import numpy as np

from ..common_source import exact_output as exact
from ..common_source import quantum_native as source


SHAPE = (512, 2)


def _contract(law, mode):
    law, mode = source._law(law), source._mode(mode)
    if mode == 'full_carrier':
        raise ValueError('an encoded protected or passive logical qubit is required')
    return law, mode


def _first_maps(law, mode):
    maps = {}
    if mode == 'protected_qubit':
        for key, weight in source._code_weights(law, (0, 1)).items():
            maps[key] = {(source.circuit.encode(state), column): weight/2
                         for column, supports in enumerate(source._logical_supports())
                         for state in supports}
    else:
        counts = defaultdict(lambda: defaultdict(int))
        for coin in (0, 1):
            for column, supports in enumerate(source._logical_supports()):
                for state in supports:
                    end, current, contrast = source._endpoint(state, (0, 1), coin)
                    counts[coin, source._raw(law, current, contrast)][
                        source.circuit.encode(end), column] += 1
        maps = {key: {entry: F(count*count, 4) for entry, count in row.items()}
                for key, row in counts.items()}
    return maps


@lru_cache(maxsize=4)
def _tree(law, mode):
    """Private cache accepts only the already validated canonical strings."""
    maps = {}
    if mode == 'protected_qubit':
        for (c1, raw1), first in source._code_weights(law, (0, 1)).items():
            for selected in (0, 1):
                seam = (2, 3) if selected == 0 else (0, 1)
                for (c2, raw2), last in source._code_weights(law, seam).items():
                    maps[c1, raw1, selected, c2, raw2] = {
                        (source.circuit.encode(state), column): first*last/2
                        for column, supports in enumerate(source._logical_supports())
                        for state in supports}
    else:
        counts = defaultdict(lambda: defaultdict(int))
        for c1, c2, selected in product((0, 1), repeat=3):
            seam = (2, 3) if selected == 0 else (0, 1)
            for column, supports in enumerate(source._logical_supports()):
                for state in supports:
                    middle, j1, d1 = source._endpoint(state, (0, 1), c1)
                    end, j2, d2 = source._endpoint(middle, seam, c2)
                    key = c1, source._raw(law, j1, d1), selected, c2, source._raw(law, j2, d2)
                    counts[key][source.circuit.encode(end), column] += 1
        # Each path has amplitude (1/sqrt(2))**3. Add coincident paths
        # before squaring, retaining their constructive interference.
        maps = {key: {entry: F(count*count, 8) for entry, count in row.items()}
                for key, row in counts.items()}
    return tuple((key, tuple(sorted(row.items()))) for key, row in sorted(maps.items()))


def ideal_tree(law, mode):
    """Fresh exact nonzero targets keyed by both native records and selector.

    Missing keys are zero maps. ``ideal_map`` additionally validates the full
    retained history menu and handles the explicit early-abort records.
    """
    law, mode = _contract(law, mode)
    return {key: dict(row) for key, row in _tree(law, mode)}


@lru_cache(maxsize=2)
def _menu(law):
    return frozenset(source._menu(law))


def ideal_map(law, mode, record):
    """Return the exact ideal map for one complete native/work record.

    All four smooth bits must be present on continued records, even though
    only the third bit selects the final seam. Invalid raw first records
    terminate with ``(coin, raw, (), None, None)`` and have the zero target.
    """
    law, mode = _contract(law, mode)
    if (type(record) is not tuple or len(record) != 5
            or type(record[2]) is not tuple
            or any(type(value) is not int for value in (*record[:2], *record[2]))
            or any(value is not None and type(value) is not int for value in record[3:])
            or record not in _menu(law)):
        raise ValueError('a complete canonical native/work record, including aborts, is required')
    c1, raw1, bits, c2, raw2 = record
    if not bits:
        return {}
    key = c1, raw1, bits[2], c2, raw2
    return dict(next((row for candidate, row in _tree(law, mode) if candidate == key), ()))


def _sparse(value):
    if type(value) is not dict:
        raise ValueError('a sparse dict of exact squared ideal amplitudes is required')
    result = {}
    for key, weight in value.items():
        if (type(key) is not tuple or len(key) != 2
                or any(type(index) is not int or not 0 <= index < bound
                       for index, bound in zip(key, SHAPE, strict=True))
                or type(weight) is not F or weight < 0):
            raise ValueError('valid indices and nonnegative Fraction squared amplitudes are required')
        if weight:
            result[key] = weight
    return result


def compare_map(array, ideal):
    """Rational upper bound on half the diamond distance of single CP maps.

    For D=N-M and arbitrary operator X, including any reference system,
    N X N* - M X M* = D X N* + M X D*. Schatten submultiplicativity gives
    ||difference||_1 <= ||D||op (||N||op+||M||op) ||X||_1. Thus eps >= ||D||F
    and m >= ||M||F imply half-diamond <= eps*(2*m+eps)/2. Neither map
    must be trace preserving, and the reported map may be complex.

    Exact binary components and directed rational square-root intervals
    enclose ||N-M||F. No eigensolver, rounded ideal map, probability cap or
    normalization is used. This is an upper bound, not a distance estimator;
    even a harmless global phase can make it loose.
    """
    # Subclasses can override copy or array metadata and replace the supplied
    # entries before exact conversion. Reject them before any array dispatch.
    if type(array) is not np.ndarray:
        raise ValueError('an original base ndarray is required')
    array, ideal = exact._array(array, SHAPE), _sparse(ideal)
    entries = tuple(tuple(exact._component(value) for value in row) for row in array)
    reported = sum((real*real+imag*imag for row in entries for real, imag in row), F(0))
    target = sum(ideal.values(), F(0))
    terms = defaultdict(F)
    for (row, column), weight in ideal.items():
        terms[weight] += entries[row][column][0]
    cross_lower, cross_upper = exact._radical_sum(terms)
    squared_lower, squared_upper = exact._nonnegative(
        reported+target-2*cross_upper, reported+target-2*cross_lower)
    epsilon_lower = exact._sqrt_interval(squared_lower)[0]
    epsilon = exact._sqrt_interval(squared_upper)[1]
    target_norm = exact._sqrt_interval(target)[1]
    codes = frozenset(source.CODES)
    padding = sum((real*real+imag*imag for index, row in enumerate(entries)
                   if index not in codes for real, imag in row), F(0))
    return {
        'distance_upper': epsilon*(2*target_norm+epsilon)/2,
        'frobenius_error_lower': epsilon_lower,
        'frobenius_error_upper': epsilon,
        'frobenius_error_squared_lower': squared_lower,
        'frobenius_error_squared_upper': squared_upper,
        'ideal_frobenius_squared': target,
        'reported_frobenius_squared': reported,
        'padding_frobenius_squared': padding,
    }


def _effect(matrix):
    """Compute M*M exactly for these sparse native maps."""
    rows = defaultdict(dict)
    for (row, column), weight in matrix.items():
        rows[row][column] = weight
    result = {(i, j): F(0) for i, j in product(range(2), repeat=2)}
    for row in rows.values():
        for i, left in row.items():
            for j, right in row.items():
                value = left*right
                numerator, denominator = math.isqrt(value.numerator), math.isqrt(value.denominator)
                if numerator*numerator != value.numerator or denominator*denominator != value.denominator:
                    raise ArithmeticError('native Gram entry is outside the exact rational contract')
                result[i, j] += F(numerator, denominator)
    return result


def completeness(law, mode):
    """Verify adaptive final instruments on each first-record image exactly.

    For every first outcome a and selected bit s, sum_b M[a,s,b]*M[a,s,b]
    equals F[a]*F[a] = q[a] I. Therefore an independently prepared work
    instrument selecting s can be inserted without losing completeness.
    The first-record current weights are exactly (-1,0,+1)=(1/4,1/2,1/4).
    Zero maps supply every absent or abort raw outcome without postselection.
    """
    law, mode = _contract(law, mode)
    tree, first = ideal_tree(law, mode), _first_maps(law, mode)
    weights = source._code_weights(law, (0, 1))
    identity = {(i, j): F(i == j) for i, j in product(range(2), repeat=2)}
    current_weights = dict.fromkeys((-1, 0, 1), F(0))
    conditional = {}
    for prefix, matrix in first.items():
        effect = _effect(matrix)
        if effect != {entry: weights[prefix]*value for entry, value in identity.items()}:
            raise ArithmeticError('first native record does not have its exact scalar code effect')
        current = source._record_current(law, *prefix)
        if current is None:
            raise ArithmeticError('a nonzero ideal first native map was marked abort')
        current_weights[current] += weights[prefix]
        for selected in (0, 1):
            total = dict.fromkeys(identity, F(0))
            for key, row in tree.items():
                if key[:3] == (*prefix, selected):
                    for entry, value in _effect(row).items():
                        total[entry] += value
            if total != effect:
                raise ArithmeticError('adaptive final native instrument is not complete on its first image')
            conditional[(*prefix, selected)] = total
    if set(first) != set(weights):
        raise ArithmeticError('first native outcome support changed')
    for selected in (0, 1):
        total = dict.fromkeys(identity, F(0))
        for key, row in tree.items():
            if key[:2] not in first:
                raise ArithmeticError('final native map has no first outcome')
            if key[2] == selected:
                for entry, value in _effect(row).items():
                    total[entry] += value
        if total != identity:
            raise ArithmeticError('complete selected native instrument is not trace preserving')
    if current_weights != {-1: F(1, 4), 0: F(1, 2), 1: F(1, 4)}:
        raise ArithmeticError('native-current code weights changed')
    return {'current_weights': current_weights, 'first_record_weights': dict(weights),
            'conditional_final_effects': conditional,
            'selected_complete_effect': dict(identity),
            'all_native_output_rows': 512,
            'arbitrary_finite_reference': True}

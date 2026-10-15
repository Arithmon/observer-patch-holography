"""Exact fixed-seam completion attachment for the native binary C/E channels.

The channel forgets every intermediate coin and repair record. It is not an
implementation of a retained-record history or an adaptive seam schedule.
Sparse exact operators use flattened system-first tensor indices and Gaussian
rational entries: a rational, or a pair (real rational, imaginary rational).
The returned entries are always pairs of fractions. NumPy Kraus adapters are
separate floating diagnostics; the matrix-unit audit uses exact arithmetic.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from fractions import Fraction as F
import math

import numpy as np

from ..observer import native


SEAM = (0, 1)
STATES = native.states(binary=True)
DIMENSION = len(STATES)
INDEX = {state: i for i, state in enumerate(STATES)}
EQUAL = tuple(state[0] == state[1] for state in STATES)
ZERO = (F(0), F(0))


def _contract(seam, retain_records, adaptive_seams):
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(value) is not int for value in seam) or seam != SEAM):
        raise ValueError('this attachment requires the fixed seam (0, 1)')
    if type(retain_records) is not bool or retain_records:
        raise ValueError('compression requires every repair record to be forgotten')
    if adaptive_seams is not None:
        raise ValueError('adaptive or multi-seam histories cannot use this compression')


def _iterations(value):
    if type(value) is not int or value < 1:
        raise ValueError('a positive integer number of fixed-seam events is required')
    return value


def _channel(value):
    if type(value) is not str or value not in ('C', 'E'):
        raise ValueError('native channel C or E required')
    return value


def _rational(value):
    if type(value) not in (int, F):
        raise ValueError('exact integer or Fraction entries required')
    return F(value)


def _gaussian(value):
    if type(value) is tuple and len(value) == 2:
        return tuple(_rational(part) for part in value)
    return (_rational(value), F(0))


def _operator(value, reference_dimension):
    if type(reference_dimension) is not int or reference_dimension < 1:
        raise ValueError('a positive integer reference dimension is required')
    if not isinstance(value, Mapping):
        raise ValueError('a sparse mapping of matrix units to exact entries is required')
    limit = DIMENSION * reference_dimension
    result = {}
    for key, coefficient in value.items():
        if (type(key) is not tuple or len(key) != 2
                or any(type(i) is not int or not 0 <= i < limit for i in key)):
            raise ValueError('matrix-unit indices lie outside the declared tensor carrier')
        entry = _gaussian(coefficient)
        if entry != ZERO:
            result[key] = entry
    return result


def _native_unit(channel, i, j, *, states=STATES, index=INDEX):
    """Trace out actual equal-record branches, with exact coin weight 1/2."""
    output = defaultdict(F)
    left, right = states[i], states[j]
    for coin in (0, 1):
        a, ar = native.native(left, SEAM, coin)
        b, br = native.native(right, SEAM, coin)
        same_record = (ar == br if channel == 'C' else
                       left[0] - left[1] == right[0] - right[1])
        if same_record:
            output[index[a], index[b]] += F(1, 2)
    return dict(output)


def _completion_unit(i, j):
    """Independent block partial trace: M46 + (I2 tensor M10)."""
    if EQUAL[i] and EQUAL[j]:
        return {(i, j): F(1)}
    if EQUAL[i] or EQUAL[j] or STATES[i][0] != STATES[j][0]:
        return {}
    output = {}
    for orientation in (0, 1):
        a, b = list(STATES[i]), list(STATES[j])
        a[:2] = b[:2] = (orientation, 1-orientation)
        output[INDEX[tuple(a)], INDEX[tuple(b)]] = F(1, 2)
    return output


def _cross_unit(i, j):
    return {(i, j): F(1)} if EQUAL[i] != EQUAL[j] else {}


def _real_apply(matrix, unit_action):
    result = defaultdict(F)
    for (i, j), value in matrix.items():
        for key, coefficient in unit_action(i, j).items():
            result[key] += value * coefficient
    return {key: value for key, value in result.items() if value}


def _plus(left, right, scale=F(1)):
    result = defaultdict(F, left)
    for key, value in right.items():
        result[key] += scale * value
    return {key: value for key, value in result.items() if value}


def apply_compressed(operator, iterations, *, reference_dimension=1, seam=SEAM,
                     retain_records=False, adaptive_seams=None):
    """Apply exact C**n, retaining the complete system/reference operator.

    Inputs need not be normalized, positive or Hermitian: this is the linear
    extension on the entire declared algebra. No record register is returned.
    State entries can be complex Gaussian rationals, with arbitrary finite
    reference dimension and system-first index ``system*r + reference``.
    """
    _contract(seam, retain_records, adaptive_seams)
    n = _iterations(iterations)
    matrix = _operator(operator, reference_dimension)
    scale = F(1, 2**n)
    output = {}
    for (i, j), value in matrix.items():
        a, r = divmod(i, reference_dimension)
        b, s = divmod(j, reference_dimension)
        image = _plus(_completion_unit(a, b), _cross_unit(a, b), scale)
        for (x, y), coefficient in image.items():
            key = (x*reference_dimension+r, y*reference_dimension+s)
            previous = output.get(key, ZERO)
            output[key] = (previous[0]+coefficient*value[0],
                           previous[1]+coefficient*value[1])
    return {key: value for key, value in output.items() if value != ZERO}


def binary_kraus(channel='C', *, seam=SEAM, retain_records=False, adaptive_seams=None):
    """Fresh floating Kraus matrices assembled from native transitions.

    Branch labels are discarded by contract. Completeness and the ideal
    algebraic meaning follow from the exact native matrix-unit attachment;
    the floating matrices are not a roundoff or physical service certificate.
    """
    _contract(seam, retain_records, adaptive_seams)
    channel = _channel(channel)
    groups = {}
    for i, state in enumerate(STATES):
        for coin in (0, 1):
            output, record = native.native(state, SEAM, coin)
            key = record if channel == 'C' else (coin, state[0]-state[1])
            groups.setdefault(key, np.zeros((DIMENSION, DIMENSION), dtype=complex))[
                INDEX[output], i] = 1/math.sqrt(2)
    return tuple(groups[key] for key in sorted(groups))


def compressed_kraus(iterations, *, seam=SEAM, retain_records=False, adaptive_seams=None):
    """C**n as the positive mixture (1-2**(1-n))*E + 2**(1-n)*C.

    This executable adapter compresses a forgotten-record channel. It neither
    generates the eliminated event records nor reduces physical elapsed time.
    """
    _contract(seam, retain_records, adaptive_seams)
    n = _iterations(iterations)
    exact = F(1, 2**(n-1))
    weight = float(exact)
    if not math.isfinite(weight) or weight <= 0:
        raise ValueError('positive mixture weight is unresolved in binary64')
    result = [math.sqrt(weight)*branch for branch in binary_kraus('C')]
    if n > 1:
        result.extend(math.sqrt(float(1-exact))*branch for branch in binary_kraus('E'))
    return tuple(result)


def half_diamond_distance(iterations):
    """Exact half-diamond distance from C**n to the completion E."""
    return F(1, 2**(_iterations(iterations)+1))


def witness_indices():
    """Equal-block preparation a, unequal u and its opposite orientation v."""
    return tuple(INDEX[tuple(int(k in pair) for k in range(12))]
                 for pair in ((2, 3), (0, 2), (1, 2)))


def _witness():
    a, u, _ = witness_indices()
    state = {(a, a): F(1, 2), (u, u): F(1, 2),
             (a, u): F(1, 2), (u, a): F(1, 2)}
    completed = _real_apply(state, _completion_unit)
    plus = lambda matrix: sum(matrix.get(key, F(0))/2 for key in state)
    reference_plus = plus(completed)
    if reference_plus != F(3, 8):
        raise ArithmeticError('completion witness has the wrong complete binary law')
    rows = []
    for n in (1, 2, 3, 8):
        actual = apply_compressed(state, n)
        if any(value[1] for value in actual.values()):
            raise ArithmeticError('real witness acquired imaginary entries')
        probability = plus({key: value[0] for key, value in actual.items()})
        expected = reference_plus + half_diamond_distance(n)
        if probability != expected or not 0 <= probability <= 1:
            raise ArithmeticError('complete binary measurement fails to attain the bound')
        rows.append({'events': n, 'C_power_binary_law': [str(probability), str(1-probability)],
                     'E_binary_law': ['3/8', '5/8'],
                     'total_variation': str(half_diamond_distance(n))})
    return {'state': '(|a>+|u>)/sqrt(2)', 'basis_indices_a_u_v': list(witness_indices()),
            'plus_effect': '|a+u><a+u|/2', 'minus_effect': 'I - plus_effect',
            'postselection': False, 'runs': rows}


def _full_carrier_counterexample():
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    image = defaultdict(F)
    for i in range(len(states)):
        for key, coefficient in _native_unit('E', i, i, states=states, index=index).items():
            image[key] += coefficient
    changes = []
    for i, state in enumerate(states):
        diagonal = image[i, i]
        if diagonal != 1:
            changes.append({'loads': list(state), 'E_of_identity_diagonal': str(diagonal)})
    if (len(states) != 78 or len(changes) != 3
            or sorted(F(row['E_of_identity_diagonal']) for row in changes) != [0, 0, 3]
            or sum(image.values()) != len(states)):
        raise ArithmeticError('full-carrier nonunitality counterexample failed')
    return {'dimension': 78, 'unital': False, 'changed_diagonal_entries': changes,
            'trace_preserved': True, 'binary_tracial_completion_transfers': False}


def audit():
    """Verify exact source attachment and return a small JSON-ready summary."""
    if DIMENSION != 66 or sum(EQUAL) != 46:
        raise ArithmeticError('the native binary-shell decomposition changed')
    traces = {'E': F(0), 'Q': F(0)}
    units = 0
    for i in range(DIMENSION):
        for j in range(DIMENSION):
            c, e = (_native_unit(law, i, j) for law in ('C', 'E'))
            q = _cross_unit(i, j)
            if e != _completion_unit(i, j) or c != _plus(e, q, F(1, 2)):
                raise ArithmeticError('native C/E differs from the exact completion decomposition')
            if (_real_apply(e, _completion_unit) != e
                    or _real_apply(q, _cross_unit) != q
                    or _real_apply(q, _completion_unit)
                    or _real_apply(e, _cross_unit)):
                raise ArithmeticError('completion and cross-block projections fail their identities')
            c_action = lambda a, b: _native_unit('C', a, b)
            if _real_apply(c, c_action) != _plus(e, q, F(1, 4)):
                raise ArithmeticError('the second actual native event disagrees with compression')
            for key, coefficient in e.items():
                if _completion_unit(*key).get((i, j), F(0)) != coefficient:
                    raise ArithmeticError('completion is not Hilbert-Schmidt self-adjoint')
            traces['E'] += e.get((i, j), F(0))
            traces['Q'] += q.get((i, j), F(0))
            units += 1
    if traces != {'E': F(2216), 'Q': F(1840)}:
        raise ArithmeticError('exact spectral multiplicities changed')
    return {
        'schema': 'oph.primitive_observer.fixed_seam_completion.v1',
        'all_pass': True,
        'native_matrix_units_checked': units,
        'fixed_seam': list(SEAM),
        'binary_dimension': DIMENSION,
        'retained_algebra': 'M46 direct-sum (I2 tensor M10)',
        'reference': 'I66/66',
        'identities': ['E^2=E', 'Q^2=Q', 'EQ=QE=0', 'C=E+Q/2',
                       'C^n=E+2^(-n)Q for integer n>=1'],
        'C_eigenvalue_multiplicities': {'1': 2216, '1/2': 1840, '0': 300},
        'exact_half_diamond_distance': '2^(-n-1) for integer n>=1',
        'arbitrary_reference_bound': True,
        'saturating_complete_binary_read': _witness(),
        'full_carrier_counterexample': _full_carrier_counterexample(),
        'scope': {
            'parameter': 'number of fixed-seam native events; not physical duration',
            'all_repair_records_forgotten': True,
            'retained_record_history_preserved': False,
            'adaptive_seams_supported': False,
            'retained_algebra_and_reference_supplied': True,
            'finite_step_law_selected_by_completion': False,
            'physical_prediction_or_empirical_validation': False,
            'native_pulse_compilation': False,
        },
    }

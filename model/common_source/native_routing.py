"""Execute native port routing and recorded repair through source words.

The nine-qubit encoding is unchanged. Swapping physical port labels 1 and 5
conjugates the existing (0,1) repair into the (0,5) repair. One actually
executed nine-bit controlled flip is reused by fixed wire remapping and
source X conjugation. All 512 matrix columns and padding amplitudes remain.
This is a finite encoded construction under the compiler's declared services,
not a new physical load/export or memory implementation.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import math

import numpy as np

from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_discrimination import experiment as recorded
from ..observer.source_scalar import compiler
from .experiment import COUNT_KEYS


QUBITS = circuit.SOURCE_QUBITS
SIZE = 1 << QUBITS
SEAMS = ((0, 1), (0, 5))


def _immutable(value):
    value = np.ascontiguousarray(value, dtype=np.complex128)
    return np.frombuffer(value.tobytes(), dtype=np.complex128).reshape(value.shape)


def _seam(seam):
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(port) is not int for port in seam) or seam not in SEAMS):
        raise ValueError('the implemented native seam is (0,1) or (0,5)')
    return seam


def _candidate(candidate):
    if type(candidate) is not str or candidate not in ('C', 'E'):
        raise ValueError('the recorded native candidate must be C or E')
    return candidate


def _original_scalars(value):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError('native amplitudes cannot be Boolean, missing or masked')
    if isinstance(value, (tuple, list)):
        for item in value:
            _original_scalars(item)
    elif isinstance(value, np.ndarray) and value.dtype.kind not in 'iufc':
        raise ValueError('native amplitudes require an original numeric dtype')


def _state(value):
    _original_scalars(value)
    raw = np.asarray(value)
    if (raw.ndim not in (1, 2) or raw.shape[0] != SIZE
            or (raw.ndim == 2 and raw.shape[1] < 1)
            or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('finite 512-component native amplitudes or reference columns required')
    result = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(result).all():
        raise ValueError('native amplitudes must remain finite in the execution representation')
    return result


def routing_pairs():
    """Disjoint code transpositions for the port-label swap; padding is fixed."""
    seen, pairs = set(), []
    for state in native.states():
        routed = list(state)
        routed[1], routed[5] = routed[5], routed[1]
        routed = tuple(routed)
        if state != routed and state not in seen:
            pairs.append((circuit.encode(state), circuit.encode(routed)))
            seen.update((state, routed))
    return tuple(pairs)


def routing_edges(*, inverse=False):
    """Chronological Gray-neighbor swaps, each (target_bit, control_pattern)."""
    if type(inverse) is not bool:
        raise ValueError('inverse must be an explicit Boolean')
    edges = []
    for left, right in routing_pairs():
        path = [left]
        for bit in range(QUBITS):
            if (left ^ right) >> bit & 1:
                path.append(path[-1] ^ (1 << bit))
        adjacent = tuple(zip(path[:-1], path[1:], strict=True))
        for a, b in (*adjacent, *adjacent[-2::-1]):
            target = (a ^ b).bit_length()-1
            edges.append((target, a & ~(1 << target)))
    return tuple(reversed(edges)) if inverse else tuple(edges)


def _canonical_word(inverse=False):
    word = circuit.controlled_flip_word(QUBITS, QUBITS-1, 0)
    return compiler.inverse(word) if inverse else word


@lru_cache(maxsize=2)
def _executed_macro(inverse):
    # Execute the inverse's negative-angle reversed word independently. An
    # adjoint of a floating matrix is not substituted for this source replay.
    raw = compiler.execute(_canonical_word(inverse), QUBITS)
    matrix, guard = compiler.stabilize_executed_unitary(raw)
    return _immutable(matrix), guard


def _chart(matrix, target):
    # This is the tensor chart of compiler.remap(word, wires), not a physical
    # port permutation supplied for free. Every native request uses these wires.
    wires = [bit for bit in range(QUBITS) if bit != target]+[target]
    indices = np.arange(SIZE)
    source_indices = sum(((indices >> wire) & 1) << local
                         for local, wire in enumerate(wires))
    return matrix[np.ix_(source_indices, source_indices)]


def _conjugation_word(pattern):
    return [compiler.pulse('X', math.pi/2, [bit])
            for bit in range(QUBITS) if pattern >> bit & 1]


@lru_cache(maxsize=2)
def _executed_route(inverse):
    macro, macro_guard = _executed_macro(inverse)
    charts = {target: _chart(macro, target) for target in range(QUBITS)}
    result = np.eye(SIZE, dtype=complex)
    for target, pattern in routing_edges(inverse=inverse):
        before = _conjugation_word(pattern)
        result = compiler.execute(before, QUBITS, result)
        result = charts[target]@result
        result = compiler.execute(compiler.inverse(before), QUBITS, result)
    # Both guards act on the WHOLE 512-dimensional matrix, never its code
    # compression. Numerical padding leakage is neither clipped nor erased.
    result, route_guard = compiler.stabilize_executed_unitary(result)
    return _immutable(result), {'macro': macro_guard, 'route': route_guard}


def routing_unitary(*, inverse=False):
    """Immutable full source-executed matrix, including every padding column."""
    if type(inverse) is not bool:
        raise ValueError('inverse must be an explicit Boolean')
    return _executed_route(inverse)[0]


def routing_word_summary():
    """Compact exact schedule and source-word resource counts for each route.

    Matrix caches reduce host work only. Every repeated native pulse, load
    and export remains charged in each physical source-word invocation.
    """
    word = _canonical_word()
    macro = compiler.resources(word)
    edges = routing_edges()
    conjugations = 2*sum(pattern.bit_count() for _, pattern in edges)
    one_x = compiler.resources([compiler.pulse('X', math.pi/2, [0])])
    counts = {key: len(edges)*macro[key]+conjugations*one_x[key] for key in COUNT_KEYS}
    serialized = json.dumps(word, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return {'physical_port_permutation': [1, 5], 'native_carrier_dimension': 78,
            'executed_register_dimension': SIZE, 'native_code_qubits': QUBITS,
            'padding_columns_retained': SIZE-78, 'disjoint_transpositions': len(routing_pairs()),
            'controlled_flip_invocations_per_direction': len(edges),
            'canonical_controlled_flip_native_pulses': macro['native_pulses'],
            'canonical_controlled_flip_word_sha256': hashlib.sha256(serialized.encode()).hexdigest(),
            'conjugating_source_x_pulses_per_direction': conjugations,
            'resources_per_direction': counts,
            'inverse_source_word_executed_independently': True,
            'tensor_chart_is_fixed_wire_addressing': True,
            'physical_services_or_roundoff_certified': False}


def routing_execution_diagnostics():
    """Guard diagnostics from the complete forward and inverse executions."""
    forward, first = _executed_route(False)
    inverse, second = _executed_route(True)
    return {'forward': deepcopy(first), 'inverse': deepcopy(second),
            'inverse_composition_operator_error': float(np.linalg.norm(inverse@forward-np.eye(SIZE), 2)),
            'full_512_columns_executed': True, 'padding_projection_applied': False,
            'errors_are_numerical_diagnostics_not_certificates': True}


def routed_repair_branches(candidate, seam, state):
    """Complete source-recorded repair on the two implemented physical seams.

    Input/output arrays have 512 rows and may have arbitrary reference columns;
    unnormalized inputs remain unnormalized. Every raw record, including an
    abort code, is retained. ``abort`` means an invalid RECORD code: quantum
    amplitudes in the native-code complement remain in the full output array.
    The native interpretation concerns the 78-dimensional input code only.
    """
    candidate, seam, value = _candidate(candidate), _seam(seam), _state(state)
    if seam == (0, 5):
        value = routing_unitary()@value
    branches = recorded.repair_branches(candidate, value)
    for branch in branches:
        amplitude = branch['amplitude']
        if seam == (0, 5):
            amplitude = routing_unitary(inverse=True)@amplitude
        branch['amplitude'] = _immutable(amplitude)
        branch['seam'] = seam
    return branches


def routed_repair_resources(candidate, seam):
    """One complete attempted repair path, including its recorded source coin.

    Counts cover source words, not physical initialization/storage/readout costs.
    The input source is already supplied; no new native input preparation is
    charged. All raw record codes are retained even where their norm is zero.
    """
    candidate, seam = _candidate(candidate), _seam(seam)
    coins = [compiler.resources(circuit.repair_word(candidate, coin)) for coin in (0, 1)]
    if any(coins[0][key] != coins[1][key] for key in COUNT_KEYS):
        raise ArithmeticError('native repair branches no longer have common path costs')
    coin = compiler.resources([compiler.pulse('Y', -math.pi/4, [0])])
    route = routing_word_summary()['resources_per_direction'] if seam == (0, 5) else dict.fromkeys(COUNT_KEYS, 0)
    counts = {key: coins[0][key]+coin[key]+2*route[key] for key in COUNT_KEYS}
    return {**counts, 'route_directions_per_path': 2 if seam == (0, 5) else 0,
            'retained_classical_bits': 1+circuit.RECORD_BITS[candidate],
            'blank_register_initializations': 2, 'classical_readout_registers': 2,
            'native_input_preparations': 0,
            'peak_live_code_qubits': QUBITS+circuit.RECORD_BITS[candidate],
            'six_level_processor_additional': True,
            'physical_initialization_storage_readout_costs_included': False}

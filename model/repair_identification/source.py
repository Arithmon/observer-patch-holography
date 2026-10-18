"""Emitted source probes for the complete native record-law coordinates.

The nine-bit total-two encoding is unchanged. Preparations and detectors
are routed by actual three-bit source words, not target permutation gates.
All raw outcomes and all 512 output amplitudes survive. Coherent proper-code
loading, blank preparation and classical readout remain declared services.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import math
from typing import NamedTuple

import numpy as np

from . import core
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from ..observer.scalar_response.experiment import canonical


CANDIDATES = ('C', 'E', 'C_phase_quarter', 'E_collision_half')
PROBE_NAMES = tuple(probe.name for probe in core.PROBES)
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS,
    'model/repair_identification/source.py',
    'model/observer/source_discrimination/circuit.py',
    'model/observer/source_discrimination/experiment.py',
    'model/observer/source_scalar/compiler.py',
    'model/observer/scalar_response/experiment.py',
    compiler.RESPONSE, compiler.CONTROL_RECEIPT))))
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'zero_angle_requests', 'identity_phase_pulses',
              'proper_code_load_events', 'proper_code_export_events',
              'proper_code_service_events')


class Setting(NamedTuple):
    name: str
    parameter: str
    phase: int
    current: int
    input_pair: tuple
    output_pair: tuple


class Branch(NamedTuple):
    # A repair record is (coin, raw_code, auxiliary_tuple). An observed
    # branch appends pointer +1/-1; aborts append None and retain their state.
    record: tuple
    current: int | None
    contrast: int | None
    abort: bool
    amplitude: np.ndarray
    resources: tuple


class ProbeExecution(NamedTuple):
    candidate: str
    probe: str
    reference_dimension: int
    branches: tuple
    word_hashes: tuple


class HistoryBranch(NamedTuple):
    record: tuple
    currents: tuple
    abort: bool
    amplitude: np.ndarray
    resources: tuple


class HistoryExecution(NamedTuple):
    candidate: str
    probe: str
    reference_dimension: int
    feedback: bool
    branches: tuple
    word_hashes: tuple


def _candidate(value):
    if type(value) is not str or value not in CANDIDATES:
        raise ValueError('an original source candidate name is required')
    return value


def settings(name):
    if type(name) is not str or name not in PROBE_NAMES:
        raise ValueError('an original named quadrature probe is required')
    probe = next(probe for probe in core.PROBES if probe.name == name)
    pairs = tuple(tuple(circuit.SOURCE_CODES[i] for i in pair)
                  for pair in (probe.input_pair, probe.output_pair))
    if any(not 0 <= i < 6 for pair in pairs for i in pair):
        raise ArithmeticError('the retained probes left the declared spectator-zero chart')
    return Setting(name, probe.parameter, probe.phase, probe.current, *pairs)


def _flags(reference=False, feedback=False):
    if type(reference) is not bool or type(feedback) is not bool:
        raise ValueError('original Boolean reference and feedback choices required')


def _state(value):
    if (type(value) is not np.ndarray or value.ndim not in (1, 2)
            or value.shape[0] != 512 or (value.ndim == 2 and not 1 <= value.shape[1] <= 512)
            or value.dtype.kind not in 'iufc'
            or (value.dtype.kind == 'f' and value.dtype.itemsize > 8)
            or (value.dtype.kind == 'c' and value.dtype.itemsize > 16)
            or not np.isfinite(value).all()):
        raise ValueError('an original finite 512-amplitude array or reference-column array is required')
    return np.array(value, dtype=np.complex128, copy=True)


def _immutable(value):
    value = np.ascontiguousarray(value, dtype=np.complex128)
    return np.frombuffer(value.tobytes(), dtype=np.complex128).reshape(value.shape)


def _sum_resources(*rows):
    result = {}
    for row in rows:
        for key, count in row:
            result[key] = result.get(key, 0)+count
    return tuple(sorted(result.items()))


def _resources(word=(), **services):
    counts = compiler.resources(word)
    return tuple(sorted((*((key, counts[key]) for key in COUNT_KEYS), *services.items())))


def _permutation_word(pair):
    images = list(range(8))
    for source, target in enumerate(pair):
        other = images.index(target)
        images[source], images[other] = images[other], images[source]
    seen, word = set(), []
    for start in range(8):
        if start in seen:
            continue
        cycle, current = [], start
        while current not in seen:
            seen.add(current)
            cycle.append(current)
            current = images[current]
        for other in cycle[1:]:
            word.extend(circuit.basis_swap_word(start, other, 3))
    return word


def _preparation_word(probe, reference=False):
    setting = settings(probe)
    word = circuit.preparation_word('phase_pi_over_2' if setting.phase else 'phase_0')
    if reference:
        word += compiler.cnot(0, 9)
    return word+_permutation_word(setting.input_pair)


def _quarter_word():
    return compiler.compile_diagonal([math.pi/2 if i in (1, 2) else 0. for i in range(8)])


def feedback_word():
    """Phase-exact Z on source bit zero; invoked only from common records."""
    return [compiler.pulse('Z', math.pi/2, [0]), compiler.pulse('I', -math.pi/2, [0])]


def probe_words(probe, reference=False):
    """Fresh emitted words, including the independently executed inverse route."""
    _flags(reference)
    setting = settings(probe)
    route = _permutation_word(setting.output_pair)
    return {'preparation': _preparation_word(probe, reference),
            'read_route_inverse': compiler.inverse(route),
            'pointer_read': circuit.readout_word(), 'read_route_forward': route,
            'native_fair_coin': [compiler.pulse('Y', -math.pi/4, [0])],
            'quarter_waiting_phase': _quarter_word(), 'feedback': feedback_word(),
            **{f'{candidate}_repair_coin{coin}': circuit.repair_word(candidate, coin)
               for candidate in ('C', 'E') for coin in (0, 1)},
            **{f'parity_copy_coin{coin}': compiler.cnot(coin, 9) for coin in (0, 1)}}


@lru_cache(maxsize=18)
def _word_hashes(probe, reference):
    return tuple((name, hashlib.sha256(canonical(word).encode()).hexdigest())
                 for name, word in sorted(probe_words(probe, reference).items()))


def word_hashes(probe, reference=False):
    settings(probe); _flags(reference)
    return _word_hashes(probe, reference)


@lru_cache(maxsize=20)
def _route(pair, inverse):
    word = _permutation_word(pair)
    return _immutable(compiler.execute(compiler.inverse(word) if inverse else word, 3))


def _apply_local(matrix, value):
    vector = value.ndim == 1
    columns = value.reshape(512, -1)
    local = columns.reshape(64, 8, -1).transpose(1, 0, 2)
    result = (matrix@local.reshape(8, -1)).reshape(local.shape).transpose(1, 0, 2).reshape(columns.shape)
    return result[:, 0] if vector else result


def preparation(probe, reference=False):
    _flags(reference)
    word = _preparation_word(probe, reference)
    blank = np.zeros(1024 if reference else 512, complex)
    blank[0] = 1
    result = compiler.execute(word, 10 if reference else 9, blank)
    return _immutable(result.reshape(2, 512).T if reference else result)


def preparation_resources(probe, reference=False):
    _flags(reference)
    return _resources(_preparation_word(probe, reference), declared_blank_data_preparations=1,
                      joint_reference_preparations=int(reference))


@lru_cache(maxsize=1)
def _extra_coin():
    return _immutable(compiler.execute([compiler.pulse('Y', -math.pi/4, [0])],
                                      1, np.array([1., 0.])))


@lru_cache(maxsize=4)
def _repair_resources(base, coin):
    coin_word = [compiler.pulse('Y', -math.pi/4, [0])]
    return _resources(coin_word+circuit.repair_word(base, coin), native_repairs=1,
        blank_coin_preparations=1, coin_measurements=1,
        blank_record_bits=circuit.RECORD_BITS[base], record_bit_measurements=circuit.RECORD_BITS[base])


def _parity_branches(coin, state):
    vector = state.ndim == 1
    columns = state.reshape(512, -1)
    blank = np.concatenate((columns, np.zeros_like(columns)))
    output = compiler.execute(compiler.cnot(coin, 9), 10, blank)
    return tuple(block[:, 0] if vector else block for block in (output[:512], output[512:]))


def repair_branches(candidate, state):
    """Actual native/source branches, including private refinement and aborts."""
    candidate, state = _candidate(candidate), _state(state)
    base = 'E' if candidate == 'E' else 'C'
    coin_word = [compiler.pulse('Y', -math.pi/4, [0])]
    output = []
    for branch in recorded.repair_branches(base, state):
        coin, raw, current = branch['coin'], branch['record_code'], branch['current']
        value = branch['amplitude']
        costs = _repair_resources(base, coin)
        if candidate == 'C_phase_quarter' and current == 0:
            word = _quarter_word()
            value = compiler.execute(word, 9, value)
            costs = _sum_resources(costs, _resources(word, waiting_phase_controls=1))
        alternatives = [((), value, costs)]
        if candidate == 'E_collision_half' and not branch['abort']:
            alternatives = []
            choices = ((None, 1.),) if current == 0 else tuple(enumerate(_extra_coin()))
            for choice, factor in choices:
                chosen_cost = costs if choice is None else _sum_resources(costs,
                    _resources(coin_word, blank_coin_preparations=1, coin_measurements=1))
                if choice == 0:
                    alternatives.append(((choice, None), factor*value, chosen_cost))
                else:
                    parity_cost = _sum_resources(chosen_cost, _resources(compiler.cnot(coin, 9),
                        blank_record_bits=1, record_bit_measurements=1))
                    alternatives.extend(((choice, parity), factor*projected, parity_cost)
                                        for parity, projected in enumerate(_parity_branches(coin, value)))
        output.extend(Branch((coin, raw, auxiliary), current, branch['contrast'], branch['abort'],
                             _immutable(value), cost) for auxiliary, value, cost in alternatives)
    return tuple(output)


@lru_cache(maxsize=1)
def _full_pointer_maps():
    # A continuing numerical state can contain tiny spectator amplitudes.
    # Execute every column once rather than erasing those amplitudes or
    # substituting the ideal identity on the detector's complement.
    return tuple((pointer, _immutable(value)) for pointer, value in
                 recorded.readout_branches(np.eye(512, dtype=complex)).items())


def readout_branches(probe, state):
    """Complete Lueders plus/complement read, restoring native coordinates."""
    setting, state = settings(probe), _state(state)
    routed = _apply_local(_route(setting.output_pair, True), state)
    measured = ({pointer: matrix@routed for pointer, matrix in _full_pointer_maps()}
                if np.any(routed[8:]) else recorded.readout_branches(routed))
    return {pointer: _immutable(_apply_local(_route(setting.output_pair, False), value))
            for pointer, value in measured.items()}


@lru_cache(maxsize=9)
def _readout_resources(probe):
    route = _permutation_word(settings(probe).output_pair)
    return _resources(compiler.inverse(route)+circuit.readout_word()+route,
                      blank_pointer_preparations=1, pointer_measurements=1)


def readout_resources(probe):
    settings(probe)
    return _readout_resources(probe)


def _event(candidate, probe, state):
    output = []
    for branch in repair_branches(candidate, state):
        if branch.abort:
            output.append(branch._replace(record=(*branch.record, None)))
        else:
            costs = _sum_resources(branch.resources, readout_resources(probe))
            output.extend(branch._replace(record=(*branch.record, pointer), amplitude=value,
                                          resources=costs)
                          for pointer, value in readout_branches(probe, branch.amplitude).items())
    return tuple(output)


def event_matrices(candidate, probe):
    """Every emitted event Kraus map on all eight active data columns.

    The 504 spectator-complement rows are retained, not clipped. This is
    source execution; no target instrument or normalized branch is inserted.
    """
    _candidate(candidate); settings(probe)
    inputs = np.zeros((512, 8), complex)
    inputs[:8] = np.eye(8)
    return _event(candidate, probe, inputs)


def execute_probe(candidate, probe, reference=False):
    _candidate(candidate); settings(probe); _flags(reference)
    prepared = preparation(probe, reference)
    costs = preparation_resources(probe, reference)
    branches = tuple(branch._replace(resources=_sum_resources(costs, branch.resources))
                     for branch in _event(candidate, probe, prepared))
    return ProbeExecution(candidate, probe, 2 if reference else 1, branches,
                          word_hashes(probe, reference))


def continue_history(candidate, probe, *, feedback=True, reference=False):
    """Two repairs/readouts on one continuing source; no second preparation.

    Feedback sees only first current and pointer, never the native coin,
    contrast, auxiliary refinement record or hidden source state.
    Canonical runtime checkpointing is provided by the experiment adapter.
    """
    _flags(reference, feedback)
    first = execute_probe(candidate, probe, reference)
    output = []
    for branch in first.branches:
        if branch.abort:
            output.append(HistoryBranch((branch.record, None), (branch.current, None),
                                        True, branch.amplitude, branch.resources))
            continue
        state, costs = branch.amplitude, branch.resources
        if feedback and branch.current == 0 and branch.record[-1] == 1:
            state = compiler.execute(feedback_word(), 9, state)
            costs = _sum_resources(costs, _resources(feedback_word(), feedback_phase_controls=1))
        for final in _event(candidate, probe, state):
            output.append(HistoryBranch((branch.record, final.record),
                (branch.current, final.current), final.abort, final.amplitude,
                _sum_resources(costs, final.resources)))
    return HistoryExecution(candidate, probe, first.reference_dimension, feedback,
                            tuple(output), first.word_hashes)

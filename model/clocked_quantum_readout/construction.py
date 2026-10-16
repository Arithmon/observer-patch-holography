"""Complete native preparation, stopping and phase-read instrument.

The ideal balanced service is executed numerically on its invariant native
sector. No normalized successful state or target POVM enters execution.
Source compilation is a separately accounted implementation of that service.
"""
from __future__ import annotations

from fractions import Fraction as F
from dataclasses import dataclass

import numpy as np

from ..native_interface import construction as interface
from ..native_relational_readout import construction as relational
from ..primitive_observer import runtime
from ..clocked_common_observer import clock, construction as work


YIELD = 2*interface.SUCCESS
CAPS = (4, 16, 64)


def preparation_program(law):
    law = clock.quantum.candidate(law)
    herald, _ = interface.herald_and_pointer(law)
    identity = runtime.Instrument(kraus=((np.eye(66),),))
    balanced = runtime.Instrument(kraus=((interface.mixer(),),),
        resources=(('balanced_mixer_calls', 1),),
        binding={'ideal_service_surrogate': True,
                 'source_implementation_error_charged_separately': True})
    return runtime.Program('heralded_native_clock_preparation', 66, np.eye(66)/66,
        (runtime.Event('fixed_seam_current_herald', instruments=(herald,)),
         runtime.Event('balanced_if_heralded', instruments=(identity, balanced, balanced), control=0)),
        preparation_resources=(('native_preparations', 1),),
        binding={'stationary_diagonal_native_input_supplied': True,
                 'herald_flags': ['failed', 'plus_phase', 'minus_phase'],
                 'no_conditionally_normalized_execution': True})


def records(cap):
    clock._horizon(cap)
    return ((0, 0, 0, -1),)+tuple(
        (herald, n, current, bit)
        for herald in (1, 2) for n in range(1, cap+1)
        for current in (-1, 1) for bit in (0, 1)) + (
        (1, cap, 0, -1), (2, cap, 0, -1))


@dataclass(frozen=True)
class NativeExecution:
    law: str
    cap: int
    preparation: runtime.History
    # Frozen ordered pairs; every numerical density has its original entries.
    outputs: tuple


def execute_native(law, cap):
    law = clock.quantum.candidate(law)
    clock._horizon(cap)
    prepared = runtime.execute(preparation_program(law))
    initial = {branch.record[0]: branch.density for branch in prepared.branches}
    outputs = {(0, 0, 0, -1): initial[0]}
    _, pointer = interface.herald_and_pointer(law)
    balanced = runtime.Instrument(kraus=((interface.mixer(),),))
    source = clock.source(law).instrument
    for herald in (1, 2):
        state = initial[herald]
        for n in range(1, cap+1):
            negative, state, positive = runtime._apply(source, state, 1)
            for current, stopped in ((-1, negative), (1, positive)):
                analyzed = runtime._apply(balanced, stopped, 1)[0]
                phases = runtime._apply(pointer, analyzed, 1)
                for bit in (0, 1):
                    # Match/opposite decoding uses the committed herald sign.
                    outputs[herald, n, current, bit] = runtime._immutable(
                        phases[bit if herald == 1 else 1-bit])
        outputs[herald, cap, 0, -1] = runtime._immutable(state)
    return NativeExecution(law, cap, prepared,
                           tuple((r, outputs[r]) for r in records(cap)))


def _balanced_exact(state):
    """Exact rational conjugation on this source's invariant support.

    A general cross term between the active corner and its complement would
    require irrational coefficients. It is refused, not silently discarded.
    Every actual stopped state in this preparation has block-diagonal support.
    """
    columns = {}
    for left, right in interface.phase_pairs():
        columns[left] = ((left, 1), (right, 1))
        columns[right] = ((left, 1), (right, -1))
    result = {}
    for (i, j), value in state.items():
        active = i in columns
        if active != (j in columns):
            if value:
                raise ValueError('cross-corner coherence needs an irrational exact domain')
            continue
        scale = F(1, 2) if active else F(1)
        for a, x in columns.get(i, ((i, 1),)):
            for b, y in columns.get(j, ((j, 1),)):
                result[a, b] = result.get((a, b), F(0))+value*x*y*scale
    return {key: value for key, value in result.items() if value}


def _pointer_exact(state, law):
    outputs = ({}, {})
    for branch in interface.fixed_maps(law):
        bit = int(2*branch.current+2*branch.coin-1 < 0)
        mapping = dict(branch.mapping)
        for (i, j), value in state.items():
            if i in mapping and j in mapping:
                key = mapping[i], mapping[j]
                outputs[bit][key] = outputs[bit].get(key, F(0))+value/2
    return tuple({key: value for key, value in state.items() if value} for state in outputs)


def exact_outputs(law, cap):
    law = clock.quantum.candidate(law)
    clock._horizon(cap)
    failed = {}
    for branch in interface.fixed_maps(law):
        if branch.current == 0:
            for _, output in branch.mapping:
                failed[output, output] = failed.get((output, output), F(0))+F(1, 132)
    outputs = {(0, 0, 0, -1): failed}
    for herald, sign in ((1, 1), (2, -1)):
        stopped = clock.exact_outputs(law, cap, sign=sign)
        for (n, current), state in stopped.items():
            state = {key: interface.SUCCESS*value for key, value in state.items()}
            if current:
                phases = _pointer_exact(_balanced_exact(state), law)
                for bit in (0, 1):
                    outputs[herald, n, current, bit] = phases[bit if herald == 1 else 1-bit]
            else:
                outputs[herald, cap, 0, -1] = state
    ordered = {r: outputs[r] for r in records(cap)}
    if sum(map(clock._trace, ordered.values()), F(0)) != 1:
        raise ArithmeticError('complete herald/stop/phase instrument lost mass')
    return ordered


def validate_native(execution):
    if type(execution) is not NativeExecution:
        raise ValueError('a complete native execution required')
    law = clock.quantum.candidate(execution.law)
    cap = clock._horizon(execution.cap)
    if (type(execution.outputs) is not tuple
            or any(type(row) is not tuple or len(row) != 2 for row in execution.outputs)
            or any(type(r) is not tuple or len(r) != 4
                   or any(type(x) is not int for x in r) for r, _ in execution.outputs)
            or tuple(r for r, _ in execution.outputs) != records(cap)):
        raise ValueError('every ordered native outcome including failure and timeout is required')
    prep = runtime.validate_history(execution.preparation)
    expected = preparation_program(law)
    if (prep.program.digest != expected.digest or prep.cursor != 2
            or prep.preparation_mode != 'declared' or prep.reference_dimension != 1):
        raise ValueError('the complete declared diagonal-input preparation is required')
    snapshots = []
    for record, array in execution.outputs:
        array = clock._array(array, (66, 66))
        snapshots.append((record, np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)))
    return NativeExecution(law, cap, prep, tuple(snapshots))


def certify_native(execution):
    execution = validate_native(execution)
    targets = exact_outputs(execution.law, execution.cap)
    results = {}
    for record, array in execution.outputs:
        target = targets[record]
        support = set(target) | {(int(i), int(j)) for i, j in zip(*array.nonzero(), strict=True)}
        squared = F(0)
        for i, j in support:
            real, imag = clock.exact._component(array[i, j])
            squared += (real-target.get((i, j), F(0)))**2+imag**2
        error = clock.exact._sqrt_interval(66*squared)[1]/2
        results[record] = {'error': error, 'mass': clock._trace(target)}
    return {'error': sum((r['error'] for r in results.values()), F(0)), 'per_record': results}


def output_rows(execution, work_execution):
    """Complete native/work factors, without a floating tensor multiplication."""
    execution, w = validate_native(execution), work.validate_work(work_execution)
    if execution.cap > w.horizon:
        raise ValueError('the work history has not executed the requested cap')
    histories = w.successful()
    rows = []
    for (h, n, j, bit), density in execution.outputs:
        selected = histories[n, j] if j else w.prefixes[n]
        for branch in selected.branches:
            rows.append(((h, n, j, branch.record, bit), density, branch.density))
    if len(rows) != 128*execution.cap+3:
        raise ArithmeticError('a complete joint outcome was omitted')
    return iter(tuple(rows))

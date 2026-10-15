"""Execute a native first-current clock and the retained interacting work.

The clock's first signed current schedules the existing four-read continuation.
The same work state advances at every slot before the latch, even native-idle
slots in the idling control. Its declared
Hamiltonian coordinate per slot is not inferred from the native count. Native
and work preparation are independent; their joint output is stored in factors,
without forming a prohibitively large tensor density or selecting successes.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np

from tools import strict_json
from ..common_source import exact_output, protocol
from ..primitive_observer import common, runtime
from . import clock


DIMENSION = 2
STEPS = 128
DELTA = F(1, 4)
CAPS = (4, 16, 64, 128)


def _parameters(horizon, steps, delta):
    clock._horizon(horizon)
    if type(steps) is not int or not 1 <= steps <= 4096:
        raise ValueError('an integer source split count from one through 4096 is required')
    if type(delta) is not F or not 0 < delta <= 1:
        raise ValueError('an exact positive model duration at most one is required')


def compile_wait(horizon, *, steps=STEPS, delta=DELTA):
    """One preparation and repeated source evolution; no intervening reset."""
    _parameters(horizon, steps, delta)
    canonical = protocol.continuation_program(DIMENSION, steps, native_current=0)
    prepared = common.compile_program(canonical)
    maps, cost = protocol._operation_maps(canonical, protocol.Operation('evolve', delta), None, 'source')
    if len(maps) != 1:
        raise ArithmeticError('one source-contracted waiting unitary is required')
    operation = runtime.Instrument(kraus=((maps[0],),), resources=tuple(cost),
        binding={'kind': 'existing_source_interacting_evolution',
                 'model_duration_per_controller_slot': str(delta), 'split_count': steps})
    return runtime.Program('work_waiting_for_first_native_current', canonical.work_dimension,
        prepared.initial_density,
        tuple(runtime.Event(f'wait_{n}', instruments=(operation,)) for n in range(1, horizon+1)),
        prepared.preparation_resources,
        {'occupation_dimension': DIMENSION, 'attempt_cap': horizon,
         'model_duration_per_controller_slot': str(delta),
         'time_is_a_supplied_model_coordinate': True,
         'one_initial_work_preparation': True})


@dataclass(frozen=True)
class WorkExecution:
    horizon: int
    steps: int
    delta: F
    prefixes: tuple[runtime.History, ...]
    # Ordered n=1,...,N, current=-1,+1. All 16 suffix records are retained.
    successes: tuple[tuple[int, int, runtime.History], ...]

    def successful(self):
        return {(n, current): history for n, current, history in self.successes}


def execute_work(horizon, *, steps=STEPS, delta=DELTA, initial=None):
    """Shared-runtime execution; a host prefix fork is not physical copying.

    Every slot extends the preceding unread state. The alternative suffixes
    are the mutually exclusive first-event branches of the complete process.
    Supplying work/reference input replaces, and does not charge, preparation.
    """
    waiting = compile_wait(horizon, steps=steps, delta=delta)
    suffixes = {j: common.compile_program(protocol.continuation_program(
        DIMENSION, steps, native_current=j)) for j in (-1, 1)}
    state = runtime.execute(waiting, initial=initial, stop_after=0)
    prefixes, successes = [state], []
    for n in range(1, horizon+1):
        state = runtime.execute(waiting, checkpoint=state, stop_after=n)
        prefixes.append(state)
        for current in (-1, 1):
            retained = runtime.branch_experiment(state, suffixes[current].events,
                name='native_current_triggered_interacting_observer',
                binding={'first_current_slot': n, 'signed_current': current,
                         'work_preparation_repeated': False,
                         'suffix': suffixes[current].specification()})
            successes.append((n, current, retained))
    return WorkExecution(horizon, steps, delta, tuple(prefixes), tuple(successes))


def validate_work(work):
    """Detach a complete, source-bound in-process execution family.

    Seals check integrity, not authentication of an entirely fabricated and
    resealed object. Saved untrusted histories must pass the runtime's replay
    importer. Resource menus are reconstructed independently of their seals.
    """
    if type(work) is not WorkExecution:
        raise ValueError('a typed executed work family is required')
    _parameters(work.horizon, work.steps, work.delta)
    expected_keys = tuple((n, j) for n in range(1, work.horizon+1) for j in (-1, 1))
    if (type(work.prefixes) is not tuple or len(work.prefixes) != work.horizon+1
            or type(work.successes) is not tuple or len(work.successes) != len(expected_keys)
            or any(type(row) is not tuple or len(row) != 3
                   or type(row[0]) is not int or type(row[1]) is not int
                   for row in work.successes)
            or tuple((n, j) for n, j, _ in work.successes) != expected_keys):
        raise ValueError('every ordered waiting prefix and signed successful work history is required')
    waiting = compile_wait(work.horizon, steps=work.steps, delta=work.delta)
    wait_digest = waiting.digest
    prefixes = tuple(runtime.validate_history(history) for history in work.prefixes)
    initial = prefixes[0]
    preparation = waiting.preparation_resources if initial.preparation_mode == 'declared' else ()
    wait_cost = waiting.events[0].instruments[0].resources
    cost = preparation
    for index, prefix in enumerate(prefixes):
        if (prefix.program.digest != wait_digest or prefix.cursor != index
                or prefix.reference_dimension != initial.reference_dimension
                or prefix.preparation_mode != initial.preparation_mode
                or not np.array_equal(prefix.initial_density, initial.initial_density)
                or prefix.records != ((),)):
            raise ValueError('waiting prefixes must retain the same source program and initial work/reference state')
        if (prefix.branches[0].resources != cost or prefix.expanded_tree_resources != cost):
            raise ValueError('waiting resource counts must charge one preparation and every executed slot')
        cost = runtime.add_resources(cost, wait_cost)
    suffixes = {j: common.compile_program(protocol.continuation_program(
        DIMENSION, work.steps, native_current=j)) for j in (-1, 1)}
    checked = []
    for n, current, value in work.successes:
        history = runtime.validate_history(value)
        suffix = suffixes[current]
        expected = runtime.Program(
            'native_current_triggered_interacting_observer', waiting.work_dimension,
            waiting.initial_density, (*waiting.events[:n], *suffix.events),
            waiting.preparation_resources,
            {'prefix_declaration': strict_json.loads(waiting.binding),
             'alternative_continuation': {
                 'first_current_slot': n, 'signed_current': current,
                 'work_preparation_repeated': False, 'suffix': suffix.specification()},
             'host_prefix_reuse_is_physical_copying': False})
        if (history.program.digest != expected.digest or history.cursor != len(expected.events)
                or history.reference_dimension != initial.reference_dimension
                or history.preparation_mode != initial.preparation_mode
                or not np.array_equal(history.initial_density, initial.initial_density)):
            raise ValueError('successful work must continue its waiting prefix through the declared signed source suffix')
        costs = {(): prefixes[n].branches[0].resources}
        expanded = prefixes[n].expanded_tree_resources
        for event in suffix.events:
            following = {}
            for record, previous in costs.items():
                instrument = event.instruments[0 if event.control is None else record[event.control]]
                expanded = runtime.add_resources(expanded, instrument.resources)
                for outcome in range(instrument.outcomes):
                    output = (*record, outcome) if event.visible else record
                    following[output] = runtime.add_resources(previous, instrument.resources)
            costs = following
        if (history.records != tuple(costs)
                or history.expanded_tree_resources != expanded
                or any(branch.resources != costs[branch.record] for branch in history.branches)):
            raise ValueError('complete suffix resources must retain every outcome and charge no new preparation')
        checked.append((n, current, history))
    return WorkExecution(work.horizon, work.steps, work.delta, prefixes, tuple(checked))


def output_rows(work, native_outputs, *, cap=None):
    """Stream complete labelled density factors without a normalized selection.

    Success labels are (slot, signed current, four work bits), timeout is
    (cap, 0, ()). Native seam/coin labels and hidden repair labels are summed
    by the declared signed-current instrument, not exposed as archived data.
    """
    work = validate_work(work)
    cap = work.horizon if cap is None else cap
    clock._horizon(cap)
    if cap > work.horizon:
        raise ValueError('the requested cap has not been executed')
    expected = {(n, j) for n in range(1, cap+1) for j in (-1, 1)} | {(cap, 0)}
    if (type(native_outputs) is not dict
            or any(type(record) is not tuple or len(record) != 2
                   or any(type(value) is not int for value in record) for record in native_outputs)
            or set(native_outputs) != expected):
        raise ValueError('all signed stopping outcomes and the timeout are required')
    native = {}
    for record, value in native_outputs.items():
        snapshot = exact_output._array(value, (66, 66))
        # Check the numerical density domain without using the validator's
        # cast as output: the original finite entries remain unchanged.
        runtime._density(snapshot, 66, normalized=False)
        native[record] = np.frombuffer(snapshot.tobytes(), dtype=snapshot.dtype).reshape(snapshot.shape)
    mass = sum((np.trace(value) for value in native.values()), 0j)
    if not np.isfinite(mass) or abs(mass-1) > runtime.TOLERANCE:
        raise ValueError('the complete native output must retain its total mass')
    rows = []
    for n, current, history in work.successes:
        if n <= cap:
            for branch in history.branches:
                rows.append(((n, current, branch.record), native[n, current], branch.density))
    rows.append(((cap, 0, ()), native[cap, 0], work.prefixes[cap].branches[0].density))
    # Validate and snapshot before exposing a lazy iterator. Later caller
    # mutations cannot alter a quantum factor or remove a retained record.
    return iter(tuple(rows))


def observations(work, native_outputs, *, cap=None):
    """Unconditional, complete-law probabilities; timeout is not discarded."""
    rows = []
    for record, native_density, work_density in output_rows(work, native_outputs, cap=cap):
        probability = float((np.trace(native_density)*np.trace(work_density)).real)
        if not math.isfinite(probability) or probability < -runtime.TOLERANCE:
            raise ArithmeticError('nonfinite or negative complete output probability')
        rows.append((record, probability))
    total = math.fsum(p for _, p in rows)
    if abs(total-1) > runtime.TOLERANCE:
        raise ArithmeticError('native-clocked complete output lost mass')
    return {'mass': total, 'outcomes': len(rows),
            'timeout_probability': math.fsum(p for (_, j, _), p in rows if j == 0),
            'final_gauge_one_and_success': math.fsum(p for (_, j, bits), p in rows if j and bits[3] == 1),
            'different_smooth_bits_and_success': math.fsum(p for (_, j, bits), p in rows if j and bits[1] != bits[2]),
            'first_slot_gauge_one': math.fsum(p for (n, j, bits), p in rows if n == 1 and j and bits[3] == 1)}

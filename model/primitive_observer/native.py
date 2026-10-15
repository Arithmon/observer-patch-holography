"""Close the native repair/work/record/native repair loop on the shared runtime.

The first executed native current selects the work preparation kick. The
second smooth work record selects the final native seam. Protected and passive
logical memories use the existing executed encoder and recovery words. Every
raw record, abort and padding output is retained, including zero-mass branches.

Work is independently prepared. A native input may be entangled with its
untouched reference, but an initially correlated native/work input is outside
this executable contract. Factorized output avoids a dense 512*work density.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import math

import numpy as np

from ..common_source import encoded_memory as memory
from ..common_source import native_routing as routing
from ..common_source import open_environment as environment
from ..common_source import open_native as original
from ..common_source import protocol, quantum_native as source
from . import common, runtime


@dataclass(frozen=True)
class NativeHistory:
    law: str
    dimension: int
    steps: int
    native_mode: str
    environment: environment.Environment
    work_histories: tuple[runtime.History, ...]
    branches: tuple[original.OpenNativeBranch, ...]
    seal: str = ''


def _rows(law, mode, workers):
    by_current = dict(zip((-1, 0, 1), workers, strict=True))
    prepared = workers[1].initial_density
    rows, snapshots = [], {}
    def matrix(value):
        if id(value) not in snapshots:
            snapshots[id(value)] = source._immutable(source._array(
                value, (512, source.input_dimension(mode)), 'compiled native map'))
        return snapshots[id(value)]
    for c1, raw1, current, abort, first, endings in source._native_tree(law, mode):
        if (type(c1) is not int or type(raw1) is not int or type(abort) is not bool
                or current != source._record_current(law, c1, raw1)
                or abort != (current is None)):
            raise ValueError('native source tree has an invalid current/abort record')
        if abort:
            rows.append(original.OpenNativeBranch((c1, raw1, (), None, None), True,
                                                   matrix(first), prepared))
            continue
        for middle in by_current[current].branches:
            for c2, raw2, stopped, complete in endings[middle.record[2]]:
                if (type(c2) is not int or type(raw2) is not int or type(stopped) is not bool
                        or stopped != (source._record_current(law, c2, raw2) is None)):
                    raise ValueError('native continuation has an invalid final raw outcome')
                rows.append(original.OpenNativeBranch((c1, raw1, middle.record, c2, raw2),
                                                       stopped, matrix(complete), middle.density))
    expected = source._menu(law)
    if len(rows) != len(expected) or {row.record: row.abort for row in rows} != expected:
        raise ValueError('native source tree omitted or relabelled a complete raw output')
    return tuple(rows)


def _seal(history):
    header = {'law': history.law, 'dimension': history.dimension, 'steps': history.steps,
              'native_mode': history.native_mode, 'environment': history.environment.specification(),
              'work_seals': [worker.seal for worker in history.work_histories],
              'records': [(branch.record, branch.abort) for branch in history.branches]}
    digest = hashlib.sha256(runtime._json(header).encode())
    cache = {}
    for branch in history.branches:
        for array in (branch.native_map, branch.work_density):
            if id(array) not in cache:
                cache[id(array)] = runtime._array_digest(array).encode()
            digest.update(cache[id(array)])
    return digest.hexdigest()


def validate_history(history):
    """Validate and detach a trusted in-memory execution; no work is reset."""
    if type(history) is not NativeHistory:
        raise ValueError('a typed native observer history required')
    law, mode = source._law(history.law), source._mode(history.native_mode)
    base = protocol.continuation_program(history.dimension, history.steps)
    env = environment.validated(base, history.environment)
    if type(history.work_histories) is not tuple or len(history.work_histories) != 3:
        raise ValueError('all three native-current work continuations required')
    workers = tuple(runtime.validate_history(worker) for worker in history.work_histories)
    prepared = workers[1].initial_density
    for current, worker in zip((-1, 0, 1), workers, strict=True):
        expected = common.compile_program(protocol.continuation_program(
            history.dimension, history.steps, native_current=current), env)
        if (worker.program.digest != expected.digest or worker.cursor != len(expected.events)
                or worker.reference_dimension != 1 or worker.preparation_mode != 'declared'
                or not np.array_equal(worker.initial_density, prepared)):
            raise ValueError('native attachment requires the complete independently prepared canonical work')
    expected_rows = _rows(law, mode, workers)
    if type(history.branches) is not tuple or len(history.branches) != len(expected_rows):
        raise ValueError('all complete native/work outputs required')
    checked = set()
    for actual, expected in zip(history.branches, expected_rows, strict=True):
        if (type(actual) is not original.OpenNativeBranch or actual.record != expected.record
                or type(actual.record) is not tuple or type(actual.record[2]) is not tuple
                or any(type(x) is not int for x in (*actual.record[:2], *actual.record[2]))
                or any(x is not None and type(x) is not int for x in actual.record[3:])
                or type(actual.abort) is not bool or actual.abort != expected.abort):
            raise ValueError('native/raw records and abort outputs must remain complete')
        for left, right in ((actual.native_map, expected.native_map),
                            (actual.work_density, expected.work_density)):
            key = (id(left), id(right))
            if key not in checked:
                if not isinstance(left, np.ndarray) or not np.array_equal(left, right):
                    raise ValueError('native map or mixed work output detached from its executed history')
                checked.add(key)
    if not history.seal or _seal(history) != history.seal:
        raise ValueError('native observer integrity digest changed')
    return replace(history, environment=env, work_histories=workers, branches=expected_rows)


def execute(env=None, law='C', dimension=2, steps=8, *, native_mode='protected_qubit'):
    """Execute native input -> current -> shared work runtime -> native seam."""
    law, native_mode = source._law(law), source._mode(native_mode)
    base = protocol.continuation_program(dimension, steps)
    env = environment.validated(base, environment.Environment(()) if env is None else env)
    workers = tuple(runtime.execute(common.compile_program(protocol.continuation_program(
        dimension, steps, native_current=current), env)) for current in (-1, 0, 1))
    history = NativeHistory(law, dimension, steps, native_mode, env, workers,
                            _rows(law, native_mode, workers))
    return replace(history, seal=_seal(history))


def output(history, initial):
    """Stream (record, native/reference amplitude, mixed work density).

    A row denotes |amplitude><amplitude| tensor work_density; neither factor
    is normalized conditionally. Bind the complete output before returning a
    lazy iterator, so later mutation of caller envelopes cannot change it.
    """
    history = validate_history(history)
    initial = source._input(initial, source.input_dimension(history.native_mode))
    rows = tuple((branch.record, branch.native_map, branch.work_density) for branch in history.branches)
    return ((record, source._immutable(matrix@initial), density) for record, matrix, density in rows)


def audit(history):
    """All-input native completeness and independent endpoint-map diagnostics."""
    history = validate_history(history)
    size = source.input_dimension(history.native_mode)
    targets = (source._ideal_native_tree(history.law) if history.native_mode == 'full_carrier'
               else source._ideal_encoded_tree(history.law, history.native_mode))
    zero = np.zeros((512, size), complex)
    padding = np.array([i for i in range(512) if i not in set(source.CODES)])
    effects, residuals = {}, {}
    effect = np.zeros((size, size), complex)
    maximum = abort_mass = padding_mass = 0.
    for branch in history.branches:
        c1, raw1, bits, c2, raw2 = branch.record
        target = zero if not bits else targets.get((c1, raw1, bits[2], c2, raw2), zero)
        key = id(branch.native_map), id(target)
        if key not in residuals:
            residuals[key] = float(np.linalg.norm(branch.native_map-target, 'fro'))
        maximum = max(maximum, residuals[key])
        key = id(branch.native_map)
        if key not in effects:
            gram = branch.native_map.conj().T@branch.native_map
            effects[key] = (gram, float(np.trace(gram).real)/size,
                            float(np.linalg.norm(branch.native_map[padding], 'fro')**2)/size)
        gram, trace, pad = effects[key]
        weight = float(np.trace(branch.work_density).real)
        effect += weight*gram
        abort_mass += weight*trace*branch.abort
        padding_mass += weight*pad
    error = float(np.linalg.norm(effect-np.eye(size), 2))
    if any(not math.isfinite(value) or value > source.TOLERANCE
           for value in (maximum, error, abort_mass, padding_mass)):
        raise ArithmeticError('native/runtime observer failed its complete independent endpoint contract')
    return {'raw_outcome_count_including_abort': len(history.branches),
            'first_abort_outcome_count': sum(not b.record[2] for b in history.branches),
            'final_abort_outcome_count': sum(b.abort and bool(b.record[2]) for b in history.branches),
            'maximum_native_map_frobenius_error_diagnostic': maximum,
            'prepared_work_native_input_completeness_error_diagnostic': error,
            'maximally_mixed_native_abort_mass_diagnostic': abort_mass,
            'maximally_mixed_native_padding_mass_diagnostic': padding_mass,
            'all_native_output_rows': 512, 'work_output_is_mixed_density': True,
            'native_current_selects_work_kick': True, 'second_smooth_record_selects_native_seam': True,
            'arbitrary_initial_native_work_correlations_executed': False,
            'floating_or_physical_error_certified': False}


def _phase_distance(history):
    flip = np.array([[0., 1.], [1., 0.]])
    cache, terms = {}, []
    for branch in history.branches:
        key = id(branch.native_map)
        if key not in cache:
            _, triangular = np.linalg.qr(branch.native_map, mode='reduced')
            delta = triangular@flip@triangular.conj().T
            cache[key] = float(np.abs(np.linalg.eigvalsh(delta)).sum())/2
        terms.append(float(np.trace(branch.work_density).real)*cache[key])
    return math.fsum(terms)


def memory_witness(active, passive):
    """Instantiate logical memory using the actual noisy work-record weights."""
    active, passive = validate_history(active), validate_history(passive)
    if (active.native_mode != 'protected_qubit' or passive.native_mode != 'passive_qubit'
            or any(getattr(active, key) != getattr(passive, key)
                   for key in ('law', 'dimension', 'steps', 'environment'))):
        raise ValueError('matched active and passive encoded observers required')
    probabilities = [math.fsum(branch.probability for branch in worker.branches if branch.record[2] == 0)
                     for worker in active.work_histories]
    qbar = (probabilities[0]+2*probabilities[1]+probabilities[2])/4
    distances = {'active': _phase_distance(active), 'passive': _phase_distance(passive)}
    bell = np.eye(2, dtype=complex)/math.sqrt(2)
    target = source._logical_code()@bell
    fidelity = math.fsum(float(np.trace(branch.work_density).real)*abs(
        np.vdot(target, branch.native_map@bell))**2 for branch in active.branches)
    if (any(not math.isfinite(value) for value in (*probabilities, qbar, fidelity, *distances.values()))
            or abs(distances['active']-1) > source.TOLERANCE
            or abs(distances['passive']-(1-qbar)) > source.TOLERANCE
            or abs(fidelity-1) > source.TOLERANCE):
        raise ArithmeticError('shared work history failed its attached encoded-memory identity')
    return {'second_smooth_plus_probability_diagnostic': qbar,
            'pre_read_phase_distances_diagnostic': distances,
            'active_encoded_bell_fidelity_diagnostic': fidelity,
            'passive_distance_formula': '1 - qbar',
            'passive_any_decoder_error_lower_diagnostic': qbar/2,
            'existing_work_only_CPTP_memory_theorem_instantiated': True,
            'native_storage_or_record_fault_tolerance_certified': False,
            'terminal_destructive_read_executed': False}


def resources(history):
    """Count each compatible realized path, including native code services."""
    history = validate_history(history)
    keys = tuple(name for name, _ in history.work_histories[0].branches[0].resources)
    empty = dict.fromkeys(keys, 0)
    encoder = memory.encoding_resources() if history.native_mode != 'full_carrier' else {}
    first = routing.routed_repair_resources(history.law, (0, 1))
    fixed = {key: encoder.get(key, 0)+first.get(key, 0) for key in keys}
    work = {current: {b.record: dict(b.resources) for b in worker.branches}
            for current, worker in zip((-1, 0, 1), history.work_histories, strict=True)}
    prepared = dict(history.work_histories[1].program.preparation_resources)
    final = {bit: routing.routed_repair_resources(history.law, source.final_seam(history.native_mode, bit))
             for bit in (0, 1)}
    cache = {}
    def recovery(seam, coin, raw):
        if history.native_mode != 'protected_qubit':
            return empty
        key = seam, coin, raw
        if key not in cache:
            cache[key] = memory.recovery_resources(history.law, seam, coin, raw)
        return cache[key]
    maxima, witness = dict.fromkeys(keys, 0), None
    for branch in history.branches:
        c1, raw1, bits, c2, raw2 = branch.record
        if not bits:
            costs = {key: fixed[key]+prepared[key] for key in keys}
        else:
            current = source._record_current(history.law, c1, raw1)
            seam = source.final_seam(history.native_mode, bits[2])
            additions = (recovery((0, 1), c1, raw1), final[bits[2]],
                         empty if branch.abort else recovery(seam, c2, raw2))
            costs = {key: fixed[key]+work[current][bits][key]+sum(s.get(key, 0) for s in additions)
                     for key in keys}
        if witness is None or costs['native_pulses'] > maxima['native_pulses']:
            witness = branch.record
        maxima = {key: max(maxima[key], costs[key]) for key in keys}
    return {'maximum_realized_path_resources': maxima,
            'native_pulse_maximizing_record': [witness[0], witness[1], list(witness[2]), witness[3], witness[4]],
            'maximum_resources_are_componentwise_over_compatible_paths': True,
            'native_input_qubits': 9,
            'supplied_unknown_logical_qubits': int(history.native_mode != 'full_carrier'),
            'native_encoder_count': int(history.native_mode != 'full_carrier'),
            'maximum_native_coin_initializations': 2, 'maximum_native_record_initializations': 2,
            'retained_classical_bits_including_two_status_bits': 2*(1+source.circuit.RECORD_BITS[history.law])+4+2,
            'native_logical_storage_through_work_required': history.native_mode != 'full_carrier',
            'diagnostic_input_preparation_and_terminal_read_pulses': 0,
            'environment_draws_are_observer_records': False,
            'source_counts_include_every_repeated_invocation': True,
            'host_tree_storage_or_physical_service_errors_counted': False}


def compatibility_report(history, initial=None):
    """Compare every full factorized output with the older independent engine."""
    history = validate_history(history)
    expected = original.execute(history.environment, history.law, history.dimension, history.steps,
                                native_mode=history.native_mode)
    size = source.input_dimension(history.native_mode)
    if initial is None:
        initial = np.eye(2)/math.sqrt(2) if size == 2 else np.eye(size)[:, 0]
    initial = source._input(initial, size)
    if len(history.branches) != len(expected.branches):
        raise ArithmeticError('native observer comparison omitted raw branches')
    cache, terms, maximum = {}, [], 0.
    for left, right in zip(history.branches, expected.branches, strict=True):
        if left.record != right.record or left.abort != right.abort or not np.array_equal(left.native_map, right.native_map):
            raise ArithmeticError('native observer comparison changed a native operation or raw outcome')
        key = id(left.work_density), id(right.work_density)
        if key not in cache:
            delta = left.work_density-right.work_density
            cache[key] = float(np.abs(np.linalg.eigvalsh((delta+delta.conj().T)/2)).sum())/2
        maximum = max(maximum, cache[key])
        amplitude = left.native_map@initial
        terms.append(float(np.vdot(amplitude, amplitude).real)*cache[key])
    distance = math.fsum(terms)
    if not math.isfinite(distance) or distance > runtime.TOLERANCE:
        raise ArithmeticError('shared runtime changed a complete native/reference/work output')
    if resources(history) != original.resources(expected):
        raise ArithmeticError('shared native observer changed a compatible path resource ledger')
    return {'raw_outcome_count_including_abort': len(history.branches),
            'native_input_dimension': size, 'reference_dimension': initial.shape[1],
            'maximum_work_block_trace_distance': maximum,
            'complete_reference_output_trace_distance': distance,
            'every_native_map_and_raw_record_agrees': True,
            'source_resources_agree': True,
            'comparison': 'independent_existing_open_native_interpreter_on_the_same_source_operations',
            'prepared_work_independent_of_native_reference_input': True}

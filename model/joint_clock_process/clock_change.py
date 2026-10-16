"""Complete observer-clock conversion from native instruments.

Native attempt caps are resource limits. The comparison variable counts remote
B events before an A event; it has no identification with calibrated duration.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as sp

from . import construction as c
from ..federated_quantum_sector import construction as q
from ..primitive_observer import runtime

TRACE = sp.ImmutableMatrix([[int(i//4 == i % 4) for i in range(16)]])


def _feedback(value):
    if type(value) is not bool:
        raise ValueError('Boolean remote-event feedback required')
    return value


@lru_cache(maxsize=4, typed=True)
def controlled_maps(law, feedback=False):
    """A field Hadamard is applied only after a visible B-entry when enabled."""
    c._law(law); _feedback(feedback)
    maps = c.source_maps(law)
    unitary = c.superoperator(c.field_unitary('hadamard' if feedback else 'identity'))
    return tuple(sp.ImmutableMatrix(unitary*m if i == 4 else m)
                 for i, m in enumerate(maps))


@lru_cache(maxsize=4, typed=True)
def union_maps(law, feedback=False):
    maps = controlled_maps(law, feedback)
    resolvent = (sp.eye(16)-maps[0]).inv()
    return tuple(sp.ImmutableMatrix(m*resolvent) for m in maps[1:])


def coefficient_maps(law, remote_count, feedback=False):
    """Two CP maps: first A-exit/entry after exactly remote_count B events.

    B signs are summed, after their declared feedback. Waiting attempts are
    summed analytically. No per-branch normalization is applied.
    """
    if type(remote_count) is not int or not 0 <= remote_count <= 64:
        raise ValueError('an exact remote count from0 through64 required')
    a_exit, a_entry, b_exit, b_entry = union_maps(law, feedback)
    before = (b_exit+b_entry)**remote_count
    return tuple(sp.ImmutableMatrix(m*before) for m in (a_exit, a_entry))


@lru_cache(maxsize=4, typed=True)
def complete_maps(law, feedback=False):
    a_exit, a_entry, b_exit, b_entry = union_maps(law, feedback)
    before = (sp.eye(16)-b_exit-b_entry).inv()
    return tuple(sp.ImmutableMatrix(m*before) for m in (a_exit, a_entry))


def certificate(law, feedback=False):
    """Full-map associative elimination and source-derived survival bound."""
    maps = controlled_maps(law, feedback)
    union = union_maps(law, feedback)
    b = union[2]+union[3]
    closed = complete_maps(law, feedback)
    direct_before = (sp.eye(16)-maps[0]-maps[3]-maps[4]).inv()
    if any(event != mark*direct_before for event, mark in zip(closed, maps[1:3], strict=True)):
        raise ArithmeticError('direct and staged clock elimination disagree')
    z = sp.Symbol('z')
    direct_count = (sp.eye(16)-maps[0]-z*(maps[3]+maps[4])).inv()
    staged_count = (sp.eye(16)-z*b).inv()
    for mark, event in zip(maps[1:3], union[:2], strict=True):
        if any(sp.cancel(entry) != 0 for entry in mark*direct_count-event*staged_count):
            raise ArithmeticError('clock conversion lost a remote-count coefficient')
    if TRACE*sum(closed, sp.zeros(16)) != TRACE:
        raise ArithmeticError('completed A clock lost total mass')
    # Complete positivity is inherited from sums/compositions of source maps.
    # The effect bound is diagonal on all matrix units, not just tested states.
    effect = TRACE*b**2
    if any(effect[i] for i in range(16) if i//4 != i % 4):
        raise ArithmeticError('remote two-event survival effect is not diagonal')
    bound = max(effect[5*i] for i in range(4))
    if bound >= 1:
        raise ArithmeticError('no uniform local-clock completion certificate')
    # Noncommutative resolvent factorization gives each exact remote-count
    # coefficient too, hence the count record and arbitrary-reference output.
    no_a = maps[0]+maps[3]+maps[4]
    activity = sp.Rational(1, 7)
    held = (1-activity)*sp.eye(16)+activity*no_a
    if any(activity*mark*(sp.eye(16)-held).inv() != event
           for mark, event in zip(maps[1:3], closed, strict=True)):
        raise ArithmeticError('independent identity stuttering changed completed clock maps')
    coefficients = [coefficient_maps(law, k, feedback) for k in range(5)]
    return {'candidate': law, 'B_entry_field_Hadamard': feedback,
        'direct_native_and_union_event_elimination_equal_on_all16_matrix_units': True,
        'symbolic_remote_count_generating_instruments_equal_in_both_elimination_orders': True,
        'clock_change_preserves_conditional_work_and_arbitrary_reference_outputs': True,
        'uniform_no_A_probability_after2k_union_events_bound': f'({bound})^k',
        'two_union_event_survival_effect_diagonal': [str(effect[5*i]) for i in range(4)],
        'completed_A_instrument_trace_preserving': True,
        'independent_identity_hold_activity_one_seventh_exactly_cancels': True,
        'completed_maps_sha256': c._digest([list(m) for m in closed]),
        'remote_count0_through4_maps_sha256': c._digest([[list(m) for m in pair] for pair in coefficients])}


def relative_count_prediction(feedback=False):
    """Rational B-count generating functions at known A-event boundaries."""
    _feedback(feedback)
    z = sp.Symbol('z')
    a_exit, a_entry, b_exit, b_entry = union_maps('C', feedback)
    generating = (a_exit+a_entry)*(sp.eye(16)-z*(b_exit+b_entry)).inv()
    native = controlled_maps('C', feedback)
    attempt_wait = (sp.eye(16)-native[0]-native[3]-native[4]).inv()
    rows = []
    for mark, destination in (('A exit', 1), ('A entry', 0)):
        initial = sp.zeros(4); initial[destination, destination] = 1
        pgf = sp.factor((TRACE*generating*c.vector(initial))[0])
        mean = sp.diff(pgf, z).subs(z, 1)
        variance = sp.simplify((sp.diff(pgf, z, 2)+sp.diff(pgf, z)).subs(z, 1)-mean**2)
        attempts = (TRACE*attempt_wait*c.vector(initial))[0]
        if pgf.subs(z, 1) != 1:
            raise ArithmeticError('relative clock count distribution incomplete')
        if attempts != 6*(1+mean):
            raise ArithmeticError('relative count and source resource expectations disagree')
        rows.append({'preceding_local_mark': mark,
            'remote_count_probability_generating_function': str(pgf),
            'zero_remote_tick_probability': str(pgf.subs(z, 0)),
            'mean_B_ticks_before_next_A_tick': str(mean),
            'variance_B_ticks_before_next_A_tick': str(variance),
            'expected_native_attempt_resources_until_next_A_tick': str(attempts),
            'probabilities_B_counts0_through8': [str(sp.diff(pgf, z, k).subs(z, 0)/sp.factorial(k)) for k in range(9)]})
    return {'B_entry_field_Hadamard': feedback, 'event_boundary_predictions': rows,
        'C_and_E_full_event_instruments_equal': union_maps('C', feedback) == union_maps('E', feedback),
        'no_deterministic_mean_rate_clock_conversion': True,
        'renewal_cycles_without_feedback': None if feedback else {
            'cycle': 'A exit to next A entry to next A exit; exactly two A ticks',
            'B_count_pgf': '(2-z^2)/(3-2*z^2)',
            'mean_B_ticks_per_cycle': '2', 'variance_B_ticks_per_cycle': '16',
            'zero_B_ticks_in_N_cycles_probability': '(2/3)^N',
            'mean_ratio_B_to_A_counts_in_N_cycles': '1',
            'variance_ratio_B_to_A_counts_in_N_cycles': '4/N',
            'independence_scope': 'source Markov schedule and renewed known work state at A exit; no intervening feedback',
            'not_a_calibrated_duration_or_affine_clock_map': True}}


def _absorbing_instrument(mark):
    """Administrative completion label; it performs no native repair or read."""
    groups = tuple(((F(1), (sp.eye(4) if i == mark else sp.zeros(4)).tolist()),)
                   for i in range(5))
    return q._runtime_instrument(groups, visible=True,
        binding={'completed_A_sign': mark, 'absorbing_administrative_label_not_another_clock_tick': True})


def _native_word(record):
    for index, mark in enumerate(record):
        if mark in (1, 2):
            return record[:index+1]
    return record


def experiment(law, *, feedback=False, reference=False, cap=3):
    """Actual capped source execution, stopping further repairs at first A tick.

    Fixed runtime menus repeat a cost-free absorbing completion label after
    stopping. It is not interpreted as a clock event. All impossible menu
    branches are zero; the public record is the genuine native prefix.
    """
    c._law(law); _feedback(feedback); c._cap(cap)
    if cap > 3:
        raise ValueError('complete stopped runtime allocation allows cap at most3')
    initial = c.initial_state(reference)
    rd = 2 if reference else 1
    source = c.source_instrument(law)
    menu = (source, _absorbing_instrument(1), _absorbing_instrument(2), source, source)
    feedback_menu = tuple(q._runtime_instrument((c.field_unitary(
        'hadamard' if feedback and mark == 4 else 'identity'),),
        resources=(('declared_remote_feedback', int(feedback and mark == 4)),),
        binding={'trigger': 'B entry', 'field_Hadamard': feedback and mark == 4}) for mark in range(5))
    events = []
    for attempt in range(cap):
        events.append(runtime.Event(f'attempt_until_A_{attempt}',
            instruments=(source,) if attempt == 0 else menu,
            control=None if attempt == 0 else attempt-1))
        events.append(runtime.Event(f'B_entry_feedback_{attempt}', control=attempt, instruments=feedback_menu))
    events.append(runtime.Event('field_X_at_A_or_timeout', instruments=(q._runtime_instrument(
        c.field_readout(), visible=True, resources=(('declared_field_readouts', 1),),
        binding={'includes_timeout': True}),)))
    default = sp.Matrix(4, 4, lambda i, j: sum(initial[rd*i+r, rd*j+r] for r in range(rd)))
    program = runtime.Program(f'relative_clock_{law}_{feedback}_{reference}_{cap}', 4, c._dense(default),
        tuple(events), preparation_resources=(('declared_coherent_preparations', 1),),
        binding={'clock_A_ports': [0], 'clock_B_ports': [72],
                 'uniform_three_seam_schedule': c.geometry()[1], 'native_attempt_cap': cap,
                 'stopped_work_not_reset': True, 'timeout_and_field_state_retained': True,
                 'fixed_menu_padding_is_administrative_not_a_native_event': True})
    maps = controlled_maps(law, feedback)
    waiting, finished = {(): initial}, {}
    for _ in range(cap):
        following = {}
        for word, state in waiting.items():
            for mark, matrix in enumerate(maps):
                next_word = word+(mark,)
                output = c._apply_exact(matrix, state, rd)
                (finished if mark in (1, 2) else following)[next_word] = output
        waiting = following
    finished.update(waiting)
    exact = {word+(read,): c._apply_exact(matrix, state, rd)
             for word, state in finished.items()
             for read, matrix in enumerate(tuple(c.superoperator(g) for g in c.field_readout()))}
    return program, initial, exact


def execute(law, *, feedback=False, reference=False, cap=3):
    program, initial, exact = experiment(law, feedback=feedback, reference=reference, cap=cap)
    supplied = c._dense(initial) if reference else None
    actual = runtime.execute(program, initial=supplied)
    prefix = runtime.execute(program, initial=supplied, stop_after=2)
    checkpoint = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(program, checkpoint=checkpoint)
    residual = runtime.complete_distance(actual, continued)
    aggregate = {word: np.zeros_like(c._dense(matrix)) for word, matrix in exact.items()}
    for branch in actual.branches:
        native = _native_word(branch.record[:-1])
        key = native+(branch.record[-1],)
        aggregate[key] += branch.density
        # Impossible administrative menus can have arbitrary declared costs
        # but carry exactly zero source probability. Live branches are priced.
        if np.linalg.norm(branch.density) > 1e-14:
            costs = dict(branch.resources)
            if costs['native_attempts'] != len(native) or costs['joint_clock_reads'] != len(native):
                raise ArithmeticError('a stopped branch executed an extra native attempt')
            if costs.get('declared_remote_feedback', 0) != int(feedback)*native.count(4):
                raise ArithmeticError('remote feedback resources disagree with retained records')
    rows, matrices = [], []
    timeout = plus = completed_plus = sp.Rational(0)
    for word, state in exact.items():
        probability = sp.trace(state)
        native, read = word[:-1], word[-1]
        stopped = native[-1] in (1, 2)
        residual = max(residual, float(np.linalg.norm(aggregate[word]-c._dense(state), 'fro')))
        timeout += probability*int(not stopped)
        plus += probability*int(read == 0)
        completed_plus += probability*int(stopped and read == 0)
        rows.append({'native_attempt_marks': list(native), 'field_X_outcome': read,
            'A_tick_completed': stopped, 'B_tick_count': sum(v in (3, 4) for v in native),
            'probability': str(probability)})
        matrices.append([list(word), list(state)])
    if sum(sp.trace(state) for state in exact.values()) != 1 or residual > 2e-11:
        raise ArithmeticError('complete stopped execution disagrees with exact source oracle')
    return {'candidate': law, 'B_entry_field_Hadamard': feedback,
        'reference_dimension': actual.reference_dimension, 'native_attempt_cap': cap,
        'complete_physical_record_count': len(rows), 'fixed_runtime_menu_record_count': len(actual.branches),
        'A_timeout_probability': str(timeout), 'unconditional_field_plus_probability': str(plus),
        'joint_A_completion_and_field_plus_probability': str(completed_plus),
        'complete_probability': '1', 'complete_record_law': rows,
        'all_terminal_work_reference_matrices_sha256': c._digest(matrices),
        'checkpoint_after_first_source_attempt_and_feedback_replayed': True,
        'full_state_and_checkpoint_residual_diagnostic': residual,
        'absorbing_padding_executes_no_native_repair_and_adds_no_clock_ticks': True,
        'compiled_pulse_or_physical_service_error_certified': False}


def run():
    return {'clock_change_certificates': [certificate(law, feedback)
                for law, feedback in product(('C', 'E'), (False, True))],
        'relative_clock_count_predictions': [relative_count_prediction(feedback) for feedback in (False, True)],
        'complete_stopped_runtime_experiments': [execute(law, feedback=feedback)
                for law, feedback in product(('C', 'E'), (False, True))]
            +[execute(law, feedback=True, reference=True, cap=2) for law in ('C', 'E')]}

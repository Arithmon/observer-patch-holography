"""The existing interacting observer expressed through the shared CQ runtime.

Compilation calls the canonical preparation, native kick, complete pointer
maps and finite environment branches. The older mixed interpreter is used
only as an independent compatibility check, never as an inserted trajectory.
"""
from __future__ import annotations

import math

import numpy as np

from ..common_source import protocol
from ..common_source import open_environment as environment
from ..common_source import open_protocol as original
from . import runtime


def _noise(program, event):
    schur, metadata = environment._snapshot(program, event)
    indices = np.arange(program.work_dimension)
    numbers = (indices//4) % program.dimension if event.site == 0 else indices//(4*program.dimension)
    work_schur = schur[numbers[:, None], numbers[None, :]]
    cost = runtime.add_resources(tuple(metadata['resources_per_realization'].items()),
                                 (('independent_environment_draws', 1),))
    instrument = runtime.Instrument(schur=work_schur, resources=cost, binding={
        'kind': 'canonical_source_phase_environment', 'site': event.site,
        'variance': str(event.variance), 'order': event.order,
        'source_branch_word_hashes': [row['source']['word_sha256'] for row in metadata['branches']],
        'environment_outcomes_discarded_without_selection': True})
    return runtime.Event(f'phase_noise_after_{event.after}_site_{event.site}', instruments=(instrument,))


def compile_program(program=None, env=None):
    """Compile canonical controls; preserve supplied action and noise premises."""
    program = protocol.validate_program(protocol.continuation_program(2, 8) if program is None else program)
    env = environment.validated(program, environment.Environment(()) if env is None else env)
    preparation = protocol.execute(program, backend='source', stop_after=0)
    state = preparation.initial_state.ravel()
    initial = np.outer(state, state.conj())
    kick, kick_cost = protocol._kick(program.dimension, program.native_current, 'source')
    events = [runtime.Event('native_current_kick', instruments=(runtime.Instrument(
        kraus=((kick,),), resources=original._cost(kick_cost), binding={
            'kind': 'canonical_source_native_current_kick', 'native_current': program.native_current}),))]
    events.extend(_noise(program, event) for event in env.events if event.after == 0)
    for cursor, operation in enumerate(program.operations, 1):
        alternatives = []
        for selected in ((0, 1) if operation.kind == 'feedback' else (None,)):
            maps, cost = protocol._operation_maps(program, operation, selected, 'source')
            alternatives.append(runtime.Instrument(
                kraus=tuple((matrix,) for matrix in maps), visible=operation.kind in protocol.READS,
                resources=original._cost(cost), binding={'kind': 'canonical_source_operation', 'operation': operation.kind,
                    'model_coordinate_duration': str(operation.time), 'selected_record_value': selected,
                    'enabled': operation.enabled}))
        events.append(runtime.Event(f'{cursor}_{operation.kind}', instruments=tuple(alternatives),
                                    control=operation.record))
        events.extend(_noise(program, event) for event in env.events if event.after == cursor)
    return runtime.Program('common_interacting_source', program.work_dimension, initial, tuple(events),
                           preparation_resources=original._cost(preparation.preparation_resources), binding={
        'canonical_program': program.specification(), 'canonical_environment': env.specification(),
        'prediction_backend': 'contracted_native_source_words',
        'time_parameter': 'supplied dimensionless model coordinate',
        'source_selection_or_physical_service_derivation': False})


def execute(program=None, env=None, **kwargs):
    return runtime.execute(compile_program(program, env), **kwargs)


def compatibility_report(program=None, env=None, *, initial=None):
    """Compare complete work/reference outputs, not just visible probabilities."""
    program = protocol.validate_program(protocol.continuation_program(2, 8) if program is None else program)
    env = environment.validated(program, environment.Environment(()) if env is None else env)
    actual = execute(program, env, initial=initial)
    expected = original.execute(program, env, initial=initial)
    if actual.records != expected.records or actual.reference_dimension != expected.reference_dimension:
        raise ArithmeticError('the shared runtime changed complete canonical output labels')
    differences = []
    for left, right in zip(actual.branches, expected.branches, strict=True):
        difference = left.density-right.density
        differences.append(float(np.abs(np.linalg.eigvalsh((difference+difference.conj().T)/2)).sum())/2)
        if dict(left.resources) != dict(right.resources):
            raise ArithmeticError('the shared runtime changed canonical per-path source resources')
    distance = math.fsum(differences)
    if not math.isfinite(distance) or distance > runtime.TOLERANCE:
        raise ArithmeticError('the shared runtime changed a complete canonical quantum output')
    if dict(actual.expanded_tree_resources) != dict(expected.expanded_tree_resources):
        raise ArithmeticError('the shared runtime changed canonical expanded-tree source resources')
    return {'complete_cq_trace_distance': distance, 'visible_records': len(actual.branches),
            'reference_dimension': actual.reference_dimension,
            'probability_mass': actual.probability_mass,
            'source_resources_agree': True,
            'all_zero_probability_outcomes_retained': True,
            'canonical_program_digest': program.digest,
            'comparison': 'independent_existing_mixed_interpreter_on_the_same_canonical_maps'}

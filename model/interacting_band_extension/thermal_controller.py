"""Audit the existing homogeneous controller on the cold work code.

This is an audit of ProductControllerInstrument, not a second controller.
The caller supplies its already executed preparation and numerical target.
Every native tape and timeout prefix remains part of the complete output.
"""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np

from ..occupation_control import intervals as iv
from . import construction as source
from .thermal_source import _amplitude, _source_budget, pure_record_distance


def audit(raw, target, *, occupation=8, sign=1, interval=F(1, 16), steps=64):
    """All-input comparison, then the retained actual work/reference state.

    The complete controller is enumerated once. Identical work operators
    may be reused in reductions only after exact array equality is checked
    on every actual tape. This does not merge or postselect record labels.
    """
    program = source.ProductControllerInstrument(occupation=occupation, sign=sign,
                                                   interval=interval, steps=steps)
    dimension = program.lattice.dimension
    for value in (raw, target):
        if _amplitude(value).shape != (dimension, dimension):
            raise ValueError('one complete work/reference preparation matrix is required')
    preparation_error = F(iv.frobenius_upper(*iv.complex_sub(
        iv.complex_exact(raw), iv.complex_exact(target))))
    if preparation_error > F(1, 10**9):
        raise ArithmeticError('the supplied source preparation exceeds its numerical guard')
    identity = np.eye(dimension)
    direct = source._product_clock_prefixes(program.lattice, identity, source=True,
                         interval=program.interval, steps=program.steps, sign=sign)
    ideal = source._product_clock_prefixes(program.lattice, target, source=False,
                         interval=program.interval, steps=program.steps, sign=sign)
    expected_tapes = {(first, second): final for first, second, final in source.product_incident_tapes()}
    expected_labels = {(tape, bits) for tape in expected_tapes
                       for bits in product((0, 1), repeat=sum(record[3] for record in tape))}
    seen_labels = set()
    groups, tape_counts, lengths = {}, {}, {}
    maximum_error = 0.
    completed = program.completed()
    for state, operation, weight, count in completed:
        source.product_controller_transition(state, 'idle', steps=steps)
        _amplitude(operation)
        if (type(weight) is not F or weight != F(1, 3600)
                or state.instruction != 'halt' or state.tape not in expected_tapes
                or state.loads != expected_tapes[state.tape]):
            raise ArithmeticError('the complete native controller alphabet or weight changed')
        # One initial kick, two native attempts and one halt transition;
        # each incident contributes its palindrome and a read. Only the
        # first read can append the positive feedback instruction.
        expected_count = (4+len(state.bits)*(2*steps*len(program.lattice.terms)+1)
                          +(state.bits[0] if state.bits else 0))
        if type(count) is not int or count != expected_count:
            raise ArithmeticError('the complete controller instruction budget changed')
        label = state.tape, state.bits
        if label in seen_labels:
            raise ArithmeticError('a complete native/field output label was duplicated')
        seen_labels.add(label)
        key = len(state.bits), state.bits
        residual = float(np.linalg.norm(operation-direct[key[0]][key[1]], 'fro'))
        if not math.isfinite(residual):
            raise ArithmeticError('the all-input controller residual is nonfinite')
        maximum_error = max(maximum_error, residual)
        if key in groups:
            previous, mass = groups[key]
            if not np.array_equal(previous, operation):
                raise ArithmeticError('native-tape grouping changed a work operator')
            groups[key] = previous, mass+weight
        else:
            groups[key] = operation, weight
        tape_counts[state.tape] = tape_counts.get(state.tape, 0)+1
        lengths[count] = lengths.get(count, 0)+1
    if (len(completed) != 4900 or seen_labels != expected_labels
            or set(tape_counts) != set(expected_tapes)
            or set(groups) != {(tick, bits) for tick, rows in enumerate(direct) for bits in rows}
            or maximum_error > 2e-12):
        raise ArithmeticError('the complete all-input controller changed its independent history')
    gram = sum((float(weight)*(operation.conj().T@operation)
                for operation, weight in groups.values()), np.zeros((dimension, dimension), complex))
    completeness = float(np.linalg.norm(gram-identity, 2))
    if not math.isfinite(completeness) or completeness > 2e-11:
        raise ArithmeticError('the complete all-input controller is not normalized within its guard')
    rows, complete_error = [], 0.
    for (tick, bits), (operation, weight) in sorted(groups.items()):
        actual = operation@raw
        _amplitude(actual)
        probability = float(weight)*float(np.linalg.norm(actual)**2)
        error = pure_record_distance([actual], [ideal[tick][bits]])
        if not math.isfinite(probability) or probability < 0:
            raise ArithmeticError('the retained source-prefix probability is nonfinite')
        complete_error += float(weight)*error
        rows.append({'incident_ticks': tick, 'timeout': tick < 2, 'field_bits': list(bits),
                     'native_tape_weight': weight, 'probability': probability})
    if abs(sum(row['probability'] for row in rows)-1) > 1e-8:
        raise ArithmeticError('the retained source-prefix output lost probability')
    split = source.splitting_bound(program.lattice, program.interval, program.steps)
    # Charge the declared, checked preparation guard, rather than serializing
    # a platform-dependent measured norm as an exact rational requirement.
    budget = _source_budget(target, F(1, 10**9), 4*steps*len(program.lattice.terms)+4,
                            split['two_interval_complete_output_trace_distance_upper'])
    if (not math.isfinite(complete_error) or complete_error
            > float(budget['target_scaled_splitting_error_upper']+budget['source_numerical_guard'])):
        raise ArithmeticError('the retained controller output exceeds its conditional finite budget')
    # Idle and the existing full six-level service exports preserve the
    # committed native prefix; failure payloads have their own dimensions.
    state = source.product_controller_transition(source.ProductControllerState(), steps=steps)
    edge = next(i for i, seam in enumerate(source.native_seams()) if 10 in seam)
    written = source.product_controller_transition(state, (edge, 0), steps=steps)
    program._checked_branches(written, ((source.product_controller_transition(
        written, 'idle', steps=steps), None, F(1)),))
    for code_dimension in (2, 4):
        for failed, _, label in source.product_controller_export(written, code_dimension, steps=steps):
            if label != 'code':
                program._checked_branches(written, ((failed, None, F(1)),))
                if program.branches(failed) != ((failed, None, F(1)),):
                    raise ArithmeticError('a committed failure output is not absorbing')
    return {'occupation': occupation, 'initial_kick_sign': sign,
        'initial_kick_sign_is_fixed_program_setting': True,
        'feedback_kick': F(1, 4), 'retained_work_input_columns': dimension,
        'retained_preparation_reference_columns': dimension,
        'complete_native_tapes': len(tape_counts), 'complete_native_field_outputs': len(completed),
        'maximum_active_instructions': max(lengths),
        'every_terminal_instruction_count_independently_checked': True,
        'all_input_source_Kraus_comparison_residual_numerical': maximum_error,
        'all_input_channel_completeness_residual_numerical': completeness,
        'complete_prepared_work_reference_prefix_distance_numerical': complete_error,
        'complete_clock_field_prefix_law': rows,
        'timeout_probability': F(35, 36),
        'postcommit_prefix_transition_classes_checked': len(program.prefix_transition_classes),
        'postcommit_public_projector_fixed_point_argument_dimension_independent': True,
        'native_transitions_and_prefixes_changed': False,
        'source_preparation_executed_again': False,
        'new_d8_quantum_access_algebra_rank_claimed': False,
        'all_terminal_native_timeout_work_reference_outputs_retained': True,
        'full_six_level_export_complement_labels_retained': True,
        'valid_code_service_failure_probability_exact': F(0),
        'nonzero_service_failure_on_actual_preparation_executed': False,
        'arbitrary_unfinished_mid_service_output_bound_claimed': False,
        'physical_controller_or_service_fidelity_certified': False,
        'finite_source_budget_under_numerical_guards': budget}

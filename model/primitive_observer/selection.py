"""Attach certain-state instrument selection to the compiled spin pointer.

The existing exact theorem selects a projective outcome map once its effect
and invariance on every outcome-certain state are required. Here every matrix
unit of the actual native CNOT pointer is checked against that exact map.
These numerical compiler checks do not themselves produce Lean certificates
or remove the supplied pointer basis, blank, commitment or physical clock.
"""
from __future__ import annotations

import math

import numpy as np

from . import anharmonic, runtime


DIMENSION = 4
TOLERANCE = 2e-10
THEOREM_PATH = 'Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean'
THEOREM = 'EventAlgebra.instrument_component_fixing_certain_states_is_lueders'


def _metric(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ArithmeticError(f'{name} must be a finite nonnegative diagnostic')
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ArithmeticError(f'{name} must be a finite nonnegative diagnostic')
    return value


def _map_residual(matrices, targets):
    """Compare all linear-map basis inputs, including off-diagonal matrices."""
    size = targets[0].shape[0]
    maximum = 0.
    for matrix, target in zip(matrices, targets, strict=True):
        for i in range(size):
            for j in range(size):
                # K |i><j| K* is the outer product of columns i and j.
                actual = np.outer(matrix[:, i], matrix[:, j].conj())
                expected = np.outer(target[:, i], target[:, j].conj())
                maximum = max(maximum, _metric(np.linalg.norm(actual-expected, 'fro'),
                                               'matrix-unit outcome map'))
    return maximum


def _instrument(matrices, metadata):
    return runtime.Instrument(
        kraus=tuple((matrix,) for matrix in matrices), visible=True,
        resources=anharmonic._resources(metadata, fresh_pointer_blanks=1, committed_classical_bits=1),
        binding={'kind': 'native_CNOT_spin_pointer', 'source': metadata,
                 'basis_blank_and_commitment_supplied': True})


def _reference_report(matrices, targets, metadata):
    size, reference = 2*DIMENSION, 3
    state = np.zeros((size, reference), complex)
    state[0, 0], state[1, 1], state[DIMENSION, 2], state[-1, 0] = 1., 2j, 3-1j, .5
    state /= np.linalg.norm(state)
    density = np.outer(state.ravel(), state.ravel().conj())
    event = runtime.Event('complete_spin_pointer', instruments=(_instrument(matrices, metadata),))
    program = runtime.Program('native_projective_pointer_selection', size, np.eye(size)/size,
                              (event,), binding={'reference_is_untouched': True,
                                                 'time_parameter': 'one source instrument event'})
    history = runtime.execute(program, initial=density)
    ideal = runtime.Program('exact_projector_comparison', size, np.eye(size)/size,
                            (runtime.Event('complete_spin_pointer', instruments=(
                                runtime.Instrument(kraus=tuple((p,) for p in targets), visible=True),)),))
    expected = runtime.execute(ideal, initial=density)
    phases = (np.exp(.37j), np.exp(-1.1j))
    changed = tuple(phase*matrix for phase, matrix in zip(phases, matrices, strict=True))
    phased_program = runtime.Program('committed_pointer_phase_control', size, np.eye(size)/size,
                                    (runtime.Event('complete_spin_pointer', instruments=(
                                        _instrument(changed, metadata),)),))
    phased = runtime.execute(phased_program, initial=density)
    marginal = sum(np.trace(branch.density.reshape(size, reference, size, reference),
                            axis1=0, axis2=2) for branch in history.branches)
    if history.records != ((0,), (1,)) or history.reference_dimension != reference:
        raise ArithmeticError('complete pointer/reference output was reduced')
    return {
        'reference_dimension': reference,
        'quantum_dimension_per_outcome': size*reference,
        'complete_outcome_probabilities': [branch.probability for branch in history.branches],
        'projector_complete_cq_distance': _metric(runtime.complete_distance(history, expected), 'reference map'),
        'outcome_phase_complete_cq_distance': _metric(runtime.complete_distance(history, phased), 'outcome phase'),
        'outcome_phase_matrix_unit_residual': _map_residual(changed, matrices),
        'reference_marginal_residual': _metric(np.linalg.norm(marginal-state.T@state.conj(), 'fro'),
                                               'reference marginal'),
        'both_outcomes_and_disturbed_quantum_states_retained': True,
        'matrix_unit_identity_extends_algebraically_to_every_finite_reference': True,
        'numerical_checks_are_exact_compiler_or_diamond_error_certificates': False,
        'uncommitted_coherent_pointer_interference_is_covered': False,
    }


def run():
    """Produce one source-bound conditional selection and its countercontrol."""
    matrices, metadata = anharmonic._spin_record(DIMENSION)
    size = 2*DIMENSION
    targets = tuple(np.diag([float(i//DIMENSION == spin) for i in range(size)]) for spin in (0, 1))
    # Runtime validation checks that the two actual output maps form one
    # complete CP instrument before the theorem attachment is evaluated.
    _instrument(matrices, metadata)
    effect_error = _metric(max(np.linalg.norm(matrix.conj().T@matrix-target, 'fro')
                              for matrix, target in zip(matrices, targets, strict=True)), 'projective effect')
    map_error = _map_residual(matrices, targets)
    certain_error = 0.
    certain_units = 0
    for spin, matrix in enumerate(matrices):
        for i in range(spin*DIMENSION, (spin+1)*DIMENSION):
            for j in range(spin*DIMENSION, (spin+1)*DIMENSION):
                unit = np.zeros((size, size), complex)
                unit[i, j] = 1.
                certain_error = max(certain_error, _metric(
                    np.linalg.norm(matrix@unit@matrix.conj().T-unit, 'fro'), 'certain-state invariance'))
                certain_units += 1

    # An intra-outcome unitary preserves every outcome effect and repeatable
    # spin label, but changes a state whose outcome is already certain.
    twist = np.eye(size, dtype=complex)
    twist[[0, 1]] = twist[[1, 0]]
    twisted = tuple(twist@matrix for matrix in matrices)
    _instrument(twisted, metadata)
    negative_effect = _metric(max(np.linalg.norm(a.conj().T@a-b.conj().T@b, 'fro')
                                  for a, b in zip(twisted, matrices, strict=True)), 'twisted effect')
    certain = np.zeros((size, size), complex)
    certain[0, 0] = 1.
    changed = twisted[0]@certain@twisted[0].conj().T
    delta = changed-certain
    disturbance = _metric(np.abs(np.linalg.eigvalsh((delta+delta.conj().T)/2)).sum()/2,
                          'twisted certain-state disturbance')
    repeatability = _metric(np.linalg.norm(targets[0]@changed@targets[0]-changed, 'fro'),
                            'twisted spin certainty')
    reference = _reference_report(matrices, targets, metadata)
    checks = {
        'actual_source_effects_are_declared_projectors': effect_error < TOLERANCE,
        'all_matrix_units_match_selected_outcome_maps': map_error < TOLERANCE,
        'whole_certain_state_corners_are_fixed': certain_error < TOLERANCE,
        'complete_reference_outputs_match': reference['projector_complete_cq_distance'] < TOLERANCE,
        'committed_outcome_phases_cancel': reference['outcome_phase_complete_cq_distance'] < TOLERANCE
                                         and reference['outcome_phase_matrix_unit_residual'] < TOLERANCE,
        'spectator_reference_marginal_preserved': reference['reference_marginal_residual'] < TOLERANCE,
        'same_effect_control_changes_a_certain_state': negative_effect < TOLERANCE and disturbance > .99,
        'repeatability_alone_does_not_select_the_map': repeatability < TOLERANCE and disturbance > .99,
    }
    if not all(checks.values()):
        raise ArithmeticError('compiled spin-pointer selection attachment failed')
    return {
        'schema': 'oph.primitive_observer.instrument_selection.v1',
        'scope': 'conditional selection of the native projective spin-pointer outcome map',
        'theorem_path': THEOREM_PATH, 'theorem': THEOREM,
        'theorem_conditions': ['finite Kraus-form outcome', 'specified projective effect',
                               'invariance of every state certain of that outcome'],
        'source': metadata,
        'occupation_dimension': DIMENSION, 'spin_work_dimension': size,
        'outcome_count': 2, 'outcome_matrix_units_checked': len(matrices)*size**2,
        'certain_corner_matrix_units_checked': certain_units,
        'projective_effect_residual': effect_error,
        'outcome_map_matrix_unit_residual': map_error,
        'certain_corner_identity_residual': certain_error,
        'reference_and_phase_checks': reference,
        'same_effect_negative_control': {
            'intervention': 'swap occupation levels zero and one inside spin-zero outcome',
            'effect_residual': negative_effect,
            'certain_state_trace_distance': disturbance,
            'normalized_output_spin_certainty_residual': repeatability,
            'outcome_zero_probability_on_witness': float(np.trace(changed).real),
        },
        'independent_postmeasurement_map_choice_removed_under_theorem_conditions': True,
        'pointer_basis_blank_commitment_or_physical_clock_selected': False,
        'complete_spin_plus_smooth_four_outcome_detector_selected': False,
        'all_OPH_observer_instruments_selected': False,
        'new_experimental_data_consumed': False,
        'compiled_matrices_formally_certified': False,
        'checks': checks, 'all_pass': all(checks.values()),
    }

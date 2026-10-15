"""A real source pointer satisfies the selecting condition; effects alone do not."""
import json

import numpy as np
import pytest

from . import anharmonic, selection


@pytest.fixture(scope='module')
def result():
    return selection.run()


def test_source_pointer_attaches_every_matrix_unit_to_the_existing_selection_theorem(result):
    assert result['all_pass'] and all(result['checks'].values())
    assert result['outcome_count'] == 2
    assert result['outcome_matrix_units_checked'] == 128
    assert result['certain_corner_matrix_units_checked'] == 32
    assert result['projective_effect_residual'] < 2e-12
    assert result['outcome_map_matrix_unit_residual'] < 2e-12
    assert result['certain_corner_identity_residual'] < 2e-12
    assert result['theorem'] == 'EventAlgebra.instrument_component_fixing_certain_states_is_lueders'
    json.dumps(result, allow_nan=False)


def test_complete_reference_and_committed_outcome_phase_laws_are_preserved(result):
    reference = result['reference_and_phase_checks']
    assert reference['reference_dimension'] == 3
    assert reference['quantum_dimension_per_outcome'] == 24
    assert sum(reference['complete_outcome_probabilities']) == pytest.approx(1., abs=2e-12)
    assert all(0 < mass < 1 for mass in reference['complete_outcome_probabilities'])
    assert reference['projector_complete_cq_distance'] < 2e-12
    assert reference['outcome_phase_complete_cq_distance'] < 2e-12
    assert reference['outcome_phase_matrix_unit_residual'] < 2e-12
    assert reference['reference_marginal_residual'] < 2e-12
    assert not reference['uncommitted_coherent_pointer_interference_is_covered']


def test_same_effect_repeatable_control_fails_the_selecting_condition(result):
    control = result['same_effect_negative_control']
    assert control['effect_residual'] < 2e-12
    assert control['certain_state_trace_distance'] == pytest.approx(1., abs=2e-12)
    assert control['normalized_output_spin_certainty_residual'] < 2e-12
    assert control['outcome_zero_probability_on_witness'] == pytest.approx(1., abs=2e-12)


def test_the_attachment_does_not_select_unprovided_basis_preparation_or_smooth_instrument(result):
    assert result['independent_postmeasurement_map_choice_removed_under_theorem_conditions']
    assert not result['pointer_basis_blank_commitment_or_physical_clock_selected']
    assert not result['complete_spin_plus_smooth_four_outcome_detector_selected']
    assert not result['all_OPH_observer_instruments_selected']
    assert not result['new_experimental_data_consumed']
    assert not result['compiled_matrices_formally_certified']


def test_a_source_change_with_identical_effects_cannot_pass_the_attachment(monkeypatch):
    original = anharmonic._spin_record
    def altered(dimension):
        matrices, metadata = original(dimension)
        twist = np.eye(2*dimension, dtype=complex)
        twist[[0, 1]] = twist[[1, 0]]
        return tuple(twist@matrix for matrix in matrices), metadata
    monkeypatch.setattr(anharmonic, '_spin_record', altered)
    with pytest.raises(ArithmeticError, match='selection attachment failed'):
        selection.run()


def test_committed_kraus_phase_is_not_a_new_selectable_parameter(monkeypatch):
    original = anharmonic._spin_record
    def phased(dimension):
        matrices, metadata = original(dimension)
        return tuple(np.exp(1j*phase)*matrix for phase, matrix in zip((.2, -.7), matrices, strict=True)), metadata
    monkeypatch.setattr(anharmonic, '_spin_record', phased)
    assert selection.run()['all_pass']


def test_nonfinite_diagnostics_fail_closed():
    for value in (float('nan'), float('inf'), -1., True):
        with pytest.raises(ArithmeticError, match='finite nonnegative'):
            selection._metric(value, 'invalid')

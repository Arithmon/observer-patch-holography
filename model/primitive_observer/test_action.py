"""Exact native response attachment, ambiguity and failed-hypothesis controls."""
from fractions import Fraction as F
import json

import pytest

from . import action


def test_actual_native_response_recovers_relative_action_and_preserves_scope():
    result = action.run()
    assert result['all_pass']
    assert result['response_recovery']['full_native_inputs_checked'] == 78
    relative = result['relative_action']
    assert relative['normalized_capacities'] == ['1/12']*12
    assert len(relative['seam_conductances']) == 30
    assert {row[2] for row in relative['seam_conductances']} == {'1/720'}
    assert relative['recovered_seam_rates'] == ['1/30']*30
    assert relative['unit_port_reading_energy'] == '1/288'
    assert not relative['absolute_action_scale_selected']
    assert not result['history_register_boundary']['vacuum_covariance_theorem_instantiated']
    assert not result['scope']['physical_field_action_or_interaction_coefficients_selected']
    assert result['diagonal_quantum_extension']['recorded_branch_addresses_checked'] == 9360
    assert result['diagonal_quantum_extension']['mean_response_extends_to_all_quantum_and_reference_inputs']
    spectrum = result['recovered_response_spectrum']
    assert spectrum['slow_response_dimension'] == 3
    assert spectrum['characteristic_polynomial_coefficients'][:3] == [1, -60, 1620]
    assert not spectrum['physical_spatial_dimension_selected']
    json.dumps(result, allow_nan=False)


def test_reciprocity_solver_handles_nonuniform_capacities_without_using_native_formula():
    # A three-node path with capacities proportional to 1,2,3 and conductances 2,5.
    response = ((F(2), F(-2), F(0)),
                (F(-1), F(7, 2), F(-5, 2)),
                (F(0), F(-5, 3), F(5, 3)))
    result = action.relative_action(response)
    assert result['capacities'] == (F(1, 6), F(1, 3), F(1, 2))
    assert result['seams'] == ((0, 1, F(1, 3)), (1, 2, F(5, 6)))
    assert action.energy((2, -1, 3), result) == F(49, 6)
    assert action._weighted_repair_check(response, result) == [F(3), F(25, 6)]


def test_response_scale_changes_energy_rates_but_not_relative_capacity():
    response = ((1, -1), (-2, 2))
    baseline = action.relative_action(response)
    scaled = action.relative_action(tuple(tuple(7*x for x in row) for row in response))
    assert baseline['capacities'] == scaled['capacities'] == (F(2, 3), F(1, 3))
    assert action.energy((1, 0), scaled) == 7*action.energy((1, 0), baseline)


@pytest.mark.parametrize('matrix', [
    ((0, 0), (0, 0)),  # disconnected, no unique normalized capacity
    ((1, -1), (0, 0)),  # one-way support
    ((-1, 1), (1, -1)),  # anti-diffusive signs
    ((2, -1), (-1, 1)),  # grounding violates the conservative declared sector
    ((2, -1, -1), (-2, 3, -1), (-1, -1, 2)),  # inconsistent cycle ratios
    ((True, -1), (-1, 1)), ((1., -1.), (-1., 1.)), (), ((1, 2),),
])
def test_invalid_or_nonidentifiable_relative_action_fails_closed(matrix):
    with pytest.raises(ValueError):
        action.relative_action(matrix)


def test_native_rule_mutation_cannot_hide_behind_supplied_laplacian(monkeypatch):
    original = action.native.native
    def always_first(state, seam, coin):
        return original(state, seam, 0)
    monkeypatch.setattr(action.native, 'native', always_first)
    with pytest.raises(ArithmeticError, match='not linear'):
        action.recover_response()


def test_quantum_extension_rejects_record_erasure_that_merges_input_amplitudes(monkeypatch):
    original = action.native.native
    seams = tuple(action.slow_geometry.source().seams())
    def erased_current(state, seam, coin):
        output, _ = original(state, seam, coin)
        return output, (coin, 0)
    monkeypatch.setattr(action.native, 'native', erased_current)
    with pytest.raises(ArithmeticError, match='not injective'):
        action._diagonal_quantum_extension(seams)


def test_actual_native_noise_violates_Gaussian_fourth_moment():
    result = action._noise_control()
    assert result['variance'] == '1/4'
    assert result['centered_fourth_moment'] == '1/16'
    assert result['Gaussian_fourth_moment_at_same_variance'] == '3/16'


def test_same_scheduler_mean_action_does_not_select_quantum_instrument():
    seams = tuple(action.slow_geometry.source().seams())
    result = action._coherent_control(seams)
    laws = result['complete_binary_probabilities']
    assert laws['C'] != laws['E']
    assert all(sum(map(F, law)) == 1 for law in laws.values())
    assert result['same_classical_population_transition']


def test_native_stationary_covariance_fails_ground_relation_after_removing_constant_mode():
    _, response, seams = action.recover_response()
    result = action._stationary_vacuum_control(action.relative_action(response), seams)
    assert result['stationarity_branches_checked'] == 3960
    assert result['covariance_diagonal'] == '5/36'
    assert result['covariance_off_diagonal'] == '-5/396'
    assert result['hbar_squared_required_by_adjacent_entry'] == '5/39204'
    assert result['hbar_squared_required_by_nonadjacent_entry'] == '0'
    assert result['conserved_constant_mode_removed']
    assert not result['one_common_hbar_matches_ground_relation']
    assert not result['excludes_other_observer_fields_or_states']


def test_energy_requires_exact_readings_of_the_declared_dimension():
    attachment = action.relative_action(((1, -1), (-1, 1)))
    for bad in ((1,), (1, 0, 1), (True, 0), (.5, 0)):
        with pytest.raises(ValueError):
            action.energy(bad, attachment)

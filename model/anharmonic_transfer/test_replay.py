"""Exact local seals, independently validated certificates, portable floats."""
from copy import deepcopy
from fractions import Fraction
import math

import pytest

from . import bounds, experiment, replay
from .test_experiment import build_prediction


@pytest.fixture(scope='module')
def receipt():
    with pytest.MonkeyPatch.context() as monkeypatch:
        prediction, _ = build_prediction(monkeypatch)
    return experiment.evaluate(prediction)


def reseal(receipt):
    receipt[replay.SEAL] = experiment.digest(receipt['prediction'])


def _law_row():
    return {'acceptance_probability': .75, 'rejection_probability': .25,
            'predicted_detector_probabilities': [.6, .4],
            'spin_detector_joint_probabilities': [[.45, .3], [.125, .125]],
            'predicted_spin_mean': .2}


@pytest.fixture(scope='module')
def zero_rejection_receipt(receipt):
    other = deepcopy(receipt)
    row = other['prediction']['runs'][0]['zero_time_predictions'][0]
    row.update(acceptance_probability=1., rejection_probability=0.,
               predicted_detector_probabilities=[.5, .5],
               spin_detector_joint_probabilities=[[.5, .5], [0., 0.]],
               predicted_spin_mean=0., predicted_position_estimator=0.)
    budget = row['ideal_finite_model_error']
    budget['selected_readout'] = bounds.conditional_readout_bound(
        Fraction(budget['exact_complete_trace_distance_upper']), 1.)
    return experiment.evaluate(other['prediction'])


@pytest.mark.parametrize('side', ['expected', 'actual', 'both'])
def test_resealed_negative_zero_outcome_cannot_hide_in_portable_tolerance(
        zero_rejection_receipt, side):
    other = deepcopy(zero_rejection_receipt)
    row = other['prediction']['runs'][0]['zero_time_predictions'][0]
    row['spin_detector_joint_probabilities'][1][0] = -1e-10
    reseal(other)
    expected = other if side in ('expected', 'both') else zero_rejection_receipt
    actual = other if side in ('actual', 'both') else zero_rejection_receipt
    with pytest.raises(ValueError, match='raw boundary roundoff'):
        replay.compare_receipt(expected, actual)


def test_literal_tiny_signed_roundoff_and_upper_boundary_overshoot_are_preserved(
        zero_rejection_receipt):
    other = deepcopy(zero_rejection_receipt)
    row = other['prediction']['runs'][0]['zero_time_predictions'][0]
    # Magnitudes present in the retained canonical source execution.
    row['spin_detector_joint_probabilities'][1][0] = -4.216224026370642e-17
    row['acceptance_probability'] = 1.0000000000000029
    reseal(other)
    before = deepcopy(other)
    replay.compare_receipt(other, zero_rejection_receipt)
    replay.compare_receipt(zero_rejection_receipt, other)
    assert other == before
    assert row['spin_detector_joint_probabilities'][1][0] < 0
    assert row['acceptance_probability'] > 1


@pytest.mark.parametrize('field', ['acceptance_probability', 'rejection_probability',
                                  'selected', 'joint'])
@pytest.mark.parametrize('bad', [True, 0, Fraction(1, 2), float('nan'), float('inf'),
                               -1e-10, 1+1e-10])
def test_probability_values_require_original_finite_floats_in_the_raw_domain(field, bad):
    row = _law_row()
    if field == 'selected':
        row['predicted_detector_probabilities'][0] = bad
    elif field == 'joint':
        row['spin_detector_joint_probabilities'][1][0] = bad
    else:
        row[field] = bad
    with pytest.raises(ValueError, match='Finite original probabilities'):
        replay._validate_probability_law(row)


@pytest.mark.parametrize('mutation', ['selected_tuple', 'missing_selected', 'joint_tuple',
                                     'missing_spin_row', 'missing_detector', 'extra_detector'])
def test_probability_law_requires_every_spin_and_detector_outcome(mutation):
    row = _law_row()
    if mutation == 'selected_tuple':
        row['predicted_detector_probabilities'] = (.6, .4)
    elif mutation == 'missing_selected':
        del row['predicted_detector_probabilities']
    elif mutation == 'joint_tuple':
        row['spin_detector_joint_probabilities'] = ((.45, .3), (.125, .125))
    elif mutation == 'missing_spin_row':
        row['spin_detector_joint_probabilities'].pop()
    elif mutation == 'missing_detector':
        row['spin_detector_joint_probabilities'][1].pop()
    else:
        row['predicted_detector_probabilities'].append(0.)
    with pytest.raises(ValueError, match='complete two-spin by two-detector'):
        replay._validate_probability_law(row)


@pytest.mark.parametrize('mutation', ['aggregate_mass', 'selected_mass', 'spin_marginal',
                                     'conditional', 'mean'])
def test_many_portable_float_changes_cannot_forge_a_complete_distribution(mutation):
    original = _law_row()
    changed = deepcopy(original)
    if mutation == 'aggregate_mass':
        for row in changed['spin_detector_joint_probabilities']:
            row[0] += 6e-11
            row[1] += 6e-11
    elif mutation == 'selected_mass':
        changed['predicted_detector_probabilities'][0] += 1e-11
    elif mutation == 'spin_marginal':
        changed['rejection_probability'] += 1e-11
    elif mutation == 'conditional':
        changed['spin_detector_joint_probabilities'][0][0] += 1e-11
        changed['spin_detector_joint_probabilities'][0][1] -= 1e-11
    else:
        changed['predicted_spin_mean'] += 1e-11
    # Every individual change is allowed by portable closeness alone.
    replay._compare(original, changed)
    with pytest.raises(ValueError, match='normalization or marginal consistency'):
        replay._validate_probability_law(changed)


def test_consistent_accumulated_source_mass_roundoff_is_not_renormalized():
    row = _law_row()
    factor = 1+3e-11
    row['acceptance_probability'] *= factor
    row['rejection_probability'] *= factor
    row['spin_detector_joint_probabilities'] = [
        [value*factor for value in part] for part in row['spin_detector_joint_probabilities']]
    before = deepcopy(row)
    replay._validate_probability_law(row)
    assert row == before
    assert row['acceptance_probability']+row['rejection_probability'] > 1


def test_identical_receipt_replays_without_mutating_inputs(receipt):
    before = deepcopy(receipt)
    replay.compare_receipt(receipt, deepcopy(receipt))
    assert receipt == before


def test_portable_acceptance_revalidates_derived_exact_rationals(receipt):
    other = deepcopy(receipt)
    row = other['prediction']['runs'][0]['predictions'][0]
    row['acceptance_probability'] += 1e-13
    budget = row['ideal_finite_model_error']
    budget['selected_readout'] = bounds.conditional_readout_bound(
        Fraction(budget['exact_complete_trace_distance_upper']), row['acceptance_probability'])
    reseal(other)
    assert other[replay.SEAL] != receipt[replay.SEAL]
    assert budget['selected_readout']['exact_conditional_detector_tv_upper'] != receipt[
        'prediction']['runs'][0]['predictions'][0]['ideal_finite_model_error'][
            'selected_readout']['exact_conditional_detector_tv_upper']
    replay.compare_receipt(receipt, other)


def test_numerical_prediction_roundoff_can_change_local_seal(receipt):
    other = deepcopy(receipt)
    row = other['prediction']['runs'][0]['predictions'][0]
    _shift_mean(row, 1e-13)
    other = experiment.evaluate(other['prediction'])
    replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('change, accepted', [
    (1.9e-10, True), (2.1e-10, False), (5e-9, False), (5e-8, False)])
def test_portable_tolerance_remains_bounded_after_source_arithmetic_repair(receipt, change, accepted):
    other = deepcopy(receipt)
    _shift_mean(other['prediction']['runs'][0]['predictions'][0], change)
    other = experiment.evaluate(other['prediction'])
    if accepted:
        replay.compare_receipt(receipt, other)
    else:
        with pytest.raises(ValueError, match='Floating replay mismatch'):
            replay.compare_receipt(receipt, other)


def test_refinement_summary_propagates_individual_prediction_tolerance(receipt):
    other = deepcopy(receipt)
    for run, change in ((other['prediction']['runs'][1], -9e-11),
                        (other['prediction']['runs'][2], 9e-11)):
        _shift_mean(run['predictions'][0], change)
    other = experiment.evaluate(other['prediction'])
    assert other['numerical_refinements']['occupation']['maximum_position_estimator_change'] > 2e-10
    replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('field', ['maximum_binary_probability_change', 'maximum_position_estimator_change'])
def test_refinement_summary_cannot_be_changed_independently(receipt, field):
    other = deepcopy(receipt)
    other['numerical_refinements']['occupation'][field] += 1e-13
    with pytest.raises(ValueError, match='Derived receipt.*numerical_refinements'):
        replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('mutation', ['stale_seal', 'conditional_fraction', 'unconditional_fraction',
                                    'word_digest', 'integer', 'large_float', 'type', 'nonfinite'])
def test_replay_rejects_integrity_and_scientific_contract_changes(receipt, mutation):
    other = deepcopy(receipt)
    row = other['prediction']['runs'][0]['predictions'][0]
    if mutation == 'stale_seal':
        row['predicted_spin_mean'] += 1e-13
    elif mutation == 'conditional_fraction':
        row['ideal_finite_model_error']['selected_readout']['exact_conditional_detector_tv_upper'] = '1/2'
        reseal(other)
    elif mutation == 'unconditional_fraction':
        row['ideal_finite_model_error']['splitting']['exact_half_diamond_error_upper'] = '1/2'
        reseal(other)
    elif mutation == 'word_digest':
        key = next(iter(other['prediction']['source_bindings']))
        other['prediction']['source_bindings'][key] = '0'*64
        reseal(other)
    elif mutation == 'integer':
        row['steps'] += 1
        reseal(other)
    elif mutation == 'large_float':
        _shift_mean(row, 1e-6)
        other = experiment.evaluate(other['prediction'])
    elif mutation == 'type':
        row['steps'] = float(row['steps'])
        reseal(other)
    else:
        other['scores']['primary']['metrics']['combined']['position_estimator_rmse'] = float('nan')
    with pytest.raises(ValueError):
        replay.compare_receipt(receipt, other)


def test_changed_acceptance_requires_recomputed_local_conditional_certificate(receipt):
    other = deepcopy(receipt)
    other['prediction']['runs'][0]['predictions'][0]['acceptance_probability'] += 1e-13
    reseal(other)
    with pytest.raises(ValueError, match='Conditional exact'):
        replay.compare_receipt(receipt, other)


def _shift_mean(row, change):
    row['predicted_spin_mean'] += change
    row['predicted_position_estimator'] = row['predicted_spin_mean']/(math.sqrt(2)*.4)
    for bit, sign in ((0, 1), (1, -1)):
        row['predicted_detector_probabilities'][bit] += sign*change/2
        row['spin_detector_joint_probabilities'][0][bit] += sign*change*row['acceptance_probability']/2


def test_honest_primitive_neighbor_recomputes_amplified_estimator_scores_and_refinements(receipt):
    other = deepcopy(receipt)
    for index, change in ((0, -1.9e-10), (1, 1.9e-10), (2, -1.9e-10)):
        _shift_mean(other['prediction']['runs'][index]['predictions'][0], change)
    other = experiment.evaluate(other['prediction'])
    old = receipt['prediction']['runs'][0]['predictions'][0]
    new = other['prediction']['runs'][0]['predictions'][0]
    assert abs(old['predicted_spin_mean']-new['predicted_spin_mean']) < 2e-10
    assert abs(old['predicted_position_estimator']-new['predicted_position_estimator']) > 2e-10
    assert abs(receipt['numerical_refinements']['occupation']['maximum_position_estimator_change']
               - other['numerical_refinements']['occupation']['maximum_position_estimator_change']) > 2e-10
    before = deepcopy(other)
    replay.compare_receipt(receipt, other)
    replay.compare_receipt(other, receipt)
    assert other == before


@pytest.mark.parametrize('side', ['expected', 'actual', 'both'])
@pytest.mark.parametrize('field', ['estimator', 'score', 'refinement'])
def test_adjacent_float_forgery_fails_own_input_validation_even_when_both_match(receipt, side, field):
    corrupt = deepcopy(receipt)
    if field == 'estimator':
        target = corrupt['prediction']['runs'][0]['predictions'][0]
        key = 'predicted_position_estimator'
    elif field == 'score':
        target = corrupt['scores']['primary']['metrics']['combined']
        key = 'position_estimator_rmse'
    else:
        target = corrupt['numerical_refinements']['occupation']
        key = 'maximum_position_estimator_change'
    target[key] = math.nextafter(target[key], math.inf)
    reseal(corrupt)
    expected = corrupt if side in ('expected', 'both') else receipt
    actual = corrupt if side in ('actual', 'both') else receipt
    with pytest.raises(ValueError, match='Derived receipt'):
        replay.compare_receipt(expected, actual)


@pytest.mark.parametrize('mutation', ['nonfinite', 'type', 'omission', 'count', 'flag',
                                     'missing_estimator', 'missing_scores', 'missing_refinements'])
def test_derived_structure_and_metadata_are_not_exempted(receipt, mutation):
    other = deepcopy(receipt)
    if mutation == 'nonfinite':
        other['numerical_refinements']['step']['maximum_position_estimator_change'] = float('nan')
    elif mutation == 'type':
        other['numerical_refinements']['step']['maximum_binary_probability_change'] = 0
    elif mutation == 'omission':
        other['scores']['primary']['rows'].pop()
    elif mutation == 'count':
        other['scores']['primary']['metrics']['combined']['observations'] = 124
    elif mutation == 'flag':
        other['numerical_refinements']['occupation']['rigorous_occupation_error_bound'] = True
    elif mutation == 'missing_estimator':
        del other['prediction']['runs'][0]['predictions'][0]['predicted_position_estimator']
        reseal(other)
    elif mutation == 'missing_scores':
        del other['scores']
    else:
        del other['numerical_refinements']
    with pytest.raises(ValueError, match='Derived receipt|Nonfinite derived'):
        replay.compare_receipt(other, other)


def test_primitive_spin_guard_is_not_relaxed_by_dependent_recalculation(receipt):
    other = deepcopy(receipt)
    _shift_mean(other['prediction']['runs'][0]['predictions'][0], 2.1e-10)
    other = experiment.evaluate(other['prediction'])
    with pytest.raises(ValueError, match='Floating replay mismatch.*predicted_spin_mean'):
        replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('run_index,dimension', [(0, 32), (2, 64), (6, 32)])
@pytest.mark.parametrize('normalized_change,accepted', [(1.9e-10, True), (2.1e-10, False)])
def test_independent_position_observable_remains_compared_in_operator_norm_units(
        receipt, run_index, dimension, normalized_change, accepted):
    other = deepcopy(receipt)
    row = other['prediction']['runs'][run_index]['predictions'][0]
    change = normalized_change*math.sqrt(2*(dimension-1))
    assert change > 2e-10
    row['exact_position_mean_diagnostic'] += change
    reseal(other)
    if accepted:
        replay.compare_receipt(receipt, other)
    else:
        with pytest.raises(ValueError, match='Floating replay mismatch.*exact_position_mean_diagnostic'):
            replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), True, 100.])
def test_independent_position_observable_cannot_hide_invalid_values(receipt, bad):
    other = deepcopy(receipt)
    other['prediction']['runs'][0]['predictions'][0]['exact_position_mean_diagnostic'] = bad
    if math.isfinite(bad):
        reseal(other)
    with pytest.raises(ValueError):
        replay.compare_receipt(other, other)


def test_frozen_diagnostic_cannot_claim_a_larger_normalization_dimension(receipt):
    other = deepcopy(receipt)
    frozen = other['prediction']['runs'][-1]
    frozen['dimension'] = 1000000
    frozen['predictions'][0]['exact_position_mean_diagnostic'] = 100.
    reseal(other)
    with pytest.raises(ValueError, match='operator norm bound'):
        replay.compare_receipt(other, other)

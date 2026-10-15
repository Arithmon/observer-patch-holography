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

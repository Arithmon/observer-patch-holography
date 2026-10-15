"""Exact local seals, independently validated certificates, portable floats."""
from copy import deepcopy
from fractions import Fraction
import math

import pytest

from . import bounds, experiment, replay
from .test_experiment import build_prediction, synthetic_observations


@pytest.fixture(scope='module')
def receipt():
    with pytest.MonkeyPatch.context() as monkeypatch:
        prediction, _ = build_prediction(monkeypatch)
        monkeypatch.setattr(experiment, 'load_validation', synthetic_observations)
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


def test_numerical_prediction_roundoff_can_change_local_seal(receipt, monkeypatch):
    other = deepcopy(receipt)
    other['prediction']['runs'][0]['predictions'][0]['predicted_spin_mean'] += 1e-13
    monkeypatch.setattr(experiment, 'load_validation', synthetic_observations)
    other = experiment.evaluate(other['prediction'])
    replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('change, accepted', [(5e-9, True), (5e-8, False)])
def test_portable_tolerance_covers_measured_platform_drift_only(receipt, monkeypatch,
                                                                change, accepted):
    other = deepcopy(receipt)
    other['prediction']['runs'][0]['predictions'][0]['predicted_spin_mean'] += change
    monkeypatch.setattr(experiment, 'load_validation', synthetic_observations)
    other = experiment.evaluate(other['prediction'])
    if accepted:
        replay.compare_receipt(receipt, other)
    else:
        with pytest.raises(ValueError, match='Floating replay mismatch'):
            replay.compare_receipt(receipt, other)


def test_refinement_summary_propagates_individual_prediction_tolerance(receipt, monkeypatch):
    other = deepcopy(receipt)
    for run, change in ((other['prediction']['runs'][1], -9e-11),
                        (other['prediction']['runs'][2], 9e-11)):
        row = run['predictions'][0]
        row['predicted_spin_mean'] += change
        row['predicted_position_estimator'] += change/(math.sqrt(2)*.4)
        row['predicted_detector_probabilities'][0] += change/2
        row['predicted_detector_probabilities'][1] -= change/2
        row['spin_detector_joint_probabilities'][0][0] += .75*change/2
        row['spin_detector_joint_probabilities'][0][1] -= .75*change/2
    monkeypatch.setattr(experiment, 'load_validation', synthetic_observations)
    other = experiment.evaluate(other['prediction'])
    assert other['numerical_refinements']['occupation']['maximum_position_estimator_change'] > 2e-10
    replay.compare_receipt(receipt, other)


@pytest.mark.parametrize('field', ['maximum_binary_probability_change', 'maximum_position_estimator_change'])
def test_refinement_summary_cannot_be_changed_independently(receipt, field):
    other = deepcopy(receipt)
    other['numerical_refinements']['occupation'][field] += 1e-13
    with pytest.raises(ValueError, match='Refinement summary'):
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
        row['predicted_spin_mean'] += 1e-6
        reseal(other)
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

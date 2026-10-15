"""Prediction isolation and receipt integrity without expensive trajectories."""
from copy import deepcopy
import math
from types import SimpleNamespace

import numpy as np
import pytest

from . import experiment as e
from . import __main__ as cli


def _synthetic_history(apparatus, alpha, step):
    mean = .1*math.tanh(alpha)+step/500+apparatus.phi/100
    accepted = .75
    probabilities = np.array([[accepted*(1+mean)/2, accepted*(1-mean)/2],
                              [(1-accepted)/2, (1-accepted)/2]])
    reading = {
        'position_estimator': mean/(math.sqrt(2)*apparatus.h),
        'selected_mean': mean,
        'spin_detector_joint_probabilities': probabilities.tolist(),
        'acceptance_probability': accepted,
        'rejection_probability': 1-accepted,
    }
    branches = tuple(SimpleNamespace(density=np.diag([p, 0.]),
                                    resources=(('native_pulses', 11+step),))
                     for p in probabilities.ravel())
    return SimpleNamespace(reading=reading, branches=branches,
                           probability_mass=float(probabilities.sum()),
                           work_dimension=2*apparatus.dimension)


@pytest.fixture(scope='module')
def synthetic_prediction():
    """The producer gets apparatus inputs but cannot open later measurements."""
    with pytest.MonkeyPatch.context() as patch:
        calls = []

        def snapshots(apparatus, alpha, thermal, steps):
            calls.append((apparatus, alpha, thermal, steps))
            for step in steps:
                yield step, _synthetic_history(apparatus, alpha, step)

        patch.setattr(e.anharmonic, 'snapshots', snapshots)
        patch.setattr(e.anharmonic, 'readout', lambda history: history.reading)
        patch.setattr(e.calibration, 'transferred_parameters', lambda _: {
            'alpha_real': -.9, 'thermal_occupation': .1})
        patch.setattr(e.data, 'load_validation', lambda: pytest.fail(
            'The prediction producer opened later experimental values'))
        patch.setattr(e.strict_json, 'load', lambda _: pytest.fail(
            'The prediction producer opened an earlier prediction receipt'))
        prediction = e.predict_laboratory()
        return prediction, calls


def _observations():
    return [{'case': case, 'time_ms': step/5,
             'observed_position_estimator': .25}
            for case in 'abcde' for step in e.PRIMARY_STEPS[1:]]


def _oracle(prediction):
    return {'prediction': {'runs': [{'name': 'primary',
                                    'predictions': deepcopy(prediction['predictions']),
                                    'zero_time_predictions': deepcopy(prediction['zero_time_predictions'])}]}}


def _scoring_inputs(monkeypatch, prediction):
    monkeypatch.setattr(e.data, 'load_validation', _observations)
    monkeypatch.setattr(e.data, 'verify_provenance', lambda: {'synthetic': True})
    monkeypatch.setattr(e.strict_json, 'load', lambda _: _oracle(prediction))


def test_prediction_executes_all_endpoints_without_later_data_or_old_trajectories(synthetic_prediction):
    prediction, calls = synthetic_prediction
    assert len(calls) == 5
    assert [call[1] for call in calls] == [-.9, -.9, .9, 0., -.9]
    assert [call[0].phi/math.pi for call in calls] == [0., .05, .05, .05, .1]
    assert all(call[2:] == (.1, e.PRIMARY_STEPS) for call in calls)
    assert all((call[0].dimension, call[0].subdivisions) == (32, 128) for call in calls)
    assert len(prediction['predictions']) == 125
    assert len(prediction['zero_time_predictions']) == 5
    assert prediction['positive_time_outcomes_used'] is False
    assert prediction['earlier_endpoint_measurements_applied_to_later_trials'] is False
    assert next(row for row in prediction['predictions'] if row['steps'] == 12)['time_ms'] == float('2.4')


def test_prediction_is_sealed_before_observations_and_old_executor_open(synthetic_prediction, monkeypatch):
    prediction = deepcopy(synthetic_prediction[0])
    digest = e.digest
    seals, opened = [], []

    def seal(value):
        result = digest(value)
        seals.append(result)
        return result

    def observations():
        assert digest(prediction) in seals
        opened.append('measurements')
        return _observations()

    def oracle(_):
        assert digest(prediction) in seals
        assert opened == ['measurements']
        opened.append('old_executor')
        return _oracle(prediction)

    _scoring_inputs(monkeypatch, prediction)
    monkeypatch.setattr(e, 'digest', seal)
    monkeypatch.setattr(e.data, 'load_validation', observations)
    monkeypatch.setattr(e.strict_json, 'load', oracle)
    result = e.evaluate_laboratory(prediction)
    assert opened == ['measurements', 'old_executor']
    assert result['prediction_sha256_before_scoring'] == digest(prediction)
    assert result['scores']['metrics']['combined']['observations'] == 125
    assert result['independent_existing_executor_comparison']['maximum_joint_probability_residual'] == 0.
    assert result['empirical_exact_agreement_established'] is False
    assert result['physical_discrimination_from_conventional_quantum_mechanics'] is False


@pytest.mark.parametrize('mutation', [
    'missing_zero', 'duplicate_zero', 'missing_positive', 'duplicate_positive',
    'wrong_protocol', 'wrong_dimension', 'wrong_subdivisions', 'leak_flag',
    'negative_probability', 'probability_shape', 'nonfinite_probability',
    'bad_mass', 'bad_estimator', 'bad_acceptance', 'bad_purity', 'wrong_clock',
])
def test_invalid_prediction_rejected_before_opening_outcomes(synthetic_prediction, monkeypatch, mutation):
    prediction = deepcopy(synthetic_prediction[0])
    row = prediction['predictions'][0]
    if mutation == 'missing_zero':
        prediction['zero_time_predictions'].pop()
    elif mutation == 'duplicate_zero':
        prediction['zero_time_predictions'][-1] = deepcopy(prediction['zero_time_predictions'][0])
    elif mutation == 'missing_positive':
        prediction['predictions'].pop()
    elif mutation == 'duplicate_positive':
        prediction['predictions'][-1] = deepcopy(row)
    elif mutation == 'wrong_protocol':
        prediction['protocol_sha256'] = '0'*64
    elif mutation == 'wrong_dimension':
        prediction['dimension'] = 64
    elif mutation == 'wrong_subdivisions':
        prediction['subdivisions'] = True
    elif mutation == 'leak_flag':
        prediction['positive_time_outcomes_used'] = True
    elif mutation == 'negative_probability':
        row['spin_detector_joint_probabilities'][0][0] = -.1
    elif mutation == 'probability_shape':
        row['spin_detector_joint_probabilities'] = [.2, .3, .2, .3]
    elif mutation == 'nonfinite_probability':
        row['spin_detector_joint_probabilities'][0][0] = float('nan')
    elif mutation == 'bad_mass':
        row['total_probability_mass'] = .5
    elif mutation == 'bad_estimator':
        row['predicted_position_estimator'] += .1
    elif mutation == 'bad_acceptance':
        row['acceptance_probability'] = .1
    elif mutation == 'bad_purity':
        row['unnormalized_quantum_branch_purities'][0] = 2.
    else:
        row['drive_elapsed_seconds'] += .001
    monkeypatch.setattr(e.data, 'load_validation', lambda: pytest.fail(
        'An invalid prediction opened later experimental values'))
    with pytest.raises((ValueError, ArithmeticError)):
        e.evaluate_laboratory(prediction)


def test_scoring_cannot_mutate_the_sealed_prediction(synthetic_prediction, monkeypatch):
    prediction = deepcopy(synthetic_prediction[0])
    _scoring_inputs(monkeypatch, prediction)
    score = e.prior.score

    def mutating_score(rows, observations):
        result = score(rows, observations)
        rows[0]['unexecuted_annotation'] = 'changed after seal'
        return result

    monkeypatch.setattr(e.prior, 'score', mutating_score)
    with pytest.raises(ArithmeticError, match='changed during scoring'):
        e.evaluate_laboratory(prediction)


@pytest.mark.parametrize('rows_key', ['predictions', 'zero_time_predictions'])
def test_old_executor_comparison_covers_rejected_outcomes_at_every_endpoint(
        synthetic_prediction, monkeypatch, rows_key):
    prediction = deepcopy(synthetic_prediction[0])
    _scoring_inputs(monkeypatch, prediction)
    oracle = _oracle(prediction)
    rejected = oracle['prediction']['runs'][0][rows_key][0]['spin_detector_joint_probabilities'][1]
    # The scored accepted estimator is identical; this changes only the
    # unselected detector law that a comparison of means would miss.
    rejected[0] += 1e-4
    rejected[1] -= 1e-4
    monkeypatch.setattr(e.strict_json, 'load', lambda _: oracle)
    with pytest.raises(ArithmeticError, match='changed the frozen apparatus'):
        e.evaluate_laboratory(prediction)


def test_changed_observations_change_scores_without_changing_prediction(synthetic_prediction, monkeypatch):
    prediction = deepcopy(synthetic_prediction[0])
    _scoring_inputs(monkeypatch, prediction)
    first = e.evaluate_laboratory(prediction)
    observations = _observations()
    for row in observations:
        row['observed_position_estimator'] += 10.
    monkeypatch.setattr(e.data, 'load_validation', lambda: observations)
    second = e.evaluate_laboratory(prediction)
    assert first['prediction'] == second['prediction'] == prediction
    assert first['prediction_sha256_before_scoring'] == second['prediction_sha256_before_scoring']
    assert first['scores'] != second['scores']


def _minimal_result():
    prediction = {'a_primitive_probability': .5}
    return {'all_pass': True, 'checks': {'finite_contract': True},
            'laboratory': {'prediction': prediction,
                           'prediction_sha256_before_scoring': e.digest(prediction),
                           'scores': {'metrics': {'combined': {
                               'observations': 125, 'position_estimator_rmse': .28}}}},
            'scope': {'all_universe_observations_reproduced': False}}


def test_portable_receipt_accepts_small_float_changes_only_with_its_own_seal():
    expected = _minimal_result()
    actual = deepcopy(expected)
    actual['laboratory']['prediction']['a_primitive_probability'] += 1e-12
    with pytest.raises(ValueError, match='seal'):
        cli.compare(expected, actual)
    actual['laboratory']['prediction_sha256_before_scoring'] = e.digest(actual['laboratory']['prediction'])
    cli.compare(expected, actual)
    actual['laboratory']['prediction']['a_primitive_probability'] += 1e-7
    actual['laboratory']['prediction_sha256_before_scoring'] = e.digest(actual['laboratory']['prediction'])
    with pytest.raises(ValueError, match='replay mismatch'):
        cli.compare(expected, actual)


@pytest.mark.parametrize('failure', ['result', 'check'])
def test_cli_cannot_publish_a_failed_contract(monkeypatch, tmp_path, failure):
    result = _minimal_result()
    if failure == 'result':
        result['all_pass'] = False
    else:
        result['checks']['finite_contract'] = False
    destination = tmp_path/'must_not_exist.json'
    monkeypatch.setattr(cli, 'run', lambda **_: result)
    monkeypatch.setattr('sys.argv', ['primitive_observer', '--write', str(destination)])
    with pytest.raises((ValueError, ArithmeticError, SystemExit)):
        cli.main()
    assert not destination.exists()


def test_actual_producer_bindings_match_retained_portable_paths(monkeypatch):
    from . import action, completion_report, observer_report, selection

    retained = e.strict_json.load(e.FOLDER/'results.json')
    # Isolate expensive numerical components, but exercise run's real source
    # discovery, path conversion, deduplication and file hashing. On Windows
    # the former str(relative_path) produced backslashes and duplicate aliases.
    for module, key in ((action, 'native_action'), (selection, 'pointer_selection'),
                        (completion_report, 'completion'), (observer_report, 'observer')):
        monkeypatch.setattr(module, 'run', lambda key=key: deepcopy(retained[key]))
    monkeypatch.setattr(e, 'predict_laboratory', lambda **_: None)
    monkeypatch.setattr(e, 'evaluate_laboratory', lambda _: deepcopy(retained['laboratory']))
    actual = e.run()['source_bindings']
    assert all('\\' not in path for path in actual)
    assert actual == retained['source_bindings']

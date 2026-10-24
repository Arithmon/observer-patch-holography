"""Prediction isolation and complete scoring with mocked source execution."""
from copy import deepcopy
import hashlib
import json
import math

import pytest

from . import bounds, experiment as e


def fake_execution(dimension, alpha, thermal, steps, settings, **kwargs):
    """Consistent synthetic joint laws, never empirical response values."""
    records = []
    for step in steps:
        mean = .1*math.tanh(alpha)+step/500
        probabilities = [(1+mean)/2, (1-mean)/2]
        records.append({'steps': step, 'simulated_duration_seconds': step*.0002,
                        'drive_elapsed_duration_seconds': step*.000203,
                        'spin_probabilities': [.75, .25], 'selected_position_mean': .123,
                        'execution_error_budget': bounds.execution_error_budget(
                            dimension, step, kwargs['subdivisions'], .75),
                        'reads': [{**settings[0], 'selected_mean': mean,
                                   'selected_detector_probabilities': probabilities,
                                   'spin_detector_joint_probabilities': [
                                       [p*.75 for p in probabilities], [.125, .125]]}]})
    return {'records': records, 'dimension': dimension, 'alpha': alpha,
            'thermal': thermal, 'requested': kwargs}


def build_prediction(monkeypatch):
    monkeypatch.setattr(e, 'load_metadata', lambda: {'labels': list('abcde'),
                        'times_ms': [index*2/5 for index in range(26)], 'provenance': {'synthetic': True}})
    monkeypatch.setattr(e, 'load_zero', lambda: dict(zip('abcde', (-1., -.8, .8, 0., -1.1))))
    # Keep the authenticated prior preparation so the independent occupation
    # certificate describes the same apparatus; only source execution is fake.
    calls = []
    def execute(*args, **kwargs):
        calls.append((args, kwargs))
        return fake_execution(*args, **kwargs)
    monkeypatch.setattr(e.source, 'run', execute)
    monkeypatch.setattr(e, 'load_validation', lambda: pytest.fail('Later outcomes entered prediction'))
    return e.predict(), calls


def synthetic_observations():
    return [{'case': case, 'time_ms': step/5, 'observed_position_estimator': .25,
             'observed_position_estimator_decimal': '0.25'}
            for step in range(2, 51, 2) for case in 'abcde']


@pytest.fixture(scope='module')
def fixed_prediction():
    with pytest.MonkeyPatch.context() as monkeypatch:
        return build_prediction(monkeypatch)


def test_prediction_never_opens_later_values_and_executes_every_predeclared_run(fixed_prediction):
    prediction, calls = fixed_prediction
    assert not prediction['positive_time_outcomes_used']
    assert len(calls) == 30
    assert [run['name'] for run in prediction['runs']] == list(e.RUN_NAMES)
    assert all(len(run['predictions']) == 125 and len(run['zero_time_predictions']) == 5
               for run in prediction['runs'])
    assert prediction['protocol_sha256'] == hashlib.sha256((e.FOLDER/'protocol.json').read_bytes()).hexdigest()
    for path, expected in prediction['source_bindings'].items():
        assert expected == hashlib.sha256((e.ROOT/path).read_bytes()).hexdigest()


def test_primary_transfer_changes_only_declared_preparation_sign_and_centre(fixed_prediction):
    prediction, calls = fixed_prediction
    alpha = prediction['prior_parameters']['alpha_real']
    thermal = prediction['prior_parameters']['thermal_occupation']
    assert [args[1] for args, _ in calls[:5]] == [alpha, alpha, -alpha, 0., alpha]
    assert all(args[2] == thermal for args, _ in calls)
    assert [kwargs['phi']/math.pi for _, kwargs in calls[:5]] == [0., .05, .05, .05, .1]
    adjusted = [args[1] for args, _ in calls[15:20]]
    assert adjusted == [prediction['preparations'][case]['zero_time_adjusted_alpha'] for case in 'abcde']
    assert adjusted != [args[1] for args, _ in calls[:5]]


def test_refinements_vary_one_axis_and_controls_keep_declared_settings(fixed_prediction):
    _, calls = fixed_prediction
    assert [(calls[index][0][0], calls[index][1]['subdivisions']) for index in (0, 5, 10)] == [
        (32, 128), (32, 256), (64, 256)]
    assert all(kwargs['finite_spin_duration'] is False for _, kwargs in calls[20:25])
    assert all(kwargs['omega'] == 0 and kwargs['theta'] == 0 for _, kwargs in calls[25:30])


@pytest.mark.parametrize('section,key,value', [
    ('apparatus', 'theta', .81), ('apparatus', 'schedule', ['wrong chronology']),
    ('readout', 'h', .45), ('readout', 'estimator', 'exact position'),
    ('numerics', 'primary_dimension', 64), ('numerics', 'quadrature_order', True),
    ('preparation', 'prior_receipt', 'another.json'), ('preparation', 'secondary', 'fit all later values'),
])
def test_unsupported_protocol_changes_rejected(section, key, value):
    protocol = json.loads((e.FOLDER/'protocol.json').read_text())
    protocol[section][key] = value
    with pytest.raises(ValueError, match='Unsupported'):
        e._configuration(protocol)


def test_complete_scoring_uses_literal_decimal_coordinates(fixed_prediction, monkeypatch):
    prediction, _ = fixed_prediction
    # 12*.2 is 2.4000000000000004; the CSV coordinate is float('2.4').
    assert next(row for row in prediction['runs'][0]['predictions'] if row['steps'] == 12)['time_ms'] == float('2.4')
    seal_before = e.digest(prediction)
    def open_validation():
        assert e.digest(prediction) == seal_before
        return synthetic_observations()
    monkeypatch.setattr(e, 'load_validation', open_validation)
    result = e.evaluate(prediction)
    assert result['all_pass']
    assert result['prediction_payload_sha256_before_scoring'] == seal_before
    assert e.digest(prediction) == seal_before
    assert set(result['scores']) == set(e.RUN_NAMES)
    for score in result['scores'].values():
        assert score['metrics']['combined']['observations'] == 125
        assert all(score['metrics'][case]['observations'] == 25 for case in 'abcde')
        assert not score['statistical_significance_claimed']
    assert result['numerical_refinements']['occupation']['runs'] == ['step_refinement', 'occupation_refinement']


@pytest.mark.parametrize('mutation', ['protocol_hash', 'protocol_content', 'run_inventory',
                                    'target_missing', 'target_duplicate', 'dimension', 'subdivisions', 'leak_flag'])
def test_prediction_contract_fails_before_later_values_open(fixed_prediction, monkeypatch, mutation):
    prediction = deepcopy(fixed_prediction[0])
    monkeypatch.setattr(e, 'load_validation', lambda: pytest.fail('Invalid predictions opened later data'))
    if mutation == 'protocol_hash':
        prediction['protocol_sha256'] = '0'*64
    elif mutation == 'protocol_content':
        prediction['protocol']['readout']['h'] = .5
    elif mutation == 'run_inventory':
        prediction['runs'].pop()
    elif mutation == 'target_missing':
        prediction['runs'][0]['predictions'].pop()
    elif mutation == 'target_duplicate':
        prediction['runs'][0]['predictions'][0] = prediction['runs'][0]['predictions'][1]
    elif mutation in ('dimension', 'subdivisions'):
        prediction['runs'][1][mutation] = 64
    else:
        prediction['positive_time_outcomes_used'] = True
    with pytest.raises(ValueError):
        e.evaluate(prediction)


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'nonfinite', 'unrelated_time'])
def test_scoring_rejects_missing_extra_or_invalid_observations(fixed_prediction, mutation):
    observations = synthetic_observations()
    if mutation == 'missing':
        observations.pop()
    elif mutation == 'duplicate':
        observations.append(observations[0])
    elif mutation == 'nonfinite':
        observations[0]['observed_position_estimator'] = float('nan')
    else:
        observations[0]['time_ms'] += 1e-6
    with pytest.raises(ValueError):
        e.score(fixed_prediction[0]['runs'][0]['predictions'], observations)


def test_frozen_control_reuses_zero_time_law_and_keeps_requested_coordinate(fixed_prediction):
    prediction, _ = fixed_prediction
    primary, frozen = prediction['runs'][0], prediction['runs'][-1]
    initial = {row['case']: row for row in primary['zero_time_predictions']}
    for row in frozen['predictions']:
        assert row['steps'] == 0 and row['requested_steps'] > 0
        assert row['time_ms'] == row['requested_steps']/5
        assert row['predicted_detector_probabilities'] == initial[row['case']]['predicted_detector_probabilities']
        assert row['drive_elapsed_seconds'] == 0

"""Scoring and prediction isolation on synthetic observations only."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest

from . import experiment


def _observations():
    characteristic_beta = np.linspace(0., 5., 50)
    position_beta = np.linspace(-.4, .4, 13)/math.sqrt(2)
    return {
        'characteristic': {
            'times_ms': np.array([2., 4.]),
            'probe_coordinates': characteristic_beta,
            'spin_real': np.full((2, 50), .2),
            'spin_imag': np.full((2, 50), -.3)},
        'position': {
            'times_ms': np.arange(1, 40)*.4,
            'probe_coordinates': position_beta,
            'spin_imag': np.tile(.4*position_beta+.02, (39, 1))}}


def _predictions(observations):
    rows = []
    for family, components in (('characteristic', ('real', 'imag')), ('position', ('imag',))):
        data = observations[family]
        for index, time in enumerate(data['times_ms']):
            for probe, beta in enumerate(data['probe_coordinates']):
                for component in components:
                    residual = .2 if family == 'characteristic' else -.1
                    mean = float(data['spin_'+component][index, probe]+residual)
                    probabilities = [(1+mean)/2, (1-mean)/2]
                    rows.append({'family': family, 'time_ms': float(time),
                                 'beta_imag': float(beta), 'component': component,
                                 'predicted_mean': mean,
                                 'predicted_probabilities': probabilities,
                                 'spin_detector_joint_probabilities': [
                                     [.71*value for value in probabilities], [.145, .145]],
                                 'acceptance_probability': .71, 'rejection_probability': .29})
    return rows


def _fake_run(dimension, alpha, thermal, steps, settings, **kwargs):
    records = []
    for step in steps:
        acceptance = .85-.001*step
        reads = []
        for setting in settings:
            mean = .1*math.sin(setting['v'])+.001*step
            probabilities = [(1+mean)/2, (1-mean)/2]
            reads.append({**setting, 'selected_mean': mean,
                          'selected_detector_probabilities': probabilities,
                          'spin_detector_joint_probabilities': [
                              [acceptance*p for p in probabilities],
                              [(1-acceptance)/2, (1-acceptance)/2]]})
        records.append({'steps': step, 'spin_probabilities': [acceptance, 1-acceptance],
                        'reads': reads})
    return {'dimension': dimension, 'records': records, 'settings': kwargs,
            'preparation': {'alpha': alpha, 'thermal': thermal}}


@pytest.fixture
def isolated_prediction(monkeypatch, tmp_path):
    observations = _observations()
    metadata = deepcopy(observations)
    zero = deepcopy(observations)
    for family in ('characteristic', 'position'):
        metadata[family] = {'times_ms': np.r_[0., observations[family]['times_ms']],
                            'probe_coordinates': observations[family]['probe_coordinates']}
        zero[family]['times_ms'] = np.array([0.])
        for component in ('real', 'imag'):
            key = 'spin_'+component
            if key in zero[family]:
                zero[family][key] = zero[family][key][:1].copy()
    metadata['provenance'] = {'test_only': True}
    protocol = json.loads((experiment.FOLDER/'protocol.json').read_text())
    (tmp_path/'protocol.json').write_text(json.dumps(protocol))
    monkeypatch.setattr(experiment, 'FOLDER', tmp_path)
    monkeypatch.setattr(experiment, 'load_metadata', lambda: deepcopy(metadata))
    monkeypatch.setattr(experiment, 'load_calibration', lambda: deepcopy(zero))
    calls = []

    def fake_calibrate(beta, real, imag):
        np.testing.assert_array_equal(beta, zero['characteristic']['probe_coordinates'])
        np.testing.assert_array_equal(real, zero['characteristic']['spin_real'][0])
        np.testing.assert_array_equal(imag, zero['characteristic']['spin_imag'][0])
        return {'parameters': {'alpha_real': -.7, 'thermal_occupation': .12},
                'test_only': True}

    def fake_run(*args, **kwargs):
        calls.append((args, kwargs))
        return _fake_run(*args, **kwargs)

    def forbidden_validation():
        raise AssertionError('positive-time observations were read during prediction')

    monkeypatch.setattr(experiment, 'calibrate', fake_calibrate)
    monkeypatch.setattr(experiment.source, 'run', fake_run)
    monkeypatch.setattr(experiment, 'load_validation', forbidden_validation)
    return observations, protocol, calls, tmp_path


def test_every_707_mean_is_scored_once_with_equal_weight():
    observations = _observations()
    comparison = experiment.score(_predictions(observations), observations)
    assert comparison['characteristic']['scalar_means'] == 200
    assert comparison['position']['scalar_means'] == 507
    assert comparison['combined']['scalar_means'] == 707
    assert comparison['characteristic']['spin_mean_rmse'] == pytest.approx(.2)
    assert comparison['position']['spin_mean_rmse'] == pytest.approx(.1)
    expected = math.sqrt((200*.2**2+507*.1**2)/707)
    assert comparison['combined']['spin_mean_rmse'] == pytest.approx(expected)
    assert comparison['combined']['binary_probability_rmse'] == pytest.approx(expected/2)
    assert len(comparison['derived_position']) == 39
    for row in comparison['derived_position']:
        assert row['observed_position'] == pytest.approx(.4/math.sqrt(2))
        assert row['predicted_position'] == pytest.approx(.4/math.sqrt(2))
    assert comparison['statistical_significance_claimed'] is False


@pytest.mark.parametrize('corruption', ['missing', 'duplicate', 'wrong_setting', 'wrong_component'])
def test_scoring_rejects_incomplete_or_duplicate_prediction_inventory(corruption):
    observations = _observations()
    rows = _predictions(observations)
    if corruption == 'missing':
        rows.pop()
    elif corruption == 'duplicate':
        rows.append(deepcopy(rows[0]))
    elif corruption == 'wrong_setting':
        rows[0]['beta_imag'] = 17.
    else:
        rows[0]['component'] = 'unmeasured'
    with pytest.raises(ValueError, match='exactly once'):
        experiment.score(rows, observations)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), 1.1])
def test_scoring_rejects_nonphysical_observed_means(value):
    observations = _observations()
    rows = _predictions(observations)
    observations['position']['spin_imag'][0, 0] = value
    with pytest.raises(ValueError, match='binary spin means'):
        experiment.score(rows, observations)


def _selected_example():
    targets = [{'family': 'position', 'time_ms': .4, 'beta_imag': .2, 'component': 'imag'}]
    execution = {'records': [{'steps': 2, 'spin_probabilities': [.25, .75], 'reads': [
        {'u': 0., 'v': .2, 'theta': -math.pi/2, 'selected_mean': .6,
         'selected_detector_probabilities': [.8, .2],
         'spin_detector_joint_probabilities': [[.2, .05], [.3, .45]]}]}]}
    return execution, targets


def test_extraction_uses_conditional_mean_and_retains_both_joint_branches():
    execution, targets = _selected_example()
    row = experiment._extract(execution, targets)[0]
    assert row['predicted_mean'] == .6
    assert row['acceptance_probability'] == .25
    assert row['rejection_probability'] == .75
    assert row['predicted_probabilities'] == [.8, .2]
    assert row['spin_detector_joint_probabilities'] == [[.2, .05], [.3, .45]]
    assert np.sum(row['spin_detector_joint_probabilities'], axis=0).tolist() == [.5, .5]


@pytest.mark.parametrize('corruption', ['normalization', 'conditional', 'joint', 'mean', 'spin_mass'])
def test_extraction_rejects_inconsistent_probability_laws(corruption):
    execution, targets = _selected_example()
    snapshot = execution['records'][0]
    read = snapshot['reads'][0]
    if corruption == 'normalization':
        read['selected_detector_probabilities'][0] = .7
    elif corruption == 'conditional':
        read['selected_detector_probabilities'] = [.5, .5]
        read['selected_mean'] = 0.
    elif corruption == 'joint':
        read['spin_detector_joint_probabilities'][1][1] = .4
    elif corruption == 'mean':
        read['selected_mean'] = .3
    else:
        snapshot['spin_probabilities'] = [.3, .7]
    with pytest.raises(ArithmeticError, match='selected-detector law'):
        experiment._extract(execution, targets)


def test_prediction_never_reads_validation_and_maps_all_refinements(isolated_prediction):
    _, protocol, calls, _ = isolated_prediction
    prediction = experiment.predict()
    assert prediction['positive_duration_values_used'] is False
    assert prediction['protocol'] == protocol
    assert len(calls) == 5
    assert [args[0] for args, _ in calls] == [32, 64, 32, 32, 32]
    assert [kwargs['subdivisions'] for _, kwargs in calls] == [8, 8, 16, 8, 8]
    assert all(args[1:3] == (-.7, .12) for args, _ in calls)
    assert all(args[3] == list(range(0, 79, 2)) for args, _ in calls)
    assert calls[3][1]['omega'] == calls[3][1]['theta'] == 0.
    assert calls[4][1]['dephasing_rate'] == 0.
    assert all(len(run['predictions']) == 707 for run in prediction['runs'])
    assert len(prediction['zero_duration_position_predictions']) == 13


def test_frozen_control_preserves_the_prepared_readout_at_every_time(isolated_prediction):
    prediction = experiment.predict()
    frozen = next(run for run in prediction['runs'] if run['name'] == 'frozen_control')
    assert frozen['source_execution']['evolution_executed'] is False
    for row in frozen['predictions']:
        assert row['predicted_mean'] == pytest.approx(.1*math.sin(row['beta_imag']))
        assert row['acceptance_probability'] == .85
    primary = next(run for run in prediction['runs'] if run['name'] == 'primary')
    assert any(a['predicted_mean'] != b['predicted_mean']
               for a, b in zip(primary['predictions'], frozen['predictions'], strict=True))


def test_scoring_does_not_mutate_frozen_prediction_payload(isolated_prediction, monkeypatch):
    observations, _, _, _ = isolated_prediction
    prediction = experiment.predict()
    original = deepcopy(prediction)
    monkeypatch.setattr(experiment, 'load_validation', lambda: deepcopy(observations))
    result = experiment.evaluate(prediction)
    assert prediction == original
    assert result['prediction_payload_sha256_before_scoring'] == experiment._digest(original)
    assert result['checks']['all_707_later_statistics_retained'] is True
    assert result['all_pass'] is True
    from .replay import prediction_payload, validate_prediction_seal
    assert prediction_payload(result) == original
    assert validate_prediction_seal(result) == result['prediction_payload_sha256_before_scoring']
    for refinement in ('occupation_refinement', 'dephasing_refinement'):
        assert result['numerical_refinements'][refinement]['maximum_spin_mean_change'] == 0.
        assert result['numerical_refinements'][refinement]['rigorous_error_bound'] is False


@pytest.mark.parametrize('axis', ['times_ms', 'probe_coordinates'])
def test_scoring_rejects_duplicate_observation_coordinates(axis):
    observations = _observations()
    rows = _predictions(observations)
    observations['characteristic'][axis][1] = observations['characteristic'][axis][0]
    with pytest.raises(ValueError, match='duplicate|unique'):
        experiment.score(rows, observations)


@pytest.mark.parametrize('corruption', ['snapshot', 'detector', 'real_coordinate', 'phase',
                                       'nonfinite_mean', 'zero_acceptance'])
def test_extraction_rejects_ambiguous_or_undefined_source_records(corruption):
    execution, targets = _selected_example()
    snapshot = execution['records'][0]
    read = snapshot['reads'][0]
    if corruption == 'snapshot':
        execution['records'].append(deepcopy(snapshot))
    elif corruption == 'detector':
        snapshot['reads'].append(deepcopy(read))
    elif corruption == 'real_coordinate':
        read['u'] = .3
    elif corruption == 'phase':
        read['theta'] = .7
    elif corruption == 'nonfinite_mean':
        read['selected_mean'] = float('nan')
    else:
        snapshot['spin_probabilities'] = [0., 1.]
        read['spin_detector_joint_probabilities'] = [[0., 0.], [.4, .6]]
    with pytest.raises((ValueError, ArithmeticError)):
        experiment._extract(execution, targets)


@pytest.mark.parametrize(('key', 'value'), [
    ('delta_radians_per_second', '2*pi*501'),
    ('omega_radians_per_second', 'pi/0.000151'),
    ('pulse_duration_seconds', .00006),
    ('step_duration_seconds', .0003),
    ('spin_rotation_angle', .9),
    ('asymmetry_phase', .1),
    ('dephasing_per_second', 19.),
    ('dephasing_operator', 'sqrt(36) * N'),
])
def test_unsupported_apparatus_protocol_changes_fail_before_execution(
        isolated_prediction, key, value):
    _, protocol, calls, path = isolated_prediction
    protocol['apparatus'][key] = value
    (path/'protocol.json').write_text(json.dumps(protocol))
    with pytest.raises(ValueError, match='apparatus|protocol'):
        experiment.predict()
    assert calls == []


def test_verified_apparatus_values_are_explicit_source_arguments(isolated_prediction):
    _, _, calls, _ = isolated_prediction
    experiment.predict()
    expected = {'delta': 2*math.pi*500, 'omega': math.pi/.000150,
                'pulse_duration': .00005, 'theta': .8, 'phi': 0., 'dephasing_rate': 18.}
    assert {key: calls[0][1][key] for key in expected} == expected
    assert {key: calls[1][1][key] for key in expected} == expected
    assert {key: calls[2][1][key] for key in expected} == expected


def test_numeric_protocol_values_drive_all_run_configurations(isolated_prediction):
    _, protocol, calls, path = isolated_prediction
    protocol['numerics'].update(primary_occupation_dimension=4, occupation_check_dimension=8,
                                dephasing_subdivisions_per_50us=2,
                                dephasing_refinement_subdivisions_per_50us=4)
    (path/'protocol.json').write_text(json.dumps(protocol))
    experiment.predict()
    assert [args[0] for args, _ in calls] == [4, 8, 4, 4, 4]
    assert [kwargs['subdivisions'] for _, kwargs in calls] == [2, 2, 4, 2, 2]


@pytest.mark.parametrize('corruption', ['missing', 'duplicate', 'unknown'])
def test_evaluation_rejects_run_inventory_changes_before_unsealing(isolated_prediction, corruption):
    prediction = experiment.predict()
    if corruption == 'missing':
        prediction['runs'].pop()
    elif corruption == 'duplicate':
        prediction['runs'].append(deepcopy(prediction['runs'][0]))
    else:
        prediction['runs'][0]['name'] = 'not_the_primary'
    with pytest.raises(ValueError, match='run|inventory'):
        experiment.evaluate(prediction)


def test_refinement_uses_primary_name_even_if_runs_are_reordered(isolated_prediction, monkeypatch):
    observations, _, _, _ = isolated_prediction
    prediction = experiment.predict()
    refinement = next(run for run in prediction['runs'] if run['name'] == 'occupation_refinement')
    for row in refinement['predictions']:
        row['predicted_mean'] += .02
        row['predicted_probabilities'] = [(1+row['predicted_mean'])/2, (1-row['predicted_mean'])/2]
        row['spin_detector_joint_probabilities'][0] = [
            row['acceptance_probability']*p for p in row['predicted_probabilities']]
    prediction['runs'].remove(refinement)
    prediction['runs'].insert(0, refinement)
    monkeypatch.setattr(experiment, 'load_validation', lambda: deepcopy(observations))
    result = experiment.evaluate(prediction)
    assert result['numerical_refinements']['occupation_refinement']['maximum_spin_mean_change'] == pytest.approx(.02)


def test_compacted_operation_identity_does_not_hash_float_diagnostics():
    first = {'name': 'example', 'word_sha256': 'a'*64, 'register_qubits': 4,
             'full_phase_target_residual': 1e-13}
    second = dict(first, full_phase_target_residual=2e-13)
    first_words, second_words = {}, {}
    first_reference = experiment._compact_operations(first, first_words)
    second_reference = experiment._compact_operations(second, second_words)
    assert first_reference == second_reference
    assert first_words != second_words
    assert list(first_words) == list(second_words)


def test_conflicting_metadata_cannot_replace_a_contracted_word_within_a_run():
    words = {}
    operation = {'name': 'example', 'word_sha256': 'a'*64, 'register_qubits': 4,
                 'full_phase_target_residual': 1e-13}
    first_reference = experiment._compact_operations(operation, words)
    assert experiment._compact_operations(operation, words) == first_reference
    with pytest.raises(ArithmeticError, match='Conflicting metadata'):
        experiment._compact_operations(dict(operation, full_phase_target_residual=2e-13), words)


def test_compacted_words_keep_register_and_role_distinctions():
    words = {}
    operation = {'name': 'example', 'word_sha256': 'a'*64, 'register_qubits': 4}
    references = [experiment._compact_operations(value, words)
                  for value in (operation, dict(operation, name='other'),
                                dict(operation, register_qubits=5))]
    assert len({row['operation_metadata_sha256'] for row in references}) == 3

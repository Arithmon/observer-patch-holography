"""Fresh end-to-end audit of inference access and coupled held-out histories."""
from copy import deepcopy
from fractions import Fraction
import hashlib
from itertools import product
import json

import numpy as np
import pytest

from . import execution, experiment, recovery, source
from .prediction import RecoveredField


def canonical(value):
    return json.dumps(experiment.serializable(value), sort_keys=True, allow_nan=False)


def probabilities(rows):
    assert [row['bits'] for row in rows] == [[0, 0], [0, 1], [1, 0], [1, 1]]
    return np.array([row['probability'] for row in rows])


@pytest.fixture(scope='module')
def executed():
    """Execute once; wrappers inspect calls and never substitute their values."""
    events, calibrations, coupled, inversions, predictions = [], [], [], [], []
    original_calibration = source.calibration_records
    original_coupled = execution.coupled_history
    original_recover = recovery.recover
    original_history = RecoveredField.history
    original_audit = execution.target_audit

    def forbidden(*args, **kwargs):
        raise AssertionError('inverse or prediction accessed a hidden source/target')

    def block_backend(patch):
        patch.setattr(source, 'supplied_target', forbidden)
        patch.setattr(source, 'calibration_records', forbidden)
        patch.setattr(source, 'onsite_probe', forbidden)
        patch.setattr(source, 'bond_probe', forbidden)
        patch.setattr(source.c, 'execute', forbidden)
        patch.setattr(execution, 'coupled_history', forbidden)
        patch.setattr(execution, 'target_audit', forbidden)

    def calibration(duration, *args, **kwargs):
        events.append(('calibration', duration))
        value = original_calibration(duration, *args, **kwargs)
        calibrations.append(deepcopy(value))
        return value

    def infer(*args, **kwargs):
        events.append(('inverse', kwargs.get('probability_error', 0.)))
        before = canonical(args)
        with pytest.MonkeyPatch.context() as patch:
            block_backend(patch)
            value = original_recover(*args, **kwargs)
        assert canonical(args) == before
        inversions.append({'arguments': deepcopy(args), 'kwargs': deepcopy(kwargs)})
        return value

    def predict(field, duration, **kwargs):
        record = {'dimensions': field.dimensions, 'duration': duration,
                  'sign': kwargs.get('sign', -1),
                  'feedback_phase': kwargs.get('feedback_phase', np.pi/2),
                  'quartic_deleted': bool(np.all(field.quartic == 0))}
        events.append(('prediction', *record.values()))
        with pytest.MonkeyPatch.context() as patch:
            block_backend(patch)
            value = original_history(field, duration, **kwargs)
        predictions.append({**record, 'rows': deepcopy(value)})
        return value

    def evolve(dimensions, duration, **kwargs):
        record = {'dimensions': tuple(dimensions), 'duration': duration,
                  'sign': kwargs.get('sign', -1),
                  'feedback_phase': kwargs.get('feedback_phase', np.pi/2)}
        events.append(('coupled', *record.values()))
        rows, metadata = original_coupled(dimensions, duration, **kwargs)
        coupled.append({**record, 'rows': deepcopy(rows), 'metadata': deepcopy(metadata)})
        return rows, metadata

    def audit():
        events.append(('audit',))
        return original_audit()

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(source, 'calibration_records', calibration)
        patch.setattr(recovery, 'recover', infer)
        patch.setattr(RecoveredField, 'history', predict)
        patch.setattr(execution, 'coupled_history', evolve)
        patch.setattr(execution, 'target_audit', audit)
        result = experiment.run()
    return {'result': result, 'events': events, 'calibrations': calibrations,
            'coupled': coupled, 'inversions': inversions, 'predictions': predictions}


def test_only_retained_calibration_and_clock_records_enter_inverse(executed):
    result = executed['result']
    assert len(executed['calibrations']) == 2
    assert len(executed['inversions']) == 2  # exact-law inference plus sensitivity
    for inverse in executed['inversions']:
        assert len(inverse['arguments']) == 4
        for actual, retained in zip(inverse['arguments'][:2], result['calibration_records'], strict=True):
            assert canonical(actual) == canonical(retained)
        for actual, retained in zip(inverse['arguments'][2:],
                                     result['calibration_clock_probabilities'], strict=True):
            np.testing.assert_array_equal(actual, np.asarray(retained)[:, 0])
        assert inverse['kwargs']['omega'] == 2.
        assert inverse['kwargs']['norm_bound_onsite'] == 16.
        assert inverse['kwargs']['norm_bound_bond'] == 2.
    assert executed['inversions'][0]['kwargs'].get('probability_error', 0.) == 0.
    assert executed['inversions'][1]['kwargs']['probability_error'] == 1e-7
    for raw, retained in zip(executed['calibrations'], result['calibration_records'], strict=True):
        assert canonical(raw) == canonical(retained)
        assert retained['fresh_trials'] == 88
        assert len(retained['pairs']) == 28
        for key, shape in [('onsite_probabilities', (8, 2, 2, 2)),
                           ('bond_probabilities', (28, 2, 2))]:
            laws = np.asarray(retained[key])
            assert laws.shape == shape
            assert laws.min() > 0
            np.testing.assert_allclose(laws.sum(axis=-1), 1., atol=4e-13, rtol=0)
        resources = retained['resources']
        assert resources['native_pulses'] > 0
        assert resources['proper_code_service_events'] == 2*resources['native_pulses']


def test_predictions_precede_withheld_executions_and_final_target_audit(executed):
    events = executed['events']
    assert events[:3] == [('calibration', t) for t in experiment.CALIBRATION_DURATIONS]+[('inverse', 0.)]
    cursor = 3
    for dimensions in experiment.LEVELS:
        for sign in (-1, 1):
            for deleted in (False, True):
                event = events[cursor]
                assert event[0] == 'prediction'
                assert event[1] == dimensions
                assert event[2] == pytest.approx(float(experiment.INTERVAL), abs=3e-13)
                assert event[3:] == (sign, np.pi/2, deleted)
                cursor += 1
        for sign in (-1, 1):
            assert events[cursor] == ('coupled', dimensions, experiment.INTERVAL, sign, np.pi/2)
            cursor += 1
    prediction = events[cursor]
    assert prediction[0] == 'prediction'
    assert prediction[1] == experiment.LEVELS[0]
    assert prediction[3:] == (-1, 0., False)
    assert events[cursor+1] == ('coupled', experiment.LEVELS[0], experiment.INTERVAL, -1, 0.)
    assert events[cursor+2:] == [('inverse', 1e-7), ('audit',)]
    assert len(executed['predictions']) == 9
    assert len(executed['coupled']) == 5


def test_source_histories_and_resources_survive_receipt_without_rewriting(executed):
    result = executed['result']
    for level_index, level in enumerate(result['levels']):
        assert level['tensor_dimension'] == int(np.prod(level['local_dimensions']))
        for phase in (0, 1):
            call = executed['coupled'][2*level_index+phase]
            assert call['rows'] == level['executed_histories_by_native_phase'][phase]
            assert canonical(call['metadata']) == canonical(level['execution_by_phase'][phase])
            laws = probabilities(call['rows'])
            assert laws.min() > 0
            assert laws.sum() == pytest.approx(1., abs=2e-9)
            resources = call['metadata']['worst_branch_resources']
            assert resources['native_pulses'] > 1_000_000
            assert resources['proper_code_service_events'] == 2*resources['native_pulses']
            assert call['metadata']['all_outcomes_retained']
            assert not call['metadata']['reset_between_reads']
            assert call['metadata']['steps_per_interval'] == experiment.STEPS
            assert call['metadata']['maximum_local_target_residual'] < 2e-10
            predicted = executed['predictions'][4*level_index+2*phase]
            deleted = executed['predictions'][4*level_index+2*phase+1]
            assert predicted['rows'] == level['predicted_histories_by_native_phase'][phase]
            assert deleted['rows'] == level['quadratic_deletion_predictions_by_native_phase'][phase]
    assert result['controls']['no_feedback_execution'] == executed['coupled'][-1]['rows']
    assert result['controls']['no_feedback_prediction'] == executed['predictions'][-1]['rows']


def test_native_joint_history_retains_each_quantum_record_and_phase(executed):
    for level in executed['result']['levels']:
        assert level['native_history_count'] == 1056
        for prefix in ('predicted', 'executed'):
            rows = level[f'native_{prefix}_history']
            assert len(rows) == 1056
            assert sum(row['probability'] for row in rows) == pytest.approx(1., abs=2e-9)
            assert len({(row['initial'], *row['records']) for row in rows}) == 1056
            histories = level[f'{prefix}_histories_by_native_phase']
            for row in rows:
                first_coin, first_current, first_bit, last_bit, last_coin, _ = row['records']
                assert first_coin in (0, 1) and last_coin in (0, 1)
                assert 0 <= row['initial'] < 66
                phase = int(first_current > 0)
                expected = histories[phase][2*first_bit+last_bit]['probability']/264
                assert row['probability'] == expected
        actual_tv = sum(abs(a['probability']-b['probability']) for a, b in zip(
            level['native_predicted_history'], level['native_executed_history'], strict=True))/2
        assert level['native_prediction_tv'] == actual_tv
        assert actual_tv <= max(level['prediction_history_tv_by_phase'])+2e-12


def test_recovered_local_action_and_geometry_match_independent_audit(executed):
    result = executed['result']
    recovered = result['recovered']
    # Direct cube coordinates are used only here, after the completed run.
    positions = np.asarray(list(product((0, 1), repeat=3)))[list(source.PERMUTATION)]
    squared_distances = np.sum((positions[:, None]-positions[None, :])**2, axis=-1)
    adjacency = (squared_distances == 1).astype(int)
    np.testing.assert_array_equal(recovered['adjacency'], adjacency)
    np.testing.assert_array_equal(recovered['graph_distances'], squared_distances)
    np.testing.assert_allclose(recovered['gram'], (positions-.5)@(positions-.5).T, atol=2e-13, rtol=0)
    np.testing.assert_allclose(recovered['pair_weights'], adjacency, atol=2e-8, rtol=0)
    np.testing.assert_allclose(recovered['onsite_quadratic'], 4., atol=2e-8, rtol=0)
    np.testing.assert_allclose(recovered['onsite_quartic'], 1/24, atol=2e-8, rtol=0)
    assert recovered['bare_mass_squared'] == pytest.approx(1., abs=2e-8)
    assert recovered['quartic'] == pytest.approx(1/24, abs=2e-8)
    assert recovered['bond_kappa'] == pytest.approx(1., abs=2e-8)
    assert recovered['spatial_dimension'] == 3
    assert recovered['bare_mass_is_not_measured_interacting_gap']
    assert recovered['positive_quartic_interaction_resolved']
    for key, truth in [('bare_mass_squared', 1.), ('quartic', 1/24), ('edge_weight', 1.)]:
        interval = recovered['homogeneous_intervals'][key]
        assert interval['lower']-2e-11 <= truth <= interval['upper']+2e-11
    assert recovered['homogeneous_intervals']['quartic']['lower'] > 0


def test_coupled_predictions_beat_quadratic_deletion_with_declared_error_budget(executed):
    result = executed['result']
    for level in result['levels']:
        budget = level['ideal_inferred_to_source_history_tv_upper']
        parameter = level['parameter_error_certificate']['complete_history_tv_upper']
        splitting = float(Fraction(level['ideal_splitting_history_tv_upper_exact']))
        assert budget == min(1., parameter+splitting)
        assert parameter > 0 and splitting > 0
        for index in (0, 1):
            expected = probabilities(level['predicted_histories_by_native_phase'][index])
            actual = probabilities(level['executed_histories_by_native_phase'][index])
            deleted = probabilities(level['quadratic_deletion_predictions_by_native_phase'][index])
            difference = float(np.sum(abs(expected-actual))/2)
            separation = float(np.sum(abs(deleted-actual))/2)
            assert level['prediction_history_tv_by_phase'][index] == pytest.approx(difference, abs=1e-18)
            assert level['quadratic_deletion_history_tv_by_phase'][index] == pytest.approx(separation, abs=1e-18)
            assert difference < 2e-8
            assert difference < budget+2e-8
            assert separation > budget+2e-8
        assert level['ideal_source_to_untruncated_tv_upper'] < .002
        assert level['ideal_inferred_continuous_to_untruncated_tv_upper'] < .002
        assert all(level['checks'].values()) and level['all_pass']
    controls = result['controls']
    free = controls['executed_free_onsite_parameters']
    assert abs(free['quartic']) <= free['quartic_radius']
    actual = probabilities(controls['no_feedback_execution'])
    prediction = probabilities(controls['no_feedback_prediction'])
    assert .5*np.sum(abs(actual-prediction)) < 2e-8
    original = probabilities(result['levels'][0]['executed_histories_by_native_phase'][0])
    assert .5*np.sum(abs(original-actual)) == pytest.approx(controls['feedback_deletion_history_tv'], abs=1e-18)
    assert controls['feedback_deletion_history_tv'] > 1e-3


def test_source_bindings_and_scope_flags_remain_explicit(executed):
    result = executed['result']
    assert set(result['source_binding_sha256']) == set(experiment.SOURCE_PATHS)
    for filename, digest in result['source_binding_sha256'].items():
        assert digest == hashlib.sha256((experiment.ROOT/filename).read_bytes()).hexdigest()
    protocol = result['protocol']
    assert protocol['record_laws_are_exact_probability_executions_not_samples']
    assert protocol['isolation_controls_and_clock_calibration_supplied']
    assert not protocol['inference_arithmetic_compiled_into_source']
    sensitivity = result['finite_sampling_sensitivity']
    assert not sensitivity['finite_sampling_was_executed']
    assert sensitivity['clock_laws_treated_as_exact']
    shots = sensitivity['independent_Bernoulli_budget']
    assert shots['probability_laws'] == 176
    assert shots['total_fresh_trials'] == 176*shots['trials_per_law']
    assert not shots['sampling_executed']
    assert not result['floating_or_physical_service_error_certified']
    assert all(result['checks'].values()) and result['all_pass']
    json.dumps(result, allow_nan=False)

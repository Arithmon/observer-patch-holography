"""Orchestration checks separate source generation from observed outcomes."""
from copy import deepcopy
from fractions import Fraction as F
import inspect
import math

import pytest

from . import executable as e


AXES = {(i, 0): .5j for i in range(11)} | {(0, j): .5 for j in range(11)}
KEYS = {(i, j) for i in range(1, 11) for j in range(1, 11)}


def fake_calibration(axis_means, *, extent):
    assert axis_means == AXES
    assert len(axis_means) == 21
    return {'parameters': {'alpha_real': .7, 'thermal_occupation': .1},
            'coordinate_extent': extent}


def fake_source(dimension, alpha, thermal, settings):
    assert dimension == 4
    assert alpha == .7 and thermal == .1
    records = []
    for setting in settings:
        # Deliberately unrelated to the Gaussian reference. Scoring must
        # nevertheless use these source-emitted values without replacement.
        mean = (.2+.01*setting['u'] if setting['theta'] == 0
                else -.1+.02*setting['v'])
        records.append({**setting, 'mean': mean,
                        'probabilities': [(1+mean)/2, (1-mean)/2],
                        'probability_mass': 1.,
                        'both_detector_outcomes_retained': True})
    return {'records': records, 'analytic_characteristic_used_for_prediction': False}


def fake_certificate(dimension, **parameters):
    assert dimension == 4
    assert all(isinstance(value, F) for value in parameters.values())
    assert set(parameters) == {'thermal', 'alpha_real', 'beta_real', 'beta_imag'}
    return {'probability_error_upper': F(1, 100)}


def install_fakes(monkeypatch):
    monkeypatch.setattr(e, 'calibrate_state', fake_calibration)
    monkeypatch.setattr(e, 'source_run', fake_source)
    monkeypatch.setattr(e, 'cutoff_certificate', fake_certificate)


@pytest.mark.parametrize('extent', [4, 5])
def test_all_one_hundred_predictions_are_emitted_by_both_source_readouts(monkeypatch, extent):
    install_fakes(monkeypatch)
    calls = []

    def traced_source(*args):
        calls.append(deepcopy(args))
        return fake_source(*args)

    monkeypatch.setattr(e, 'source_run', traced_source)
    prediction = e.predict(AXES, extent=extent, dimension=4)
    assert len(prediction['predictions']) == 100
    assert {tuple(point['index']) for point in prediction['predictions']} == KEYS
    assert len(calls) == 1 and len(calls[0][3]) == 200
    for point in prediction['predictions']:
        i, j = point['index']
        assert point['u'] == extent*i/10
        assert point['v'] == extent*j/10
        assert point['source_real'] == pytest.approx(.2+.01*extent*i/10)
        assert point['source_imag'] == pytest.approx(-.1+.02*extent*j/10)
        assert point['source_real_probabilities'][0] == pytest.approx((1+point['source_real'])/2)
        assert point['source_imag_probabilities'][1] == pytest.approx((1-point['source_imag'])/2)
        assert point['occupation_certificate']['probability_error_upper'] == '1/100'
    settings = calls[0][3]
    assert [setting['theta'] for setting in settings] == [0., -math.pi/2]*100
    assert set(inspect.signature(e.predict).parameters) == {'axis_means', 'extent', 'dimension'}


def test_scoring_uses_source_values_and_validation_cannot_change_prediction(monkeypatch):
    install_fakes(monkeypatch)
    prediction = e.predict(AXES, extent=4, dimension=4)
    before = deepcopy(prediction)
    observed = {tuple(row['index']): complex(row['source_real'], row['source_imag'])
                for row in prediction['predictions']}
    assert e.score(prediction, observed)['spin_mean_rmse'] == 0
    observed.update({key: 0j for key in observed})
    comparison = e.score(prediction, observed)
    direct = sum(row['source_real']**2+row['source_imag']**2
                 for row in prediction['predictions'])/200
    assert comparison['spin_mean_rmse'] == pytest.approx(math.sqrt(direct))
    assert comparison['zero_response_spin_mean_rmse'] == 0
    assert prediction == before
    assert e.predict(AXES, extent=4, dimension=4) == before


def example_prediction():
    return {'predictions': [
        {'index': list(key), 'source_real': .1, 'source_imag': .6,
         'source_real_probabilities': [.55, .45],
         'source_imag_probabilities': [.8, .2]}
        for key in sorted(KEYS)]}


def test_independent_binary_probability_and_spin_score():
    comparison = e.score(example_prediction(), {key: .3+.4j for key in KEYS})
    assert comparison['withheld_complex_points'] == 100
    assert comparison['withheld_real_statistics'] == 200
    assert comparison['spin_mean_rmse'] == pytest.approx(.2)
    assert comparison['binary_probability_rmse'] == pytest.approx(.1)
    assert comparison['maximum_component_residual'] == pytest.approx(.2)
    assert comparison['zero_response_spin_mean_rmse'] == pytest.approx(math.sqrt(.125))
    assert not comparison['parameter_uncertainty_propagated']
    assert not comparison['sampling_or_apparatus_error_certified']
    assert not comparison['statistical_significance_claimed']


@pytest.mark.parametrize('mutation', ['missing', 'extra'])
def test_all_observed_coordinates_are_required(mutation):
    observations = {key: 0j for key in KEYS}
    if mutation == 'missing':
        observations.pop((1, 1))
    else:
        observations[0, 0] = 1
    with pytest.raises(ValueError):
        e.score(example_prediction(), observations)


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'off_grid'])
def test_predictions_must_cover_exactly_the_same_one_hundred_coordinates(mutation):
    prediction = example_prediction()
    if mutation == 'missing':
        prediction['predictions'].pop()
    elif mutation == 'duplicate':
        prediction['predictions'][-1] = deepcopy(prediction['predictions'][0])
    else:
        prediction['predictions'][0]['index'] = [0, 0]
    with pytest.raises(ValueError):
        e.score(prediction, {key: 0j for key in KEYS})


@pytest.mark.parametrize('value', [float('nan'), float('inf'), 1.1])
def test_invalid_predicted_spin_component_is_not_a_scientific_comparison(value):
    prediction = example_prediction()
    prediction['predictions'][0]['source_real'] = value
    with pytest.raises(ValueError):
        e.score(prediction, {key: 0j for key in KEYS})


@pytest.mark.parametrize('value', [float('nan'), 1.1, -1.1j])
def test_invalid_observed_component_rejected(value):
    observations = {key: 0j for key in KEYS}
    observations[1, 1] = value
    with pytest.raises(ValueError):
        e.score(example_prediction(), observations)


@pytest.mark.parametrize('mutation', ['missing', 'order', 'phase', 'mass', 'mean', 'range', 'nan'])
def test_source_records_must_match_every_requested_setting(monkeypatch, mutation):
    install_fakes(monkeypatch)

    def corrupted(*args):
        result = fake_source(*args)
        records = result['records']
        if mutation == 'missing':
            records.pop()
        elif mutation == 'order':
            records[0], records[2] = records[2], records[0]
        elif mutation == 'phase':
            records[1]['theta'] = 0.
        elif mutation == 'mass':
            records[0]['probabilities'] = [.2, .2]
        elif mutation == 'mean':
            records[0]['mean'] = -.8
        elif mutation == 'range':
            records[0]['probabilities'] = [1.2, -.2]
        else:
            records[0]['probabilities'][0] = float('nan')
        return result

    monkeypatch.setattr(e, 'source_run', corrupted)
    with pytest.raises(ArithmeticError):
        e.predict(AXES, extent=4, dimension=4)


def test_run_retains_both_coordinate_extents_without_validation_selection(monkeypatch):
    # Use the production file/binding inventory, changing only the protocol's
    # small dimension in memory so integration does not repeat the long run.
    install_fakes(monkeypatch)
    observations = {'calibration': AXES, 'validation': {key: .3+.4j for key in KEYS},
                    'provenance': {'fixture': True}, 'orientation': {'fixture': True}}
    monkeypatch.setattr(e, 'load_observations', lambda: deepcopy(observations))
    original_loads = e.json.loads

    def protocol_loads(raw):
        result = original_loads(raw)
        result['occupation_dimension'] = 4
        return result

    monkeypatch.setattr(e.json, 'loads', protocol_loads)
    monkeypatch.setattr(e, '_audit', lambda prediction: {
        'maximum_probability_mass_residual': 0.,
        'maximum_full_phase_operation_residual': 0.,
        'cutoff_consistent_with_independent_numeric_reference': True})
    receipt = e.run()
    assert receipt['all_pass']
    assert [run['coordinate_extent'] for run in receipt['runs']] == [4, 5]
    assert all(run['comparison']['withheld_real_statistics'] == 200 for run in receipt['runs'])
    assert all(run['calibration']['parameters'] == {'alpha_real': .7, 'thermal_occupation': .1}
               for run in receipt['runs'])
    observations['validation'] = {key: -.7+.2j for key in KEYS}
    modified = e.run()
    assert [run['coordinate_extent'] for run in modified['runs']] == [4, 5]
    for before, after in zip(receipt['runs'], modified['runs'], strict=True):
        assert before['predictions'] == after['predictions']
        assert before['comparison'] != after['comparison']

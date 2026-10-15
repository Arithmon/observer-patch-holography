"""Exact per-run seals coexist with narrowly numerical replay tolerance."""
from copy import deepcopy
import hashlib
import json

import pytest

from .replay import (SEAL_KEY, compare_receipt, prediction_payload,
                     validate_prediction_seal)


def _seal(receipt):
    receipt[SEAL_KEY] = hashlib.sha256(json.dumps(
        prediction_payload(receipt), sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    return receipt


def _receipt():
    rows = {family: {'rows': [dict(family=family, time_ms=2., beta_imag=.2,
                                 component='imag', predicted_mean=.4,
                                 predicted_probabilities=[.7, .3],
                                 observed_mean=.3, observed_probability_zero=.65,
                                 residual=.1)]}
            for family in ('characteristic', 'position')}
    return _seal({
        'schema': 'oph.anharmonic_observations.results.v1',
        'protocol_sha256': 'a'*64, 'protocol': {'duration': .0002},
        'calibration': {'alpha_real': -.8},
        'source_binding_sha256': {'source.py': 'b'*64},
        'provenance': {'file': {'sha256': 'c'*64}},
        'zero_duration_position_predictions': [],
        'positive_duration_values_used': False,
        'runs': [{'name': 'primary', 'source_execution': {'native_pulses': 123,
                                                        'word_sha256': 'd'*64},
                  'comparison': rows}],
        'numerical_refinements': {'maximum_spin_mean_change': 1e-5},
        'independent_zero_duration_position_check': {'scalar_means': 13},
        'checks': {'test': True}, 'all_pass': True,
        'all_pass_meaning': 'software only'})


def test_reconstruction_removes_only_scoring_fields_and_preserves_row_order():
    receipt = _receipt()
    prediction = prediction_payload(receipt)
    assert prediction['schema'] == 'oph.anharmonic_observations.predictions.v1'
    assert SEAL_KEY not in prediction
    assert 'checks' not in prediction
    assert prediction['source_binding_sha256'] == receipt['source_binding_sha256']
    assert [row['family'] for row in prediction['runs'][0]['predictions']] == ['characteristic', 'position']
    assert all(not {'observed_mean', 'observed_probability_zero', 'residual'} & row.keys()
               for row in prediction['runs'][0]['predictions'])
    assert validate_prediction_seal(receipt) == receipt[SEAL_KEY]


def test_replay_accepts_last_bit_difference_only_after_own_seal_validation():
    expected = _receipt()
    actual = deepcopy(expected)
    row = actual['runs'][0]['comparison']['position']['rows'][0]
    row['predicted_mean'] += 1e-12
    with pytest.raises(ValueError, match='seal'):
        compare_receipt(expected, actual)
    _seal(actual)
    assert expected[SEAL_KEY] != actual[SEAL_KEY]
    compare_receipt(expected, actual)


@pytest.mark.parametrize('side', ['expected', 'actual'])
def test_both_receipts_must_have_a_valid_pre_scoring_seal(side):
    receipts = {'expected': _receipt(), 'actual': _receipt()}
    receipts[side][SEAL_KEY] = '0'*64
    with pytest.raises(ValueError, match='seal'):
        compare_receipt(receipts['expected'], receipts['actual'])


def test_valid_new_seal_does_not_excuse_a_changed_prediction():
    expected = _receipt()
    actual = deepcopy(expected)
    actual['runs'][0]['comparison']['position']['rows'][0]['predicted_mean'] += .001
    _seal(actual)
    with pytest.raises(ValueError, match='Floating replay'):
        compare_receipt(expected, actual)


@pytest.mark.parametrize('field', ['source_hash', 'word_hash', 'integer_count', 'boolean'])
def test_valid_new_seal_does_not_excuse_discrete_binding_changes(field):
    expected = _receipt()
    actual = deepcopy(expected)
    if field == 'source_hash':
        actual['source_binding_sha256']['source.py'] = 'e'*64
    elif field == 'word_hash':
        actual['runs'][0]['source_execution']['word_sha256'] = 'e'*64
    elif field == 'integer_count':
        actual['runs'][0]['source_execution']['native_pulses'] += 1
    else:
        actual['positive_duration_values_used'] = True
    _seal(actual)
    with pytest.raises(ValueError, match='Receipt mismatch'):
        compare_receipt(expected, actual)


def test_observation_and_scoring_changes_are_still_compared():
    expected = _receipt()
    actual = deepcopy(expected)
    actual['runs'][0]['comparison']['position']['rows'][0]['observed_mean'] += .01
    assert validate_prediction_seal(actual) == actual[SEAL_KEY]
    with pytest.raises(ValueError, match='Floating replay'):
        compare_receipt(expected, actual)


def test_only_top_level_seal_is_exempted_from_exact_string_replay():
    expected = _receipt()
    expected['runs'][0]['source_execution'][SEAL_KEY] = 'nested declaration'
    _seal(expected)
    actual = deepcopy(expected)
    actual['runs'][0]['source_execution'][SEAL_KEY] = 'different declaration'
    _seal(actual)
    with pytest.raises(ValueError, match='Receipt mismatch'):
        compare_receipt(expected, actual)


def test_nonfinite_prediction_cannot_be_sealed_or_replayed():
    expected = _receipt()
    actual = deepcopy(expected)
    actual['runs'][0]['comparison']['position']['rows'][0]['predicted_mean'] = float('nan')
    with pytest.raises(ValueError):
        compare_receipt(expected, actual)

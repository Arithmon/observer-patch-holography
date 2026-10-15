"""Validate each prediction seal before ordinary numerical receipt replay.

The seal binds one exact pre-scoring floating payload. Different valid runs
may have different seals because their numerical diagnostics differ in the
last bits. Each seal is independently checked; only that top-level seal is
then removed for the numerical comparison. Source word hashes, protocol/data
bindings, integer counts, and all other strings remain exact.

Fresh Linux and Windows executions of the macOS receipt differ from it by
about 1.6e-10 on one runner and by more than 2e-10 at the same derived
position on another, so floats compare within ``PORTABLE_REPLAY_TOLERANCE``.
It is not a physical or floating-point error certificate.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math


PORTABLE_REPLAY_TOLERANCE = 2e-9
SEAL_KEY = 'prediction_payload_sha256_before_scoring'
SCORING_FIELDS = frozenset({SEAL_KEY, 'numerical_refinements',
                           'independent_zero_duration_position_check', 'checks',
                           'all_pass_meaning', 'all_pass'})
OBSERVATION_FIELDS = frozenset({'observed_mean', 'observed_probability_zero', 'residual'})


def prediction_payload(receipt):
    """Recover the exact prediction object from its scored, ordered rows.

    Prediction generation fixes the row order to characteristic then position;
    scoring preserves the within-family order. No observed field belongs to
    the recovered prediction object.
    """
    if receipt.get('schema') != 'oph.anharmonic_observations.results.v1':
        raise ValueError('Expected an anharmonic result receipt')
    result = {key: deepcopy(value) for key, value in receipt.items()
              if key not in SCORING_FIELDS and key != 'runs'}
    result['schema'] = 'oph.anharmonic_observations.predictions.v1'
    result['runs'] = []
    try:
        for run in receipt['runs']:
            recovered = {key: deepcopy(value) for key, value in run.items() if key != 'comparison'}
            recovered['predictions'] = [
                {key: deepcopy(value) for key, value in row.items() if key not in OBSERVATION_FIELDS}
                for family in ('characteristic', 'position')
                for row in run['comparison'][family]['rows']]
            result['runs'].append(recovered)
    except (KeyError, TypeError) as error:
        raise ValueError('Cannot reconstruct complete pre-scoring prediction payload') from error
    return result


def _compare(expected, actual, path='root'):
    if type(expected) is not type(actual):
        raise ValueError(f'Receipt type mismatch at {path}')
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError(f'Receipt keys mismatch at {path}')
        for key in expected:
            _compare(expected[key], actual[key], f'{path}.{key}')
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f'Receipt length mismatch at {path}')
        for i, (a, b) in enumerate(zip(expected, actual)):
            _compare(a, b, f'{path}[{i}]')
    elif isinstance(expected, float):
        if (not math.isfinite(expected) or not math.isfinite(actual)
                or abs(expected-actual) > PORTABLE_REPLAY_TOLERANCE):
            raise ValueError(f'Floating replay mismatch at {path}')
    elif expected != actual:
        raise ValueError(f'Receipt mismatch at {path}')


def validate_prediction_seal(receipt):
    """Check the stored pre-score seal against this receipt's own predictions."""
    payload = json.dumps(prediction_payload(receipt), sort_keys=True,
                         separators=(',', ':'), allow_nan=False).encode()
    digest = hashlib.sha256(payload).hexdigest()
    if receipt.get(SEAL_KEY) != digest:
        raise ValueError('Pre-scoring prediction seal does not match this receipt')
    return digest


def compare_receipt(expected, actual):
    """Self-validate both exact seals, then compare all remaining receipt data."""
    validate_prediction_seal(expected)
    validate_prediction_seal(actual)
    expected = {key: value for key, value in expected.items() if key != SEAL_KEY}
    actual = {key: value for key, value in actual.items() if key != SEAL_KEY}
    _compare(expected, actual)

"""Validate each prediction seal before ordinary numerical receipt replay.

The seal binds one exact pre-scoring floating payload. Different valid runs
may have different seals because their numerical diagnostics differ in the
last bits. Each seal is independently checked; only that top-level seal is
then removed for the existing 2e-10 numerical comparison. Source word hashes,
protocol/data bindings, integer counts, and all other strings remain exact.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json

from ..observer.scalar_response.experiment import compare_receipt as _compare


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

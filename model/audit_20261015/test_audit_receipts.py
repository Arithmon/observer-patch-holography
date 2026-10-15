"""Replay audit evidence while keeping provenance and structures exact."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import pytest

from model.audit_20261015 import architecture_counterexamples, empirical_probe


FOLDER = Path(__file__).resolve().parent
MODEL = FOLDER.parent
ABSOLUTE_TOLERANCE = 2e-9


def _compare(actual, expected, path='root'):
    """Only finite floats may vary; record types and all other values are exact."""
    assert type(actual) is type(expected), (
        f'{path}: type {type(actual).__name__} != {type(expected).__name__}')
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), f'{path}: dictionary keys changed'
        for key in expected:
            _compare(actual[key], expected[key], f'{path}.{key}')
    elif isinstance(expected, list):
        assert len(actual) == len(expected), f'{path}: list length changed'
        for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
            _compare(left, right, f'{path}[{index}]')
    elif isinstance(expected, float):
        assert math.isfinite(actual) and math.isfinite(expected), (
            f'{path}: receipt floats must be finite')
        assert math.isclose(actual, expected, rel_tol=0.,
                            abs_tol=ABSOLUTE_TOLERANCE), (
            f'{path}: {actual!r} != {expected!r}')
    else:
        assert actual == expected, f'{path}: exact value changed'


def _read(name):
    return json.loads((FOLDER/name).read_text(encoding='utf-8'))


def _receipt(value):
    # Compare the actual JSON receipt schema, not Python implementation
    # scalar classes: the empirical producer retains some numpy.float64
    # intermediates which its public JSON output serializes as JSON floats.
    # Nonfinite values are rejected before any tolerance comparison.
    return json.loads(json.dumps(value, allow_nan=False))


def test_architecture_receipt_replays_exact_controls_and_reference_contrasts():
    _compare(_receipt(architecture_counterexamples.run()), _read('architecture_results.json'))


def test_empirical_receipt_binds_probe_and_all_external_inputs():
    paths = {
        'audit_probe': FOLDER/'empirical_probe.py',
        'transfer_data': MODEL/'anharmonic_transfer/data/asymmetric_wells_expectation_x.csv',
        'transfer_receipt': MODEL/'anharmonic_transfer/results.json',
    }
    expected = _read('empirical_probe.json')['input_sha256']
    assert expected.keys() == paths.keys()
    actual = {name: hashlib.sha256(path.read_bytes()).hexdigest()
              for name, path in paths.items()}
    _compare(actual, expected, 'input_sha256')


def test_empirical_receipt_replays_independent_lindblad_comparison():
    _compare(_receipt(empirical_probe.run()), _read('empirical_probe.json'))


@pytest.mark.parametrize(('actual', 'expected'), [
    ({'flag': 1}, {'flag': True}),
    ({'count': 3.0}, {'count': 3}),
    ({'hash': 'a'*63+'b'}, {'hash': 'a'*64}),
    ({'words': ['one', 'changed']}, {'words': ['one', 'two']}),
    ({'number': float('nan')}, {'number': float('nan')}),
    ({'number': float('inf')}, {'number': float('inf')}),
    ({'number': 1.}, {'number': float('-inf')}),
    ({'number': 1.+4e-9}, {'number': 1.}),
    ({'missing': 1}, {'missing': 1, 'required': 2}),
    ([1, 2], [1]),
])
def test_replay_rejects_corruption_and_nonfinite_values(actual, expected):
    with pytest.raises(AssertionError):
        _compare(actual, expected)


def test_replay_allows_only_declared_float_roundoff():
    expected = {'probability': .25, 'dimension': 32, 'normalized': True,
                'input_sha256': {'source': 'a'*64}}
    actual = deepcopy(expected)
    actual['probability'] += 1e-9
    _compare(actual, expected)

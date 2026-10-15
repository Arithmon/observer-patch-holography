"""Preparation transfer and the finite readout's principal inverse."""
from copy import deepcopy
import gzip
import hashlib
import json
import math

import pytest

from . import calibration, data, experiment
from .test_data import synthetic_csv


@pytest.mark.parametrize('alpha', [-1.7, -.9, 0., .9, 1.7])
@pytest.mark.parametrize('thermal', [0., .13, .7])
def test_finite_estimator_inverse_recovers_principal_preparation(alpha, thermal):
    measured = calibration.estimator(alpha, thermal)
    assert calibration.zero_time_displacement(measured, thermal) == pytest.approx(alpha, abs=2e-15)


def test_finite_estimator_is_not_exact_position_expectation():
    alpha, thermal = .96, .13
    measured = calibration.estimator(alpha, thermal)
    assert abs(measured/math.sqrt(2)-alpha) > .1
    assert calibration.zero_time_displacement(measured, thermal) == pytest.approx(alpha)


def test_principal_branch_does_not_claim_global_identifiability():
    alpha = 2.4
    reconstructed = calibration.zero_time_displacement(calibration.estimator(alpha, .1), .1)
    assert abs(reconstructed) <= math.pi/(4*.4)
    assert reconstructed != pytest.approx(alpha)
    assert calibration.estimator(reconstructed, .1) == pytest.approx(calibration.estimator(alpha, .1))


@pytest.mark.parametrize('position,thermal,h', [(3., .1, .4), (0., -.1, .4), (0., .1, 0.),
                                               (True, .1, .4), (float('nan'), .1, .4),
                                               (0., float('inf'), .4), (0., .1, -.4)])
def test_invalid_preparation_domain_rejected(position, thermal, h):
    with pytest.raises(ValueError):
        calibration.zero_time_displacement(position, thermal, h)


def test_zero_time_fit_ignores_poisoned_later_responses(tmp_path):
    path = synthetic_csv(tmp_path/'poison.csv', later='LATER DATA FORBIDDEN')
    estimates = {case: calibration.zero_time_displacement(value, .1)
                 for case, value in data.parse_zero(path).items()}
    assert len(estimates) == 5


@pytest.mark.parametrize('compressed', [False, True])
def test_prior_receipt_uses_exact_digest_and_only_transfers_parameters(tmp_path, monkeypatch, compressed):
    parameters = {'alpha_real': -.91, 'thermal_occupation': .12}
    receipt = {'calibration': {'parameters': parameters}, 'later_outcomes': 'UNUSED'}
    raw = json.dumps(receipt).encode()
    name = 'old.json.gz' if compressed else 'old.json'
    pack = (lambda value: gzip.compress(value, mtime=0)) if compressed else (lambda value: value)
    (tmp_path/name).write_bytes(pack(raw))
    protocol = {'preparation': {'prior_receipt': name,
                               'prior_receipt_sha256': hashlib.sha256(raw).hexdigest()}}
    monkeypatch.setattr(calibration, 'ROOT', tmp_path)
    assert calibration.transferred_parameters(protocol) == parameters
    (tmp_path/name).write_bytes(pack(raw+b' '))
    with pytest.raises(ValueError, match='changed'):
        calibration.transferred_parameters(protocol)


@pytest.mark.parametrize('corruption', ['header', 'truncated', 'payload'])
def test_corrupt_frozen_compression_is_rejected(tmp_path, monkeypatch, corruption):
    raw = b'{"calibration": {"parameters": {"alpha_real": -0.91}}}'
    compressed = bytearray(gzip.compress(raw, mtime=0))
    if corruption == 'header':
        compressed[0] ^= 1
    elif corruption == 'truncated':
        compressed = compressed[:-4]
    else:
        compressed[12] ^= 0xff
    (tmp_path/'old.json.gz').write_bytes(compressed)
    monkeypatch.setattr(calibration, 'ROOT', tmp_path)
    protocol = {'preparation': {'prior_receipt': 'old.json.gz',
                               'prior_receipt_sha256': hashlib.sha256(raw).hexdigest()}}
    with pytest.raises(ValueError):
        calibration.transferred_parameters(protocol)


def test_transfer_uses_frozen_input_even_when_live_receipt_changes(tmp_path, monkeypatch):
    protocol = json.loads((experiment.FOLDER/'protocol.json').read_text())
    relative = protocol['preparation']['prior_receipt']
    target = tmp_path/relative
    target.parent.mkdir(parents=True)
    target.write_bytes((calibration.ROOT/relative).read_bytes())
    expected = calibration.transferred_parameters(protocol)
    live = tmp_path/'model/anharmonic_observations/results.json'
    live.parent.mkdir(parents=True)
    live.write_text('A later live source receipt must not recalibrate the frozen transfer.')
    monkeypatch.setattr(calibration, 'ROOT', tmp_path)
    assert calibration.transferred_parameters(protocol) == expected


def test_actual_frozen_prior_digest_is_supported_without_reading_new_outcomes():
    protocol = json.loads((experiment.FOLDER/'protocol.json').read_text())
    transferred = calibration.transferred_parameters(protocol)
    assert abs(transferred['alpha_real']) > 0
    assert transferred['thermal_occupation'] >= 0
    altered = deepcopy(protocol)
    altered['preparation']['prior_receipt_sha256'] = '0'*64
    with pytest.raises(ValueError):
        experiment._configuration(altered)

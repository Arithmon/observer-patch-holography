"""Exact residual floors, common apparatus identity and complete data coverage."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math
from pathlib import Path
import shutil

import pytest

from . import certified_comparison as c, data, experiment as e, replay
from .calibration import transferred_parameters


@pytest.fixture
def inputs():
    protocol = json.loads((e.FOLDER/'protocol.json').read_text())
    return protocol, transferred_parameters(protocol), e._configuration(protocol)


@pytest.fixture
def attachment(inputs):
    return c.certified_inputs(*inputs)


@pytest.mark.parametrize('observed', [F(-5), F(-1), F(0), F(2), F(7, 3)])
@pytest.mark.parametrize('radius', [F(0), F(1, 10), F(2)])
def test_allowance_is_attained_and_bounds_every_compatible_model(observed, radius):
    centre = F(1, 3)
    floor = c.unavoidable_allowance(observed, centre, radius)
    lower, upper = centre-radius, centre+radius
    nearest = min(upper, max(lower, observed))
    assert floor == abs(observed-nearest)
    for index in range(11):
        truth = lower+F(index, 10)*(upper-lower)
        assert abs(observed-truth) >= floor
    assert c.unavoidable_allowance(observed, centre, radius+1) <= floor


@pytest.mark.parametrize('values', [(True, F(0), F(1)), (F(0), 0., F(1)),
                                    (F(0), F(0), -F(1)), (F(0), F(0), float('nan'))])
def test_allowance_rejects_inexact_or_negative_bounds(values):
    with pytest.raises(ValueError):
        c.unavoidable_allowance(*values)


@pytest.mark.parametrize('token', [None, True, .1, 'nan', 'inf', '-inf', '1/3', 'not a number'])
def test_deposited_values_require_finite_decimal_tokens(token):
    with pytest.raises(ValueError):
        c.decimal_token(token)


def test_decimal_custody_does_not_substitute_float_rounding():
    assert c.decimal_token('0.047169811') == F(47169811, 10**9)
    assert c.decimal_token('0.047169811') != F(float('0.047169811'))
    values = data.load_validation()
    assert all(float(c.decimal_token(row['observed_position_estimator_decimal'])) ==
               row['observed_position_estimator'] for row in values)


@pytest.mark.parametrize('section, field', [('prior', 'alpha_real'), ('prior', 'thermal_occupation'),
                                          ('apparatus', 'delta'), ('apparatus', 'omega'),
                                          ('apparatus', 'pulse_duration'), ('apparatus', 'theta'),
                                          ('apparatus', 'dephasing_rate'), ('apparatus', 'omega_carrier'),
                                          ('readout', 'h')])
def test_certificate_cannot_attach_to_a_different_preparation_or_apparatus(inputs, section, field):
    protocol, prior, apparatus = inputs
    target = prior if section == 'prior' else apparatus if section == 'apparatus' else protocol['readout']
    target[field] = math.nextafter(target[field], math.inf)
    with pytest.raises(ValueError, match='parameters differ'):
        c.certified_inputs(protocol, prior, apparatus)


@pytest.mark.parametrize('name', list(c.CERTIFICATE_PINS))
def test_changed_certificate_or_proof_source_is_not_silently_rebound(inputs, tmp_path, name):
    for relative in c.CERTIFICATE_PINS:
        target = tmp_path/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(c.ROOT/relative, target)
    target = tmp_path/name
    target.write_bytes(target.read_bytes()+b' ')
    with pytest.raises(ValueError, match='source or artifact changed'):
        c.certified_inputs(*inputs, root=tmp_path)


def test_every_held_out_row_and_both_refinements_remain_visible(attachment):
    result = c.compare(attachment, data.load_validation(), data.load_zero_tokens())
    assert [run['dimension'] for run in result['runs']] == [32, 64]
    for run, certified in zip(result['runs'], (2, 3), strict=True):
        assert len(run['held_out_rows']) == run['held_out_count'] == 125
        assert run['certified_overlap_count'] == certified
        assert run['unavailable_count'] == 125-certified
        assert {row['case'] for row in run['held_out_rows']} == set('abcde')
        assert run['certificate_endpoints_without_observations'] == [{
            'case': 'a', 'time_ms': '5',
            'reason': 'no observation at this exact published coordinate; no interpolation'}]
        assert len(run['zero_time_diagnostics_not_held_out']) == 1
    low, high = result['runs']
    assert low['minimum_common_uncovered_allowance_on_certified_held_out_rows'] == '0'
    assert F(high['minimum_common_uncovered_allowance_on_certified_held_out_rows']) > F('0.17431748')
    available = {F(row['time_ms']): row for row in high['held_out_rows'] if row['status'] == 'certified_interval'}
    assert set(available) == {F('0.4'), F(2), F(10)}
    assert F(available[F('0.4')]['minimum_uncovered_absolute_allowance']) > F('0.04690045')
    assert available[F(10)]['minimum_uncovered_absolute_allowance'] == '0'
    assert F(high['zero_time_diagnostics_not_held_out'][0]['minimum_uncovered_absolute_allowance']) > F('0.34054967')
    assert result['experimental_uncertainty_applied'] is None
    assert result['accepted_counts_or_raw_trials_authenticated'] is False
    assert result['statistical_rejection_claimed'] is False
    assert result['full_uncertainty_agreement_established'] is False


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'foreign_time', 'token_mismatch',
                                    'nonfinite', 'boolean', 'missing_zero'])
def test_comparison_rejects_invalid_data_or_silent_selection(attachment, mutation):
    observations, zero = data.load_validation(), data.load_zero_tokens()
    if mutation == 'missing':
        observations.pop()
    elif mutation == 'duplicate':
        observations.append(observations[0])
    elif mutation == 'foreign_time':
        observations[0]['time_ms'] = .41
    elif mutation == 'token_mismatch':
        observations[0]['observed_position_estimator_decimal'] = '0'
    elif mutation == 'nonfinite':
        observations[0]['observed_position_estimator'] = float('nan')
    elif mutation == 'boolean':
        observations[0]['observed_position_estimator'] = True
    else:
        zero.pop('e')
    with pytest.raises(ValueError):
        c.compare(attachment, observations, zero)


def test_forged_radius_is_rejected_before_later_data_are_opened(monkeypatch):
    from .test_experiment import build_prediction
    prediction, _ = build_prediction(monkeypatch)
    prediction['certified_unrestricted_apparatus']['runs'][1]['records'][2]['radius_exact_binary'] = '1000'
    with pytest.raises(ValueError, match='Sealed occupation certificate'):
        e.evaluate(prediction)


@pytest.mark.parametrize('mutation', ['floor', 'omit_unavailable', 'claim_significance'])
def test_receipt_cannot_rewrite_its_exact_empirical_consequences(monkeypatch, mutation):
    from .test_experiment import build_prediction
    prediction, _ = build_prediction(monkeypatch)
    monkeypatch.setattr(e, 'load_validation', data.load_validation)
    result = e.evaluate(prediction)
    changed = deepcopy(result)
    comparison = changed['certified_unrestricted_apparatus_comparison']
    if mutation == 'floor':
        row = next(r for r in comparison['runs'][1]['held_out_rows'] if r['time_ms'] == '2' and r['case'] == 'a')
        row['minimum_uncovered_absolute_allowance'] = '0'
    elif mutation == 'omit_unavailable':
        comparison['runs'][0]['held_out_rows'].pop()
    else:
        comparison['statistical_rejection_claimed'] = True
    with pytest.raises(ValueError, match='Derived receipt'):
        replay.compare_receipt(changed, changed)

"""Adversarial counts, complete records and exact predeclared error thresholds."""
from copy import deepcopy
from fractions import Fraction as F
import json

import pytest

from . import protocol as p
from .__main__ import read_json


@pytest.fixture(scope='module')
def plan():
    return p.build_protocol()


@pytest.fixture(scope='module')
def c_counts(plan):
    return p.synthetic_counts('C', 1078, plan)


def test_fixed_plan_has_exact_simultaneous_error_bound(plan):
    certificate = plan['concentration']
    assert certificate['simultaneous_statistics'] == 56
    assert F(certificate['exponential_argument']) == 16
    assert 0 < F(certificate['familywise_error_upper']) < F(1, 1000)
    # Sum of 33 positive Taylor terms is below e^16; multiplying by
    # exp(-16) bounds every statistic simultaneously, without data fitting.
    assert F(certificate['exp_lower_partial_sum']) > 112_000


@pytest.mark.parametrize('candidate', ('C', 'E'))
def test_supplied_complete_software_data_separate_candidates(plan, candidate):
    data = p.synthetic_counts(candidate, 1091, plan)
    result = p.evaluate(data, plan)
    assert result['status'] == 'comparison_complete'
    assert result['candidate_verdicts'][candidate]['verdict'] == 'not_rejected'
    other = 'E' if candidate == 'C' else 'C'
    assert result['candidate_verdicts'][other]['verdict'] == 'rejected_under_contract'
    assert result['empirical_confirmation'] is False
    assert result['provenance_verified'] is False


def test_no_data_and_unfilled_template_never_validate_physics(plan):
    result = p.evaluate(None, plan)
    assert result['status'] == 'no_data'
    assert result['candidate_verdicts'] == {}
    assert result['empirical_confirmation'] is False
    with pytest.raises(ValueError, match='integer count'):
        p.evaluate(p.data_template(plan), plan)


@pytest.mark.parametrize('bad', (-1, 0.0, True, None, '0', float('nan')))
def test_noninteger_or_invalid_counts_are_rejected(plan, c_counts, bad):
    data = deepcopy(c_counts)
    data['settings'][0]['events'][0]['count'] = bad
    with pytest.raises(ValueError, match='integer count'):
        p.evaluate(data, plan)


def test_missing_trials_cannot_be_renormalized(plan, c_counts):
    data = deepcopy(c_counts)
    event = next(event for event in data['settings'][0]['events'] if event['count'])
    event['count'] -= 1
    with pytest.raises(ValueError, match='fixed trial count'):
        p.evaluate(data, plan)


def test_aborts_are_retained_and_budget_failure_is_separate(plan, c_counts):
    data = deepcopy(c_counts)
    for row in data['settings']:
        row['aborts'] = p.TRIALS
        for event in row['events']:
            event['count'] = 0
    result = p.evaluate(data, plan)
    assert result['status'] == 'apparatus_budget_exceeded'
    assert set(result['apparatus_budget_failures']) == set(p.SETTINGS)
    assert result['empirical_confirmation'] is False


def test_joint_record_tampering_is_detected_even_if_detector_fringe_unchanged(plan, c_counts):
    data = deepcopy(c_counts)
    for row in data['settings']:
        plus = sum(e['count'] for e in row['events'] if e['outcome'] == 1)
        for event in row['events']:
            event['count'] = 0
            # Force all events into impossible old-current +1 records.
            if event['coin'] == 0 and event['current'] == 1:
                event['count'] = plus if event['outcome'] == 1 else p.TRIALS-plus
    result = p.evaluate(data, plan)
    assert result['candidate_verdicts']['C']['verdict'] == 'rejected_under_contract'
    assert any(row['statistic'] not in ('plus', 'abort')
               for row in result['candidate_verdicts']['C']['separating_statistics'])


@pytest.mark.parametrize('mutation', ('duplicate_setting', 'duplicate_event', 'missing_setting',
                                      'missing_event', 'different_protocol', 'pretend_verified'))
def test_invalid_data_contracts(plan, c_counts, mutation):
    data = deepcopy(c_counts)
    if mutation == 'duplicate_setting':
        data['settings'][1] = deepcopy(data['settings'][0])
    elif mutation == 'duplicate_event':
        data['settings'][0]['events'][1] = deepcopy(data['settings'][0]['events'][0])
    elif mutation == 'missing_setting':
        data['settings'].pop()
    elif mutation == 'missing_event':
        data['settings'][0]['events'].pop()
    elif mutation == 'different_protocol':
        data['protocol_sha256'] = '0'*64
    else:
        data['data_kind'] = 'verified_experiment'
    with pytest.raises(ValueError):
        p.evaluate(data, plan)


def test_changing_prediction_after_observing_data_is_rejected(plan, c_counts):
    modified = deepcopy(plan)
    modified['predictions']['C']['phase_0']['plus'] = '3/8'
    with pytest.raises(ValueError, match='fixed source-bound'):
        p.evaluate(c_counts, modified)


def test_external_label_does_not_authenticate_records(plan, c_counts):
    data = deepcopy(c_counts)
    data['data_kind'] = 'external_unverified'
    data['provenance'] = 'An externally supplied acquisition identifier.'
    result = p.evaluate(data, plan)
    assert result['data_kind'] == 'external_unverified'
    assert result['candidate_verdicts']['C']['verdict'] == 'not_rejected'
    assert result['empirical_confirmation'] is False
    assert result['provenance_verified'] is False


@pytest.mark.parametrize('text', ('{"a": 1, "a": 2}', '{"a": NaN}', '{"a": Infinity}'))
def test_cli_refuses_ambiguous_or_nonfinite_json(tmp_path, text):
    path = tmp_path/'counts.json'
    path.write_text(text)
    with pytest.raises(ValueError):
        read_json(path)

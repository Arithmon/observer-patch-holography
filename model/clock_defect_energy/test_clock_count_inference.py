"""Independent binary-loss, marginal-control, and source-custody checks."""
from collections import defaultdict
from fractions import Fraction as F
import hashlib

import pytest

from tools import strict_json
from . import clock_count_inference as inference, joint_clock


@pytest.mark.parametrize('law', ('C', 'E'))
def test_exact_event_prior_and_both_readouts_have_real_information(law):
    rows = joint_clock.exact_clock_rows(law)
    prior = defaultdict(lambda: (F(), F(), F()))
    for (_, _, _, _, n), coefficients in rows.items():
        prior[n] = tuple(a+b for a,b in zip(prior[n], coefficients, strict=True))
    assert prior == {1: (F(307,360), F(), F()), 2: (F(53,360), F(), F())}
    result = inference.exact_law(law)
    risks = {name: tuple(map(F, pair)) for name,pair in
             result['optimal_squared_error_risk_count_units'].items()}
    assert risks['joint'][1] < risks['clock_only'][0]
    assert risks['clock_only'][1] < risks['field_only'][0]
    assert risks['field_only'][1] < risks['prior'][0]
    assert F(27, 10000) < risks['joint_gain_over_clock'][0] < F(28, 10000)
    # The same current has opposite strength at the two field signs, so the
    # joint readout is more than a renamed clock output.
    odds = {key: tuple(map(F, value)) for key,value in
            result['posterior_probability_N2_by_readout'].items()}
    assert odds['-1,-1'][0] > F(87, 100)
    assert odds['-1,1'][1] < F(62, 100)
    assert odds['1,1'] == (0, 0)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_matched_marginal_field_resampling_removes_extra_information(law):
    result = inference.exact_law(law)
    risks = result['optimal_squared_error_risk_count_units']
    assert result['matched_marginal_independent_field_control_risk'] == risks['clock_only']
    # Independent field resampling keeps both one-variable marginals exactly.
    rows = joint_clock.exact_clock_rows(law)
    count_mass = defaultdict(F)
    for (_, _, _, _, n), coefficients in rows.items():
        count_mass[n] += coefficients[0]
    plus = sum((mass*joint_clock.field_probabilities(n)[0]
                for n,mass in count_mass.items()), F())
    assert plus == F(307,360)*F(21,22)+F(53,360)*F(91,110)
    assert 0 < plus < 1
    assert count_mass == {1:F(307,360), 2:F(53,360)}
    # In the control P(field sign | N,current) = P(field sign), so its
    # posterior P(N | current,field sign) equals P(N | current).
    assert (F(1)-plus)*plus > 0


def test_retained_source_receipts_and_complete_failure_alphabet():
    for name in ('joint_clock_results.json','ramsey_service_results.json',
                 '../observer_time_fields/wave_reference_results.json'):
        inference._pred(name)
    result = strict_json.load(inference.OUTPUT)
    assert result['source_sha256'] == inference.source_bindings(result['source_sha256'])
    assert result['physical_seconds_calibrated'] is False
    assert result['contract_sha256'] == hashlib.sha256(inference.CAMPAIGN.read_bytes()).hexdigest()
    for row in result['compiled']:
        assert row['source_clock_export_abort_slots'] > 0
        assert row['source_read_slots_including_both_field_signs'] > 0
        assert set(row['finite_preparation_flag_probabilities_diagnostic_intervals']) == {
            'herald_failure','warmup_timeout','completion_timeout'}
        assert F(row['source_complete_mass_diagnostic_interval'][0]) <= 1 <= F(row['source_complete_mass_diagnostic_interval'][1])
    for row in result['exact']:
        assert F(row['source_robust_joint_gain_lower']) > F(27,10000)


def test_contract_predecessor_result_hashes():
    contract = strict_json.load(inference.CAMPAIGN)
    for key, name in (('joint_result_sha256','joint_clock_results.json'),
                      ('ramsey_result_sha256','ramsey_service_results.json'),
                      ('wave_reference_result_sha256','../observer_time_fields/wave_reference_results.json')):
        assert contract[key] == hashlib.sha256((inference.HERE/name).read_bytes()).hexdigest()

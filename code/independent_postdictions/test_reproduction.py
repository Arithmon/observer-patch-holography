"""Independent replay and hostile controls; no synthetic physical evidence."""
import ast
from copy import deepcopy
from decimal import Decimal
import inspect
import json
from pathlib import Path
import subprocess
import sys

import pytest

from independent_postdictions import admission, build, comparison, replay, verify


@pytest.fixture(scope='module')
def packet():
    return verify.load(verify.HERE/'receipt.json')


@pytest.fixture(scope='module')
def fresh():
    return build.calculations()


def test_committed_packet_reconstructs_every_equation(packet):
    verify.verify()


def test_fresh_independent_producer_matches_receipt(packet, fresh):
    # Compare the entire output, not a success flag or selected favorable row.
    verify.same(fresh, packet['calculations'])


def test_implementation_independence():
    producer = ast.parse(inspect.getsource(replay))
    checker = ast.parse(inspect.getsource(verify))
    for tree in (producer, checker):
        imports = [n.module or '' for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
        imports += [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
        assert not any('P_derivation' in n or 'particles' in n for n in imports)
    imports = [n.module for n in ast.walk(checker) if isinstance(n, ast.ImportFrom)]
    assert not {'replay', 'build', 'comparison'} & set(imports)
    assert not any(isinstance(n, ast.Call) and
                   ((isinstance(n.func, ast.Name) and n.func.id == 'open') or
                    (isinstance(n.func, ast.Attribute) and n.func.attr in ('read_text', 'read_bytes', 'load')))
                   for n in ast.walk(producer))
    assert set(inspect.signature(replay.tau_replay).parameters) == {'electron', 'muon', 'ordered'}


@pytest.mark.parametrize('field', ['P', 'alpha_inverse', 'alpha_U', 'mZ_over_v',
                                   'anchor_inverse', 'lepton_transport', 'unscreened_quark_transport'])
def test_false_alpha_numbers_are_rejected_even_with_valid_other_fields(packet, field):
    row = deepcopy(packet['calculations']['alpha'])
    row['rows'][0][field] = str(Decimal(row['rows'][0][field])+Decimal('.00001'))
    with pytest.raises(ValueError):
        verify.verify_alpha(row)


@pytest.mark.parametrize('mutation', ['omit_map', 'swap_mode', 'invent_mode', 'coupling',
                                     'omit_selector', 'invent_exponent', 'selector_score',
                                     'mass_ratio', 'boolean', 'tolerance'])
def test_alpha_selector_and_schema_attacks(packet, mutation):
    data = deepcopy(packet['calculations']['alpha'])
    if mutation == 'omit_map': data['rows'].pop()
    elif mutation == 'swap_mode': data['rows'][0], data['rows'][1] = data['rows'][1], data['rows'][0]
    elif mutation == 'invent_mode': data['rows'][0]['mode'] = 'best_fit'
    elif mutation == 'coupling': data['rows'][0]['couplings'][0] = '0.1'
    elif mutation == 'omit_selector': data['electron_menu'].pop()
    elif mutation == 'invent_exponent': data['electron_menu'][0]['exponent'] = 1
    elif mutation == 'selector_score': data['electron_menu'][0]['score'] = '0'
    elif mutation == 'mass_ratio': data['masses_over_v']['e'] = '0.5'
    elif mutation == 'boolean': data['electron_menu'][0]['exponent'] = True
    elif mutation == 'tolerance': data['rows'][0]['printed_residual_tolerance'] = '1e-3'
    with pytest.raises(ValueError):
        verify.verify_alpha(data)


@pytest.mark.parametrize('field', ['central_mev', 'excluded_small_root_mev', 'measured_Q',
                                   'first_order_sigma_max_mev'])
def test_tau_scalar_corruption(packet, field):
    data = deepcopy(packet['calculations']['tau'])
    data[field] = '1'
    with pytest.raises(ValueError):
        verify.verify_tau(data, verify.load(verify.HERE/'inputs.json'))


@pytest.mark.parametrize('field', ['corners_mev', 'measured_Q_corner_range', 'outward_mev',
                                   'jacobian', 'residual_in_sigma_with_arbitrary_correlations'])
@pytest.mark.parametrize('attack', ['change', 'omit', 'append', 'boolean'])
def test_tau_enclosure_census_and_sensitivity_corruption(packet, field, attack):
    data = deepcopy(packet['calculations']['tau'])
    if attack == 'change': data[field][0] = '0'
    elif attack == 'omit': data[field].pop()
    elif attack == 'append': data[field].append(data[field][0])
    else: data[field][0] = True
    with pytest.raises(ValueError):
        verify.verify_tau(data, verify.load(verify.HERE/'inputs.json'))


@pytest.mark.parametrize('section,field,value', [
    ('approximate_trunk', 'relative_pair_defect', '0'),
    ('approximate_trunk', 'printed_minus_converged_asymptotic', '0'),
    ('measured_pixel_diagnostic', 'P_measured', '1.63'),
    ('measured_pixel_diagnostic', 'mixed_inverse', '137.035999177'),
    ('historical_tau', 'frozen_center_verdict', 'COMPATIBLE'),
    ('historical_tau', 'window_robust_verdict', 'FAIL'),
    ('normal_error_reference_diagnostic', 'oph_to_koide', '1'),
    ('normal_error_reference_diagnostic', 'twice_log_free_mass_to_balanced', '0'),
    ('conclusions', 'new_natural_outcomes', True),
    ('conclusions', 'global_chance_probability', '0.00001'),
    ('conclusions', 'corners_are_joint_68_percent_interval', True),
    ('conclusions', 'tau_independently_tests_P', True),
    ('conclusions', 'natural_comparison', 'PASS')])
def test_report_promotions_fail_semantic_replay(packet, monkeypatch, section, field, value):
    # The alpha checker is exercised separately above; isolate report semantics.
    monkeypatch.setattr(verify, 'verify_alpha', lambda _: replay.context())
    data = deepcopy(packet['calculations'])
    data[section][field] = value
    with pytest.raises(ValueError):
        verify.verify_calculations(data)


@pytest.mark.parametrize('offset,sigma,expected', [
    ('0', '.045', 'COMPATIBLE'), ('.09', '.045', 'COMPATIBLE'),
    ('.090000000000000000001', '.045', 'INCONCLUSIVE'),
    ('.135', '.045', 'INCONCLUSIVE'), ('.135000000000000000001', '.045', 'FAIL'),
    ('0', '.045000000000000001', 'INCONCLUSIVE'),
    ('.301', '.1', 'FAIL'), ('0', '.09', 'INCONCLUSIVE'),
    ('-.135000000000000000001', '.045', 'FAIL')])
def test_exact_frozen_rule_edges(offset, sigma, expected):
    value = str(Decimal('1776.969027')+Decimal(offset))
    assert comparison.compare(dict(value=value, sigma=str(Decimal(sigma)), unit='MeV'))['frozen_center_verdict'] == expected


def test_rounded_center_cannot_silently_kill_entire_window():
    result = comparison.compare(dict(value='1776.968991', sigma='0.000000001', unit='MeV'))
    assert result['frozen_center_verdict'] == 'FAIL'
    assert result['window_robust_verdict'] == 'INCONCLUSIVE'
    result = comparison.compare(dict(value='1776.969027', sigma='0.045', unit='MeV'))
    assert result['window_robust_verdict'] == 'COMPATIBLE'
    result = comparison.compare(dict(value='1776', sigma='0.045', unit='MeV'))
    assert result['window_robust_verdict'] == 'FAIL'


@pytest.mark.parametrize('field,value', [
    ('unit', 'GeV'), ('unit', None), ('value', 'NaN'), ('value', 'Infinity'),
    ('value', float('nan')), ('value', True), ('value', '0'), ('value', '-1'),
    ('value', '1e999'), ('value', '9'*91), ('sigma', '0'), ('sigma', '-1'),
    ('sigma', '.045'), ('sigma', '1e-999'), ('sigma', '101'), ('sigma', 0.045)])
def test_invalid_observation_fails_closed(field, value):
    # Leading-dot strings are deliberately outside the canonical decimal grammar.
    payload = dict(value='1776.969027', sigma='0.045', unit='MeV')
    payload[field] = value
    with pytest.raises(ValueError): comparison.compare(payload)


def test_no_undeclared_fit_or_nuisance_fields():
    for payload in ({}, [], {'value':'1776.9', 'sigma':'0.04', 'unit':'MeV', 'offset':'0.1'}):
        with pytest.raises(ValueError): comparison.compare(payload)


@pytest.mark.parametrize('y', ['1700', '1776.969027', '1800'])
def test_identical_reference_likelihoods_for_all_synthetic_observations(y):
    result = comparison.reference_log_likelihood_ratios(y, '0.045', '1776.969027')
    assert result['oph_to_koide'] == '0'
    assert Decimal(result['twice_log_free_mass_to_balanced']) >= 0


def release(name='arxiv/2710.00001', date='2027-10-01T12:00:00Z', **changes):
    return dict(id=name, first_public_utc=date, kind='dedicated_measurement',
                observable='charged_tau_mass', exposed_by=[]) | changes


START = '2027-09-01T00:00:00Z'


def test_first_release_selection_is_order_independent_and_not_admission():
    a, b = release(), release('arxiv/2710.00002')
    for rows in ([a,b], [b,a]):
        result = admission.nominate(rows, START)
        assert result['selected'] == a['id']
        assert result['state'] == 'NOMINATED_NOT_ADMITTED'


def test_never_skip_exposed_first_release():
    a, b = release(exposed_by=['analyst']), release('arxiv/2711.00001', '2027-11-01T00:00:00Z')
    result = admission.nominate([b,a], START)
    assert result['selected'] == a['id'] and result['state'] == 'NOT_READY'


@pytest.mark.parametrize('change', [dict(first_public_utc=START), dict(kind='theory'),
                                   dict(kind='world_average'), dict(kind='revision'), dict(observable='other')])
def test_predeclared_metadata_exclusions(change):
    result = admission.nominate([release(**change)], START)
    assert result['state'] == 'NOT_READY' and result['selected'] is None


@pytest.mark.parametrize('change', [dict(value='1776.96'), dict(sigma='0.01'),
                                   dict(first_public_utc='2027-13-01T00:00:00Z'),
                                   dict(first_public_utc=True), dict(kind='favorable'),
                                   dict(observable='tau_like'), dict(exposed_by=False), dict(id='../mass')])
def test_metadata_rejects_outcomes_and_invalid_records(change):
    with pytest.raises(ValueError): admission.nominate([release(**change)], START)


def test_metadata_inventory_bounds_and_duplicate_releases():
    for rows in ([release(),release()], [release()]*101, {}, None):
        with pytest.raises(ValueError): admission.nominate(rows, START)


@pytest.mark.parametrize('name,field,value', [
    ('protocol.json', 'state', 'READY'), ('protocol.json', 'dataset', 'favorable-data'),
    ('protocol.json', 'stopping', 'keep trying'), ('protocol.json', 'nuisance_policy', 'fit offset'),
    ('inputs.json', 'covariance', 'independent Gaussian'), ('inputs.json', 'scope', 'held-out')])
def test_rebuilt_custody_cannot_loosen_reviewed_policy(monkeypatch, name, field, value):
    original = verify.load
    def modified(path):
        data = original(path)
        if Path(path).name == name: data[field] = value
        return data
    monkeypatch.setattr(verify, 'load', modified)
    with pytest.raises(ValueError): verify.verify_policy()


@pytest.mark.parametrize('attack', ['outcome_fields', 'hide_failed_query', 'drop_positive_control', 'add_candidate'])
def test_discovery_does_not_accept_outcomes_or_launder_failed_search(monkeypatch, attack):
    original = verify.load
    def modified(path):
        data = original(path)
        if Path(path).name == 'catalogue_discovery.json':
            if attack == 'outcome_fields': data['queries'][1]['rows'][0]['abstracts'] = ['synthetic outcome']
            elif attack == 'hide_failed_query': data['queries'].pop(0)
            elif attack == 'drop_positive_control':
                data['queries'][3]['rows'] = [r for r in data['queries'][3]['rows'] if r['id'] != '2663717']
            else: data['queries'][1]['rows'][0]['id'] = '123'
        return data
    monkeypatch.setattr(verify, 'load', modified)
    with pytest.raises(ValueError): verify.verify_policy()


@pytest.mark.parametrize('raw', [b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":Infinity}',
                                b'{"x":1e999}', b'{"x":'+b'9'*51+b'}', b'\xff',
                                b'['*2000+b']'*2000, b' '*200_001],
                         ids=['duplicate', 'nan', 'infinity', 'overflow', 'huge_integer',
                              'encoding', 'nesting', 'oversized'])
def test_hostile_json(tmp_path, raw):
    target = tmp_path/'bad.json'
    target.write_bytes(raw)
    with pytest.raises(ValueError): verify.load(target)


def test_false_math_with_rebuilt_source_pins_fails_under_optimized_python(tmp_path, packet):
    data = deepcopy(packet)
    data['sources'] = verify.pins()
    data['calculations']['alpha']['rows'][0]['alpha_inverse'] = '137.035999177'
    target = tmp_path/'false.json'
    target.write_text(json.dumps(data), encoding='utf-8')
    result = subprocess.run([sys.executable, '-O', '-m', 'independent_postdictions.verify', str(target)],
                            cwd=verify.ROOT, capture_output=True, text=True, timeout=120)
    assert result.returncode != 0 and 'independent numerical replay' in result.stderr


def test_positive_optimized_python_control():
    result = subprocess.run([sys.executable, '-O', '-m', 'independent_postdictions.verify'],
                            cwd=verify.ROOT, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr

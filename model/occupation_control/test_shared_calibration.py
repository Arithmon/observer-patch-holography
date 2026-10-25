"""Exact all-endpoint calibration obstruction and hostile receipt checks."""
from copy import deepcopy
from fractions import Fraction as F
import json

import pytest

from tools import strict_json
from . import shared_calibration as s


@pytest.fixture(scope='module')
def retained():
    return strict_json.load(s.ROOT/'model/occupation_control/shared_calibration_results.json')


@pytest.fixture(scope='module')
def rows():
    return s.source_rows()


def test_retained_exact_global_optimum_and_witness(retained, rows):
    optimum = s.verify(retained)
    assert F('0.49161779') < optimum < F('0.49161781')
    certificate = retained['certificate']
    assert [term['constraint'] for term in certificate['dual']] == [
        'a@6/5:lo', 'a@16/5:lo', 'a@42/5:hi']
    assert len(rows) == 125
    assert certificate['input_binding']['endpoint_count'] == 125
    assert len(certificate['input_binding']['constraint_order']) == 252
    assert certificate == s.solve_exact(rows)
    # The exact certificate is a relaxation: all 125 intervals can choose their
    # own interior value. It proves a lower bound for the actual shared dynamics.
    point = certificate['primal']
    g, b, e = (F(point[k]) for k in ('gain', 'offset', 'allowance'))
    assert g > 0
    for row in rows:
        lower, upper, observed = (F(row[k]) for k in ('lower', 'upper', 'observed'))
        assert max(F(0), g*lower+b-observed, observed-g*upper-b) <= e


def test_verifier_never_invokes_optimization_or_producer(retained, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Exact verifier must not solve or trust a producer')
    monkeypatch.setattr(s, 'solve_exact', forbidden)
    monkeypatch.setattr(s, '_solve_square', forbidden)
    monkeypatch.setattr(s, 'produce', forbidden)
    monkeypatch.setattr(s, 'independent_case_certificate', forbidden)
    assert s.verify(retained) > 0


@pytest.mark.parametrize('mutation', [
    'omit', 'duplicate', 'reorder', 'interval', 'observation', 'reverse', 'nonfinite', 'extra',
])
def test_endpoint_mutations_cannot_keep_certificate(retained, rows, mutation):
    changed = deepcopy(rows)
    if mutation == 'omit':
        changed.pop()
    elif mutation == 'duplicate':
        changed[1] = deepcopy(changed[0])
    elif mutation == 'reorder':
        changed[1], changed[2] = changed[2], changed[1]
    elif mutation == 'interval':
        changed[-1]['upper'] = str(F(changed[-1]['upper'])+1)
    elif mutation == 'observation':
        changed[-1]['observed'] = '0'
    elif mutation == 'reverse':
        changed[0]['lower'], changed[0]['upper'] = changed[0]['upper'], changed[0]['lower']
    elif mutation == 'nonfinite':
        changed[-1]['observed'] = 'NaN'
    else:
        changed[-1]['trash'] = 1
    with pytest.raises(ValueError):
        s.verify_exact(changed, retained['certificate'])


@pytest.mark.parametrize('mutation', [
    'omit_constraint', 'duplicate_constraint', 'reorder_constraint', 'bool_count',
    'lower_objective', 'negative_gain', 'bool_gain', 'nan_gain', 'noncanonical',
    'dual_sign', 'dual_weight', 'dual_duplicate', 'dual_reorder', 'foreign_dual',
    'extra', 'source', 'claim', 'float_nan',
])
def test_adversarial_receipts_fail_closed(retained, rows, monkeypatch, mutation):
    monkeypatch.setattr(s, 'source_rows', lambda: deepcopy(rows))
    changed = deepcopy(retained)
    cert = changed['certificate']
    order = cert['input_binding']['constraint_order']
    primal = cert['primal']
    dual = cert['dual']
    if mutation == 'omit_constraint':
        order.pop()
    elif mutation == 'duplicate_constraint':
        order[1] = order[0]
    elif mutation == 'reorder_constraint':
        order[0], order[1] = order[1], order[0]
    elif mutation == 'bool_count':
        cert['input_binding']['endpoint_count'] = True
    elif mutation == 'lower_objective':
        primal['allowance'] = str(F(primal['allowance'])-F(1, 10**12))
    elif mutation == 'negative_gain':
        primal['gain'] = '-1'
    elif mutation == 'bool_gain':
        primal['gain'] = True
    elif mutation == 'nan_gain':
        primal['gain'] = 'NaN'
    elif mutation == 'noncanonical':
        primal['gain'] = '2/2'
    elif mutation == 'dual_sign':
        dual[0]['weight'] = str(-F(dual[0]['weight']))
    elif mutation == 'dual_weight':
        dual[0]['weight'] = str(F(dual[0]['weight'])+F(1, 10**12))
    elif mutation == 'dual_duplicate':
        dual.append(deepcopy(dual[0]))
    elif mutation == 'dual_reorder':
        dual.reverse()
    elif mutation == 'foreign_dual':
        dual[0]['constraint'] = 'a@7/5:lo'
    elif mutation == 'extra':
        cert['trash'] = 'ignored?'
    elif mutation == 'source':
        changed['source']['sha256'] = '0'*64
    elif mutation == 'claim':
        changed['contract']['statistical_rejection_claimed'] = True
    else:
        changed['contract']['experimental_uncertainty_applied'] = float('nan')
    with pytest.raises(ValueError):
        s.verify(changed)


def test_gain_zero_boundary_is_certified_with_its_dual_constraint():
    # Decreasing point values cannot be followed by a nonnegative gain.
    rows = [{'id': 'first', 'lower': '0', 'upper': '0', 'observed': '1'},
            {'id': 'second', 'lower': '1', 'upper': '1', 'observed': '0'}]
    certificate = s.solve_exact(rows)
    assert certificate['primal'] == {'gain': '0', 'offset': '1/2', 'allowance': '1/2'}
    assert any(term['constraint'] == 'gain_nonnegative' for term in certificate['dual'])
    assert s.verify_exact(rows, certificate) == F(1, 2)


def test_zero_remaining_allowance_boundary_can_be_verified():
    rows = [{'id': 'first', 'lower': '0', 'upper': '0', 'observed': '0'},
            {'id': 'second', 'lower': '1', 'upper': '1', 'observed': '1'}]
    certificate = s.solve_exact(rows)
    assert s.verify_exact(rows, certificate) == 0


def test_duplicate_json_keys_rejected_before_verification(retained):
    text = json.dumps(retained)
    text = text.replace('"primal": {', '"primal": {}, "primal": {', 1)
    with pytest.raises(ValueError, match='duplicate'):
        strict_json.loads(text)


def test_artifact_bytes_are_pinned_before_parsing(monkeypatch):
    monkeypatch.setattr(s.Path, 'read_bytes', lambda self: b'{"not": "the retained certificate"}')
    def forbidden(*args, **kwargs):
        raise AssertionError('A changed artifact must fail before receipt parsing')
    monkeypatch.setattr(s.transfer, 'read_receipt', forbidden)
    with pytest.raises(ValueError, match='Pinned transfer'):
        s.source_rows()


def test_independent_case_optimum_reuses_both_exact_witnesses(retained, rows, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('A second optimization is unnecessary')
    monkeypatch.setattr(s, 'solve_exact', forbidden)
    common, extension = retained['certificate'], retained['independent_cases']
    assert extension == s.independent_case_certificate(rows, common)
    optimum = s.verify_independent_cases(rows, common, extension)
    # Independent reconstruction of the restricted dual; every positive term
    # is an endpoint of case a, so the other four parameter pairs cancel out.
    selected = [row for row in rows if row['id'].startswith('a@')]
    certificate = deepcopy(common)
    certificate['input_binding'] = s._binding(selected)
    assert len(selected) == 25
    assert s.verify_exact(selected, certificate) == optimum
    point = common['primal']
    gain, offset = F(point['gain']), F(point['offset'])
    for case in 'abcde':
        own = [row for row in rows if row['id'].startswith(case+'@')]
        assert len(own) == 25
        assert all(gain*F(row['lower'])+offset-optimum <= F(row['observed'])
                   <= gain*F(row['upper'])+offset+optimum for row in own)


@pytest.mark.parametrize('gain,offset', [('0', '0'), ('1', '-1/3'), ('17/5', '2')])
def test_heterogeneous_allowances_obey_the_exact_weighted_bound(retained, rows, gain, offset):
    allowance = retained['independent_cases']['heterogeneous_allowance']
    endpoints = {row['id']: row for row in rows}
    weights = {term['endpoint']: F(term['weight']) for term in allowance['terms']}
    assert list(weights) == ['a@6/5', 'a@16/5', 'a@42/5']
    assert sum(weights.values()) == 1
    gain, offset = F(gain), F(offset)
    # Exact heterogeneous allowances sufficient for this chosen calibration.
    required = {}
    for name in weights:
        row = endpoints[name]
        required[name] = max(F(), gain*F(row['lower'])+offset-F(row['observed']),
                             F(row['observed'])-gain*F(row['upper'])-offset)
    assert sum(weights[name]*required[name] for name in weights) >= F(allowance['lower_bound'])
    # Check cancellation directly without the producer's constraint helpers.
    first, middle, last = [endpoints[name] for name in weights]
    w0, w1, w2 = weights.values()
    assert w0+w1-w2 == 0
    assert w0*F(first['lower'])+w1*F(middle['lower'])-w2*F(last['upper']) == 0
    assert w2*F(last['observed'])-w0*F(first['observed'])-w1*F(middle['observed']) == F(allowance['lower_bound'])


@pytest.mark.parametrize('mutation', [
    'case', 'foreign_case', 'bool_case', 'omit_case', 'duplicate_case', 'reorder_case',
    'case_hash', 'case_count', 'case_extra', 'objective', 'time', 'other_case_time',
    'weight', 'negative_weight', 'bool_weight', 'noncanonical_weight', 'duplicate_term',
    'reorder_terms', 'lower_bound', 'extra', 'missing_extension', 'old_schema',
])
def test_independent_case_receipt_mutations_fail_closed(retained, rows, monkeypatch, mutation):
    monkeypatch.setattr(s, 'source_rows', lambda: deepcopy(rows))
    changed = deepcopy(retained)
    extension = changed['independent_cases']
    bindings = extension['case_bindings']
    allowance = extension['heterogeneous_allowance']
    if mutation == 'case':
        extension['witness_case'] = 'b'
    elif mutation == 'foreign_case':
        extension['witness_case'] = 'f'
    elif mutation == 'bool_case':
        extension['witness_case'] = True
    elif mutation == 'omit_case':
        bindings.pop()
    elif mutation == 'duplicate_case':
        bindings[-1] = deepcopy(bindings[0])
    elif mutation == 'reorder_case':
        bindings.reverse()
    elif mutation == 'case_hash':
        bindings[-1]['endpoints_sha256'] = '0'*64
    elif mutation == 'case_count':
        bindings[0]['endpoint_count'] = True
    elif mutation == 'case_extra':
        bindings[0]['extra'] = 'ignored?'
    elif mutation == 'objective':
        extension['uniform_allowance'] = '0'
    elif mutation == 'time':
        allowance['terms'][0]['endpoint'] = 'a@7/5'
    elif mutation == 'other_case_time':
        allowance['terms'][0]['endpoint'] = 'b@6/5'
    elif mutation == 'weight':
        allowance['terms'][0]['weight'] = str(F(allowance['terms'][0]['weight'])+F(1, 10**12))
    elif mutation == 'negative_weight':
        allowance['terms'][0]['weight'] = str(-F(allowance['terms'][0]['weight']))
    elif mutation == 'bool_weight':
        allowance['terms'][0]['weight'] = True
    elif mutation == 'noncanonical_weight':
        allowance['terms'][2]['weight'] = '2/4'
    elif mutation == 'duplicate_term':
        allowance['terms'].append(deepcopy(allowance['terms'][0]))
    elif mutation == 'reorder_terms':
        allowance['terms'].reverse()
    elif mutation == 'lower_bound':
        allowance['lower_bound'] = '0'
    elif mutation == 'extra':
        extension['verified'] = True
    elif mutation == 'missing_extension':
        del changed['independent_cases']
    else:
        changed['schema'] = 'oph.occupation_control.shared_calibration.v1'
    with pytest.raises(ValueError):
        s.verify(changed)


@pytest.mark.parametrize('new_id', ['a@7/5', 'a@1.2', 'b@6/5', 'f@6/5'])
def test_rebinding_cannot_relabel_a_witness_endpoint(retained, rows, new_id):
    changed = deepcopy(rows)
    changed[2]['id'] = new_id
    common = deepcopy(retained['certificate'])
    common['dual'][0]['constraint'] = new_id+':lo'
    # Rebinding a general LP certificate does not authenticate changed time or
    # case labels. The fixed independent-case contract must reject it itself.
    with pytest.raises(ValueError):
        common['input_binding'] = s._binding(changed)
        s.verify_independent_cases(changed, common, retained['independent_cases'])

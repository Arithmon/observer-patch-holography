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

"""Independent source-block, elimination and approximation controls."""
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as sp

from ..observer.scalar_response.experiment import compare_receipt
from . import source_green as g


def test_exact_two_carrier_elimination_against_full_port_inverse():
    _, local = g.passage.recovered_carrier()
    internal = sp.diag(local, local)
    d = sp.zeros(24, 1); d[0], d[12] = 1, -1
    l = internal+d*d.T
    v = sp.diag(sp.ones(12, 1), sp.ones(12, 1))
    g0 = sp.diag((local+sp.ones(12)/12).inv()-sp.ones(12)/12,
                 (local+sp.ones(12)/12).inv()-sp.ones(12)/12)
    r = d.T*v
    s = r.T*(sp.eye(1)+d.T*g0*d).inv()*r/12
    actual = v.T*((l+sp.ones(24)/24).inv()-sp.ones(24)/24)*v/12
    assert actual == (s+sp.ones(2)/2).inv()-sp.ones(2)/2
    a = r.T*r/12
    assert s != a  # Simply replacing each carrier by its mean misses relaxation.


@pytest.mark.parametrize('level', [0, 1, 2])
def test_actual_source_matching_blocks_and_native_green(level):
    result = g.diagnostic(level)
    # Live output must have the same scalar types as its retained JSON.
    compare_receipt(json.loads(json.dumps(result)), result)
    assert result['carriers'] == 20*4**level
    assert result['ports'] == 240*4**level
    assert result['all_intercarrier_endpoints_distinct']
    assert abs(result['cell_contrast_native_Green_diagnostic']-
               result['cell_contrast_Schur_Green_diagnostic']) < 2e-9
    errors = result['controlled_Schur_approximations']
    assert errors[-1]['relative_Green_error_diagnostic'] < errors[0]['relative_Green_error_diagnostic']
    for row in errors:
        assert row['relative_Green_error_diagnostic'] <= float(F(row['relative_Green_error_bound']))


def test_matching_hypothesis_is_checked_before_using_factor_two(monkeypatch):
    original = g.source.federation(0)
    edges = list(original.inter_seams)
    first, other = edges[0][0], edges[1][1]
    edges[1] = tuple(sorted((first, other)))
    monkeypatch.setattr(g.source, 'federation', lambda level: replace(original, inter_seams=tuple(edges)))
    with pytest.raises(ArithmeticError, match='matching'):
        g.matrices(0)


def test_polynomial_error_is_uniform_and_does_not_require_commutation_of_graph_blocks():
    # F and the two-cell edge geometry need not be simultaneously diagonal.
    f = np.array([[1.2, .1], [.1, 1.6]])
    r = np.array([[1., -1., 0.], [0., 1., -1.]])
    exact = g._green(r.T@np.linalg.solve(f, r)/12)
    for order in (1, 4, 8):
        approximate = g._green(g.schur_series(f, r, order))
        difference = approximate-exact
        bound = float(g.approximation_relative_bound(order))*exact-difference
        assert np.linalg.eigvalsh(difference)[0] > -1e-12
        assert np.linalg.eigvalsh(bound)[0] > -1e-12
    assert g.approximation_relative_bound(16) == F(1, 131071)


@pytest.mark.parametrize('order', [True, 0, -1, 33, 1.0])
def test_approximation_order_requires_declared_integer(order):
    with pytest.raises(ValueError):
        g.approximation_relative_bound(order)

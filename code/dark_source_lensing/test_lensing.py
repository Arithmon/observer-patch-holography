"""Scientific and hostile-input tests; no acceptance solely by regenerated hashes."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
from fractions import Fraction as F
import mpmath
import pytest
from . import check, model, verify
from .format import decimal, load, rational


@pytest.fixture(scope='module')
def evidence():
    return model.build()


def test_committed_receipt_and_independent_replay(evidence):
    row = verify.verify()
    assert row['evidence'] == evidence


@pytest.mark.parametrize('path,value', [
    (('interpretation', 'physical_promotion'), True),
    (('interpretation', 'full_issue_751_closed'), True),
    (('interpretation', 'natural_data_used'), True),
    (('interpretation', 'numeric_controls'), 'rigorous_intervals'),
    (('annulus', 0, 'u'), '1/100'),
    (('annulus', 0, 'R_over_r0'), True),
    (('annulus', 0, 'branch'), 'zero_radial'),
    (('annulus', 0, 'kappa'), '1/2'),
    (('annulus', 0, 'radial_ratio'), '0'),
    (('annulus', 0, 'tangential_ratio'), '0'),
    (('annulus', 0, 'alpha'), '0'),
    (('annulus', 0, 'sweep'), '3.141592653589793'),
    (('annulus', 0, 'inferred_kappa'), '0'),
    (('annulus', 26, 'alpha'), '0'),
    (('cold_bounds', 0, 'lower'), '0'),
    (('cold_bounds', 8, 'upper'), '1/2'),
    (('global_samples', 4, 'm'), '0'),
    (('global_samples', 4, 'u'), '1'),
    (('global_samples', 5, 'density_scaled'), '0'),
    (('global_samples', 5, 'radial_scaled'), '0'),
    (('global_samples', 5, 'tangential_scaled'), '0'),
    (('global_bounds', 'density_lower'), '0'),
    (('global_bounds', 'radial_upper'), '1'),
    (('global_bounds', 'tangential_upper'), '0'),
    (('global_bounds', 'compactness_upper'), '0'),
    (('global_bounds', 'bending_gap_per_epsilon'), '0'),
    (('global_bounds', 'bending_gap_at_epsilon_001'), '1/100000'),
    (('global_ray_enclosure', 'integral_lower'), '0'),
    (('global_ray_enclosure', 'integral_upper'), '1'),
    (('global_ray_enclosure', 'rows', 1, 'lower'), '0'),
    (('global_ray_enclosure', 'rows', 0, 'upper'), '1'),
    (('abel', 0, 'Q_minus_one'), '0'),
    (('abel', 0, 'alpha'), '0'),
    (('abel', 14, 'alpha'), '0'),
    (('abel', 14, 'terms'), [[1, '1/3']]),
])
def test_resealed_semantic_corruption_rejected(evidence, path, value):
    packet = copy.deepcopy(evidence)
    target = packet
    for name in path[:-1]:
        target = target[name]
    target[path[-1]] = value
    # Custody is intentionally bypassed: this must fail for its mathematics.
    with pytest.raises(ValueError):
        check.verify(packet)


@pytest.mark.parametrize('field', ['annulus', 'cold_bounds', 'global_samples', 'abel'])
@pytest.mark.parametrize('mode', ['omit', 'duplicate', 'reverse'])
def test_complete_case_catalogues(evidence, field, mode):
    packet = copy.deepcopy(evidence)
    if mode == 'omit':
        packet[field].pop()
    elif mode == 'duplicate':
        packet[field][-1] = copy.deepcopy(packet[field][0])
    else:
        packet[field].reverse()
    with pytest.raises(ValueError):
        check.verify(packet)


def test_wrong_side_tiny_enclosure_rejected(evidence):
    row = copy.deepcopy(evidence['global_ray_enclosure'])
    row['integral_lower'] = str(F(row['integral_lower'])+F(1, 10**18))
    with pytest.raises(ValueError, match='outward'):
        check.verify_ray_enclosure(row)


def test_slightly_wrong_bending_is_not_rounded_away(evidence):
    packet = copy.deepcopy(evidence)
    ctx = mpmath.mp.clone()
    ctx.dps = 70
    value = ctx.mpf(packet['annulus'][0]['alpha'])*(1+ctx.mpf('1e-40'))
    packet['annulus'][0]['alpha'] = ctx.nstr(value, 60)
    with pytest.raises(ValueError, match='bending'):
        check.verify(packet)


def test_pressure_endpoint_classification_exhaustive_rational_controls():
    # Direct energy inequalities versus claimed interval over a two-dimensional grid.
    for u in (F(i, 31) for i in range(1, 31)):
        endpoints = (u/(2+2*u), u/(1+2*u), u*u/(2+2*u*u))
        kappas = [F(j, 202) for j in range(1, 101)]+list(endpoints)
        for k in kappas:
            radial = (1-2*k)*u-k
            tangential = (1-2*k)*u*u/2
            dec = abs(radial) <= k and abs(tangential) <= k
            assert dec == (u/(2+2*u) <= k < F(1, 2))
            assert (dec and radial >= 0 and tangential >= 0) == (endpoints[0] <= k <= endpoints[1])


def test_cold_limit_needs_tangential_support():
    for u in (F(1, 10**6), F(1, 100), F(1, 5)):
        k = u/(1+2*u)
        assert (1-2*k)*u-k == 0
        assert (1-2*k)*u*u/(2*k) == u/2 > 0
        # Newtonian identification retains tension instead of silently setting it to zero.
        assert ((1-2*u)*u-u)/u == -2*u


@pytest.mark.parametrize('value', [True, 1, 0.1, None, [], 'NaN', 'Inf', '1/0', '2/4', '+1', '00', '1e999999'])
def test_noncanonical_rationals_rejected(value):
    with pytest.raises(ValueError):
        rational(value)


@pytest.mark.parametrize('value', [True, 1, 0.1, None, [], 'nan', 'inf', '1/2', ' 1', '1e99999'])
def test_nonfinite_and_wrong_numeric_types_rejected(value):
    ctx = mpmath.mp.clone()
    with pytest.raises(ValueError):
        decimal(ctx, value)


@pytest.mark.parametrize('raw', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":1.0}', '{"x":Infinity}'])
def test_strict_json(tmp_path, raw):
    path = tmp_path/'bad.json'
    path.write_text(raw, encoding='ascii')
    with pytest.raises(ValueError):
        load(path)


def test_extra_fields_and_empty_packet_rejected(evidence):
    for packet in ({}, {**evidence, 'passed':True}, None, []):
        with pytest.raises(ValueError):
            check.verify(packet)


def test_source_custody_rejected(tmp_path):
    row = load(verify.HERE/'receipt.json')
    row['sources'].pop(next(iter(row['sources'])))
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(row), encoding='ascii')
    with pytest.raises(ValueError, match='custody'):
        verify.verify(path)


def test_global_mpmath_context_is_not_an_input(evidence):
    previous = mpmath.mp.dps
    try:
        mpmath.mp.dps = 5
        assert model.build() == evidence
        check.verify(evidence)
    finally:
        mpmath.mp.dps = previous


def test_optimized_cli_without_producer_imports(tmp_path):
    program = '''
import importlib.abc, runpy, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in ('dark_source_lensing.model', 'dark_source_lensing.build'):
            raise RuntimeError('producer import forbidden')
sys.meta_path.insert(0, Block())
sys.argv = ['verify', sys.argv[1]]
runpy.run_module('dark_source_lensing.verify', run_name='__main__')
'''
    env = {**os.environ, 'PYTHONPATH':str(verify.ROOT/'code')}
    valid = subprocess.run([sys.executable, '-O', '-c', program, str(verify.HERE/'receipt.json')],
                           env=env, capture_output=True, text=True, timeout=120)
    assert valid.returncode == 0, valid.stderr
    row = load(verify.HERE/'receipt.json')
    row['evidence']['annulus'][0]['alpha'] = '0'
    path = tmp_path/'resealed_bad.json'
    path.write_text(json.dumps(row), encoding='ascii')
    invalid = subprocess.run([sys.executable, '-O', '-c', program, str(path)],
                             env=env, capture_output=True, text=True, timeout=120)
    assert invalid.returncode != 0
    assert 'ValueError: finite-endpoint bending' in invalid.stderr
    assert 'producer import forbidden' not in invalid.stderr

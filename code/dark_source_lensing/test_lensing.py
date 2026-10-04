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
from . import check, global_ray, model, verify
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
    (('annulus', 0, 'R_over_r0'), 1),
    (('annulus', 0, 'branch'), 'zero_radial'),
    (('annulus', 0, 'kappa'), '1/2'),
    (('annulus', 0, 'radial_ratio'), '0'),
    (('annulus', 0, 'tangential_ratio'), '0'),
    (('annulus', 0, 'alpha'), '0'),
    (('annulus', 0, 'sweep'), '3.141592653589793'),
    (('annulus', 0, 'inferred_kappa'), '0'),
    (('annulus', 26, 'alpha'), '0'),
    (('winding_controls', 1, 'turns'), 0),
    (('winding_controls', 1, 'sweep_over_pi'), '1'),
    (('winding_controls', 1, 'alpha_over_pi'), '1/3'),
    (('winding_controls', 1, 'reduced_sweep_over_pi'), '3'),
    (('winding_controls', 1, 'kappa'), '11/72'),
    (('cold_bounds', 0, 'lower'), '0'),
    (('cold_bounds', 8, 'upper'), '1/2'),
    (('kinetic_bounds', 0, 'speed_squared_cap'), '0'),
    (('kinetic_bounds', 0, 'lower_kappa'), '0'),
    (('kinetic_bounds', 8, 'bending_ratio_lower'), '1'),
    (('kinetic_bounds', 8, 'bending_ratio_upper'), '0'),
    (('kinetic_bounds', 8, 'sqrt_low_kappa'), ['1', '1']),
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


@pytest.mark.parametrize('field', ['annulus', 'winding_controls', 'cold_bounds', 'kinetic_bounds', 'global_samples', 'abel'])
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


def test_unknown_winding_retains_distinct_mass_solutions(evidence):
    rows = evidence['winding_controls']
    check.verify_winding(rows)
    assert len({r['kappa'] for r in rows}) == 3
    assert len({r['endpoint_angle_over_pi'] for r in rows}) == 1
    assert len({r['reduced_sweep_over_pi'] for r in rows}) == 1
    assert [F(r['sweep_over_pi']) % 2 for r in rows] == [1, 1, 1]


@pytest.mark.parametrize('epsilon', [F(1, 100), F(1, 10**6)])
def test_full_source_ray_precision_refinement(epsilon):
    ctx = mpmath.mp.clone()
    ctx.dps = 60
    coarse = ctx.mpf(global_ray.bending_difference(epsilon, 32))
    fine = ctx.mpf(global_ray.bending_difference(epsilon, 48))
    assert abs(coarse-fine) < abs(fine)*ctx.mpf('1e-22')


def test_full_ray_control_rejects_wrong_geometric_enclosure(evidence):
    # Bypass the exact-bound checker to challenge this independent integral control.
    rows = copy.deepcopy(evidence['global_ray_enclosure']['rows'])
    rows[0]['upper'] = str(F(rows[0]['lower'])+F(1, 10**20))
    with pytest.raises(ValueError, match='full-source null integral'):
        global_ray.verify_enclosures(rows)


@pytest.mark.parametrize('epsilon', [F(0), F(-1), F(1), 0.01, True, None])
def test_ray_integrator_rejects_out_of_domain_strength(epsilon):
    with pytest.raises(ValueError, match='strength domain'):
        global_ray.bending_difference(epsilon)


def test_cached_integrator_does_not_accept_wrong_precision_type():
    global_ray.bending_difference(F(1, 100), 40)
    with pytest.raises(ValueError, match='precision'):
        global_ray.bending_difference(F(1, 100), 40.0)


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


def test_kinetic_threshold_and_endpoint_classification():
    for u in (F(i, 19) for i in range(1, 19)):
        for cap in (F(0), u/2, u, (1+u)/2):
            lower = u*(1+u)/(1+2*u+2*u*u+cap)
            upper = u/(1+2*u)
            assert (lower <= upper) == (cap >= u)
            for k in (lower, upper, (lower+upper)/2, lower/2, (upper+F(1, 2))/2):
                radial = u-(1+2*u)*k
                tangential = (1-2*k)*u*u/2
                allowed = radial >= 0 and tangential >= 0 and radial+2*tangential <= cap*k
                assert allowed == (lower <= k <= upper)


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


def test_oversized_receipt_is_rejected(tmp_path):
    path = tmp_path/'oversized.json'
    path.write_bytes(b' '*250_001)
    with pytest.raises(ValueError, match='size limit'):
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

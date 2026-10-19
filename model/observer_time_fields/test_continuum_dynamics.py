"""Common native continuum generator, time parameters and wave-path controls."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from model.observer.scalar_response.experiment import compare_receipt
from . import continuum_dynamics as d


@pytest.mark.parametrize('level', (0, 1, 3, 8))
def test_native_integer_step_and_coherent_clock_errors_use_actual_normalizations(level):
    bounds = d.propagation_bounds(level, band=16, horizon=4)
    q, n = 630*16**level, 1260*16**level
    for a in (F(1, 10), F(1), F(7), F(16)):
        for parameter in (F(1, 13), F(7, 10), F(4)):
            with mp.workdps(60):
                a0, t = mp.mpf(a.numerator)/a.denominator, mp.mpf(parameter.numerator)/parameter.denominator
                attempts = int(mp.floor(t*n))
                heat_error = abs((1-a0/n)**attempts-mp.exp(-t*a0))
                steps = int(mp.floor(t*mp.sqrt(q)))
                theta = 2*mp.asin(mp.sqrt(a0/(4*q)))
                wave_error = abs(mp.cos(steps*theta)-mp.cos(t*mp.sqrt(a0)))
            assert heat_error <= float(bounds['integer_attempt_heat_error_on_band'])
            assert wave_error <= float(bounds['integer_coherent_wave_error_on_band'])
    assert bounds['spatial_subsequence_convergence_rate_supplied'] is False


@pytest.mark.parametrize('level', (0, 2, 5, 10))
def test_actual_positive_length_polynomial_freezes_at_fixed_cap(level):
    cap, band = 8, F(16)
    value = d.fixed_cap_collapse(level, cap, band)
    coefficients = d.wave_reference.coefficients(cap)
    with mp.workdps(80):
        q = 630*mp.mpf(16)**level
        x = 1-mp.mpf(16)/(2*q)
        weights = tuple(mp.mpf(c.numerator)/c.denominator for c in coefficients)
        theta = -mp.log(sum(c*x**i for i, c in enumerate(weights))/sum(weights))
        scaled = mp.sqrt(q)*theta
        ideal = mp.sqrt(q)*mp.acos(x)
    assert 0 < scaled <= float(value['scaled_modular_frequency_upper'])
    assert ideal >= 4-mp.mpf('1e-60')
    if level >= 5:
        assert scaled < .01


def test_growing_cap_budget_has_finite_scaled_error_and_distinct_clock_parameter():
    old = None
    for level in range(6):
        value = d.propagation_bounds(level)
        assert value['sufficient_native_return_cap'] == 64*256**level
        if old is not None:
            assert value['joint_cap_modular_cosine_error_upper'] <= old
        old = value['joint_cap_modular_cosine_error_upper']
    # The passive attempt and coherent step units diverge under refinement.
    assert F(1260*16**5, 1260*16**4) == 16
    assert F(4**5, 4**4) == 4


def test_native_green_wave_covariance_has_identity_velocity_not_quantum_vacuum():
    value = d.exact_wave_synthesis()
    assert value['constant_mode_exactly_removed']
    _, _, _, _, _, source = d.dynamics._base_blocks()
    eigen, vectors = np.linalg.eigh(np.asarray(source, dtype=float))
    positive = eigen[1:]
    covariance = 1/positive-1/2520
    assert np.all(covariance > 0)
    np.testing.assert_allclose(positive*covariance, 1-positive/2520, atol=1e-15)
    assert np.linalg.norm(1/positive-1/(2*np.sqrt(positive))) > 1
    for t, u in ((0., .3), (.5, 1.2), (2., 2.)):
        actual = covariance*(np.cos(t*np.sqrt(positive))*np.cos(u*np.sqrt(positive))+
                             np.sin(t*np.sqrt(positive))*np.sin(u*np.sqrt(positive)))
        expected = covariance*np.cos((t-u)*np.sqrt(positive))
        np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=2e-14)


def test_reusing_one_gaussian_block_breaks_stationary_wave_covariance():
    # Two disjoint CLT increments are essential for the declared Gaussian law.
    a, t = 3., .4
    one_block = (np.cos(t*np.sqrt(a))+np.sin(t*np.sqrt(a)))**2/a
    two_blocks = (np.cos(t*np.sqrt(a))**2+np.sin(t*np.sqrt(a))**2)/a
    assert abs(two_blocks-1/a) < 1e-15
    assert abs(one_block-1/a) > .1


def test_weak_compressed_static_data_alone_do_not_determine_the_resolvent():
    # In a growing basis the second coordinate can escape every fixed dictionary.
    # G_n=I+(|1><n|+|n><1|)/2 converges weakly to I, but not strongly.
    block = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 1]])
    assert block[0, 0] == 1
    actual = (block*(sp.eye(2)+block).inv())[0, 0]
    assert actual == sp.Rational(7, 15)
    assert actual != sp.Rational(1, 2)
    # Source compactness excludes the escaping off-diagonal mass in this control.
    assert (block-sp.eye(2))[:, 0].dot((block-sp.eye(2))[:, 0]) == sp.Rational(1, 4)


def test_actual_full_native_coherent_oracle_keeps_both_Ramsey_outcomes():
    value = d.native_propagator_control()
    assert value['source_ports'] == 240
    assert value['full_native_seams'] == 630
    rows = value['readouts']
    for row in rows:
        assert abs(row['Ramsey_plus_diagnostic']+row['Ramsey_minus_diagnostic']-1) < 1e-15
        assert row['spectral_oracle_residual_diagnostic'] < 3e-12
    assert rows[0]['Ramsey_plus_diagnostic'] == 1
    assert rows[-1]['Ramsey_minus_diagnostic'] > .01


def test_wave_path_tail_controls_the_supremum_and_tends_to_zero():
    for cutoff in (0, 10, 100, 1000):
        value = d.wave_path_bounds(cutoff, F(3, 2))
        expected = 4*(d.spatial.GREEN_BOUND+F(9, 4))/(cutoff+1)**2
        assert value['expected_squared_supremum_tail_upper'] == expected
        assert value['expected_fourth_increment_upper_per_dt_fourth'] == 12
    assert d.wave_path_bounds(1000)['expected_squared_supremum_tail_upper'] < F(1, 100)


@pytest.mark.parametrize('bad', (True, -1, 1., '1'))
def test_ambiguous_levels_caps_and_cutoffs_rejected(bad):
    with pytest.raises(ValueError):
        d.propagation_bounds(bad)
    with pytest.raises(ValueError):
        d.fixed_cap_collapse(0, bad)
    with pytest.raises(ValueError):
        d.wave_path_bounds(bad)


def test_retained_continuum_receipt_replays_and_rejects_physical_promotion():
    result = d.run()
    path = Path(__file__).with_name('continuum_dynamics_results.json')
    compare_receipt(json.loads(path.read_text()), result)
    changed = copy.deepcopy(result)
    changed['scope']['observed_one_plus_three_spacetime_identified'] = True
    with pytest.raises(ValueError):
        compare_receipt(changed, result)


def test_fresh_runtime_source_closure_is_pinned():
    command = '''
import sys
from pathlib import Path
from model.observer_time_fields import continuum_dynamics as d
d.run()
missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path is None:continue
    try:name=Path(path).resolve().relative_to(d.ROOT).as_posix()
    except ValueError:continue
    if name.endswith('.py') and name not in d.SOURCE_PATHS:missing.append(name)
assert not missing,missing
'''
    subprocess.run([sys.executable, '-c', command], cwd=d.ROOT, check=True)

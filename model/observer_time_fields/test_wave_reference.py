"""Positive duration law, full native closure, modular cap and wrong-state controls."""
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
from . import wave_reference as w


def test_exact_recurrence_matches_the_generating_function_differential_equation():
    x = sp.Symbol('x')
    values = w.coefficients(24)
    polynomial = sum(sp.Rational(value)*x**n for n, value in enumerate(values))
    residual = sp.Poly((1-x*x)*sp.diff(polynomial,x,2)-x*sp.diff(polynomial,x)-polynomial,x)
    assert all(residual.nth(n) == 0 for n in range(23))
    assert values[:5] == (1,1,F(1,2),F(1,3),F(5,24))
    assert all(value > 0 for value in values)
    assert sp.exp(sp.asin(x)).series(x,0,9).removeO() == sum(sp.Rational(values[n])*x**n for n in range(9))


def test_rational_Machin_exponential_enclosure_contains_independent_high_precision_value():
    value = w.normalization()
    with mp.workdps(100):
        def number(x):return mp.mpf(x.numerator)/x.denominator
        lo,hi = value['pi_interval']
        assert number(lo) < mp.pi < number(hi)
        lo,hi = value['c0_interval']
        target = mp.exp(-mp.pi/2)
        assert number(lo) < target < number(hi)
        assert number(value['implemented_c0']) <= number(lo)
        assert hi-lo < F(1,10**27)


@pytest.mark.parametrize('cap',(0,1,2,8,64,512))
def test_finite_sampler_keeps_all_length_and_cap_probabilities(cap):
    value = w.sampler(cap)
    assert len(value['weights']) == cap+1
    assert sum(value['weights'])+value['cap'] == 1
    assert sum(value['uniform_integer_counts']) == value['uniform_denominator']
    assert all(x > 0 for x in value['uniform_integer_counts'])
    assert 0 <= value['expected_native_attempts'] <= cap
    for n,weight in enumerate(value['weights']):
        assert weight == w.normalization()['implemented_c0']*w.coefficients(cap)[n]


@pytest.mark.parametrize('cap',(0,1,2,3))
def test_every_native_length_and_endpoint_matches_independent_repair_enumeration(cap):
    assert w.complete_law(cap) == w.endpoint_oracle(cap)
    value = w.reference(cap)
    assert value['return_probability']+value['different_probability']+value['cap_probability'] == 1
    assert value['return_probability'] >= w.normalization()['implemented_c0']
    if cap:
        assert value['different_probability'] > 0


@pytest.fixture(scope='module')
def spectrum():
    kernel = np.asarray(w.dynamics.native_density()['counts'],dtype=float)/60
    values,vectors = np.linalg.eigh(kernel)
    values[-1] = 1 # exact stochastic ground eigenvalue, independent of roundoff
    assert values[0] > 0
    theta = np.arccos(values)
    weights = np.exp(-theta)
    return values,vectors,theta,weights/weights.sum()


@pytest.mark.parametrize('cap',(0,2,4,8))
def test_closed_reference_is_positive_and_has_uniform_wave_log_bound(cap,spectrum):
    eigen,vectors,theta,target = spectrum
    value = w.reference(cap)
    polynomial = np.array([sum(float(b)*x**n for n,b in enumerate(w.coefficients(cap))) for x in eigen])
    probability = polynomial/polynomial.sum()
    expected = (vectors*probability)@vectors.T
    np.testing.assert_allclose(np.asarray(value['density'],dtype=float),expected,atol=2e-15,rtol=3e-13)
    truncated_theta = -np.log(polynomial/float(sum(w.coefficients(cap))))
    allowance = w.cap_bounds(cap)
    assert min(theta-truncated_theta) > -2e-14
    assert max(theta-truncated_theta) <= float(allowance['wave_modular_generator_norm_error_upper'])
    assert np.abs(target-probability).sum()/2 <= float(allowance['reconstructed_reference_trace_distance_upper'])
    assert sp.trace(value['density']) == 1


def test_zero_mode_constant_and_positive_tail_bound_do_not_require_a_source_gap():
    cap = 8
    b = w.coefficients(cap)
    bound = w.cap_bounds(cap)
    for x in (0.,.01,.25,.5,.9,1.):
        polynomial = sum(float(value)*x**n for n,value in enumerate(b))
        normalized = polynomial/float(sum(b))
        exact = np.exp(-np.arccos(x))
        assert normalized >= exact-2e-16
        assert np.log(normalized/exact) <= float(bound['wave_modular_generator_norm_error_upper'])
    # At x=1 both generators vanish, while x=0 attains the uniform log bound.
    assert sum(b)/sum(b) == 1


def test_return_conditioning_favors_short_lengths_but_does_not_make_mean_finite(spectrum):
    eigen,_,_,_ = spectrum
    full_return = np.exp(-np.arccos(eigen)).mean()
    for cap in (0,2,8):
        polynomial = sum(float(b)*eigen**n for n,b in enumerate(w.coefficients(cap)))
        retained = np.exp(-np.pi/2)*polynomial.mean()
        tail = 1-retained/full_return
        assert 0 <= tail <= float(w.cap_bounds(cap)['conditional_unforced_return_history_TV_upper'])
    # Every power keeps the constant eigenvector: return chance >=1/66.
    for n in (0,1,8,64,512):
        assert np.mean(eigen**n) >= 1/66
    assert w.cap_bounds(512)['ideal_unforced_return_conditioned_expected_length_is_infinite']


def test_a_native_return_prepares_a_different_state_from_the_reconstructed_reference():
    value = w.reference(2)
    assert value['density'] != value['native_return_state']
    assert any(value['density'][i,j] for i in range(66) for j in range(66) if i!=j)
    assert all(value['native_return_state'][i,j] == 0 for i in range(66) for j in range(66) if i!=j)
    assert w.reference(0)['density'] == w.reference(0)['native_return_state']


def test_cap_errors_decrease_but_native_work_and_sampler_resources_grow():
    values = [w.cap_bounds(h) for h in (8,16,64,256,512)]
    assert all(a['complete_implemented_cap_probability'] > b['complete_implemented_cap_probability']
               for a,b in zip(values,values[1:]))
    assert all(a['wave_modular_generator_norm_error_upper'] > b['wave_modular_generator_norm_error_upper']
               for a,b in zip(values,values[1:]))
    assert values[-1]['wave_modular_generator_norm_error_upper'] < F(36,1000)
    assert values[-1]['reconstructed_reference_trace_distance_upper'] < F(353,10000)
    assert values[-1]['expected_executed_native_attempts_interval'][0] > 17
    assert values[-1]['uniform_categorical_denominator_bits'] > 64


def test_explicit_cap_growth_controls_the_rescaled_modular_generator():
    for h in (1,8,64,512):
        tail = w.cap_bounds(h)['ideal_length_tail_interval'][1]
        assert tail*tail <= F(16,h+1)
    for level in (0,1,4,20):
        value = w.refinement_budget(level)
        scale = 4**level
        delta = value['reference_trace_distance_upper']
        assert F(16,value['sufficient_native_length_cap']+1) < delta**2
        assert delta/(1-delta) <= value['unscaled_modular_generator_error_upper']
        assert value['scaled_wave_generator_error_upper'] == F(26,scale)
        assert not value['cap_sampler_at_this_budget_executed']
    # A fixed cap gives a fixed nonzero all-source log error, so cannot
    # silently be promoted to a vanishing bound after multiplication bys.
    assert 4**20*w.cap_bounds(512)['wave_modular_generator_norm_error_upper'] > 1


def test_normalization_cache_cannot_be_corrupted_through_returned_Fractions():
    before = w.normalization()['implemented_c0']
    exposed = w.normalization()
    exposed['implemented_c0']._numerator = 0
    exposed['c0_interval'][0]._denominator = 1
    assert w.normalization()['implemented_c0'] == before
    assert w.normalization()['c0_interval'][0] < F(21,100)


@pytest.mark.parametrize('value',(True,-1,1.,'2',513))
def test_invalid_length_domain_rejected(value):
    with pytest.raises(ValueError):w.coefficients(value)
    with pytest.raises(ValueError):w.sampler(value)
    with pytest.raises(ValueError):w.cap_bounds(value)


def test_receipt_replays_and_rejects_coherent_preparation_promotion():
    expected = json.loads(Path(__file__).with_name('wave_reference_results.json').read_text())
    actual = w.run(); compare_receipt(expected,actual)
    changed = copy.deepcopy(actual)
    changed['scope']['reconstructed_density_physically_prepared_on_native_working_source'] = True
    with pytest.raises(ValueError):compare_receipt(changed,actual)


def test_fresh_runtime_source_dependencies_are_pinned():
    command='''
import sys
from pathlib import Path
from model.observer_time_fields import wave_reference as w
w.run()
missing=[]
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name is None:continue
    try:path=Path(name).resolve().relative_to(w.ROOT).as_posix()
    except ValueError:continue
    if path.startswith('.'):continue
    if path.endswith('.py') and path not in w.SOURCE_PATHS:missing.append(path)
assert not missing,missing
'''
    subprocess.run([sys.executable,'-c',command],cwd=w.ROOT,check=True)

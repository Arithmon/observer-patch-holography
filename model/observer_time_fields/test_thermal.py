"""Exact native closure, thermal boundary and wrong-reference controls."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from model.observer.scalar_response.experiment import compare_receipt
from . import core
from . import thermal as t


@pytest.fixture(scope='module')
def source_transfer():
    data=core.source()
    a=np.array(data.remote_counts,dtype=float)/60
    b=np.array(data.incident_counts,dtype=float)/60
    d=np.eye(66)-a
    eigen,vectors=np.linalg.eigh(d)
    root=(vectors/np.sqrt(eigen))@vectors.T
    transfer=root@b@root
    values, modes=np.linalg.eigh(transfer)
    return transfer,values,modes


@pytest.mark.parametrize('cap',(1,2,3))
@pytest.mark.parametrize('feedback',(False,True))
def test_all_endpoint_current_histories_match_independent_oracle(cap,feedback):
    actual=t.loop_law(cap=cap,feedback=feedback)
    independent=t.native_oracle(cap=cap,feedback=feedback)
    assert actual==independent
    assert len(actual.records)==7
    assert sum(actual.probabilities().values())==1
    for table in actual.counts:
        assert len(table)==66 and all(len(row)==66 for row in table)
        assert all(type(x) is int and x>=0 for row in table for x in row)


def test_waiting_timeout_is_recorded_and_zero_current_is_not_discarded():
    law=t.loop_law(cap=1)
    p=law.probabilities()
    assert sum(value for (kind,_),value in p.items() if kind=='timeout')==F(91,100)
    assert sum(value for (kind,_),value in p.items() if kind=='return')==F(659,9900)
    # A one-attempt return is often the real zero-current native branch.
    one=t.loop_law((0,),cap=1)
    assert sum(value for (kind,_),value in one.probabilities().items() if kind=='return')>0


def test_capped_complete_error_retains_every_failure_and_timeout():
    exact=t.loop_law().probabilities()
    for cap in (1,2,3):
        censored=t.loop_law(cap=cap).probabilities()
        tv=sum((abs(exact.get(k,F())-censored.get(k,F())) for k in set(exact)|set(censored)),F())/2
        assert tv==t.cap_bound(2,cap)['complete_output_error']
        p=sum(v for (s,_),v in exact.items() if s=='return')
        q=sum(v for (s,_),v in censored.items() if s=='return')
        conditional=sum((abs(v/p-censored[key]/q) for key,v in exact.items() if key[0]=='return'),F())/2
        assert conditional<=(p-q)/p<=t.cap_bound(2,cap)['conditional_return_error']
    assert t.cap_bound(2,64)['conditional_return_error']<F(1,10**8)


def test_ground_reference_has_a_quantitative_all_positive_beta_KMS_obstruction(source_transfer):
    transfer,values,modes=source_transfer
    vector=np.array([state[0]-1/6 for state in core.source().states])/np.sqrt(66)
    spectral=(modes.T@vector)**2
    moments=t.open_moments()
    assert moments[:2]==(F(5,36),F(115,1188))
    assert moments[0]-moments[1]==F(25,594)
    for beta in (F(1,4),F(1,2),1,2,4):
        defect=float(np.sum(spectral*(1-values**float(beta))))
        assert defect+2e-14>=float(t.kms_defect_lower(beta))>0
    assert np.linalg.norm(transfer@np.ones(66)-np.ones(66))<1e-13
    # The original source density is I/66; Omega is a transfer amplitude.
    assert not np.allclose(np.ones((66,66))/66,np.eye(66)/66)


def test_reflection_positive_native_half_history_grams_have_exact_certificates():
    result=t.reflection_certificate()
    assert result['ground_vector_fixed']
    assert not result['native_full_configuration_kernel_is_classically_reversible']
    for gram in result['grams']:
        assert gram['principal_minor_count']==15
        assert gram['strictly_positive_minor_count']==15
        assert np.linalg.eigvalsh(np.array(gram['matrix'],dtype=float)).min()>0


@pytest.mark.parametrize('length',(1,2,3,4))
def test_actual_closed_words_match_positive_Gibbs_spectrum_and_KMS_wrap(length,source_transfer):
    _,values,modes=source_transfer
    result=t.thermal_cycle(length)
    rho=(modes*(values**length))@modes.T
    rho/=np.trace(rho)
    field=np.diag([state[0] for state in core.source().states])
    assert abs(np.trace(rho@field)-float(result['local_occupation_mean']))<2e-14
    assert result['two_point_cycle']==tuple(reversed(result['two_point_cycle']))
    assert result['mean_differs_from_open']
    eigen_field=modes.T@field@modes
    weights=values**length/sum(values**length)
    # Forward/reverse transition weights satisfy the finite KMS ratio.
    spectral=weights[:,None]*abs(eigen_field)**2
    ratio=(values[:,None]/values[None,:])**length
    np.testing.assert_allclose(spectral,ratio*spectral.T,atol=2e-15,rtol=2e-13)
    if length<=2:
        assert t.boundary_certificate(length)['exact_two_sided_transfer_inverse_checked']


def test_unconditioned_and_local_return_ensembles_are_not_the_thermal_circle():
    wrong=t.local_closure_control()
    assert wrong['local_return_probability']>wrong['full_return_probability']
    assert wrong['conditional_field_distance']>0
    assert t.thermal_cycle(2)['local_occupation_mean']!=F(1,6)
    ordinary=t.loop_law().probabilities(); feedback=t.loop_law(feedback=True).probabilities()
    assert sum(abs(ordinary[k]-feedback[k]) for k in ordinary)>0
    assert t.loop_law().denominator!=t.loop_law(feedback=True).denominator


def test_native_gap_controls_the_closed_Gibbs_ground_limit(source_transfer):
    _,values,_=source_transfer
    assert values[-2]<float(t.TRANSFER_RATE_UPPER)
    for length in (1,2,16,147):
        mass=sum(values[:-1]**length)
        assert mass/(1+mass)<=float(t.ground_limit_bound(length))+1e-14
    assert t.ground_limit_bound(146)>F(1,1000)
    assert t.ground_limit_bound(147)<=F(1,1000)


@pytest.mark.parametrize('kwargs',({'ports':()}, {'ports':(2,)}, {'ports':(True,1)},
                                  {'cap':True}, {'cap':0}, {'cap':1.0}, {'feedback':1}))
def test_bad_protocols_reject_even_after_legal_cache_entries(kwargs):
    t.loop_law((1,1),cap=1)
    with pytest.raises(ValueError):t.loop_law(**kwargs)


@pytest.mark.parametrize('value',(True,0,-1,1.0,'1'))
def test_inverse_temperature_must_be_positive_exact(value):
    with pytest.raises(ValueError):t.kms_defect_lower(value)


@pytest.fixture(scope='module')
def produced():return t.run()


def test_receipt_recomputes_exact_science_and_rejects_changed_claim(produced):
    compare_receipt(json.loads(json.dumps(produced)),produced)
    hostile=copy.deepcopy(produced)
    hostile['open_ground_state']['finite_beta_KMS_on_the_dynamical_field_algebra']=True
    with pytest.raises((ValueError,AssertionError)):compare_receipt(hostile,produced)
    assert produced['scope']['closure_failure_and_wait_timeouts_retained']
    assert not produced['scope']['local_records_uniquely_identify_the_full_Gibbs_density']


def test_fresh_runtime_imports_are_covered_by_source_pins():
    code="""from pathlib import Path
import sys
from model.observer_time_fields import thermal as t
t.run()
from model.observer.scalar_response.experiment import compare_receipt
missing=[]
for module in tuple(sys.modules.values()):
 p=getattr(module,'__file__',None)
 if p and Path(p).resolve().is_relative_to(t.ROOT):
  name=Path(p).resolve().relative_to(t.ROOT).as_posix()
  if name.endswith('.py') and name not in t.SOURCE_PATHS:missing.append(name)
assert not missing, sorted(set(missing))
"""
    subprocess.run([sys.executable,'-c',code],check=True,capture_output=True,text=True)

@pytest.mark.parametrize('length',(1,2,3))
def test_same_native_three_class_Gibbs_reference_has_a_complete_herald(length):
    report=t.class_cycle(length)
    law=t.class_loop_law(length)
    assert report['initial_filter_success']==F(3,11)
    assert report['initial_filter_failure']==F(8,11)
    assert sum(law.values())==1
    assert law==t.class_native_oracle(length,selection='first')
    assert law==t.class_native_oracle(length,selection='last')
    assert report['combined_return_probability']==report['partition']/11
    assert report['missing_herald_conditional_distance']>0
    assert report['combined_return_probability']>=F(1,11)
    if length==1:
        assert report['partition']==F(8,3)
        assert report['combined_return_probability']==F(8,33)
        assert report['conditional_Gibbs_class_diagonal']==(F(7,20),F(27,80),F(5,16))


def test_class_reference_has_the_same_exact_spectrum_as_the_compiled_probe():
    import sympy as sp
    _, source=t.class_source()
    matrix=sp.Matrix(source['matrix'])
    x=sp.Symbol('x')
    assert sp.expand(matrix.charpoly(x).as_expr()-(x-1)*(x*x-sp.Rational(5,3)*x+sp.Rational(311,450)))==0
    pi=source['stationary']
    assert all(pi[i]*source['matrix'][i][j]==pi[j]*source['matrix'][j][i] for i in range(3) for j in range(3))

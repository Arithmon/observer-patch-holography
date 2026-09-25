import importlib.util
from pathlib import Path
import numpy as np
import pytest
import sympy as s

path=Path(__file__).with_name('run.py')
spec=importlib.util.spec_from_file_location('incident_event_run',path)
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)

@pytest.mark.parametrize('n,raised,edges',[
    (2,1,[(0,1,1)]),
    (4,2,[(0,1,1),(1,2,1),(2,3,1),(3,0,1)]),
    (4,2,[(0,1,1),(1,2,2),(2,3,3)]),
    (5,1,[(0,1,2),(1,2,3),(2,3,1),(3,4,4)])])
def test_exact_bracket_and_compensator(n,raised,edges):
    m=r.model(n,raised,edges)
    mean=s.zeros(len(m['states']),n)
    for i,j,rate,reward in m['events']:
        for a in range(n):mean[i,a]+=rate*reward[a]
    assert mean==m['Z']
    assert m['Q']*m['H']==-m['Z']
    assert m['bracket']==m['closed']
    assert min(np.linalg.eigvalsh(np.array(m['bracket'],float)))>-1e-12


def test_shot_noise_is_not_independent_of_future_load():
    m=r.model(4,2,[(0,1,1),(1,2,1),(2,3,1),(3,0,1)])
    independent_guess=4*m['kappa']*m['K']+m['J']
    assert m['bracket']!=independent_guess
    assert m['bracket']-independent_guess==-m['kappa']*(m['Pi']*m['Dinv']+m['Dinv']*m['Pi'])


def test_two_ports_iid_marked_attempts_exact_variance_all_times():
    m=r.model(2,1,[(0,1,1)])
    v=np.array([1.,-1.])/np.sqrt(2)
    for time in [.01,.3,5.,100.]:
        assert r.moment_variance(m,v,time)==pytest.approx(time/2,rel=1e-11)
    # A continuous occupation record on the same chain has twice this limit.
    assert float(v@np.array(4*m['kappa']*m['K'],float)@v)==pytest.approx(1.)


def test_irregular_graph_tilted_generator_converges_to_bracket():
    m=r.model(4,2,[(0,1,1),(1,2,2),(2,3,3)])
    eig,vec=np.linalg.eigh(np.array(m['L'],float));v=vec[:,1]
    target=float(v@np.array(m['bracket'],float)@v)
    errors=[]
    for u in [10,100,1000]:
        time=2*u/eig[1];value=r.moment_variance(m,v,time)/time
        errors.append(abs(value-target))
    assert errors[1]<.101*errors[0]
    assert errors[2]<.101*errors[1]
    assert errors[2]/target<.001


def test_regular_graph_local_correction_formula():
    n=4;d=2
    m=r.model(n,2,[(0,1,1),(1,2,1),(2,3,1),(3,0,1)])
    predicted=4*m['kappa']*m['K']-m['kappa']*s.Rational(n+2,n*d)*m['Pi']+m['kappa']*m['L']/(n*d*d)
    assert m['Pi']*m['bracket']*m['Pi']==predicted


def test_graph_domain_guards():
    with pytest.raises(ValueError):r.model(4,2,[(0,1,1),(2,3,1)])
    with pytest.raises(ValueError):r.model(2,1,[(0,1,0)])
    with pytest.raises(ValueError):r.model(2,0,[(0,1,1)])


def test_precentered_marks_differ_from_raw_marks_minus_mean_time():
    m=r.model(2,1,[(0,1,1)])
    centered_total_bracket=0
    raw_total_bracket=0
    for i,j,rate,reward in m['events']:
        centered=sum(reward)
        raw=centered+1  # two endpoints, each background 1/2, unit attempt rates
        centered_total_bracket+=rate*centered**2/len(m['states'])
        raw_total_bracket+=rate*raw**2/len(m['states'])
    assert centered_total_bracket==0
    assert raw_total_bracket==1  # total raw count minus elapsed time is Poisson noise

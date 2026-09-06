"""Falsifiers and independent arithmetic/action controls for the time tube."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

def module(name):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(name+'.py'))
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


producer=module('whitney_charged_enclosure')
verifier=module('verify_whitney_charged_enclosure')


@pytest.fixture(scope="module")
def packet():
    return verifier.load()


def test_complete_independent_replay(packet):
    result=verifier.verify(packet)
    assert result['accepted'] is True
    assert result['whole_time_interval'] is True
    assert Q(result['exact_diagnostics']['uniform_error_replay_upper'])<=Q(1,10**20)
    assert Q(result['exact_diagnostics']['historical_error_replay_upper'])<=Q(1,10**10)
    assert result['physical_clock_selected'] is False


@pytest.mark.parametrize('pair',[(-7,-3),(-2,5),(0,7),(1,1),(3,8)])
def test_directed_interval_arithmetic(pair):
    a,b=pair
    for c,d in [(-9,-2),(-1,4),(2,7)]:
        x=producer.Interval(a*producer.SCALE,b*producer.SCALE,raw=True)
        y=producer.Interval(c*producer.SCALE,d*producer.SCALE,raw=True)
        u=verifier.Box(x.lo,x.hi);v=verifier.Box(y.lo,y.hi)
        for px,vx in [(x+y,u+v),(x-y,u-v),(x*y,u*v)]:
            assert (px.lo,px.hi)==(vx.low,vx.high)
        product=u*v
        for z in [Q(a*c),Q(a*d),Q(b*c),Q(b*d)]:
            assert Q(product.low,verifier.UNIT)<=z<=Q(product.high,verifier.UNIT)
        if c*d>0:
            inverse=v.inverse()
            assert Q(inverse.low,verifier.UNIT)<=Q(1,d)<=Q(inverse.high,verifier.UNIT)
            assert Q(inverse.low,verifier.UNIT)<=Q(1,c)<=Q(inverse.high,verifier.UNIT)


def test_non_dyadic_rounding_and_zero_division():
    for q in [Q(1,3),Q(-1,3),Q(1,10**90),Q(-1,10**90)]:
        x=verifier.Box(q)
        assert Q(x.low,verifier.UNIT)<=q<=Q(x.high,verifier.UNIT)
    with pytest.raises(ValueError,match='zero'):
        verifier.Box(-1,1).inverse()


@pytest.mark.parametrize('j',range(5))
def test_exact_simplex_algebra_reconstructs_rational_rhs(j):
    y=[Q((-1)**(j+k)*(j+k+1),17) for k in range(9)]
    assert producer.vector_field(y,Q(7,16))==verifier.field(y,Q(7,16))


def test_rational_rhs_matches_original_full_action_at_distinct_times():
    old=json.loads((verifier.ROOT/verifier.PARENT).read_bytes())
    mass=np.kron(np.array([[.2,.3],[.3,1.2]]),np.eye(2))
    e=.25
    r=6/(7+3*np.sqrt(5))
    def cross(c,b):
        return e/10*np.array([c.imag+b.imag,-c.real-b.real,c.imag+2*b.imag,-c.real-2*b.real])
    for index in [0,17,48,80]:
        s=old['samples'][index]
        a,cr,ci,br,bi=s['q_reduced']
        av,cvr,cvi,bvr,bvi=s['v_reduced']
        aa,car,cai,bar,bai=s['acceleration_reduced']
        c,b,cv=cr+1j*ci,br+1j*bi,cvr+1j*cvi
        rotation=np.exp(1j*e*a)
        big=rotation*c
        bigv=rotation*(cv+1j*e*av*c)
        biga=rotation*(car+1j*cai+2j*e*av*cv+1j*e*aa*c-e*e*av*av*c)
        xdot=np.array([bigv.real,bigv.imag,bvr,bvi])
        momentum=mass@xdot+cross(big,b)*av
        pdot=mass@np.array([biga.real,biga.imag,bar,bai])+cross(bigv,bvr+1j*bvi)*av+cross(big,b)*aa
        actual=np.array(producer.vector_field([a,big.real,big.imag,br,bi,*momentum],r),float)
        assert max(abs(actual-np.r_[av,xdot,pdot]))<1e-12


@pytest.mark.parametrize('mutation',[
    'scope','quantum','physical_clock','missing_step','wrong_order','bool_steps',
    'wrong_step','wrong_dimension','missing_pin','bad_pin','false_bound','unknown_field',
    'wrong_initial','collapsed_picard','collapsed_endpoint','polynomial','leading_zero',
    'negative_zero','float_endpoint','oversized_integer','reversed_interval','missing_polynomial'])
def test_semantic_mutations_fail(packet,mutation):
    p=deepcopy(packet)
    if mutation=='scope':p['scope']='PHYSICAL_COMPLETION'
    elif mutation=='quantum':p['interpretation']['quantum_history']=True
    elif mutation=='physical_clock':p['interpretation']['physical_clock_selected']=True
    elif mutation=='missing_step':p['states'].pop()
    elif mutation=='wrong_order':p['method']['taylor_order']=31
    elif mutation=='bool_steps':p['method']['steps']=True
    elif mutation=='wrong_step':p['method']['step']='1/20'
    elif mutation=='wrong_dimension':p['method']['coordinate_order'].pop()
    elif mutation=='missing_pin':p['source_pins'].pop(next(iter(p['source_pins'])))
    elif mutation=='bad_pin':p['source_pins'][next(iter(p['source_pins']))]='0'*64
    elif mutation=='false_bound':p['bounds']['uniform_canonical_polynomial_error_upper']='0'
    elif mutation=='unknown_field':p['measured_history']=[1]
    elif mutation=='wrong_initial':p['states'][0]['state'][0]=['1','1']
    elif mutation=='collapsed_picard':p['states'][0]['picard_box']=deepcopy(p['states'][0]['state'])
    elif mutation=='collapsed_endpoint':
        pair=p['states'][1]['state'][0];mid=(int(pair[0])+int(pair[1]))//2
        p['states'][1]['state'][0]=[str(mid),str(mid)]
    elif mutation=='polynomial':p['states'][0]['polynomial'][0][0]=str(1<<80)
    elif mutation=='leading_zero':p['states'][0]['picard_box'][0][0]='01'
    elif mutation=='negative_zero':p['states'][0]['picard_box'][0][0]='-0'
    elif mutation=='float_endpoint':p['states'][0]['picard_box'][0][0]=0.0
    elif mutation=='oversized_integer':p['states'][0]['picard_box'][0][0]='1'*101
    elif mutation=='reversed_interval':p['states'][0]['picard_box'][0]=['2','1']
    elif mutation=='missing_polynomial':p['states'][0]['polynomial'].pop()
    with pytest.raises(ValueError):verifier.verify(p)


@pytest.mark.parametrize('text',[
    '{"x":1,"x":2}', '{"nested":{"x":1,"x":2}}',
    '{"x":NaN}', '{"x":Infinity}', '{"x":1e10000}', '{"x":0.0}'])
def test_strict_json_loader(tmp_path,text):
    p=tmp_path/'bad.json';p.write_text(text,encoding='utf-8')
    with pytest.raises(ValueError):verifier.load(p)


def test_fresh_source_byte_mutation_fails(packet,tmp_path,monkeypatch):
    for rel in verifier.PINS:
        p=tmp_path/rel;p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes((verifier.ROOT/rel).read_bytes())
    altered=tmp_path/'paper/tex_fragments/WHITNEY_CHARGED_MATTER.tex'
    altered.write_bytes(altered.read_bytes()+b'\n% changed source\n')
    monkeypatch.setattr(verifier,'ROOT',tmp_path)
    with pytest.raises(ValueError,match='stale source pin'):verifier.verify(packet)


def test_historical_pin_remains_immutable():
    data=(verifier.ROOT/verifier.PARENT).read_bytes()
    assert hashlib.sha256(data).hexdigest()==verifier.PARENT_SHA256


def test_standalone_verifier_file_spec_import():
    spec=importlib.util.spec_from_file_location('isolated_charged_enclosure',Path(verifier.__file__))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    assert module.ALGEBRA['inverse']==((Q(8),Q(-2)),(Q(-2),Q(4,3)))


def test_continuous_error_detects_interior_only_polynomial_corruption(packet):
    # Add c*u*(h-u): zero at both endpoints, but too large in the interior.
    p=deepcopy(packet)
    coeff=p['states'][0]['polynomial'][0]
    amplitude=40*(1<<85)
    coeff[1]=str(int(coeff[1])+amplitude//40)
    coeff[2]=str(int(coeff[2])-amplitude)
    with pytest.raises(ValueError,match='continuous-time'):verifier.verify(p)

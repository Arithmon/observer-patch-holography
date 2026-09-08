#!/usr/bin/env python3
"""Independent analytic-derivative replay of the implicit scale certificate.

No producer or interval-dual implementation is imported. Interval images of
successive boxes retain the exact fixed point by the contraction theorem.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
from mpmath import iv

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OUTPUT=HERE.parent/'certificates/implicit_scale_attachment_receipt.json'
P_DEC='1.630968209403959324879279847782648941'
A_DEC=['0.041123336195630494','0.041125336195630496']


def _object(pairs):
    result={}
    for k,v in pairs:
        if k in result:
            raise ValueError('duplicate JSON key')
        result[k]=v
    return result


def load(path=OUTPUT):
    return json.loads(Path(path).read_text(encoding='utf-8'),object_pairs_hook=_object,
                      parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)),
                      parse_float=lambda s: (_ for _ in ()).throw(ValueError('unexpected float: '+s)))


def ingredients(a,y):
    P=iv.mpf(P_DEC);pi=iv.pi
    MU=iv.exp(-2*pi)*iv.exp(iv.log(P)/6)
    v=iv.exp(-iv.log(P)/2-pi/(2*a))
    logratio=iv.log(MU/y)
    alphas=[]
    for b in [iv.mpf(33)/5,iv.mpf(1),iv.mpf(-3)]:
        den=1/a+b*logratio/(2*pi)
        if not den.a>0:
            raise ValueError('nonpositive running denominator')
        alphas.append(1/den)
    a1,a2,_=alphas
    S=a2+iv.mpf(3)*a1/5
    B=a2*a2+iv.mpf(99)*a1*a1/25
    F=v*iv.sqrt(pi*S)
    Fy=F*B/(4*pi*y*S)
    Fa=F*(pi/(2*a*a)+(a2*a2+iv.mpf(3)*a1*a1/5)/(2*a*a*S))
    return F,Fy,Fa,alphas


def _residual(a,y,with_derivative):
    _,Fy,Fa,alphas=ingredients(a,y)
    yp=Fa/(1-Fy)
    result=-iv.mpf(P_DEC)/4
    derivative=iv.mpf(0)
    for kind,alpha,beta in [(2,alphas[1],iv.mpf(1)),(3,alphas[2],iv.mpf(-3))]:
        t=4*iv.pi**2*alpha
        dt=4*iv.pi**2*alpha**2*(1/a**2+beta*yp/(2*iv.pi*y))
        Z=S=C=CL=iv.mpf(0)
        labels=((n+1,iv.mpf(n*(n+2))/4) for n in range(129)) if kind==2 else (
            (iv.mpf((p+1)*(q+1)*(p+q+2))/2,iv.mpf(p*p+q*q+p*q+3*p+3*q)/3)
            for p in range(65) for q in range(65))
        for d,c in labels:
            ld=iv.log(d);weight=d*iv.exp(-t*c)
            Z+=weight;S+=weight*ld
            if with_derivative:
                C+=weight*c;CL+=weight*c*ld
        result+=S/Z
        if with_derivative:
            derivative-=dt*(CL/Z-C*S/(Z*Z))
    return result,derivative


def contained(inner,outer):
    if (not isinstance(outer,list) or len(outer)!=2 or
        any(type(s) is not str for s in outer)):
        raise ValueError('interval requires two decimal strings')
    # Exact parser rejects NaN/infinities before interval conversion.
    lo,hi=map(Fraction,outer)
    if lo>hi:
        raise ValueError('reversed interval')
    # Decimal bound values are themselves enclosed: use inward comparisons
    # so accepted computed boxes are inside the exact declared endpoints.
    l=iv.mpf(outer[0]);h=iv.mpf(outer[1])
    if not (l.b<=inner.a and inner.b<=h.a):
        raise ValueError('recomputed interval is outside the certificate')


def _verify(r):
    fields={'schema','pixel_decimal','coupling_interval','scale_rectangle','scale_image_outer',
            'scale_partial_mu_outer','contraction_bound','endpoint_residual_outer',
            'implicit_residual_derivative_outer','finite_representation_cutoffs',
            'interval_image_refinements','scope','source_pins'}
    if not isinstance(r,dict) or set(r)!=fields:
        raise ValueError('receipt key inventory')
    exact={'schema':'oph-hierarchy-implicit-scale-v1','pixel_decimal':P_DEC,'coupling_interval':A_DEC,
           'scale_rectangle':['7e-18','8e-18'],'contraction_bound':'1/100',
           'finite_representation_cutoffs':{'su2_max_n':128,'su3_max_p_and_q':64},
           'interval_image_refinements':{'whole_coupling_interval':8,'each_endpoint':16},
           'scope':{'true_implicit_scale':True,'unique_scale_on_declared_rectangle':True,
                    'unique_hierarchy_root_on_declared_interval':True,
                    'infinite_representation_sums':False,'physical_hierarchy_attachment':False,
                    'pixel_target_independent':False}}
    for key,value in exact.items():
        if json.dumps(r[key],sort_keys=True)!=json.dumps(value,sort_keys=True):
            raise ValueError('changed declared input or scope: '+key)
    sources=[HERE.parent/'computations/hierarchy_implicit_scale_witness.py',
             HERE.parent/'computations/hierarchy_interval_witness.py',Path(__file__).resolve(),
             HERE.parent/'test_interval_witness.py',ROOT/'code/interval_decimal.py']
    expected={p.relative_to(ROOT).as_posix():{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources}
    if json.dumps(r['source_pins'],sort_keys=True)!=json.dumps(expected,sort_keys=True):
        raise ValueError('source pin mismatch')
    a=iv.mpf(A_DEC);y=iv.mpf(r['scale_rectangle'])
    F,Fy,_,_=ingredients(a,y)
    contained(F,r['scale_image_outer']);contained(Fy,r['scale_partial_mu_outer'])
    if not (y.a<F.a and F.b<y.b and Fy.a>=0 and Fy.b<iv.mpf('0.01').a):
        raise ValueError('self-map/contraction failure')
    for _ in range(8):
        y=ingredients(a,y)[0]
    _,derivative=_residual(a,y,True)
    contained(derivative,r['implicit_residual_derivative_outer'])
    if not Fraction(r['implicit_residual_derivative_outer'][1])<0:
        raise ValueError('derivative sign not certified')
    endpoints=r['endpoint_residual_outer']
    if not isinstance(endpoints,list) or len(endpoints)!=2:
        raise ValueError('endpoint shape')
    for i,point in enumerate(A_DEC):
        aa=iv.mpf(point);yy=iv.mpf(r['scale_rectangle'])
        for _ in range(16):
            yy=ingredients(aa,yy)[0]
        value,_=_residual(aa,yy,False)
        contained(value,endpoints[i])
    if not (Fraction(endpoints[0][0])>0 and Fraction(endpoints[1][1])<0):
        raise ValueError('endpoint signs not certified')
    return {'true_implicit_scale':True,'unique_hierarchy_root':True,'contraction_bound':'1/100',
            'finite_representation_cutoffs':[128,64],'derivative_outer':r['implicit_residual_derivative_outer'],
            'infinite_representation_sums':False,'physical_hierarchy_attachment':False}


def verify(r):
    previous=iv.dps
    try:
        iv.dps=60
        return _verify(r)
    finally:
        iv.dps=previous


if __name__=='__main__':
    print(json.dumps(verify(load()),sort_keys=True))

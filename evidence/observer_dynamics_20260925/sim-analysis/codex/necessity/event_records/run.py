"""Incident attempted-event record: exact martingale bracket and independent moments."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
import sympy as sp

HERE=Path(__file__).resolve().parent

def model(n, raised, edges):
    if not 0 < raised < n: raise ValueError('Require a fluctuating fixed-total class')
    states=[tuple(int(i in c) for i in range(n)) for c in itertools.combinations(range(n),raised)]
    where={z:i for i,z in enumerate(states)}
    L=sp.zeros(n); d=[sp.Rational(0) for _ in range(n)]
    for a,b,w in edges:
        w=sp.Rational(w)
        if a==b or w<=0: raise ValueError('Positive rates on distinct endpoints required')
        L[a,a]+=w;L[b,b]+=w;L[a,b]-=w;L[b,a]-=w;d[a]+=w;d[b]+=w
    if L.rank()!=n-1: raise ValueError('Connected graph required')
    Pi=sp.eye(n)-sp.ones(n)/n
    K=(L+sp.ones(n)/n).inv()-sp.ones(n)/n
    Z=sp.Matrix(states)-sp.ones(len(states),n)*sp.Rational(raised,n)
    H=2*Z*K
    Q=sp.zeros(len(states)); J=sp.zeros(n); bracket=sp.zeros(n)
    events=[]
    for i,z in enumerate(states):
        for a,b,w in edges:
            reward=sp.zeros(n,1)
            for v in (a,b): reward[v]=Z[i,v]/d[v]
            for swap in (False,True):
                dest=list(z)
                if swap: dest[a],dest[b]=dest[b],dest[a]
                j=where[tuple(dest)]; rate=sp.Rational(w,2)
                Q[i,j]+=rate;Q[i,i]-=rate
                delta=H.row(j).T-H.row(i).T
                bracket+=rate*(reward+delta)*(reward+delta).T/len(states)
                J+=rate*reward*reward.T/len(states)
                events.append((i,j,rate,reward))
    kappa=sp.Rational(raised*(n-raised),n*(n-1)); Dinv=sp.diag(*[1/v for v in d])
    closed=4*kappa*K-kappa*(Pi*Dinv+Dinv*Pi)+J
    return dict(n=n,raised=raised,edges=edges,states=states,L=L,Pi=Pi,K=K,Z=Z,H=H,Q=Q,
                events=events,kappa=kappa,Dinv=Dinv,J=J,bracket=bracket,closed=closed)

def moment_variance(model, vector, time):
    """Different method: exact-size matrix exponential of tilted-generator moment ODE."""
    N=len(model['states']); q=np.array(model['Q'],float)
    R1=np.zeros_like(q);R2=np.zeros_like(q)
    vector=np.asarray(vector,float)
    for i,j,rate,reward in model['events']:
        value=float(vector @ np.array(reward,float).ravel())
        R1[i,j]+=float(rate)*value;R2[i,j]+=float(rate)*value*value
    A=np.zeros((3*N,3*N))
    for order in range(3): A[order*N:(order+1)*N,order*N:(order+1)*N]=q
    A[N:2*N,:N]=R1;A[2*N:,N:2*N]=2*R1;A[2*N:,:N]=R2
    initial=np.r_[np.ones(N),np.zeros(2*N)]
    evolved=la.expm(time*A)@initial
    return float(np.mean(evolved[2*N:])-np.mean(evolved[N:2*N])**2)

def encode(M): return [[str(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]

def execute():
    spec=json.loads((HERE/'spec.json').read_text());rows=[]
    for case in spec['cases']:
        m=model(case['ports'],case['raised'],case['edges'])
        assert m['bracket']==m['closed']
        assert m['Q']*m['H']==-m['Z']
        assert m['Q']*m['Z']==-m['Z']*m['L']/2
        assert m['Q']==m['Q'].T
        assert m['Z'].T*m['Z']/len(m['states'])==m['kappa']*m['Pi']
        values,vectors=la.eigh(np.array(m['L'],float));v=vectors[:,1]
        limiting=float(v @ np.array(m['bracket'],float) @ v)
        relaxed=[]
        for u in spec['windows_in_slowest_relaxation_times']:
            time=u*2/values[1]
            measured=moment_variance(m,v,time)/time
            relaxed.append({'aT':u,'time':time,'variance_per_time':measured,'limit':limiting,
                            'relative_error':measured/limiting-1})
        projected=m['Pi']*m['bracket']*m['Pi']
        rows.append({'name':case['name'],'configurations':len(m['states']),
                     'Q':encode(m['Q']),'L':encode(m['L']),
                     'incident_record_covariance':encode(m['bracket']),
                     'projected_incident_covariance':encode(projected),
                     'occupation_covariance':encode(4*m['kappa']*m['K']),
                     'reward_shot_noise':encode(m['J']),
                     'poisson_equation_exact':True,'bracket_identity_exact':True,
                     'slow_mode_moments':relaxed})
    return {'schema':'oph.native-incident-event-receipt.v1',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'spec_sha256':hashlib.sha256((HERE/'spec.json').read_bytes()).hexdigest(),
            'scope':'exact rational finite-chain identities and floating-point moment checks; no physical clock identification',
            'cases':rows}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=execute();target=HERE/'receipt.json'
    if args.verify:
        previous=json.loads(target.read_text())
        # Compare exact identities/hashes literally, numerical moments within platform roundoff.
        for old,new in zip(previous['cases'],result['cases'],strict=True):
            for a,b in zip(old.pop('slow_mode_moments'),new.pop('slow_mode_moments'),strict=True):
                for key in a: np.testing.assert_allclose(a[key],b[key],rtol=1e-8,atol=1e-8)
        assert previous==result
        print('Verified exact incident-event identities and independent finite-window moments')
    else:
        target.write_text(json.dumps(result,indent=2)+'\n')
        print('Wrote incident-event receipt for',len(result['cases']),'chains')

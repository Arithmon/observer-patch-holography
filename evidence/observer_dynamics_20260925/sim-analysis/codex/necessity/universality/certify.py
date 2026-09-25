"""Exact finite-chain certificates for weighted native repair readbacks.

Fractions are used for the generator, stationary law, pseudoinverses and
asymptotic covariance identities. Floating checks are explicitly separate.
"""
from __future__ import annotations

import ast
from fractions import Fraction as F
import hashlib
import itertools
import json
from math import isqrt
from pathlib import Path
import subprocess

import numpy as np
import scipy.linalg as la
from scipy.integrate import quad_vec

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[3]


def zeros(rows,cols):
    return [[F(0) for _ in range(cols)] for _ in range(rows)]


def eye(size):
    return [[F(i==j) for j in range(size)] for i in range(size)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def add(a,b):
    return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def scale(c,a):
    return [[c*x for x in row] for row in a]


def inverse(a):
    n=len(a)
    block=[list(row)+identity for row,identity in zip(a,eye(n))]
    for col in range(n):
        pivot=next(row for row in range(col,n) if block[row][col])
        block[col],block[pivot]=block[pivot],block[col]
        norm=block[col][col]
        block[col]=[x/norm for x in block[col]]
        for row in range(n):
            if row!=col:
                multiple=block[row][col]
                block[row]=[x-multiple*y for x,y in zip(block[row],block[col])]
    return [row[n:] for row in block]


def pseudoinverse_connected_laplacian(a):
    n=len(a)
    constant=[[F(1,n) for _ in range(n)] for _ in range(n)]
    result=add(inverse(add(a,constant)),scale(F(-1),constant))
    projection=add(eye(n),scale(F(-1),constant))
    assert multiply(a,result)==projection
    assert multiply(result,a)==projection
    assert result==transpose(result)
    return result


def congruence(b,a):
    return multiply(multiply(b,a),transpose(b))


def rational_sqrt(value):
    a,b=isqrt(value.numerator),isqrt(value.denominator)
    if a*a!=value.numerator or b*b!=value.denominator:
        raise ValueError("Certificate incidence weights must have rational square roots")
    return F(a,b)


def production_rule(commit):
    path="oph_exact/carrier.py"
    source=subprocess.check_output(["git","-C",str(WORKSPACE/"oph-physics-sim"),"show",f"{commit}:{path}"])
    definition=next(x for x in ast.parse(source).body
                    if isinstance(x,ast.FunctionDef) and x.name=="integer_nearest_agreement")
    namespace={}
    exec(compile(ast.Module(body=[definition],type_ignores=[]),f"{commit}:{path}","exec"),namespace)
    return namespace[definition.name],{"commit":commit,"path":path,"sha256":hashlib.sha256(source).hexdigest()}


def configuration_data(n,raised,edges,rule=None):
    states=[tuple(int(i in selected) for i in range(n)) for selected in itertools.combinations(range(n),raised)]
    where={state:i for i,state in enumerate(states)}
    q=zeros(len(states),len(states))
    microscopic=zeros(len(states),len(states))
    for i,state in enumerate(states):
        for a,b,w in edges:
            changed=list(state)
            changed[a],changed[b]=changed[b],changed[a]
            q[i][where[tuple(changed)]]+=w/2
            q[i][i]-=w/2
            if rule is not None:
                for coin in [False,True]:
                    changed=list(state)
                    first,second=rule(7+state[a],7+state[b],ceiling_to_first=coin)
                    changed[a],changed[b]=first-7,second-7
                    microscopic[i][where[tuple(changed)]]+=w/2
                    microscopic[i][i]-=w/2
    if rule is not None:
        assert q==microscopic
    lap=zeros(n,n)
    incidence=zeros(len(edges),n)
    for index,(a,b,w) in enumerate(edges):
        lap[a][a]+=w; lap[b][b]+=w; lap[a][b]-=w; lap[b][a]-=w
        incidence[index][a]=rational_sqrt(w)
        incidence[index][b]=-rational_sqrt(w)
    centered=[[F(x)-F(raised,n) for x in state] for state in states]
    kappa=F(raised*(n-raised),n*(n-1))
    pi=add(eye(n),scale(F(-1,n),[[F(1) for _ in range(n)] for _ in range(n)]))
    assert q==transpose(q) and all(sum(row)==0 for row in q)
    assert multiply(q,centered)==scale(F(-1,2),multiply(centered,lap))
    assert scale(F(1,len(states)),multiply(transpose(centered),centered))==scale(kappa,pi)
    assert multiply(transpose(incidence),incidence)==lap
    return states,q,lap,incidence,centered,kappa,pi


def strings(a):
    return [[str(x) for x in row] for row in a]


def response_fraction(u):
    if u<1e-3:
        return u/2-u*u/6+u**3/24-u**4/120+u**5/720
    return 1+np.expm1(-u)/u


def finite_checks(q,lap,centered,kappa,records,windows):
    q=np.array(q,dtype=float)
    lap=np.array(lap,dtype=float)
    z=np.array(centered,dtype=float)
    eigen,vectors=la.eigh(lap)
    positive=eigen[1:]
    v=vectors[:,1:]
    answer=[]
    for name,b in records.items():
        b=np.array(b,dtype=float)
        observables=z@b.T
        limit=4*float(kappa)*(b@v/positive)@(b@v).T
        for window in windows:
            factors=np.array([response_fraction(lam*window/2) for lam in positive])
            expected=4*float(kappa)*(b@v*(factors/positive))@(b@v).T
            integrand=lambda t:2*(window-t)/window*(observables.T@la.expm(q*t)@observables)/len(q)
            actual,error=quad_vec(integrand,0,window,epsabs=1e-10,epsrel=1e-10)
            lower=response_fraction(positive[0]*window/2)*limit
            answer.append({"record":name,"window":window,
                           "quadrature_relative_matrix_error":float(la.norm(actual-expected)/la.norm(expected)),
                           "quadrature_error_estimate":float(error),
                           "lower_bound_min_eigenvalue":float(la.eigvalsh(expected-lower)[0]),
                           "upper_bound_min_eigenvalue":float(la.eigvalsh(limit-expected)[0])})
    # Independent moment ODE for one deterministic initial configuration.
    f=z@vectors[:,1]
    size=len(q)
    block=np.zeros((3*size,3*size))
    for degree in range(3):
        block[degree*size:(degree+1)*size,degree*size:(degree+1)*size]=q
        if degree:
            block[degree*size:(degree+1)*size,(degree-1)*size:degree*size]=degree*np.diag(f)
    initial=np.zeros(3*size); initial[:size]=1
    nonstationary=[]
    for window in windows:
        moments=la.expm(block*window)@initial
        mean=moments[size]/np.sqrt(window)
        variance=moments[2*size]/window-mean*mean
        a=positive[0]/2
        predicted_mean=f[0]*(-np.expm1(-a*window))/(a*np.sqrt(window))
        nonstationary.append({"window":window,"mean_moment_ode":float(mean),
                              "mean_graph_formula":float(predicted_mean),
                              "variance_deterministic_initial":float(variance),
                              "variance_stationary_formula":float(2*float(kappa)/a*response_fraction(a*window))})
    return answer,nonstationary


def graph_certificate(description,rule,spec):
    n,raised=description["ports"],description["raised"]
    edges=[(a,b,F(w)) for a,b,w in description["edges"]]
    states,q,lap,incidence,z,kappa,pi=configuration_data(n,raised,edges,rule)
    lp=pseudoinverse_connected_laplacian(lap)
    gp=pseudoinverse_connected_laplacian(scale(F(-1),q))
    degree=max(lap[i][i] for i in range(n))
    records={"load":eye(n),"drive":lap,"seam":incidence,
             "lazy_average":add(eye(n),scale(-F(1)/(2*degree),lap))}
    w_total=sum(w for _,_,w in edges)
    record_results={}
    for name,b in records.items():
        observable=multiply(z,transpose(b))
        config=scale(F(2,len(states)),multiply(multiply(transpose(observable),gp),observable))
        vertex=scale(4*kappa,congruence(b,lp))
        assert config==vertex
        static=scale(F(1,len(states)),multiply(transpose(observable),observable))
        # (I-P)^+ = W*(-G)^+ for the fixed weighted-attempt chain.
        attempt_kernel=add(scale(2*w_total,gp),scale(F(-1),eye(len(states))))
        discrete_config=scale(F(1)/(w_total*len(states)),congruence(transpose(observable),attempt_kernel))
        discrete_graph=add(vertex,scale(-kappa/w_total,congruence(b,pi)))
        assert discrete_config==discrete_graph
        record_results[name]={"B":strings(b),"stationary_covariance":strings(static),
                              "continuous_occupation_limit":strings(vertex),
                              "discrete_attempt_occupation_limit":strings(discrete_graph),
                              "configuration_and_graph_covariance_equal_exactly":True}
    assert [[F(x) for x in row] for row in record_results["drive"]["continuous_occupation_limit"]]==scale(4*kappa,lap)
    cut=congruence(incidence,lp)
    assert multiply(cut,cut)==cut and cut==transpose(cut)
    c=F(3,2)
    assert pseudoinverse_connected_laplacian(scale(c,lap))==scale(1/c,lp)
    # Successful state-change clock biases the configuration sampling measure.
    exit_rates=[-q[i][i] for i in range(len(q))]
    jump_weights=[rate/sum(exit_rates) for rate in exit_rates]
    jump=zeros(len(q),len(q))
    for i in range(len(q)):
        for j in range(len(q)):
            if i!=j:
                jump[i][j]=q[i][j]/exit_rates[i]
    assert multiply([jump_weights],jump)==[jump_weights]
    jump_mean=[sum(weight*F(state[j]) for weight,state in zip(jump_weights,states)) for j in range(n)]
    jump_centered=[[F(x)-jump_mean[j] for j,x in enumerate(state)] for state in states]
    weighted=[[weight*x for x in row] for weight,row in zip(jump_weights,jump_centered)]
    jump_cov=multiply(transpose(jump_centered),weighted)
    symmetry=None
    if description["name"]=="cycle5":
        perm=zeros(n,n)
        for i in range(n):perm[(i+1)%n][i]=1
        covariance=scale(4*kappa,lp)
        assert congruence(perm,lap)==lap and congruence(perm,covariance)==covariance
        gain=eye(n);gain[0][0]=2
        changed=congruence(gain,covariance)
        assert congruence(perm,changed)!=changed
        symmetry={"cycle_rotation_preserves_load_and_drive_laws":True,
                  "unequal_local_record_gain_breaks_covariance_symmetry":True}
    finite,nonstationary=finite_checks(q,lap,z,kappa,records,spec["numerical_checks"]["finite_windows"])
    return {"name":description["name"],"ports":n,"raised":raised,
            "edge_rates":description["edges"],"total_attempt_rate":str(w_total),"kappa":str(kappa),
            "states":states,"configuration_generator":strings(q),"laplacian":strings(lap),
            "laplacian_pseudoinverse":strings(lp),"records":record_results,
            "pinned_production_rule_equals_swap_generator":True,
            "uniform_rate_scaling_check_exact":True,"cut_projector_check_exact":True,
            "successful_event_clock":{"exit_rates":list(map(str,exit_rates)),
                "stationary_configuration_weights":list(map(str,jump_weights)),
                "stationary_load_mean":list(map(str,jump_mean)),"stationary_load_covariance":strings(jump_cov),
                "differs_from_attempt_clock_covariance":jump_cov!=scale(kappa,pi)},
            "symmetry":symmetry,"finite_window_checks":finite,"nonstationary_initial_checks":nonstationary}


def build():
    spec=json.loads((HERE/"spec.json").read_text())
    rule,pin=production_rule(spec["source_commit"])
    refinement=[]
    control=spec["refinement_countercontrol"]
    for size in control["sizes"]:
        gap=4*np.sin(np.pi/size)**2
        refinement.append({"cycle_size":size,"gap":float(gap),
            "fixed_window":control["fixed_window"],
            "fixed_window_gff_fraction":float(response_fraction(control["fixed_window"]*gap/2)),
            "scaled_window":float(control["window_times_gap"]/gap),
            "scaled_window_gff_fraction":float(response_fraction(control["window_times_gap"]/2))})
    return {"schema":"oph.weighted-native-readback-universality.receipt.v1",
            "inputs_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ["spec.json","certify.py"]},
            "production_rule_source":pin,
            "graphs":[graph_certificate(row,rule,spec) for row in spec["graphs"]],
            "refinement_countercontrol":refinement,
            "scope":spec["scope"]}


if __name__=="__main__":
    result=build()
    (HERE/"receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("Exact weighted readback certificates passed; wrote receipt.json")

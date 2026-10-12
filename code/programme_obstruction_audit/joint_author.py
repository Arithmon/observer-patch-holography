"""Bounded exact attachment controls for OPH issue 1078; standard library only.

Target 0510578e317de5da0d485d09792dc8bcad2e2f3c. No search/optimization.
Reads the target's response.json; all arithmetic below is rational or Q(sqrt(5)).
Stop after the fixed checks. These checks do not establish full A1 model membership.
"""
import argparse
import hashlib
import itertools
import json
import platform
from fractions import Fraction as F
from pathlib import Path

SHA = "0510578e317de5da0d485d09792dc8bcad2e2f3c"
ZERO = (F(0), F(0))
ONE = (F(1), F(0))
ROOT = (F(0), F(1))

def q(x):
    return (F(x), F(0))
def add(x, y):
    return (x[0]+y[0], x[1]+y[1])
def neg(x):
    return (-x[0], -x[1])
def sub(x,y):
    return add(x,neg(y))
def mul(x, y):
    return (x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def div(x,y):
    den=y[0]*y[0]-5*y[1]*y[1]
    assert den
    return mul(x,(y[0]/den,-y[1]/den))
def conj(x):
    return (x[0],-x[1])
def summ(xs):
    z=ZERO
    for x in xs:z=add(z,x)
    return z
def fmt(x):
    return f"{x[0]} + ({x[1]})sqrt(5)"
def mm(a,b):
    return [[summ(mul(x,y) for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def scale(a,s):
    return [[mul(x,s) for x in row] for row in a]
def order(p):
    r=list(range(len(p)))
    for k in range(1,61):
        r=[p[j] for j in r]
        if r==list(range(len(p))):return k
    raise AssertionError("permutation order above bound")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("response",type=Path)
    args=parser.parse_args()
    response=json.loads(args.response.read_text())
    checks=[]
    def check(name, condition):
        assert condition, name
        checks.append(name)

    # v recovered from the real skew part of the first source-response block.
    v=[]
    for blocks in response["generators"]:
        a=blocks[0]
        v.append([mul(q(2),(F(a[i][j][0]),F(a[i][j][1])))
                  for i,j in [(2,1),(0,2),(1,0)]])
    dot=lambda x,y:summ(mul(a,b) for a,b in zip(x,y))
    radius=dot(v[0],v[0])
    check("all twelve recovered vectors have common norm",all(dot(x,x)==radius for x in v))
    gram=[[div(dot(x,y),radius) for y in v] for x in v]
    adjacency=[[int(gram[i][j]==div(ROOT,q(5))) for j in range(12)] for i in range(12)]
    check("icosahedral graph has 30 edges and degree five",sum(map(sum,adjacency))==60 and all(sum(row)==5 for row in adjacency))
    check("slow Gram G squared equals 4G",mm(gram,gram)==scale(gram,q(4)))
    check("slow Gram trace twelve",summ(gram[i][i] for i in range(12))==q(12))
    lap=[[q(5*int(i==j)-adjacency[i][j]) for j in range(12)] for i in range(12)]
    check("slow Gram lies in 5-sqrt(5) band",mm(lap,gram)==scale(gram,sub(q(5),ROOT)))
    mirror=[[conj(x) for x in row] for row in gram]
    check("Galois Gram lies in 5+sqrt(5) band",mm(lap,mirror)==scale(mirror,add(q(5),ROOT)))
    p=next(p for p in response["closed_paths"] if order(p["port_action"])==5)
    traces=[summ((F(b[i][i][0]),F(b[i][i][1])) for i in range(3)) for b in p["observed_response"]]
    check("same order-five action has inequivalent triplet traces",traces[0]!=traces[1] and traces[1]==conj(traces[0]))
    check("order-five matrix action permutes recovered slow frame",all(
        [summ(mul((F(p["observed_response"][0][i][j][0]),F(p["observed_response"][0][i][j][1])),v[k][j]) for j in range(3)) for i in range(3)]==v[p["port_action"][k]]
        for k in range(12)))

    shell=list(itertools.combinations(range(12),2))
    check("binary shell dimension 66",len(shell)==66)
    odd=F(sum((0 in x)!=(1 in x) for x in shell),len(shell))
    kappa=odd/2
    bcoef=kappa/30
    nu=3*bcoef
    check("uniform native move fraction",kappa==F(5,33))
    check("uniform native response covariance",bcoef==F(1,198))
    check("triangle circulation variance",nu==F(1,66))
    ratio=mul(q(kappa/nu/F(60)),sub(q(5),ROOT))
    check("native common-rate-free ratio",ratio==div(sub(q(5),ROOT),q(6)))
    unit_scale=F(2)
    rescaled_nu=unit_scale**2*nu
    rescaled_ratio=mul(q(kappa/rescaled_nu/F(60)),sub(q(5),ROOT))
    check("dimensionless kappa ratio requires atomic current normalization",rescaled_ratio==div(ratio,q(4)))
    check("record Laplacian is independent of current units",3*unit_scale**2*bcoef/rescaled_nu==1)
    # Clock/source-law control: insert independently retained idle attempts.
    # With eta=1/2, the working transition has same eigenvectors and stationary
    # state; kappa, mean relaxation and circulation variance each gain eta.
    eta=F(1,2)
    idle_ratio=mul(q(eta*kappa/(eta*nu)*eta/F(60)),sub(q(5),ROOT))
    check("idle attempts change the kappa-weighted ratio",idle_ratio==mul(q(eta),ratio))

    # M66 occupation projectors fail A1's central orthogonal partition.
    x=(0,1);y=(1,2)
    check("occupation N0 has nonzero matrix-unit commutator",int(0 in x)-int(0 in y)==1)
    check("occupations are not pairwise orthogonal",0 in x and 1 in x)
    check("occupation sum is twice identity",all(sum(int(i in x) for i in range(12))==2 for x in shell))
    # Analytic centre proof, recorded as a witness rather than a fake computed
    # classification: commutation with Eii kills all off-diagonal entries;
    # commutation with Eij equates all diagonal entries, hence centre=C I.

    # Finite oscillator cutoff: [a,a*]=I-d|d-1><d-1|, trace zero.
    ccr={}
    for d in [2,3,4,6]:
        comm=[F(1)]*(d-1)+[F(1-d)]
        check(f"dimension {d} finite CCR boundary defect",sum(comm)==0 and comm!=[F(1)]*d)
        ccr[str(d)]=[str(x) for x in comm]

    # Exact three-setting clock control; use rational pairs for real/imaginary.
    clock=[]
    for re,im in [(F(1,2),F(0)),(F(0),F(-1,2)),(F(-1,2),F(0))]:
        rho=[[(F(1,2),F(0)),(re,im)],[(re,-im),(F(1,2),F(0))]]
        e=[[F(1,2),F(1,2)],[F(1,2),F(1,2)]]
        coherent=sum(e[j][i]*rho[i][j][0] for i in range(2) for j in range(2))
        dephased=sum(e[i][i]*rho[i][i][0] for i in range(2))
        clock.append([str(coherent),str(dephased)])
    check("complete dephasing erases the three-setting clock",clock==[['1','1/2'],['1/2','1/2'],['0','1/2']])
    out={"target_revision":SHA,"python":platform.python_version(),
         "budget":"fixed twelve-vector algebra; 60 source paths maximum; 66 binary states; cutoffs 2,3,4,6; no numerical optimization",
         "response_input":"data/response.json","response_sha256":hashlib.sha256(args.response.read_bytes()).hexdigest(),"check_count":len(checks),"checks":checks,
         "order_five_port_action":p["port_action"],"slow_trace":fmt(traces[0]),"su2_mirror_trace":fmt(traces[1]),
         "slow_laplacian":"5-sqrt(5)","su2_mirror_laplacian":"5+sqrt(5)",
         "kappa":str(kappa),"B_L_coefficient":str(bcoef),"nu":str(nu),
         "native_ratio":fmt(ratio),"idle_eta":str(eta),"idle_ratio":fmt(idle_ratio),
         "current_unit_scale":str(unit_scale),"rescaled_dimensionless_kappa_ratio":fmt(rescaled_ratio),
         "central_port_commutator":"[N0, |01><12|] = |01><12|",
         "central_port_analytic_proof":"Z(M66)=C I by Eii/Eij commutation; only nonzero central projection is I",
         "truncated_oscillator_commutator_diagonals":ccr,
         "coherent_and_dephased_clock_probabilities":clock,
         "scope":"attachment controls; no all-model OPH no-go; response source-selection and full model membership not certified"}
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()

"""Producer-free Fraction replay of the protected-record memory witness.

The kernel is reconstructed from conditional probability cells. The proposed
resolvent is checked by both Poisson inverse identities, not recomputed by
the producer's Gaussian elimination. Every correlation, tail and distribution
is independently iterated. This verifies this finite fixture; the general
finite-fibre argument is analytic, not a newly formalized Lean theorem.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUTPUT=HERE/"runtime"/"protected_memory_receipt.json"
SOURCES=("protected_memory.py","verify_protected_memory.py","test_protected_memory.py")
PARENTS=("Lean/Thermodynamics/GreenKubo.lean","Lean/Thermodynamics/FiniteConditionalRepair.lean",
         "code/thermodynamics/common_reference_obstruction/verify_common_reference_obstruction.py",
         "paper/tex_fragments/PROTECTED_RECORD_MEMORY.tex")


def need(value,label):
    if not value:raise ValueError(label)


def exact(actual,expected,label):
    need(type(actual) is type(expected),label+" type")
    if isinstance(expected,dict):
        need(set(actual)==set(expected),label+" fields")
        for k,v in expected.items():exact(actual[k],v,label+"."+k)
    elif isinstance(expected,list):
        need(len(actual)==len(expected),label+" length")
        for a,b in zip(actual,expected):exact(a,b,label+" entry")
    else:need(actual==expected,label+" value")


def load(path=OUTPUT):
    def pairs(rows):
        d={}
        for k,v in rows:
            need(k not in d,"duplicate JSON key");d[k]=v
        return d
    def invalid(x):raise ValueError("floating JSON values forbidden: "+x)
    return json.loads(Path(path).read_text(encoding="utf-8"),object_pairs_hook=pairs,
                      parse_constant=invalid,parse_float=invalid)


def rational(x):
    need(type(x) is str,"rational must be canonical string")
    try:f=F(x)
    except (ValueError,ZeroDivisionError) as exc:raise ValueError("invalid rational") from exc
    need(str(f)==x,"noncanonical rational")
    return f


def vector(xs,n):
    need(type(xs) is list and len(xs)==n,"vector shape")
    return [rational(x) for x in xs]


def matrix(rows,n,m):
    need(type(rows) is list and len(rows)==n,"matrix shape")
    return [vector(row,m) for row in rows]


def pin(path):
    raw=path.read_bytes()
    return {"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}


def act(A,v):
    return [sum((a*b for a,b in zip(row,v)),F(0)) for row in A]


def product(A,B):
    columns=list(zip(*B))
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in columns] for row in A]


def verify(p):
    need(type(p) is dict,"root mapping")
    exact(sorted(p),sorted(["schema","scope","model","reference","coordinate_kernels","equilibrium_projection",
          "transition","stationary_pair_coupling","two_step_minorization","poisson_resolvent","currents",
          "poisson_solutions","green_kubo","global_centering_counterexample","distribution_history",
          "source_pins","parent_pins"]),"receipt inventory")
    exact(p["schema"],"oph.protected_record_memory.v1","schema")
    exact(p["source_pins"],{x:pin(HERE/x) for x in SOURCES},"source pins")
    exact(p["parent_pins"],{x:pin(ROOT/x) for x in PARENTS},"parent pins")
    exact(p["scope"],{"new_supplied_local_stochastic_law":True,"same_reference_equilibrium_and_dynamics":True,
          "nonconstant_record_preserved":True,"exact_rational_replay":True,"native_source_attachment":False,
          "historical_obstruction_overturned":False,"global_ergodicity":False,"physical_clock_or_conductivity":False,
          "continuum_limit":False,"empirical_evidence":False,"new_Lean_formalization":False},"scope")
    states=[[i//4,(i//2)%2,i%2] for i in range(8)]
    exact(p["model"],{"states":states,"protected_coordinate":"c","fibre_masses":["1/3","2/3"],
          "conditional_reference_weights":[[1,2,3,4],[4,1,2,3]],"reference_denominator":10,
          "step_unit":"one declared stochastic update; no laboratory calibration",
          "operations":[{"name":"hold","probability":"1/2","reads":[],"writes":[]},
                        {"name":"resample_u","probability":"1/4","reads":["c","v"],"writes":["u"]},
                        {"name":"resample_v","probability":"1/4","reads":["c","u"],"writes":["v"]}]},"declared source law")
    pi=vector(p["reference"],8)
    need(pi==[F(1,30),F(2,30),F(3,30),F(4,30),F(8,30),F(2,30),F(4,30),F(6,30)],"reference law")
    need(sum(pi)==1 and all(x>0 for x in pi),"faithful reference")
    kernels={}
    # Group each read-port value; this differs from pair-by-pair construction.
    for name,coordinate in (("u",1),("v",2),("P",None)):
        groups={}
        for i,x in enumerate(states):
            key=tuple(x[k] for k in range(3) if k!=coordinate) if coordinate else (x[0],)
            groups.setdefault(key,[]).append(i)
        K=[[F(0)]*8 for _ in range(8)]
        for cell in groups.values():
            mass=sum((pi[j] for j in cell),F(0))
            for i in cell:
                for j in cell:K[i][j]=pi[j]/mass
        kernels[name]=K
    exact(sorted(p["coordinate_kernels"]),["u","v"],"coordinate kernel names")
    for name in ("u","v"):
        need(matrix(p["coordinate_kernels"][name],8,8)==kernels[name],"conditional "+name+" kernel")
    P=kernels["P"];need(matrix(p["equilibrium_projection"],8,8)==P,"equilibrium projection")
    I=[[F(i==j) for j in range(8)] for i in range(8)]
    T=matrix(p["transition"],8,8)
    need(T==[[I[i][j]/2+(kernels["u"][i][j]+kernels["v"][i][j])/4 for j in range(8)] for i in range(8)],"lazy local transition")
    for K in (P,T,*[kernels[n] for n in ("u","v")]):
        need(all(sum(row)==1 and min(row)>=0 for row in K),"stochastic rows")
        need(all(pi[i]*K[i][j]==pi[j]*K[j][i] for i in range(8) for j in range(8)),"detailed balance")
        need(all(not K[i][j] for i in range(8) for j in range(8) if states[i][0]!=states[j][0]),"protected-record support")
    need(product(P,P)==P and product(T,P)==P and product(P,T)==P,"equilibrium/dynamics composition")
    for name in ("u","v"):need(product(kernels[name],kernels[name])==kernels[name],"coordinate projection")
    T2=product(T,T);need(T2!=T,"nonidempotent dynamics")
    need(matrix(p["stationary_pair_coupling"],8,8)==[[pi[i]*T[i][j] for j in range(8)] for i in range(8)],"stationary stochastic coupling")
    block=p["two_step_minorization"]
    exact(sorted(block),sorted(["epsilon","remainder","gap_lower","L2_decay"]),"minorization fields")
    epsilon=rational(block["epsilon"])
    need(epsilon==min(T2[i][j]/P[i][j] for i in range(8) for j in range(8) if P[i][j]) and 0<epsilon<1,"strict maximal two-step minorization")
    remainder=matrix(block["remainder"],8,8)
    need(remainder==[[T2[i][j]-epsilon*P[i][j] for j in range(8)] for i in range(8)],"minorization remainder")
    need(all(min(row)>=0 and sum(row)==1-epsilon for row in remainder),"stationary Markov contraction remainder")
    need(rational(block["gap_lower"])==epsilon/2,"spectral-gap lower bound")
    exact(block["L2_decay"],"(1-epsilon)^floor(n/2) on ker(P)","decay domain")
    R=matrix(p["poisson_resolvent"],8,8)
    D=[[I[i][j]-T[i][j] for j in range(8)] for i in range(8)]
    centre=[[I[i][j]-P[i][j] for j in range(8)] for i in range(8)]
    zero=[[F(0)]*8 for _ in range(8)]
    need(product(D,R)==centre and product(R,D)==centre,"two-sided Poisson inverse")
    need(product(P,R)==zero and product(R,P)==zero,"fibre-centred Poisson range")
    need(all(pi[i]*R[i][j]==pi[j]*R[j][i] for i in range(8) for j in range(8)),"selfadjoint resolvent")
    currents=matrix(p["currents"],2,8);solutions=matrix(p["poisson_solutions"],2,8)
    expected=[]
    for k in (1,2):
        raw=[F(s[k]) for s in states]
        expected.append([a-b for a,b in zip(raw,act(P,raw))])
    need(currents==expected,"conditional rather than global current centering")
    need(solutions==[act(R,f) for f in currents],"Poisson current solutions")
    pair=lambda f,g:sum((a*b*w for a,b,w in zip(f,g,pi)),F(0))
    gram=[[pair(f,g) for g in solutions] for f in currents]
    need(gram[0][1]==gram[1][0] and gram[0][0]>0 and gram[1][1]>0 and gram[0][0]*gram[1][1]>gram[0][1]**2,"positive definite Green-Kubo matrix")
    gk=p["green_kubo"]
    exact(sorted(gk),sorted(["convention","matrix","equal_time","correlations","cutoff","partial_sum","exact_remainder","tail_upper"]),"Green-Kubo fields")
    exact(gk["convention"],"sum n=0 to infinity <f,T^n g>_pi; no physical conductivity factor","Green-Kubo convention")
    need(matrix(gk["matrix"],2,2)==gram,"Green-Kubo value")
    equal=[[pair(f,g) for g in currents] for f in currents]
    need(matrix(gk["equal_time"],2,2)==equal,"equal-time covariance")
    need(type(gk["correlations"]) is list and len(gk["correlations"])==65,"correlation count")
    exact(gk["cutoff"],32,"cutoff")
    transported=[list(f) for f in currents];partial=[[F(0)]*2 for _ in range(2)]
    for n,row in enumerate(gk["correlations"]):
        corr=[[pair(f,g) for g in transported] for f in currents]
        need(matrix(row,2,2)==corr,"time-lag correlation")
        if n<=32:
            partial=[[partial[i][j]+corr[i][j] for j in range(2)] for i in range(2)]
        transported=[act(T,f) for f in transported]
    need(rational(gk["correlations"][1][0][0])>0,"nonzero positive-lag memory")
    need(matrix(gk["partial_sum"],2,2)==partial,"finite Green-Kubo sum")
    tail=matrix(gk["exact_remainder"],2,2)
    need(tail==[[gram[i][j]-partial[i][j] for j in range(2)] for i in range(2)],"exact remainder subtraction")
    propagated=[list(g) for g in solutions]
    for _ in range(33):propagated=[act(T,g) for g in propagated]
    need(tail==[[pair(f,g) for g in propagated] for f in currents],"propagated Poisson remainder")
    upper=matrix(gk["tail_upper"],2,2)
    need(upper==[[(equal[i][i]+equal[j][j])*(1-epsilon)**16/epsilon for j in range(2)] for i in range(2)],"geometric infinite tail bound")
    need(all(abs(tail[i][j])<=upper[i][j] for i in range(2) for j in range(2)),"tail containment")
    counter=p["global_centering_counterexample"]
    exact(sorted(counter),sorted(["current","global_mean","persistent_correlation"]),"centering counterexample")
    protected=[F(s[0])-F(2,3) for s in states]
    need(vector(counter["current"],8)==protected and rational(counter["global_mean"])==0,"globally centred protected current")
    need(act(T,protected)==protected and act(P,protected)==protected,"persistent protected mode")
    need(rational(counter["persistent_correlation"])==pair(protected,protected)==F(2,9),"nonzero persistent correlation")
    need(type(p["distribution_history"]) is list and len(p["distribution_history"])==8,"history count")
    initial=[pi[i]*(1+currents[0][i]/2) for i in range(8)];distribution=initial
    index=0;tv=[]
    for n in range(65):
        if n in (0,1,2,4,8,16,32,64):
            row=p["distribution_history"][index];index+=1
            exact(sorted(row),sorted(["step","distribution","total_variation_to_own_fibre_equilibrium"]),"history fields")
            exact(row["step"],n,"history step")
            need(vector(row["distribution"],8)==distribution,"source-law distribution evolution")
            need(sum(distribution)==1 and min(distribution)>0,"faithful distribution")
            need([sum(distribution[:4]),sum(distribution[4:])]==[F(1,3),F(2,3)],"protected record mass")
            error=sum((abs(a-b) for a,b in zip(distribution,pi)),F(0))/2
            need(rational(row["total_variation_to_own_fibre_equilibrium"])==error,"equilibrium distance")
            tv.append(error)
        distribution=[sum((distribution[i]*T[i][j] for i in range(8)),F(0)) for j in range(8)]
    need(all(a>b>0 for a,b in zip(tv,tv[1:])) and tv[-1]<F(1,100000),"nontrivial convergent distribution history")
    return {"accepted":True,"scope":"EXACT_NEW_LOCAL_STOCHASTIC_LAW__PROTECTED_RECORD_AND_DECAYING_MEMORY",
            "states":8,"protected_fibres":2,"current_count":2,"correlation_lags":65,"distribution_samples":8,
            "two_step_minorization":str(epsilon),"gap_lower":str(epsilon/2),
            "green_kubo_matrix":[[str(x) for x in row] for row in gram],
            "native_source_attachment":False,"historical_obstruction_overturned":False,
            "physical_clock_or_conductivity":False,"new_Lean_formalization":False}


if __name__=="__main__":print(json.dumps(verify(load()),sort_keys=True,indent=2))

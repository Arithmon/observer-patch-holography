"""Exact protected-record dynamics with dissipative memory.

This is a NEW supplied stochastic law, not the historical source-counted
transition table.  The full visible-fibre conditional expectation P is its
equilibrium projection; the lazy coordinate heat bath T is its dynamics.

For a faithful law on finite product fibres, coordinate conditional
expectations are orthogonal projections in L2(pi). Their positive lazy
mixture is reversible, positive and fixes the protected fibre. A sequence
updating every coordinate has positive probability to reach every state of
that fibre. Thus T**m >= epsilon P for some epsilon>0. For the two-coordinate
example, m=2 and every entry is checked exactly. Q=(T**2-epsilon P)/(1-epsilon)
is another stationary stochastic kernel, hence an L2(pi) contraction by
Jensen. On ker(P), ||T**n|| <= (1-epsilon)**floor(n/2); the spectral gap is at
least epsilon/2. The limit of an arbitrary distribution is its P image,
not necessarily pi, because the visible record remains conserved.

R=(I-T+P)**(-1)-P is the unique fibre-centred Poisson solver. Detailed balance
makes <f,Rg> symmetric, and positivity follows from the Dirichlet form.
It equals sum(n>=0)<f,T**n g> for fibre-centred currents. The exact remainder
after n=N is <f,T**(N+1)Rg>. The geometric tail is bounded without square
roots by (||f||²+||g||²)*(1-epsilon)**floor((N+1)/2)/epsilon.

The additional source law supplies bounded patches (c,u,v), visible record c,
coordinate ports/readback, conditional repair probabilities and a discrete
update counter. No physical clock, native attachment or continuum is inferred.
Run this module to regenerate only its adjacent receipt.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "runtime" / "protected_memory_receipt.json"
SOURCES = ("protected_memory.py", "verify_protected_memory.py", "test_protected_memory.py")
PARENTS = ("Lean/Thermodynamics/GreenKubo.lean", "Lean/Thermodynamics/FiniteConditionalRepair.lean",
           "code/thermodynamics/common_reference_obstruction/verify_common_reference_obstruction.py",
           "paper/tex_fragments/PROTECTED_RECORD_MEMORY.tex")
STATES = [(c,u,v) for c in range(2) for u in range(2) for v in range(2)]


def identity(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]


def multiply(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def apply(a,v):
    return [sum((x*y for x,y in zip(row,v)),F(0)) for row in a]


def inverse(a):
    n=len(a);rows=[list(row)+eye for row,eye in zip(a,identity(n))]
    for j in range(n):
        pivot=next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[pivot]=rows[pivot],rows[j]
        scale=rows[j][j];rows[j]=[x/scale for x in rows[j]]
        for i in range(n):
            if i!=j:
                scale=rows[i][j]
                rows[i]=[x-scale*y for x,y in zip(rows[i],rows[j])]
    return [row[n:] for row in rows]


def encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,(tuple,list)):return [encode(x) for x in value]
    if isinstance(value,dict):return {k:encode(x) for k,x in value.items()}
    return value


def pin(path):
    data=path.read_bytes()
    return {"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}


def produce():
    weights=((1,2,3,4),(4,1,2,3));masses=(F(1,3),F(2,3))
    pi=[masses[c]*F(weights[c][2*u+v],10) for c,u,v in STATES]
    kernels=[]
    for coordinate in (1,2):
        keep=[k for k in range(3) if k!=coordinate]
        rows=[]
        for x in STATES:
            targets=[j for j,y in enumerate(STATES) if all(x[k]==y[k] for k in keep)]
            mass=sum((pi[j] for j in targets),F(0))
            rows.append([pi[j]/mass if j in targets else F(0) for j in range(8)])
        kernels.append(rows)
    U,V=kernels;I=identity(8)
    P=[[pi[j]/masses[x[0]] if x[0]==y[0] else F(0)
        for j,y in enumerate(STATES)] for x in STATES]
    T=[[I[i][j]/2+(U[i][j]+V[i][j])/4 for j in range(8)] for i in range(8)]
    T2=multiply(T,T)
    epsilon=min(T2[i][j]/P[i][j] for i in range(8) for j in range(8) if P[i][j])
    minor=[[T2[i][j]-epsilon*P[i][j] for j in range(8)] for i in range(8)]
    M=[[I[i][j]-T[i][j]+P[i][j] for j in range(8)] for i in range(8)]
    inv=inverse(M);R=[[inv[i][j]-P[i][j] for j in range(8)] for i in range(8)]
    raw=[[F(x[k]) for x in STATES] for k in (1,2)]
    currents=[[x-y for x,y in zip(f,apply(P,f))] for f in raw]
    solutions=[apply(R,f) for f in currents]
    pair=lambda f,g:sum((w*x*y for w,x,y in zip(pi,f,g)),F(0))
    norms=[pair(f,f) for f in currents]
    gram=[[pair(f,g) for g in solutions] for f in currents]
    covariance=[[pair(f,g) for g in currents] for f in currents]
    C=[];power=I
    for n in range(65):
        C.append([[pair(f,apply(power,g)) for g in currents] for f in currents])
        power=multiply(power,T)
    cutoff=32;partial=[[sum((C[n][i][j] for n in range(cutoff+1)),F(0))
                        for j in range(2)] for i in range(2)]
    remainder=[[gram[i][j]-partial[i][j] for j in range(2)] for i in range(2)]
    tails=[[(norms[i]+norms[j])*(1-epsilon)**((cutoff+1)//2)/epsilon
            for j in range(2)] for i in range(2)]
    p=[pi[i]*(1+currents[0][i]/2) for i in range(8)]
    histories=[]
    for n in range(65):
        if n in (0,1,2,4,8,16,32,64):
            histories.append({"step":n,"distribution":p,
                              "total_variation_to_own_fibre_equilibrium":sum((abs(x-y) for x,y in zip(p,pi)),F(0))/2})
        p=[sum((p[i]*T[i][j] for i in range(8)),F(0)) for j in range(8)]
    protected=[F(x[0])-masses[1] for x in STATES]
    return encode({
        "schema":"oph.protected_record_memory.v1",
        "scope":{"new_supplied_local_stochastic_law":True,"same_reference_equilibrium_and_dynamics":True,
                 "nonconstant_record_preserved":True,"exact_rational_replay":True,
                 "native_source_attachment":False,"historical_obstruction_overturned":False,
                 "global_ergodicity":False,"physical_clock_or_conductivity":False,
                 "continuum_limit":False,"empirical_evidence":False,"new_Lean_formalization":False},
        "model":{"states":STATES,"protected_coordinate":"c","fibre_masses":masses,
                 "conditional_reference_weights":weights,"reference_denominator":10,
                 "step_unit":"one declared stochastic update; no laboratory calibration",
                 "operations":[{"name":"hold","probability":F(1,2),"reads":[],"writes":[]},
                               {"name":"resample_u","probability":F(1,4),"reads":["c","v"],"writes":["u"]},
                               {"name":"resample_v","probability":F(1,4),"reads":["c","u"],"writes":["v"]}]},
        "reference":pi,"coordinate_kernels":{"u":U,"v":V},"equilibrium_projection":P,
        "transition":T,"stationary_pair_coupling":[[pi[i]*T[i][j] for j in range(8)] for i in range(8)],
        "two_step_minorization":{"epsilon":epsilon,"remainder":minor,"gap_lower":epsilon/2,
                                 "L2_decay":"(1-epsilon)^floor(n/2) on ker(P)"},
        "poisson_resolvent":R,"currents":currents,"poisson_solutions":solutions,
        "green_kubo":{"convention":"sum n=0 to infinity <f,T^n g>_pi; no physical conductivity factor",
                      "matrix":gram,"equal_time":covariance,"correlations":C,
                      "cutoff":cutoff,"partial_sum":partial,"exact_remainder":remainder,"tail_upper":tails},
        "global_centering_counterexample":{"current":protected,"global_mean":pair(protected,[F(1)]*8),
                                           "persistent_correlation":pair(protected,protected)},
        "distribution_history":histories,
        "source_pins":{name:pin(HERE/name) for name in SOURCES},
        "parent_pins":{name:pin(ROOT/name) for name in PARENTS}})


if __name__=="__main__":
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(json.dumps(produce(),sort_keys=True,separators=(",",":"),allow_nan=False)+"\n",encoding="utf-8")
    print(OUTPUT)

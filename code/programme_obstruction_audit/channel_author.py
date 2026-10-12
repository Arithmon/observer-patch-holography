"""Exact small-matrix audit at OPH target 0510578e317de5da0d485d09792dc8bcad2e2f3c.

No repository imports: Fraction arithmetic checks the binary seam channel
described by docs/MODEL_SEARCH.md: odd-block depolarization, unchanged even
block, and halved even/odd coherence. The private factor is one toy qubit;
this is an escape-route witness in a broader channel class, not an A1 model.
"""
from fractions import Fraction as F
import json

def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]

def eye(n):
    a=zero(n)
    for i in range(n): a[i][i]=F(1)
    return a

def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]

def trans(a):
    return [list(x) for x in zip(*a)]

def plus(a,b):
    return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]

def scale(r,a):
    return [[r*x for x in ar] for ar in a]

def reduced(x):
    return [[x[2*i][2*j]+x[2*i+1][2*j+1] for j in range(4)]
            for i in range(4)]

def C(x):
    # Load basis 00,10,01,11; private basis 0,1.
    y=zero(8)
    for i in range(4):
        for j in range(4):
            for p in range(2):
                for q in range(2):
                    if i in (1,2) and j in (1,2):
                        y[2*i+p][2*j+q]=((x[2+p][2+q]+x[4+p][4+q])/2
                                           if i==j else F(0))
                    else:
                        factor=F(1,2) if ((i in (1,2)) != (j in (1,2))) else F(1)
                        y[2*i+p][2*j+q]=factor*x[2*i+p][2*j+q]
    return y

def E(x):
    y=C(x)
    for i in range(4):
        for j in range(4):
            if (i in (1,2)) != (j in (1,2)):
                for p in range(2):
                    for q in range(2): y[2*i+p][2*j+q]=F(0)
    return y

def native_kraus():
    # L0=I_even+I_odd/2; remaining Lr are odd X, real iY, Z over 2.
    load=[]
    k=zero(4);k[0][0]=k[3][3]=F(1);k[1][1]=k[2][2]=F(1,2);load.append(k)
    k=zero(4);k[1][2]=k[2][1]=F(1,2);load.append(k)
    k=zero(4);k[1][2]=F(1,2);k[2][1]=F(-1,2);load.append(k)
    k=zero(4);k[1][1]=F(1,2);k[2][2]=F(-1,2);load.append(k)
    out=[]
    for k in load:
        lifted=zero(8)
        for i in range(4):
            for j in range(4):
                for p in range(2): lifted[2*i+p][2*j+p]=k[i][j]
        out.append(lifted)
    return out

ks=native_kraus()
complete=zero(8)
for k in ks: complete=plus(complete,mm(trans(k),k))
assert complete==eye(8)
for i in range(8):
    for j in range(8):
        e=zero(8);e[i][j]=F(1)
        got=zero(8)
        for k in ks:got=plus(got,mm(mm(k,e),trans(k)))
        assert got==C(e)

def CW(x,w):
    return mm(mm(w,C(mm(mm(trans(w),x),w))),trans(w))

def rotation(w,i,j,c,s):
    w[i][i]=w[j][j]=c
    w[i][j]=s
    w[j][i]=-s

def check_every_matrix_unit(w):
    assert mm(w,trans(w))==eye(8)
    for i in range(8):
        for j in range(8):
            e=zero(8);e[i][j]=F(1)
            assert reduced(CW(e,w))==reduced(C(e))
    return 64

# Independent rational private-controlled load rotations. Such a W is one
# product exponential after diagonalizing the private generator. Different
# rational rotations exercise correlated inputs through all matrix units.
w_controlled=eye(8)
rotation(w_controlled,2,4,F(3,5),F(4,5))
rotation(w_controlled,3,5,F(5,13),F(12,13))
checks=check_every_matrix_unit(w_controlled)

# Exchange W only mixes |10>_load|1>_private and |01>_load|0>_private.
# Its generator is s_+ tensor |1><0| - s_- tensor |0><1|;
# the private amplitude is nonnormal. This is a two-quadrature exchange,
# outside the single Hermitian-load tensor skew-Hermitian-private class.
w_exchange=eye(8)
rotation(w_exchange,3,4,F(3,5),F(4,5))
assert mm(w_exchange,trans(w_exchange))==eye(8)
probs=[]
for p in range(2):
    x=zero(8)
    x[2+p][2+p]=x[4+p][4+p]=F(1,2)
    y=CW(x,w_exchange)
    assert sum(y[i][i] for i in range(8))==1
    probs.append(sum(y[2+p2][2+p2] for p2 in range(2)))
assert probs[0]!=probs[1]

# All matrix units also validate inherited single-seam completion at n=2,3.
for i in range(8):
    for j in range(8):
        e=zero(8);e[i][j]=F(1)
        ew=mm(mm(w_exchange,E(mm(mm(trans(w_exchange),e),w_exchange))),trans(w_exchange))
        c1=CW(e,w_exchange)
        cn=c1
        for n in (2,3):
            cn=CW(cn,w_exchange)
            want=plus(ew,scale(F(1,2**(n-1)),plus(c1,scale(-1,ew))))
            assert cn==want

packet={
    "target_revision":"0510578e317de5da0d485d09792dc8bcad2e2f3c",
    "channel":"binary seam C, load dimension 4, toy private dimension 2",
    "arithmetic":"exact rational, no floats, no repository implementation imports",
    "controlled_rotation_matrix_units_checked":checks,
    "CPTP_Kraus_completeness_and_channel_matrix_units_checked":64,
    "exchange_completion_matrix_units_checked_each_at_n2_n3":64,
    "exchange_rotation_cos":"3/5",
    "exchange_rotation_sin":"4/5",
    "fixed_load_preparation":"I_odd/2",
    "private_preparations":["|0><0|","|1><1|"],
    "P_original_endpoint_occupied":[str(p) for p in probs],
    "difference":str(probs[0]-probs[1]),
    "scope":"independent toy escape witness; neither complete A1-A3 nor OPH exclusion",
    "review_status":"author exact check; independent review pending",
}
print(json.dumps(packet,indent=2))

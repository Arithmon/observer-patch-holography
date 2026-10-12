"""Independent Fraction replay in private-first tensor order; no author imports."""
from fractions import Fraction as F
from hashlib import sha256
import json
import platform
from pathlib import Path

def zero(n):return [[F(0) for _ in range(n)] for _ in range(n)]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def transpose(a):return list(map(list,zip(*a)))
def multiply(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def times(a,t):return [[t*x for x in r] for r in a]
def conjugate(a,x):return multiply(multiply(a,x),transpose(a))
def kron(a,b):return [[a[i//len(b)][j//len(b)]*b[i%len(b)][j%len(b)] for j in range(len(a)*len(b))] for i in range(len(a)*len(b))]
def projector(indices):
    a=zero(4)
    for i in indices:a[i][i]=F(1)
    return a
def rotation(i,j,c,s):
    a=eye(8);a[i][i]=a[j][j]=c;a[i][j]=s;a[j][i]=-s
    return a
def partial(x):return [[sum(x[4*p+i][4*p+j] for p in range(2)) for j in range(4)] for i in range(4)]
def channel(ks,x):
    a=zero(8)
    for w,k in ks:a=add(a,times(conjugate(k,x),w))
    return a
def rank(a):
    a=[r[:] for r in a];rr=0
    for c in range(len(a[0])):
        p=next((i for i in range(rr,len(a)) if a[i][c]),None)
        if p is None:continue
        a[p],a[rr]=a[rr],a[p];z=a[rr][c];a[rr]=[x/z for x in a[rr]]
        for i in range(len(a)):
            if i!=rr and a[i][c]:
                z=a[i][c];a[i]=[x-z*y for x,y in zip(a[i],a[rr])]
        rr+=1
    return rr

Pe=projector([0,3]);Po=projector([1,2])
ks=[]
k=add(Pe,times(Po,F(1,2)));ks.append((F(1),kron(eye(2),k)))
for typ in ['X','iY','Z']:
    k=zero(4)
    if typ=='X':k[1][2]=k[2][1]=F(1,2)
    if typ=='iY':k[1][2]=F(1,2);k[2][1]=F(-1,2)
    if typ=='Z':k[1][1]=F(1,2);k[2][2]=F(-1,2)
    ks.append((F(1),kron(eye(2),k)))
es=[(F(1),kron(eye(2),Pe))]
for i in [1,2]:
    for j in [1,2]:
        k=zero(4);k[i][j]=F(1);es.append((F(1,2),kron(eye(2),k)))
complete=zero(8)
for weight,k in ks:complete=add(complete,times(multiply(transpose(k),k),weight))
assert complete==eye(8)
completeE=zero(8)
for weight,k in es:completeE=add(completeE,times(multiply(transpose(k),k),weight))
assert completeE==eye(8)

W=rotation(5,2,F(3,5),F(4,5))
assert multiply(W,transpose(W))==eye(8)
N=kron(eye(2),[[F(i==j)*F([0,1,1,2][i]) for j in range(4)] for i in range(4)])
assert multiply(W,N)==multiply(N,W)
kw=[(weight,multiply(multiply(W,k),transpose(W))) for weight,k in ks]
ew=[(weight,multiply(multiply(W,k),transpose(W))) for weight,k in es]
completeW=zero(8)
for weight,k in kw:completeW=add(completeW,times(multiply(transpose(k),k),weight))
assert completeW==eye(8)
U=eye(8)
for i,j,c,s in [(1,2,F(3,5),F(4,5)),(5,6,F(5,13),F(12,13))]:
    U=multiply(rotation(i,j,c,s),U)
ku=[(weight,multiply(multiply(U,k),transpose(U))) for weight,k in ks]
reduced_failures=0
for i in range(8):
    for j in range(8):
        x=zero(8);x[i][j]=F(1)
        # Operator block definition independently of author C code.
        block=zero(8)
        for p in range(2):
            for q in range(2):
                for a in range(4):
                    for b in range(4):
                        if a in [1,2] and b in [1,2]:
                            block[4*p+a][4*q+b]=(x[4*p+1][4*q+1]+x[4*p+2][4*q+2])/2 if a==b else F(0)
                        else:block[4*p+a][4*q+b]=x[4*p+a][4*q+b]*(F(1,2) if ((a in [1,2])!=(b in [1,2])) else F(1))
        assert channel(ks,x)==block
        assert partial(channel(ku,x))==partial(channel(ks,x))
        e=channel(ew,x);c=channel(kw,x)
        assert channel(ew,e)==e
        assert channel(ew,c)==e and channel(kw,e)==e
        for n in [2,3]:
            c=channel(kw,c)
            assert c==add(e,times(add(channel(kw,x),times(e,F(-1))),F(1,2**(n-1))))
        if partial(channel(kw,x))!=partial(channel(ks,x)):reduced_failures+=1

# Private commutant: [I_load tensor G, W]=0 in private-first layout.
columns=[]
for i in range(2):
    for j in range(2):
        g=zero(2);g[i][j]=F(1);g=kron(g,eye(4))
        comm=add(multiply(g,W),times(multiply(W,g),F(-1)))
        columns.append([x for row in comm for x in row])
commutant_rank=rank(list(map(list,zip(*columns))))
assert commutant_rank==3 # 4 parameters, identity-only kernel over R or C.
probs=[]
for p in [0,1]:
    state=zero(8);state[4*p+1][4*p+1]=state[4*p+2][4*p+2]=F(1,2)
    out=channel(kw,state)
    probs.append(sum(out[4*q+1][4*q+1] for q in [0,1]))
assert probs==[F(769,1250),F(481,1250)]
assert probs[0]-probs[1]==F(144,625)
assert probs[0]-probs[1]==F(3,5)**2*F(4,5)**2

packet={'target_revision':'0510578e317de5da0d485d09792dc8bcad2e2f3c',
    'python':platform.python_version(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
    'tensor_order':'private-first; independent of author load-first implementation',
    'C_CPTP':True,'E_CPTP':True,'conjugated_C_CPTP':True,'W_number_conserving':True,
    'native_block_channel_matrix_units':64,'controlled_reduced_channel_matrix_units':64,
    'conjugated_completion_matrix_units_each_n2_n3':64,
    'native_reduced_channel_failures_for_exchange_matrix_units':reduced_failures,
    'private_commutant_constraint_rank':commutant_rank,'private_commutant_dimension':1,
    'feedback_probabilities':list(map(str,probs)),'difference':str(probs[0]-probs[1]),
    'symbolic_difference':'cos(theta)^2 sin(theta)^2',
    'scope':'finite toy channel; all operator matrix units, no full A1-A3 or physical force'}
print(json.dumps(packet,indent=2))

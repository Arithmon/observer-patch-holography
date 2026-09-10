"""Exact local gauge/Grassmann jet compiler on a declared flat Spin chart.

All connection and scalar test jets are exact Gaussian rationals. Fermion values
and first derivatives remain independent anticommuting generators. Physical
couplings, Higgs mass/quartic and Yukawa entries are symbolic sector coefficients.
"""
from fractions import Fraction as F
from dataclasses import dataclass
import hashlib,json

@dataclass(frozen=True)
class C:
 r:F=F(0)
 i:F=F(0)
 def __post_init__(self):
  if any(type(x) not in (int,F) for x in (self.r,self.i)):raise TypeError('exact integer/Fraction coefficients required')
  object.__setattr__(self,'r',F(self.r));object.__setattr__(self,'i',F(self.i))
 def __add__(self,x):
  if not isinstance(x,(C,F,int)):return NotImplemented
  x=c(x);return C(self.r+x.r,self.i+x.i)
 __radd__=__add__
 def __neg__(self):return C(-self.r,-self.i)
 def __sub__(self,x):
  if not isinstance(x,(C,F,int)):return NotImplemented
  return self+-c(x)
 def __rsub__(self,x):return c(x)+-self
 def __mul__(self,x):
  if not isinstance(x,(C,F,int)):return NotImplemented
  x=c(x);return C(self.r*x.r-self.i*x.i,self.r*x.i+self.i*x.r)
 __rmul__=__mul__
 def __truediv__(self,x):
  x=c(x);d=x.r*x.r+x.i*x.i
  if not d:raise ZeroDivisionError('Gaussian rational division')
  return self*x.conj()*C(1/d)
 def conj(self):return C(self.r,-self.i)
 def __bool__(self):return bool(self.r or self.i)
 def encode(self):return [str(self.r),str(self.i)]
def c(x):return x if isinstance(x,C) else C(x)
I=C(0,1)

class P:
 """Sparse exterior algebra; monomials are strictly increasing generator IDs."""
 def __init__(self,terms=None):self.t={k:c(v) for k,v in (terms or {}).items() if c(v)}
 @staticmethod
 def scalar(x):return P({():c(x)})
 @staticmethod
 def generator(n):
  if type(n) is not int or n<0:raise ValueError('nonnegative integer generator')
  return P({(n,):C(1)})
 def __add__(self,x):
  if isinstance(x,D):return NotImplemented
  x=p(x);out=dict(self.t)
  for k,v in x.t.items():
   out[k]=out.get(k,C())+v
   if not out[k]:del out[k]
  return P(out)
 __radd__=__add__
 def __neg__(self):return P({k:-v for k,v in self.t.items()})
 def __sub__(self,x):return self+-p(x)
 def __rsub__(self,x):return p(x)+-self
 def __mul__(self,x):
  if isinstance(x,D):return NotImplemented
  x=p(x);out={}
  for a,u in self.t.items():
   for b,v in x.t.items():
    if set(a).intersection(b):continue
    sign=-1 if sum(i>j for i in a for j in b)%2 else 1;k=tuple(sorted(a+b))
    out[k]=out.get(k,C())+sign*u*v
    if not out[k]:del out[k]
  return P(out)
 __rmul__=__mul__
 def adjoint(self,partners):
  out=P()
  for word,value in self.t.items():
   term=P.scalar(value.conj())
   for i in reversed(word):term=term*P.generator(partners[i])
   out=out+term
  return out
 def encode(self):return [[list(k),v.encode()] for k,v in sorted(self.t.items())]
 def __eq__(self,x):return self.t==p(x).t
 def __bool__(self):return bool(self.t)
def p(x):return x if isinstance(x,P) else P.scalar(x)

@dataclass
class D:
 """A first variation; the infinitesimal parameter is commuting and squares to0."""
 value:P
 delta:P
 def __init__(self,value=0,delta=0):self.value=p(value);self.delta=p(delta)
 def __add__(self,x):
  x=dual(x);return D(self.value+x.value,self.delta+x.delta)
 __radd__=__add__
 def __neg__(self):return D(-self.value,-self.delta)
 def __sub__(self,x):return self+-dual(x)
 def __rsub__(self,x):return dual(x)+-self
 def __mul__(self,x):
  x=dual(x);return D(self.value*x.value,self.delta*x.value+self.value*x.delta)
 __rmul__=__mul__
 def adjoint(self,partners):return D(self.value.adjoint(partners),self.delta.adjoint(partners))
def dual(x):return x if isinstance(x,D) else D(x)
def zero(n):return [[C() for _ in range(n)] for _ in range(n)]
def madd(A,B):return [[a+b for a,b in zip(x,y)] for x,y in zip(A,B)]
def mscale(a,A):return [[a*x for x in row] for row in A]
def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),C()) for j in range(len(B[0]))] for i in range(len(A))]
def comm(A,B):return madd(mm(A,B),mscale(-1,mm(B,A)))
def transpose(A):return list(map(list,zip(*A)))
def tr(A):return sum((A[i][i] for i in range(len(A))),0)
def basis(n):
 out=[]
 for a in range(n):
  for b in range(a+1,n):
   x=zero(n);x[a][b]=C(1);x[b][a]=C(-1);out.append(x)
   x=zero(n);x[a][b]=x[b][a]=I;out.append(x)
 for a in range(n-1):
  x=zero(n);x[a][a]=I;x[a+1][a+1]=-I;out.append(x)
 return out
BASES=[*[(0,t) for t in basis(3)],*[(1,t) for t in basis(2)],(2,[[I]])]
FIELDS=[('Q',3,2,F(1,6)),('u_c',-3,1,F(-2,3)),('d_c',-3,1,F(1,3)),('L',1,2,F(-1,2)),('e_c',1,1,F(1))]

class Algebra:
 def __init__(self,families=3):
  if type(families) is not int or not 1<=families<=3:raise ValueError('families must be1..3')
  self.families=families;self.labels=[];self.ids={};self.partners={}
  for family in range(families):
   for name,color,weak,_ in FIELDS:
    for k in range(abs(color)*weak):
     for spin in range(2):
      for jet in ('value','d0','d1','d2','d3'):
       i=len(self.labels)
       for bar in (False,True):
        label=(name,family,k,spin,jet,bar);self.ids[label]=len(self.labels);self.labels.append(label)
       self.partners[i]=i+1;self.partners[i+1]=i
 def variable(self,name,fam,k,s,jet):return P.generator(self.ids[name,fam,k,s,jet,False])

def representation(connections,color,weak,Y):
 """Tensor product of the registered fundamental/dual/singlet carriers."""
 col=zero(1) if abs(color)==1 else connections[0]
 if color<0:col=mscale(-1,transpose(col))
 w=connections[1] if weak==2 else zero(1);n=abs(color)*weak;A=zero(n)
 for ca in range(abs(color)):
  for a in range(weak):
   for cb in range(abs(color)):
    for b in range(weak):A[ca*weak+a][cb*weak+b]=(col[ca][cb] if a==b else C())+(w[a][b] if ca==cb else C())+(Y*connections[2][0][0] if ca==cb and a==b else C())
 return A

def fixture(seed=0):
 """Small exact test jets; these are not coupling values or measured fields."""
 if type(seed) is not int or not 0<=seed<=3:raise ValueError('fixture seed0..3')
 def triple(values):
  out=[zero(3),zero(2),zero(1)]
  for value,(sector,T) in zip(values,BASES):out[sector]=madd(out[sector],mscale(value,T))
  return out
 A=[triple([F(((13*mu+7*a+3*seed)%11)-5,7) for a in range(12)]) for mu in range(4)]
 dA=[[triple([F(((17*mu+19*nu+5*a+seed)%13)-6,11) for a in range(12)]) for nu in range(4)] for mu in range(4)]
 H=[C(F(2+seed,3),F(1,5)),C(F(-1,4),F(3+seed,7))]
 dH=[[C(F(mu+a+1,5),F(2*mu-a+seed,9)) for a in range(2)] for mu in range(4)]
 return {'A':A,'dA':dA,'H':H,'dH':dH}

def direction(kind='theta',generator=0,mu=0,nu=0,variant='correct'):
 if kind not in ('theta','dtheta','d2theta','connection'):raise ValueError('direction kind')
 if type(generator) is not int or not 0<=generator<12 or any(type(k) is not int or not 0<=k<4 for k in (mu,nu)):raise ValueError('direction index')
 z=lambda:[zero(3),zero(2),zero(1)]
 theta=z();dt=[z() for _ in range(4)];d2=[[z() for _ in range(4)] for _ in range(4)];external=[z() for _ in range(4)]
 sec,T=BASES[generator]
 if kind=='theta':theta[sec]=T
 elif kind=='dtheta':dt[mu][sec]=T
 elif kind=='d2theta':d2[mu][nu][sec]=T;d2[nu][mu][sec]=T
 else:external[mu][sec]=T
 return theta,dt,d2,external

def compile_action(algebra,jets,gauge_direction=None,variant='correct'):
 """Return coefficient-indexed local density and its exact first variation."""
 if variant not in ('correct','drop_dtheta_connection','drop_A_commutator','wrong_dual_color','mixed_up_hypercharge'):raise ValueError('unknown compiler variant')
 z=lambda:[zero(3),zero(2),zero(1)]
 theta,dt,d2,external=gauge_direction or (z(),[z() for _ in range(4)],[[z() for _ in range(4)] for _ in range(4)],[z() for _ in range(4)])
 A=jets['A'];dA=jets['dA'];deltA=[];deltdA=[]
 for mu in range(4):
  deltA.append([madd(madd(mscale(0 if variant=='drop_dtheta_connection' else -1,dt[mu][s]),comm(theta[s],A[mu][s])),external[mu][s]) for s in range(3)])
  deltdA.append([[madd(mscale(-1,d2[mu][nu][s]),madd(comm(dt[mu][s],A[nu][s]),comm(theta[s],dA[mu][nu][s]))) for s in range(3)] for nu in range(4)])
 conns=[[[[D(a,d) for a,d in zip(ar,dr)] for ar,dr in zip(aa,dd)] for aa,dd in zip(A[mu],deltA[mu])] for mu in range(4)]
 first=[[[[[D(a,d) for a,d in zip(ar,dr)] for ar,dr in zip(aa,dd)] for aa,dd in zip(dA[mu][nu],deltdA[mu][nu])] for nu in range(4)] for mu in range(4)]
 eta=[1,-1,-1,-1];out={}
 for s,name in enumerate(('SU3','SU2','U1')):
  density=D()
  for mu in range(4):
   for nu in range(mu+1,4):
    curvature=madd(madd(first[mu][nu][s],mscale(-1,first[nu][mu][s])),mscale(0 if variant=='drop_A_commutator' else 1,comm(conns[mu][s],conns[nu][s])))
    density=density+eta[mu]*eta[nu]*tr(mm(curvature,curvature))
  out['kappa_'+name]=density
 Th=representation(theta,1,2,F(1,2));dTh=[representation(t,1,2,F(1,2)) for t in dt];H=jets['H'];dH=jets['dH']
 h=[D(H[a],sum((Th[a][b]*H[b] for b in range(2)),C())) for a in range(2)]
 dh=[[D(dH[mu][a],sum((dTh[mu][a][b]*H[b]+Th[a][b]*dH[mu][b] for b in range(2)),C())) for a in range(2)] for mu in range(4)]
 scalar=D()
 for mu in range(4):
  conn=representation(conns[mu],1,2,F(1,2));cov=[dh[mu][a]+sum((conn[a][b]*h[b] for b in range(2)),D()) for a in range(2)]
  scalar=scalar+eta[mu]*sum((v.adjoint(algebra.partners)*v for v in cov),D())
 norm=sum((v.adjoint(algebra.partners)*v for v in h),D());out['scalar_kinetic']=scalar;out['m_squared']=-norm;out['lambda']=-norm*norm
 sigma=[[[C(1),C()],[C(),C(1)]],[[C(),C(-1)],[C(-1),C()]],[[C(),I],[-I,C()]],[[C(-1),C()],[C(),C(1)]]]
 psi={}
 for name,color,weak,Y in FIELDS:
  if variant=='wrong_dual_color' and color<0:color=-color
  T=representation(theta,color,weak,Y);Td=[representation(t,color,weak,Y) for t in dt];conn=[representation(t,color,weak,Y) for t in conns];dim=abs(color)*weak;kin=D()
  for family in range(algebra.families):
   raw=[[algebra.variable(name,family,k,a,'value') for a in range(2)] for k in range(dim)]
   spinor=[[D(raw[k][a],sum((T[k][l]*raw[l][a] for l in range(dim)),P())) for a in range(2)] for k in range(dim)]
   psi[name,family]=spinor
   for mu in range(4):
    draw=[[algebra.variable(name,family,k,a,'d'+str(mu)) for a in range(2)] for k in range(dim)]
    dspinor=[[D(draw[k][a],sum((Td[mu][k][l]*raw[l][a]+T[k][l]*draw[l][a] for l in range(dim)),P())) for a in range(2)] for k in range(dim)]
    cov=[[dspinor[k][a]+sum((conn[mu][k][l]*spinor[l][a] for l in range(dim)),D()) for a in range(2)] for k in range(dim)]
    form=sum((sigma[mu][a][b]*spinor[k][a].adjoint(algebra.partners)*cov[k][b] for k in range(dim) for a in range(2) for b in range(2)),D())
    kin=kin+D(I/2)*(form-form.adjoint(algebra.partners))
  out['kinetic_'+name]=kin
 eps=[[0,1],[-1,0]]
 for family in range(algebra.families):
  for other in range(algebra.families):
   channels={'up':D(),'down':D(),'lepton':D()}
   for a in range(2):
    for b in range(2):
     if not eps[a][b]:continue
     for color in range(3):
      for weak in range(2):
       Q=psi['Q',family][2*color+weak][a]
       for hw in range(2):
        hs=h[hw].adjoint(algebra.partners) if variant=='mixed_up_hypercharge' else h[hw]
        channels['up']=channels['up']+eps[a][b]*eps[weak][hw]*Q*hs*psi['u_c',other][color][b]
       channels['down']=channels['down']+eps[a][b]*Q*h[weak].adjoint(algebra.partners)*psi['d_c',other][color][b]
     for weak in range(2):channels['lepton']=channels['lepton']+eps[a][b]*psi['L',family][weak][a]*h[weak].adjoint(algebra.partners)*psi['e_c',other][0][b]
   for channel,T in channels.items():
    out[f'Y_{channel}_{family}_{other}_re']=-(T+T.adjoint(algebra.partners));out[f'Y_{channel}_{family}_{other}_im']=-D(I)*(T-T.adjoint(algebra.partners))
 return out

def encode_action(action,part='value'):
 return {k:getattr(v,part).encode() for k,v in sorted(action.items())}
def canonical(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()

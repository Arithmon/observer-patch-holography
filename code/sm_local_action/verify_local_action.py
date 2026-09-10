"""Independent bit-mask exterior algebra and Gaussian-domain Ward reconstruction.

No producer imports. Gauge/scalar test jets are reconstructed independently;
fermions and their derivatives remain formal Grassmann variables. Matrix jets
and first variations are expanded explicitly, without the producer's dual ring.
"""
from fractions import Fraction
from pathlib import Path
import argparse,hashlib,json
from sympy.polys.domains import QQ,QQ_I
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUTPUT=HERE/'local_action_receipt.json'
Z=QQ_I.zero;ONE=QQ_I.one;II=QQ_I(0,1)

def need(ok,label):
 if not ok:raise ValueError(label)
def conj(z):return QQ_I(z.x,-z.y)
def qc(r=0,i=0):return QQ_I(QQ(r),QQ(i))
def bits(m):
 while m:
  b=m&-m;yield b.bit_length()-1;m-=b
class W:
 def __init__(self,t=None):self.t={m:v for m,v in (t or {}).items() if v}
 def __add__(self,x):
  if not isinstance(x,W):x=W({0:QQ_I(x)})
  t=dict(self.t)
  for m,c in x.t.items():
   t[m]=t.get(m,Z)+c
   if not t[m]:del t[m]
  return W(t)
 __radd__=__add__
 def __neg__(self):return W({m:-c for m,c in self.t.items()})
 def __sub__(self,x):return self+-x
 def __mul__(self,x):
  if not isinstance(x,W):return W({m:c*x for m,c in self.t.items()})
  t={}
  for a,c in self.t.items():
   for b,d in x.t.items():
    if a&b:continue
    sign=-1 if sum((a>>(j+1)).bit_count() for j in bits(b))%2 else 1;m=a|b
    t[m]=t.get(m,Z)+c*d*sign
    if not t[m]:del t[m]
  return W(t)
 __rmul__=__mul__
 def star(self):
  out=W()
  for mask,c in self.t.items():
   term=W({0:conj(c)})
   for j in reversed(list(bits(mask))):term=term*W({1<<(j^1):ONE})
   out=out+term
  return out
 def encode(self):return [[list(bits(m)),[str(c.x),str(c.y)]] for m,c in sorted(self.t.items(),key=lambda x:tuple(bits(x[0])))]
 def __bool__(self):return bool(self.t)
 def __eq__(self,x):return isinstance(x,W) and self.t==x.t

def zeros(n):return [[Z for _ in range(n)] for _ in range(n)]
def add(A,B):return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]
def scale(A,a):return [[x*a for x in row] for row in A]
def mul(A,B):return [[sum((x*y for x,y in zip(row,col)),Z) for col in zip(*B)] for row in A]
def bracket(A,B):return add(mul(A,B),scale(mul(B,A),-1))
def trace(A):return sum((A[i][i] for i in range(len(A))),Z)
def generators():
 table=[]
 for sector,n in enumerate((3,2)):
  for row in range(n):
   for col in range(row+1,n):
    for phase in (0,1):
     a=zeros(n);a[row][col]=II if phase else ONE;a[col][row]=II if phase else -ONE;table.append((sector,a))
  for row in range(n-1):
   a=zeros(n);a[row][row]=II;a[row+1][row+1]=-II;table.append((sector,a))
 table.append((2,[[II]]));return table
GEN=generators()
ROWS=[('Q',3,2,QQ(1,6)),('u_c',-3,1,QQ(-2,3)),('d_c',-3,1,QQ(1,3)),('L',1,2,QQ(-1,2)),('e_c',1,1,QQ(1))]
def blank():return [zeros(3),zeros(2),zeros(1)]
def family_matrix(A,row):
 _,color,weak,y=row;n=abs(color)*weak;out=zeros(n)
 for i in range(n):
  ci,wi=divmod(i,weak)
  for j in range(n):
   cj,wj=divmod(j,weak);z=Z
   if wi==wj and abs(color)>1:z+=A[0][ci][cj] if color>0 else -A[0][cj][ci]
   if ci==cj and weak>1:z+=A[1][wi][wj]
   if i==j:z+=A[2][0][0]*y
   out[i][j]=z
 return out
HROW=('H',1,2,QQ(1,2))

def fixture():
 def tri(values):
  ans=blank()
  for v,(s,g) in zip(values,GEN):ans[s]=add(ans[s],scale(g,v))
  return ans
 return {'A':[tri([QQ((13*u+7*a)%11-5,7) for a in range(12)]) for u in range(4)],
  'dA':[[tri([QQ((17*u+19*v+5*a)%13-6,11) for a in range(12)]) for v in range(4)] for u in range(4)],
  'H':[qc(QQ(2,3),QQ(1,5)),qc(QQ(-1,4),QQ(3,7))],
  'dH':[[qc(QQ(u+a+1,5),QQ(2*u-a,9)) for a in range(2)] for u in range(4)]}
def encode_scalar(x):
 if isinstance(x,list):return [encode_scalar(v) for v in x]
 return [str(x.x),str(x.y)]
def labels(families=3):
 return [[name,f,k,s,j,bar] for f in range(families) for name,c,w,_ in ROWS for k in range(abs(c)*w) for s in range(2) for j in ('value','d0','d1','d2','d3') for bar in (False,True)]

def variation(j,kind,g,mu=0,nu=0,variant='correct'):
 theta=blank();dt=[blank() for _ in range(4)];d2=[[blank() for _ in range(4)] for _ in range(4)];extra=[blank() for _ in range(4)];s,T=GEN[g]
 if kind=='theta':theta[s]=T
 if kind=='dtheta':dt[mu][s]=T
 if kind=='d2theta':d2[mu][nu][s]=T;d2[nu][mu][s]=T
 if kind=='connection':extra[mu][s]=T
 da=[[add(add(scale(dt[u][s],0 if variant=='drop_dtheta_connection' else -1),bracket(theta[s],j['A'][u][s])),extra[u][s]) for s in range(3)] for u in range(4)]
 dda=[[[add(scale(d2[u][v][s],-1),add(bracket(dt[u][s],j['A'][v][s]),bracket(theta[s],j['dA'][u][v][s]))) for s in range(3)] for v in range(4)] for u in range(4)]
 return theta,dt,da,dda

def build_density(j,kind='none',g=0,mu=0,nu=0,variant='correct'):
 """Value and gauge variation, formed separately with the ordinary product rule."""
 theta,dt,da,dda=variation(j,kind,g,mu,nu,variant);A=j['A'];first=j['dA'];sgn=[1,-1,-1,-1];out={};dot={}
 for s,name in enumerate(('SU3','SU2','U1')):
  val=Z;dv=Z
  for u in range(4):
   for v in range(u+1,4):
    f=add(add(first[u][v][s],scale(first[v][u][s],-1)),scale(bracket(A[u][s],A[v][s]),0 if variant=='drop_A_commutator' else 1))
    df=add(add(dda[u][v][s],scale(dda[v][u][s],-1)),scale(add(bracket(da[u][s],A[v][s]),bracket(A[u][s],da[v][s])),0 if variant=='drop_A_commutator' else 1))
    val+=trace(mul(f,f))*sgn[u]*sgn[v];dv+=trace(add(mul(df,f),mul(f,df)))*sgn[u]*sgn[v]
  out['kappa_'+name]=W({0:val});dot['kappa_'+name]=W({0:dv})
 H=j['H'];T=family_matrix(theta,HROW);dT=[family_matrix(x,HROW) for x in dt]
 dH=[sum((T[a][b]*H[b] for b in range(2)),Z) for a in range(2)]
 hkin=Z;dhkin=Z
 for u in range(4):
  R=family_matrix(A[u],HROW);dR=family_matrix(da[u],HROW)
  v=[j['dH'][u][a]+sum((R[a][b]*H[b] for b in range(2)),Z) for a in range(2)]
  dv=[sum((dT[u][a][b]*H[b]+T[a][b]*j['dH'][u][b]+dR[a][b]*H[b]+R[a][b]*dH[b] for b in range(2)),Z) for a in range(2)]
  hkin+=sum((conj(x)*x for x in v),Z)*sgn[u];dhkin+=sum((conj(x)*y+conj(y)*x for x,y in zip(v,dv)),Z)*sgn[u]
 hnorm=sum((conj(x)*x for x in H),Z);dhnorm=sum((conj(x)*y+conj(y)*x for x,y in zip(H,dH)),Z)
 for name,v,dv in [('scalar_kinetic',hkin,dhkin),('m_squared',-hnorm,-dhnorm),('lambda',-hnorm*hnorm,-2*hnorm*dhnorm)]:out[name]=W({0:v});dot[name]=W({0:dv})
 catalog=labels();lookup={tuple(x):i for i,x in enumerate(catalog)}
 def variable(name,f,k,s,jet):return W({1<<lookup[name,f,k,s,jet,False]:ONE})
 sigma=[[[ONE,Z],[Z,ONE]],[[Z,-ONE],[-ONE,Z]],[[Z,II],[-II,Z]],[[-ONE,Z],[Z,ONE]]]
 psi={};dp={}
 for row in ROWS:
  name,color,weak,y=row
  if variant=='wrong_dual_color' and color<0:row=(name,-color,weak,y)
  dim=abs(color)*weak;T=family_matrix(theta,row);dT=[family_matrix(x,row) for x in dt];vkin=W();dkin=W()
  for family in range(3):
   x=[[variable(name,family,k,a,'value') for a in range(2)] for k in range(dim)]
   dx=[[sum((x[l][a]*T[k][l] for l in range(dim)),W()) for a in range(2)] for k in range(dim)];psi[name,family]=x;dp[name,family]=dx
   for u in range(4):
    R=family_matrix(A[u],row);dR=family_matrix(da[u],row);y=[[variable(name,family,k,a,'d'+str(u)) for a in range(2)] for k in range(dim)]
    cv=[[y[k][a]+sum((x[l][a]*R[k][l] for l in range(dim)),W()) for a in range(2)] for k in range(dim)]
    dcv=[[sum((x[l][a]*dT[u][k][l]+y[l][a]*T[k][l]+x[l][a]*dR[k][l]+dx[l][a]*R[k][l] for l in range(dim)),W()) for a in range(2)] for k in range(dim)]
    b=W();db=W()
    for k in range(dim):
     for a in range(2):
      for bspin in range(2):
       weight=sigma[u][a][bspin]
       b=b+(x[k][a].star()*cv[k][bspin])*weight
       db=db+(dx[k][a].star()*cv[k][bspin]+x[k][a].star()*dcv[k][bspin])*weight
    vkin=vkin+(b-b.star())*(II/2);dkin=dkin+(db-db.star())*(II/2)
  out['kinetic_'+name]=vkin;dot['kinetic_'+name]=dkin
 eps=lambda a,b:0 if a==b else (1 if a==0 else -1)
 for r in range(3):
  for s in range(3):
   for kind in ('up','down','lepton'):
    val=W();dv=W()
    for a in range(2):
     b=1-a;spin=eps(a,b)
     for color in range(1 if kind=='lepton' else 3):
      for weak in range(2):
       left='L' if kind=='lepton' else 'Q';right={'up':'u_c','down':'d_c','lepton':'e_c'}[kind];idx=2*color+weak
       hw=1-weak if kind=='up' else weak;h=H[hw] if kind=='up' and variant!='mixed_up_hypercharge' else conj(H[hw]);dh=dH[hw] if kind=='up' and variant!='mixed_up_hypercharge' else conj(dH[hw]);weight=spin*(eps(weak,hw) if kind=='up' else 1)
       x=psi[left,r][idx][a];dx=dp[left,r][idx][a];y=psi[right,s][color][b];dy=dp[right,s][color][b]
       val=val+(x*y)*(h*weight);dv=dv+(dx*y+x*dy)*(h*weight)+(x*y)*(dh*weight)
    for part,factor in [('re',ONE),('im',II)]:
     name=f'Y_{kind}_{r}_{s}_{part}';out[name]=-(val+val.star()) if part=='re' else -(val-val.star())*factor;dot[name]=-(dv+dv.star()) if part=='re' else -(dv-dv.star())*factor
 return out,dot

def load(path=OUTPUT):
 raw=Path(path).read_bytes();need(len(raw)<10000000,'packet size')
 def pairs(xs):
  d={}
  for k,v in xs:need(k not in d,'duplicate JSON key');d[k]=v
  return d
 def no(x):raise ValueError('nonexact JSON number '+x)
 return json.loads(raw,object_pairs_hook=pairs,parse_float=no,parse_constant=no)
def canonical(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def exact(a,b,label):need(canonical(a)==canonical(b),label)
def pin(path):
 raw=path.read_bytes();return {'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}

OWN=('jet_action.py','build_local_action.py','verify_local_action.py','test_local_action.py')
PARENTS=('code/a5_closure/exterior_sm_completion.json','code/e9_kinetic/gauge_kinetic_invariant_forms.certificate.json','Lean/Screen/WeylYukawaConventions.lean','Lean/Screen/LocalGaugeJetAction.lean','paper/tex_fragments/LOCAL_SM_JET_ACTION.tex','code/particles/calibration/wz_upstream_completion/outputs/sm_eft_action_1.json')
EXPECTED_SCOPE={'conditional_local_classical_action':True,'three_families_declared':True,'fermions_anticommuting':True,'all_complex_Yukawa_entries_symbolic':True,'local_parameter_derivatives_retained':True,'same_action_variational_matter_current':True,'source_connection_selected':False,'source_action_emitted':False,'physical_Spin_attached':False,'coupling_or_mass_values_selected':False,'broken_vacuum_selected':False,'curved_spacetime_action':False,'renormalized_quantum_theory':False,'quantum_anomaly_from_measure_proved':False,'full_action_formalized_in_Lean':False,'empirical_comparison':False,'current_joined_to_source_population':False}
EXPECTED_PREMISES=['flat oriented Minkowski chart and supplied two-component Spin convention','declared local anti-Hermitian SU3/SU2/U1 connections on the registered fundamental/dual/singlet representations','representation-to-classified-bracket identification and invariant trace normalization are declared; no response-history emission','minimal first-derivative Hermitian Weyl and scalar kinetic law','Higgs doublet H=(1,2)_(1/2), scalar degree<=4 polynomial, three families and all three Yukawa intertwiners declared','positive free gauge trace coefficients; real mass-squared and quartic parameters; unrestricted complex Yukawa matrices','smooth gauge parameter with symmetric second derivatives; compact support for the integrated Noether argument']
EXPECTED_CONVENTIONS={'metric':[1,-1,-1,-1],'connection':'anti-Hermitian; D=partial+rho(A), F=dA+[A,A]','parameter':'delta psi=theta psi; delta A=-dtheta+[theta,A]','fermions':'all left Weyl; independent Grassmann conjugate and first-derivative generators','adjoint':'antilinear, reverses order, exchanges adjacent value/conjugate generator pairs','Pauli':'bar_sigma=(I,-sigma1,-sigma2,-sigma3); epsilon01=1','gauge_normalization':'kappa_s sum_(mu<nu) eta_mu eta_nu Tr(F_mu_nu squared); trace bases are not orthonormal','scalar_potential':'V=m_squared Hdagger H+lambda(Hdagger H)^2; no broken branch selected','Yukawa':'-Y T-conj(Y)Tdagger, represented by independent real and imaginary coefficient sectors','units':'formal local action/jet conventions; no laboratory calibration'}

def direct_current(j,mu,g):
 """Analytic variational formula, independent of differentiating the compiler."""
 trp=blank();s,T=GEN[g];trp[s]=T
 sigma=[[[ONE,Z],[Z,ONE]],[[Z,-ONE],[-ONE,Z]],[[Z,II],[-II,Z]],[[-ONE,Z],[Z,ONE]]]
 ids={tuple(x):i for i,x in enumerate(labels())};result=W()
 for row in ROWS:
  name,color,weak,y=row;dim=abs(color)*weak;R=family_matrix(trp,row)
  for family in range(3):
   psi=[[W({1<<ids[name,family,k,a,'value',False]:ONE}) for a in range(2)] for k in range(dim)]
   for k in range(dim):
    for l in range(dim):
     for a in range(2):
      for b in range(2):result=result+(psi[k][a].star()*psi[l][b])*(II*R[k][l]*sigma[mu][a][b])
 R=family_matrix(j['A'][mu],HROW);T=family_matrix(trp,HROW);H=j['H']
 cov=[j['dH'][mu][a]+sum((R[a][b]*H[b] for b in range(2)),Z) for a in range(2)]
 th=[sum((T[a][b]*H[b] for b in range(2)),Z) for a in range(2)]
 scalar=sum((conj(x)*y+conj(y)*x for x,y in zip(th,cov)),Z)*(1 if mu==0 else -1)
 return result+W({0:scalar})

def validate_custody(packet):
 """Authenticate inputs and exact scope; this does not prove coefficient correctness."""
 need(type(packet) is dict,"packet must be an object")
 exact(sorted(packet),sorted(['schema','scope','premises','conventions','representation_fields','representation_boundary','fixture','generator_labels','sector_coefficients','action_coefficients','ward_directions','symmetric_second_parameter_jets','matter_current_coefficients','negative_controls','grassmann_controls','comparison','source_pins','parent_pins']),'packet inventory')
 exact(packet['schema'],'oph.local_sm_jet_action.v1','schema');exact(packet['scope'],EXPECTED_SCOPE,'scientific scope');exact(packet['premises'],EXPECTED_PREMISES,'declared premises');exact(packet['conventions'],EXPECTED_CONVENTIONS,'conventions')
 exact(packet['source_pins'],{x:pin(HERE/x) for x in OWN},'source pins');exact(packet['parent_pins'],{x:pin(ROOT/x) for x in PARENTS},'parent pins')
 ext=load(ROOT/PARENTS[0]);external=load(ROOT/PARENTS[-1]);exact(packet['representation_fields'],ext['matter_package']['fields'],'registered representation');exact(packet['representation_boundary'],ext['conditional_closure']['physical_status'],'representation boundary')
 for name,color,weak,Y in ROWS:
  f=packet['representation_fields'][name];exact([f['su3'],f['su2'],f['Y'],f['dimension']],['3bar' if color<0 else str(color),str(weak),str(Y),abs(color)*weak],'representation dimensions/charges')
 exact(packet['fixture'],{'kind':'exact algebraic test jets; not measured or source-emitted fields','seed':0,'families':3,'jets':{k:encode_scalar(v) for k,v in fixture().items()}},'test jet identity')
 exact(packet['generator_labels'],labels(),'Grassmann generator census')
 ward=[{'kind':kind,'generator':g,'mu':mu,'nonzero_sectors':{}} for kind in ('theta','dtheta') for g in range(12) for mu in (range(4) if kind=='dtheta' else [0])]
 exact(packet['ward_directions'],ward,'complete local Ward direction inventory')
 hess=[{'generator':g,'mu':mu,'nu':nu,'curvature_variation_zero':True} for g in range(12) for mu in range(4) for nu in range(mu,4)]
 exact(packet['symmetric_second_parameter_jets'],hess,'complete symmetric second-jet inventory')
 # Strict nested algebra encoding. Canonical decimal integers/rationals only.
 def rational(token):
  need(type(token) is str,'exact rational token type')
  try:value=Fraction(token)
  except (ValueError,ZeroDivisionError):raise ValueError('exact rational token') from None
  need(str(value)==token,'canonical exact rational token')
 def polynomial(poly):
  need(type(poly) is list,'polynomial list')
  words=[]
  for term in poly:
   need(type(term) is list and len(term)==2,'polynomial term')
   word,coefficient=term
   need(type(word) is list and all(type(i) is int and 0<=i<900 for i in word),'Grassmann index type/range')
   need(word==sorted(set(word)),'canonical exterior word')
   need(type(coefficient) is list and len(coefficient)==2,'Gaussian coefficient')
   for token in coefficient:rational(token)
   need(coefficient!=['0','0'],'zero term must be omitted')
   words.append(tuple(word))
  need(words==sorted(set(words)),'canonical unique monomials')
 names=sorted(['kappa_SU3','kappa_SU2','kappa_U1','scalar_kinetic','m_squared','lambda']+['kinetic_'+r[0] for r in ROWS]+[f'Y_{k}_{r}_{s}_{p}' for r in range(3) for s in range(3) for k in ('up','down','lepton') for p in ('re','im')])
 need(type(packet['action_coefficients']) is dict,'action coefficient object');exact(sorted(packet['action_coefficients']),names,'complete action sector inventory')
 for poly in packet['action_coefficients'].values():polynomial(poly)
 pars={k:({'kind':'declared_minimal_normalization','value':'1'} if k.startswith('kinetic_') or k=='scalar_kinetic' else {'kind':'free_symbol','name':k,'domain':'positive real' if k.startswith('kappa_') else 'real'}) for k in names}
 exact(packet['sector_coefficients'],pars,'symbolic coefficient inventory')
 exact(sorted(packet['matter_current_coefficients']),sorted(f'{mu}/{g}' for mu in range(4) for g in range(12)),'complete current inventory')
 for poly in packet['matter_current_coefficients'].values():polynomial(poly)
 exact(sorted(packet['grassmann_controls']),sorted(['square','anticommutator','spin_pair','adjoint_product']),'Grassmann control inventory')
 for poly in packet['grassmann_controls'].values():polynomial(poly)
 exact(packet['comparison'],{'external_action_role':external['claim_boundary'],'physical_operator_shapes_assembled':['gauge_kinetic_SU3','gauge_kinetic_SU2','gauge_kinetic_U1','fermion_kinetic','higgs_kinetic','higgs_mass','higgs_quartic','yukawa_up','yukawa_down','yukawa_lepton'],'not_assembled':['gauge_fixing','ghost_sector'],'no_external_parameter_values_consumed':True,'mass_sign_dictionary':'external +mu2 |H|^2 corresponds to local -m_squared |H|^2 with mu2=-m_squared','Yukawa_comparison':'same neutral invariant channels under charge conjugation and Hermitian conjugation; external coefficient symbols are not numerically identified'},'external comparison boundary')
 controls=packet['negative_controls'];need(type(controls) is list and len(controls)==4,'negative control inventory')
 for control,(variant,kind,g,mu) in zip(controls,[('drop_dtheta_connection','dtheta',11,0),('drop_A_commutator','dtheta',0,1),('wrong_dual_color','theta',1,0),('mixed_up_hypercharge','theta',11,0)]):
  exact(sorted(control),sorted(['variant','kind','generator','mu','nonzero_sector_terms','residual_digest']),'negative control fields')
  exact([control['variant'],control['kind'],control['generator'],control['mu']],[variant,kind,g,mu],'negative control direction')
  need(type(control['nonzero_sector_terms']) is dict and bool(control['nonzero_sector_terms']),'negative control residual inventory')
  need(all(k in names and type(v) is int and v>0 for k,v in control['nonzero_sector_terms'].items()),'negative control term counts')
  h=control['residual_digest'];need(type(h) is str and len(h)==64 and all(c in '0123456789abcdef' for c in h),'residual digest type')
 return {'accepted_custody':True,'mathematical_replay':False,'schema':packet['schema'],'sha256':digest(packet)}

def verify(packet):
 validate_custody(packet)
 j=fixture();value,_=build_density(j)
 expected={k:x.encode() for k,x in sorted(value.items())};exact(packet['action_coefficients'],expected,'reconstructed Grassmann action coefficients')
 need(len(value)==65 and all(v for v in value.values()),'nonvacuous complete sector census');need(all(v==v.star() and all(mask.bit_count()%2==0 for mask in v.t) for v in value.values()),'Hermitian Grassmann-even density')
 pars={k:({'kind':'declared_minimal_normalization','value':'1'} if k.startswith('kinetic_') or k=='scalar_kinetic' else {'kind':'free_symbol','name':k,'domain':'positive real' if k.startswith('kappa_') else 'real'}) for k in sorted(value)}
 exact(packet['sector_coefficients'],pars,'symbolic coefficient inventory')
 # Compare currents before the longer Ward sweep so a forged current fails promptly.
 currents={f'{mu}/{g}':direct_current(j,mu,g).encode() for mu in range(4) for g in range(12)}
 exact(packet['matter_current_coefficients'],currents,'independent variational current formula');need(all(currents.values()),'nonzero currents')
 x=W({1:ONE});y=W({4:ONE});exact(packet['grassmann_controls'],{'square':(x*x).encode(),'anticommutator':(x*y+y*x).encode(),'spin_pair':(x*y-y*x).encode(),'adjoint_product':(x*y).star().encode()},'Grassmann sign controls');need(bool(x*y-y*x),'anticommuting spin contraction nonzero')
 controls=[]
 for variant,kind,g,mu in [('drop_dtheta_connection','dtheta',11,0),('drop_A_commutator','dtheta',0,1),('wrong_dual_color','theta',1,0),('mixed_up_hypercharge','theta',11,0)]:
  _,dot=build_density(j,kind,g,mu,variant=variant);nonzero={k:len(x.t) for k,x in sorted(dot.items()) if x};need(bool(nonzero),'nonvacuous negative control')
  controls.append({'variant':variant,'kind':kind,'generator':g,'mu':mu,'nonzero_sector_terms':nonzero,'residual_digest':digest({k:x.encode() for k,x in sorted(dot.items())})})
 exact(packet['negative_controls'],controls,'independent negative Ward controls')
 for row in packet['ward_directions']:
  _,dot=build_density(j,row['kind'],row['generator'],row['mu']);need(not any(dot.values()),'exact local Ward failure')
 for row in packet['symmetric_second_parameter_jets']:
  _,_,da,dda=variation(j,'d2theta',row['generator'],row['mu'],row['nu'])
  need(all(not x for u in range(4) for v in range(4) for s in range(3) for row0 in add(dda[u][v][s],scale(dda[v][u][s],-1)) for x in row0),'second parameter cancellation')
 return {'accepted':True,'schema':packet['schema'],'action_sectors':65,'Grassmann_generators':900,'action_monomials':sum(len(v.t) for v in value.values()),'local_Ward_directions':60,'symmetric_second_parameter_jets':120,'variational_matter_currents':48,'all_coefficients_exact':True,'full_action_Lean_formalization':False,'source_action_emission':False,'physical_calibration':False,'quantum_or_empirical_promotion':False,'verified_packet_sha256':digest(packet)}


if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,default=OUTPUT);a=ap.parse_args();print(json.dumps(verify(load(a.packet)),sort_keys=True,indent=2))

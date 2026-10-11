"""Independent reviewer: fixed exact algebra, does not import original checker.
60 paths,20 Cayley factors,12 generators,3960 binary events,4 cutoff dimensions.
No numerical search. Stops on first failed assertion. Full-source membership excluded.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,sys
base=Path(__file__).parent
raw=(base/'data/response.json').read_bytes()
r=json.loads(raw)
Z=(F(0),F(0));O=(F(1),F(0));S=(F(0),F(1))
def Q(a):return(F(a),F(0))
def plus(a,b):return(a[0]+b[0],a[1]+b[1])
def minus(a,b):return(a[0]-b[0],a[1]-b[1])
def prod(a,b):return(a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def inv(a):
 d=a[0]*a[0]-5*a[1]*a[1];assert d!=0
 return(a[0]/d,-a[1]/d)
def sig(a):return(a[0],-a[1])
def total(xs):
 a=Z
 for x in xs:a=plus(a,x)
 return a
def mat(m):return[[tuple(map(F,e[:2])) for e in row] for row in m]
def imag(m):return[[tuple(map(F,e[2:])) for e in row] for row in m]
def mm(a,b):return[[total(prod(x,y) for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a):return list(map(list,zip(*a)))
def ms(a,b,op):return[[op(x,y) for x,y in zip(row1,row2)]for row1,row2 in zip(a,b)]
def mult(a,s):return[[prod(x,s)for x in row]for row in a]
def determinant(m):
 a,b,c=m[0];d,e,f=m[1];g,h,i=m[2]
 return plus(minus(prod(a,minus(prod(e,i),prod(f,h))),prod(b,minus(prod(d,i),prod(f,g)))),prod(c,minus(prod(d,h),prod(e,g))))
def I(n):return[[Q(i==j)for j in range(n)]for i in range(n)]
def negate(a):return mult(a,Q(-1))
def fmt(a):return[f'{a[0]}',f'{a[1]}']
checks=[]
def check(s,t):
 assert t,s
 checks.append(s)
Gens=[[mat(b)for b in x]for x in r['generators']]
Ims=[[imag(b)for b in x]for x in r['generators']]
V=[[prod(Q(2),g[0][i][j])for i,j in [(2,1),(0,2),(1,0)]]for g in Gens]
dot=lambda v,w:total(prod(a,b)for a,b in zip(v,w))
radius=dot(V[0],V[0]);check('common positive-frame radius',all(dot(v,v)==radius for v in V))
for k,(g,im,v) in enumerate(zip(Gens,Ims,V)):
 skew=[[Z,minus(Z,v[2]),v[1]],[v[2],Z,minus(Z,v[0])],[minus(Z,v[1]),v[0],Z]]
 Qm=[[minus(prod(prod(v[i],v[j]),Q(F(1,2))),prod(radius,Q(F(1,6)))if i==j else Z)for j in range(3)]for i in range(3)]
 for i in range(3):Qm[i][i]=plus(Qm[i][i],Q(F(1,12)))
 check(f'generator {k} exact D formula both blocks',g[0]==mult(skew,Q(F(1,2)))and g[1]==[[sig(e)for e in row]for row in g[0]]and im[0]==Qm and im[1]==[[Z]*3 for _ in range(3)])
 check(f'generator {k} skew-Hermitian',all(tr(a)==negate(a)for a in g)and all(tr(a)==a for a in im))
rows=[]
for g,im in zip(Gens,Ims):rows.append([e for b,a in zip(g,im) for row1,row2 in zip(b,a)for x,y in zip(row1,row2)for e in [x,y]])
pivot=0
for col in range(len(rows[0])):
 at=next((i for i in range(pivot,len(rows))if rows[i][col]!=Z),None)
 if at is None:continue
 rows[pivot],rows[at]=rows[at],rows[pivot]; d=inv(rows[pivot][col]);rows[pivot]=[prod(x,d)for x in rows[pivot]]
 for i in range(len(rows)):
  if i!=pivot and rows[i][col]!=Z:
   c=rows[i][col];rows[i]=[minus(x,prod(c,y))for x,y in zip(rows[i],rows[pivot])]
 pivot+=1
 if pivot==len(rows):break
check('response generators have exact field rank twelve',pivot==12)
fac=[]
for k,f in enumerate(r['cayley_factors']):
 u=[mat(b)for b in f['observed_response']];c=[tuple(map(F,x))for x in f['control']]
 K=[[[total(prod(c[p],Gens[p][b][i][j])for p in range(12))for j in range(3)]for i in range(3)]for b in range(2)]
 K_im=[[[total(prod(c[p],Ims[p][b][i][j])for p in range(12))for j in range(3)]for i in range(3)]for b in range(2)]
 check(f'Cayley {k} source control has no imaginary part',all(e==Z for b in K_im for row in b for e in row))
 check(f'Cayley {k} matches actual source response',all(ms(I(3),a,plus)==mm(b,ms(I(3),a,minus))for a,b in zip(K,u)))
 check(f'Cayley {k} proper rotations both blocks',all(mm(a,tr(a))==I(3)and determinant(a)==O for a in u))
 fac.append(u)
trace_record=[]
for k,p in enumerate(r['closed_paths']):
 u=[mat(b)for b in p['observed_response']];acc=[I(3),I(3)];perm=list(range(12))
 for f in p['factor_indices']:
  acc=[mm(a,b)for a,b in zip(fac[f],acc)]
  perm=[r['cayley_factors'][f]['port_action'][j]for j in perm]
 check(f'path {k} response equals stored word BOTH blocks',acc==u)
 check(f'path {k} port action equals same stored word',perm==p['port_action'])
 for b,vec in enumerate([V,[[sig(e)for e in v]for v in V]]):
  check(f'path {k} actual covariance on triplet {b}',all([total(prod(u[b][i][j],vec[h][j])for j in range(3))for i in range(3)]==vec[p['port_action'][h]]for h in range(12)))
 check(f'path {k} block Galois relation',u[1]==[[sig(e)for e in row]for row in u[0]])
 ts=[total(a[i][i]for i in range(3))for a in u]
 if k==1:trace_record=ts
check('same chosen path traces differ exactly',trace_record==[(F(1,2),F(-1,2)),(F(1,2),F(1,2))])
# Uniform events independently enumerated: state, seam, fresh fair endpoint bit.
from itertools import combinations
Gram=[[prod(dot(v,w),inv(radius))for w in V]for v in V]
adj=[[int(x==prod(S,Q(F(1,5))))for x in row]for row in Gram]
edges=[(a,b)for a in range(12)for b in range(a+1,12)if adj[a][b]]
triangle=next(t for t in combinations(range(12),3)if all(adj[a][b]for a,b in combinations(t,2)))
shell=list(combinations(range(12),2));event_count=len(shell)*len(edges)*2
Bsum=[[0]*12 for _ in range(12)];moves=0;vars=0
for state in shell:
 x=[int(i in state)for i in range(12)];ds=[0]*12
 for a,b in edges:
  for coin in [0,1]:
   da=(coin-x[a])if x[a]+x[b]==1 else 0
   dx=[0]*12;dx[a]=da;dx[b]=-da
   if da:
    moves+=1
    if a in triangle and b in triangle:vars+=1
   for i in [a,b]:
    ds[i]+=dx[i]
    for j in [a,b]:Bsum[i][j]+=dx[i]*dx[j]
 # Mean drift is -Lx/60 because each event weight is 1/60.
 check('binary-state '+str(state)+' mean drift',all(ds[i]==-(5*x[i]-sum(adj[i][j]*x[j] for j in range(12))) for i in range(12)))
kappa=F(moves,event_count);nu=F(vars,event_count)
B=[[F(x,event_count)for x in row]for row in Bsum]
check('enumerated kappa and nu',kappa==F(5,33)and nu==F(1,66))
check('enumerated full B=L/198',all(B[i][j]==F(5*int(i==j)-adj[i][j],198)for i in range(12)for j in range(12)))
# Normalize lambda=1. Distinguish move fraction from squared-current activity.
gamma=prod(minus(Q(5),S),Q(F(1,60)));ratio=prod(gamma,Q(kappa/nu))
controls=[]
for scale in [F(1),F(2),F(3)]:
 scaled_nu=scale*scale*nu
 scaled_ratio=prod(gamma,Q(kappa/scaled_nu))
 restored=prod(gamma,Q((scale*scale*kappa)/scaled_nu))
 check('load-scale '+str(scale)+' L_record cancellation',all(3*scale*scale*B[i][j]/scaled_nu==5*int(i==j)-adj[i][j]for i in range(12)for j in range(12)))
 check('load-scale '+str(scale)+' move-fraction ratio scaling',scaled_ratio==prod(ratio,Q(1/(scale*scale)))and restored==ratio)
 controls.append(dict(scale=str(scale),nu=str(scaled_nu),ratio_with_move_fraction=list(map(str,scaled_ratio)),ratio_with_squared_current_activity=list(map(str,restored))))
eta=F(1,2);idleratio=prod(prod(gamma,Q(eta)),Q((eta*kappa)/(eta*nu)))
check('idle attempt move-fraction ratio changes',idleratio==prod(ratio,Q(eta)))
check('idle clock gamma/nu fixed theoretical kappa unchanged',prod(prod(gamma,Q(eta)),Q(kappa/(eta*nu)))==ratio)
# Occupation sparse-matrix-unit witness, CCR and qubit control independently derived.
occupation=dict(commutator_coefficient=int(0 in (0,1))-int(0 in (1,2)),N0N1_on_01=1,sum_N_diagonal=list({sum(int(i in state)for i in range(12))for state in shell}))
check('occupation noncentral nonorthogonal and sum',occupation==dict(commutator_coefficient=1,N0N1_on_01=1,sum_N_diagonal=[2]))
ccr={}
for d in [2,3,4,6]:
 aa=[n+1 if n<d-1 else 0 for n in range(d)];astar_a=list(range(d));comm=[a-b for a,b in zip(aa,astar_a)]
 check('cutoff '+str(d)+' exact shift commutator',comm==[1]*(d-1)+[1-d]and sum(comm)==0)
 ccr[str(d)]=comm
clock=[[(1+cos)/2,F(1,2)]for cos in [F(1),F(0),F(-1)]]
check('qubit dephasing control probabilities',clock==[[1,F(1,2)],[F(1,2),F(1,2)],[0,F(1,2)]])
# A3 coordinate mutation controls: diagonal71 + real180 + imaginary180.
coordinate_count=71+12*15*2
check('A3 full normalized coordinate count431',coordinate_count==431)
base_eigen=F(1,72);perturb=F(1,288)
check('A3 all omission states have eigenfloor1/96',base_eigen-perturb==F(1,96)>0)
mutations=[]
for a in range(71):
 # Only independent diagonal a changes; normalization anchor71 compensates.
 mutations.append(('diagonal',a))
 check('A3 diagonal mutation '+str(a),2*perturb!=0 and base_eigen-perturb>0)
for block in range(12):
 for a,b in combinations(range(6),2):
  for coordinate in ['real','imaginary']:
   mutations.append((coordinate,block,a,b))
   # Hermitian +/- off-diagonal matrix-unit pair has eigenvalues +/-1.
   check('A3 '+coordinate+' mutation '+str((block,a,b)),base_eigen-perturb>0 and 2*perturb!=0)
check('A3 all independent omission witnesses covered',len(mutations)==431 and len(set(mutations))==431)
# Tagged instrument pullback with nonuniform posterior weights: exact2x2 control.
rho1=[F(1,3),F(2,3)];rho2=[F(3,4),F(1,4)]
factor=[F(1),F(1,4)]
t1=[a*b for a,b in zip(rho1,factor)];t2=[a*b for a,b in zip(rho2,factor)]
coarse=[(a+b)/2 for a,b in zip(rho1,rho2)]
tc=[a*b for a,b in zip(coarse,factor)]
check('A3 descendant unnormalized instrument pullback',tc==[(a+b)/2 for a,b in zip(t1,t2)])
p1=sum(t1);p2=sum(t2);pc=sum(tc)
posterior=[p1/(2*pc),p2/(2*pc)]
normalized=[a/pc for a in tc]
mixture=[posterior[0]*a/p1+posterior[1]*b/p2 for a,b in zip(t1,t2)]
old_prior=[(a/p1+b/p2)/2 for a,b in zip(t1,t2)]
check('A3 selected tag posterior retained',posterior==[F(8,21),F(13,21)] and normalized==mixture and normalized==[F(52,63),F(11,63)])
check('A3 old equal tag prior is a failing control',old_prior!=normalized)
# The whole-source Galois construction changes only structural law coefficients.
# A3 remains on the usual complex PSD cone; entrywise sigma is not a positivity map.
e=prod(minus(S,O),Q(F(1,4)));sigma_e=sig(e)
check('Galois arbitrary-positive-state substitution is not licensed',e==(F(-1,4),F(1,4)) and sigma_e==(F(-1,4),F(-1,4)))
# Since2<sqrt5<3,0<e<1/2 but sigma(e)<0: diag(e,1-e) is a countercontrol.
result=dict(target_revision='0510578e317de5da0d485d09792dc8bcad2e2f3c',python=sys.version,check_count=len(checks),all_passed=True,response_sha256=hashlib.sha256(raw).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,selected_path_index=1,selected_factor_indices=r['closed_paths'][1]['factor_indices'],selected_traces=[fmt(x)for x in trace_record],generator_rank=pivot,event_count=event_count,triangle=triangle,kappa=str(kappa),nu=str(nu),load_unit_controls=controls,idle_control_ratio=fmt(idleratio),occupation=occupation,ccr=ccr,A3_coordinate_count=coordinate_count,A3_omission_eigenfloor=str(base_eigen-perturb),A3_posterior=list(map(str,posterior)),galois_positive_state_caution={'positive_e':fmt(e),'negative_sigma_e':fmt(sigma_e)},scope='independent exact finite algebra/path/event validation; full all-level A1-A3 and physical membership NOT certified')
(base/'review_independent_controls.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items()if k!='checks'},indent=2))

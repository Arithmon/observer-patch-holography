"""Exact 3-site control; finite mechanics only, not a full A1--A3 model.
Target: revision 0510578e317de5da0d485d09792dc8bcad2e2f3c, MODEL_SEARCH 1331--1412.
Fixed budget: 8 cases, 3x3 rational products, no data, no numerical search.
Stop after prescribed cases; any failed assertion stops with nonzero exit.
Premises: M=I; A=I+path3 Laplacian; canonical symplectic reading;
H_e=(p^T(I+eA)p+q^TAq)/2, fixed position reading and common clock.
The same-capacity inertia premise is deliberately absent (violated if e!=0).
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
A=[[F(2),F(-1),F(0)],[F(-1),F(3),F(-1)],[F(0),F(-1),F(2)]]
I=[[F(i==j) for j in range(3)] for i in range(3)]
A2=[[sum(A[i][k]*A[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
assert A2[0][2]==1
# Exact orthogonal eigendirections verify the stiffness spectrum 1,2,4.
for l,v in [(1,[1,1,1]),(2,[1,0,-1]),(4,[1,-2,1])]:
 assert all(sum(A[i][j]*v[j] for j in range(3))==l*v[i] for i in range(3))
rows=[]
for e in [F(0),F(1,16),F(-1,16),F(-1,4)]:
 for g in [F(1),F(3,2)]:
  eigB=[1+e*l for l in [1,2,4]]
  positive=all(x>0 for x in eigB)
  assert positive==(e!=F(-1,4))
  Y=[[g*(A[i][j]+e*A2[i][j]) for j in range(3)] for i in range(3)]
  assert Y[0][2]==g*e
  vals=[g*(l+e*l*l) for l in [1,2,4]]
  beta=vals[1]-vals[0]; alpha=vals[0]-beta
  residual=vals[2]-(alpha+4*beta)
  assert residual==6*g*e
  assert (residual==0)==(e==0)
  # The distance-two nonzero response cannot occur for any alpha I+beta A.
  assert (Y[0][2]==0)==(e==0)
  rows.append(dict(epsilon=str(e),clock_multiplier=str(g),kinetic_eigenvalues=list(map(str,eigB)),positive_kinetic=positive,distance_two_acceleration=str(Y[0][2]),fitted_alpha=str(alpha),fitted_beta=str(beta),fitted_mass_squared_on_L=str(alpha+beta),third_mode_squared_frequency_residual=str(residual)))
result=dict(target_revision='0510578e317de5da0d485d09792dc8bcad2e2f3c',scope='exact finite mechanics, NOT complete A1--A3 source membership',A=[[str(x) for x in r] for r in A],cases=rows,assertions_passed=True,controls='epsilon=0 affine baseline; epsilon=-1/4 excluded singular inertia; two common clock multipliers; fixed position/stiffness identity',evidence='rational finite algebra; normalized classical configuration Gibbs invariance follows analytically by integrating p; full phase-space reference and quantum vacuum change',independent_review='passed analytic challenge and exact independent replay; see provenance/review_status.json',stopping_rule='all 8 prescribed cases, no adaptive expansion')
print(json.dumps(result,indent=2))

"""Bounded exact population/record experiment; no physical selection claimed."""
from __future__ import annotations
import hashlib,json
from fractions import Fraction
from itertools import product
from math import isqrt
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SPEC=HERE/'specification.json'
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
PORTS=(0,1,4,5,8,9)
def raw(d):return (json.dumps(d,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def digest(d):return hashlib.sha256(raw(d)).hexdigest()
def write(p,d):p.write_bytes(raw(d))
def floor_phi(b):return (b+isqrt(5*b*b))//2

def tables():
 edges=sorted({tuple(sorted((f[k],f[(k+1)%3]))) for f in FACES for k in range(3)})
 dist=np.full((12,12),99,dtype=np.int64);np.fill_diagonal(dist,0)
 for a,b in edges:dist[a,b]=dist[b,a]=1
 for k in range(12):dist=np.minimum(dist,dist[:,k,None]+dist[None,k,:])
 anti=[int(np.flatnonzero(dist[i]==3)[0]) for i in range(12)]
 return edges,anti

def address(loads,anti):
 z=[loads[p]-loads[anti[p]] for p in PORTS]
 return [[(z[1]-z[0])/2,(z[4]+z[5])/2],[(z[3]-z[2])/2,(z[0]+z[1])/2],[(z[5]-z[4])/2,(z[2]+z[3])/2]]

def phi_positive(a,b):
 # A+B phi = (2A+B+B sqrt(5))/2. Integer products stay below2^63 here.
 c=2*a+b
 return ((c>=0)&(b>=0)&((c!=0)|(b!=0))) | ((c>0)&(b<0)&(c*c>5*b*b)) | ((c<0)&(b>0)&(5*b*b>c*c))

def relation(points,q):
 # Coordinates encoded by integer coefficients of 1/4 and phi/4.
 n=len(points);answer=np.zeros((n,n),dtype=bool);points=np.asarray(points,dtype=np.int64)
 for start in range(0,n,128):
  dif=points[start:start+128,None,:,:]-points[None,:,:,:]
  aa=dif[:,:,:,0];bb=dif[:,:,:,1]
  A=((aa*aa+bb*bb).sum(axis=2))*q-16;B=((2*aa*bb+bb*bb).sum(axis=2))*q
  answer[start:start+128]=~phi_positive(A,B)
 return answer

def make():
 spec=json.loads(SPEC.read_text());assert spec['q']==13
 edges,anti=tables();rows=[];events=[];previous='0'*64
 def event(payload):
  nonlocal previous
  item=dict(id=len(events),previous_hash=previous,**payload);item['event_hash']=digest(item);previous=item['event_hash'];events.append(item);return item['id']
 old=[];new=[]
 for site,b in enumerate(product(range(13),repeat=3)):
  m=[-floor_phi(i) for i in b]
  z=[b[1]-m[0],b[1]+m[0],b[2]-m[1],b[2]+m[1],b[0]-m[2],b[0]+m[2]]
  loads=[Fraction(0) for _ in range(12)]
  for k,p in enumerate(PORTS):loads[p]=Fraction(max(z[k],0));loads[anti[p]]=Fraction(max(-z[k],0))
  before=address(loads,anti);assert before==[[Fraction(m[i]),Fraction(b[i])] for i in range(3)]
  init=event({'op':'initialize','site':site,'parents':[],'loads':[str(x) for x in loads],'protected_address':[[str(v) for v in x] for x in before]})
  seam=next(((a,c) for a,c in edges if loads[a]!=loads[c]),None);after=loads.copy();repair=None
  if seam:
   a,c=seam;v=(loads[a]+loads[c])/2;after[a]=after[c]=v
   repair=event({'op':'pair_mean','site':site,'parents':[init],'seam':[a,c],'reads':{str(a):str(loads[a]),str(c):str(loads[c])},'writes':{str(a):str(v),str(c):str(v)},'quadratic_decrement':str((loads[a]-loads[c])**2/2)})
  live=address(after,anti)
  assert sum(loads)==sum(after)
  rows.append({'site':site,'b':list(b),'z':z,'init_event':init,'repair_event':repair,'live_address_after':[[str(v) for v in x] for x in live],'protected_address_after':[[str(v) for v in x] for x in before]})
  def encode(x):
   r=[[v*4 for v in pair] for pair in x];assert all(v.denominator==1 for pair in r for v in pair);return [[int(v) for v in pair] for pair in r]
  old.append(encode(before));new.append(encode(live))
 original=relation(old,13);updated=relation(new,13)
 changed=np.argwhere(np.triu(original!=updated,k=1));samples=[]
 for i,j in changed[:12]:
  samples.append({'sites':[int(i),int(j)],'before':bool(original[i,j]),'after':bool(updated[i,j])})
 pin='evidence/source_net_causal_poset/build_causal_poset.py'
 packet={'schema':'oph.source_population.v1','producer_reference':{'path':Path(__file__).relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},'specification_sha256':hashlib.sha256(SPEC.read_bytes()).hexdigest(),'source_reference':{'path':pin,'sha256':hashlib.sha256((ROOT/pin).read_bytes()).hexdigest()},'q':13,'seams':edges,'antipodes':anti,'events':events,'sites':rows,'final_event_hash':previous,'summary':{'sites':len(rows),'operations':len(events),'active_mean_repairs':sum(r['repair_event'] is not None for r in rows),'live_addresses_changed':sum(a!=b for a,b in zip(old,new)),'protected_addresses_changed':0,'distinct_live_addresses_after':len({tuple(v for pair in x for v in pair) for x in new}),'original_directed_metric_reads_with_self':int(original.sum()),'updated_directed_metric_reads_with_self':int(updated.sum()),'changed_unordered_metric_relations':len(changed)},'changed_relation_witnesses':samples,'scope':{'initial_population':'declared golden preparation, not derived by the pair-mean dynamics','protected_address':'immutable initialized memory; preservation follows from the specified write footprint','live_position':'readout of the mutable twelve-port load','repair_class':'one legal scalar pair mean per nonconstant carrier; no inter-carrier routing, drive, memory erasure or alternative-law exclusion','clock':'finite operation index, not physical time','physical_premise_discharged':False,'empirical_measurement':False}}
 return packet
if __name__=='__main__':
 p=make();out=HERE/'runtime/population_receipt.json';write(out,p);print(json.dumps(p['summary'],indent=2))

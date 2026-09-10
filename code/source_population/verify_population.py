"""Independent exact register replay and rational-bracket metric census."""
from __future__ import annotations
import hashlib,json,types
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'evidence/source_net_causal_poset/build_causal_poset.py'
SOURCE_SHA256='d82750cf95e368db9fee7fafd7905c4d235e62173701184a4a56bc1f33495aeb'
SPEC_SHA256='3223f8a1eb0a67024a1a0e9294e62fa80f9f3b3ff38b7c939f70ee3094a4d765'
PRODUCER=HERE/'build_population.py'
PRODUCER_SHA256='ac3b355a160613b66a61d9a40e96bf05928c9dfda5295973777d2dca516d5dd7'

def raw(d):return (json.dumps(d,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def need(ok,msg):
 if not ok:raise ValueError(msg)
def strict_pairs(pairs):
 d={}
 for k,v in pairs:need(k not in d,'duplicate key');d[k]=v
 return d

def load(path):
 return json.loads(path.read_text(encoding='utf-8'),object_pairs_hook=strict_pairs,parse_float=lambda v:(_ for _ in ()).throw(ValueError('floating number')),parse_constant=lambda v:(_ for _ in ()).throw(ValueError(v)))

def adjacency(points,q):
 """Certify every distance sign using rational enclosures of sqrt(5)."""
 points=np.asarray(points,dtype=np.int64);n=len(points);ans=np.empty((n,n),dtype=bool)
 D=2**40;lo=isqrt(5*D*D);hi=lo+1
 need(lo*lo<5*D*D<hi*hi,'irrational root bracket')
 for start in range(0,n,64):
  difference=points[start:start+64,None]-points[None]
  a=difference[...,0];b=difference[...,1]
  # Each coordinate is (a+b phi)/4, so use sqrt(5) basis directly.
  c=2*a+b
  C=q*((c*c+5*b*b).sum(axis=2))-64
  B=2*q*(c*b).sum(axis=2)
  need(int(abs(C).max())*D+int(abs(B).max())*hi<2**63,'integer multiplication bound')
  lower=C*D+np.where(B>=0,B*lo,B*hi)
  upper=C*D+np.where(B>=0,B*hi,B*lo)
  exact_zero=(C==0)&(B==0)
  need(np.all((lower>0)|(upper<0)|exact_zero),'distance sign unresolved')
  ans[start:start+64]=(upper<0)|exact_zero
 return ans

def check_shape(p):
 need(type(p)is dict and set(p)=={'schema','producer_reference','specification_sha256','source_reference','q','seams','antipodes','events','sites','final_event_hash','summary','changed_relation_witnesses','scope'},'root schema')
 need(p['schema']=='oph.source_population.v1','schema version')
 def walk(x,path=()):
  if type(x)is dict:
   for k,v in x.items():need(type(k)is str,'JSON key type');walk(v,path+(k,))
  elif type(x)is list:
   for i,v in enumerate(x):walk(v,path+(i,))
  elif type(x)is bool:
   need((path[:1]==('scope',) and path[-1] in {'physical_premise_discharged','empirical_measurement'}) or (path[:1]==('changed_relation_witnesses',) and path[-1] in {'before','after'}),'unexpected boolean')
  else:need(type(x) in (str,int,type(None)),'JSON value type')
 walk(p)
 for k in ('physical_premise_discharged','empirical_measurement'):need(type(p['scope'][k])is bool,'scope boolean type')
 for r in p['changed_relation_witnesses']:
  need(set(r)=={'sites','before','after'} and type(r['before'])is bool and type(r['after'])is bool,'witness schema')
 for r in p['sites']:
  need(set(r)=={'site','b','z','init_event','repair_event','live_address_after','protected_address_after'},'site schema')
 for e in p['events']:
  base={'id','previous_hash','event_hash','op','site','parents'}
  extra={'loads','protected_address'} if e['op']=='initialize' else {'seam','reads','writes','quadratic_decrement'}
  need(set(e)==base|extra,'event schema')
  if e['op']=='initialize':need(all(type(v)is str for v in e['loads']),'rational string type')

def check(p):
 check_shape(p)
 need(hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA256,'fixed source custody')
 need(hashlib.sha256((HERE/'specification.json').read_bytes()).hexdigest()==SPEC_SHA256,'fixed specification custody')
 need(hashlib.sha256(PRODUCER.read_bytes()).hexdigest()==PRODUCER_SHA256,'fixed producer custody')
 need(p['producer_reference']=={'path':PRODUCER.relative_to(ROOT).as_posix(),'sha256':PRODUCER_SHA256},'producer custody')
 spec=load(HERE/'specification.json')
 need(p['specification_sha256']==hashlib.sha256((HERE/'specification.json').read_bytes()).hexdigest(),'specification hash')
 need(p['q']==spec['q']==13,'fixed population size')
 need(p['source_reference']=={'path':SOURCE.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest()},'source custody')
 source=types.ModuleType('source_archive_for_population_check');source.__file__=str(SOURCE);exec(compile(SOURCE.read_bytes(),str(SOURCE),'exec'),source.__dict__)
 table=source.carrier_tables();edges=table['seams'];anti=table['antipode'];site_labels=source.site_coordinates(13,3)
 z,_,_,_=source.source_records(13,site_labels);initial=source.place_loads(z,anti)
 need(p['seams']==[list(x) for x in edges] and p['antipodes']==anti,'canonical incidence')
 need(len(p['sites'])==2197 and len(p['events'])==4393,'complete site/event census')
 # Coefficients of twice the source-axis readout, independently assembled
 # from its six antipodal basis differences. Apply the linear matrix to N.
 Z=np.zeros((6,12),dtype=int)
 for j,k in enumerate(source.POSITIVE_PORTS):Z[j,k]=1;Z[j,anti[k]]=-1
 C=np.zeros((3,2,12),dtype=int)
 C[0,0]=Z[1]-Z[0];C[0,1]=Z[4]+Z[5]
 C[1,0]=Z[3]-Z[2];C[1,1]=Z[0]+Z[1]
 C[2,0]=Z[5]-Z[4];C[2,1]=Z[2]+Z[3]
 def coord(N):return [[sum(F(int(C[i,j,k]))*N[k] for k in range(12))/2 for j in range(2)] for i in range(3)]
 records={};protected={};last={};seen={};previous='0'*64;active=0
 for i,e in enumerate(p['events']):
  need(e['id']==i and e['previous_hash']==previous,'serial history')
  payload={k:v for k,v in e.items() if k!='event_hash'}
  need(e['event_hash']==hashlib.sha256(raw(payload)).hexdigest(),'event hash');previous=e['event_hash']
  site=e['site'];need(type(site)is int and 0<=site<2197,'site ID')
  if e['op']=='initialize':
   need(site not in records and e['parents']==[],'unique initialization')
   N=list(map(F,e['loads']));need(N==list(map(F,initial[site])),'actual source preparation')
   address=coord(N);need(e['protected_address']==[[str(v) for v in x] for x in address],'initial exact address')
   records[site]=N;protected[site]=address;last[site]=i;seen[site]=[i,None]
  elif e['op']=='pair_mean':
   need(site in records and seen[site][1] is None,'one pair repair per site')
   N=records[site];a,b=e['seam'];need((a,b) in edges and N[a]!=N[b],'active allowed seam')
   need((a,b)==next((s,t) for s,t in edges if N[s]!=N[t]),'fixed seam choice')
   need(e['parents']==[last[site]],'read-from parent')
   need(e['reads']=={str(a):str(N[a]),str(b):str(N[b])},'read footprint')
   mean=(N[a]+N[b])/2;need(e['writes']=={str(a):str(mean),str(b):str(mean)},'pair mean')
   before_sum=sum(N);before_V=sum(x*x for x in N);N=N.copy();N[a]=N[b]=mean
   need(sum(N)==before_sum,'sum conservation');decrement=before_V-sum(x*x for x in N)
   need(decrement>0 and e['quadratic_decrement']==str(decrement),'strict potential descent')
   need(coord(N)!=protected[site],'live position changed')
   records[site]=N;last[site]=i;seen[site][1]=i;active+=1
  else:raise ValueError('unexpected operation')
 need(previous==p['final_event_hash'],'terminal authentication')
 old=[];new=[]
 def ints(x):
  values=[[v*4 for v in pair] for pair in x];need(all(v.denominator==1 for pair in values for v in pair),'address encoding denominator');return [[int(v) for v in pair] for pair in values]
 for i,r in enumerate(p['sites']):
  need(r['site']==i and r['b']==site_labels[i].tolist() and r['z']==z[i].tolist(),'source site identity')
  need([r['init_event'],r['repair_event']]==seen[i],'site history references')
  live=coord(records[i]);need(r['live_address_after']==[[str(v) for v in pair] for pair in live],'live readout join')
  need(r['protected_address_after']==[[str(v) for v in pair] for pair in protected[i]],'protected write footprint')
  if seen[i][1] is None:need(len(set(records[i]))==1,'only fixed state may omit repair')
  old.append(ints(protected[i]));new.append(ints(live))
 A=adjacency(old,13);B=adjacency(new,13);diff=np.argwhere(np.triu(A!=B,k=1))
 expected={'sites':2197,'operations':4393,'active_mean_repairs':active,'live_addresses_changed':sum(a!=b for a,b in zip(old,new)),'protected_addresses_changed':0,'distinct_live_addresses_after':len({tuple(v for pair in row for v in pair) for row in new}),'original_directed_metric_reads_with_self':int(A.sum()),'updated_directed_metric_reads_with_self':int(B.sum()),'changed_unordered_metric_relations':len(diff)}
 need(p['summary']==expected,'complete exact census')
 witnesses=[{'sites':[int(i),int(j)],'before':bool(A[i,j]),'after':bool(B[i,j])} for i,j in diff[:12]]
 need(p['changed_relation_witnesses']==witnesses,'metric witnesses')
 need(p['scope']=={'initial_population':'declared golden preparation, not derived by the pair-mean dynamics','protected_address':'immutable initialized memory; preservation follows from the specified write footprint','live_position':'readout of the mutable twelve-port load','repair_class':'one legal scalar pair mean per nonconstant carrier; no inter-carrier routing, drive, memory erasure or alternative-law exclusion','clock':'finite operation index, not physical time','physical_premise_discharged':False,'empirical_measurement':False},'scope boundary')
 return {'accepted':True,**expected,'metric_check':'all pairs, exact rational sqrt(5) enclosure; no producer import','source_population_produced':False,'persistent_memory_supplied':True}
if __name__=='__main__':print(json.dumps(check(load(HERE/'runtime/population_receipt.json')),indent=2))

#!/usr/bin/env python3
"""Independent graph and invariant-space verifier; no producer imports."""
from __future__ import annotations
import json
import hashlib
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
MANIFEST=ROOT/"code/a5_closure/manifests/echosahedral_federation_reference.json"
RECEIPTS={"a3":ROOT/"code/a5_closure/manifests/a3_scheduler_kernel_reference.json","directed":ROOT/"code/a5_closure/manifests/directed_seam_repair_reference.json","repair":ROOT/"code/a5_closure/manifests/record_counting_mechanism_reference.json"}
PINNED={"a3":"sha256:4afffa0760ff59649d1f2548d8f755f1efe4c5be14594bbde80f22d7ed3d536a","directed":"sha256:29d77ed315cba2e98a1f1bece56191edf8759a472e59799200716ee70163e4b8","repair":"sha256:efad022da8fdb52b58e6e8356ef945227832467c162cb512d0261812f286a655"}

def ck(f):
 a,b,c=f;return min((a,b,c),(b,c,a),(c,a,b))

def rank(rows):
 a=[[Fraction(x) for x in r] for r in rows];m=len(a);n=len(a[0]) if m else 0;r=0
 for c in range(n):
  q=next((i for i in range(r,m) if a[i][c]),None)
  if q is None:continue
  a[r],a[q]=a[q],a[r];p=a[r][c];a[r]=[x/p for x in a[r]]
  for i in range(m):
   if i!=r and a[i][c]:
    p=a[i][c];a[i]=[x-p*y for x,y in zip(a[i],a[r])]
  r+=1
 return r

def weighted_laplacian(n,edges,weights):
 L=[[Fraction(0) for _ in range(n)] for _ in range(n)]
 for (u,v),w in zip(edges,weights):
  if w<=0:raise AssertionError("proposal rates must be strictly positive")
  L[u][u]+=w;L[v][v]+=w;L[u][v]-=w;L[v][u]-=w
 return L

def proportional(a,b):
 ratio=None
 for i in range(len(a)):
  for j in range(len(a)):
   if b[i][j]:
    q=a[i][j]/b[i][j]
    if ratio is None:ratio=q
    elif ratio!=q:return False
   elif a[i][j]:return False
 return ratio is not None

def verify_receipts(receipts=None):
 docs={k:(json.loads(v.read_text()) if receipts is None else receipts[k]) for k,v in RECEIPTS.items()}
 expected={"a3":"oph.a3_scheduler_kernel_certificate.v1","directed":"oph.directed_seam_repair_certificate.v1","repair":"oph.record_counting_mechanism_certificate.v2"}
 def digest(d):return "sha256:"+hashlib.sha256(json.dumps(d,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()
 for k,d in docs.items():
  if d.get("schema")!=expected[k] or d.get("manifest_sha256")!=digest({x:y for x,y in d.items() if x!="manifest_sha256"}):raise AssertionError("receipt schema or self digest mismatch: "+k)
  if digest(d)!=PINNED[k]:raise AssertionError("receipt differs from audited baseline pin: "+k)
 a,b,r=docs["a3"],docs["directed"],docs["repair"]
 if a.get("issue")!=614 or r.get("issue")!=628:raise AssertionError("issue identity mismatch")
 sel=a["a3_selection"];dk=sel["deck_invariance"]
 if not (a["verdict"].get("kernel_status")=="A3_selected_exact" and a["move_simplex"].get("seams")==30 and sel["reference"]=="uniform_on_moves" and sel["feasible_set"]=="the full move simplex; the complete constraint grammar contributes normalization only" and sel["selected_kernel_probability_per_seam"]=="1/30" and sel["positive_weights"] and sel["kernel_sums_to_one"] and dk["edge_transitive"] and dk["seam_orbit_size"]==30):raise AssertionError("#614 selection receipt changed")
 pins=b["upstream_pins"]
 for key,receipt,suffix in (("undirected_scheduler",a,"a3_scheduler_kernel_reference.json"),("integer_record_counting_mechanism",r,"record_counting_mechanism_reference.json")):
  pin=pins[key]
  if pin["path"]!="manifests/"+suffix or pin["schema"]!=receipt["schema"] or pin["canonical_json_sha256"]!=digest(receipt):raise AssertionError("directed parent pin mismatch: "+key)
 s1=b["uniform_s1_channel"];v=b["verdict"];scope=b["scope"]["claim_boundary"]
 if s1.get("sector")!="S=1 total integer working load" or v.get("uniform_s1_channel")!="exact_named_realization":raise AssertionError("S=1 receipt boundary changed")
 if v.get("directed_orbit_schedule_as_physical_selection")!="open" or v.get("full_self_readback_and_universe_selection")!="open":raise AssertionError("source boundary flags changed")
 if "does not prove an IID" not in scope or "path law" not in scope:raise AssertionError("directed claim boundary changed")
 mech=r.get("mechanism",{})
 if "proposal_law" in mech or "proposal_weights" in mech:raise AssertionError("#628 acquired proposal-law data; reassess source conclusion")
 return {k:{"schema":d["schema"],"raw_sha256":hashlib.sha256(RECEIPTS[k].read_bytes()).hexdigest(),"canonical_json_sha256":digest(d)} for k,d in docs.items()}

def audit(doc=None):
 c=(json.loads(MANIFEST.read_text()) if doc is None else doc)["carrier"]; lab=c["ports"];ix={x:i for i,x in enumerate(lab)}
 if len(lab)!=12 or len(ix)!=12:raise AssertionError("twelve distinct port labels required")
 fs=[tuple(ix[x] for x in f) for f in c["oriented_faces"]]
 es={tuple(sorted(z)) for a,b,d in fs for z in ((a,b),(b,d),(d,a))};n=len(lab)
 # Both incidence representations belong to the input certificate. Rebuilding
 # a valid graph from faces must not hide corruption of the serialized seams.
 serialized=c["edges"]
 if not isinstance(serialized,list) or any(not isinstance(e,list) or len(e)!=2 for e in serialized):raise AssertionError("serialized seam shape mismatch")
 seams=[tuple(sorted((ix[a],ix[b]))) for a,b in serialized]
 if len(seams)!=30 or len(set(seams))!=30 or set(seams)!=es:raise AssertionError("serialized seams differ from oriented faces")
 A=[set() for _ in lab]
 for a,b in es:A[a].add(b);A[b].add(a)
 face={ck(f) for f in fs}; perms=[];m={};used=set()
 # Dynamic minimum-remaining-values vertex choice, with induced-edge checks.
 dom={v:{w for w in range(n) if len(A[v])==len(A[w])} for v in range(n)}
 def go():
  if len(m)==n:
   p=tuple(m[i] for i in range(n))
   if all(ck(tuple(p[x] for x in f)) in face for f in fs):perms.append(p)
   return
  v=min((x for x in range(n) if x not in m),key=lambda x:(len(dom[x]-used),x))
  for w in sorted(dom[v]-used):
   if all(((v in A[x])==(w in A[m[x]])) for x in m):
    m[v]=w;used.add(w);go();used.remove(w);del m[v]
 go()
 deg=[len(x) for x in A]; seen={0};todo=[0]
 while todo:
  u=todo.pop()
  for v in A[u]-seen:seen.add(v);todo.append(v)
 ei={e:j for j,e in enumerate(sorted(es))}; equations=[]
 for p in perms:
  for e,j in ei.items():
   q=tuple(sorted((p[e[0]],p[e[1]])))
   if e!=q:
    row=[0]*len(es);row[j]=1;row[ei[q]]=-1;equations.append(row)
 dim=len(es)-rank(equations)
 vorb={frozenset(p[v] for p in perms) for v in range(n)}
 rem=set(es);eorbs=[]
 while rem:
  e=min(rem);o={tuple(sorted((p[e[0]],p[e[1]]))) for p in perms};eorbs.append(o);rem-=o
 if len(fs)!=20 or len({ck(f) for f in fs})!=20:raise AssertionError("twenty oriented face census failed")
 result={"vertices":n,"edges":len(es),"degrees":sorted(set(deg)),"connected":len(seen)==n,"proper_group_order":len(perms),"vertex_orbit_sizes":sorted(map(len,vorb)),"edge_orbit_sizes":sorted(map(len,eorbs)),"invariant_weight_dimension":dim}
 if result!={"vertices":12,"edges":30,"degrees":[5],"connected":True,"proper_group_order":60,"vertex_orbit_sizes":[12],"edge_orbit_sizes":[30],"invariant_weight_dimension":1}:raise AssertionError(result)
 # Independently reconstruct both proposal generators from this graph.
 edges=sorted(es);uniform=[Fraction(1)]*len(edges);biased=[Fraction(2 if i==0 else 1) for i in range(len(edges))]
 Lu=weighted_laplacian(n,edges,uniform);Lb=weighted_laplacian(n,edges,biased)
 if any(sum(row) for row in Lu+Lb):raise AssertionError("generator does not conserve constants")
 support_u={(i,j) for i in range(n) for j in range(n) if i!=j and Lu[i][j]}
 support_b={(i,j) for i in range(n) for j in range(n) if i!=j and Lb[i][j]}
 expected_support={(u,v) for u,v in edges}|{(v,u) for u,v in edges}
 if support_u!=expected_support or support_b!=expected_support:raise AssertionError("weighted generator support differs from seams")
 if Lu==Lb or proportional(Lb,Lu):raise AssertionError("biased law did not produce a nonproportional generator")
 if not all(w>0 for w in uniform+biased):raise AssertionError("positive-law check failed")
 result["uniform_and_biased_laws"]={"both_positive":True,"both_row_sums_zero":True,"same_30_seam_support":True,"distinct":True,"nonproportional":True}
 result["receipt_sha256"]=verify_receipts()
 return result

if __name__=="__main__":print(json.dumps(audit(),sort_keys=True))

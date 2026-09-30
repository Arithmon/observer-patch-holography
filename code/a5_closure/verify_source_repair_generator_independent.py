#!/usr/bin/env python3
"""Independent graph and invariant-space verifier; no producer imports."""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
MANIFEST=ROOT/"code/a5_closure/manifests/echosahedral_federation_reference.json"

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

def audit():
 c=json.loads(MANIFEST.read_text())["carrier"]; lab=c["ports"];ix={x:i for i,x in enumerate(lab)}
 fs=[tuple(ix[x] for x in f) for f in c["oriented_faces"]]
 es={tuple(sorted(z)) for a,b,d in fs for z in ((a,b),(b,d),(d,a))};n=len(lab)
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
 result={"vertices":n,"edges":len(es),"degrees":sorted(set(deg)),"connected":len(seen)==n,"proper_group_order":len(perms),"vertex_orbit_sizes":sorted(map(len,vorb)),"edge_orbit_sizes":sorted(map(len,eorbs)),"invariant_weight_dimension":dim}
 if result!={"vertices":12,"edges":30,"degrees":[5],"connected":True,"proper_group_order":60,"vertex_orbit_sizes":[12],"edge_orbit_sizes":[30],"invariant_weight_dimension":1}:raise AssertionError(result)
 # Invariance equations have nullspace spanned by the all-ones assignment.
 return result

if __name__=="__main__":print(json.dumps(audit(),sort_keys=True))

#!/usr/bin/env python3
"""Exact finite audit of the native twelve-port seam proposal boundary.

The graph is derived from oriented face incidence in the pinned carrier
manifest. The relation layer admits any positive edge-rate vector. A5
invariance is then solved as equality constraints on the 30 edge variables.
"""
from __future__ import annotations
import argparse, hashlib, json
from collections import deque
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "code/a5_closure/manifests/echosahedral_federation_reference.json"
A3_PATH = ROOT / "code/a5_closure/manifests/a3_scheduler_kernel_reference.json"
DIRECTED_PATH = ROOT / "code/a5_closure/manifests/directed_seam_repair_reference.json"
REPAIR_PATH = ROOT / "code/a5_closure/manifests/record_counting_mechanism_reference.json"
PINNED_RECEIPT_DIGESTS = {
    "a3": "sha256:4afffa0760ff59649d1f2548d8f755f1efe4c5be14594bbde80f22d7ed3d536a",
    "directed": "sha256:29d77ed315cba2e98a1f1bece56191edf8759a472e59799200716ee70163e4b8",
    "repair": "sha256:efad022da8fdb52b58e6e8356ef945227832467c162cb512d0261812f286a655",
}

def need(ok: bool, why: str) -> None:
    if not ok: raise ValueError(why)

def data_from_faces(doc: dict):
    c=doc["carrier"]; ports=c["ports"]; ix={p:i for i,p in enumerate(ports)}
    faces=[tuple(ix[p] for p in f) for f in c["oriented_faces"]]
    edges={tuple(sorted(e)) for a,b,c in faces for e in ((a,b),(b,c),(c,a))}
    return ports,faces,edges

def ckey(f):
    a,b,c=f; return min((a,b,c),(b,c,a),(c,a,b))

def proper_autos(n, edges, faces):
    adj=[set() for _ in range(n)]
    for u,v in edges: adj[u].add(v); adj[v].add(u)
    face_set={ckey(f) for f in faces}; result=[]; image=[-1]*n; used=set()
    # Assign vertices in a fixed rare-signature-first order and prune using
    # the already mapped induced subgraph.
    order=sorted(range(n),key=lambda x:(-len(adj[x]),x))
    def visit(k):
        if k==n:
            p=tuple(image)
            if all(ckey(tuple(p[x] for x in f)) in face_set for f in faces): result.append(p)
            return
        u=order[k]
        for v in range(n):
            if v in used or len(adj[u])!=len(adj[v]): continue
            if any(((z in adj[u]) != (image[z] in adj[v])) for z in order[:k]): continue
            image[u]=v; used.add(v); visit(k+1); used.remove(v); image[u]=-1
    visit(0); return result

def rank_q(rows):
    a=[[Fraction(x) for x in r] for r in rows]; m=len(a); n=len(a[0]) if m else 0; rank=0
    for col in range(n):
        piv=next((i for i in range(rank,m) if a[i][col]),None)
        if piv is None: continue
        a[rank],a[piv]=a[piv],a[rank]; q=a[rank][col]; a[rank]=[x/q for x in a[rank]]
        for i in range(m):
            if i!=rank and a[i][col]:
                q=a[i][col]; a[i]=[x-q*y for x,y in zip(a[i],a[rank])]
        rank+=1
    return rank

def classify(doc=None):
    doc=json.loads(MANIFEST.read_text()) if doc is None else doc
    ports,faces,edges=data_from_faces(doc); n=len(ports)
    need(len(faces)==20 and len({ckey(f) for f in faces})==20,"expected twenty unique oriented faces")
    pinned={tuple(sorted((ports.index(a),ports.index(b)))) for a,b in doc["carrier"]["edges"]}
    need(len(doc["carrier"]["edges"])==30 and len(pinned)==30,"expected thirty unique serialized seams")
    need(edges==pinned,"face-derived incidence differs from carrier receipt")
    deg=[sum(i in e for e in edges) for i in range(n)]
    adj=[set() for _ in range(n)]
    for u,v in edges: adj[u].add(v);adj[v].add(u)
    seen={0};todo=[0]
    while todo:
        u=todo.pop()
        for v in adj[u]-seen: seen.add(v);todo.append(v)
    group=proper_autos(n,edges,faces)
    v_orbits={frozenset(g[i] for g in group) for i in range(n)}
    edge_ix={e:i for i,e in enumerate(sorted(edges))}; unseen=set(edges); e_orbits=[]; rows=[]
    for g in group:
        for e,j in edge_ix.items():
            z=tuple(sorted((g[e[0]],g[e[1]])))
            if z != e:
                row=[0]*len(edges);row[j]=1;row[edge_ix[z]]=-1;rows.append(row)
    while unseen:
        e=min(unseen); orb={tuple(sorted((g[e[0]],g[e[1]]))) for g in group}
        need(orb<=edges,"proper action failed to preserve seams")
        e_orbits.append(orb);unseen-=orb
    dim=len(edges)-rank_q(rows)
    actual=(n,len(edges),set(deg),len(seen),len(group),len(v_orbits),sorted(map(len,e_orbits)),dim)
    need(actual==(12,30,{5},12,60,1,[30],1),f"incidence/A5/invariant-rank failure: {actual}")
    return {"n":n,"edges":sorted(edges),"faces":faces,"group":group,"degrees":deg,"edge_orbits":[sorted(x) for x in e_orbits],"invariant_dimension":dim}

def canonical_digest(value):
    raw=json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()
    return "sha256:"+hashlib.sha256(raw).hexdigest()

def read_receipts(receipts=None):
    if receipts is None:
        receipts={"a3":json.loads(A3_PATH.read_text()),"directed":json.loads(DIRECTED_PATH.read_text()),"repair":json.loads(REPAIR_PATH.read_text())}
    a3,directed,repair=(receipts[k] for k in ("a3","directed","repair"))
    expected={"a3":"oph.a3_scheduler_kernel_certificate.v1","directed":"oph.directed_seam_repair_certificate.v1","repair":"oph.record_counting_mechanism_certificate.v2"}
    for key,value in (("a3",a3),("directed",directed),("repair",repair)):
        need(value.get("schema")==expected[key],f"wrong {key} receipt schema")
        body={k:v for k,v in value.items() if k!="manifest_sha256"}
        need(value.get("manifest_sha256")==canonical_digest(body),f"{key} receipt self-digest mismatch")
        need(canonical_digest(value)==PINNED_RECEIPT_DIGESTS[key],f"{key} receipt differs from audited baseline pin")
    need(a3.get("issue")==614 and repair.get("issue")==628,"wrong issue pins")
    pins=directed.get("upstream_pins",{})
    for key,value,path,schema in (
        ("undirected_scheduler",a3,"manifests/a3_scheduler_kernel_reference.json",expected["a3"]),
        ("integer_record_counting_mechanism",repair,"manifests/record_counting_mechanism_reference.json",expected["repair"])):
        pin=pins.get(key,{})
        need(pin.get("path")==path and pin.get("schema")==schema,"directed receipt missing expected upstream pin")
        need(pin.get("canonical_json_sha256")==canonical_digest(value),"directed receipt upstream digest mismatch")
    selection=a3.get("a3_selection",{});deck=selection.get("deck_invariance",{})
    selected=(a3.get("verdict",{}).get("kernel_status")=="A3_selected_exact" and a3.get("move_simplex",{}).get("seams")==30 and selection.get("reference")=="uniform_on_moves" and selection.get("feasible_set")=="the full move simplex; the complete constraint grammar contributes normalization only" and selection.get("selected_kernel_probability_per_seam")=="1/30" and selection.get("positive_weights") is True and selection.get("kernel_sums_to_one") is True and deck.get("edge_transitive") is True and deck.get("seam_orbit_size")==30)
    need(selected,"#614 no longer certifies its uniform A3 move law")
    s1=directed.get("uniform_s1_channel",{});verdict=directed.get("verdict",{});scope=directed.get("scope",{})
    need(s1.get("sector")=="S=1 total integer working load","directed receipt sector changed")
    need(verdict.get("uniform_s1_channel")=="exact_named_realization","directed S=1 status changed")
    need(verdict.get("directed_orbit_schedule_as_physical_selection")=="open","directed physical-selection boundary changed")
    need(verdict.get("full_self_readback_and_universe_selection")=="open","full self-readback boundary changed")
    claim=scope.get("claim_boundary","")
    need("does not prove an IID" in claim and "path law" in claim,"directed scope boundary changed")
    # The #628 receipt types moves and settling but does not carry a proposal
    # distribution; #614's uniform reference is an external input to it.
    mech=repair.get("mechanism",{})
    relation_only=isinstance(mech,dict) and "proposal_law" not in mech and "proposal_weights" not in mech
    need(relation_only,"#628 mechanism unexpectedly acquired proposal fields; reassess verdict")
    hashes={}
    for k,path,value in (("a3",A3_PATH,a3),("directed",DIRECTED_PATH,directed),("repair",REPAIR_PATH,repair)):
        hashes[k]={"path":str(path.relative_to(ROOT)),"raw_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"canonical_json_sha256":canonical_digest(value),"schema":value["schema"]}
    return {"issue_614_selects_uniform_on_declared_uniform_move_reference":True,
      "reference_source_derived_from_issue_628_relation":False,
      "issue_614_law_domain":"thirty-seam scheduler on named reference carrier",
      "issue_628_relation_alone_selects_weights":False,
      "directed_receipt_channel_domain":s1["sector"],
      "full_integer_state_generator_binding":"not established",
      "receipts":hashes}

def laplacian(n, edges, weights):
    L=[[Fraction() for _ in range(n)] for _ in range(n)]
    for (u,v),w in zip(edges,weights):
        need(w>=0,"negative rate");L[u][u]+=w;L[v][v]+=w;L[u][v]-=w;L[v][u]-=w
    return L

def proportional(a,b):
    ratio=None
    for i in range(len(a)):
        for j in range(len(a)):
            if b[i][j]:
                q=a[i][j]/b[i][j]
                if ratio is None: ratio=q
                elif q!=ratio: return False
            elif a[i][j]: return False
    return ratio is not None

def certificate(receipts=None):
    g=classify(); edges=g["edges"]
    wu=[Fraction(1)]*30; wb=[Fraction(2 if i==0 else 1) for i in range(30)]
    Lu=laplacian(12,edges,wu); Lb=laplacian(12,edges,wb)
    for L in (Lu,Lb): need(all(sum(r)==0 for r in L),"conservation failure")
    support={(i,j) for i in range(12) for j in range(12) if i!=j and Lu[i][j]}
    need(support=={(u,v) for u,v in edges}|{(v,u) for u,v in edges},"support failure")
    need(Lu!=Lb,"positive biased countermodel collapsed to uniform")
    need(not proportional(Lb,Lu),"biased generator is unexpectedly proportional to uniform")
    need(any(Lb[u][v]!=Lu[u][v] for u,v in edges),"no edge weight distinction")
    return {"graph":{"vertices":12,"edges":30,"degree":5,"connected":True,"proper_A5_order":60,"vertex_orbits":[12],"unoriented_edge_orbits":[30]},
      "invariant_edge_weight_dimension":g["invariant_dimension"],"invariant_basis":"all-ones vector on the 30-edge orbit",
      "general_symmetric_conservative_local_generator":"(L_w f)(i)=sum_{j~i} w_ij(f(i)-f(j)); w_ij>=0; L=D_w-W",
      "uniform_law":"w_e=1; L=5I-A",
      "biased_law":"w_e=2 on lexicographically first face-derived seam, 1 on the other 29 seams",
      "countermodel":{"same_relation":True,"strictly_positive":True,"both_conservative":True,"same_30_seam_support":True,"distinct_generators":True,"nonproportional_generators":True,"biased_not_A5_invariant":True},
      "proposal_source_audit":read_receipts(receipts),
      "spectrum_of_uniform_generator":{"P3":"5-sqrt(5)","P5":"6","P3prime":"5+sqrt(5)","strict_order":"5-sqrt(5)<6<5+sqrt(5), since sqrt(5)>1","basis":"conditional consequence; the exact band identities are already certified in PortGramRepairBand.lean; generic epsilon range not replayed here"},
      "galois_consequence":"under the declared linear repair mean, P3 is the unique slow nonconstant band and its frozen Gram is G_plus=4P3; this packet does not bind that mean to the full nonlinear integer repair chain"}

def main():
    p=argparse.ArgumentParser();p.add_argument("--json",action="store_true");a=p.parse_args();out=certificate()
    if a.json: print(json.dumps(out,indent=2,sort_keys=True))
    else: print("A5_INVARIANT_REPAIR_LAW_FORCES_LOW_BAND (conditional)\nREPAIR_RELATION_DOES_NOT_SELECT_PROPOSAL_LAW\nlocal invariant edge-weight dimension=1; positive biased relation-preserving countermodel exists")
if __name__=="__main__": main()

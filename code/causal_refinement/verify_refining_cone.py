"""Independent integer path checks and complete threaded-provenance replay."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUTPUT=HERE/"refining_cone_receipt.json"

def need(ok,message):
    if not ok:
        raise ValueError(message)

def same(a,b,label):
    need(type(a) is type(b),label+" type")
    if isinstance(b,dict):
        need(set(a)==set(b),label+" fields")
        for key in b:
            same(a[key],b[key],label+"."+key)
    elif isinstance(b,list):
        need(len(a)==len(b),label+" length")
        for x,y in zip(a,b):
            same(x,y,label+" entry")
    else:
        need(a==b,label+" mismatch")

def load(path=OUTPUT):
    def pairs(items):
        out={}
        for k,v in items:
            need(k not in out,"duplicate key")
            out[k]=v
        return out
    def invalid(s):
        raise ValueError("noninteger JSON number: "+s)
    return json.loads(Path(path).read_text(encoding="utf-8"),object_pairs_hook=pairs,
                      parse_float=invalid,parse_constant=invalid)

def vector(value,size):
    need(type(value) is list and len(value)==size and all(type(x) is int for x in value),"integer vector")
    return tuple(value)

def sq(v):
    return sum(x*x for x in v)

def replay(rows):
    need(type(rows) is list and len(rows)==81,"complete finite history")
    sites=set(product(range(-1,2),repeat=3))
    states={}
    edges=set()
    for row in rows:
        need(type(row) is dict and set(row)=={"event","reads","write_value"},"commit fields")
        event=vector(row["event"],4)
        layer,*site=event
        need(layer in range(3) and tuple(site) in sites and event not in states,"fresh event")
        required=set() if layer==0 else {(layer-1,*old) for old in sites
             if sq([x-y for x,y in zip(old,site)])<=4}
        need(type(row["reads"]) is list,"read list")
        consumed=set()
        total=1
        for read in row["reads"]:
            need(type(read) is dict and set(read)=={"register","writer","value"},"read fields")
            register=vector(read["register"],4)
            need(register in required and register not in consumed,"complete distinct causal support")
            need(register in states,"read from committed prefix")
            same(read["writer"],list(register),"actual last writer")
            same(read["value"],states[register],"actual written value")
            consumed.add(register)
            total+=read["value"]
            edges.add((register,event))
        need(consumed==required,"all move-law parents consumed")
        same(row["write_value"],total,"state depends on consumed records")
        states[event]=total
    need(set(states)=={(j,*x) for j in range(3) for x in sites},"carrier coverage")
    # Same-layer antichain (all edges increase time), and 27 vertical chains
    # cover the carrier because zero displacement is allowed: width exactly27.
    need(all(((j,*x),(j+1,*x)) in edges for j in range(2) for x in sites),"vertical chain cover")
    return states,edges

def verify(packet):
    need(type(packet) is dict and set(packet)=={"schema","scope","cone_paths","history","reordered_history","fixed_axial_control"},"packet fields")
    same(packet["schema"],"oph.refining-causal-cone.v1","schema")
    same(packet["scope"],{"spatial_dimension_supplied":3,"growing_move_menu_supplied":True,
         "physical_clock_selected":False,"native_OPH_manifold":False,
         "poisson_sprinkling":False,"quantum_field_limit":False},"scope")
    need(type(packet["cone_paths"]) is list and len(packet["cone_paths"])==16,"probe census")
    maximum_step_squared=0
    for index,row in enumerate(packet["cone_paths"]):
        need(type(row) is dict and set(row)=={"radius","steps","h","endpoint","path"},"path fields")
        r=(4,8,16,32)[index//4]
        same(row["radius"],r,"radius")
        same(row["steps"],r,"steps")
        same(row["h"],str(F(1,r*r)),"spatial scale")
        reach=r*(r-3)
        expected=((reach,0,0),(reach//2,reach//2,reach//2),(-reach//2,reach//3,reach//4),(0,0,0))[index%4]
        same(row["endpoint"],list(expected),"endpoint")
        d=vector(row["endpoint"],3)
        need(sq(d)<=reach**2,"inner-cone premise")
        need(type(row["path"]) is list and len(row["path"])==r+1,"path length")
        points=[vector(x,3) for x in row["path"]]
        need(points[0]==(0,0,0) and points[-1]==d,"path endpoints")
        for j,x in enumerate(points):
            # Integer inequalities characterize floor without calling proposer.
            need(all(r*a<=j*b<r*(a+1) for a,b in zip(x,d)),"rounded path")
        for a,b in zip(points,points[1:]):
            displacement=tuple(y-x for x,y in zip(a,b))
            need(sq(displacement)<=r*r,"local move ball")
            maximum_step_squared=max(maximum_step_squared,sq(displacement))
    first=replay(packet["history"])
    second=replay(packet["reordered_history"])
    need(first==second,"independent same-layer reorder invariance")
    same(packet["fixed_axial_control"],{"steps":3,"endpoint":[2,2,0],
         "euclidean_norm_squared":8,"minimum_axial_moves":4},"axial control")
    # Independent exhaustive finite reachability excludes a Euclidean-inside pair.
    moves=[(0,0,0)]+[tuple(s*int(i==a) for i in range(3)) for a in range(3) for s in (-1,1)]
    reachable={(0,0,0)}
    for _ in range(3):
        reachable={tuple(a+b for a,b in zip(x,v)) for x in reachable for v in moves}
    need((2,2,0) not in reachable and 8<9,"fixed-menu cone discrepancy")
    return {"verdict":"PASS_SUPPLIED_LAW_CAUSAL_REFINEMENT_CONTROL",
         "path_witnesses":16,"largest_radius":32,"largest_inner_speed":"29/32",
         "history_events":81,"authenticated_edges":len(first[1]),"exact_history_width":27,
         "native_physical_spacetime_selected":False,
         "continuum_claim":"analytic inner/outer bound; supplied grid, clock and density"}

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path",nargs="?",type=Path,default=OUTPUT)
    args=parser.parse_args()
    print(json.dumps(verify(load(args.path)),indent=2,sort_keys=True))

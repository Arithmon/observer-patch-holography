"""Produce finite controls for a supplied growing-menu causal lattice."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent

def paths():
    rows=[]
    # Source menu and endpoint probes are explicit, with no fitted statistic.
    for radius in (4,8,16,32):
        steps=radius
        reach=steps*(radius-3)
        endpoints=((reach,0,0),(reach//2,reach//2,reach//2),
                   (-reach//2,reach//3,reach//4),(0,0,0))
        for endpoint in endpoints:
            rows.append({"radius":radius,"steps":steps,"h":str(F(1,radius**2)),
                "endpoint":list(endpoint),
                "path":[[j*x//steps for x in endpoint] for j in range(steps+1)]})
    return rows

def history(reverse=False):
    spatial=list(product(range(-1,2),repeat=3))
    state={}
    rows=[]
    radius=2
    for layer in range(3):
        for site in (list(reversed(spatial)) if reverse else spatial):
            event=(layer,*site)
            parents=[] if layer==0 else [(layer-1,*old) for old in spatial
                if sum((x-y)**2 for x,y in zip(site,old))<=radius**2]
            parents.sort()
            reads=[{"register":list(p),"writer":list(p),"value":state[p]} for p in parents]
            value=1+sum(x["value"] for x in reads)
            rows.append({"event":list(event),"reads":reads,"write_value":value})
            state[event]=value
    return rows

def produce():
    return {"schema":"oph.refining-causal-cone.v1",
        "scope":{"spatial_dimension_supplied":3,"growing_move_menu_supplied":True,
                 "physical_clock_selected":False,"native_OPH_manifold":False,
                 "poisson_sprinkling":False,"quantum_field_limit":False},
        "cone_paths":paths(),"history":history(),"reordered_history":history(True),
        "fixed_axial_control":{"steps":3,"endpoint":[2,2,0],
                               "euclidean_norm_squared":8,"minimum_axial_moves":4}}

if __name__=="__main__":
    (HERE/"refining_cone_receipt.json").write_text(json.dumps(produce(),sort_keys=True,separators=(",",":"))+"\n",encoding="utf-8")

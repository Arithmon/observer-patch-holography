"""Produce compact candidate evidence; verification is a separate entry point."""

import json
from itertools import groupby
from pathlib import Path

from . import model, balanced, balanced_check
from .audit import consumer_audit, downstream_audit
from .check import digest
from .verify import pins, route_targets

HERE = Path(__file__).resolve().parent


def graph_summary(n, family, intervention=None):
    sites, parents, states = model.graph(n,family,intervention=intervention)
    return {"states":[digest([state[x] for x in sites]) for state in states],
            "events":len(sites)*len(states),
            "reads":3*sum(map(len,parents.values()))}


def build():
    cuts = {}
    for n in (1,2,3,4,6,8,12,16,20):
        cuts[str(n)] = {"stencil_sha256":digest(model.sparse(n)), "cuts":model.cuts(n),
                        "actions":{family:model.action(n,family) for family in
                                   ("sparse","axial")+(("dense",) if n<=4 else ())}}
    routes = []
    for n,d in route_targets():
        steps = model.route(n,d)
        instructions = [[len(list(group)),*v] for v,group in groupby(steps)]
        place,points = [0,0,0],[[0,0,0]]
        for v in steps:
            place = [x+y for x,y in zip(place,v)]
            points.append(place)
        routes.append({"n":n,"displacement":list(d),"instructions":instructions,
                       "steps":len(steps),"points_sha256":digest(points)})
    graphs = {}
    for n in (1,2):
        _,q = model.parameters(n)
        for family in ("sparse","axial","dense"):
            graphs[f"{n}:{family}"] = {
                "baseline":graph_summary(n,family),
                "initial":graph_summary(n,family,(0,[q//2]*3)),
                "intermediate":graph_summary(n,family,(1,[q//2]*3))}
    balanced_routes=[]
    for t,d in balanced_check.route_targets():
        steps=balanced.route(t,d)
        point,points=[0,0,0],[[0,0,0]]
        for v in steps:
            point=[x+y for x,y in zip(point,v)];points.append(point)
        balanced_routes.append(dict(t=t,displacement=list(d),steps=len(steps),
                                    points_sha256=digest(points),
                                    instructions=[[len(list(g)),*v] for v,g in groupby(steps)]))
    return {"schema":"oph-m1-necessity-v1","sources":pins(),"consumers":consumer_audit(),
            "cuts":cuts,"routes":routes,"graphs":graphs,"downstream":downstream_audit(),
            "balanced":{"cuts":{str(t):balanced.certificate(t) for t in (1,2,3)},
                        "routes":balanced_routes}}


def render(packet):
    # One line per route at either level; no expanded event or path tapes.
    def encode(obj,indent=0):
        if not isinstance(obj,dict):
            return json.dumps(obj,indent=2,sort_keys=True)
        parts=[]
        for key,value in sorted(obj.items()):
            if key=='routes':
                encoded='[\n'+',\n'.join(' '*(indent+4)+json.dumps(row,sort_keys=True) for row in value)+'\n'+' '*(indent+2)+']'
            elif key=='balanced':
                encoded=encode(value,indent+2)
            else:
                encoded=json.dumps(value,indent=2,sort_keys=True).replace('\n','\n'+' '*(indent+2))
            parts.append(' '*(indent+2)+json.dumps(key)+': '+encoded)
        return '{\n'+',\n'.join(parts)+'\n'+' '*indent+'}'
    return encode(packet)+'\n'


if __name__ == "__main__":
    result = build()
    (HERE/"receipt.json").write_text(render(result),encoding="utf-8",newline="\n")
    print(f"built {len(result['routes'])} routes, {len(result['balanced']['routes'])} balanced routes, "
          f"{len(result['cuts'])} original cut levels and {len(result['graphs'])} full graphs")

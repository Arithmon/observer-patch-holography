"""Untrusted deterministic paired-route producer; witnesses require check_routes."""
from collections import deque
from heapq import heappop, heappush
import json
from pathlib import Path
import sys
import time

from . import codec

ROOT = codec.ROOT


def graph(level):
    path = ("code/source_routing/support_w12_l3.json" if level == 3 else
            f"code/source_read_routing/support_w12_l{level}.json")
    data = json.loads((ROOT/path).read_text(encoding="utf-8"))
    adjacency = [set() for _ in range(12*data["carriers"])]
    for carrier in range(data["carriers"]):
        for a,b in data["intra_carrier_seams"]:
            u,v = 12*carrier+a,12*carrier+b
            adjacency[u].add(v)
            adjacency[v].add(u)
    for c,a,d,b in data["glued_pairs"]:
        u,v = 12*c+a,12*d+b
        adjacency[u].add(v)
        adjacency[v].add(u)
    return [tuple(sorted(row)) for row in adjacency]


def distances(adjacency, target, forbidden):
    result = [-1]*len(adjacency)
    result[target] = 0
    queue = deque([target])
    while queue:
        u = queue.popleft()
        for v in adjacency[u]:
            if result[v]<0 and v not in forbidden:
                result[v] = result[u]+1
                queue.append(v)
    return result


def shortest(adjacency, source, target, forbidden, order=0, estimates=None):
    estimate = estimates[target] if estimates is not None else None
    def lower(v):
        if estimate is None:
            return 0
        if estimate[v]>=0:
            return estimate[v]
        available = [estimate[u]+1 for u in adjacency[v] if estimate[u]>=0]
        return min(available,default=len(adjacency))
    parents = {source:None}
    best = {source:0}
    queue = [(lower(source),0,source)]
    while queue:
        _,cost,u = heappop(queue)
        if best[u] != cost:
            continue
        if u == target:
            path = [u]
            while parents[u] is not None:
                u = parents[u]
                path.append(u)
            return path[::-1]
        neighbours = adjacency[u]
        if order:
            if order%2:
                neighbours = neighbours[::-1]
            offset = (order//2)%len(neighbours)
            neighbours = neighbours[offset:]+neighbours[:offset]
        for v in neighbours:
            if v not in forbidden and cost+1<best.get(v,len(adjacency)):
                parents[v] = u
                best[v] = cost+1
                heappush(queue,(cost+1+lower(v),cost+1,v))
    raise ValueError("no disjoint route")


def paired(adjacency, source, target, forbidden, estimates=None):
    for order,first,flip in ((order,first,flip) for order in range(12)
                             for first in (0,1) for flip in (False,True)):
        end = target[::-1] if flip else target
        second = 1-first
        try:
            p = shortest(adjacency,source[first],end[first],
                         forbidden|{source[second],end[second]},order,estimates)
            q = shortest(adjacency,source[second],end[second],forbidden|set(p),order,estimates)
        except ValueError:
            continue
        paths = [None,None]
        paths[first],paths[second] = p,q
        def grow(index, extra):
            used = forbidden|set(paths[0])|set(paths[1])
            path = paths[index]
            for at,(u,v) in enumerate(zip(path,path[1:])):
                if extra == 1:
                    choices = sorted((set(adjacency[u])&set(adjacency[v]))-used)
                    if choices:
                        path.insert(at+1,choices[0])
                        return True
                else:
                    for a in adjacency[u]:
                        if a in used:
                            continue
                        choices = sorted((set(adjacency[a])&set(adjacency[v]))-used-{a})
                        if choices:
                            path[at+1:at+1] = [a,choices[0]]
                            return True
            return False

        for _ in range(64):
            if len(paths[0]) == len(paths[1]):
                break
            short = 0 if len(paths[0])<len(paths[1]) else 1
            difference = abs(len(paths[0])-len(paths[1]))
            if grow(short,1):
                continue
            if difference>=2 and grow(short,2):
                continue
            if difference==1 and grow(1-short,1) and grow(short,2):
                continue
            else:
                break
        if len(paths[0]) == len(paths[1]):
            assert len(set(paths[0]+paths[1])) == 2*len(paths[0])
            return [list(p) for p in zip(*paths)]
    raise ValueError("cannot equalize disjoint paths")


def generate(n,level):
    adjacency = graph(level)
    protected = {12*c+r for c in range(1,2*n+1) for r in (0,1)}
    read_blocked = protected|{2,1,4,9,6,7}
    write_blocked = protected|{2,1,6,7,3,5}
    read_estimates = {p:distances(adjacency,p,read_blocked) for p in (3,5)}
    write_estimates = {p:distances(adjacency,p,write_blocked) for p in (4,9)}
    if not 1 <= n or 2*n+1 > len(adjacency)//12:
        raise ValueError("two-bank capacity")
    for carrier in range(1,2*n+1):
        record = (12*carrier,12*carrier+1)
        scratch = (12*carrier+5,12*carrier+9)
        for direction in ("read","write"):
            if direction == "read":
                source,target = scratch,(3,5)
                forbidden = read_blocked
                estimates = read_estimates
            else:
                source,target = record,(4,9)
                forbidden = write_blocked-set(record)
                estimates = write_estimates
            route = paired(adjacency,source,target,forbidden,estimates)
            if direction == "write":
                route.reverse()
                if route[0] != [4,9]:
                    route = [pair[::-1] for pair in route]
            yield {"carrier":carrier,"direction":direction,"pairs":route}


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites",type=int,required=True)
    parser.add_argument("--level",type=int,required=True,choices=(3,4,5))
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    with args.output.open("wb") as stream:
        for row in generate(args.sites,args.level):
            stream.write((codec.compact(row)+"\n").encode("ascii"))


if __name__ == "__main__":
    main()

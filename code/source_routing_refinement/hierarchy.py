"""Exact shared-vertex incidence; no floating-point geometry or fitted rates."""
from __future__ import annotations

from collections import defaultdict, deque

BASE = ((0,11,5), (0,5,1), (0,1,7), (0,7,10), (0,10,11),
        (1,5,9), (5,11,4), (11,10,2), (10,7,6), (7,1,8),
        (3,9,4), (3,4,2), (3,2,6), (3,6,8), (3,8,9),
        (4,9,5), (2,4,11), (6,2,10), (8,6,7), (9,8,1))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def refine(faces):
    next_vertex = 1 + max(max(face) for face in faces)
    midpoints = {}
    result = []

    def midpoint(a, b):
        key = tuple(sorted((a, b)))
        if key not in midpoints:
            midpoints[key] = next_vertex + len(midpoints)
        return midpoints[key]

    for a, b, c in faces:
        ab, bc, ca = midpoint(a,b), midpoint(b,c), midpoint(c,a)
        result.extend(((a,ab,ca), (b,bc,ab), (c,ca,bc), (ab,bc,ca)))
    return tuple(result)


def graph(faces):
    incidence = defaultdict(list)
    for i, face in enumerate(faces):
        require(len(set(face)) == 3, "degenerate triangle")
        for v in face:
            incidence[v].append(i)
    neighbours = [set() for _ in faces]
    for cells in incidence.values():
        for i in cells:
            neighbours[i].update(j for j in cells if i != j)
    return tuple(tuple(sorted(row)) for row in neighbours)


def paths(adjacency, root):
    parents = {root: root}
    queue = deque([root])
    while queue:
        u = queue.popleft()
        for v in adjacency[u]:
            if v not in parents:
                parents[v] = u
                queue.append(v)
    require(len(parents) == len(adjacency), "disconnected graph")
    result = []
    for target in range(len(adjacency)):
        path = [target]
        while path[-1] != root:
            path.append(parents[path[-1]])
        result.append(tuple(reversed(path)))
    return result


def refinement_certificate(coarse, fine):
    """Check the two combinatorial premises of the kernel path-lifting theorem."""
    require(tuple(fine) == refine(coarse), "wrong subdivision stencil")
    require(len(fine) == 4*len(coarse), "child count")
    for parent in range(len(coarse)):
        children = [set(face) for face in fine[4*parent:4*parent+4]]
        require(all(a & b for a in children for b in children), "disconnected child fibre")
        for v in coarse[parent]:
            require(any(v in child for child in children), "lost old vertex")
    # Each coarse edge shares an old vertex. Both retained copies are fine
    # vertices with the same integer ID, so they give the required bridge.
    adjacency = graph(coarse)
    bridges = 0
    for a, row in enumerate(adjacency):
        for b in row:
            if a >= b:
                continue
            shared = set(coarse[a]) & set(coarse[b])
            require(bool(shared), "coarse edge without shared vertex")
            v = min(shared)
            u = next(i for i in range(4*a,4*a+4) if v in fine[i])
            w = next(i for i in range(4*b,4*b+4) if v in fine[i])
            require(u//4 == a and w//4 == b and set(fine[u]) & set(fine[w]), "missing lifted edge")
            bridges += 1
    return {"parents": len(coarse), "children": len(fine), "lifted_coarse_edges": bridges}


def envelope(level, base_diameter=3):
    return (base_diameter+1)*2**level - 1


def minimal_level(sites):
    require(sites > 0, "empty population")
    level = 0
    while 20*4**level < sites:
        level += 1
    return level

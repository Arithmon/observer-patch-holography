"""Exact combinatorial support tower; its mesh bound is proved analytically."""
from itertools import combinations


def seed():
    # Same unlabeled boundary, in a convenient orientation chart.
    return [(0, 1, 2), (0, 3, 1), (0, 2, 4), (0, 6, 3), (0, 4, 6),
            (1, 5, 2), (1, 3, 7), (1, 7, 5), (2, 8, 4), (2, 5, 8),
            (3, 6, 9), (3, 9, 7), (4, 10, 6), (4, 8, 10), (5, 7, 11),
            (5, 11, 8), (6, 10, 9), (7, 9, 11), (8, 11, 10), (9, 10, 11)]


def subdivide(vertices, faces):
    edges = sorted({tuple(sorted(e)) for f in faces for e in combinations(f, 2)})
    midpoint = {e: vertices+i for i, e in enumerate(edges)}
    children = []
    for a, b, c in faces:
        ab = midpoint[tuple(sorted((a, b)))]
        bc = midpoint[tuple(sorted((b, c)))]
        ca = midpoint[tuple(sorted((c, a)))]
        children.extend([(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)])
    # This is a declared coarsening, transported with the presentation. It is
    # not claimed to be a canonical choice from an unlabeled edge.
    coarsen = list(range(vertices)) + [max(e) for e in edges]
    return vertices+len(edges), children, coarsen


def tower(levels=4):
    count, faces = 12, seed()
    out = [{"level": 0, "vertices": count, "faces": faces, "coarsen": None}]
    for level in range(1, levels):
        count, faces, coarsen = subdivide(count, faces)
        out.append({"level": level, "vertices": count, "faces": faces, "coarsen": coarsen})
    return out


def oriented_key(face):
    a, b, c = face
    cyclic = ((a, b, c), (b, c, a), (c, a, b))
    reverse = ((a, c, b), (c, b, a), (b, a, c))
    smallest = min(cyclic+reverse)
    return smallest, 1 if smallest in cyclic else -1


def check_tower(rows):
    """Check chain maps directly, not against the subdivision routine."""
    from collections import Counter
    require = lambda c, m: None if c else (_ for _ in ()).throw(ValueError(m))
    summaries = []
    previous = None
    for level, row in enumerate(rows):
        require(type(row) is dict and set(row) == {"level", "vertices", "faces", "coarsen"},
                "tower row schema")
        count, faces = row["vertices"], row["faces"]
        require(type(count) is int and count == 10*4**level+2, "population")
        require(type(row["level"]) is int and row["level"] == level, "level")
        require(type(faces) is list and len(faces) == 20*4**level, "face count")
        directed = Counter()
        chains = Counter()
        for face in faces:
            require(len(face) == 3 and len(set(face)) == 3 and
                    all(type(v) is int and 0 <= v < count for v in face), "face")
            key, sign = oriented_key(face)
            require(key not in chains, "duplicate face")
            chains[key] = sign
            a, b, c = face
            directed.update(((a, b), (b, c), (c, a)))
        require(all(n == 1 and directed[(b, a)] == 1 for (a, b), n in directed.items()),
                "oriented two-face edges")
        edges = {tuple(sorted(e)) for e in directed}
        require(len(edges) == 30*4**level and count-len(edges)+len(faces) == 2, "Euler census")
        # Connected vertex links establish a closed triangulated surface, not
        # merely an Euler characteristic equal to two.
        for vertex in range(count):
            link_edges = [set(f)-{vertex} for f in faces if vertex in f]
            link = set().union(*link_edges)
            require(all(sum(v in e for e in link_edges) == 2 for v in link), "link degrees")
            visited = {next(iter(link))}
            while True:
                following = visited | set().union(*(e for e in link_edges if e & visited))
                if following == visited:
                    break
                visited = following
            require(visited == link, "link connectivity")
        visited = {0}
        while True:
            following = visited | {b for a, b in directed if a in visited}
            if following == visited:
                break
            visited = following
        require(len(visited) == count, "support connectedness")
        if previous is not None:
            coarse = row["coarsen"]
            require(type(coarse) is list and len(coarse) == count and
                    all(type(v) is int and 0 <= v < previous["vertices"] for v in coarse),
                    "coarsening type")
            require(coarse[:previous["vertices"]] == list(range(previous["vertices"])),
                    "old vertex lineages")
            target = Counter()
            coarse_simplices = {frozenset(s) for f in previous["faces"]
                                for k in (1, 2, 3) for s in combinations(f, k)}
            for f in faces:
                image = tuple(coarse[v] for v in f)
                require(frozenset(image) in coarse_simplices, "simplicial map")
                if len(set(image)) == 3:
                    k, s = oriented_key(image)
                    target[k] += s
            expected = Counter(dict(oriented_key(f) for f in previous["faces"]))
            require(dict(target) == dict(expected), "degree-one fundamental chain")
            fibers = Counter(coarse)
            require(set(fibers) == set(range(previous["vertices"])), "nonempty refinement fibers")
            members = {i: {j for j, parent in enumerate(coarse) if parent == i}
                       for i in fibers}
            neighbors = {i: {b for a, b in directed if a == i} for i in range(count)}
            for parent, children in members.items():
                require(children-{parent} <= neighbors[parent], "connected star fibers")
            coarse_edges = {tuple(sorted(e)) for f in previous["faces"] for e in combinations(f, 2)}
            routes = 0
            for a, b in coarse_edges:
                middle = neighbors[a] & neighbors[b]
                require(len(middle) == 1, "coarse edge lift")
                m = next(iter(middle))
                for start in members[a]:
                    for end in members[b]:
                        path = [start, a, m, b, end]
                        require(all(x == y or (x, y) in directed for x, y in zip(path, path[1:])),
                                "tagged seam routing")
                        routes += 1
            fiber_sizes = sorted(Counter(fibers.values()).items())
        else:
            require(row["coarsen"] is None, "seed coarsening")
            fiber_sizes = []
            routes = 0
        summaries.append({"level": level, "carriers": count, "seams": len(edges),
                          "triple_overlaps": len(faces), "fiber_size_census": fiber_sizes,
                          "tagged_scalar_seam_routes": routes})
        previous = row
    require(len(rows) == 4, "declared tower coverage")
    return summaries


if __name__ == "__main__":
    import json
    print(json.dumps(check_tower(tower()), sort_keys=True))

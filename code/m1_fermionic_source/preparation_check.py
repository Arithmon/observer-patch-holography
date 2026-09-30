"""Independent local-cycle enumeration and binary Gaussian elimination."""

from collections import Counter
import itertools

from .check import exact, need
from .geometry_check import digest, geometry, loop_word


def certificate(q):
    edges, _ = geometry(q)
    size = 32*q**3
    index = {e: k for k, e in enumerate(edges)}
    def address(x):
        return 32*((x[0]*q+x[1])*q+x[2])
    def coords(v):
        c = v//32
        return (c//(q*q), (c//q) % q, c % q)
    def root_path(v):
        base, local = divmod(v, 32)
        path = [v]
        for bit in (16, 8, 4):
            if local & bit:
                local -= bit
                path.append(32*base+local)
        while local:
            local -= 1
            path.append(32*base+local)
        return path
    internal = {tuple(sorted((v, root_path(v)[1]))) for v in range(size) if v % 32}
    grid = {e for e in edges if e[0] % 32 == e[1] % 32 == 0}
    tree = set(internal)
    for e in grid:
        x, y = coords(e[0]), coords(e[1])
        axis = next(k for k in range(3) if x[k] != y[k])
        if axis == 2 or (axis == 1 and x[2] == 0) or (axis == 0 and x[1:] == (0, 0)):
            tree.add(e)
    need(len(tree) == size-1, 'axial spanning tree size')
    cycles = []
    for axes, starts in [((0, 1), itertools.product(range(q-1), range(q-1), [0])),
                         ((0, 2), itertools.product(range(q-1), range(q), range(q-1))),
                         ((1, 2), itertools.product(range(q), range(q-1), range(q-1)))]:
        for start in starts:
            path, current = [address(start)], list(start)
            for axis, delta in [(axes[0], 1), (axes[1], 1), (axes[0], -1), (axes[1], -1)]:
                current[axis] += delta
                path.append(address(current))
            cycles.append(path)
    for e in edges:
        if e in internal or e in grid:
            continue
        left, right = root_path(e[0]), root_path(e[1])
        x, y = list(coords(e[0])), coords(e[1])
        if x == list(y):
            meeting = next(v for v in left if v in right)
            path = left[:left.index(meeting)+1]+right[:right.index(meeting)][::-1]
        else:
            path = left.copy()
            for axis in range(3):
                while x[axis] != y[axis]:
                    x[axis] += 1 if y[axis] > x[axis] else -1
                    path.append(address(x))
            path += right[:-1][::-1]
        cycles.append(path+[e[0]])
    rank = len(edges)-size+1
    need(len(cycles) == rank, 'complete independent local loop basis')
    incident = [[] for _ in range(size)]
    for k, (i, j) in enumerate(edges):
        incident[i].append(k)
        incident[j].append(k)
    words, weights, anchors = [], [], []
    for path in cycles:
        need(path[0] == path[-1] and len(set(path[:-1])) == len(path)-1, 'simple local loop')
        raw, weight = loop_word(path, edges, incident)
        words.append(raw)
        weights.append(weight)
        anchor = tuple(min(coords(v)[k] for v in path) for k in range(3))
        anchors.append(anchor)
        support = int(raw[0]) | int(raw[1])
        for e, (i, j) in enumerate(edges):
            if support >> e & 1:
                midpoint = [(a+b)/2 for a, b in zip(coords(i), coords(j))]
                need(all(-.5 <= t-a <= 1.5 for t, a in zip(midpoint, anchor)), 'local measurement owner box')
    free = [k for k, e in enumerate(edges) if e not in tree]
    # Invert the cycle/edge incidence matrix with all axial-tree variables
    # fixed to zero. This uses none of the producer's prefix recurrences.
    rows = []
    for k, raw in enumerate(words):
        x = int(raw[0])
        lhs = sum(((x >> e) & 1) << j for j, e in enumerate(free))
        rows.append([lhs, 1 << k])
    for col in range(rank):
        pivot = next((j for j in range(col, rank) if rows[j][0] >> col & 1), None)
        need(pivot is not None, 'local loop independence')
        rows[col], rows[pivot] = rows[pivot], rows[col]
        for j in range(rank):
            if j != col and rows[j][0] >> col & 1:
                rows[j][0] ^= rows[col][0]
                rows[j][1] ^= rows[col][1]
    decoder = [0]*len(edges)
    for j, e in enumerate(free):
        need(rows[j][0] == 1 << j, 'binary inverse')
        decoder[e] = rows[j][1]
    occupancy = Counter(anchors)
    need(max(occupancy.values()) <= 48 and max(weights) <= 368 and max(map(len, cycles)) <= 17,
         'uniform local preparation bounds')
    return dict(side=q, cycle_count=rank, maximum_length=max(map(len, cycles))-1,
                maximum_weight=max(weights), maximum_slots_per_anchor=max(occupancy.values()),
                cycle_paths_sha256=digest(cycles), signed_loops_sha256=digest(words),
                decoder_sha256=digest([str(v) for v in decoder]), prefix_rounds=2*(q-1),
                local_rounds=4, colors=64, slots_upper=48, loop_length_upper=16, loop_weight_upper=368)


def verify_preparation(rows):
    need(type(rows) is list and len(rows) == 4, 'all local preparation cutoffs')
    for row, q in zip(rows, (1, 2, 3, 4)):
        exact(row, certificate(q))

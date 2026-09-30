"""Local loop basis and a linear, causal decoder for every preparation record."""

from collections import defaultdict
import itertools

from .model import Encoding
from .spatial import digest, graph


def local_parent(vertex):
    cell, local = divmod(vertex, 32)
    if local >= 16:
        return 32*cell+local-16
    if local >= 8:
        return vertex-8
    if local >= 4:
        return vertex-4
    return vertex-1 if local else None


def plan(q):
    vertices, edges, _ = graph(q)
    edge_id = {e: k for k, e in enumerate(edges)}
    ids = {v: i for i, v in enumerate(vertices)}
    def hub(x):
        return ids[(tuple(x), 0, 0)]
    def index(i, j):
        return edge_id[tuple(sorted((i, j)))]
    def grid_edge(x, axis):
        y = list(x)
        y[axis] += 1
        return index(hub(x), hub(y))
    internal = {index(v, local_parent(v)) for v in range(len(vertices)) if local_parent(v) is not None}
    skeleton = internal | {k for k, (i, j) in enumerate(edges) if i % 32 == j % 32 == 0}
    cycles, labels = [], []
    def plaquette(x, a, b, label):
        y, z, w = list(x), list(x), list(x)
        y[a] += 1
        z[a] += 1
        z[b] += 1
        w[b] += 1
        cycles.append([hub(x), hub(y), hub(z), hub(w), hub(x)])
        labels.append(label)
    for x, y in itertools.product(range(q-1), repeat=2):
        plaquette((x, y, 0), 0, 1, ('xy', x, y, 0))
    for x, y, z in itertools.product(range(q-1), range(q), range(q-1)):
        plaquette((x, y, z), 0, 2, ('xz', x, y, z))
    for x, y, z in itertools.product(range(q), range(q-1), range(q-1)):
        plaquette((x, y, z), 1, 2, ('yz', x, y, z))
    for e, (i, j) in enumerate(edges):
        if e in skeleton:
            continue
        left, right = [i], [j]
        while local_parent(left[-1]) is not None:
            left.append(local_parent(left[-1]))
        while local_parent(right[-1]) is not None:
            right.append(local_parent(right[-1]))
        start = list(vertices[i][0])
        target = vertices[j][0]
        path = left.copy()
        for axis in range(3):
            while start[axis] != target[axis]:
                start[axis] += 1 if target[axis] > start[axis] else -1
                path.append(hub(start))
        path += right[-2::-1]
        reduced = []
        for v in path:
            if len(reduced) > 1 and v == reduced[-2]:
                reduced.pop()
            else:
                reduced.append(v)
        cycles.append(reduced+[i])
        labels.append(('extra', e))
    # z_e is represented as an exact linear form in all syndrome bits.
    # Keeping every such form checks all 2**rank outcome branches at once.
    syndrome = {label: 1 << k for k, label in enumerate(labels)}
    correction = [0]*len(edges)
    for x in range(q-1):
        for y in range(q-1):
            correction[grid_edge((x, y+1, 0), 0)] = (correction[grid_edge((x, y, 0), 0)]
                                                       ^ syndrome['xy', x, y, 0])
    for z in range(q-1):
        for x, y in itertools.product(range(q-1), range(q)):
            correction[grid_edge((x, y, z+1), 0)] = (correction[grid_edge((x, y, z), 0)]
                                                       ^ syndrome['xz', x, y, z])
        for x, y in itertools.product(range(q), range(q-1)):
            correction[grid_edge((x, y, z+1), 1)] = (correction[grid_edge((x, y, z), 1)]
                                                       ^ syndrome['yz', x, y, z])
    for k, (cycle, label) in enumerate(zip(cycles, labels)):
        if label[0] == 'extra':
            value = 1 << k
            for a, b in zip(cycle[:-2], cycle[1:-1]):
                value ^= correction[index(a, b)]
            correction[label[1]] = value
    anchors = [tuple(min(vertices[v][0][axis] for v in cycle) for axis in range(3)) for cycle in cycles]
    return vertices, edges, cycles, labels, correction, anchors


def evidence():
    rows = []
    for q in (1, 2, 3, 4):
        vertices, edges, cycles, labels, corrections, anchors = plan(q)
        enc = Encoding(len(vertices), edges)
        loops = [enc.loop(path) for path in cycles]
        occupancy = defaultdict(int)
        for anchor in anchors:
            occupancy[anchor] += 1
        # Validate the producer's symbolic decoder before publishing it.
        for k, p in enumerate(loops):
            value = 0
            for e, form in enumerate(corrections):
                if p.x >> e & 1:
                    value ^= form
            if value != 1 << k:
                raise ValueError('local syndrome decoder is not a right inverse')
        rows.append(dict(side=q, cycle_count=len(cycles), maximum_length=max(len(p)-1 for p in cycles),
                         maximum_weight=max(p.weight() for p in loops),
                         maximum_slots_per_anchor=max(occupancy.values()),
                         cycle_paths_sha256=digest(cycles), signed_loops_sha256=digest([p.row() for p in loops]),
                         decoder_sha256=digest([str(x) for x in corrections]),
                         prefix_rounds=2*(q-1), local_rounds=4, colors=64, slots_upper=48,
                         loop_length_upper=16, loop_weight_upper=368))
    return rows

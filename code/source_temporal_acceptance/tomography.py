"""Producer of topology-only native completion plans and small exact controls."""
from collections import Counter, deque
from fractions import Fraction as F
import hashlib

from . import codec


def plan(graph, root):
    parent = {root: None}
    queue = deque([root])
    order = []
    while queue:
        vertex = queue.popleft()
        order.append(vertex)
        for neighbor in sorted(graph[vertex]):
            if neighbor not in parent:
                parent[neighbor] = vertex
                queue.append(neighbor)
    if set(parent) != set(graph):
        raise ValueError("disconnected completion domain")
    paths = []
    for vertex in order[1:]:
        path = [vertex]
        while path[-1] != root:
            path.append(parent[path[-1]])
        paths.append(path)
    return order, parent, paths


def summary(graph, root):
    order, parents, paths = plan(graph, root)
    stream = hashlib.sha256()
    for path in paths:
        stream.update(codec.canonical(path))
    lengths = [len(path)-1 for path in paths]
    work = sum(lengths)
    return {"vertices": len(order), "root": root, "phases": len(paths),
            "tree_sha256": codec.digest([[v, parents[v]] for v in order]),
            "paths_sha256": stream.hexdigest(), "native_means": work,
            "native_reads": 2*work, "native_writes": 2*work,
            "preparation_writes": len(order), "receiver_samples": len(order),
            "path_length_histogram": [[d, count] for d, count in sorted(Counter(lengths).items())],
            "max_path_length": max(lengths, default=0),
            "max_phase_endpoint_gain": str(2**max(lengths, default=0)),
            "all_ports_covered": True, "payload_replay": False}


def capture_graph():
    source = codec.load(codec.ROOT/codec.SUPPORT)
    graph = {v: set() for v in range(12*source["carriers"])}
    pairs = [(12*c+a, 12*c+b) for c in range(source["carriers"])
             for a, b in source["intra_carrier_seams"]]
    pairs += [(12*c+a, 12*d+b) for c, a, d, b in source["glued_pairs"]]
    for a, b in pairs:
        graph[a].add(b)
        graph[b].add(a)
    return graph


def scalar_run(paths, root, initial):
    state = dict(initial)
    samples = [state[root]]
    tape = []
    for path in paths:
        for a, b in zip(path, path[1:]):
            value = (state[a]+state[b])/2
            tape.append([a, b, str(state[a]), str(state[b]), str(value)])
            state[a] = state[b] = value
        samples.append(state[root])
    return samples, state, codec.digest(tape)


def build():
    # This connected path is a literal subset of the captured W12 seams.
    edges = [(0, 1), (1, 14), (14, 23), (23, 45)]
    graph = {v: set() for e in edges for v in e}
    for a, b in edges:
        graph[a].add(b)
        graph[b].add(a)
    order, _, paths = plan(graph, 45)
    vertices = sorted(graph)
    payloads = [[int(i == j) for i in range(5)] for j in range(5)]
    payloads += [[-3, 2, -1, 4, 7], [5, -7, 11, -13, 17], [0]*5]
    cases = []
    for values in payloads:
        samples, final, tape = scalar_run(paths, 45, dict(zip(vertices, map(F, values))))
        cases.append({"initial": values, "samples": list(map(str, samples)),
                      "final": [str(final[v]) for v in vertices], "tape_sha256": tape})
    return {"scope": "finite_connected_scalar_domain_exact_native_completion",
            "captured_plan": summary(capture_graph(), 45),
            "path_control": {"vertices": vertices, "root": 45, "order": order,
                             "paths": paths, "cases": cases,
                             "native_means_per_case": sum(len(p)-1 for p in paths)},
            "metric_radius_selected": False, "schedule_selected_by_full_a3": False,
            "physical_precision_derived": False}

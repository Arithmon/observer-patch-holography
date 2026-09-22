"""Independent topology validation, scalar replay, and receiver-only decoding."""
from collections import Counter
from fractions import Fraction as F
import hashlib

from . import codec


def need(condition, label):
    if not condition:
        raise ValueError(label)


def check_paths(vertices, edges, root, paths):
    """Semantic Plan conditions; never assume the generator constructed a tree."""
    need(root in vertices, "receiver outside domain")
    known = {root}
    count = 0
    for path in paths:
        need(type(path) is list and len(path) >= 2, "empty calibration phase")
        need(all(type(p) is int and p in vertices for p in path), "path domain")
        need(len(set(path)) == len(path), "repeated calibration vertex")
        need(path[0] not in known, "phase does not introduce a new record")
        need(path[-1] == root, "remote sample masquerades as receiver sample")
        need(set(path[1:]) <= known, "uncalibrated relay")
        need(all(frozenset((a, b)) in edges for a, b in zip(path, path[1:])), "unsupported hop")
        known.add(path[0])
        count += len(path)-1
    need(known == set(vertices), "incomplete continuation domain")
    return count


def recover(root, paths, samples):
    """Only topology and receiver samples are inputs; no remote initial values."""
    need(len(samples) == len(paths)+1, "receiver sample count")
    current = {root: F(samples[0])}
    original = dict(current)
    for path, output in zip(paths, samples[1:]):
        value = F(output)
        for port in reversed(path[1:]):
            value = 2*value-current[port]
        original[path[0]] = value
        current[path[0]] = value
        for a, b in zip(path, path[1:]):
            current[a] = current[b] = (current[a]+current[b])/2
        need(current[root] == F(output), "decoded phase disagrees with sample")
    return original, current


def full_plan(edges):
    """Reconstruct breadth layers, independently of the producer's queue."""
    vertices = set(range(15360))
    graph = {v: set() for v in vertices}
    for edge in edges:
        a, b = sorted(edge)
        graph[a].add(b)
        graph[b].add(a)
    order, parent, frontier = [45], {45: None}, [45]
    while frontier:
        candidates = {}
        for rank, p in enumerate(frontier):
            for v in graph[p]-parent.keys():
                candidates.setdefault(v, (rank, p))
        new = sorted(candidates, key=lambda v: (candidates[v][0], v))
        for v in new:
            parent[v] = candidates[v][1]
        order.extend(new)
        frontier = new
    need(set(order) == vertices, "captured support disconnected")
    paths = []
    for v in order[1:]:
        path = [v]
        while path[-1] != 45:
            need(len(path) <= len(vertices), "parent cycle")
            path.append(parent[path[-1]])
        paths.append(path)
    work = check_paths(vertices, edges, 45, paths)
    lengths = Counter(len(p)-1 for p in paths)
    stream = hashlib.sha256()
    for path in paths:
        stream.update(codec.canonical(path))
    return {"vertices": len(vertices), "root": 45, "phases": len(paths),
            "tree_sha256": codec.digest([[v, parent[v]] for v in order]),
            "paths_sha256": stream.hexdigest(), "native_means": work,
            "native_reads": 2*work, "native_writes": 2*work,
            "preparation_writes": len(vertices), "receiver_samples": len(vertices),
            "path_length_histogram": [[d, n] for d, n in sorted(lengths.items())],
            "max_path_length": max(lengths), "max_phase_endpoint_gain": str(2**max(lengths)),
            "all_ports_covered": True, "payload_replay": False}


def check(control, edges):
    c = control["path_control"]
    vertices = [0, 1, 14, 23, 45]
    paths = [[23, 45], [14, 23, 45], [1, 14, 23, 45], [0, 1, 14, 23, 45]]
    work = check_paths(set(vertices), edges, 45, c["paths"])
    payloads = [[int(i == j) for i in range(5)] for j in range(5)]
    payloads += [[-3, 2, -1, 4, 7], [5, -7, 11, -13, 17], [0]*5]
    expected_cases = []
    for values in payloads:
        state = dict(zip(vertices, map(F, values)))
        samples = [state[45]]
        tape = []
        for path in paths:
            for i in range(len(path)-1):
                a, b = path[i:i+2]
                before_a, before_b = state[a], state[b]
                state[a] = (before_a+before_b)*F(1, 2)
                state[b] = state[a]
                tape.append([a, b, str(before_a), str(before_b), str(state[a])])
            samples.append(state[45])
        decoded, evolved = recover(45, paths, samples)
        need([decoded[v] for v in vertices] == values, "receiver fails original records")
        need(evolved == state, "receiver fails evolved calibration")
        expected_cases.append({"initial": values, "samples": list(map(str, samples)),
                               "final": [str(state[v]) for v in vertices],
                               "tape_sha256": codec.digest(tape)})
    expected = {"scope": "finite_connected_scalar_domain_exact_native_completion",
                "captured_plan": full_plan(edges),
                "path_control": {"vertices": vertices, "root": 45, "order": [45, 23, 14, 1, 0],
                                 "paths": paths, "cases": expected_cases,
                                 "native_means_per_case": work},
                "metric_radius_selected": False, "schedule_selected_by_full_a3": False,
                "physical_precision_derived": False}
    codec.equal(control, expected, "native completion plan, state, work or scope")
    return {"captured_topology_vertices": 15360,
            "captured_plan_means": expected["captured_plan"]["native_means"], "replayed_histories": len(payloads),
            "replayed_means": len(payloads)*work}

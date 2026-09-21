"""Compare two declared finite twelve-port gluing graphs.

An edge is ``[carrier, port, carrier, port]``. Carrier identifiers form a
contiguous nonempty range and every carrier occurs in an edge. The input is
the complete declared edge list; this module neither recovers that list from
events nor asserts a canonical geometric realization.
"""

from __future__ import annotations

from collections import Counter, deque
import hashlib
import json
from typing import Iterable, Sequence


PORTS_PER_CARRIER = 12
REMOVED_EDGES = ((120, 5, 123, 4), (132, 9, 135, 11))
ADDED_EDGES = ((120, 5, 135, 11), (123, 4, 132, 9))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True) + "\n").encode("ascii")


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _graph(edges: Iterable[Sequence[int]]) -> tuple:
    rows = []
    pairs, occupied, vertices = set(), set(), set()
    for raw in edges:
        _require(isinstance(raw, (list, tuple)) and len(raw) == 4,
                 "edge must have four integer fields")
        _require(all(type(value) is int for value in raw), "noninteger edge field")
        a, p, b, q = raw
        _require(0 <= a < b, "edge endpoints must be distinct and increasing")
        _require(0 <= p < PORTS_PER_CARRIER and 0 <= q < PORTS_PER_CARRIER,
                 "edge port outside the twelve-port interface")
        _require((a, b) not in pairs, "duplicate carrier edge")
        _require((a, p) not in occupied and (b, q) not in occupied,
                 "carrier port occurs in more than one edge")
        pairs.add((a, b))
        occupied.update(((a, p), (b, q)))
        vertices.update((a, b))
        rows.append((a, p, b, q))
    _require(bool(rows), "empty graph")
    count = max(vertices) + 1
    # Check the cardinality before allocating a range indexed by input IDs.
    _require(count == len(vertices), "carrier identifiers are not contiguous")
    adjacency = [set() for _ in range(count)]
    for a, _, b, _ in rows:
        adjacency[a].add(b)
        adjacency[b].add(a)
    return rows, adjacency, occupied


def _bfs(adjacency: list[set[int]], root: int) -> tuple[list[int], list[int]]:
    """Ascending-neighbour BFS, retaining all discoveries and queue positions."""
    parents = [-1] * len(adjacency)
    parents[root] = root
    pending, order = deque([root]), []
    while pending:
        vertex = pending.popleft()
        order.append(vertex)
        for neighbour in sorted(adjacency[vertex]):
            if parents[neighbour] == -1:
                parents[neighbour] = vertex
                pending.append(neighbour)
    return parents, order


def _invariants(rows: list[tuple], adjacency: list[set[int]], occupied: set) -> dict:
    degrees = list(map(len, adjacency))
    # Each triangle contributes once at each of its three unordered edges.
    incidences = sum(len(adjacency[a] & adjacency[b])
                     for a, neighbours in enumerate(adjacency)
                     for b in neighbours if a < b)
    _require(incidences % 3 == 0, "triangle incidence arithmetic")
    return {
        "carriers": len(adjacency),
        "glued_edges": len(rows),
        "triangles": incidences // 3,
        "degree_histogram": [list(row) for row in sorted(Counter(degrees).items())],
        "degree_vector_sha256": _digest(degrees),
        "occupied_ports_sha256": _digest(sorted(occupied)),
        "canonical_edges_sha256": _digest(sorted(rows)),
    }


def make_alternative(edges: Iterable[Sequence[int]]) -> list[list[int]]:
    """Apply the declared witness's two-switch without changing its port roster.

    The input is copied. All other edge rows retain their original positions.
    A different graph lacking either original edge is rejected.
    """
    rows, _, _ = _graph(edges)
    alternative = list(rows)
    for old, new in zip(REMOVED_EDGES, ADDED_EDGES, strict=True):
        _require(old in rows, "two-switch source edge missing")
        alternative[rows.index(old)] = new
    _graph(alternative)
    return [list(row) for row in alternative]


def analyze_completions(original_edges: Iterable[Sequence[int]],
                        alternate_edges: Iterable[Sequence[int]],
                        hosts: Iterable[int]) -> dict:
    """Return exact graph and BFS comparisons, with explicit success flags.

    Malformed graph/interface data raises ``ValueError``. A well-formed pair
    that fails a witness condition returns false for that condition and for
    ``witness_valid``; callers must require that conjunction. The finite
    triangle count distinguishes abstract graph isomorphism classes, while
    degree/port and BFS equalities compare the declared carrier labels.
    """
    old_rows, old_graph, old_ports = _graph(original_edges)
    new_rows, new_graph, new_ports = _graph(alternate_edges)
    roots = list(hosts)
    _require(bool(roots) and all(type(root) is int for root in roots),
             "nonempty integer host list required")
    _require(len(roots) == len(set(roots)), "duplicate source host")
    _require(all(0 <= root < min(len(old_graph), len(new_graph)) for root in roots),
             "source host outside a graph")
    roots.sort()
    old = _invariants(old_rows, old_graph, old_ports)
    new = _invariants(new_rows, new_graph, new_ports)
    old_hash, new_hash = hashlib.sha256(), hashlib.sha256()
    parent_equal = queue_equal = True
    old_connected = new_connected = True
    first_difference = None
    for root in roots:
        old_parents, old_queue = _bfs(old_graph, root)
        new_parents, new_queue = _bfs(new_graph, root)
        old_connected &= len(old_queue) == len(old_graph)
        new_connected &= len(new_queue) == len(new_graph)
        parent_equal &= old_parents == new_parents
        queue_equal &= old_queue == new_queue
        for digest, parents, queue in ((old_hash, old_parents, old_queue),
                                      (new_hash, new_parents, new_queue)):
            digest.update(_canonical({"root": root, "parents": parents, "queue": queue}))
        if first_difference is None and (old_parents != new_parents or old_queue != new_queue):
            first_difference = {
                "root": root,
                "parent_maps_equal": old_parents == new_parents,
                "queue_orders_equal": old_queue == new_queue,
            }
    comparison = {
        "carrier_counts_equal": len(old_graph) == len(new_graph),
        "glued_edge_counts_equal": len(old_rows) == len(new_rows),
        "per_carrier_degrees_equal": list(map(len, old_graph)) == list(map(len, new_graph)),
        "occupied_carrier_ports_equal": old_ports == new_ports,
        "original_connected": old_connected,
        "alternate_connected": new_connected,
        "all_host_BFS_parent_maps_equal": parent_equal,
        "all_host_BFS_queue_orders_equal": queue_equal,
        "nonisomorphic_by_triangle_count": old["triangles"] != new["triangles"],
    }
    return {
        "schema": "oph.federation_recovery.topology_comparison.v1",
        "original": old,
        "alternate": new,
        "hosts": roots,
        "comparison": comparison,
        "original_BFS_sha256": old_hash.hexdigest(),
        "alternate_BFS_sha256": new_hash.hexdigest(),
        "first_BFS_difference": first_difference,
        "witness_valid": all(comparison.values()),
        "scope": "Declared finite twelve-port gluing graphs; no event-log recovery, "
                 "native execution, or canonical geometric realization is certified here.",
    }

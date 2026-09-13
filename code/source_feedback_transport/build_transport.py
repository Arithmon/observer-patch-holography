"""Reusable local classical record transport over canonical W12 seams.

Copy/reset feedback is an explicit additional local operation. It is never
called a consequence of pair averaging or a quantum channel.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def routes(support, lengths):
    adjacency = [set() for _ in range(12 * support["carriers"])]
    for carrier in range(support["carriers"]):
        for a, b in support["intra_carrier_seams"]:
            u, v = 12 * carrier + a, 12 * carrier + b
            adjacency[u].add(v)
            adjacency[v].add(u)
    for c, a, e, b in support["glued_pairs"]:
        u, v = 12 * c + a, 12 * e + b
        adjacency[u].add(v)
        adjacency[v].add(u)
    c, a, _, _ = support["glued_pairs"][0]
    start = 12 * c + a
    previous, distance, queue = {start: None}, {start: 0}, deque([start])
    answer = {}
    while queue:
        u = queue.popleft()
        if distance[u] in lengths and distance[u] not in answer:
            path = [u]
            while previous[path[-1]] is not None:
                path.append(previous[path[-1]])
            answer[distance[u]] = list(reversed(path))
        for v in sorted(adjacency[u]):
            if v not in previous:
                previous[v], distance[v] = u, distance[u] + 1
                queue.append(v)
    return [answer[d] for d in lengths]


class Recorder:
    def __init__(self):
        self.cells, self.events = {}, []

    def emit(self, op, owner, reads, writes, label=None, **arguments):
        before = [{"register": key, **self.cells[key]} for key in reads]
        result = []
        for key, value, immutable in writes:
            if key in self.cells and self.cells[key]["immutable"]:
                raise ValueError("attempt to overwrite an immutable archive")
            version = self.cells[key]["version"] + 1 if key in self.cells else 0
            cell = {"owner": owner if op != "mean" else int(key.split("/")[1]),
                    "value": str(Q(value)), "version": version,
                    "writer": len(self.events), "immutable": immutable}
            result.append({"register": key, **cell})
        event = {"id": len(self.events), "op": op, "owner": owner,
                 "label": label, "arguments": arguments,
                 "reads": before, "writes": result,
                 "parents": sorted({row["writer"] for row in before}),
                 "previous_hash": self.events[-1]["event_hash"] if self.events else "0" * 64}
        event["event_hash"] = digest(event)
        self.events.append(event)
        for row in result:
            self.cells[row["register"]] = {k: v for k, v in row.items() if k != "register"}

    def value(self, key):
        return Q(self.cells[key]["value"])

    def hop(self, message, index, u, v, archive, origin):
        incoming = f"buffer/{message}/{index}"
        self.emit("export", u, [archive], [(f"port/{u}", self.value(archive), False)],
                  origin, message=message, index=index)
        self.emit("reset", v, [f"zero/{v}"], [(f"port/{v}", 0, False)],
                  message=message, index=index, phase="prepare")
        mean = (self.value(f"port/{u}") + self.value(f"port/{v}")) / 2
        self.emit("mean", None, [f"port/{u}", f"port/{v}"],
                  [(f"port/{u}", mean, False), (f"port/{v}", mean, False)],
                  origin, message=message, index=index, seam=[u, v])
        self.emit("capture", v, [f"port/{v}"], [(incoming, 2 * mean, True)],
                  origin, message=message, index=index)
        for endpoint in (u, v):
            self.emit("reset", endpoint, [f"zero/{endpoint}"],
                      [(f"port/{endpoint}", 0, False)],
                      message=message, index=index, phase="restore")
        return incoming

    def transport(self, name, route, archive, origin):
        first = len(self.events)
        for index, (u, v) in enumerate(zip(route, route[1:])):
            archive = self.hop(name, index, u, v, archive, origin)
        return {"message": name, "route": route, "origin": origin,
                "first_event": first, "event_count": len(self.events) - first,
                "delivered_register": archive}


def episode(path, variant, spec):
    rec = Recorder()
    for position, port in enumerate(path):
        rec.emit("initialize_zero", port, [], [(f"zero/{port}", 0, True)])
        nuisance = Q(100 + 3 * position, 17)
        if variant == "scratch_intervention":
            nuisance += Q(position + 1, 19)
        rec.emit("initialize_scratch", port, [], [(f"port/{port}", nuisance, False)])
    locations = {"A": path[0], "B": path[-1], "C": path[max(1, (len(path)-1)//2)]}
    logical, requests = [], []
    for name, value in spec["archive_values"].items():
        value = Q(value)
        if name == "A0" and variant == "source_intervention":
            value += Q(spec["source_intervention"])
        rec.emit("publish", locations[name[0]], [], [(f"archive/{name}", value, True)], name, name=name)
        logical.append({"name": name, "owner": locations[name[0]], "parents": [],
                        "event_id": len(rec.events)-1, "value": str(value)})

    recipes = [
        ("to_B1", "A0", path, "B1", "B0"),
        ("to_A1", "B1", list(reversed(path)), "A1", "A0"),
        ("to_C1", "A1", path[:path.index(locations["C"])+1], "C1", "C0"),
        ("old_to_B2", "A0", path, "B2", "B1"),
        ("old_to_B3", "A0", path, "B3", "B2"),
    ]
    for message, origin, route, name, local in recipes:
        request = rec.transport(message, route, f"archive/{origin}", origin)
        requests.append(request)
        reads = [request["delivered_register"], f"archive/{local}"]
        value = 1 + sum(rec.value(key) for key in reads)
        delta = Q(spec["branch_intervention"]) if name == "A1" and variant == "branch_intervention" else Q(0)
        value += delta
        rec.emit("commit", route[-1], reads, [(f"archive/{name}", value, True)],
                  name, name=name, intervention=str(delta))
        logical.append({"name": name, "owner": route[-1], "parents": sorted([origin, local]),
                        "event_id": len(rec.events)-1, "value": str(value)})
    hops = sum(len(row["route"])-1 for row in requests)
    costs = {"hops": hops, "means": hops, "exports": hops, "captures": hops,
             "resets": 3*hops, "transport_events": 6*hops,
             "initialization_events": 2*len(path)+3, "logical_commits": 5,
             "total_events": len(rec.events),
             "register_reads": sum(len(e["reads"]) for e in rec.events),
             "register_writes": sum(len(e["writes"]) for e in rec.events),
             "protected_scalar_registers": sum(c["immutable"] for c in rec.cells.values()),
             "scratch_scalar_registers": len(path),
             "largest_written_scalar_bits": max(abs(Q(w["value"]).numerator).bit_length()+Q(w["value"]).denominator.bit_length()+1 for e in rec.events for w in e["writes"])}
    return {"path": path, "variant": variant, "requests": requests, "logical_events": logical,
            "events": rec.events, "event_root": rec.events[-1]["event_hash"], "costs": costs}


def target_costs():
    parent = ROOT / "evidence/source_net_causal_poset/carrier_source_net_receipt.json"
    source = json.loads(parent.read_text())
    rows = []
    for row in source["levels"]:
        if row["q"] in (13, 21):
            n, rounds, edges = row["site_count"], row["rounds"], row["neighbours"]["undirected_spatial_edges"]
            rows.append({"q": row["q"], "sites": n, "rounds": rounds, "undirected_metric_edges": edges,
                         "logical_reads": (2*edges+n)*rounds,
                         "nonself_reads": 2*edges*rounds,
                         "six_event_unicast_transport_lower_bound": 12*edges*rounds,
                         "scope": "Declared compiler target count only; assumes distinct source/destination scratch vertices; excludes routing distance, local reads and commits. Not a full execution or a lower bound for every multicast/compiler design."})
    return {"source_path": str(parent.relative_to(ROOT)), "source_sha256": hashlib.sha256(parent.read_bytes()).hexdigest(), "rows": rows}


def build():
    spec = json.loads((HERE / "specification.json").read_text())
    support_path = ROOT / spec["support"]
    support = json.loads(support_path.read_text())
    runs = [episode(path, variant, spec) for path in routes(support, spec["path_lengths"]) for variant in spec["variants"]]
    return {"schema": "oph.source_feedback_transport.receipt.v1", "specification": spec,
            "support_sha256": hashlib.sha256(support_path.read_bytes()).hexdigest(),
            "scope": "CLASSICAL_REUSABLE_LOCAL_FEEDBACK_TRANSPORT__SUPPLIED_COPY_RESET_LAW__NO_PHYSICAL_OR_QUANTUM_IDENTIFICATION",
            "noise_contract": {"per_hop": "export + reset + 2*mean + 2*read + decode", "depth_d": "initial_archive_error + d*(export + reset + 2*mean + 2*read + decode)",
                               "scope": "Separately bounded additive errors; export includes archive retention/read/write error relative to its retained payload. Mean error bounds the receiver coordinate. Exact deterministic schedule; no measured physical noise limits."},
            "compiler_targets_not_executed": target_costs(), "episodes": runs,
            "summary": {"episodes": len(runs), "events": sum(len(run["events"]) for run in runs), "hops": sum(run["costs"]["hops"] for run in runs)}}


if __name__ == "__main__":
    output = HERE / "transport_receipt.json"
    output.write_bytes(canonical(build()))
    print(output)

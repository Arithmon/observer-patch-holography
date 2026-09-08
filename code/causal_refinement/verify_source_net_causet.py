"""Independent source-coded causet replay in Q(sqrt(5)).

No producer or sibling repository is imported. Source records are checked
by actual twelve-port incidence followed by the primitive six-axis Gram
readback. Every finite graph and every read/write trace is reconstructed.
This certifies the declared finite model, not native physical selection.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "source_net_causet_receipt.json"
PIN_PATHS = (
    "code/causal_refinement/source_net_causet.py",
    "code/causal_refinement/verify_source_net_causet.py",
    "code/causal_refinement/test_source_net_causet.py",
    "Lean/Screen/PrimitivePortFrameQuotient.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/Screen/SeamCurrentHomogeneousAction.lean",
    "Lean/Screen/PortFrameGram.lean",
)


def require(c, label):
    if not c:
        raise ValueError(label)


def exact(x, y, label="receipt"):
    require(type(x) is type(y), label+": wrong type")
    if isinstance(y, dict):
        require(x.keys() == y.keys(), label+": wrong fields")
        for k in y:
            exact(x[k], y[k], label+"/"+k)
    elif isinstance(y, list):
        require(len(x) == len(y), label+": wrong length")
        for i, (a, b) in enumerate(zip(x, y)):
            exact(a, b, label+"/"+str(i))
    else:
        require(x == y, label+": wrong value")


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def forbidden_token(token):
    raise ValueError("noninteger JSON numeric token")


def load(path=OUTPUT):
    data = Path(path).read_bytes()
    require(len(data) <= 2_000_000, "receipt size limit")
    return json.loads(data.decode("utf-8"), object_pairs_hook=pairs,
                      parse_float=forbidden_token, parse_constant=forbidden_token)


def packed(x):
    return (json.dumps(x, ensure_ascii=True, allow_nan=False, sort_keys=True,
                       separators=(",", ":"))+"\n").encode("ascii")


def hashed(x):
    return sha256(packed(x)).hexdigest()


def plus(x, y):
    return (x[0]+y[0], x[1]+y[1])


def minus(x, y):
    return (x[0]-y[0], x[1]-y[1])


def mult(x, y):
    return (x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0])


def times(c, x):
    return (c*x[0], c*x[1])


def sgn(x):
    a, b = x
    if b == 0:
        return (a > 0)-(a < 0)
    if a == 0:
        return (b > 0)-(b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    d = a*a-5*b*b
    return ((d > 0)-(d < 0))*(1 if a > 0 else -1)


def encoded_phi(x):
    return [str(F(x[0]-x[1])), str(F(2*x[1]))]


def table(relative, name):
    text = (ROOT/relative).read_text(encoding="utf-8")
    match = re.search(r"def "+name+r"\s*:[^=]*:=\s*!\[(.*?)\]", text, re.S)
    require(match is not None, "missing primitive source table "+name)
    return [int(v.strip()) for v in match[1].split(",")]


def source_basis():
    seam_file = "Lean/Screen/SeamCurrentCarrierQuotient.lean"
    frame_file = "Lean/Screen/PrimitivePortFrameQuotient.lean"
    left, right = table(seam_file, "seamLeft"), table(seam_file, "seamRight")
    relabel = table(frame_file, "sourceToRERPort")
    positive = [relabel[i] for i in table(frame_file, "positiveSourcePort")]
    axes = (((-2, 0), (1, 1), (0, 0)), ((2, 0), (1, 1), (0, 0)),
            ((0, 0), (-2, 0), (1, 1)), ((0, 0), (2, 0), (1, 1)),
            ((1, 1), (0, 0), (-2, 0)), ((1, 1), (0, 0), (2, 0)))
    gram_text = (ROOT/"Lean/Screen/PortFrameGram.lean").read_text(encoding="utf-8")
    adjacency = {int(i): [int(v) for v in row.split(",")]
                 for i, row in re.findall(r"\| (\d+) => \[([\d, ]+)\]", gram_text)}
    require(len(adjacency) == 12 and len(left) == len(right) == 30, "source table dimensions")
    for i in range(6):
        for j in range(6):
            u, v = positive[i], positive[j]
            g5 = (5, 0) if u == v else ((-5, 0) if u+v == 11 else
                                       ((0, 1) if v in adjacency[u] else (0, -1)))
            dot = (0, 0)
            for a, b in zip(axes[i], axes[j]):
                dot = plus(dot, mult(a, b))
            require(times(F(1, 4), dot) == mult(times(F(1, 5), g5), (F(5, 2), F(1, 2))),
                    "primitive axes disagree with source Gram")
    return left, right, positive, axes


def graph_and_source(n):
    require(type(n) is int and 4 <= n <= 7, "finite regulator domain")
    fib = [0, 1]
    for _ in range(n):
        fib.append(fib[-1]+fib[-2])
    q, p = fib[n], fib[n+1]
    floors = []
    for b in range(q):
        floor = 0
        while sgn((b-2*(floor+1), b)) >= 0:
            floor += 1
        floors.append(floor)
        require(floor == b*p//q, "golden rational cell mismatch")
    # Twice each orbit coordinate, in the independent sqrt(5) basis.
    two_xi = [(b-2*floors[b], b) for b in range(q)]
    sites = list(product(range(q), repeat=3))
    pair_sq4 = [[mult(minus(x, y), minus(x, y)) for y in two_xi] for x in two_xi]
    neighbors = [[i] for i in range(q**3)]
    # Independent unordered all-pairs graph; no producer neighborhood pruning.
    for i, a in enumerate(sites):
        for j in range(i):
            b = sites[j]
            d2 = plus(plus(pair_sq4[a[0]][b[0]], pair_sq4[a[1]][b[1]]), pair_sq4[a[2]][b[2]])
            if sgn(minus(times(q, d2), (4, 0))) <= 0:
                neighbors[i].append(j)
                neighbors[j].append(i)
    for row in neighbors:
        row.sort()
    left, right, positive, axes = source_basis()
    records, costs = [], []
    for b in sites:
        a = [-floors[i] for i in b]
        # Explicit current section, followed by twelve-port incidence; z is
        # obtained from primitive source records rather than control(a,b).
        currents = [-a[2]-b[0], a[0]+b[0]+b[2], b[0]+b[1]-b[2],
                    -a[0]-b[0]+b[2], a[2]+b[1]-b[2], a[1]+b[2]]
        load = [0]*12
        for seam, value in zip((5, 8, 12, 13, 22, 24), currents):
            load[left[seam]] -= value
            load[right[seam]] += value
        require(sum(load) == 0, "seam incidence did not conserve load")
        z = [load[v]-load[11-v] for v in positive]
        require(sum(z) % 2 == 0, "nonconservative D6 record")
        for k in range(3):
            coordinate = (0, 0)
            for value, axis in zip(z, axes):
                coordinate = plus(coordinate, times(value, axis[k]))
            require(coordinate == times(2, two_xi[b[k]]), "record Gram readback mismatch")
        records.append([list(b), z, currents])
        costs.append(sum(map(abs, currents)))
    require(max(costs) <= 27*(q-1), "source word cost bound")
    return q, p, two_xi, sites, pair_sq4, neighbors, records, costs


def execute(neighbors, layers, reversed_order=False, intervention=None):
    count = len(neighbors)
    registers = {}
    chain = bytes(32)
    layer_hashes, sums = [], []
    for layer in range(layers+1):
        input_registers = dict(registers)
        order = range(count-1, -1, -1) if reversed_order else range(count)
        for site in order:
            reads = []
            if layer:
                for parent in neighbors[site]:
                    version, writer, value = input_registers[parent]
                    require(version == layer and writer == [layer-1, parent], "actual writer/version mismatch")
                    reads.append([parent, version, writer, value])
            value = site+1+(site == intervention) if layer == 0 else 1+sum(r[-1] for r in reads)
            write = [site, layer+1, [layer, site], value]
            registers[site] = (write[1], write[2], value)
            chain = sha256(chain+packed([[layer, site], reads, write])).digest()
        values = [registers[i][2] for i in range(count)]
        layer_hashes.append(hashed(values))
        sums.append(sum(values))
    return {"audit_trace_sha256": chain.hex(), "layer_value_sha256": layer_hashes,
            "layer_value_sums": sums, "event_count": count*(layers+1),
            "authenticated_read_count": layers*sum(len(r) for r in neighbors)}


def breadth(neighbors, start):
    seen, frontier = {start: 0}, {start}
    while frontier:
        level = seen[next(iter(frontier))]+1
        future = {j for i in frontier for j in neighbors[i]}-seen.keys()
        seen.update({j: level for j in future})
        frontier = future
    require(len(seen) == len(neighbors), "disconnected finite graph")
    return [seen[i] for i in range(len(neighbors))]


def verify_level(row):
    require(type(row) is dict, "level object")
    n = row.get("fibonacci_index")
    q, p, xi2, sites, ds4, graph, records, costs = graph_and_source(n)
    orbit = [times(F(1, 2), x) for x in xi2]
    permutation = [b*p % q for b in range(q)]
    require(sorted(permutation) == list(range(q)), "golden permutation")
    order = sorted(range(q), key=lambda i: permutation[i])
    for i in range(q):
        error = minus(orbit[i], (F(permutation[i], q), 0))
        require(sgn(minus((F(1, q*q), 0), mult(error, error))) > 0, "quadrature cell displacement")
    radius = minus((1, 0), orbit[order[-1]])
    for a, b in zip(order, order[1:]):
        half_gap = times(F(1, 2), minus(orbit[b], orbit[a]))
        if sgn(minus(half_gap, radius)) > 0:
            radius = half_gap
    h2 = times(3, mult(radius, radius))
    h_a2 = times(q, h2)
    ratio = row.get("h_over_a_upper")
    require(type(ratio) is str, "rational fill ratio")
    ratio = F(ratio)
    require(ratio >= 0 and (ratio*10**6).denominator == 1, "fill ratio grid")
    require(sgn(minus((ratio*ratio, 0), h_a2)) >= 0, "unsafe fill upper bound")
    require(ratio == 0 or sgn(minus(((ratio-F(1, 10**6))**2, 0), h_a2)) < 0,
            "noncanonical outward fill bound")
    inner = max(F(0), 1-2*ratio)
    center_axis = 0
    for b in range(1, q):
        d, old = minus(orbit[b], (F(1, 2), 0)), minus(orbit[center_axis], (F(1, 2), 0))
        if sgn(minus(mult(d, d), mult(old, old))) < 0:
            center_axis = b
    center = sites.index((center_axis,)*3)
    layers = isqrt(q)+(isqrt(q)**2 < q)
    diagnostics = []
    for start in sorted({0, center, q**3-1}):
        steps = breadth(graph, start)
        for k in range(1, layers+1):
            reachable, missed = [], []
            for i, b in enumerate(sites):
                sq4 = (0, 0)
                for axis in range(3):
                    sq4 = plus(sq4, ds4[sites[start][axis]][b[axis]])
                outer = sgn(minus(times(q, sq4), (4*k*k, 0))) <= 0
                small = sgn(minus(times(q, sq4), (4*k*k*inner*inner, 0))) <= 0
                if steps[i] <= k:
                    require(outer, "reached outside outer cone")
                    reachable.append(i)
                else:
                    require(not small, "inner cone unreachable")
                    if outer:
                        missed.append(i)
            diagnostics.append({"start": start, "layers": k,
                "reachable_count": len(reachable), "reachable_ids_sha256": hashed(reachable),
                "outer_cone_violations": 0, "certified_inner_cone_misses": 0,
                "exact_finite_cone_missing_count": len(missed), "missing_ids_sha256": hashed(missed)})
    hops = breadth(graph, center)
    clearance = orbit[center_axis]
    if sgn(minus(minus((1, 0), clearance), clearance)) < 0:
        clearance = minus((1, 0), clearance)
    intervals = []
    for k in range(1, layers+1):
        counts = [sum(d <= j and d <= k-j for d in hops) for j in range(k+1)]
        intervals.append({"layers": k, "inclusive_event_count": sum(counts), "counts_by_layer": counts,
            "continuum_diamond_inside_cube": sgn(minus(times(4*q, mult(clearance, clearance)), (k*k, 0))) >= 0,
            "count_over_L4_times_sqrt_q": str(F(sum(counts), q**3)),
            "diamond_volume_over_pi_L4": str(F(k**4, 24*q*q)),
            "comparison_role": "finite_diagnostic_not_dimension_fit_or_convergence_test"})
    forward = execute(graph, layers)
    reverse = execute(graph, layers, reversed_order=True)
    intervention = execute(graph, layers, intervention=center)
    require(forward["layer_value_sha256"] == reverse["layer_value_sha256"], "schedule dependence")
    require(forward["audit_trace_sha256"] != reverse["audit_trace_sha256"], "lost actual serial audit order")
    delta = {center: 1}
    response = []
    for k in range(layers+1):
        support = sorted(delta)
        require(support == [i for i, d in enumerate(hops) if d <= k], "intervention support mismatch")
        require(sum(delta.values()) == intervention["layer_value_sums"][k]-forward["layer_value_sums"][k],
                "intervention amplitude mismatch")
        response.append({"layer": k, "support_count": len(support),
            "support_ids_sha256": hashed(support), "positive_integer_delta_sum": sum(delta.values()),
            "positive_integer_delta_maximum": max(delta.values())})
        future = {}
        for source, amplitude in delta.items():
            for target in graph[source]:
                future[target] = future.get(target, 0)+amplitude
        delta = future
    expected = {"fibonacci_index": n, "q": q, "p": p, "layer_steps": layers,
        "orbit_Qphi": [encoded_phi(x) for x in orbit], "grid_permutation": permutation,
        "record_count": q**3, "source_records_sha256": hashed(records),
        "source_record_examples": [records[i] for i in sorted({0, center, q**3-1})],
        "maximum_word_length": max(costs), "sum_word_lengths": sum(costs), "word_length_bound": 27*(q-1),
        "whole_cube_fill_h_squared_over_L2_Qphi": encoded_phi(h2),
        "one_dimensional_fill_over_L_Qphi": encoded_phi(radius),
        "quadrature_displacement_H_squared_over_L2": str(F(12, q*q)),
        "h_over_a_upper": str(ratio), "certified_inner_speed_lower": str(inner),
        "positive_inner_cone": sgn(minus((1, 0), times(4*q, h2))) > 0,
        "radius_squared_over_L2": str(F(1, q)), "layer_time_equals_radius": True,
        "neighbors_including_wait_sha256": hashed(graph),
        "undirected_spatial_edges": (sum(map(len, graph))-q**3)//2,
        "minimum_neighbor_count_including_wait": min(map(len, graph)),
        "maximum_neighbor_count_including_wait": max(map(len, graph)),
        "exact_width": q**3, "exact_height_in_events": layers+1,
        "chain_cover": "one wait chain per source record; layers are antichains",
        "forward_execution": forward, "within_layer_reversed_execution": reverse,
        "center_plus_one_intervention": intervention, "intervention_source_id": center,
        "intervention_response": response,
        "reachability_probes": diagnostics, "center_intervals": intervals}
    exact(row, expected, "level"+str(q))
    return {"q": q, "source_records": q**3, "width": q**3, "height": layers+1,
            "reads_per_full_trace": forward["authenticated_read_count"],
            "positive_inner_cone": expected["positive_inner_cone"], "inner_speed_lower": str(inner)}


def verify(packet):
    require(type(packet) is dict, "packet object")
    require(type(packet.get("source_pins")) is dict, "source pins object")
    exact(packet["source_pins"], {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in PIN_PATHS}, "source_pins")
    expected_scope = {
        "conservative_source_record_positions": True, "exact_gram_metric_edge_decisions": True,
        "full_finite_read_write_traces_executed": True, "full_traces_stored": False,
        "intervention_replayed": True, "cube_product_population_and_neighbor_law_supplied": True,
        "model_clock_supplied": True, "native_repair_selected": False,
        "primitive_source_words_executed_as_repairs": False, "wave_action_derived": False,
        "physical_clock_or_spacetime_identified": False, "dimension_statistic_used_as_acceptance": False,
        "poisson_sprinkling": False, "raw_audit_chain_is_semantic_causal_order": False,
        "finite_energy_signal_or_quantum_field_established": False,
        "finite_runs_demonstrate_fixed_time_asymptotic_limit": False}
    expected = {"schema": "oph.source-coded-golden-causet.v1", "scope": expected_scope,
        "source_scale_L_squared_Qphi": ["12/5", "-4/5"],
        "window": "[0,L]^3 in the existing source-Gram coordinate witness",
        "source_current_seams": [5, 8, 12, 13, 22, 24],
        "event_law": "q_i(0)=i+1; q_i(j)=1+sum(q_r(j-1) for every metric neighbor including i)",
        "event_order": "closure of consecutive-layer authenticated local reads; audit digests are metadata",
        "sampling_measure": "L^3/q^3 per site, layer weight a_q=L/sqrt(q)",
        "finite_schedule": "K_q=ceil(sqrt(q)); T_q=K_q L/sqrt(q) tends to L; four traces alone do not establish the limit",
        "execution_format": "all finite rows executed and hash-chained; full value/read/write trace replayed from templates"}
    exact({k:v for k,v in packet.items() if k not in ("levels", "source_pins")}, expected)
    levels = packet.get("levels")
    require(type(levels) is list and len(levels) == 4, "four frozen levels")
    exact([r.get("fibonacci_index") if type(r) is dict else None for r in levels], [4, 5, 6, 7], "level census")
    rows = [verify_level(r) for r in levels]
    return {"accepted": True, "scope": "CONDITIONAL_SOURCE_CODED_CAUSETS__SUPPLIED_POPULATION_READ_LAW_AND_CLOCK",
            "native_physical_spacetime_selected": False, "levels": rows,
            "packet_sha256": sha256(packed(packet)).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), indent=2, sort_keys=True))

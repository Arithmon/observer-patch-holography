"""Conditional causets on golden-orbit conservative source records.

Positions are the existing Gram readback of actual D6 controls, not lattice
nodes renamed as records. The cube window, product population, all-neighbour
layer law, integer signal, and clock a=L/sqrt(q) are supplied. They are not
selected by the native repair law or derived from a wave action. Every layer
of every finite run is executed; full read/write traces are digested rather
than stored. See the independently replayed geometric and source-word data.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "source_net_causet_receipt.json"
LEVELS = (4, 5, 6, 7)
SEAMS = (5, 8, 12, 13, 22, 24)
PINS = (
    "code/causal_refinement/source_net_causet.py",
    "code/causal_refinement/verify_source_net_causet.py",
    "code/causal_refinement/test_source_net_causet.py",
    "Lean/Screen/PrimitivePortFrameQuotient.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/Screen/SeamCurrentHomogeneousAction.lean",
    "Lean/Screen/PortFrameGram.lean",
)
SCHEMA = "oph.source-coded-golden-causet.v1"
SCOPE = {
    "conservative_source_record_positions": True,
    "exact_gram_metric_edge_decisions": True,
    "full_finite_read_write_traces_executed": True,
    "full_traces_stored": False,
    "intervention_replayed": True,
    "cube_product_population_and_neighbor_law_supplied": True,
    "model_clock_supplied": True,
    "native_repair_selected": False,
    "primitive_source_words_executed_as_repairs": False,
    "wave_action_derived": False,
    "physical_clock_or_spacetime_identified": False,
    "dimension_statistic_used_as_acceptance": False,
    "poisson_sprinkling": False,
    "raw_audit_chain_is_semantic_causal_order": False,
    "finite_energy_signal_or_quantum_field_established": False,
    "finite_runs_demonstrate_fixed_time_asymptotic_limit": False,
}


def canonical(x):
    return (json.dumps(x, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def digest(x):
    return sha256(canonical(x)).hexdigest()


def sign(x):
    """Exact sign of a+b phi, also for rational coefficients."""
    a, b = x
    c = 2*a+b
    if b == 0:
        return (c > 0)-(c < 0)
    if c >= 0 and b > 0:
        return 1
    if c <= 0 and b < 0:
        return -1
    d = c*c-5*b*b
    return ((d > 0)-(d < 0))*(1 if c > 0 else -1)


def add(x, y):
    return (x[0]+y[0], x[1]+y[1])


def sub(x, y):
    return (x[0]-y[0], x[1]-y[1])


def scale(c, x):
    return (c*x[0], c*x[1])


def square(x):
    a, b = x
    return (a*a+b*b, 2*a*b+b*b)


def encode(x):
    return [str(F(v)) for v in x]


def fibonacci(n):
    if type(n) is not int or not 4 <= n <= 7:
        raise ValueError("finite packet levels are integers 4 through 7")
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a, b


def orbit(n):
    q, p = fibonacci(n)
    return q, p, [(-((b+isqrt(5*b*b))//2), b) for b in range(q)]


def source_control(a, b):
    if (not isinstance(a, (list, tuple)) or not isinstance(b, (list, tuple))
            or len(a) != 3 or len(b) != 3 or any(type(v) is not int for v in (*a, *b))):
        raise ValueError("source controls require two integer triples")
    return (b[1]-a[0], b[1]+a[0], b[2]-a[1], b[2]+a[1],
            b[0]-a[2], b[0]+a[2])


def source_currents(z):
    if (not isinstance(z, (list, tuple)) or len(z) != 6
            or any(type(v) is not int for v in z) or sum(z) % 2):
        raise ValueError("source current section requires an even-sum integer six-vector")
    h = sum(z)//2
    return (-z[5], h-z[0], h-z[2]-z[3], h-z[1]-z[4]-z[5],
            h-z[2]-z[3]-z[4], z[3])


def geometry(n):
    q, p, values = orbit(n)
    sites = list(product(range(q), repeat=3))
    deltas = [[square(sub(x, y)) for y in values] for x in values]
    candidates = [[j for j in range(q) if sign(sub(scale(q, deltas[i][j]), (1, 0))) <= 0]
                  for i in range(q)]
    neighbors = []
    for x, y, z in sites:
        row = []
        for a in candidates[x]:
            for b in candidates[y]:
                first = add(deltas[x][a], deltas[y][b])
                for c in candidates[z]:
                    if sign(sub(scale(q, add(first, deltas[z][c])), (1, 0))) <= 0:
                        row.append((a*q+b)*q+c)
        neighbors.append(row)
    return q, p, values, sites, deltas, neighbors


def trace(neighbors, layers, *, reverse=False, intervention=None):
    """All reads consume version j of a register written at layer j-1.

    The audit digest is separate metadata, absent from semantic input support.
    Initial register writes have explicit source seeds and no event parent.
    Every positive coefficient makes each declared input influential.
    """
    count = len(neighbors)
    previous_values = []
    chain = "0"*64
    layer_hashes, layer_sums = [], []
    order = list(range(count))
    if reverse:
        order.reverse()
    for j in range(layers+1):
        current = [None]*count
        for i in order:
            reads = [] if j == 0 else [[r, j, [j-1, r], previous_values[r]]
                                       for r in neighbors[i]]
            value = 1+i+(i == intervention) if j == 0 else 1+sum(r[3] for r in reads)
            current[i] = value
            material = [[j, i], reads, [i, j+1, [j, i], value]]
            chain = sha256(bytes.fromhex(chain)+canonical(material)).hexdigest()
        layer_hashes.append(digest(current))
        layer_sums.append(sum(current))
        previous_values = current
    return {"audit_trace_sha256": chain, "layer_value_sha256": layer_hashes,
            "layer_value_sums": layer_sums, "event_count": (layers+1)*count,
            "authenticated_read_count": layers*sum(map(len, neighbors))}


def distances(neighbors, start):
    result = [-1]*len(neighbors)
    result[start] = 0
    queue = [start]
    for i in queue:
        for j in neighbors[i]:
            if result[j] == -1:
                result[j] = result[i]+1
                queue.append(j)
    return result


def root_upper(x, denominator=10**6):
    """Rational upper bound for sqrt(a+b phi), proved by squaring."""
    if sign(x) < 0:
        raise ValueError("negative radicand")
    lo, hi = 0, denominator
    while sign(sub((F(hi, denominator)**2, 0), x)) < 0:
        hi *= 2
    while hi-lo > 1:
        mid = (hi+lo)//2
        if sign(sub((F(mid, denominator)**2, 0), x)) >= 0:
            hi = mid
        else:
            lo = mid
    return F(0) if sign(x) == 0 else F(hi, denominator)


def build_level(n):
    q, p, values, sites, ds, neighbors = geometry(n)
    ordered = sorted(range(q), key=lambda b: (b*p) % q)
    for b in range(q):
        assert values[b][0] == -(b*p//q)
        error = sub(values[b], (F(b*p % q, q), 0))
        assert sign(sub((F(1, q*q), 0), square(error))) > 0
    radii = [sub((1, 0), values[ordered[-1]])]
    radii += [scale(F(1, 2), sub(values[b], values[a])) for a, b in zip(ordered, ordered[1:])]
    radius = radii[0]
    for r in radii[1:]:
        if sign(sub(r, radius)) > 0:
            radius = r
    h_squared = scale(3, square(radius))  # h²/L², not quadrature H²/L².
    ratio_upper = root_upper(scale(q, h_squared))
    inner_speed = max(F(0), 1-2*ratio_upper)
    # Exact nearest-to-half selection, with smallest source label breaking ties.
    center_axis = 0
    for b in range(1, q):
        if sign(sub(square(sub(values[b], (F(1, 2), 0))),
                    square(sub(values[center_axis], (F(1, 2), 0))))) < 0:
            center_axis = b
    center = (center_axis*q+center_axis)*q+center_axis
    layers = isqrt(q)+(isqrt(q)**2 < q)
    records = []
    word_lengths = []
    for b in sites:
        z = source_control([values[i][0] for i in b], b)
        currents = source_currents(z)
        length = sum(map(abs, currents))
        assert sum(z) % 2 == 0 and length <= 27*(q-1)
        records.append([list(b), list(z), list(currents)])
        word_lengths.append(length)
    forward = trace(neighbors, layers)
    reverse = trace(neighbors, layers, reverse=True)
    intervention = trace(neighbors, layers, intervention=center)
    assert forward["layer_value_sha256"] == reverse["layer_value_sha256"]
    reach_diagnostics = []
    for start in sorted({0, center, len(sites)-1}):
        distance = distances(neighbors, start)
        assert min(distance) >= 0
        for k in range(1, layers+1):
            reachable = [i for i, d in enumerate(distance) if d <= k]
            outer, inner, misses = 0, 0, []
            for i, b in enumerate(sites):
                d2 = (0, 0)
                for j in range(3):
                    d2 = add(d2, ds[sites[start][j]][b[j]])
                cone = sign(sub(scale(q, d2), (k*k, 0))) <= 0
                inner_cone = sign(sub(scale(q, d2), (k*k*inner_speed**2, 0))) <= 0
                outer += distance[i] <= k and not cone
                inner += inner_cone and distance[i] > k
                if cone and distance[i] > k:
                    misses.append(i)
            assert outer == inner == 0
            reach_diagnostics.append({"start": start, "layers": k,
                "reachable_count": len(reachable), "reachable_ids_sha256": digest(reachable),
                "outer_cone_violations": outer, "certified_inner_cone_misses": inner,
                "exact_finite_cone_missing_count": len(misses), "missing_ids_sha256": digest(misses)})
    distance = distances(neighbors, center)
    delta = [int(i == center) for i in range(len(sites))]
    response = []
    for k in range(layers+1):
        support = [i for i, v in enumerate(delta) if v]
        assert support == [i for i, d in enumerate(distance) if d <= k]
        assert sum(delta) == intervention["layer_value_sums"][k]-forward["layer_value_sums"][k]
        response.append({"layer": k, "support_count": len(support),
            "support_ids_sha256": digest(support), "positive_integer_delta_sum": sum(delta),
            "positive_integer_delta_maximum": max(delta)})
        delta = [sum(delta[j] for j in row) for row in neighbors]
    intervals = []
    clearance = values[center_axis]
    if sign(sub(sub((1, 0), clearance), clearance)) < 0:
        clearance = sub((1, 0), clearance)
    for k in range(1, layers+1):
        per_layer = [sum(d <= min(j, k-j) for d in distance) for j in range(k+1)]
        inside = sign(sub(scale(4*q, square(clearance)), (k*k, 0))) >= 0
        intervals.append({"layers": k, "inclusive_event_count": sum(per_layer),
            "counts_by_layer": per_layer, "continuum_diamond_inside_cube": inside,
            "count_over_L4_times_sqrt_q": str(F(sum(per_layer), q**3)),
            "diamond_volume_over_pi_L4": str(F(k**4, 24*q*q)),
            "comparison_role": "finite_diagnostic_not_dimension_fit_or_convergence_test"})
    return {"fibonacci_index": n, "q": q, "p": p, "layer_steps": layers,
        "orbit_Qphi": [encode(v) for v in values], "grid_permutation": [(b*p) % q for b in range(q)],
        "record_count": len(sites), "source_records_sha256": digest(records),
        "source_record_examples": [records[i] for i in sorted({0, center, len(sites)-1})],
        "maximum_word_length": max(word_lengths), "sum_word_lengths": sum(word_lengths),
        "word_length_bound": 27*(q-1), "whole_cube_fill_h_squared_over_L2_Qphi": encode(h_squared),
        "one_dimensional_fill_over_L_Qphi": encode(radius),
        "quadrature_displacement_H_squared_over_L2": str(F(12, q*q)),
        "h_over_a_upper": str(ratio_upper), "certified_inner_speed_lower": str(inner_speed),
        "positive_inner_cone": sign(sub((1, 0), scale(4*q, h_squared))) > 0,
        "radius_squared_over_L2": str(F(1, q)), "layer_time_equals_radius": True,
        "neighbors_including_wait_sha256": digest(neighbors),
        "undirected_spatial_edges": (sum(map(len, neighbors))-len(sites))//2,
        "minimum_neighbor_count_including_wait": min(map(len, neighbors)),
        "maximum_neighbor_count_including_wait": max(map(len, neighbors)),
        "exact_width": len(sites), "exact_height_in_events": layers+1,
        "chain_cover": "one wait chain per source record; layers are antichains",
        "forward_execution": forward, "within_layer_reversed_execution": reverse,
        "center_plus_one_intervention": intervention, "intervention_source_id": center,
        "intervention_response": response,
        "reachability_probes": reach_diagnostics, "center_intervals": intervals}


def build():
    return {"schema": SCHEMA, "scope": SCOPE,
        "source_scale_L_squared_Qphi": ["12/5", "-4/5"],
        "window": "[0,L]^3 in the existing source-Gram coordinate witness",
        "source_current_seams": list(SEAMS),
        "event_law": "q_i(0)=i+1; q_i(j)=1+sum(q_r(j-1) for every metric neighbor including i)",
        "event_order": "closure of consecutive-layer authenticated local reads; audit digests are metadata",
        "sampling_measure": "L^3/q^3 per site, layer weight a_q=L/sqrt(q)",
        "finite_schedule": "K_q=ceil(sqrt(q)); T_q=K_q L/sqrt(q) tends to L; four traces alone do not establish the limit",
        "execution_format": "all finite rows executed and hash-chained; full value/read/write trace replayed from templates",
        "levels": [build_level(n) for n in LEVELS],
        "source_pins": {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = canonical(build())
    if args.write:
        OUTPUT.write_bytes(data)
    if args.check and OUTPUT.read_bytes() != data:
        raise ValueError("source-coded causet receipt is stale")
    print("SOURCE_CODED_CAUSETS_REPLAYED", len(data), sha256(data).hexdigest())


if __name__ == "__main__":
    main()

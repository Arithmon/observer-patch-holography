"""Exact causal cones for fixed schedules of the actual native seam primitive.

The schedule and coin word are supplied and shared in each paired comparison.
No adaptive controller or record broadcast is silently treated as a local
operation. This is graph locality in seam-event order, not a speed in metres
per second or a causality theorem for all quantum services.
"""
from __future__ import annotations

import argparse
from collections import deque
import hashlib
from itertools import product
import json
from pathlib import Path

from ..observer import native
from . import slow_geometry as geometry


def edges():
    return tuple(geometry.source().seams())


def _schedule(schedule):
    allowed = set(edges())
    result = []
    for seam in schedule:
        if (not isinstance(seam, (list, tuple)) or len(seam) != 2
                or any(type(port) is not int for port in seam)
                or tuple(sorted(seam)) not in allowed):
            raise ValueError("each scheduled event must be an actual native seam")
        result.append(tuple(seam))
    return tuple(result)


def _port(target):
    if type(target) is not int or not 0 <= target < 12:
        raise ValueError("a native port from zero through eleven is required")
    return target


def distances():
    adjacency = [set() for _ in range(12)]
    for a, b in edges():
        adjacency[a].add(b)
        adjacency[b].add(a)
    result = []
    for origin in range(12):
        row, queue = {origin: 0}, deque([origin])
        while queue:
            port = queue.popleft()
            for neighbor in sorted(adjacency[port]):
                if neighbor not in row:
                    row[neighbor] = row[port] + 1
                    queue.append(neighbor)
        result.append(tuple(row[target] for target in range(12)))
    return tuple(result)


def shortest_path(origin, target):
    origin, target = _port(origin), _port(target)
    distance = distances()
    path = [origin]
    while path[-1] != target:
        current = path[-1]
        neighbors = {b if a == current else a for a, b in edges() if current in (a, b)}
        path.append(min(port for port in neighbors if distance[port][target] == distance[current][target] - 1))
    return tuple(path)


def backward_cone(schedule, target):
    """Sufficient initial-load and local-coin dependencies of one output port."""
    schedule, target = _schedule(schedule), _port(target)
    cone, coin_events = {target}, []
    for index in range(len(schedule) - 1, -1, -1):
        a, b = schedule[index]
        if a in cone or b in cone:
            cone.update((a, b))
            coin_events.append(index)
    return {"initial_ports": sorted(cone), "coin_event_indices": sorted(coin_events)}


def causal_hops(schedule):
    """Minimum seam transmissions along time-ordered dependency paths.

    Waiting at a port is free in this hop count. Crossing an activated seam
    adds one. These are structural upper cones, not claims that every formal
    input can be varied independently within the fixed-total-mass sector.
    """
    dependencies = [{port: 0} for port in range(12)]
    for a, b in _schedule(schedule):
        old_a, old_b = dependencies[a], dependencies[b]
        new_a, new_b = dict(old_a), dict(old_b)
        for origin, hops in old_b.items():
            new_a[origin] = min(new_a.get(origin, hops + 1), hops + 1)
        for origin, hops in old_a.items():
            new_b[origin] = min(new_b.get(origin, hops + 1), hops + 1)
        dependencies[a], dependencies[b] = new_a, new_b
    return tuple(dependencies)


def execute(loads, schedule, coins):
    schedule = _schedule(schedule)
    if (not isinstance(loads, (list, tuple)) or len(loads) != 12
            or any(type(value) is not int or value < 0 for value in loads) or sum(loads) != 2):
        raise ValueError("an actual nonnegative total-two native state is required")
    coins = tuple(coins)
    if len(coins) != len(schedule) or any(type(coin) is not int or coin not in (0, 1) for coin in coins):
        raise ValueError("one fixed binary placement coin per seam event is required")
    state, trace = tuple(loads), [tuple(loads)]
    records = []
    for seam, coin in zip(schedule, coins, strict=True):
        state, record = native.native(state, seam, coin)
        trace.append(state)
        records.append(record)
    return {"final": state, "states": tuple(trace), "records": tuple(records)}


def check_fixed_schedule(schedule, coins):
    """Exhaust all 78 inputs, grouping only states that agree on each cone."""
    schedule, coins = _schedule(schedule), tuple(coins)
    states = native.states()
    outputs = [execute(state, schedule, coins)["final"] for state in states]
    dependencies, metric = causal_hops(schedule), distances()
    rows = []
    for target in range(12):
        cone = backward_cone(schedule, target)
        groups, counts = {}, {}
        for state, output in zip(states, outputs, strict=True):
            key = tuple(state[port] for port in cone["initial_ports"])
            if key in groups and groups[key] != output[target]:
                raise AssertionError("native output depends on a load outside its event cone")
            groups[key] = output[target]
            counts[key] = counts.get(key, 0) + 1
        if set(dependencies[target]) != set(cone["initial_ports"]):
            raise AssertionError("forward and backward dependency constructions disagree")
        if any(metric[origin][target] > hops or hops > len(schedule)
               for origin, hops in dependencies[target].items()):
            raise AssertionError("time-ordered dependency outruns the native graph metric")
        rows.append({"output_port": target, **cone,
                     "minimum_causal_hops": {str(origin): hops for origin, hops in sorted(dependencies[target].items())},
                     "equal_cone_input_pairs_checked": sum(n * (n - 1) // 2 for n in counts.values()),
                     "distinct_cone_input_groups": len(groups)})
    return {"schedule": [list(seam) for seam in schedule], "coins": list(coins),
            "valid_initial_states_checked": len(states), "outputs": rows,
            "native_rule_evaluations": len(states) * len(schedule), "all_pass": True}


def _state(ports):
    loads = [0] * 12
    for port in ports:
        loads[port] += 1
    return tuple(loads)


def delayed_state_witness():
    """Two admissible binary states first influence the receiver after two edges."""
    metric, target = distances(), 0
    first = next(port for port in range(12) if metric[port][target] == 2)
    other = next(port for port in range(12) if port != first and metric[port][target] >= 2)
    path = shortest_path(first, target)
    spectator = next(port for port in range(12) if port not in {*path, other})
    initial = (_state((first, spectator)), _state((other, spectator)))
    schedule = tuple(zip(path, path[1:]))
    forward = [execute(state, schedule, (0, 0)) for state in initial]
    reverse = [execute(state, schedule[::-1], (0, 0)) for state in initial]
    support = [port for port in range(12) if initial[0][port] != initial[1][port]]
    readings = [[state[target] for state in branch["states"]] for branch in forward]
    reversed_readings = [[state[target] for state in branch["states"]] for branch in reverse]
    minimum = min(metric[port][target] for port in support)
    if minimum != 2 or readings[0][:2] != readings[1][:2] or readings[0][2] == readings[1][2]:
        raise AssertionError("two-edge delayed influence witness failed")
    if reversed_readings[0] != reversed_readings[1]:
        raise AssertionError("reversed event order transmitted the witness too soon")
    return {"receiver_port": target, "initial_states": [list(state) for state in initial],
            "initial_difference_support": support, "difference_support_graph_distance": minimum,
            "initial_states_in_binary_total_two_shell": all(max(state) == 1 for state in initial),
            "ordered_path": list(path), "forward_schedule": [list(seam) for seam in schedule],
            "coins_shared_between_states": [0, 0], "receiver_trace": readings,
            "reversed_schedule": [list(seam) for seam in schedule[::-1]],
            "reversed_receiver_trace": reversed_readings,
            "full_forward_load_traces": [[list(state) for state in branch["states"]] for branch in forward],
            "all_pass": True}


def delayed_coin_witness():
    """One local coin perturbation needs two further seams to reach a remote port."""
    metric, target = distances(), 0
    origin = next(port for port in range(12) if metric[port][target] == 3)
    path = shortest_path(origin, target)
    spectator = next(port for port in range(12) if port not in path)
    initial = _state((origin, spectator))
    schedule = tuple(zip(path, path[1:]))
    coins = ((0, 0, 0), (1, 0, 0))
    runs = [execute(initial, schedule, word) for word in coins]
    support = [port for port in range(12) if runs[0]["states"][1][port] != runs[1]["states"][1][port]]
    readings = [[state[target] for state in run["states"]] for run in runs]
    if min(metric[port][target] for port in support) != 2 or readings[0][:3] != readings[1][:3] or readings[0][3] == readings[1][3]:
        raise AssertionError("local coin influence did not respect the two remaining seams")
    return {"receiver_port": target, "initial_state_shared": list(initial),
            "ordered_path": list(path), "schedule_shared": [list(seam) for seam in schedule],
            "coin_words": [list(word) for word in coins], "only_changed_coin_event": 0,
            "post_first_event_difference_support": support,
            "post_first_event_difference_support_graph_distance": 2,
            "receiver_trace": readings,
            "further_ordered_edges_before_receiver_can_differ": 2,
            "full_load_traces": [[list(state) for state in run["states"]] for run in runs],
            "all_pass": True}


def report():
    state_witness, coin_witness = delayed_state_witness(), delayed_coin_witness()
    short = [(), (native.E,), tuple(map(tuple, state_witness["forward_schedule"])),
             tuple(map(tuple, state_witness["reversed_schedule"])),
             tuple(map(tuple, coin_witness["schedule_shared"]))]
    checks = [check_fixed_schedule(schedule, coins) for schedule in short
              for coins in product((0, 1), repeat=len(schedule))]
    # Long fixed sweeps have four declared paired coin words. The general
    # dependency proof covers arbitrary words; these runs do not claim to
    # enumerate all 2^30 words on a sweep.
    for schedule in (edges(), tuple(reversed(edges()))):
        words = ((0,) * len(schedule), (1,) * len(schedule),
                 tuple(i % 2 for i in range(len(schedule))),
                 tuple((i // 3) % 2 for i in range(len(schedule))))
        checks.extend(check_fixed_schedule(schedule, word) for word in words)
    metric = distances()
    root = Path(__file__).resolve().parents[2]
    native_path = Path(native.__file__).resolve()
    return {
        "schema": "oph.common_source.fixed_schedule_native_causality.v1",
        "source": {"native_rule": str(native_path.relative_to(root)),
                   "native_sha256": hashlib.sha256(native_path.read_bytes()).hexdigest(),
                   "incidence": "existing pinned native response carrier, checked by slow_geometry.source"},
        "native_ports": 12, "native_seams": [list(seam) for seam in edges()],
        "graph_distances": [list(row) for row in metric],
        "graph_diameter": max(max(row) for row in metric),
        "fixed_schedule_checks": checks,
        "valid_state_coin_schedule_runs": 78 * len(checks),
        "total_agreeing_cone_pairs_checked": sum(row["equal_cone_input_pairs_checked"] for check in checks for row in check["outputs"]),
        "delayed_initial_state_influence": state_witness,
        "delayed_local_coin_influence": coin_witness,
        "all_pass": all(check["all_pass"] for check in checks),
        "scope": {
            "event_parameter": "one ordered native seam update",
            "graph_unit_speed_bound": "at most one edge per causal layer; the retained witnesses attain it",
            "external_physical_unit_selection_required": False,
            "solver_counter_identified_with_observer_clock": False,
            "sequential_sweep_is_one_parallel_round": False,
            "coins_and_schedule_shared_in_load_comparisons": True,
            "adaptive_remote_controller_or_record_broadcast_included": False,
            "quantum_service_localization_established": False,
            "physical_length_or_time_calibrated": False,
            "speed_of_light_or_Lorentz_causality_derived": False,
            "diameter_three_implies_independent_single_port_variation_in_mass_two_sector": False,
            "meaning": "exact dependency cones and delayed native influence on a fixed supplied finite graph",
        },
    }


def main():
    from ..observer.scalar_response.experiment import canonical, compare_receipt
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = report()
    if args.write:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding="utf-8")), result)
    print(json.dumps({"all_pass": result["all_pass"], "graph_diameter": result["graph_diameter"],
                      "valid_state_coin_schedule_runs": result["valid_state_coin_schedule_runs"],
                      "scope": result["scope"]}, indent=2))


if __name__ == "__main__":
    main()

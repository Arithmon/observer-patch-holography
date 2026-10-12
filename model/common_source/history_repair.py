"""Local constraint repair of one complete, source-generated event history.

The unknowns are unnormalized quantum work vectors at twelve event nodes.
Every repair reads one parent and writes one child. The prepared boundary,
source maps and causal graph are fixed construction inputs. The repair counter
is an auxiliary solver parameter, not a clock inside the represented history.
See HISTORY_REPAIR.md for the finite-DAG theorem and its scope.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from types import MappingProxyType

import numpy as np

from ..observer.source_scalar import compiler as c
from . import experiment as e
from . import scalar_coupling as scalar


BOUNDARY = "boundary"


@dataclass(frozen=True)
class Constraint:
    parent: str
    operator: np.ndarray


@dataclass(frozen=True)
class HistoryConstraints:
    boundary: np.ndarray
    constraints: object
    order: tuple[str, ...]
    weights: object
    depths: object
    dimension: int
    steps: int
    current: int
    feedback_enabled: bool
    coupling: float


def terminal_name(a, b):
    return f"terminal:{a}:{b}"


def build_constraints(d, steps=512, current=0, feedback_enabled=True, coupling=.25):
    """Build only source maps and the fixed preparation, without a solved history."""
    qubits, _ = e._register(d)
    if type(steps) is not int or steps < 1 or type(feedback_enabled) is not bool:
        raise ValueError("positive integer split count and Boolean feedback switch required")
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError("native current -1, 0 or 1 required")
    if isinstance(coupling, (bool, np.bool_)) or not isinstance(coupling, (int, float)) or not math.isfinite(coupling):
        raise ValueError("finite real coupling required")
    initial, _ = scalar.preparation(d)
    kick_word = [c.pulse("X", -current * e.FEEDBACK_ANGLE, [0])]
    kick, _ = c.stabilize_executed_unitary(c.execute(kick_word, qubits))
    boundary = kick @ initial
    boundary.flags.writeable = False
    first, _ = scalar.source_unitary(d, steps, time=1., coupling=coupling)
    second, _ = scalar.source_unitary(d, steps, time=.5, coupling=coupling)
    read_a, _ = e.pointer(d, 0)
    read_b, _ = e.scalar_pointer(d)
    raw = [("first_evolution", BOUNDARY, first)]
    raw += [(f"first_record:{a}", "first_evolution", read_a[a]) for a in (0, 1)]
    raw += [(f"feedback:{a}", f"first_record:{a}", e.feedback(d, a, feedback_enabled)[0]) for a in (0, 1)]
    raw += [(f"second_evolution:{a}", f"feedback:{a}", second) for a in (0, 1)]
    raw += [(terminal_name(a, b), f"second_evolution:{a}", read_b[b]) for a in (0, 1) for b in (0, 1)]
    constraints, depths = {}, {BOUNDARY: 0}
    for name, parent, operator in raw:
        operator = operator.copy()
        operator.flags.writeable = False
        constraints[name] = Constraint(parent, operator)
        depths[name] = depths[parent] + 1
    weights = {name: 1 for name in (BOUNDARY, *constraints)}
    for name in reversed(constraints):
        weights[constraints[name].parent] += weights[name]
    return HistoryConstraints(boundary, MappingProxyType(constraints), tuple(constraints),
                              MappingProxyType(weights), MappingProxyType(depths),
                              d, steps, current, feedback_enabled, float(coupling))


def initialize(system, initial=None):
    """Copy arbitrary inconsistent work vectors while fixing the prepared boundary."""
    if initial is None:
        result = {name: np.zeros_like(system.boundary) for name in system.order}
    else:
        if set(initial) not in (set(system.order), {BOUNDARY, *system.order}):
            raise ValueError("initial assignment must specify every nonboundary node")
        if BOUNDARY in initial and not np.array_equal(initial[BOUNDARY], system.boundary):
            raise ValueError("the prepared history boundary is fixed")
        result = {}
        for name in system.order:
            value = np.asarray(initial[name], dtype=complex)
            if value.shape != system.boundary.shape or not np.isfinite(value).all():
                raise ValueError("finite work vector on the common carrier required")
            result[name] = value.copy()
    result[BOUNDARY] = system.boundary.copy()
    result[BOUNDARY].flags.writeable = False
    return result


def residuals(system, assignment):
    return {name: float(np.linalg.norm(assignment[name] - item.operator @ assignment[item.parent]))
            for name, item in system.constraints.items()}


def potential(system, assignment):
    """Real-valued weighted constraint residual; not an integer tower rank."""
    return math.fsum(system.weights[name] * residual for name, residual in residuals(system, assignment).items())


def repair_once(system, assignment, node):
    """Change exactly one child using its current parent and declared local map."""
    if node not in system.constraints:
        raise ValueError("repair node must be a nonboundary event")
    item = system.constraints[node]
    target = item.operator @ assignment[item.parent]
    residual = float(np.linalg.norm(assignment[node] - target))
    assignment[node] = target
    return residual


def cover_schedule(system, kind="reverse", rounds=5, seed=0):
    """Finite examples of covering schedules; fairness itself is a theorem premise."""
    if kind not in ("forward", "reverse", "shuffled") or type(rounds) is not int or rounds < 1:
        raise ValueError("a declared schedule kind and positive cover count required")
    rng = np.random.default_rng(seed)
    schedule = []
    for _ in range(rounds):
        order = list(system.order)
        if kind == "reverse":
            order.reverse()
        elif kind == "shuffled":
            rng.shuffle(order)
        schedule.extend(order)
    return schedule


def repair(system, schedule, initial=None, tolerance=2e-10):
    """Execute the supplied finite local schedule, even if it is not fair.

    No target terminal state or probability is available to this solver.
    Returned assignments include all intermediate work vectors for inspection;
    audit() emits only compact summaries. Exhausting a bad schedule does not
    certify closure: the returned residual decides that numerical check.
    """
    if isinstance(tolerance, (bool, np.bool_)) or not isinstance(tolerance, (int, float)) or not 0 <= tolerance < math.inf:
        raise ValueError("finite nonnegative residual tolerance required")
    assignment = initialize(system, initial)
    values = [potential(system, assignment)]
    updated, active, largest_defect = 0, 0, 0.
    for node in schedule:
        residual = repair_once(system, assignment, node)
        after = potential(system, assignment)
        # In ideal arithmetic Phi_after <= Phi_before - ||r_node||.
        largest_defect = max(largest_defect, after - values[-1] + residual)
        values.append(after)
        updated += 1
        active += int(residual > tolerance)
    last = residuals(system, assignment)
    return {
        "assignment": assignment,
        "terminal_vectors": {(a, b): assignment[terminal_name(a, b)]
                             for a in (0, 1) for b in (0, 1)},
        "updates": updated,
        "numerically_active_updates": active,
        "initial_potential": values[0],
        "final_potential": values[-1],
        "maximum_constraint_residual": max(last.values()),
        "maximum_potential_drop_defect": largest_defect,
        "closed_within_tolerance": max(last.values()) <= tolerance,
        "potential_trace": values,
    }


def audit(d=4, steps=512, current=0, feedback_enabled=True, coupling=.25):
    """Compare repaired outputs with execution only AFTER solving the constraints."""
    from . import operational
    system = build_constraints(d, steps, current, feedback_enabled, coupling)
    rng = np.random.default_rng(1097)
    initial = {name: rng.normal(size=len(system.boundary)) + 1j * rng.normal(size=len(system.boundary))
               for name in system.order}
    repaired = {kind: repair(system, cover_schedule(system, kind, rounds=1 if kind == "forward" else 5,
                                                  seed=1078), initial)
                for kind in ("forward", "reverse", "shuffled")}
    # Validation-only oracle: neither it nor its terminal vectors enter the
    # graph builder, the initial assignments, or any local repair operation.
    reference = e.history(d, steps, coupling=coupling, feedback_enabled=feedback_enabled, current=current)
    targets = {tuple(row["record"]): np.asarray(row["unnormalized_terminal_state"])
               for row in reference["branches"]}
    summaries, scalar_states = [], []
    for kind, result in repaired.items():
        probabilities, differences, density_differences = [], [], []
        for record, encoded in targets.items():
            state = result["assignment"][terminal_name(*record)]
            target = encoded[:, 0] + 1j * encoded[:, 1]
            probabilities.append({"record": list(record), "probability": float(np.vdot(state, state).real)})
            differences.append(float(np.linalg.norm(state - target)))
            density_differences.append(float(np.linalg.norm(np.outer(state, state.conj()) - np.outer(target, target.conj()), 2)))
        reduced = sum(operational.scalar_reduction(d, np.outer(state, state.conj()))
                      for state in result["terminal_vectors"].values())
        scalar_states.append(reduced)
        eigenvalues = np.linalg.eigvalsh(reduced)
        faithful = bool(eigenvalues[0] > operational.RANK_TOLERANCE)
        extraction = {
            "scalar0_density_real": reduced.real.tolist(),
            "scalar0_density_imaginary": reduced.imag.tolist(),
            "scalar0_density_eigenvalues": [float(value) for value in eigenvalues],
            "numerically_faithful": faithful,
            "numerically_nontrivial_modular_flow": faithful and bool(eigenvalues[-1] - eigenvalues[0] > 2e-10),
            "physical_clock_calibration_claimed": False,
        }
        if d == 2 and faithful:
            extraction["modular_frequency"] = math.log(float(eigenvalues[1] / eigenvalues[0]))
        summaries.append({"schedule": kind, **{key: value for key, value in result.items()
                                               if key not in ("assignment", "terminal_vectors", "potential_trace")},
                          "branches": probabilities,
                          "observer_extraction": extraction,
                          "maximum_terminal_vector_difference_from_source_history": max(differences),
                          "maximum_terminal_density_operator_difference_from_source_history": max(density_differences)})
    maximum_norm = max(float(np.linalg.norm(item.operator, 2)) for item in system.constraints.values())
    extraction_independent = all(np.array_equal(scalar_states[0], state) for state in scalar_states[1:])
    all_pass = maximum_norm <= 1 + 2e-10 and extraction_independent and all(row["closed_within_tolerance"] and row["maximum_potential_drop_defect"] < 2e-10
                   and row["maximum_terminal_vector_difference_from_source_history"] < 2e-10 for row in summaries)
    if not all_pass:
        raise AssertionError("finite history-repair audit failed")
    return {
        "schema": "oph.common_source.finite_history_repair.v1",
        "dimension": d, "steps_per_segment": steps, "native_current": current,
        "feedback_enabled": feedback_enabled, "coupling": float(coupling),
        "event_nodes_including_fixed_boundary": len(system.order) + 1,
        "maximum_causal_depth": max(system.depths.values()),
        "maximum_cover_rounds_for_arbitrary_initial_assignment": 5,
        "constraint_weights": dict(system.weights),
        "maximum_executed_edge_operator_norm": maximum_norm,
        "deterministic_scalar_extraction_schedule_independent": extraction_independent,
        "schedules": summaries, "all_pass": all_pass,
        "scope": {
            "solved_object": "complete unnormalized quantum event history with a fixed prepared boundary",
            "repair_parameter": "auxiliary local update count, not physical time",
            "history_evolution_parameter": "declared Hamiltonian model-time segments 1 and 1/2",
            "unique_schedule_independent_history_for_fixed_inputs": True,
            "one_realized_record_pair_selected_by_repair": False,
            "native_source_execution_of_global_history_repair_proved": False,
            "integer_consensus_tower_witness_claimed": False,
            "graph_action_or_boundary_derived": False,
            "whole_universe_or_empirical_completion": False,
        },
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))

"""Adverse schedules, arbitrary assignments and complete quantum output checks."""
import json

import numpy as np
import pytest

from . import history_repair as h
from . import experiment as e


def inconsistent(system, seed=527):
    rng = np.random.default_rng(seed)
    return {name: rng.normal(size=len(system.boundary)) + 1j * rng.normal(size=len(system.boundary))
            for name in system.order}


def test_graph_is_same_source_history_with_twelve_nodes_and_five_event_layers():
    system = h.build_constraints(2, 16)
    assert len(system.order) + 1 == 12
    assert max(system.depths.values()) == 5
    assert system.weights[h.BOUNDARY] == 12
    for name, item in system.constraints.items():
        assert system.depths[name] == system.depths[item.parent] + 1
        children = [child for child in system.order if system.constraints[child].parent == name]
        assert system.weights[name] == 1 + sum(system.weights[child] for child in children)
        assert np.linalg.norm(item.operator, 2) <= 1 + 3e-13


@pytest.mark.parametrize("d", (2, 4))
@pytest.mark.parametrize("current", (-1, 0, 1))
@pytest.mark.parametrize("enabled", (False, True))
def test_every_terminal_vector_matches_source_history_for_reverse_and_random_fair_schedules(d, current, enabled):
    system = h.build_constraints(d, 16, current, enabled)
    initial = inconsistent(system)
    results = [h.repair(system, h.cover_schedule(system, kind, 5, seed=619), initial)
               for kind in ("reverse", "shuffled")]
    oracle = e.history(d, 16, current=current, feedback_enabled=enabled)
    for result in results:
        assert result["closed_within_tolerance"]
        assert result["maximum_constraint_residual"] < 1e-14
        assert result["maximum_potential_drop_defect"] < 2e-11
        assert np.array_equal(result["assignment"][h.BOUNDARY], system.boundary)
        probability = 0.
        for row in oracle["branches"]:
            encoded = np.asarray(row["unnormalized_terminal_state"])
            expected = encoded[:, 0] + 1j * encoded[:, 1]
            actual = result["terminal_vectors"][tuple(row["record"])]
            np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0)
            probability += np.vdot(actual, actual).real
        assert probability == pytest.approx(1, abs=2e-11)
    for record in results[0]["terminal_vectors"]:
        assert np.array_equal(results[0]["terminal_vectors"][record], results[1]["terminal_vectors"][record])


def test_adverse_reverse_order_needs_five_covers_but_topological_order_needs_one():
    system = h.build_constraints(2, 16)
    initial = inconsistent(system)
    forward = h.repair(system, h.cover_schedule(system, "forward", 1), initial)
    short = h.repair(system, h.cover_schedule(system, "reverse", 4), initial)
    reverse = h.repair(system, h.cover_schedule(system, "reverse", 5), initial)
    assert forward["closed_within_tolerance"] and reverse["closed_within_tolerance"]
    assert not short["closed_within_tolerance"]
    assert forward["updates"] == 11 and reverse["updates"] == 55
    assert forward["numerically_active_updates"] < reverse["numerically_active_updates"]
    for name in system.order:
        assert np.array_equal(forward["assignment"][name], reverse["assignment"][name])


def test_local_repair_changes_only_named_child_and_has_weighted_potential_descent():
    system = h.build_constraints(2, 16)
    assignment = h.initialize(system, inconsistent(system))
    for name in h.cover_schedule(system, "shuffled", 5, seed=1974):
        previous = {node: value.copy() for node, value in assignment.items()}
        before = h.potential(system, assignment)
        residual = h.repair_once(system, assignment, name)
        assert h.potential(system, assignment) <= before - residual + 2e-11
        for other in previous:
            if other != name:
                assert np.array_equal(previous[other], assignment[other])
        item = system.constraints[name]
        assert np.array_equal(assignment[name], item.operator @ assignment[item.parent])


def test_unfair_schedule_does_not_earn_a_false_closure_certificate():
    system = h.build_constraints(2, 16)
    result = h.repair(system, [h.terminal_name(0, 0)] * 55, inconsistent(system))
    assert not result["closed_within_tolerance"]
    assert result["maximum_constraint_residual"] > 1


def test_terminal_oracle_is_unavailable_to_builder_and_solver(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("terminal history oracle was used to construct the solution")
    monkeypatch.setattr(e, "history", forbidden)
    system = h.build_constraints(2, 8, current=-1)
    result = h.repair(system, h.cover_schedule(system, "reverse", 5), inconsistent(system))
    assert result["closed_within_tolerance"]
    assert sum(np.vdot(v, v).real for v in result["terminal_vectors"].values()) == pytest.approx(1, abs=2e-11)


def test_fixed_boundary_cannot_be_repaired_or_replaced():
    system = h.build_constraints(2, 16)
    assignment = h.initialize(system)
    with pytest.raises(ValueError):
        h.repair_once(system, assignment, h.BOUNDARY)
    with pytest.raises(ValueError):
        assignment[h.BOUNDARY][0] = 10
    changed = inconsistent(system)
    changed[h.BOUNDARY] = -system.boundary
    with pytest.raises(ValueError):
        h.initialize(system, changed)


def test_fixed_point_stays_fixed_and_deterministic_state_extraction_is_schedule_independent():
    system = h.build_constraints(2, 16)
    forward = h.repair(system, h.cover_schedule(system, "forward", 1), inconsistent(system))
    reverse = h.repair(system, h.cover_schedule(system, "reverse", 5), inconsistent(system, seed=89))
    repeated = h.repair(system, h.cover_schedule(system, "shuffled", 3), forward["assignment"])
    assert repeated["numerically_active_updates"] == 0
    densities = [sum(np.outer(v, v.conj()) for v in result["terminal_vectors"].values())
                 for result in (forward, reverse, repeated)]
    assert np.array_equal(densities[0], densities[1])
    assert np.array_equal(densities[0], densities[2])


def test_compact_audit_receipt_distinguishes_repair_from_native_hardware_and_physical_time():
    report = h.audit(2, 16, current=1)
    json.dumps(report)
    assert report["all_pass"]
    assert report["deterministic_scalar_extraction_schedule_independent"]
    extractions = [row["observer_extraction"] for row in report["schedules"]]
    assert extractions[0] == extractions[1] == extractions[2]
    assert extractions[0]["numerically_nontrivial_modular_flow"]
    assert extractions[0]["modular_frequency"] > 0
    assert not extractions[0]["physical_clock_calibration_claimed"]
    assert report["scope"]["unique_schedule_independent_history_for_fixed_inputs"]
    assert not report["scope"]["native_source_execution_of_global_history_repair_proved"]
    assert not report["scope"]["integer_consensus_tower_witness_claimed"]
    assert not report["scope"]["whole_universe_or_empirical_completion"]
    assert all("assignment" not in row and "terminal_vectors" not in row for row in report["schedules"])


@pytest.mark.parametrize("kwargs", ({"current": True}, {"current": 2}, {"steps": 0},
                                    {"feedback_enabled": 1}, {"coupling": float("nan")}))
def test_rejects_undeclared_controls(kwargs):
    with pytest.raises(ValueError):
        h.build_constraints(2, **kwargs)

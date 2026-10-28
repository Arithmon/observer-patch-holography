"""Independent endpoint, first-exit, quantum-effect and receipt controls."""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import pytest

from . import localized as loc
from ..observer import adaptive_record_recovery as adaptive
from ..observer import record_recovery as records


@pytest.fixture(scope="module")
def report():
    return loc.run()


def endpoint(state, edge, coin, law):
    """Independent integer balancing and raw record formula."""
    a, b = edge
    low, parity = divmod(state[a] + state[b], 2)
    out = list(state)
    out[a], out[b] = low + coin * parity, low + (1 - coin) * parity
    current = state[a] - out[a]
    raw = {0: 0, 1: 1, -1: 2}[current]
    return tuple(out), current, raw + (4 * parity if law == "E" else 0)


@pytest.mark.parametrize("law", ["C", "E"])
def test_complete_quantum_maps_and_population_kernel_match_independent_endpoints(law):
    states = loc.basis()
    index = {x: i for i, x in enumerate(states)}
    expected = [[0] * 66 for _ in states]
    maps = [row for row in adaptive.branch_maps(binary=True)
            if adaptive.actions()[row.action][0] == law]
    for row in maps:
        _, edge = adaptive.actions()[row.action]
        pairs = []
        for i, x in enumerate(states):
            y, _, raw = endpoint(x, edge, row.coin, law)
            if raw == row.raw:
                pairs.append((i, index[y]))
                expected[i][index[y]] += 1
        assert pairs == list(zip(row.inputs, row.outputs, strict=True))
    assert tuple(map(tuple, expected)) == loc.transition_counts(law)
    assert all(sum(row) == 60 for row in expected)
    assert sum(not row.inputs for row in maps) == (120 if law == "C" else 300)
    assert all(loc.movement_activity(x) == F(60 - expected[i][i], 60)
               for i, x in enumerate(states))


@pytest.mark.parametrize("law", ["C", "E"])
def test_exit_effect_has_no_coherence_cross_terms_including_entangled_inputs(law):
    # Explicit K* P_out K on all matrix units, using the independent endpoints.
    states, cluster = loc.basis(), set(loc.cluster_indices())
    outside = {states[i] for i in range(66) if i not in cluster}
    effect = [[0] * 66 for _ in range(66)]
    for edge in records.seams():
        for coin in (0, 1):
            images = [endpoint(x, edge, coin, law) for x in states]
            for i in range(66):
                for j in range(66):
                    yi, _, ri = images[i]
                    yj, _, rj = images[j]
                    effect[i][j] += int(ri == rj and yi == yj and yi in outside)
    assert all(effect[i][j] == 0 for i in range(66) for j in range(66) if i != j)
    assert all(effect[i][i] == 4 for i in cluster)
    # Therefore Tr((P_out tensor I)Gamma(rho_AR))=1/15 for EVERY rho_AR
    # supported on A tensor R, not only a diagonal preparation.


def test_killed_chain_transport_return_and_stationary_occupancy(report):
    uniform = report["uniform"]
    assert uniform["distance_class_sizes"] == {"1": 30, "2": 30, "3": 6}
    assert uniform["distance_transition_counts"] == [[56, 4, 0], [4, 54, 2], [0, 10, 50]]
    killed = uniform["killed_transition_counts"]
    assert all(sum(row) == 56 and row[i] == 52 for i, row in enumerate(killed))
    assert all(sum(value != 0 for j, value in enumerate(row) if j != i) == 4
               for i, row in enumerate(killed))
    # Independent finite-step multiplication verifies the geometric loss.
    distribution = [F(int(i == 0)) for i in range(30)]
    for n in range(8):
        assert sum(distribution) == F(14, 15) ** n
        distribution = [sum(distribution[i] * F(killed[i][j], 60) for i in range(30))
                        for j in range(30)]
    assert F(uniform["mean_first_exit_attempts"]) == 15
    assert F(uniform["mean_return_attempts_after_first_exit"]) == 18
    assert F(15, 15 + 18) == F(uniform["stationary_cluster_mass"])
    assert F(uniform["mean_within_cluster_moves_before_exit"]) == 1


def literal_history(law, epsilon, horizon):
    initial = tuple(int(i in (0, 1)) for i in range(12))
    edges = records.seams()
    paths = [(initial, 0, False, F(1))]
    for _ in range(horizon):
        following = []
        for state, memory, exited, weight in paths:
            preferred = (0, 5) if memory else (0, 1)
            for edge in edges:
                probability = epsilon / 30 + (1 - epsilon if edge == preferred else 0)
                for coin in (0, 1):
                    out, current, _ = endpoint(state, edge, coin, law)
                    lost = tuple(i for i, v in enumerate(out) if v) not in edges
                    following.append((out, int(current > 0), exited or lost, weight * probability / 2))
        paths = following
    return {"survival": sum(w for _, _, lost, w in paths if not lost),
            "currently_clustered": sum(w for x, _, _, w in paths
                                       if tuple(i for i, v in enumerate(x) if v) in edges),
            "total": sum(w for _, _, _, w in paths)}


@pytest.mark.parametrize("law", ["C", "E"])
def test_feedback_killed_population_dp_matches_literal_unselected_history(law):
    epsilon = F(1, 16)
    literal = literal_history(law, epsilon, 2)
    actual = loc.feedback_history(law, epsilon, (2,))[0]
    assert literal["total"] == 1
    assert literal["survival"] == actual["survival"]
    assert literal["currently_clustered"] == actual["currently_clustered"]
    assert actual["returned_and_currently_clustered"] > 0
    assert actual["currently_clustered"] + actual["currently_separated"] == 1
    assert actual["survival"] + actual["ever_exited"] == 1


def test_feedback_floor_bound_and_ineligible_freeze_controls(report):
    for epsilon in (F(1, 16), F(1, 2), F(1)):
        bound = loc.policy_bounds(epsilon)
        rows, denominator = loc.feedback_counts("C", epsilon)
        cluster = {2 * i + m for i in loc.cluster_indices() for m in (0, 1)}
        assert all(F(sum(row[j] for j in range(132) if j not in cluster), denominator)
                   >= epsilon / 15 for i, row in enumerate(rows) if i in cluster)
        actual = loc.feedback_history("C", epsilon, (1, 8))
        assert all(row["survival"] <= bound["survival_factor_upper"] ** row["attempts"]
                   for row in actual)
    assert loc.policy_bounds(F(1, 16))["mean_first_exit_attempts_upper"] == 240
    assert loc.policy_bounds(F(1, 16), F(1, 16))["mean_first_exit_attempts_upper"] == 3840
    frozen = tuple(int(i in (0, 1)) for i in range(12))
    assert all(endpoint(frozen, (0, 1), c, "C")[0] == frozen for c in (0, 1))
    assert report["negative_controls"]["zero_exploration_same_occupied_seam"]["outside_positive_exploration"]
    with pytest.raises(ValueError):
        loc.policy_bounds(F(0))
    with pytest.raises(ValueError):
        loc.policy_bounds(F(1), F(0))


def test_equal_populations_do_not_supply_equal_quantum_records_or_memory(report):
    c, e = (report["complete_C_E_instruments"][law] for law in ("C", "E"))
    assert c["raw_maps"] != e["raw_maps"]
    assert F(c["maximum_cluster_pair_agreement"]) == F(5, 6)
    assert F(e["maximum_cluster_pair_agreement"]) == F(7, 10)
    assert report["feedback"]["histories"]["C"] == report["feedback"]["histories"]["E"]
    bounds = report["quantum_memory_comparison"]["bounds_at_512_actual_repairs"]
    assert F(bounds["C"]["encoded_final_recovery_trace_error_lower"]) > F(34, 100)
    assert F(bounds["E"]["encoded_final_recovery_trace_error_lower"]) > F(49, 100)
    assert not report["quantum_memory_comparison"]["finite_time_entanglement_breaking_inferred"]


def test_reference_exterior_failure_survives_uniform_location_mixture(report):
    exterior = report["exterior"]
    assert F(exterior["fixed_core_exterior_trace_distance"]) == F(65, 66)
    assert F(exterior["uniform_cluster_mixture_exterior_trace_distance"]) == F(16, 33)
    assert F(exterior["fixed_core_exterior_distance_from_shared_reference"]) == F(43, 45)
    assert not exterior["reference_like_exterior_established"]
    assert not exterior["new_localization_detector_implemented"]


@pytest.mark.parametrize("bad", [True, 0.5, F(-1), F(2)])
def test_invalid_probabilities_are_not_coerced(bad):
    with pytest.raises(ValueError):
        loc.policy_bounds(bad)


def test_retained_exact_receipt_and_mutation_rejection(report, tmp_path):
    loc.contract.check_receipt(loc.OUTPUT, report)
    assert report["campaign"]["predecessor_git_commit"] == "3d282b3d59763fa1a5171ad9a4fbc09fc0abf35a"
    assert not report["scope"]["full_CP03_accepted"]
    assert not report["scope"]["mass_or_energy_derived"]
    changed = json.loads(loc.OUTPUT.read_text())
    changed["uniform"]["mean_first_exit_attempts"] = "150"
    path = tmp_path / "mutated.json"
    path.write_text(json.dumps(changed))
    process = subprocess.run([sys.executable, "-m", "model.clock_defect_energy.localized", "--check", str(path)],
                             cwd=Path(loc.__file__).resolve().parents[2], capture_output=True, text=True)
    assert process.returncode != 0
    assert "diagnostic receipt differs from exact replay" in process.stderr


@pytest.mark.parametrize("replacement", [False, 0.0])
def test_receipt_rejects_type_changes_equal_under_python_equality(report, tmp_path, replacement):
    changed = json.loads(loc.OUTPUT.read_text())
    changed["feedback"]["initial_memory"] = replacement
    assert changed == report  # Why plain dict equality is insufficient.
    path = tmp_path / "type_change.json"
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="differs from exact replay"):
        loc.contract.check_receipt(path, report)


def test_receipt_rejects_duplicate_keys_even_with_identical_values(report, tmp_path):
    text = json.dumps(report)
    needle = '"full_CP03_accepted": false'
    assert needle in text
    path = tmp_path / "duplicate.json"
    path.write_text(text.replace(needle, needle + ", " + needle))
    with pytest.raises(ValueError, match="duplicate JSON key"):
        loc.contract.check_receipt(path, report)

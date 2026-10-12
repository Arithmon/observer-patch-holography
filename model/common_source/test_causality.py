"""Native locality, admissible delayed signals and fixed schedule receipts."""
from itertools import product
import json

import pytest

from ..observer import native
from . import causality as c


@pytest.fixture(scope="module")
def receipt():
    return c.report()


def test_actual_carrier_graph_and_shortest_paths():
    metric = c.distances()
    assert len(c.edges()) == 30
    assert max(map(max, metric)) == 3
    for start in range(12):
        assert sorted(metric[start]) == [0] + [1] * 5 + [2] * 5 + [3]
        for end in range(12):
            path = c.shortest_path(start, end)
            assert path[0] == start and path[-1] == end
            assert len(path) - 1 == metric[start][end]
            assert all(tuple(sorted(edge)) in c.edges() for edge in zip(path, path[1:]))


def test_all_initial_states_and_declared_coin_words_checked(receipt):
    assert receipt["valid_state_coin_schedule_runs"] == 78 * 27
    assert len(receipt["fixed_schedule_checks"]) == 27
    assert receipt["total_agreeing_cone_pairs_checked"] > 100000
    assert all(row["valid_initial_states_checked"] == 78 for row in receipt["fixed_schedule_checks"])
    assert receipt["all_pass"]
    json.dumps(receipt, allow_nan=False)


def test_primitive_is_local_conservative_and_cone_coin_dependencies_are_sufficient():
    schedule = ((3, 2), (2, 10), (10, 0), (0, 1))
    for state in native.states():
        for coins in product((0, 1), repeat=len(schedule)):
            run = c.execute(state, schedule, coins)
            for before, after, (a, b), coin in zip(run["states"], run["states"][1:], schedule, coins):
                assert sum(after) == 2 and min(after) >= 0
                assert all(before[port] == after[port] for port in range(12) if port not in (a, b))
                total = before[a] + before[b]
                assert after[a] == (total + coin) // 2
                assert after[b] == total - after[a]
            for target in range(12):
                cone = c.backward_cone(schedule, target)
                for index in set(range(len(schedule))) - set(cone["coin_event_indices"]):
                    changed = list(coins)
                    changed[index] ^= 1
                    assert c.execute(state, schedule, changed)["final"][target] == run["final"][target]


def test_actual_binary_state_signal_saturates_two_edge_lower_bound(receipt):
    witness = receipt["delayed_initial_state_influence"]
    assert witness["initial_difference_support"] == [2, 3]
    assert witness["difference_support_graph_distance"] == 2
    assert witness["initial_states_in_binary_total_two_shell"]
    assert all(tuple(state) in native.states(binary=True) for state in witness["initial_states"])
    assert witness["receiver_trace"] == [[0, 0, 1], [0, 0, 0]]
    assert witness["reversed_receiver_trace"] == [[0, 0, 0], [0, 0, 0]]
    forward = c.backward_cone(witness["forward_schedule"], 0)["initial_ports"]
    reverse = c.backward_cone(witness["reversed_schedule"], 0)["initial_ports"]
    assert 2 in forward and 2 not in reverse


def test_local_coin_signal_requires_two_further_seams(receipt):
    witness = receipt["delayed_local_coin_influence"]
    assert witness["ordered_path"] == [3, 2, 10, 0]
    assert witness["post_first_event_difference_support"] == [2, 3]
    assert witness["post_first_event_difference_support_graph_distance"] == 2
    assert witness["receiver_trace"] == [[0, 0, 0, 1], [0, 0, 0, 0]]
    assert sum(a != b for a, b in zip(*witness["coin_words"])) == 1


def test_mass_constraint_precludes_independent_single_port_perturbation():
    states = native.states()
    assert all(sum(a != b for a, b in zip(first, second)) >= 2
               for index, first in enumerate(states) for second in states[index + 1:])


def test_scope_does_not_identify_a_sweep_or_host_counter_with_observer_time(receipt):
    scope = receipt["scope"]
    assert not scope["sequential_sweep_is_one_parallel_round"]
    assert not scope["solver_counter_identified_with_observer_clock"]
    assert not scope["external_physical_unit_selection_required"]
    assert not scope["adaptive_remote_controller_or_record_broadcast_included"]
    assert not scope["quantum_service_localization_established"]
    assert not scope["speed_of_light_or_Lorentz_causality_derived"]


@pytest.mark.parametrize("schedule", (((0, 0),), ((0, 3),), ((True, 0),), ((0, 12),)))
def test_nonseam_schedules_rejected(schedule):
    with pytest.raises(ValueError):
        c.backward_cone(schedule, 0)


@pytest.mark.parametrize("coins", ((), (2,), (True,), (0, 1)))
def test_wrong_coin_word_rejected(coins):
    with pytest.raises(ValueError):
        c.execute(native.states()[0], (native.E,), coins)


def test_receipt_cli_round_trip_and_tampering_rejected(tmp_path, monkeypatch, capsys):
    path = tmp_path / "causality.json"
    monkeypatch.setattr("sys.argv", ["causality", "--write", str(path)])
    c.main()
    monkeypatch.setattr("sys.argv", ["causality", "--check", str(path)])
    c.main()
    result = json.loads(path.read_text())
    result["delayed_initial_state_influence"]["receiver_trace"][0][-1] = 0
    path.write_text(json.dumps(result))
    with pytest.raises(ValueError, match="Receipt mismatch"):
        c.main()
    capsys.readouterr()


def test_retained_native_rule_path_is_platform_independent(receipt):
    assert receipt["source"]["native_rule"] == "model/observer/native.py"

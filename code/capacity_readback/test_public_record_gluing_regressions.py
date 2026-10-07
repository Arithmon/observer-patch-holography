"""Reproductions against main: invalid or merged public record sections."""
from copy import deepcopy

import pytest

from correctable_public_record_capacity import evaluate_terminal, section_id
from public_record_csp import public_global_sections_csp
from test_correctable_public_record_capacity import binary_packet


def loop(right):
    return {"left_observer": "a", "right_observer": "a",
            "left_readout": {"0": "x", "1": "y"}, "right_readout": right}


def test_inconsistent_self_interface_has_no_global_sections():
    assert public_global_sections_csp({"a": ["0", "1"]},
                                     [loop({"0": "y", "1": "x"})]) == []


def test_self_interface_selects_only_its_equalizer():
    assert public_global_sections_csp({"a": ["0", "1"]},
                                     [loop({"0": "x", "1": "x"})]) == [{"a": "0"}]


def test_self_interface_cannot_manufacture_a_saturated_packet():
    packet = binary_packet()
    packet["interfaces"].append({
        "left_observer": "alice", "right_observer": "alice",
        "left_readout": {"0": "x", "1": "y"},
        "right_readout": {"0": "y", "1": "x"},
    })
    assert evaluate_terminal(packet)["status"] == "NO_PUBLIC_RECORD_REACHABILITY"


def test_distinct_sections_never_share_a_record_identifier():
    first = {"a": "x|b=y", "b": "z"}
    second = {"a": "x", "b": "y|b=z"}
    assert first != second
    assert section_id(first) != section_id(second)


def test_readouts_must_land_in_the_declared_interface_atoms():
    seam = loop({"0": "x", "1": "y"})
    seam["interface_atoms"] = ["elsewhere"]
    with pytest.raises(ValueError):
        public_global_sections_csp({"a": ["0", "1"]}, [seam])


@pytest.mark.parametrize("observers", [{"a": "01"}, {"a": [0, 1]}, {1: ["x"]}, {"a": [""]}])
def test_atom_and_observer_labels_are_explicit_strings(observers):
    with pytest.raises(ValueError):
        public_global_sections_csp(observers, [])


@pytest.mark.parametrize("value", [True, 0, None])
def test_readout_values_are_atom_labels_not_coercible_values(value):
    seam = loop({"0": value, "1": value})
    seam["left_readout"] = deepcopy(seam["right_readout"])
    with pytest.raises(ValueError):
        public_global_sections_csp({"a": ["0", "1"]}, [seam])


def test_long_singleton_diagram_does_not_depend_on_python_recursion_limit():
    observers = {f"o{i}": ["x"] for i in range(1200)}
    assert public_global_sections_csp(observers, []) == [{x: "x" for x in sorted(observers)}]


def test_unreachable_record_is_a_packet_failure_not_an_uncaught_exception():
    packet = binary_packet()
    packet["reachability_witnesses"]["not-a-public-section"] = ["write"]
    assert evaluate_terminal(packet)["status"] == "NO_PUBLIC_RECORD_REACHABILITY"


@pytest.mark.parametrize("bad", [42, None, ("a",)])
def test_invalid_reachability_identifiers_cannot_crash_packet_classification(bad):
    packet = binary_packet()
    packet["reachability_witnesses"].update({"unknown": ["write"], bad: ["write"]})
    assert evaluate_terminal(packet)["status"] == "NO_PUBLIC_RECORD_REACHABILITY"

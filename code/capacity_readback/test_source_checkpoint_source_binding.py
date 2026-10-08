"""Original-input controls for the fixed source's named checkpoint action.

Capacity and even an abstract multiplication table cannot identify the physical
operation denoted by a continuation name.  The scalar oracle below uses the
declared port coordinates, independently of the producer's action helpers.
"""
from __future__ import annotations

import copy
import hashlib
import json

import pytest

from source_derived_public_checkpoint_packet import (
    build_source_derived_packet,
    certify_source_derived_packet,
)


@pytest.fixture(scope="module")
def original_packet():
    return build_source_derived_packet()


def rehash(packet):
    packet.pop("packet_sha256", None)
    packet["packet_sha256"] = hashlib.sha256(
        json.dumps(packet, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def source_action(name, slot):
    """Independent coordinate action: reflection, rotation, antipode, flip."""
    rotation, reflected, antipodal, flipped = (
        int(part[1:]) for part in name.split("_")
    )
    port, orientation = slot.split("/")
    if port in ("north", "south"):
        if antipodal:
            port = "south" if port == "north" else "north"
    else:
        ring, index_text = port.split("_")
        index = int(index_text)
        if reflected:
            index = -index - (ring == "lower")
        index += rotation
        if antipodal:
            index += 2 if ring == "upper" else 3
            ring = "lower" if ring == "upper" else "upper"
        port = f"{ring}_{index % 5}"
    if flipped:
        orientation = "check" if orientation == "write" else "write"
    return f"{port}/{orientation}"


def supplied_actions(packet):
    return {
        channel["continuation_id"]: {
            source: next(output for output, weight in row.items() if weight == 1)
            for source, row in channel["rows"].items()
        }
        for channel in packet["global_checkpoint_kernels"]
    }


def composition_failures(packet):
    """Count exact pointwise failures from the supplied rows, without helpers."""
    actions = supplied_actions(packet)
    return sum(
        actions[left][actions[right][source]] != actions[result][source]
        for left, rights in packet["support_relation_composition"]["table"].items()
        for right, result in rights.items()
        for source in actions[left]
    )


def replace_actions(packet, actions):
    """Preserve agreement of the public and every local supplied checkpoint."""
    aliases = packet["public_section_aliases"]
    slot_by_id = {sid: slot for slot, sid in aliases.items()}
    for channel in packet["global_checkpoint_kernels"]:
        name = channel["continuation_id"]
        channel["rows"] = {source: {target: 1.0} for source, target in actions[name].items()}
        packet["local_checkpoint_packets"][name] = {
            observer: {
                slot: {f"{observer}::record::{slot_by_id[actions[name][sid]]}": 1.0}
                for slot, sid in aliases.items()
            }
            for observer in packet["observers"]
        }


def require_refusal(packet):
    rehash(packet)
    receipt = certify_source_derived_packet(packet)
    assert receipt["status"] != "PASS", "false fixed-source checkpoint certificate"


@pytest.mark.parametrize("presentation", ["canonical", "reordered", "explicit_zeros"])
def test_valid_source_presentations_are_certified(original_packet, presentation):
    packet = copy.deepcopy(original_packet)
    aliases = packet["public_section_aliases"]
    slots = {sid: slot for slot, sid in aliases.items()}
    for name, action in supplied_actions(packet).items():
        assert all(slots[target] == source_action(name, slots[source])
                   for source, target in action.items())
    assert composition_failures(packet) == 0
    if presentation == "reordered":
        packet["global_checkpoint_kernels"].reverse()
        packet["interfaces"].reverse()
        packet["public_global_sections"].reverse()
        for channel in packet["global_checkpoint_kernels"]:
            channel["rows"] = dict(reversed(list(channel["rows"].items())))
        packet["local_checkpoint_packets"] = dict(
            reversed(list(packet["local_checkpoint_packets"].items()))
        )
        packet["public_section_aliases"] = dict(reversed(list(aliases.items())))
    elif presentation == "explicit_zeros":
        channel = packet["global_checkpoint_kernels"][0]
        source, row = next(iter(channel["rows"].items()))
        target = next(iter(row))
        zero_target = next(sid for sid in slots if sid != target)
        row[zero_target] = 0
        for observer in packet["observers"]:
            packet["local_checkpoint_packets"][channel["continuation_id"]][observer][slots[source]][
                f"{observer}::record::{slots[zero_target]}"
            ] = 0
    rehash(packet)
    receipt = certify_source_derived_packet(packet)
    assert receipt["status"] == "PASS"
    assert receipt["capacity"]["exact_zero_error_capacity"] == 24
    assert receipt["checkpoint_family"]["composition"]["composition_entries_checked"] == 1600


@pytest.mark.parametrize("replacement", ["swap_identity", "inverse_all", "source_automorphism"])
def test_capacity_and_group_properties_do_not_certify_named_source(original_packet, replacement):
    packet = copy.deepcopy(original_packet)
    original = supplied_actions(packet)
    if replacement == "swap_identity":
        actions = dict(original)
        identity, rotation = "r0_s0_a0_f0", "r1_s0_a0_f0"
        actions[identity], actions[rotation] = actions[rotation], actions[identity]
        expected_failed_points = 4620
    elif replacement == "inverse_all":
        actions = {name: {target: source for source, target in action.items()}
                   for name, action in original.items()}
        expected_failed_points = 19200
    else:
        actions = {name: original[f"r{(-int(name[1])) % 5}{name[2:]}"] for name in original}
        expected_failed_points = 0
    replace_actions(packet, actions)
    assert composition_failures(packet) == expected_failed_points
    aliases = packet["public_section_aliases"]
    if replacement == "source_automorphism":
        # Every supplied composition holds, yet the named source rotation is wrong.
        assert source_action("r1_s0_a0_f0", "upper_0/write") == "upper_1/write"
        assert actions["r1_s0_a0_f0"][aliases["upper_0/write"]] == aliases["upper_4/write"]
    elif replacement == "inverse_all":
        assert all(source == target for source, target in actions["r0_s0_a0_f0"].items())
    require_refusal(packet)


@pytest.mark.parametrize("change", ["missing", "duplicate", "unknown"])
def test_named_continuation_census_is_bound_to_the_source(original_packet, change):
    packet = copy.deepcopy(original_packet)
    channels = packet["global_checkpoint_kernels"]
    if change == "missing":
        channels[:] = [channel for channel in channels
                       if channel["continuation_id"] != "r2_s0_a0_f0"]
    else:
        extra = copy.deepcopy(channels[0])
        if change == "unknown":
            name = extra["continuation_id"]
            extra["continuation_id"] = "unbound_source_operation"
            packet["local_checkpoint_packets"][extra["continuation_id"]] = copy.deepcopy(
                packet["local_checkpoint_packets"][name]
            )
        channels.append(extra)
    require_refusal(packet)


def test_source_aliases_are_bound_to_actual_section_atoms(original_packet):
    packet = copy.deepcopy(original_packet)
    aliases = packet["public_section_aliases"]
    aliases["upper_0/write"], aliases["upper_1/write"] = (
        aliases["upper_1/write"], aliases["upper_0/write"]
    )
    require_refusal(packet)


@pytest.mark.parametrize("extra_domain", ["continuation", "observer", "source"])
def test_complete_local_read_manifest_has_no_unbound_domain(original_packet, extra_domain):
    packet = copy.deepcopy(original_packet)
    local = packet["local_checkpoint_packets"]
    name = "r1_s0_a0_f0"
    if extra_domain == "continuation":
        local["unbound_source_operation"] = copy.deepcopy(local[name])
    elif extra_domain == "observer":
        local[name]["unbound_observer"] = copy.deepcopy(local[name]["north"])
    else:
        local[name]["north"]["unbound_source_slot"] = copy.deepcopy(local[name]["north"]["north/write"])
    require_refusal(packet)


@pytest.mark.parametrize("weight", [True, "1", "0.99999999999999999999"])
def test_local_probabilities_are_not_coerced_to_a_different_input(original_packet, weight):
    packet = copy.deepcopy(original_packet)
    row = packet["local_checkpoint_packets"]["r1_s0_a0_f0"]["north"]["north/write"]
    row[next(iter(row))] = weight
    require_refusal(packet)

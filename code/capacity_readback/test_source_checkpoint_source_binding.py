"""Original-input controls for the fixed source's named checkpoint action.

Capacity and even an abstract multiplication table cannot identify the physical
operation denoted by a continuation name.  The scalar oracle below uses the
declared port coordinates, independently of the producer's action helpers.
"""
from __future__ import annotations

import copy
import hashlib
import json
from decimal import Decimal
from fractions import Fraction

import pytest
import source_derived_public_checkpoint_packet as source_packet

from source_derived_public_checkpoint_packet import (
    build_source_derived_packet,
    certify_source_derived_packet,
    verify_local_marginal_consistency,
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


@pytest.mark.parametrize("presentation", [
    "canonical", "reordered", "explicit_zeros", "near_unit_weights", "equivalent_public_cuts",
])
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
        packet["observers"] = dict(reversed(list(packet["observers"].items())))
        packet["publicness_policy"][0].reverse()
        packet["interfaces"].reverse()
        packet["public_global_sections"].reverse()
        for channel in packet["global_checkpoint_kernels"]:
            channel["authorized_observers"].reverse()
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
    elif presentation == "near_unit_weights":
        # The shared global-kernel contract treats near-unit masses as relative
        # weights; preserve exact normalization instead of tightening that domain.
        row = next(iter(packet["global_checkpoint_kernels"][0]["rows"].values()))
        row[next(iter(row))] = 1 - 2**-50
    elif presentation == "equivalent_public_cuts":
        # The existing generic policy is a family of sets, not ordered lists.
        packet["publicness_policy"].append(list(reversed(packet["publicness_policy"][0])))
        packet["publicness_policy"][0].append("north")
        for channel in packet["global_checkpoint_kernels"]:
            channel["authorized_observers"].append("north")
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


@pytest.mark.parametrize("weight", [Fraction(10**20 - 1, 10**20), Decimal("0.99999999999999999999")])
def test_exact_local_deficits_survive_float_rounding(original_packet, weight):
    packet = copy.deepcopy(original_packet)
    row = packet["local_checkpoint_packets"]["r1_s0_a0_f0"]["north"]["north/write"]
    row[next(iter(row))] = weight
    assert float(weight) == 1.0 and weight != 1
    # This public checker accepts native exact scalars; they do not need an
    # unsupported Fraction/Decimal representation in the JSON packet hash.
    receipt = verify_local_marginal_consistency(packet)
    assert receipt["status"] == "LOCAL_MARGINAL_MISMATCH"
    json.dumps(receipt, allow_nan=False)


@pytest.mark.parametrize("weight", [1, 1.0, Fraction(1), Decimal(1)])
def test_exact_numeric_local_units_remain_valid(original_packet, weight):
    packet = copy.deepcopy(original_packet)
    row = packet["local_checkpoint_packets"]["r1_s0_a0_f0"]["north"]["north/write"]
    row[next(iter(row))] = weight
    assert verify_local_marginal_consistency(packet)["status"] == "PASS"


@pytest.mark.parametrize("field,value", [
    ("complete", False),
    ("complete", "true"),
    ("composition_order", "right-after-left"),
])
def test_composition_claim_is_bound_to_its_declared_semantics(original_packet, field, value):
    packet = copy.deepcopy(original_packet)
    packet["support_relation_composition"][field] = value
    require_refusal(packet)


@pytest.mark.parametrize("field,value", [
    ("continuation_family_kind", "C40"),
    ("continuation_manifest_complete", "true"),
    ("local_marginal_manifest_complete", False),
])
def test_source_family_claim_is_not_unchecked_metadata(original_packet, field, value):
    packet = copy.deepcopy(original_packet)
    packet[field] = value
    require_refusal(packet)


def test_composition_is_replayed_without_trusting_its_producer(original_packet, monkeypatch):
    packet = copy.deepcopy(original_packet)
    table = packet["support_relation_composition"]["table"]
    # Swap chronology on one genuinely noncommuting pair, keeping the source
    # kernels and complete table domain intact. A shared table-construction bug
    # must not make its own faulty result into the verifier's expected result.
    reflection, rotation = "r0_s1_a0_f0", "r1_s0_a0_f0"
    table[reflection][rotation] = table[rotation][reflection]
    assert composition_failures(packet) == 20
    monkeypatch.setattr(source_packet, "_continuation_composition_table", lambda _: table)
    assert source_packet._verify_composition(packet)["status"] == "COMPOSITION_TABLE_MISMATCH"
    require_refusal(packet)


@pytest.mark.parametrize("omission", ["empty", "outer_row", "inner_cell", "unknown_result"])
def test_every_claimed_composition_has_a_supplied_named_result(original_packet, omission):
    packet = copy.deepcopy(original_packet)
    table = packet["support_relation_composition"]["table"]
    identity, rotation = "r0_s0_a0_f0", "r1_s0_a0_f0"
    assert len(table) == 40 and all(len(row) == 40 for row in table.values())
    if omission == "empty":
        table.clear()
    elif omission == "outer_row":
        del table[rotation]
    elif omission == "inner_cell":
        del table[identity][rotation]
    else:
        table[identity][rotation] = "unbound_composition_result"
    # All source kernels, local reads and packet hash remain valid. A partial
    # or empty replay must not claim that all 1,600 compositions were checked.
    rehash(packet)
    receipt = certify_source_derived_packet(packet)
    assert receipt["status"] == "COMPOSITION_TABLE_MISMATCH"
    json.dumps(receipt, allow_nan=False)


def test_explicit_local_zero_requires_an_atom_in_the_declared_read_alphabet(original_packet):
    packet = copy.deepcopy(original_packet)
    row = packet["local_checkpoint_packets"]["r1_s0_a0_f0"]["north"]["north/write"]
    assert "nonexistent_read_atom" not in packet["observers"]["north"]
    row["nonexistent_read_atom"] = 0
    # Omission and zero mass agree only on the declared finite read alphabet;
    # filtering zero entries must not conceal an unbound output label.
    rehash(packet)
    receipt = certify_source_derived_packet(packet)
    assert receipt["status"] == "LOCAL_MARGINAL_MISMATCH"
    json.dumps(receipt, allow_nan=False)


def test_zero_first_presentation_preserves_the_actual_noise_decoder(original_packet):
    packet = copy.deepcopy(original_packet)
    channel = next(c for c in packet["global_checkpoint_kernels"]
                   if c["continuation_id"] == "r1_s0_a0_f0")
    source, row = next(iter(channel["rows"].items()))
    target = next(iter(row))
    zero = next(label for label in channel["rows"] if label != target)
    channel["rows"][source] = {zero: 0, target: 1.0}
    # Local marginals may omit the same zero: both presentations denote the
    # identical measure. The control must decode the actual positive support.
    rehash(packet)
    receipt = certify_source_derived_packet(packet)
    assert receipt["status"] == "PASS"
    noise = receipt["controls"]["full_support_noise"]
    assert noise["tv_identity"] is True
    assert noise["inverse_decoder_worst_input_success"] == pytest.approx(
        1 - noise["mixture_weight"] * Fraction(23, 24), abs=1e-15
    )


def test_actual_observer_domain_is_checked_before_source_certification(original_packet):
    packet = copy.deepcopy(original_packet)
    packet["observers"]["extra_observer"] = packet["observers"]["north"][:]
    # A local twelve-observer declaration cannot silently omit an actual
    # thirteenth observer. Refuse before attempting incompatible section maps.
    assert verify_local_marginal_consistency(packet)["status"] == "LOCAL_MARGINAL_MISMATCH"
    require_refusal(packet)


def test_source_publicness_is_bound_to_actual_authorized_observers(original_packet):
    packet = copy.deepcopy(original_packet)
    packet["publicness_policy"] = [["north"]]
    for channel in packet["global_checkpoint_kernels"]:
        channel["authorized_observers"] = ["north"]
    # This is a valid generic one-observer policy and every source operation
    # still agrees, but it is not the declared universal twelve-port policy.
    require_refusal(packet)


def test_positive_side_branch_is_not_an_explicit_zero(original_packet):
    packet = copy.deepcopy(original_packet)
    channel = packet["global_checkpoint_kernels"][0]
    row = next(iter(channel["rows"].values()))
    target = next(iter(row))
    other = next(label for label in channel["rows"] if label != target)
    row[other] = 1e-100
    # Exact normalization keeps this positive branch even though the dominant
    # weight rounds to one at ordinary float precision.
    assert source_packet._verify_composition(packet)["status"] == "SOURCE_CHECKPOINT_MISMATCH"


@pytest.mark.parametrize("bad_kernel", ["erasure", "randomized"])
def test_noise_control_cannot_certify_a_nonreversible_supplied_kernel(original_packet, bad_kernel):
    packet = copy.deepcopy(original_packet)
    channel = next(c for c in packet["global_checkpoint_kernels"]
                   if c["continuation_id"] == "r1_s0_a0_f0")
    labels = list(channel["rows"])
    row = {labels[0]: 1.0} if bad_kernel == "erasure" else {labels[0]: 0.5, labels[1]: 0.5}
    channel["rows"] = {source: dict(row) for source in labels}
    assert source_packet._full_support_noise_control(packet)["status"] == "FAIL"

"""The declared source generators must produce the certified family."""
import copy
import hashlib
from itertools import product
import json

import pytest

from source_derived_public_checkpoint_packet import (
    build_source_derived_packet,
    certify_source_derived_packet,
)


@pytest.fixture(scope="module")
def source_packet():
    return build_source_derived_packet()


def _certify(packet):
    packet.pop("packet_sha256", None)
    packet["packet_sha256"] = hashlib.sha256(json.dumps(
        packet, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    return certify_source_derived_packet(packet)


def _actual_maps(packet):
    return {channel["continuation_id"]: {
        source: next(output for output, weight in row.items() if weight > 0)
        for source, row in channel["rows"].items()
    } for channel in packet["global_checkpoint_kernels"]}


@pytest.mark.parametrize("case", ["identity_only", "empty", "unknown", "duplicate", "string"])
def test_incomplete_or_malformed_generator_declaration_is_not_certified(source_packet, case):
    packet = copy.deepcopy(source_packet)
    generators = packet["continuation_generators"]
    if case == "identity_only":
        packet["continuation_generators"] = ["r0_s0_a0_f0"]
        identity = _actual_maps(packet)["r0_s0_a0_f0"]
        # Every word in the sole supplied generator fixes every record:
        # its generated group has order one, regardless of the other rows.
        assert all(source == target for source, target in identity.items())
        assert packet["continuation_family_order"] == 40
    elif case == "empty":
        packet["continuation_generators"] = []
    elif case == "unknown":
        generators.append("unbound_geometric_operation")
    elif case == "duplicate":
        generators.append(generators[0])
    else:
        packet["continuation_generators"] = "r1_s0_a0_f0"
    assert _certify(packet)["status"] != "PASS"


@pytest.mark.parametrize("include_identity", [False, True])
def test_inverse_rotation_generates_the_same_source_with_optional_identity(source_packet, include_identity):
    packet = copy.deepcopy(source_packet)
    generators = ["r4_s0_a0_f0", "r0_s1_a0_f0", "r0_s0_a1_f0", "r0_s0_a0_f1"]
    packet["continuation_generators"] = generators + (["r0_s0_a0_f0"] if include_identity else [])
    actual = _actual_maps(packet)
    records = sorted(actual[generators[0]])
    # Independent finite normal-form oracle: R^r S^s A^a F^f, with r<5
    # and s,a,f<2, evaluated on the original supplied transition rows.
    normal_forms = set()
    for powers in product(range(5), range(2), range(2), range(2)):
        images = list(records)
        for name, count in zip(generators, powers):
            for _ in range(count):
                images = [actual[name][record] for record in images]
        normal_forms.add(tuple(images))
    assert len(normal_forms) == 40
    assert normal_forms == {tuple(action[record] for record in records) for action in actual.values()}
    receipt = _certify(packet)
    assert receipt["status"] == "PASS"
    assert receipt["checkpoint_family"]["continuation_count"] == 40
    assert receipt["capacity"]["exact_zero_error_capacity"] == 24

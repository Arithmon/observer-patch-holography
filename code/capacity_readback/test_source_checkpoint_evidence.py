"""Bind the fixed source's diagram and carrier declarations to actual evidence."""
import copy
import hashlib
import json

import pytest

from source_derived_public_checkpoint_packet import (
    build_source_derived_packet,
    certify_source_derived_packet,
)


@pytest.fixture(scope="module")
def packet():
    return build_source_derived_packet()


def certify(packet):
    packet.pop("packet_sha256", None)
    packet["packet_sha256"] = hashlib.sha256(json.dumps(
        packet, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    return certify_source_derived_packet(packet)


@pytest.mark.parametrize("change", [
    "missing_edge", "duplicate_edge", "collapsed_edge_readout",
    "empty_section_inventory", "stale_local_inventory", "stale_interface_inventory",
    "missing_observer_registry", "empty_sparse_projections", "wrong_sparse_entry",
    "wrong_basis", "false_rank_sum", "boolean_basis_index",
    "empty_source_roots", "wrong_carrier_recipe", "unverified_refinement",
])
def test_certificate_requires_the_supplied_source_evidence(packet, change):
    p = copy.deepcopy(packet)
    if change == "missing_edge":
        p["interfaces"].pop()
        p["interface_atom_sets"] = {x["interface_id"]: x["interface_atoms"] for x in p["interfaces"]}
    elif change == "duplicate_edge":
        extra = copy.deepcopy(p["interfaces"][0])
        extra["interface_id"] += "_duplicate"
        p["interfaces"].append(extra)
        p["interface_atom_sets"][extra["interface_id"]] = extra["interface_atoms"]
    elif change == "collapsed_edge_readout":
        interface = p["interfaces"][0]
        # The other 29 edges still force exactly the original public sections.
        # Capacity alone cannot establish this interface's source readout.
        for endpoint in ("left", "right"):
            row = interface[endpoint + "_readout"]
            keys = list(row)
            row[keys[0]] = row[keys[1]]
    elif change == "empty_section_inventory":
        p["public_global_sections"] = []
    elif change == "stale_local_inventory":
        p["local_record_atom_sets"] = {}
    elif change == "stale_interface_inventory":
        p["interface_atom_sets"] = {}
    elif change == "missing_observer_registry":
        p["observer_registry"] = p["observer_registry"][:-1]
    elif change == "empty_source_roots":
        p["carrier_source"]["source_roots"] = []
    elif change == "wrong_carrier_recipe":
        p["carrier_source"]["record_register"] = "one classical bit"
    elif change == "unverified_refinement":
        p["refinement_maps"] = [{"identifies_distinct_records": True}]
    else:
        manifest = p["carrier_projection_manifest"]
        sid, projection = next(iter(manifest["sparse_rank_one_projections"].items()))
        if change == "empty_sparse_projections":
            manifest["sparse_rank_one_projections"] = {}
        elif change == "wrong_sparse_entry":
            projection["entry"][-1] = 0
        elif change == "wrong_basis":
            manifest["basis_order"][0], manifest["basis_order"][1] = manifest["basis_order"][1], manifest["basis_order"][0]
        elif change == "false_rank_sum":
            manifest["rank_sum"] = 0
        else:
            assert projection["basis_index"] == 0
            projection["basis_index"] = False
            projection["entry"][0] = False
            projection["entry"][1] = False
            p["projection_supports"][sid] = [False]
    assert certify(p)["status"] != "PASS", f"unverified source evidence: {change}"


def test_equivalent_diagram_presentations_and_carrier_basis_remain_valid(packet):
    p = copy.deepcopy(packet)
    # Relabel interface atoms, reverse edge presentation, and reorder all
    # inventories without changing the source's equality readouts.
    for interface in p["interfaces"]:
        labels = {old: f"read_{i}" for i, old in enumerate(interface["interface_atoms"])}
        interface["interface_atoms"] = list(reversed(list(labels.values())))
        for side in ("left", "right"):
            field = side + "_readout"
            interface[field] = {k: labels[v] for k, v in interface[field].items()}
        interface["left_observer"], interface["right_observer"] = interface["right_observer"], interface["left_observer"]
        interface["left_readout"], interface["right_readout"] = interface["right_readout"], interface["left_readout"]
    p["interface_atom_sets"] = {x["interface_id"]: x["interface_atoms"] for x in p["interfaces"]}
    p["interfaces"].reverse()
    p["observer_registry"].reverse()
    p["public_global_sections"].reverse()
    manifest = p["carrier_projection_manifest"]
    manifest["basis_order"].reverse()
    for index, slot in enumerate(manifest["basis_order"]):
        sid = p["public_section_aliases"][slot]
        p["projection_supports"][sid] = [index]
        manifest["sparse_rank_one_projections"][sid] = {"basis_index": index, "rank": 1, "entry": [index, index, 1]}
    receipt = certify(p)
    assert receipt["status"] == "PASS"
    assert receipt["capacity"]["exact_zero_error_capacity"] == 24


@pytest.mark.parametrize("change", ["missing_channel_id", "null_channel", "bad_atom_readout"])
def test_malformed_source_is_a_defined_public_refusal(packet, change):
    p = copy.deepcopy(packet)
    if change == "missing_channel_id":
        p["global_checkpoint_kernels"][0].pop("continuation_id")
    elif change == "null_channel":
        p["global_checkpoint_kernels"].append(None)
    else:
        p["observers"]["north"].append("invalid_atom")
    receipt = certify(p)
    assert receipt["status"] != "PASS"
    json.dumps(receipt, allow_nan=False)


@pytest.mark.parametrize("field,value", [
    ("measured_lambda_used", True), ("ew_bridge_target_used", True),
    ("rho_used_as_capacity_producer", True), ("lambda_used", 0),
])
def test_every_declared_target_use_flag_is_checked_without_truth_coercion(packet, field, value):
    p = copy.deepcopy(packet)
    p[field] = value
    assert certify(p)["status"] == "TARGET_TAINTED"


def test_certificate_boundary_is_computed_not_imported_from_free_text(packet):
    p = copy.deepcopy(packet)
    p["claim_boundary"] = "This input claims a hardware and cosmic realization."
    receipt = certify(p)
    assert receipt["status"] == "PASS"
    assert receipt["claim_boundary"] == packet["claim_boundary"]

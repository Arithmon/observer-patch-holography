"""Coherent custody mutations and independent noise/noninterference checks."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import random

import pytest

import verify_transport as verifier

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def packet():
    return verifier.load(HERE/"transport_receipt.json")


def rehash(run):
    previous = "0"*64
    for event in run["events"]:
        event["previous_hash"] = previous
        event["event_hash"] = hashlib.sha256(verifier.raw({k:v for k,v in event.items() if k != "event_hash"})).hexdigest()
        previous = event["event_hash"]
    run["event_root"] = previous


def test_independent_complete_replay(packet):
    assert verifier.verify(packet) == {"episodes":16,"events":4928,"hops":740}


@pytest.mark.parametrize("mutation",["wrong_mean","remote_capture","stale_read","wrong_version","archive_mutable","fake_parent","missing_reset","extra_reset_dependency","wrong_owner","old_version_relabel"])
def test_coherently_rehashed_primitive_mutations_rejected(packet,mutation):
    candidate = deepcopy(packet)
    run = candidate["episodes"][-4]
    if mutation == "wrong_mean":
        event = next(e for e in run["events"] if e["op"] == "mean")
        for write in event["writes"]:
            write["value"] = str(F(write["value"])+1)
    elif mutation == "remote_capture":
        event = next(e for e in run["events"] if e["op"] == "capture")
        previous_mean = run["events"][event["id"]-1]
        event["reads"] = [dict(previous_mean["writes"][0])]
    elif mutation == "stale_read":
        event = next(e for e in run["events"] if e["op"] == "export" and e["arguments"]["message"] == "old_to_B2")
        newer = next(e for e in run["events"] if e["op"] == "commit" and e["label"] == "A1")
        event["reads"] = [dict(newer["writes"][0])]
        event["parents"] = [newer["id"]]
    elif mutation == "wrong_version":
        event = next(e for e in run["events"] if e["op"] == "mean")
        event["reads"][0]["version"] += 1
    elif mutation == "archive_mutable":
        event = next(e for e in run["events"] if e["op"] == "publish")
        event["writes"][0]["immutable"] = False
    elif mutation == "fake_parent":
        event = next(e for e in run["events"] if e["op"] == "capture")
        event["parents"] = []
    elif mutation == "missing_reset":
        index = next(e["id"] for e in run["events"] if e["op"] == "reset" and e["arguments"]["phase"] == "restore")
        del run["events"][index]
    elif mutation == "extra_reset_dependency":
        event = next(e for e in run["events"] if e["op"] == "reset" and e["arguments"]["message"] == "old_to_B2")
        branch = next(e for e in run["events"] if e["op"] == "commit" and e["label"] == "A1")
        event["reads"].append(dict(branch["writes"][0]))
        event["parents"] = sorted(event["parents"]+[branch["id"]])
    elif mutation == "wrong_owner":
        event = next(e for e in run["events"] if e["op"] == "capture")
        event["owner"] = run["path"][0]
    elif mutation == "old_version_relabel":
        event = next(e for e in run["events"] if e["op"] == "export" and e["arguments"]["message"] == "old_to_B2")
        event["label"] = "A1"
    rehash(run)
    with pytest.raises(ValueError):
        verifier.verify(candidate)


@pytest.mark.parametrize("mutation",["capacity","bool_count","missing_variant","fake_order","fake_noise","physical_scope","full_execution","wrong_support"])
def test_claim_and_cost_mutations_rejected(packet,mutation):
    candidate = deepcopy(packet)
    if mutation == "capacity": candidate["episodes"][0]["costs"]["protected_scalar_registers"] -= 1
    elif mutation == "bool_count": candidate["episodes"][0]["costs"]["means"] = True
    elif mutation == "missing_variant": del candidate["episodes"][1]
    elif mutation == "fake_order": candidate["episodes"][0]["logical_events"][-2]["parents"].append("A1")
    elif mutation == "fake_noise": candidate["noise_contract"]["per_hop"] = "mean + read"
    elif mutation == "physical_scope": candidate["scope"] = "PHYSICAL_SOURCE_DERIVED_CHANNEL"
    elif mutation == "full_execution": candidate["specification"]["full_q13_q21_routing_executed"] = True
    elif mutation == "wrong_support": candidate["support_sha256"] = "0"*64
    with pytest.raises(ValueError): verifier.verify(candidate)


def test_sharp_noise_corners_and_unit_gain():
    # Reimplement the primitive arithmetic, without using the producer.
    bounds = [F(1,101),F(1,103),F(1,107),F(1,109),F(1,113)]
    budget = verifier.error_budget(1,F(0),*bounds)
    errors = []
    for signs in itertools.product((-1,1),repeat=5):
        e,z,m,r,a = [v*s for v,s in zip(bounds,signs)]
        x = F(17,19)
        exported = x+e
        receiver_before = z
        mean = (exported+receiver_before)/2+m
        decoded = 2*(mean+r)+a
        errors.append(decoded-x)
        assert abs(decoded-x) <= budget
    assert max(errors) == budget and min(errors) == -budget
    for depth in (0,1,2,4,12,24,100):
        x, observed = F(17,19), F(17,19)+F(1,97)
        for _ in range(depth):
            e,z,m,r,a = bounds
            observed = 2*((observed+e+z)/2+m+r)+a
        assert observed-x == verifier.error_budget(depth,F(1,97),*bounds)


def test_mixed_sign_route_noise_is_bounded():
    rng = random.Random(777)
    bounds = [F(1,17),F(1,19),F(1,23),F(1,29),F(1,31)]
    for depth in range(25):
        initial = F(rng.randrange(-100,100),37)
        x = initial
        for _ in range(depth):
            e,z,m,r,a = [bound*F(rng.randrange(-100,101),100) for bound in bounds]
            x = 2*((x+e+z)/2+m+r)+a
        assert abs(x-initial) <= verifier.error_budget(depth,F(0),*bounds)


def test_branch_and_scratch_interventions_have_required_noninterference(packet):
    values = {(len(r["path"])-1,r["variant"]):{e["name"]:F(e["value"]) for e in r["logical_events"]} for r in packet["episodes"]}
    for depth in (1,4,12,24):
        base = values[depth,"baseline"]
        assert values[depth,"scratch_intervention"] == base
        branch = values[depth,"branch_intervention"]
        assert branch["A1"]-base["A1"] == F(2,9)
        assert branch["C1"]-base["C1"] == F(2,9)
        assert branch["B2"] == base["B2"] and branch["B3"] == base["B3"]


def test_coherent_intervention_on_different_valid_route_rejected(packet):
    candidate = deepcopy(packet)
    run = candidate["episodes"][-3]
    relabel = dict(zip(run["path"],reversed(run["path"])))
    def key(name):
        if name.startswith(("port/","zero/")):
            prefix,vertex = name.split("/")
            return f"{prefix}/{relabel[int(vertex)]}"
        return name
    run["path"] = [relabel[v] for v in run["path"]]
    for request in run["requests"]:
        request["route"] = [relabel[v] for v in request["route"]]
    for event in run["events"]:
        if event["owner"] is not None: event["owner"] = relabel[event["owner"]]
        if "seam" in event["arguments"]: event["arguments"]["seam"] = [relabel[v] for v in event["arguments"]["seam"]]
        for row in event["reads"]+event["writes"]:
            row["owner"] = relabel[row["owner"]]
            row["register"] = key(row["register"])
    for event in run["logical_events"]: event["owner"] = relabel[event["owner"]]
    rehash(run)
    # All primitives and every internal result still pass in isolation.
    support = verifier.load(verifier.ROOT/candidate["specification"]["support"])
    verifier.check_run(run,verifier.verify_support(support),candidate["specification"])
    with pytest.raises(ValueError,match="fixed route"):
        verifier.verify(candidate)


def test_rebased_parent_hash_does_not_hide_changed_census(packet,tmp_path,monkeypatch):
    candidate = deepcopy(packet)
    source = verifier.ROOT/candidate["compiler_targets_not_executed"]["source_path"]
    parent = json.loads(source.read_text())
    next(row for row in parent["levels"] if row["q"] == 21)["site_count"] += 1
    target = tmp_path/candidate["compiler_targets_not_executed"]["source_path"]
    target.parent.mkdir(parents=True)
    target.write_text(json.dumps(parent))
    candidate["compiler_targets_not_executed"]["source_sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    support_source = verifier.ROOT/candidate["specification"]["support"]
    support_target = tmp_path/candidate["specification"]["support"]
    support_target.parent.mkdir(parents=True)
    support_target.write_bytes(support_source.read_bytes())
    monkeypatch.setattr(verifier,"ROOT",tmp_path)
    with pytest.raises(ValueError,match="parent census"):
        verifier.verify(candidate)


@pytest.mark.parametrize("text",['{"x":1,"x":2}','{"x":0.5}','{"x":NaN}'])
def test_nonexact_json_is_rejected(tmp_path,text):
    path = tmp_path/"bad.json"
    path.write_text(text)
    with pytest.raises(ValueError): verifier.load(path)

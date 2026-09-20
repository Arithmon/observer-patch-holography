"""Replay, scientific-boundary and false-green controls for the captured bus."""
import ast
from copy import deepcopy
from fractions import Fraction as F
import inspect
import re
import shlex

import pytest

from . import build, codec, verify


@pytest.fixture(scope="module")
def controls():
    return codec.load(codec.HERE/"controls.json")


def test_independent_replay_and_retained_receipt(controls):
    codec.equal(verify.verify(controls), codec.load(codec.HERE/"receipt.json"), "receipt")


def test_strict_cross_platform_artifact_bytes():
    assert codec.canonical(build.build()) == (codec.HERE/"controls.json").read_bytes()
    assert codec.canonical(verify.verify()) == (codec.HERE/"receipt.json").read_bytes()


def test_verifier_has_no_producer_import():
    tree = ast.parse((codec.HERE/"verify.py").read_text(encoding="utf-8"))
    imports = [node for node in ast.walk(tree) if isinstance(node, (ast.Import,ast.ImportFrom))]
    assert not any("build" in (getattr(node,"module","") or "") or
                   any("build" in name.name for name in node.names) for node in imports)


@pytest.mark.parametrize("mutation", [
    "precision", "initial_error", "per_mean_error", "readout_error", "pin", "port_alias",
    "remote_port", "requests", "round_count", "payload", "prepared_bus", "delete_mean",
    "extra_mean", "changed_edge", "changed_input", "stale_writer", "assigned_reset",
    "wrong_round", "noninteger", "boolean", "decoded", "read_version", "read_use",
    "read_position", "read_writer", "read_units", "final_load", "final_writer",
    "register_cost", "preparation_cost", "mean_cost", "read_cost", "sample_cost",
    "write_cost", "state_bits", "workspace_bits", "program_bits", "port_bits",
    "counter_bits", "version_bits", "request_bits", "receiver_workspace", "read_schedule",
    "receiver_address", "cross_carrier_cost", "missing_history", "extra_field",
])
def test_rejects_corrupt_or_optimistic_artifacts(controls, mutation):
    data = deepcopy(controls)
    case = data["executions"][0]
    if mutation == "precision": data["precision"] = 19
    elif mutation == "initial_error": data["assumed_error_bounds"]["initial"] = "0"
    elif mutation == "per_mean_error": data["assumed_error_bounds"]["extra_per_mean"] = "0"
    elif mutation == "readout_error": data["assumed_error_bounds"]["terminal_per_rail"] = "0"
    elif mutation == "pin": data["source_sha256"][codec.PIN_PATHS[0]] = "0"*64
    elif mutation == "port_alias": data["global_ports"][6] = data["global_ports"][4]
    elif mutation == "remote_port": data["global_ports"][10] = 44
    elif mutation == "requests": data["requests"] = [0,1,1]
    elif mutation == "round_count": data["cleanup_sweeps"] = 79
    elif mutation == "payload": case["payloads"][0] = "-1/2"
    elif mutation == "prepared_bus": case["initial_units"][4] += 1
    elif mutation == "delete_mean": case["tape"].pop(8)
    elif mutation == "extra_mean": case["tape"].append(case["tape"][-1])
    elif mutation == "changed_edge": case["tape"][2][0] = 0
    elif mutation == "changed_input": case["tape"][0][2] += 1
    elif mutation == "stale_writer": case["tape"][9][4] = -5
    elif mutation == "assigned_reset": case["tape"][9][-1] = 2*2**20
    elif mutation == "wrong_round": case["tape"][30][-1] += 1
    elif mutation == "noninteger": case["tape"][0][0] = 0.0
    elif mutation == "boolean": case["tape"][0][0] = False
    elif mutation == "decoded": case["reads"][2]["decoded"] = "3/2"
    elif mutation == "read_version": case["reads"][2]["record"] = 1
    elif mutation == "read_use": case["reads"][2]["use"] = 1
    elif mutation == "read_position": case["reads"][0]["after_mean"] = 6
    elif mutation == "read_writer": case["reads"][1]["writers"][0] = 7
    elif mutation == "read_units": case["reads"][0]["local_units"][0] += 1
    elif mutation == "final_load": case["final_units"][0] += 1
    elif mutation == "final_writer": case["final_writers"][0] = -1
    elif mutation == "missing_history": data["executions"].pop()
    elif mutation == "extra_field": case["remote_payload_side_channel"] = "-3/2"
    else:
        keys = {"register_cost":"scalar_registers", "preparation_cost":"preparation_writes",
                "mean_cost":"means", "read_cost":"mean_reads", "sample_cost":"receiver_scalar_samples",
                "write_cost":"writes", "state_bits":"scalar_storage_bits_bound",
                "workspace_bits":"mean_adder_bits_bound", "program_bits":"expanded_local_schedule_bits",
                "port_bits":"global_port_map_bits", "counter_bits":"program_counter_bits",
                "version_bits":"version_use_counter_bits", "request_bits":"request_record_bits",
                "receiver_workspace":"receiver_arithmetic_bits_bound",
                "read_schedule":"read_schedule_bits", "receiver_address":"receiver_address_bits",
                "cross_carrier_cost":"cross_carrier_means"}
        case["resources"][keys[mutation]] -= 1
    with pytest.raises(ValueError):
        verify.verify(data)


def test_deleted_physical_seam_is_rejected(monkeypatch, controls):
    edges, count = verify.support_edges()
    assert (23,45) in edges
    edges.remove((23,45))
    monkeypatch.setattr(verify,"support_edges",lambda: (edges,count))
    with pytest.raises(ValueError, match="non-native seam"):
        verify.layout(controls["global_ports"])


def test_uniform_coefficient_error_is_not_only_a_sample_check():
    rows, _ = verify.ideal_maps()
    assert rows[2][1] != 0
    corrupt = (rows[2][0],rows[2][1]+F(1,4))
    with pytest.raises(ValueError,match="universal ideal transfer error"):
        verify.read_bounds(2,corrupt)


def test_old_record_after_other_record_with_all_interventions(controls):
    observations = {}
    for case in controls["executions"]:
        p0,p1 = case["payloads"]
        assert [r["decoded"] for r in case["reads"]] == [p0,p1,p0]
        observations[(p0,p1)] = [r["decoded"] for r in case["reads"]]
    levels = ("-3/2","-1/2","1/2","3/2")
    for first in levels:
        assert {observations[(first,other)][2] for other in levels} == {first}
    receipt = codec.load(codec.HERE/"receipt.json")
    for case in receipt["histories"]:
        assert set(range(4)) <= set(case["reads"][2]["initial_ancestor_ports"])


def test_rounded_resource_bounds_and_identity_steps_are_retained(controls):
    for case in controls["executions"]:
        assert all(0 <= row[6] < 2**23 and row[2]+row[3] < 2**24 for row in case["tape"])
        assert any(row[2] == row[3] == row[6] for row in case["tape"])
    altered = deepcopy(controls)
    tape = altered["executions"][0]["tape"]
    index = next(i for i,row in enumerate(tape) if row[2] == row[3] == row[6])
    tape.pop(index)
    with pytest.raises(ValueError,match="missing/extra mean"):
        verify.verify(altered)


@pytest.mark.parametrize("u", range(-10,11))
def test_integer_rounding_bound_and_ties_match_theorem(u):
    for v in range(-10,11):
        result = build.round_mean(u,v)
        assert result == (u+v)//2 + int((u+v)%4 == 3)
        assert abs(2*result-u-v) <= 1
        assert min(u,v) <= result <= max(u,v)
        if (u+v)%2:
            assert result%2 == 0


def test_low_precision_collapses_the_receiver_record():
    case = build.execute((F(1,2),F(-3,2)),precision=2)
    first = case["reads"][0]
    assert first["local_units"][0] == first["local_units"][1]
    assert first["decoded"] is None
    assert verify.decode_receiver(*first["local_units"],1,4) is None
    assert F(8,2*4) > F(1,32)


def test_receiver_uses_local_samples_and_declared_counter_only():
    assert list(inspect.signature(build.decode_local).parameters) == ["positive","negative","use","grid"]
    assert list(inspect.signature(verify.decode_receiver).parameters) == ["positive","negative","use","grid"]
    rows,_ = verify.ideal_maps()
    for cycle,row in enumerate(rows):
        bound = verify.read_bounds(cycle,row)
        radius = F(bound["total_error_bound"])
        use = bound["use"]
        for a in (F(-3,2),F(-1,2),F(1,2),F(3,2)):
            for sign in (-1,1):
                contrast = a/2**(use+3)+sign*radius
                p,m = (2+contrast)*2**20,(2-contrast)*2**20
                assert verify.decode_receiver(p,m,use,2**20) == str(a)


def test_exact_midpoint_is_ambiguous():
    for midpoint in (-1,0,1):
        contrast = F(midpoint,16)
        p,m = (2+contrast)*2**20,(2-contrast)*2**20
        assert build.decode_local(p,m,1,2**20) is None
        assert verify.decode_receiver(p,m,1,2**20) is None


def test_strict_json_rejects_duplicate_keys(tmp_path):
    p=tmp_path/"duplicate.json"
    p.write_text('{"x":1,"x":2}',encoding="utf-8")
    with pytest.raises(ValueError,match="duplicate"):
        codec.load(p)


def test_every_public_theorem_is_transitively_audited():
    proof=(codec.ROOT/"Lean/Geometry/SourceReusableBus.lean").read_text(encoding="utf-8")
    audit=(codec.ROOT/"Lean/Geometry/SourceReusableBusAxiomAudit.lean").read_text(encoding="utf-8")
    names=re.findall(r"^theorem (\w+)",proof,re.M)
    assert len(names)==24
    covered=re.findall(r"^audit_reusable_bus_axioms OPH.SourceReusableBus\.(\w+)$",audit,re.M)
    assert covered==names
    assert "Lean.collectAxioms" in audit and "Lean.ofReduceBool" in audit and "sorryAx" in audit


def test_ci_enforces_audit_on_dependency_changes():
    workflow=(codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    targets=workflow.split("targets=(",1)[1].split(")",1)[0]
    words=shlex.split(targets,comments=True)
    assert "Geometry.SourceReusableBusAxiomAudit" in words
    dedicated=(codec.ROOT/".github/workflows/source-reusable-bus.yml").read_text(encoding="utf-8")
    for path in ("code/source_reusable_bus/**", "Lean/Geometry/SourceReusableBus*.lean",
                 "Lean/Geometry/SourceEncodedMemory.lean", codec.SUPPORT,
                 "Lean/ObserverPatchHolography/ScalarSeamRepair.lean"):
        assert dedicated.count(f'"{path}"')==2
    assert "python -m pytest -q code/source_reusable_bus" in dedicated
    assert "ubuntu-latest" in dedicated and "windows-latest" in dedicated

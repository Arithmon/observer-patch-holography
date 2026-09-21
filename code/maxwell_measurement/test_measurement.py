"""Synthetic conformance and adversarial controls; no physical run is generated."""
from __future__ import annotations

import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

HERE=Path(__file__).resolve().parent


def peer(filename):
    spec=importlib.util.spec_from_file_location("maxwell_measurement_"+filename,HERE/(filename+".py"))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


contract=peer("contract")
verify=peer("verify")


def canonical(value):
    return (json.dumps(value,sort_keys=True,separators=(",",":"))+"\n").encode()


def independently_list_samples(event,edges):
    """Fixture encoding follows public capture roles, without the adapter sampler."""
    out=[]
    def item(role,key,coordinate,is_read=False):
        value=event["reads"][key]["value"] if is_read else event["writes"][key]
        out.append((role,"slot"+str(coordinate),Q(value)))
    op=event["op"]
    if op in ("seed","advance"):
        for coordinate in range(42): item("loaded_"+str(coordinate),"x/"+str(coordinate),coordinate)
    if op=="baseline":
        frame,node=event["args"]
        item("input",f"x/{node}",node,True)
        item("retained_baseline",f"b/{frame}/{node}",node)
    if op=="probe":
        _,edge,side=event["args"]
        parent=edges[edge][side]
        item("before_left",f"x/{parent}",parent,True)
        item("before_right",f"x/{edge+12}",edge+12,True)
        item("after_left",f"x/{parent}",parent)
        item("after_right",f"x/{edge+12}",edge+12)
    if op=="response":
        frame,edge,side=event["args"]
        parent=edges[edge][side]
        item("response",f"x/{parent}",parent,True)
        item("retained_response",f"r/{frame}/{edge}/{side}",parent)
    if op=="feedback":
        _,edge,side=event["args"]
        parent=edges[edge][side]
        item("restored_left",f"x/{parent}",parent)
        item("restored_right",f"x/{edge+12}",edge+12)
    if op=="decode":
        for coordinate in range(42):
            item("decoded_"+str(coordinate),f"d/{event['args'][0]}/{coordinate}",coordinate)
    return out


@pytest.fixture(scope="module")
def fixture():
    p=contract.draft()
    p.update(status="FROZEN",study_id="synthetic-control",run_id="synthetic-only-001",
             capture_kind="synthetic",preregistered_utc="2026-01-02T00:00:00Z")
    p["apparatus"].update(device_id="NO_PHYSICAL_DEVICE",firmware_sha256="1"*64,
                          claimant_organization="synthetic-test-claimant")
    p["clock"].update(seconds_per_tick=["1/1000000","1/1000000"],tick_error="0",
                      reference_sha256="2"*64)
    for k,ch in enumerate(p["channels"]):
        gain=Q(k+100,100)
        if k in (3,12):
            ch["orientation"]=-1
            if k>=12: ch["oriented_ports"].reverse()
            gain=-gain
        offset=Q(k-21,30)
        ch["calibration"]={"id":f"cal{k}","gain":[str(gain)]*2,"offset":[str(offset)]*2,
            "raw_error":"1/1000000000","reference_sha256":hashlib.sha256(f"synthetic-cal{k}".encode()).hexdigest(),
            "reference_samples":[{"id":f"ref{k}-{i}","raw":str((Q(i)-offset)/gain),
                "model_value":str(i),"captured_utc":"2026-01-01T00:00:00Z",
                "source":"external_reference"} for i in (0,1)]}
    source=json.loads((contract.ROOT/contract.PARENT).read_bytes())
    edges=[tuple(ch["oriented_ports"]) for ch in contract.draft()["channels"][12:]]
    c={"schema":"oph.maxwell_measurement.capture.v1","study_id":p["study_id"],
       "run_id":p["run_id"],"capture_kind":"synthetic","device_id":"NO_PHYSICAL_DEVICE",
       "firmware_sha256":"1"*64,"captured_utc":"2026-01-03T00:00:00Z",
       "controls":[],"executions":[]}
    counter=0
    def window(samples):
        nonlocal counter
        start=counter*100;counter+=1
        result=[]
        for i,(role,channel,value) in enumerate(samples):
            cal=p["channels"][int(channel[4:])]["calibration"]
            raw=(value-Q(cal["offset"][0]))/Q(cal["gain"][0])
            result.append({"role":role,"channel":channel,"raw":str(raw),
                           "tick":str(start+i+1),"nonce":f"SYNTHETIC-{counter}-{i}"})
        return {"start_tick":str(start),"end_tick":str(start+99),"samples":result}
    for label in p["controls"]:
        values=[Q(-1,2),Q(3,2)]+([Q(-1,2),Q(3,2)] if label=="sham_no_probe" else [Q(1,2),Q(1,2)])
        roles=["before_left","before_right","after_left","after_right"]
        samples=[(roles[i],"slot0" if i%2==0 else "slot12",value) for i,value in enumerate(values)]
        c["controls"].append({"name":label,**window(samples)})
    for execution in source["executions"]:
        rows=[{"trace":copy.deepcopy(event),**window(independently_list_samples(event,edges))}
              for event in execution["events"]]
        c["executions"].append({"gauge":execution["gauge"],"events":rows})
    return p,c


def save(tmp_path,p,c):
    pp=tmp_path/"prereg.json";cp=tmp_path/"capture.json"
    raw_p,raw_c=canonical(p),canonical(c)
    pp.write_bytes(raw_p);cp.write_bytes(raw_c)
    return pp,cp,hashlib.sha256(raw_p).hexdigest(),hashlib.sha256(raw_c).hexdigest()


def run(tmp_path,fixture):
    return verify.verify(*save(tmp_path,*fixture))


def test_complete_synthetic_interval_readout(tmp_path,fixture):
    result=run(tmp_path,fixture)
    assert result["verdict"]=="SYNTHETIC_CONFORMANCE_ONLY"
    assert result["raw_samples"]==3540
    assert result["measurement_evidence_eligible_relative_to_operator"] is False
    assert result["complete_physical_episode_certified"] is False
    assert result["model_action_time_calibrated"] is False
    parent=json.loads((contract.ROOT/contract.PARENT).read_bytes())
    for data,reference in zip(result["executions"],parent["executions"],strict=True):
        public=reference["events"][-1]["writes"]
        for n,row in enumerate(data["electric"]):
            for edge,interval in enumerate(row):
                assert verify.Interval.parse(interval).contains(Q(public[f"E/{n}/{edge}"]))
        for n,row in enumerate(data["magnetic"][1:],1):
            for face,interval in enumerate(row):
                assert verify.Interval.parse(interval).contains(Q(public[f"B/{n}/{face}"]))


def test_draft_is_not_a_freeze(tmp_path,fixture):
    with pytest.raises(ValueError,match="frozen preregistration"):
        run(tmp_path,(contract.draft(),fixture[1]))


def test_hardware_capture_without_external_trust_stays_unattested(tmp_path,fixture):
    p,c=copy.deepcopy(fixture)
    p["capture_kind"]=c["capture_kind"]="hardware_emulator"
    result=run(tmp_path,(p,c))
    assert result["verdict"]=="CAPTURE_CONFORMS_UNATTESTED"
    assert not result["measurement_evidence_eligible_relative_to_operator"]


def test_synthetic_cannot_be_promoted_by_external_packet(tmp_path,fixture):
    args=save(tmp_path,*fixture)
    trust=tmp_path/"trust.json";trust.write_bytes(b"{}")
    with pytest.raises(ValueError,match="synthetic data cannot"):
        verify.verify(*args,trust,hashlib.sha256(b"{}").hexdigest())


def first_event(c,op):
    return next(e for e in c["executions"][0]["events"] if e["trace"]["op"]==op)


def mutate(label,p,c):
    if label=="missing_channel": p["channels"].pop()
    elif label=="duplicate_channel": p["channels"][1]=p["channels"][0]
    elif label=="wrong_edge_orientation": p["channels"][13]["oriented_ports"].reverse()
    elif label=="gain_sign": p["channels"][12]["calibration"]["gain"]=["1","1"]
    elif label=="wide_uncertainty": p["channels"][0]["calibration"]["raw_error"]="100"
    elif label=="raised_ceiling": p["limits"]["maximum_coordinate_width"]="100"
    elif label=="zero_separation": p["limits"]["minimum_control_separation"]="0"
    elif label=="reversed_interval": p["channels"][0]["calibration"]["gain"]=["2","1"]
    elif label=="heldout_leak": p["channels"][0]["calibration"]["reference_samples"][0]["source"]="slice_2"
    elif label=="post_capture_calibration": p["channels"][0]["calibration"]["reference_samples"][0]["captured_utc"]="2026-01-04T00:00:00Z"
    elif label=="copied_reference": p["channels"][0]["calibration"]["reference_samples"][1]=p["channels"][0]["calibration"]["reference_samples"][0]
    elif label=="post_capture_freeze": p["preregistered_utc"]="2026-01-04T00:00:00Z"
    elif label=="clock_reversed": p["clock"]["orientation"]=-1
    elif label=="vacuous_clock_scale": p["clock"]["seconds_per_tick"]=["1/1000000","1000000"]
    elif label=="clock_is_cycle_count": p["clock"]["raw_unit"]="feedback_count"
    elif label=="model_time_is_device_time": p["clock"]["model_action_time_equals_device_time"]=True
    elif label=="physical_field_promotion": p["apparatus"]["physical_electromagnetic_field_claim"]=True
    elif label=="missing_control": c["controls"].pop()
    elif label=="sham_as_feedback_control": c["controls"][-1]["samples"]=c["controls"][0]["samples"]
    elif label=="missing_event": c["executions"][0]["events"].pop(3)
    elif label=="forged_writer": first_event(c,"probe")["trace"]["reads"]["x/0"]["writer"]=999
    elif label=="synthetic_read": first_event(c,"probe")["trace"]["reads"]["manufactured"]={"writer":0,"value":"0"}
    elif label=="forged_parent": first_event(c,"probe")["trace"]["parents"].append(500)
    elif label=="no_feedback": first_event(c,"feedback")["trace"]["op"]="noop"
    elif label=="wrong_source": c["executions"][0]["events"][1]["trace"]["writes"]["J/0/0"]="0"
    elif label=="bool_event_id": c["executions"][0]["events"][0]["trace"]["id"]=False
    elif label=="bool_gauge": c["executions"][0]["gauge"]=0
    elif label=="integer_causal_false": p["execution"]["record_order_is_physical_causal_order"]=0
    elif label=="missing_probe_sample": first_event(c,"probe")["samples"].pop()
    elif label=="wrong_raw_value": first_event(c,"probe")["samples"][0]["raw"]="999"
    elif label=="wrong_sample_channel": first_event(c,"probe")["samples"][0]["channel"]="slot1"
    elif label=="missing_heldout_coordinate":
        next(e for e in c["executions"][0]["events"] if e["trace"]["op"]=="decode" and e["trace"]["args"]==[2])["samples"].pop()
    elif label=="forged_heldout_decode":
        next(e for e in c["executions"][0]["events"] if e["trace"]["op"]=="decode" and e["trace"]["args"]==[2])["samples"][12]["raw"]="999"
    elif label=="reused_nonce": c["controls"][1]["samples"][0]["nonce"]=c["controls"][0]["samples"][0]["nonce"]
    elif label=="reordered_capture": first_event(c,"probe")["samples"].reverse()
    elif label=="timestamp_outside_event": first_event(c,"probe")["samples"][0]["tick"]="0"
    elif label=="overlapping_events": c["executions"][0]["events"][1]["start_tick"]="0"
    elif label=="unknown_root_key": c["self_attested_physical"]=True
    elif label=="changed_adapter": p["adapter_sources"]["code/maxwell_measurement/verify.py"]="0"*64
    else: raise AssertionError(label)


@pytest.mark.parametrize("label",[
    "missing_channel","duplicate_channel","wrong_edge_orientation","gain_sign",
    "wide_uncertainty","raised_ceiling","zero_separation","reversed_interval","heldout_leak",
    "post_capture_calibration","copied_reference","post_capture_freeze","clock_reversed","vacuous_clock_scale",
    "clock_is_cycle_count","model_time_is_device_time","physical_field_promotion","missing_control",
    "sham_as_feedback_control","missing_event","forged_writer","synthetic_read","forged_parent",
    "no_feedback","wrong_source","bool_event_id","bool_gauge","integer_causal_false","missing_probe_sample",
    "wrong_raw_value","wrong_sample_channel","missing_heldout_coordinate","forged_heldout_decode",
    "reused_nonce","reordered_capture","timestamp_outside_event","overlapping_events","unknown_root_key","changed_adapter"])
def test_coherently_rehashed_mutations_fail(tmp_path,fixture,label):
    p,c=copy.deepcopy(fixture)
    mutate(label,p,c)
    # Recompute both file hashes: pin consistency alone must not manufacture conformance.
    with pytest.raises(ValueError): run(tmp_path,(p,c))


def test_external_file_pin_rejects_mutation(tmp_path,fixture):
    args=save(tmp_path,*fixture)
    args[1].write_bytes(args[1].read_bytes()+b" ")
    with pytest.raises(ValueError,match="file pin mismatch"): verify.verify(*args)


@pytest.mark.parametrize("raw",[b'{"x":1,"x":2}',b'{"x":0.5}',b'{"x":NaN}'])
def test_strict_json_rejects_ambiguous_numbers_and_keys(raw):
    with pytest.raises(ValueError): verify.strict_json(raw)


def test_interval_signed_arithmetic():
    box=verify.Interval(-2,3)
    assert box.square().json()==["0","9"]
    assert (box*verify.Interval(-4,-1)).json()==["-12","8"]
    gain=verify.Interval(Q(-3,2),Q(-1,2))
    assert (gain*verify.Interval(2,3)+verify.Interval(-1,1)).json()==["-11/2","0"]


def attestation_fixture(p,c):
    # An unbound policy-unit-test object, not a measurement receipt or trust root.
    return {"schema":"oph.maxwell_measurement.operator_attestation.v1",
        "study_id":p["study_id"],"run_id":p["run_id"],"preregistration_sha256":"3"*64,
        "capture_sha256":"4"*64,"device_id":c["device_id"],"firmware_sha256":c["firmware_sha256"],
        "preregistered_utc":p["preregistered_utc"],"captured_utc":c["captured_utc"],
        "operator_identity":"test operator","independent_witness_identity":"test witness",
        "scope":{key:True for key in ("physical_device_and_raw_capture_witnessed",
            "actual_read_write_consumption_witnessed","complete_run_and_control_population",
            "firmware_and_device_identity_checked","calibration_references_checked",
            "calibration_independent_of_episode_and_heldout_data","preregistration_precedes_capture",
            "clock_reference_and_capture_timestamps_checked")}}


def test_attestation_identity_does_not_require_three_institutions(fixture):
    p,c=fixture
    witness=attestation_fixture(p,c)
    verify.attest(witness,p,c,"3"*64,"4"*64)
    witness["operator_identity"]=witness["independent_witness_identity"]
    with pytest.raises(ValueError,match="distinct operator"):
        verify.attest(witness,p,c,"3"*64,"4"*64)


@pytest.mark.parametrize("mutation",["hash","scope","bool","run"])
def test_incomplete_or_misbound_attestation_rejected(fixture,mutation):
    p,c=fixture
    witness=attestation_fixture(p,c)
    if mutation=="hash": witness["capture_sha256"]="5"*64
    if mutation=="scope": witness["scope"].pop("actual_read_write_consumption_witnessed")
    if mutation=="bool": witness["scope"]["actual_read_write_consumption_witnessed"]=1
    if mutation=="run": witness["run_id"]="another run"
    with pytest.raises(ValueError): verify.attest(witness,p,c,"3"*64,"4"*64)


def test_fresh_parent_failure_after_success(tmp_path,fixture,monkeypatch):
    run(tmp_path,fixture)
    # A previous successful call must not suppress later parent verification.
    original=verify.fresh_parent
    calls=[]
    def fail():
        calls.append(True)
        raise ValueError("parent recheck failed")
    monkeypatch.setattr(verify,"fresh_parent",fail)
    with pytest.raises(ValueError,match="parent recheck failed"): run(tmp_path,fixture)
    assert calls==[True]
    monkeypatch.setattr(verify,"fresh_parent",original)


def test_changed_parent_source_rejected_after_prior_success(tmp_path,fixture,monkeypatch):
    run(tmp_path,fixture)
    original=Path.read_bytes
    changed=contract.ROOT/"Lean/Screen/SerialMaxwellReadout.lean"
    def changed_bytes(path):
        raw=original(path)
        return raw+b"\n-- synthetic stale-source mutation\n" if path==changed else raw
    # Intercept file reads only; no scientific parent is changed on disk.
    monkeypatch.setattr(Path,"read_bytes",changed_bytes)
    with pytest.raises(ValueError,match="source digest"):
        run(tmp_path,fixture)


def test_cli_drafts_but_refuses_overwrite(tmp_path):
    path=tmp_path/"draft.json"
    command=[sys.executable,str(HERE/"contract.py"),"--output",str(path)]
    first=subprocess.run(command,capture_output=True,text=True)
    assert first.returncode==0 and json.loads(path.read_bytes())["status"]=="DRAFT"
    assert subprocess.run(command,capture_output=True).returncode!=0

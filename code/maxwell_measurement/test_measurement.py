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


def chunked_rational_text(value):
    """Independent base-10 spelling without Decimal or a global digit limit."""
    def integer(number):
        sign="-" if number<0 else ""
        number=abs(number)
        chunks=[]
        while number:
            number,remainder=divmod(number,10**9)
            chunks.append(remainder)
        return sign+(str(chunks[-1])+"".join(f"{part:09d}" for part in reversed(chunks[:-1]))
                     if chunks else "0")
    numerator=integer(value.numerator)
    return numerator if value.denominator==1 else numerator+"/"+integer(value.denominator)


def scalar_uncertainty_oracle(p,c):
    """Original-input center/radius calculation for singleton affine transfers.

    This uses no verifier interval, calibration solver, report, or matrix.
    Linear radii are sums of absolute-coefficient contributions; squared
    radii use the exact closed scalar formula, including zero crossings.
    """
    import re

    edges=[tuple(ch["oriented_ports"]) for ch in contract.draft()["channels"][12:]]
    source=(contract.ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    faces=source.split("def orientedFaces",1)[1].split("def faceEdges",1)[0]
    triangles=[tuple(map(int,t)) for t in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)",faces)]
    assert len(triangles)==20
    transfers=[]
    for ch in p["channels"]:
        cal=ch["calibration"]
        assert cal["gain"][0]==cal["gain"][1] and cal["offset"][0]==cal["offset"][1]
        gain,offset,error=Q(cal["gain"][0]),Q(cal["offset"][0]),Q(cal["raw_error"])
        for ref in cal["reference_samples"]:
            assert gain*Q(ref["raw"])+offset==Q(ref["model_value"])
        transfers.append((gain,offset,abs(gain)*error))

    def measured(sample):
        gain,offset,radius=transfers[int(sample["channel"][4:])]
        return gain*Q(sample["raw"])+offset,radius

    def linear(terms):
        return (sum(coefficient*value[0] for coefficient,value in terms),
                sum(abs(coefficient)*value[1] for coefficient,value in terms))

    def squared(value):
        center,radius=value
        if abs(center)<=radius:
            half=(abs(center)+radius)**2/2
            return half,half
        return center**2+radius**2,2*abs(center)*radius

    def pair(value):
        center,radius=value
        return [chunked_rational_text(center-radius),chunked_rational_text(center+radius)]

    # Every original reading retains one exact gain/offset witness, delta=0.
    checked=0
    for control in c["controls"]:
        expected=[Q(-1,2),Q(3,2)]+([Q(-1,2),Q(3,2)] if control["name"]=="sham_no_probe" else [Q(1,2)]*2)
        assert [measured(sample)[0] for sample in control["samples"]]==expected
        checked+=4
    results=[]
    for execution in c["executions"]:
        baselines,responses,ticks={},{},{}
        for row in execution["events"]:
            event=row["trace"]
            required=independently_list_samples(event,edges)
            assert len(required)==len(row["samples"])
            for sample,(role,channel,value) in zip(row["samples"],required,strict=True):
                assert (sample["role"],sample["channel"])==(role,channel)
                assert measured(sample)[0]==value
                checked+=1
            if event["op"] in ("baseline","response"):
                values={sample["role"]:measured(sample) for sample in row["samples"]}
                target=baselines if event["op"]=="baseline" else responses
                target[tuple(event["args"])]=values["retained_"+event["op"]]
            elif event["op"]=="decode":
                ticks[event["args"][0]]=(Q(row["start_tick"]),Q(row["end_tick"]))
        phi=[[baselines[n,u] for u in range(12)] for n in range(3)]
        potentials=[[linear([(2,responses[n,e,0]),(-1,baselines[n,left])])
                     for e,(left,right) in enumerate(edges)] for n in range(3)]
        electric=[[linear([(-2,potentials[n+1][e]),(2,potentials[n][e]),
                            (-1,phi[n][right]),(1,phi[n][left])])
                   for e,(left,right) in enumerate(edges)] for n in range(2)]
        magnetic=[]
        for n in range(3):
            row=[]
            for triangle in triangles:
                terms=[]
                for left,right in zip(triangle,triangle[1:]+triangle[:1]):
                    oriented=(min(left,right),max(left,right))
                    terms.append((1 if left<right else -1,potentials[n][edges.index(oriented)]))
                row.append(linear(terms))
            magnetic.append(row)
        sources=execution["events"][1]["trace"]["writes"]
        action_terms=[]
        for n in range(2):
            action_terms.extend((Q(1,4),squared(value)) for value in electric[n])
            action_terms.extend((Q(-1,4),squared(value)) for value in magnetic[n+1])
            action_terms.extend((Q(sources[f"J/{n}/{e}"])/2,value) for e,value in enumerate(potentials[n+1]))
            action_terms.extend((Q(sources[f"rho/{n}/{u}"])/2,value) for u,value in enumerate(phi[n]))
        action=pair(linear(action_terms))
        clock=p["clock"]
        assert clock["seconds_per_tick"][0]==clock["seconds_per_tick"][1]
        scale,error=Q(clock["seconds_per_tick"][0]),Q(clock["tick_error"])
        durations=[[str((ticks[n+1][0]-ticks[n][1]-2*error)*scale),
                    str((ticks[n+1][1]-ticks[n][0]+2*error)*scale)] for n in range(2)]
        results.append({"gauge":execution["gauge"],
            "decoded_slices":[list(map(pair,phi[n]+potentials[n])) for n in range(3)],
            "electric":[list(map(pair,row)) for row in electric],
            "magnetic":[list(map(pair,row)) for row in magnetic],
            "field_source_action":action,"device_decode_durations_seconds":durations})
    assert checked==3540
    return results


@pytest.mark.parametrize("entry_point",["api","cli"])
@pytest.mark.parametrize("channel_count,digits",[(6,9),(12,60),(12,197)],
                         ids=["long_action","long_fields","full_input_limit"])
def test_short_inputs_retain_long_exact_outputs(tmp_path,fixture,entry_point,channel_count,digits):
    from sympy import nextprime

    p,c=copy.deepcopy(fixture)
    denominator=10**digits
    for ch in p["channels"][:channel_count]:
        denominator=int(nextprime(denominator))
        error=Q(1,denominator)
        assert 0<error<Q(ch["calibration"]["raw_error"])
        ch["calibration"]["raw_error"]=str(error)
        assert len(str(error))==digits+3
    expected=scalar_uncertainty_oracle(p,c)
    assert all(len(value)>200 for row in expected for value in row["field_source_action"])
    if digits>=60:
        assert any(len(value)>200 for row in expected for field in ("electric","magnetic")
                   for values in row[field] for pair in values for value in pair)
    if digits==197:
        assert any(len(part)>4300 for row in expected for value in row["field_source_action"]
                   for part in value.split("/"))
    result=assert_shared_calibration_result(tmp_path,p,c,True,entry_point)
    assert result["executions"]==expected
    # Ordinary JSON replay preserves every exact output digit.
    assert json.loads(json.dumps(result))["executions"]==expected


def test_interval_json_preserves_decimal_context_and_integer_limit():
    from decimal import Inexact, ROUND_UP, localcontext

    limit=sys.get_int_max_str_digits()
    values=[Q(-(10**5000+1),10**4500+3),Q(10**5000+7,10**4500+3)]
    expected=list(map(chunked_rational_text,values))
    with localcontext() as context:
        context.prec=1
        context.rounding=ROUND_UP
        context.Emin=-1
        context.Emax=1
        for signal in context.traps:
            context.traps[signal]=True
        context.flags[Inexact]=True
        before=context.copy()
        assert verify.Interval(*values).json()==expected
        assert verify.Interval(Q(-2),Q(0)).json()==["-2","0"]
        assert (context.prec,context.rounding,context.Emin,context.Emax)==(
            before.prec,before.rounding,before.Emin,before.Emax)
        assert context.flags==before.flags and context.traps==before.traps
    assert sys.get_int_max_str_digits()==limit


@pytest.mark.parametrize("location",["calibration","capture"])
@pytest.mark.parametrize("entry_point",["api","cli"])
def test_external_rational_length_limit_is_preserved(tmp_path,fixture,location,entry_point):
    p,c=copy.deepcopy(fixture)
    too_long=str(Q(1,10**200+1))
    assert len(too_long)>200
    if location=="calibration":
        p["channels"][0]["calibration"]["raw_error"]=too_long
    else:
        c["controls"][0]["samples"][0]["raw"]=too_long
    args=save(tmp_path,p,c)
    if entry_point=="api":
        with pytest.raises(ValueError,match="bounded rational string required"):
            verify.verify(*args)
    else:
        command=[sys.executable,"-W","error",str(HERE/"verify.py"),
                 "--preregistration",str(args[0]),"--capture",str(args[1]),
                 "--preregistration-sha256",args[2],"--capture-sha256",args[3]]
        process=subprocess.run(command,capture_output=True,text=True,check=False)
        assert process.returncode==2 and not process.stderr
        assert json.loads(process.stdout)=={"verdict":"INVALID","reason":"bounded rational string required"}


def assert_shared_calibration_result(tmp_path,p,c,feasible,entry_point):
    args=save(tmp_path,p,c)
    if entry_point=="api":
        if not feasible:
            with pytest.raises(ValueError,match="calibration"):
                verify.verify(*args)
            return
        result=verify.verify(*args)
    else:
        command=[sys.executable,"-W","error",str(HERE/"verify.py"),
                 "--preregistration",str(args[0]),"--capture",str(args[1]),
                 "--preregistration-sha256",args[2],"--capture-sha256",args[3]]
        process=subprocess.run(command,capture_output=True,text=True,check=False)
        result=json.loads(process.stdout)
        assert not process.stderr
        if not feasible:
            assert process.returncode==2 and result["verdict"]=="INVALID"
            assert "calibration" in result["reason"]
            return
        assert process.returncode==0
    assert result["verdict"]=="SYNTHETIC_CONFORMANCE_ONLY" and result["raw_samples"]==3540
    assert result["physical_truth_certified"] is False
    return result


def independently_admitted_gains(cal,raw,model):
    # Original-input scalar algebra, independent of the verifier's interval code.
    assert cal["offset"]==["-7/10","-7/10"]
    raw,error,target=Q(raw),Q(cal["raw_error"]),Q(model)+Q(7,10)
    assert raw-error>0 and target>0
    lo,hi=map(Q,cal["gain"])
    return max(lo,target/(raw+error)),min(hi,target/(raw-error))


@pytest.mark.parametrize("entry_point",["api","cli"])
@pytest.mark.parametrize("feasible",[False,True])
def test_shared_calibration_references(tmp_path,fixture,entry_point,feasible):
    p,c=copy.deepcopy(fixture)
    cal=p["channels"][0]["calibration"]
    upper=1+Q(1,10**9)
    cal.update(gain=["1",str(upper)],raw_error="0")
    first,second=cal["reference_samples"]
    second_gain=Q(1) if feasible else upper
    second["model_value"]=str(Q(second["raw"])*second_gain-Q(7,10))
    assert independently_admitted_gains(cal,first["raw"],first["model_value"])==(1,1)
    assert independently_admitted_gains(cal,second["raw"],second["model_value"])==(second_gain,second_gain)
    assert (second_gain==1)==feasible
    result=assert_shared_calibration_result(tmp_path,p,c,feasible,entry_point)
    if feasible:
        # References force gain 1, but their constraint cannot narrow the
        # declared uncertainty used to report the retained raw baseline.
        assert result["executions"][0]["decoded_slices"][0][0]==["0",str(Q(7,10)*(upper-1))]


@pytest.mark.parametrize("entry_point",["api","cli"])
@pytest.mark.parametrize("late_gauge",[0,1],ids=["cross_slice","cross_gauge"])
@pytest.mark.parametrize("feasible",[False,True])
def test_shared_calibration_episode(tmp_path,fixture,entry_point,late_gauge,feasible):
    p,c=copy.deepcopy(fixture)
    cal=p["channels"][0]["calibration"]
    epsilon=Q(1,10**9); upper=1+epsilon; error=3*epsilon
    cal.update(gain=["1",str(upper)],raw_error=str(error))
    # References and every unchanged capture row admit either gain endpoint.
    for ref in cal["reference_samples"]:
        for gain in (Q(1),upper):
            assert abs((Q(ref["model_value"])+Q(7,10))/gain-Q(ref["raw"]))<=error
    for row in c["controls"]+[row for execution in c["executions"] for row in execution["events"]]:
        for sample in row["samples"]:
            if sample["channel"]=="slot0":
                assert abs(Q(sample["raw"]))*epsilon/upper<=error
    early=c["executions"][0]["events"][3]
    late=c["executions"][late_gauge]["events"][583]
    assert early["trace"]["op"]=="baseline" and early["trace"]["args"]==[0,0]
    assert late["trace"]["op"]=="decode" and late["trace"]["args"]==[2]
    assert early["samples"][0]["role"]=="input" and early["samples"][0]["channel"]=="slot0"
    assert late["samples"][0]["role"]=="decoded_0" and late["samples"][0]["channel"]=="slot0"
    assert early["trace"]["reads"]["x/0"]["value"]==late["trace"]["writes"]["d/2/0"]=="0"
    early["samples"][0]["raw"]=str(Q(7,10)+error)  # Forces gain exactly 1.
    late_gain=Q(1) if feasible else upper
    late["samples"][0]["raw"]=str(Q(7,10)/late_gain+(error if feasible else -error))
    assert independently_admitted_gains(cal,early["samples"][0]["raw"],"0")== (1,1)
    assert independently_admitted_gains(cal,late["samples"][0]["raw"],"0")== (late_gain,late_gain)
    assert (late_gain==1)==feasible
    assert_shared_calibration_result(tmp_path,p,c,feasible,entry_point)


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


def calibration_corner_oracle(gain, offset, error, samples):
    """Enumerate boundary intersections; evaluate original inverse transfers.

    A nonempty bounded feasible polygon (including a segment or singleton)
    has a vertex at two independent boundary lines. The raw-error boundaries
    use both raw endpoints without branching on the gain sign. Candidate
    acceptance uses (model-offset)/gain, independently of the solver's
    elimination and half-plane sign formulas.
    """
    from itertools import combinations

    lines = [(Q(1), Q(0), gain.lo), (Q(1), Q(0), gain.hi),
             (Q(0), Q(1), offset.lo), (Q(0), Q(1), offset.hi)]
    for raw, model in samples:
        lines.extend((raw + delta, Q(1), model) for delta in (-error, error))
    for (a, b, c), (d, e, f) in combinations(lines, 2):
        determinant = a*e - b*d
        if determinant == 0:
            continue
        g = (c*e - b*f)/determinant
        o = (a*f - c*d)/determinant
        if not (gain.lo <= g <= gain.hi and offset.lo <= o <= offset.hi):
            continue
        if all(abs((model-o)/g - raw) <= error for raw, model in samples):
            return True
    return False


@pytest.mark.parametrize("sign", [-1, 1])
def test_common_calibration_matches_exhaustive_corner_oracle(sign):
    from itertools import combinations_with_replacement, product

    gain = verify.Interval(*sorted((Q(sign), Q(2*sign))))
    offset = verify.Interval(-1, 1)
    alphabet = [(Q(raw), Q(model)) for raw, model in product((-1, 0, 1), (-2, 0, 2))]
    outcomes = set()
    for error in (Q(0), Q(1, 2)):
        for samples in combinations_with_replacement(alphabet, 3):
            expected = calibration_corner_oracle(gain, offset, error, samples)
            actual = verify.common_calibration(gain, offset, error, samples)
            assert actual is expected, (sign, error, samples)
            outcomes.add(expected)
    assert outcomes == {False, True}


def test_common_calibration_matches_seeded_corner_oracle():
    import random

    rng = random.Random(1033)
    outcomes = set()
    for case in range(160):
        sign = rng.choice((-1, 1))
        endpoints = sorted(Q(sign*rng.randrange(1, 13), 4) for _ in range(2))
        gain = verify.Interval(endpoints[0], endpoints[0] if case % 5 == 0 else endpoints[1])
        endpoints = sorted(Q(rng.randrange(-8, 9), 4) for _ in range(2))
        offset = verify.Interval(endpoints[0], endpoints[0] if case % 7 == 0 else endpoints[1])
        error = Q(rng.randrange(4), 8)
        g = gain.lo + Q(rng.randrange(5), 4)*(gain.hi-gain.lo)
        o = offset.lo + Q(rng.randrange(5), 4)*(offset.hi-offset.lo)
        samples = []
        for _ in range(rng.randrange(7)):
            raw = Q(rng.randrange(-12, 13), 4)
            model = (g*(raw + error*Q(rng.randrange(-2, 3), 2)) + o
                     if case % 2 == 0 else Q(rng.randrange(-20, 21), 4))
            samples.append((raw, model))
        expected = calibration_corner_oracle(gain, offset, error, samples)
        if case % 2 == 0:
            assert expected, "constructed witness must satisfy the original transfer"
        assert verify.common_calibration(gain, offset, error, samples) is expected, case
        outcomes.add(expected)
    assert outcomes == {False, True}


@pytest.mark.parametrize("sign", [-1, 1])
def test_common_calibration_pairwise_feasible_triple_is_jointly_impossible(sign):
    from itertools import combinations

    delta = Q(1, 10000)
    gain = verify.Interval(*sorted((sign*(1-delta), sign*(1+delta))))
    offset = verify.Interval(-delta, delta)
    samples = [(Q(0), Q(0)), (Q(sign), Q(1)), (Q(2*sign), 2+delta)]
    for pair in combinations(samples, 2):
        assert calibration_corner_oracle(gain, offset, Q(0), pair)
        assert verify.common_calibration(gain, offset, Q(0), pair)
    assert not calibration_corner_oracle(gain, offset, Q(0), samples)
    assert not verify.common_calibration(gain, offset, Q(0), samples)


@pytest.mark.parametrize("sign", [-1, 1])
def test_common_calibration_exact_raw_error_boundary(sign):
    gain = verify.Interval(sign)
    offset = verify.Interval(-1, 1)
    samples = [(Q(0), Q(-5, 4)), (Q(sign), Q(1, 4))]
    for error, expected in ((Q(1, 4), True), (Q(1, 4)-Q(1, 10**30), False)):
        assert calibration_corner_oracle(gain, offset, error, samples) is expected
        assert verify.common_calibration(gain, offset, error, samples) is expected


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("shape",["line","corner","point"])
def test_common_calibration_retains_degenerate_feasible_sets(sign,shape):
    gain=verify.Interval(*sorted((Q(sign),Q(sign if shape=="point" else 2*sign))))
    offset=verify.Interval(-1,1) if shape=="corner" else verify.Interval(0)
    samples=([(Q(0),Q(-1)),(Q(sign),Q(0))] if shape=="corner"
             else [(Q(0),Q(0))])
    assert calibration_corner_oracle(gain,offset,Q(0),samples)
    assert verify.common_calibration(gain,offset,Q(0),samples)


@pytest.mark.parametrize("sign", [-1, 1])
def test_preregistration_rejects_pairwise_feasible_reference_triple(fixture,sign):
    from itertools import combinations

    p,_=copy.deepcopy(fixture)
    cal=p["channels"][0 if sign==1 else 3]["calibration"]
    delta=Q(1,10000)
    cal.update(gain=list(map(str,sorted((sign*(1-delta),sign*(1+delta))))),
               offset=[str(-delta),str(delta)],raw_error="0")
    references=[{**cal["reference_samples"][0],"id":f"triple-{i}",
                 "raw":str(sign*raw),"model_value":str(model)}
                for i,(raw,model) in enumerate(((0,0),(1,1),(2,2+delta)))]
    edges=[tuple(ch["oriented_ports"]) for ch in contract.draft()["channels"][12:]]
    for pair in combinations(references,2):
        cal["reference_samples"]=list(pair)
        verify.validate_prereg(p,edges)
    cal["reference_samples"]=references
    with pytest.raises(ValueError,match="calibration references have no common"):
        verify.validate_prereg(p,edges)


@pytest.mark.parametrize("feasible",[False,True])
def test_shared_calibration_includes_references_with_capture(tmp_path,fixture,feasible):
    p,c=copy.deepcopy(fixture)
    cal=p["channels"][0]["calibration"]
    upper=1+Q(1,10**9)
    cal.update(gain=["1",str(upper)],raw_error="0")
    capture_gain=Q(1) if feasible else upper
    for row in c["controls"]+[row for execution in c["executions"] for row in execution["events"]]:
        for sample in row["samples"]:
            if sample["channel"]=="slot0":
                sample["raw"]=str(Q(sample["raw"])/capture_gain)
    ref=cal["reference_samples"][0]
    assert independently_admitted_gains(cal,ref["raw"],ref["model_value"])==(1,1)
    raw=c["controls"][0]["samples"][0]["raw"]
    assert independently_admitted_gains(cal,raw,Q(-1,2))==(capture_gain,capture_gain)
    assert_shared_calibration_result(tmp_path,p,c,feasible,"api")


@pytest.mark.parametrize("feasible",[False,True])
def test_shared_calibration_includes_controls_with_episode(tmp_path,fixture,feasible):
    p,c=copy.deepcopy(fixture)
    cal=p["channels"][0]["calibration"]
    epsilon=Q(1,10**9); upper=1+epsilon; error=3*epsilon
    cal.update(gain=["1",str(upper)],raw_error=str(error))
    control=c["controls"][0]["samples"][0]
    control["raw"]=str(Q(1,5)+error)
    late=c["executions"][1]["events"][583]["samples"][0]
    late_gain=Q(1) if feasible else upper
    late["raw"]=str(Q(7,10)/late_gain+(error if feasible else -error))
    assert independently_admitted_gains(cal,control["raw"],Q(-1,2))==(1,1)
    assert independently_admitted_gains(cal,late["raw"],0)==(late_gain,late_gain)
    assert_shared_calibration_result(tmp_path,p,c,feasible,"api")


def test_shared_calibration_is_fresh_for_each_verification(tmp_path,fixture):
    upper=1+Q(1,10**9)
    for admitted_gain in (Q(1),upper,Q(1)):
        p,c=copy.deepcopy(fixture)
        cal=p["channels"][0]["calibration"]
        cal.update(gain=["1",str(upper)],raw_error="0")
        for ref in cal["reference_samples"]:
            ref["raw"]=str(Q(ref["raw"])/admitted_gain)
        for row in c["controls"]+[row for execution in c["executions"] for row in execution["events"]]:
            for sample in row["samples"]:
                if sample["channel"]=="slot0":
                    sample["raw"]=str(Q(sample["raw"])/admitted_gain)
        assert_shared_calibration_result(tmp_path,p,c,True,"api")


@pytest.mark.parametrize("gain,error",[("zero",Q(0)),("crossing",Q(0)),("positive",Q(-1))])
def test_common_calibration_rejects_invalid_solver_domain(gain,error):
    box={"zero":verify.Interval(0,1),"crossing":verify.Interval(-1,1),
         "positive":verify.Interval(1,2)}[gain]
    with pytest.raises(ValueError,match="calibration"):
        verify.common_calibration(box,verify.Interval(0),error,[])

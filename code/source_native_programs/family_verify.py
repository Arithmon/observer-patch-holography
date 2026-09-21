"""Independent full metric census, route checking and q=3 native program replay."""
import argparse
import json
from math import isqrt
from pathlib import Path

import numpy as np

from . import codec
from .check_routes import check,read_rows,edges,require
from .verify import certify_plan,replay


def metric_menu(q):
    require(type(q) is int and q in (3,13,21),"declared family")
    coords = np.array(list(np.ndindex(q,q,q)),dtype=np.int64)
    floors = np.array([(i+isqrt(5*i*i))//2 for i in range(q)],dtype=np.int64)
    integer = -floors[coords]
    rows = []
    # All pairs, including the exact boundary and self, tested independently
    # in Q(sqrt(5)). For q<=21 the squared expressions fit signed int64.
    for start in range(0,q**3,32):
        da = integer[start:start+32,None,:]-integer[None,:,:]
        db = coords[start:start+32,None,:]-coords[None,:,:]
        cross = (2*da*db+db*db).sum(axis=2)
        c = q*(2*(da*da+db*db).sum(axis=2)+cross)-2
        r = q*cross
        signs = np.where(c*r >= 0,np.sign(c+r),np.sign(c)*np.sign(c*c-5*r*r))
        rows.extend(np.flatnonzero(mask).tolist() for mask in signs <= 0)
    if q == 3:
        prior = codec.load(codec.ROOT/"code/source_read_routing/controls/q3_baseline.json")
        require(codec.digest(rows) == prior["inputs"]["neighbour_sha256"],"q3 metric custody")
    else:
        # This pinned historical receipt also contains unrelated floating-point
        # diagnostics. Its metric commitment is a string; no diagnostic enters
        # the exact integer reconstruction above. Reject duplicate keys and
        # nonfinite constants while admitting that historical numeric schema.
        from source_read_acceptance.codec import unique,reject
        prior = json.loads((codec.ROOT/"evidence/source_net_causal_poset/carrier_source_net_receipt.json").read_text(
            encoding="ascii"),object_pairs_hook=unique,parse_constant=reject)
        require(prior["schema"] == "oph.exact.carrier-source-net.v1","historical metric receipt")
        entry = next(x for x in prior["levels"] if x["q"] == q)
        require(codec.digest(rows) == entry["provenance"]["read_relation_sha256"],"complete metric custody")
    return rows


def check_bound(bound, q, menu, route_summary):
    size = q*q*q
    rounds = isqrt(q-1)+1
    fan_in = max(len(row) for row in menu)
    use_count = [0]*size
    for row in menu:
        for source in row:
            use_count[source] += 1
    fan_out = max(use_count)
    length = route_summary["max_paired_cells"]
    commits, reads = size*rounds,sum(map(len,menu))*rounds
    # Induct over layers. Every old-bank record is read at most fan_out
    # times before retirement; there is no circular within-layer dependence.
    exponent = commits
    for _ in range(rounds):
        exponent += fan_out+fan_in+2*length+1
    stages = 2*commits+reads
    k = exponent+(64*(size+1)*stages*(length-1)**2-1).bit_length()
    work = (commits*(3*k+3)+(commits+reads)*(2*length+1+(2*length-3)*(length-1)**2*k)
            +reads*(3*exponent+12))
    expected = {"q":q,"sites":size,"rounds":rounds,"complete_metric_reads":reads,
                "metric_menu_sha256":codec.digest(menu),"max_reads_per_site":fan_in,
                "max_reads_of_record":fan_out,"input_bound":size+1,
                "maximum_route_cells":length,"maximum_scale_upper":exponent,
                "cleanup_blocks_sufficient":k,"scalar_means_upper":work,
                "grid_bits_sufficient":exponent+(64*(work+1)-1).bit_length(),
                "physical_scalar_registers":12*20*4**({3:3,13:4,21:5}[q]),
                "scope":"analytic_finite_compilation_bound_not_executed_production_word"}
    codec.equal(bound,expected,"family resource bound or scope")


def summarize(case):
    committed = []
    samples = 0
    for checkpoint in case["checkpoints"]:
        samples += 2*len(checkpoint["observations"])
        if checkpoint["kind"] == "store":
            version = [checkpoint["layer"],checkpoint["site"],-1]
            values = [o["decoded"] for o in checkpoint["observations"] if o["version"] == version]
            require(len(values) == 1,"one publication per committed version")
            committed.append(values[0])
    return {"payload":case["payload"],"history_sha256":codec.digest(case),
            "scalar_means":case["scalar_means"],"evaluated_means":case["evaluated_means"],
            "stationary_means":case["stationary_means"],"stored_values":committed,
            "scalar_samples":samples,"plan_sha256":case["plan_sha256"],
            "final_state_sha256":case["final_state_sha256"]}


def verify(packet, directory):
    require(type(packet) is dict and set(packet) == {"schema","pins","source_selected","m1_derived",
            "routes","bounds","q3_plan","q3_cases","production_native_replay","q3_native_replay","placement"},
            "family artifact fields")
    codec.equal(packet["schema"],"oph-native-program-families-v1","family schema")
    for flag in ("source_selected","m1_derived","production_native_replay"):
        codec.equal(packet[flag],False,"family scope promotion")
    codec.equal(packet["q3_native_replay"],True,"missing q3 native execution")
    codec.equal(packet["placement"],"carrier_zero_core_banks_one_through_twice_sites","supplied placement")
    pins = codec.pins()
    codec.equal(packet["pins"],pins,"family custody")
    require(type(packet["routes"]) is list and len(packet["routes"]) == 3 and
            type(packet["bounds"]) is list and len(packet["bounds"]) == 3,"complete family inventory")
    for index,(q,level) in enumerate(((3,3),(13,4),(21,5))):
        summary = check(read_rows(directory/f"routes-{q**3}.jsonl"),q**3,level)
        codec.equal(packet["routes"][index],summary,"route certificate")
        menu = metric_menu(q)
        check_bound(packet["bounds"][index],q,menu,summary)
        print(f"q={q}: all metric decisions and paired routes verified",flush=True)
        if q == 3:
            _,route_map = check(read_rows(directory/"routes-27.jsonl"),27,3,True)
            plan = codec.load_artifact(directory/"q3-plan.json")
            certify_plan(plan,[menu,menu],route_map)
            require(plan["input_bound"] == 28,"q3 initial bound")
            codec.equal(packet["q3_plan"],{k:v for k,v in plan.items() if k not in ("segments","layers")},
                        "q3 exact plan metadata")
            require(type(packet["q3_cases"]) is list and len(packet["q3_cases"]) == 2,"q3 controls")
            for intervention in (0,1):
                payload = list(range(1,28))
                payload[13] += intervention
                result = summarize(replay(plan,payload,edges(3)))
                codec.equal(packet["q3_cases"][intervention],result,"q3 native replay or intervention")
                print(f"q=3 intervention={intervention}: complete native program verified",flush=True)
    codec.equal(codec.pins(),pins,"family sources changed during verification")
    return {"schema":"oph-native-program-families-receipt-v1","pins":pins,
            "families_sha256":codec.digest(packet),"route_witnesses":sum(r["routes"] for r in packet["routes"]),
            "q3_histories":2,"q3_scalar_means":sum(c["scalar_means"] for c in packet["q3_cases"]),
            "q3_evaluated_means":sum(c["evaluated_means"] for c in packet["q3_cases"]),
            "q3_stationary_means":sum(c["stationary_means"] for c in packet["q3_cases"]),
            "production_native_replay":False,"m1_derived":False,"source_selected":False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir",type=Path,required=True)
    parser.add_argument("--families",type=Path,default=codec.HERE/"families.json")
    parser.add_argument("--receipt",type=Path,default=codec.HERE/"families_receipt.json")
    parser.add_argument("--write-receipt",action="store_true")
    args = parser.parse_args()
    receipt = verify(codec.load_artifact(args.families),args.workdir)
    if args.write_receipt:
        args.receipt.write_bytes(codec.canonical(receipt))
    else:
        codec.equal(codec.load_artifact(args.receipt),receipt,"family receipt mismatch")


if __name__ == "__main__":
    main()

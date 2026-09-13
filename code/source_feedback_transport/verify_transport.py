"""Independent primitive, locality, version and causal-projection replay.

No transport producer is imported. The earlier independent checker is reused
only for strict parsing and the complete captured W12 support combinatorics.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "source_routing"))
from verify_routing import load, verify_support


def require(condition, message):
    if not condition:
        raise ValueError(message)


def raw(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("utf-8")


def equal(actual, expected, message):
    require(raw(actual) == raw(expected), message)


def fraction(value):
    require(type(value) is str, "rational must have exact string encoding")
    result = F(value)
    require(str(result) == value, "noncanonical rational")
    return result


def parent_census(path):
    # The older parent also contains decimal diagnostics. Parse those exactly
    # without importing them into this certificate's integer census.
    def pairs(items):
        result = {}
        for key,value in items:
            require(key not in result,"duplicate parent JSON key")
            result[key] = value
        return result
    def reject(value):
        raise ValueError("nonfinite parent diagnostic: "+value)
    return json.loads(path.read_text(encoding="utf-8"),object_pairs_hook=pairs,parse_float=F,parse_constant=reject)


def error_budget(depth, initial, export, reset, mean, read, decode):
    require(type(depth) is int and depth >= 0, "nonnegative integral depth")
    require(all(x >= 0 for x in (initial, export, reset, mean, read, decode)), "nonnegative noise bounds")
    return initial + depth*(export+reset+2*mean+2*read+decode)


def check_run(run, edges, specification):
    path, variant = run["path"], run["variant"]
    require(all(type(v) is int and 0 <= v < 15360 for v in path), "path domain")
    require(len(path)-1 in specification["path_lengths"] and len(path) == len(set(path)), "simple declared-depth path")
    require(all(tuple(sorted((a,b))) in edges for a,b in zip(path,path[1:])), "non-support hop")
    require(variant in specification["variants"], "undeclared variant")
    a, b, c = path[0], path[-1], path[max(1,(len(path)-1)//2)]
    owners = {"A0":a,"B0":b,"C0":c,"B1":b,"A1":a,"C1":c,"B2":b,"B3":b}
    dependencies = {"A0":[],"B0":[],"C0":[],"B1":["A0","B0"],"A1":["A0","B1"],
                    "C1":["A1","C0"],"B2":["A0","B1"],"B3":["A0","B2"]}
    logical_ancestors = {}
    for name, parents in dependencies.items():
        logical_ancestors[name] = {name}.union(*(logical_ancestors[p] for p in parents))
    # Expected event interfaces are fixed independently of payload values.
    interfaces = []
    def add(op, owner, label, reads, writes, arguments):
        interfaces.append((op,owner,label,reads,writes,arguments))
    for vertex in path:
        add("initialize_zero",vertex,None,[],[f"zero/{vertex}"],{})
        add("initialize_scratch",vertex,None,[],[f"port/{vertex}"],{})
    for name in ("A0","B0","C0"):
        add("publish",owners[name],name,[],[f"archive/{name}"],{"name":name})
    recipe = [("to_B1","A0",path,"B1","B0"),
              ("to_A1","B1",path[::-1],"A1","A0"),
              ("to_C1","A1",path[:path.index(c)+1],"C1","C0"),
              ("old_to_B2","A0",path,"B2","B1"),
              ("old_to_B3","A0",path,"B3","B2")]
    expected_requests = []
    for message,origin,route,target,local in recipe:
        start = len(interfaces)
        archive = f"archive/{origin}"
        for index,(u,v) in enumerate(zip(route,route[1:])):
            args = {"message":message,"index":index}
            add("export",u,origin,[archive],[f"port/{u}"],args)
            add("reset",v,None,[f"zero/{v}"],[f"port/{v}"],dict(args,phase="prepare"))
            add("mean",None,origin,[f"port/{u}",f"port/{v}"],[f"port/{u}",f"port/{v}"],dict(args,seam=[u,v]))
            archive = f"buffer/{message}/{index}"
            add("capture",v,origin,[f"port/{v}"],[archive],args)
            for vertex in (u,v):
                add("reset",vertex,None,[f"zero/{vertex}"],[f"port/{vertex}"],dict(args,phase="restore"))
        expected_requests.append({"message":message,"route":route,"origin":origin,
                                  "first_event":start,"event_count":len(interfaces)-start,"delivered_register":archive})
        delta = specification["branch_intervention"] if target == "A1" and variant == "branch_intervention" else "0"
        add("commit",owners[target],target,[archive,f"archive/{local}"],[f"archive/{target}"],{"name":target,"intervention":delta})
    equal(run["requests"],expected_requests,"request schedule/version mismatch")
    events = run["events"]
    require(len(events) == len(interfaces), "missing or extra primitive event")
    state, ancestors, labels, logical_rows = {}, [], [], []
    previous = "0"*64
    for i,(event, interface) in enumerate(zip(events,interfaces)):
        op,owner,label,reads,writes,args = interface
        equal({k:event[k] for k in ("id","op","owner","label","arguments")},
              {"id":i,"op":op,"owner":owner,"label":label,"arguments":args},"primitive interface mismatch")
        require(event["previous_hash"] == previous,"broken custody chain")
        require(event["event_hash"] == hashlib.sha256(raw({k:v for k,v in event.items() if k != "event_hash"})).hexdigest(),"event hash")
        previous = event["event_hash"]
        equal(event["reads"],[{"register":key,**state[key]} for key in reads],"stale, remote or missing read")
        if op != "mean":
            require(all(state[key]["owner"] == owner for key in reads), "remote local-operation read")
        parents = sorted({state[key]["writer"] for key in reads})
        equal(event["parents"],parents,"forged semantic parent")
        ancestry = 0
        for p in parents:
            require(p < i,"noncausal writer")
            ancestry |= ancestors[p] | (1 << p)
            parent_label = labels[p]
            if parent_label is not None:
                require(label is not None and parent_label in logical_ancestors[label],"spurious logical dependency")
        ancestors.append(ancestry)
        labels.append(label)
        values = [fraction(state[key]["value"]) for key in reads]
        if op == "initialize_zero" or op == "reset":
            result = [F(0)]
            if op == "reset":
                require(values == [0] and state[reads[0]]["immutable"],"reset must read local protected zero only")
        elif op == "initialize_scratch":
            position = path.index(owner)
            result = [F(100+3*position,17) + (F(position+1,19) if variant == "scratch_intervention" else 0)]
        elif op == "publish":
            result = [fraction(specification["archive_values"][label]) +
                      (fraction(specification["source_intervention"]) if label == "A0" and variant == "source_intervention" else 0)]
        elif op == "export":
            require(state[reads[0]]["immutable"],"export must consume retained version")
            result = values
        elif op == "mean":
            require(values[0] != values[1],"the retained episode unexpectedly contains an inactive mean")
            result = [sum(values)/2]*2
        elif op == "capture":
            result = [2*values[0]]
        elif op == "commit":
            require(all(state[key]["immutable"] for key in reads),"logical commit must use immutable versions")
            result = [1+sum(values)+fraction(args["intervention"])]
        else:
            raise ValueError("unknown primitive")
        output = []
        for key,value in zip(writes,result):
            require(key not in state or not state[key]["immutable"],"immutable overwrite")
            destination = int(key.split("/")[1]) if op == "mean" else owner
            immutable = op in ("initialize_zero","publish","capture","commit")
            output.append({"register":key,"owner":destination,"value":str(value),
                           "version":state[key]["version"]+1 if key in state else 0,
                           "writer":i,"immutable":immutable})
        equal(event["writes"],output,"incorrect write law, local owner, version or archive flag")
        state.update({row["register"]:{k:v for k,v in row.items() if k != "register"} for row in output})
        if op in ("publish","commit"):
            logical_rows.append({"name":label,"owner":owner,"parents":dependencies[label],"event_id":i,"value":str(result[0])})
    equal(run["logical_events"],logical_rows,"logical event projection/value mismatch")
    require(run["event_root"] == previous,"event root")
    require(all(fraction(state[f"port/{v}"]["value"]) == 0 for v in path),"scratch not restored")
    # Full induced order, not just direct edges or delivered scalar equality.
    for left in logical_rows:
        for right in logical_rows:
            operational = left["event_id"] == right["event_id"] or bool(ancestors[right["event_id"]] & (1 << left["event_id"]))
            expected = left["name"] in logical_ancestors[right["name"]]
            require(operational == expected,"logical reachability differs from induced operational order")
    for request in expected_requests:
        require(fraction(state[request["delivered_register"]]["value"]) == fraction(state[f"archive/{request['origin']}"]["value"]),"delivered wrong immutable version")
    counts = Counter(e["op"] for e in events)
    hops = sum(len(r["route"])-1 for r in expected_requests)
    expected_costs = {"hops":hops,"means":counts["mean"],"exports":counts["export"],"captures":counts["capture"],"resets":counts["reset"],
                      "transport_events":6*hops,"initialization_events":2*len(path)+3,"logical_commits":5,"total_events":len(events),
                      "register_reads":sum(len(e["reads"]) for e in events),"register_writes":sum(len(e["writes"]) for e in events),
                      "protected_scalar_registers":sum(v["immutable"] for v in state.values()),"scratch_scalar_registers":len(path),
                      "largest_written_scalar_bits":max(abs(fraction(v["value"]).numerator).bit_length()+fraction(v["value"]).denominator.bit_length()+1 for e in events for v in e["writes"])}
    equal(run["costs"],expected_costs,"cost or capacity undercount")
    return {row["name"]:fraction(row["value"]) for row in logical_rows}


def verify(packet=None):
    specification = load(HERE/"specification.json")
    packet = load(HERE/"transport_receipt.json") if packet is None else packet
    equal(packet["specification"],specification,"frozen specification")
    require(packet["schema"] == "oph.source_feedback_transport.receipt.v1","schema")
    require(packet["scope"] == "CLASSICAL_REUSABLE_LOCAL_FEEDBACK_TRANSPORT__SUPPLIED_COPY_RESET_LAW__NO_PHYSICAL_OR_QUANTUM_IDENTIFICATION","scope promotion")
    support_path = ROOT/specification["support"]
    require(packet["support_sha256"] == hashlib.sha256(support_path.read_bytes()).hexdigest(),"support pin")
    edges = verify_support(load(support_path))
    expected_keys = {(d,v) for d in specification["path_lengths"] for v in specification["variants"]}
    baselines = {len(run["path"])-1:run["path"] for run in packet["episodes"] if run["variant"] == "baseline"}
    results = {}
    for run in packet["episodes"]:
        key = (len(run["path"])-1,run["variant"])
        require(key not in results,"duplicate episode")
        equal(run["path"],baselines.get(key[0]),"intervention changed the fixed route")
        results[key] = check_run(run,edges,specification)
    require(set(results) == expected_keys,"missing declared depth/intervention")
    # Interventions are checked against a separate closed logical recurrence.
    for d in specification["path_lengths"]:
        for variant in specification["variants"]:
            s,t,u = (fraction(specification["archive_values"][k]) for k in ("A0","B0","C0"))
            if variant == "source_intervention": s += fraction(specification["source_intervention"])
            branch = fraction(specification["branch_intervention"]) if variant == "branch_intervention" else F(0)
            b1 = 1+s+t
            a1 = 1+b1+s+branch
            b2 = 1+s+b1
            expected = {"A0":s,"B0":t,"C0":u,"B1":b1,"A1":a1,"C1":1+a1+u,"B2":b2,"B3":1+s+b2}
            require(results[d,variant] == expected,"logical intervention compatibility")
        for k in ("A0","B0","C0","B1","B2","B3"):
            require(results[d,"branch_intervention"][k] == results[d,"baseline"][k],"branch contaminated unrelated version read")
        require(results[d,"scratch_intervention"] == results[d,"baseline"],"scratch preparation changes logical world")
    equal(packet["noise_contract"],{"per_hop":"export + reset + 2*mean + 2*read + decode","depth_d":"initial_archive_error + d*(export + reset + 2*mean + 2*read + decode)",
                                  "scope":"Separately bounded additive errors; export includes archive retention/read/write error relative to its retained payload. Mean error bounds the receiver coordinate. Exact deterministic schedule; no measured physical noise limits."},"noise scope/formula")
    source_path = ROOT/"evidence/source_net_causal_poset/carrier_source_net_receipt.json"
    target = packet["compiler_targets_not_executed"]
    require(target["source_path"] == source_path.relative_to(ROOT).as_posix() and target["source_sha256"] == hashlib.sha256(source_path.read_bytes()).hexdigest(),"compiler target pin")
    expected_rows = []
    parent_rows = {row["q"]:row for row in parent_census(source_path)["levels"]}
    for q,n,r,e in ((13,2197,4,145997),(21,9261,5,1451292)):
        parent = parent_rows[q]
        equal([parent["site_count"],parent["rounds"],parent["neighbours"]["undirected_spatial_edges"],parent["operation_costs"]["total_reads"]],
              [n,r,e,(2*e+n)*r],"pinned parent census disagrees with compiler target")
        expected_rows.append({"q":q,"sites":n,"rounds":r,"undirected_metric_edges":e,"logical_reads":(2*e+n)*r,"nonself_reads":2*e*r,"six_event_unicast_transport_lower_bound":12*e*r,
                              "scope":"Declared compiler target count only; assumes distinct source/destination scratch vertices; excludes routing distance, local reads and commits. Not a full execution or a lower bound for every multicast/compiler design."})
    equal(target["rows"],expected_rows,"incorrect full-family target count")
    equal(packet["summary"],{"episodes":len(expected_keys),"events":sum(len(r["events"]) for r in packet["episodes"]),"hops":sum(r["costs"]["hops"] for r in packet["episodes"])},"summary")
    return packet["summary"]


if __name__ == "__main__":
    print(json.dumps(verify(),sort_keys=True))

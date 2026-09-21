"""Reproduce the finite M1 controls and conditional production bounds."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


@contextmanager
def source_imports(**bindings):
    """Scope the inherited scripts' bare imports without polluting other tests."""
    absent = object()
    previous = {name: sys.modules.get(name, absent) for name in bindings}
    sys.modules.update(bindings)
    try:
        yield
    finally:
        for name, value in previous.items():
            if value is absent:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value


hierarchy = module("m1_hierarchy", HERE / "hierarchy.py")
storage = module("m1_storage", HERE / "storage.py")
codec = module("m1_tape_codec", ROOT / "code/source_read_routing/tape.py")
with source_imports(tape=codec):
    packing = module("m1_source_packing", ROOT / "code/source_read_routing/pack.py")
with source_imports(tape=codec, pack=packing):
    source_verifier = module("m1_source_verifier", ROOT / "code/source_read_routing/verify.py")
with source_imports(verify=source_verifier, pack=packing):
    oracle = module("m1_source_oracle", ROOT / "code/source_read_routing/check_control.py")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hierarchy_receipt():
    faces = hierarchy.BASE
    base_graph = hierarchy.graph(faces)
    all_paths = [hierarchy.paths(base_graph, root) for root in range(20)]
    diameter = max(len(path)-1 for paths in all_paths for path in paths)
    require(diameter == 3, "base diameter")
    for root, paths in enumerate(all_paths):
        for target, path in enumerate(paths):
            require(path[0] == root and path[-1] == target, "base endpoints")
            require(all(b in base_graph[a] for a,b in zip(path,path[1:])), "base edge")
    captures = {3: ROOT/"code/source_routing/support_w12_l3.json",
                4: ROOT/"code/source_read_routing/support_w12_l4.json",
                5: ROOT/"code/source_read_routing/support_w12_l5.json"}
    levels, refinements, pins = [], [], {}
    for level in range(6):
        adjacency = hierarchy.graph(faces)
        carriers = len(faces)
        require(carriers == 20*4**level, "carrier count")
        degree = max(map(len, adjacency))
        require(degree <= 12, "port capacity")
        seams = {(a,b) for a,row in enumerate(adjacency) for b in row if a < b}
        if level in captures:
            path = captures[level]
            packet = codec.load_json(path)
            require(list(map(list, faces)) == packet["faces"], "captured face mismatch")
            require(seams == {(a,b) for a,p,b,r in packet["glued_pairs"]}, "captured support mismatch")
            pins[path.relative_to(ROOT).as_posix()] = digest(path)
        levels.append({"level": level, "carriers": carriers, "glued_seams": len(seams),
                       "maximum_degree": degree, "proved_route_envelope": hierarchy.envelope(level)})
        if level < 5:
            fine = hierarchy.refine(faces)
            refinements.append(hierarchy.refinement_certificate(faces, fine))
            faces = fine
    return {"base_diameter": diameter, "levels": levels, "refinements": refinements,
            "captured_support_sha256": pins,
            "scope": "Analytic all-level bound assembled from four kernel-checked components; "
                     "no Lean theorem of the composed all-level graph identification. "
                     "Exact captured-face and seam comparison only at L3-L5."}


def control_receipt(write=False):
    source = ROOT / "code/source_read_routing/controls"
    target = HERE / "controls"
    result = []
    for variant in ("baseline", "source", "branch", "scratch"):
        path = source/f"q3_{variant}.json"
        packet = codec.load_json(path)
        # This oracle checks preparation, every original primitive, the full
        # menu, logical values, and every induced logical pair independently.
        with source_imports(verify=source_verifier):
            source_check = oracle.check(packet, source)
        rows = list(oracle.rows_for(packet, source))
        mapped, census = storage.replay(rows, packet["inputs"]["carriers"], packet["inputs"]["sites"])
        for name, field in (("events","events"), ("reads","register_reads"),
                            ("writes","register_writes"), ("logical_allocations","registers")):
            require(census[name] == packet["costs"][field], "source/target cost: "+name)
        tape = io.BytesIO()
        encoded = codec.encode(io.BytesIO(mapped), tape)
        require(encoded["decoded_sha256"] == census["mapped_decoded_sha256"], "target tape digest")
        destination = target/f"q3_{variant}.tape"
        if write:
            target.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(tape.getvalue())
        else:
            require(destination.read_bytes() == tape.getvalue(), "committed target tape drift")
        # Replay the bytes the reader receives, not just the pre-codec buffer.
        with destination.open("rb") as stream:
            retained = b"".join(codec.decode(stream))
        require(retained == mapped, "target codec round trip")
        compared = storage.check_refinement(rows, retained, packet["inputs"]["carriers"])
        physical = storage.replay_physical(retained, packet["inputs"]["carriers"])
        require(physical["events"] == census["events"] == compared, "target event count")
        require(physical["physical_slots_used"] == census["physical_slots_used"], "target slot census")
        result.append({"variant": variant, **census,
                       "source_receipt_sha256": digest(path),
                       "source_tape_decoded_sha256": packet["tape"]["decoded_sha256"],
                       "target_tape_sha256": hashlib.sha256(tape.getvalue()).hexdigest(),
                       "logical_pairs_checked": source_check["logical_pairs"],
                       "retained_events_replayed": physical["events"],
                       "source_target_events_compared": compared})
    return result


def production_bounds():
    result = []
    for q in (13,21):
        path = ROOT/f"evidence/source_net_causal_poset/routed_read_law/q{q}_baseline.json"
        packet = codec.load_json(path)
        inp, cost = packet["inputs"], packet["costs"]
        C,n,K,R,H = inp["carriers"], inp["sites"], inp["rounds"], inp["logical_reads"], cost["hops"]
        level = hierarchy.minimal_level(n)
        require(C == 20*4**level and n <= C < 4*n, "minimal support capacity")
        require(cost["events"] == 13*C+8*n+2*K*n+R+6*H, "event formula")
        require(cost["registers"] == 13*C+7*n+(K+1)*n+H, "archive formula")
        require(K*n <= R <= K*n*n and H <= K*n*(C-1), "resource envelope")
        require(cost["max_read_depth"] <= hierarchy.envelope(level), "route envelope")
        result.append({"q": q, "source_receipt_sha256": digest(path),
                       "carriers": C, "sites": n, "rounds": K,
                       "inherited_events": cost["events"], "logical_events": (K+1)*n,
                       "inherited_archived_registers": cost["registers"],
                       "derived_working_slot_upper_bound": 14*C+9*n,
                       "route_length_upper_bound": hierarchy.envelope(level),
                       "scope": "Bounds applied to pinned inherited counters; no new full production replay or measured target peak."})
    return result


def build_receipt(write=False):
    pins = [HERE/"hierarchy.py", HERE/"storage.py", Path(__file__),
            ROOT/"Lean/Geometry/SourceRoutingHierarchy.lean",
            ROOT/"Lean/Geometry/SourceRoutingStorage.lean",
            ROOT/"Lean/Geometry/SourceRoutingBudget.lean",
            ROOT/"Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean"]
    pins += [ROOT/"code/source_read_routing"/name for name in
             ("check_control.py", "verify.py", "pack.py", "tape.py", "produce.cpp")]
    return {"schema": "oph.m1_routing_refinement.v1",
            "source_sha256": {p.relative_to(ROOT).as_posix(): digest(p) for p in pins},
            "hierarchy": hierarchy_receipt(), "controls": control_receipt(write),
            "production_bound_applications": production_bounds(),
            "scope": {"M1_source_selected": False, "fixed_bit_memory_proved": False,
                      "physical_clock_identified": False, "operation_volume_measure_selected": False,
                      "arbitrary_retired_version_reads_supported": False,
                      "controller_and_audit_log_in_23_slot_bound": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate the reviewed receipt and complete small mapped tapes")
    args = parser.parse_args()
    receipt = build_receipt(args.write)
    destination = HERE/"receipt.json"
    if args.write:
        destination.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
    else:
        require(receipt == codec.load_json(destination), "receipt drift")
    print(json.dumps({"controls": len(receipt["controls"]),
                      "mapped_events": sum(r["events"] for r in receipt["controls"]),
                      "support_levels": len(receipt["hierarchy"]["levels"]), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()

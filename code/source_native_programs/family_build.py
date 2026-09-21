"""Regenerate compact family evidence; retain large witnesses only in a work directory."""
import argparse
import importlib.util
from math import isqrt
from pathlib import Path

from . import build, codec, routes


def menu(q):
    path = codec.ROOT/"evidence/source_net_causal_poset/build_causal_poset.py"
    spec = importlib.util.spec_from_file_location("native_program_golden_definition",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    points = module.site_coordinates(q,3)
    orbit = module.orbit(q)
    a,b = module.axis_tables(orbit)
    offsets,indices = module.site_graph(q,3,orbit,a,b,points)
    return [indices[offsets[i]:offsets[i+1]].tolist() for i in range(q**3)]


def resource_bound(q, rows, maximum_cells):
    n = q**3
    rounds = isqrt(q)+(isqrt(q)**2 < q)
    starts = rounds*n
    reads = rounds*sum(map(len,rows))
    incoming = [0]*n
    for row in rows:
        for s in row:
            incoming[s] += 1
    max_in,max_out = max(incoming),max(map(len,rows))
    # A complete input bank survives all reads in its layer. A source scale
    # grows at most max_in; a local accumulation grows by at most max_out.
    maximum_scale = starts+rounds*(max_in+max_out+2*maximum_cells+1)
    transfers = starts+reads
    residuals = starts+transfers
    depth = maximum_cells-1
    amplitude = n+1
    k = maximum_scale+build.ceil_log2(64*amplitude*residuals*depth**2)
    means = (starts*(3*k+3)+transfers*(2*maximum_cells+1+(2*maximum_cells-3)*depth**2*k)
             +reads*(3*maximum_scale+12))
    return {"q":q,"sites":n,"rounds":rounds,"complete_metric_reads":reads,
            "metric_menu_sha256":codec.digest(rows),"max_reads_per_site":max_out,
            "max_reads_of_record":max_in,"input_bound":amplitude,
            "maximum_route_cells":maximum_cells,"maximum_scale_upper":maximum_scale,
            "cleanup_blocks_sufficient":k,"scalar_means_upper":means,
            "grid_bits_sufficient":maximum_scale+build.ceil_log2(64*(means+1)),
            "physical_scalar_registers":12*20*4**({3:3,13:4,21:5}[q]),
            "scope":"analytic_finite_compilation_bound_not_executed_production_word"}


def case_summary(case):
    outputs = [next(o["decoded"] for o in checkpoint["observations"]
                    if o["version"] == [checkpoint["layer"],checkpoint["site"],-1])
               for checkpoint in case["checkpoints"] if checkpoint["kind"] == "store"]
    return {"payload":case["payload"],"history_sha256":codec.digest(case),
            "scalar_means":case["scalar_means"],"evaluated_means":case["evaluated_means"],
            "stationary_means":case["stationary_means"],"stored_values":outputs,
            "scalar_samples":2*sum(len(c["observations"]) for c in case["checkpoints"]),
            "plan_sha256":case["plan_sha256"],"final_state_sha256":case["final_state_sha256"]}


def generate(directory, reuse_routes=False):
    directory.mkdir(parents=True,exist_ok=True)
    bounds,route_summaries = [],[]
    from .check_routes import check,read_rows
    for q,level in ((3,3),(13,4),(21,5)):
        n = q**3
        path = directory/f"routes-{n}.jsonl"
        if not reuse_routes or not path.exists():
            with path.open("wb") as stream:
                for row in routes.generate(n,level):
                    stream.write((codec.compact(row)+"\n").encode("ascii"))
        summary = check(read_rows(path),n,level)
        route_summaries.append(summary)
        rows = menu(q)
        bounds.append(resource_bound(q,rows,summary["max_paired_cells"]))
        print(f"q={q}: complete routes and finite resource bound",flush=True)
        if q == 3:
            witnesses = list(read_rows(path))
            plan = build.compile_program([rows,rows],witnesses,28)
            (directory/"q3-plan.json").write_bytes(codec.canonical(plan))
            payloads = [list(range(1,28)),list(range(1,28))]
            payloads[1][13] += 1
            cases = [case_summary(build.execute(plan,p)) for p in payloads]
            exact_plan = {k:v for k,v in plan.items() if k not in ("segments","layers")}
    return {"schema":"oph-native-program-families-v1","pins":codec.pins(),
            "source_selected":False,"m1_derived":False,"routes":route_summaries,
            "bounds":bounds,"q3_plan":exact_plan,"q3_cases":cases,
            "production_native_replay":False,"q3_native_replay":True,
            "placement":"carrier_zero_core_banks_one_through_twice_sites"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir",type=Path,required=True)
    parser.add_argument("--output",type=Path,default=codec.HERE/"families.json")
    parser.add_argument("--reuse-routes",action="store_true")
    args = parser.parse_args()
    args.output.write_bytes(codec.canonical(generate(args.workdir,args.reuse_routes)))


if __name__ == "__main__":
    main()

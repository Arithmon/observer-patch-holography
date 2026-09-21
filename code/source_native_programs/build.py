"""Untrusted banked compiler and rounded native evaluator.

The supplied request lists determine the word without consulting payloads.
Only preparation and scalar pair means write physical state. Repeated words
may be evaluated at a proved integer fixed point; their full work is charged.
"""
from fractions import Fraction as F
import hashlib
from itertools import product
from math import ceil, floor

from . import codec, routes

LAYERS = [[[0,1],[0,0]], [[1,0],[]], [[0,1,0],[1]]]


def ceil_log2(n):
    return (n-1).bit_length()


def compile_program(layers, witnesses, input_bound=2):
    if type(input_bound) is not int or input_bound < 0:
        raise ValueError("nonnegative public integer input bound")
    if type(layers) is not list or not layers or type(layers[0]) is not list:
        raise ValueError("nonempty layered program")
    size = len(layers[0])
    if size < 1 or any(type(layer) is not list or len(layer) != size or
                       any(type(ss) is not list or any(type(s) is not int or not 0 <= s < size
                                                      for s in ss) for ss in layer) for layer in layers):
        raise ValueError("complete finite layer/request inventory")
    lookup = {(r["carrier"],r["direction"]):r["pairs"] for r in witnesses}
    bank = {i+1:{"pair":[12*(i+1),12*(i+1)+1],"scale":0,"version":[0,i,-1],"cap":input_bound}
            for i in range(size)}
    unit_scale = accumulator_scale = 0
    segments = []
    maximum_scale = 0
    cleanup_count = 0
    maximum_bus = 1

    def paired(a,b):
        return [[a[0],b[0]],[a[1],b[1]]]

    def segment(kind, layer, site, blocks, extra=()):
        nonlocal maximum_scale
        polls = [dict(row) for _,row in sorted(bank.items())]+list(extra)
        maximum_scale = max(maximum_scale,*(p["scale"] for p in polls))
        segments.append({"kind":kind,"layer":layer,"site":site,"blocks":blocks,"polls":polls})

    def block(edges, repeat=1, cleanup=False):
        return {"edges":edges,"repeat":repeat,"cleanup":cleanup}

    def shuttle(path, source, retired):
        nonlocal cleanup_count, maximum_bus
        cleanup_count += 1
        bus = path[:-1]
        maximum_bus = max(maximum_bus,len(bus))
        forward = sum((paired(a,b) for a,b in zip(path,path[1:])),[])
        sweep = [list(bus[0])]+sum((paired(a,b) for a,b in zip(bus,bus[1:])),[])
        initial = ([list(path[-1])] if retired else [])+paired(source,path[0])+forward
        return [block(initial),block(sweep,len(bus)**2,True)]

    A,W,U = [2,1],[4,9],[6,7]
    for layer, rows in enumerate(layers,1):
        old_base = 1+((layer-1)%2)*size
        destination_base = 1+(layer%2)*size
        for site, sources in enumerate(rows):
            cleanup_count += 1
            unit_scale += 1
            accumulator_scale = unit_scale
            cap = 1
            current = {"pair":A,"scale":accumulator_scale,"version":[layer,site,0],"cap":cap}
            seed = {"pair":U,"scale":unit_scale,"version":[-1,0,0],"cap":1}
            segment("start",layer,site,[block([W]+paired(A,W),1,True),
                                       block([W]+paired(U,A))],(current,seed))
            for position, source in enumerate(sources,1):
                carrier = old_base+source
                record = bank[carrier]
                path = lookup[carrier,"read"]
                if record["pair"][0] > record["pair"][1]:
                    path = [p[::-1] for p in path]
                operand = path[-1]
                operand_scale = record["scale"]+len(path)
                blocks = shuttle(path,record["pair"],False)
                record["scale"] += 1
                incoming = {"pair":operand,"scale":operand_scale,
                            "version":record["version"],"cap":record["cap"]}
                segment("read",layer,site,blocks,(current,seed,incoming))
                target = max(accumulator_scale,operand_scale+1)
                blocks = []
                for pair, repetitions in ((A,target-accumulator_scale),(operand,target-operand_scale-1)):
                    if repetitions:
                        blocks.append(block(paired(pair,W)+[W],repetitions))
                blocks.append(block(paired(operand,W)+paired(A,W)+[W]))
                cap += record["cap"]
                accumulator_scale = target+1
                current = {"pair":A,"scale":accumulator_scale,"version":[layer,site,position],"cap":cap}
                segment("add",layer,site,blocks,(current,seed))
                rectangle = paired(operand,W)+[[operand[0],W[1]],[operand[1],W[0]]]
                segment("retire_operand",layer,site,[block(rectangle)],(current,seed))
            destination = destination_base+site
            path = lookup[destination,"write"]
            blocks = shuttle(path,A,True)
            bank[destination] = {"pair":path[-1],"scale":accumulator_scale+len(path),
                                 "version":[layer,site,-1],"cap":cap}
            accumulator_scale += 1
            current = {**current,"scale":accumulator_scale}
            segment("store",layer,site,blocks,(current,seed))
    amplitude = max(1,input_bound)
    cleanup_blocks = maximum_scale+ceil_log2(64*amplitude*cleanup_count*maximum_bus**2)
    means = 0
    for stage in segments:
        for block in stage["blocks"]:
            if block["cleanup"]:
                block["repeat"] *= cleanup_blocks
            means += len(block["edges"])*block["repeat"]
    bits = maximum_scale+ceil_log2(64*(means+1))
    return {"sites":size,"input_bound":input_bound,"layers":layers,"segments":segments,"scalar_means":means,
            "maximum_scale":maximum_scale,"cleanup_stages":cleanup_count,
            "maximum_bus_cells":maximum_bus,"cleanup_blocks":cleanup_blocks,"grid_bits":bits,
            "residual_numerator":amplitude*cleanup_count*maximum_bus**2,
            "residual_denominator_exponent":cleanup_blocks}


def publish(plus,minus,exponent,cap,count,Q,residual):
    center = F((plus-minus)*(1 << exponent),2*(Q-3))
    radius = (F(4+5*count,8)+Q*residual)*(1 << exponent)/(Q-3)
    low,high = max(-cap,ceil(center-radius)),min(cap,floor(center+radius))
    return low if low == high else None


def execute(plan, payload, level=3, accelerate=True):
    n = plan["sites"]
    if len(payload) != n or any(type(z) is not int or abs(z) > plan["input_bound"] for z in payload):
        raise ValueError("prepared integer payload domain")
    Q = 1 << plan["grid_bits"]
    state = [(max(1,plan["input_bound"])+1)*Q]*(12*20*4**level)
    for pair,value in [([6,7],1)]+[([12*(i+1),12*(i+1)+1],z) for i,z in enumerate(payload)]:
        state[pair[0]] += (Q-3)*value
        state[pair[1]] -= (Q-3)*value
    preparation = codec.digest(state)
    count = evaluated = skipped = 0
    checkpoints = []
    history = hashlib.sha256()
    residual = F(plan["residual_numerator"],1 << plan["residual_denominator_exponent"])
    for stage in plan["segments"]:
        for block in stage["blocks"]:
            touched = sorted({v for edge in block["edges"] for v in edge})
            for repetition in range(block["repeat"]):
                before = [state[p] for p in touched]
                for u,v in block["edges"]:
                    total = state[u]+state[v]
                    midpoint = total//2+(total%4 == 3)
                    state[u] = state[v] = midpoint
                evaluated += len(block["edges"])
                if accelerate and before == [state[p] for p in touched]:
                    skipped += (block["repeat"]-repetition-1)*len(block["edges"])
                    break
            count += block["repeat"]*len(block["edges"])
            history.update((codec.compact([count,[[p,state[p]] for p in touched]])+"\n").encode("ascii"))
        observations = []
        for item in stage["polls"]:
            plus,minus = (state[p] for p in item["pair"])
            observations.append({**item,"samples":[plus,minus],
                "decoded":publish(plus,minus,item["scale"],item["cap"],count,Q,residual)})
        checkpoints.append({"kind":stage["kind"],"layer":stage["layer"],"site":stage["site"],
                            "scalar_means":count,"observations":observations})
    return {"payload":payload,"plan_sha256":codec.digest(plan),"preparation_sha256":preparation,
            "checkpoints":checkpoints,"scalar_means":count,"evaluated_means":evaluated,
            "stationary_means":skipped,"block_states_sha256":history.hexdigest(),
            "final_state_sha256":codec.digest(state)}


def controls():
    witnesses = list(routes.generate(2,3))
    plan = compile_program(LAYERS,witnesses)
    cases = []
    for payload in product((-2,0,2),repeat=2):
        case = execute(plan,list(payload))
        committed = [next(o["decoded"] for o in c["observations"]
                          if o["version"] == [c["layer"],c["site"],-1])
                     for c in case["checkpoints"] if c["kind"] == "store"]
        cases.append({"payload":case["payload"],"history_sha256":codec.digest(case),
                      "scalar_means":case["scalar_means"],"evaluated_means":case["evaluated_means"],
                      "stationary_means":case["stationary_means"],"stored_values":committed,
                      "final_state_sha256":case["final_state_sha256"]})
    return {"schema":"oph-native-stored-program-v1","pins":codec.pins(),
            "m1_derived":False,"source_selected":False,"routes":witnesses,
            "program":LAYERS,"plan":plan,"cases":cases,
            "scope":"finite_three_layer_two_bank_integer_program",
            "evaluation":"rounded_native_means_with_certified_stationary_repeat_acceleration"}


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=__import__("pathlib").Path,default=codec.HERE/"controls.json")
    args = parser.parse_args()
    args.output.write_bytes(codec.canonical(controls()))


if __name__ == "__main__":
    main()

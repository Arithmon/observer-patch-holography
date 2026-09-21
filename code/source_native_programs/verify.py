"""Independent schedule, lifetime, native replay and local publication verifier."""
from fractions import Fraction
import hashlib
from itertools import product

from . import codec
from .check_routes import check, edges, require


def certify_plan(plan, program, route_map):
    """Derive every block from requests and live-version metadata, not the producer.

    Bank writes are accepted only into the opposite bank. A layer reads the
    same committed input version throughout, including reads after other sites
    have published their outputs. Scale and cap calculations do not use values.
    """
    require(type(program) is list and program and type(program[0]) is list and program[0],
            "nonempty program")
    n = len(program[0])
    require(all(type(row) is list and len(row) == n and
                all(type(ss) is list and all(type(s) is int and 0 <= s < n for s in ss)
                    for ss in row) for row in program), "complete layer/request inventory")
    require(type(plan) is dict and type(plan.get("segments")) is list, "plan structure")
    k = plan.get("cleanup_blocks")
    require(type(k) is int and k >= 1, "positive cleanup precision")
    input_bound = plan.get("input_bound")
    require(type(input_bound) is int and input_bound >= 0, "public input bound")
    amplitude = max(1,input_bound)
    records = {i+1:((12*(i+1),12*(i+1)+1),0,0,i,input_bound) for i in range(n)}
    cursor = means = largest = residuals = 0
    depth = 1
    seed_exponent = 0

    def descriptor(pair, exponent, generation, site, prefix, cap):
        return {"pair":list(pair),"scale":exponent,"version":[generation,site,prefix],"cap":cap}

    def twice(left,right):
        return [[left[r],right[r]] for r in range(2)]

    def block(word, repetitions=1, cleaning=False):
        return {"edges":word,"repeat":repetitions,"cleanup":cleaning}

    def consume(kind, generation, site, blocks, additional):
        nonlocal cursor, means, largest
        observations = [descriptor(p,e,g,s,-1,c) for _,(p,e,g,s,c) in sorted(records.items())]
        observations.extend(additional)
        require(cursor < len(plan["segments"]), "truncated native program")
        expected = {"kind":kind,"layer":generation,"site":site,"blocks":blocks,"polls":observations}
        codec.equal(plan["segments"][cursor],expected,"native schedule/lifetime/scale/cap certificate")
        cursor += 1
        means += sum(len(b["edges"])*b["repeat"] for b in blocks)
        largest = max(largest,*(o["scale"] for o in observations))

    def transfer(path, source, clear_destination):
        nonlocal residuals, depth
        residuals += 1
        length = len(path)-1
        depth = max(depth,length)
        prefix = [list(path[-1])] if clear_destination else []
        prefix.extend(twice(source,path[0]))
        for i in range(1,len(path)):
            prefix.extend(twice(path[i-1],path[i]))
        sweep = [list(path[0])]
        for i in range(1,len(path)-1):
            sweep.extend(twice(path[i-1],path[i]))
        return [block(prefix),block(sweep,length*length*k,True)]

    acc, work, seed = (2,1),(4,9),(6,7)
    for generation, layer in enumerate(program,1):
        input_base = 1 if generation%2 else 1+n
        output_base = 1+n if generation%2 else 1
        for site, requests in enumerate(layer):
            residuals += 1
            seed_exponent += 1
            accumulator_exponent = seed_exponent
            cap = 1
            current = descriptor(acc,accumulator_exponent,generation,site,0,cap)
            unit = descriptor(seed,seed_exponent,-1,0,0,1)
            consume("start",generation,site,
                    [block([list(work)]+twice(acc,work),k,True),
                     block([list(work)]+twice(seed,acc))],[current,unit])
            for position, source in enumerate(requests,1):
                carrier = input_base+source
                pair,exponent,version,address,limit = records[carrier]
                require(version == generation-1 and address == source, "stale or overwritten input version")
                path = route_map[carrier,"read"]
                if pair[0] != 12*carrier:
                    path = [[b,a] for a,b in path]
                operand = path[-1]
                input_exponent = exponent+len(path)
                records[carrier] = (pair,exponent+1,version,address,limit)
                incoming = descriptor(operand,input_exponent,version,address,-1,limit)
                consume("read",generation,site,transfer(path,pair,False),[current,unit,incoming])
                common = max(accumulator_exponent,input_exponent+1)
                additions = []
                for pair,repeat in ((acc,common-accumulator_exponent),
                                    (operand,common-input_exponent-1)):
                    if repeat > 0:
                        additions.append(block(twice(pair,work)+[list(work)],repeat))
                additions.append(block(twice(operand,work)+twice(acc,work)+[list(work)]))
                accumulator_exponent,cap = common+1,cap+limit
                current = descriptor(acc,accumulator_exponent,generation,site,position,cap)
                consume("add",generation,site,additions,[current,unit])
                retirement = twice(operand,work)+[[operand[0],work[1]],[operand[1],work[0]]]
                consume("retire_operand",generation,site,[block(retirement)],[current,unit])
            destination = output_base+site
            require(destination not in range(input_base,input_base+n), "overwrite of live input bank")
            path = route_map[destination,"write"]
            records[destination] = (tuple(path[-1]),accumulator_exponent+len(path),generation,site,cap)
            current = {**current,"scale":accumulator_exponent+1}
            consume("store",generation,site,transfer(path,acc,True),[current,unit])
    require(cursor == len(plan["segments"]), "extra native program segment")
    want_k = largest+(64*amplitude*residuals*depth*depth-1).bit_length()
    require(k == want_k, "cleanup budget does not cover whole history")
    expected = {"sites":n,"input_bound":input_bound,"layers":program,"segments":plan["segments"],"scalar_means":means,
                "maximum_scale":largest,"cleanup_stages":residuals,"maximum_bus_cells":depth,
                "cleanup_blocks":k,"grid_bits":largest+(64*(means+1)-1).bit_length(),
                "residual_numerator":amplitude*residuals*depth*depth,"residual_denominator_exponent":k}
    codec.equal(plan,expected,"finite resource/precision certificate")
    return expected


def local_decode(plus,minus,shift,cap,steps,grid,residual):
    observation = Fraction(plus-minus,2*(grid-3))*(2**shift)
    allowance = (Fraction(4+5*steps,8)+grid*residual)*Fraction(2**shift,grid-3)
    left,right = observation-allowance,observation+allowance
    first = max(-cap,-((-left.numerator)//left.denominator))
    last = min(cap,right.numerator//right.denominator)
    return (first if first == last else None),observation,allowance


def replay(plan, payload, allowed, accelerate=True, preparation_error=False, level=3, disturbed=False):
    n = plan["sites"]
    require(type(payload) is list and len(payload) == n and
            all(type(z) is int and abs(z) <= plan["input_bound"] for z in payload), "payload domain")
    grid = 2**plan["grid_bits"]
    ports = 12*20*4**level
    baseline = (max(1,plan["input_bound"])+1)*grid
    cells = {p:baseline for p in range(ports)}
    cells[6],cells[7] = baseline+grid-3,baseline-grid+3
    for i,value in enumerate(payload,1):
        cells[12*i],cells[12*i+1] = baseline+(grid-3)*value,baseline-(grid-3)*value
    initial_hash = codec.digest([cells[i] for i in range(ports)])
    if preparation_error:
        cells = {p:value+Fraction(1 if p%2 else -1,4) for p,value in cells.items()}
    require(not (preparation_error and disturbed),"choose one disturbance control")
    phase = 0
    if disturbed:
        # Eighth-grid integer coordinates. A lazy alternating offset applies
        # +/-1/8 to every port after each mean, including all idle records.
        cells = {p:8*value+(2 if p%2 else -2) for p,value in cells.items()}
    # The logical oracle is confined to assertions. It never writes native state
    # and is not an argument of the local interval comparator.
    logical = [payload]
    for layer in plan["layers"]:
        previous = logical[-1]
        logical.append([1+sum(previous[source] for source in requests) for requests in layer])
    def answer(version):
        layer,site,prefix = version
        if layer == -1:
            return 1
        if prefix == -1:
            return logical[layer][site]
        return 1+sum(logical[layer-1][s] for s in plan["layers"][layer-1][site][:prefix])

    residual = Fraction(plan["residual_numerator"],2**plan["residual_denominator_exponent"])
    steps = computed = stationary = 0
    transcript = hashlib.sha256()
    checkpoints = []
    for stage in plan["segments"]:
        for group in stage["blocks"]:
            word = group["edges"]
            require(all(type(e) is list and len(e) == 2 and e[0] != e[1] and
                        tuple(sorted(e)) in allowed for e in word), "unsupported scalar operation")
            visited = sorted({p for edge in word for p in edge})
            remaining = group["repeat"]
            require(type(remaining) is int and remaining > 0, "invalid repetition count")
            while remaining:
                snapshot = {p:cells[p] for p in visited}
                old_phase = phase
                for left,right in word:
                    if disturbed:
                        dl = phase*(1 if left%2 else -1)
                        dr = phase*(1 if right%2 else -1)
                        quotient,remainder = divmod(cells[left]+dl+cells[right]+dr,16)
                        value = 8*(quotient+(remainder>8 or (remainder == 8 and quotient%2)))
                        cells[left],cells[right] = value-dl,value-dr
                        phase = 1-phase
                    else:
                        quotient,remainder = divmod(cells[left]+cells[right],2)
                        value = (int(quotient)+(remainder == 1 and int(quotient)%2)
                                 if not preparation_error else round(Fraction(cells[left]+cells[right],2)))
                        cells[left] = value
                        cells[right] = value
                computed += len(word)
                remaining -= 1
                if accelerate and phase == old_phase and all(cells[p] == snapshot[p] for p in visited):
                    stationary += remaining*len(word)
                    remaining = 0
            steps += group["repeat"]*len(word)
            transcript.update((codec.compact([steps,[[p,cells[p]] for p in visited]])+"\n").encode("ascii"))
        observations = []
        for record in stage["polls"]:
            plus,minus = [cells[p] for p in record["pair"]]
            if disturbed:
                l,r = record["pair"]
                plus = Fraction(plus+phase*(1 if l%2 else -1),8)+Fraction(1,4)
                minus = Fraction(minus+phase*(1 if r%2 else -1),8)-Fraction(1,4)
            if preparation_error:
                plus,minus = plus+Fraction(1,4),minus-Fraction(1,4)
            result,center,radius = local_decode(plus,minus,record["scale"],record["cap"],steps,grid,residual)
            expected = answer(record["version"])
            require(abs(expected) <= record["cap"], "public cap does not contain recurrence")
            require(abs(center-expected) <= radius, "native state violates global error bound")
            require(2*radius < 1, "precision cannot identify record")
            require(result == expected, "incorrect local publication or lost committed version")
            observations.append({**record,"samples":[plus,minus],"decoded":result})
        checkpoints.append({"kind":stage["kind"],"layer":stage["layer"],"site":stage["site"],
                            "scalar_means":steps,"observations":observations})
    require(steps == plan["scalar_means"] and computed+stationary == steps, "unpaid native work")
    if preparation_error or disturbed:
        return {"scalar_means":steps,"checkpoints":len(checkpoints),"published":logical[-1]}
    return {"payload":payload,"plan_sha256":codec.digest(plan),"preparation_sha256":initial_hash,
            "checkpoints":checkpoints,"scalar_means":steps,"evaluated_means":computed,
            "stationary_means":stationary,"block_states_sha256":transcript.hexdigest(),
            "final_state_sha256":codec.digest([cells[i] for i in range(ports)])}


def verify(packet):
    require(type(packet) is dict and set(packet) == {
        "schema","pins","m1_derived","source_selected","routes","program","plan","cases","scope","evaluation"},
        "artifact fields")
    codec.equal(packet["schema"],"oph-native-stored-program-v1","schema")
    codec.equal(packet["m1_derived"],False,"unproved M1 promotion")
    codec.equal(packet["source_selected"],False,"unproved source selection")
    codec.equal(packet["scope"],"finite_three_layer_two_bank_integer_program","scope")
    codec.equal(packet["evaluation"],"rounded_native_means_with_certified_stationary_repeat_acceleration","evaluation")
    pins = codec.pins()
    codec.equal(packet["pins"],pins,"source/proof custody")
    route_summary, route_map = check(packet["routes"],2,3,retain=True)
    program = [[[0,1],[0,0]], [[1,0],[]], [[0,1,0],[1]]]
    codec.equal(packet["program"],program,"complete control program")
    plan = certify_plan(packet["plan"],program,route_map)
    require(plan["input_bound"] == 2,"signed control bound")
    require(type(packet["cases"]) is list and len(packet["cases"]) == 9,"signed payload inventory")
    allowed = edges(3)
    totals = {"scalar_means":0,"evaluated_means":0,"stationary_means":0,"scalar_samples":0}
    for case,values in zip(packet["cases"],product((-2,0,2),repeat=2)):
        require(type(case) is dict,"case structure")
        codec.equal(case.get("payload"),list(values),"complete signed payload inventory")
        result = replay(plan,list(values),allowed)
        committed = []
        for checkpoint in result["checkpoints"]:
            if checkpoint["kind"] == "store":
                version = [checkpoint["layer"],checkpoint["site"],-1]
                outputs = [o["decoded"] for o in checkpoint["observations"] if o["version"] == version]
                require(len(outputs) == 1,"one committed output per site/layer")
                committed.append(outputs[0])
        expected = {"payload":result["payload"],"history_sha256":codec.digest(result),
                    "scalar_means":result["scalar_means"],"evaluated_means":result["evaluated_means"],
                    "stationary_means":result["stationary_means"],"stored_values":committed,
                    "final_state_sha256":result["final_state_sha256"]}
        codec.equal(case,expected,"native state, version, observation or work differs")
        for key in ("scalar_means","evaluated_means","stationary_means"):
            totals[key] += result[key]
        totals["scalar_samples"] += 2*sum(len(c["observations"]) for c in result["checkpoints"])
    codec.equal(codec.pins(),pins,"sources changed during verification")
    return {"schema":"oph-native-stored-program-receipt-v1","pins":pins,
            "controls_sha256":codec.digest(packet),"histories":9,"layers_per_history":3,
            "sites":2,"route_certificate":route_summary,**totals,
            "grid_bits":plan["grid_bits"],"maximum_scale":plan["maximum_scale"],
            "preparation_writes_per_history":12*1280,
            "dynamic_writes":"supported_scalar_pair_means_only",
            "source_selected":False,"m1_derived":False,
            "large_family_scope":"route_witnesses_only_no_native_production_replay",
            "supplied":["initial records and unit","placement and protected lifetime",
                        "layer requests and temporal controller","precision and readout instrument"]}


def main():
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls",type=Path,default=codec.HERE/"controls.json")
    parser.add_argument("--receipt",type=Path,default=codec.HERE/"receipt.json")
    parser.add_argument("--write-receipt",action="store_true")
    args = parser.parse_args()
    receipt = verify(codec.load_artifact(args.controls))
    if args.write_receipt:
        args.receipt.write_bytes(codec.canonical(receipt))
    else:
        codec.equal(codec.load_artifact(args.receipt),receipt,"receipt mismatch")
    print(f"Verified {receipt['histories']} stored programs; {receipt['scalar_means']} charged native means.")


if __name__ == "__main__":
    main()

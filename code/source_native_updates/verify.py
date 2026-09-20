"""Independent semantic replay; never import the producer or trust its digests."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
if __package__:
    from . import codec
else:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_native_updates import codec


def require(condition, message):
    if not condition:
        raise ValueError(message)


def operation(k):
    if k < 7:
        return ((0,6),(1,7),(2,4),(3,5),(6,4),(7,5),(6,7))[k]
    cycle, phase = divmod(k-7, 568)
    if phase < 8:
        hop, rail = divmod(phase, 2)
        source = (2,0,2,1)[cycle] if hop == 0 else hop+2
        return 2*source+rail, 2*(hop+3)+rail
    phase = (phase-8) % 7
    if phase == 0:
        return 6,7
    hop, rail = divmod(phase-1, 2)
    return 6+2*hop+rail, 8+2*hop+rail


def support():
    source = codec.load(codec.ROOT/codec.SUPPORT)
    edges = {tuple(sorted((12*c+u,12*c+v))) for c in range(source["carriers"])
             for u,v in source["intra_carrier_seams"]}
    edges.update(tuple(sorted((12*c+u,12*d+v))) for c,u,d,v in source["glued_pairs"])
    return edges, 12*source["carriers"]


def layout(ports):
    codec.equal(ports, [0,9,7,4,8,11,1,5,14,40,23,38,45,39], "changed embedding")
    edges, count = support()
    require(len(set(ports)) == 14 and all(0 <= p < count for p in ports), "aliased ports")
    for k in range(2279):
        u,v = operation(k)
        require(tuple(sorted((ports[u],ports[v]))) in edges, "unsupported mean")
    return edges, count


@lru_cache(maxsize=1)
def ideal_reads():
    rows = [(F(0),F(0)) for _ in range(14)]
    rows[:4] = [(F(1),F(0)), (F(-1),F(0)), (F(0),F(1)), (F(0),F(-1))]
    reads = []
    for k in range(2279):
        u,v = operation(k)
        rows[u] = rows[v] = tuple((x+y)/2 for x,y in zip(rows[u],rows[v]))
        if k == 6:
            codec.equal([[str(x) for x in row] for row in rows[:6]],
                        [["1/2","0"],["-1/2","0"],["0","1/2"],["0","-1/2"],
                         ["1/4","1/4"],["-1/4","-1/4"]], "native commit formula")
            require(all(row == (0,0) for row in rows[6:]), "nonblank initial bus")
        if k >= 7 and (k-7) % 568 == 7:
            reads.append(tuple((p-m)/2 for p,m in zip(rows[12],rows[13])))
        if k >= 7 and (k-7) % 568 == 567:
            require(all(F(3,2)*sum(abs(x) for x in row) <= F(3,2)*F(7,8)**80
                        for row in rows[6:]), "cleanup residual")
    return tuple(reads)


def bounds(cycle, coefficients, grid=2**20):
    require(type(grid) is int and grid > 0, "invalid grid")
    record, use, exponent = ((2,1,6),(0,1,5),(2,2,7),(1,1,5))[cycle]
    gain = F(1,2**exponent)
    target = (gain,gain) if record == 2 else ((gain,F(0)) if record == 0 else (F(0),gain))
    actual = F(3,2)*sum(abs(x-y) for x,y in zip(coefficients,target))
    cleanup = F(3,2)*F(7,8)**80
    require(actual <= cleanup, "universal transfer estimate")
    count = 15+568*cycle
    implementation = F(1,2**18)+count*(F(1,2*grid)+F(1,2**28))+F(1,2**16)
    total = cleanup+implementation
    return {"record": record, "use": use, "exponent": exponent,
            "means_before_read": count, "coefficients": list(map(str,coefficients)),
            "ideal_error": str(actual), "cleanup_bound": str(cleanup),
            "implementation_bound": str(implementation), "half_spacing": str(gain/2),
            "total_error": str(total), "margin": str(gain/2-total)}


def check_case(case, payloads, maps, retained=None):
    q = 2**20
    a,b = payloads
    state = [2*q]*14
    for j, amplitude in enumerate((a,b)):
        state[2*j] += int(amplitude*q)
        state[2*j+1] -= int(amplitude*q)
    initial = state[:]
    writers = [-j-1 for j in range(14)]
    ancestry = [1 << j for j in range(14)]
    tape, reads, audit = [], [], []
    max_rounding = F(0)
    for k in range(2279):
        u,v = operation(k)
        total = state[u]+state[v]
        floor, remainder = divmod(total,2)
        output = floor+(remainder if floor % 2 else 0)
        event = [u,v,state[u],state[v],writers[u],writers[v],output]
        tape.append(event)
        if retained is not None:
            require(len(retained) == 2279, "truncated/extra tape")
            codec.equal(retained[k], event, "changed mean/input/output/writer")
        error = abs(F(output)-F(total,2))/q
        require(error <= F(1,2*q), "rounding bound")
        max_rounding = max(max_rounding, error)
        require(0 <= output <= 4*q, "overflow")
        lineage = ancestry[u] | ancestry[v]
        state[u] = state[v] = output
        writers[u] = writers[v] = k
        ancestry[u] = ancestry[v] = lineage
        if k == 6:
            committed = state[:]
            require(F(state[4]-state[5],2*q) == (a+b)/4, "sum not natively written")
            require(F(state[0]-state[1],2*q) == a/2 and
                    F(state[2]-state[3],2*q) == b/2, "old record destroyed at commit")
            require(ancestry[4] | ancestry[5] == (1 << 8)-1, "commit consumed wrong records")
        if k >= 7 and (k-7) % 568 == 7:
            cycle = (k-7)//568
            record,use,exponent = ((2,1,6),(0,1,5),(2,2,7),(1,1,5))[cycle]
            raw = F(state[12]-state[13],2*q)
            scaled = raw*2**exponent
            codebook = [F(v,2) for v in (-3,-1,1,3)] if record < 2 else list(map(F,range(-3,4)))
            # Independent decision rule: open Voronoi cells with saturation at ends.
            cuts = [(x+y)/2 for x,y in zip(codebook,codebook[1:])]
            require(scaled not in cuts, "ambiguous receiver")
            decoded = codebook[sum(scaled > cut for cut in cuts)]
            require(decoded == (a,b,a+b)[record], "incorrect decoded version")
            ideal = sum(x*y for x,y in zip(maps[cycle],payloads))
            require(abs(raw-ideal) <= F(k+1,2*q), "composed rounding error")
            bound = bounds(cycle,maps[cycle])
            require(codec.rational(bound["margin"]) > 0, "nonpositive margin")
            reads.append({"after_mean": k, "record": record, "use": use,
                          "exponent": exponent, "local_units": state[12:14],
                          "writers": writers[12:14], "decoded": str(decoded)})
            origins = ancestry[12] | ancestry[13]
            require(origins & 15 == 15, "raw cross-record ancestry omitted")
            audit.append({"contrast": str(raw), "rounding_error": str(abs(raw-ideal)),
                          "initial_ancestors": [j for j in range(14) if origins & (1 << j)]})
    expected = {"payloads": list(map(str,payloads)), "initial_units": initial,
                "commit_units": committed, "reads": reads, "final_units": state,
                "event_count": len(tape), "tape_sha256": codec.digest(tape)}
    codec.equal(case, expected, "case does not match independently replayed history")
    return {"payloads": list(map(str,payloads)), "reads": audit,
            "max_rounding_error": str(max_rounding)}, tape


def scaling_controls():
    certificates = []
    for n in (1,2,3,8,16,32,64):
        mass = (n+1)**2
        weights = [F(mass-(n-j)**2) for j in range(n+1)]
        rate = 1-F(1,mass)
        values = weights[:]
        values[0] = 0
        for j in range(n):
            values[j] = values[j+1] = (values[j]+values[j+1])/2
        require(all(x <= rate*w for x,w in zip(values,weights)), "weighted contraction")
        require(min(weights) == 2*n+1 and max(weights) == mass, "weight condition number")
        certificates.append({"edges": n, "rate": str(rate), "weight_min": 2*n+1,
                             "weight_max": mass, "means_per_sweep": 1+2*n,
                             "sweeps_per_binary_block": mass,
                             "means_per_binary_block": (1+2*n)*mass,
                             "weighted_image": list(map(str,values))})
    return certificates


def schedule_control(edges):
    result = []
    for order in (((0,1),(1,14)), ((1,14),(0,1))):
        values = {0:F(1),1:F(3),14:F(5)}
        origins = {i:{i} for i in values}
        tape = []
        for u,v in order:
            require(tuple(sorted((u,v))) in edges, "unsupported schedule control")
            before = sum(x*x for x in values.values())
            x,y = values[u],values[v]
            values[u] = values[v] = (x+y)/2
            origins[u] = origins[v] = origins[u] | origins[v]
            after = sum(x*x for x in values.values())
            require(after < before and before-after == (x-y)**2/2, "non-strict repair")
            tape.append([u,v,str(x),str(y),str(values[u])])
        result.append({"tape": tape, "read_port": 0, "value": str(values[0]),
                       "initial_ancestors": sorted(origins[0])})
    require(result[0]["value"] == "2" and result[1]["value"] == "5/2", "schedule witness")
    return result


def verify(artifact):
    codec.equal(sorted(artifact), sorted(("schema","pins","ports","cases",
                "representative_payloads","representative_tape")), "artifact fields")
    codec.equal(artifact["schema"], "oph-native-updates-v1", "schema")
    codec.equal(artifact["pins"], codec.pins(), "source pins")
    edges, vertices = layout(artifact["ports"])
    codec.equal(artifact["representative_payloads"], ["1/2","-3/2"], "representative selection")
    require(type(artifact["representative_tape"]) is list and
            len(artifact["representative_tape"]) == 2279, "missing/truncated representative tape")
    maps = ideal_reads()
    require(type(artifact["cases"]) is list and len(artifact["cases"]) == 16, "missing/extra cases")
    audits = []
    levels = [F(v,2) for v in (-3,-1,1,3)]
    for case, payloads in zip(artifact["cases"],product(levels, repeat=2)):
        retained = artifact["representative_tape"] if payloads == (F(1,2),F(-3,2)) else None
        audits.append(check_case(case,payloads,maps,retained)[0])
    # Intervention comparisons are on decoded records, not on raw linear maps.
    for left,right in product(artifact["cases"],repeat=2):
        a,b = map(codec.rational,left["payloads"])
        c,d = map(codec.rational,right["payloads"])
        for position, unchanged in enumerate((a+b == c+d,a == c,a+b == c+d,b == d)):
            if unchanged:
                require(left["reads"][position]["decoded"] == right["reads"][position]["decoded"],
                        "decoded intervention correspondence")
    ports = artifact["ports"]
    resources = {"support_ports": vertices, "active_scalar_registers": 14,
                 "logical_records_after_commit": 3, "native_commit_means": 7,
                 "means_per_history": 2279, "histories": 16, "means_replayed": 16*2279,
                 "preparation_writes": 14, "native_scalar_reads": 2*2279,
                 "native_scalar_writes": 2*2279, "receiver_samples": 8,
                 "cross_carrier_means": sum(ports[operation(k)[0]]//12 != ports[operation(k)[1]]//12
                                              for k in range(2279)),
                 "scalar_storage_bits_bound": 14*(4*2**20).bit_length(),
                 "mean_adder_bits_bound": (8*2**20).bit_length(),
                 "receiver_scaled_difference_bits_bound": (4*2**20*2**7).bit_length()+1,
                 "expanded_schedule_bits": 2*4*2279,
                 "global_port_map_bits": 14*(vertices-1).bit_length(),
                 "program_counter_bits": (2279).bit_length(),
                 "record_scale_bits": 3*2, "read_use_counter_bits": 3*2,
                 "record_identifier_bits": 3*2, "commit_flag_bits": 1,
                 "codebook_numerator_bits": 4*3+7*4,
                 "request_schedule_bits": 4*(2+2+(2279).bit_length())}
    read_bounds = [bounds(i,row) for i,row in enumerate(maps)]
    insufficient = [bounds(i,row,2**8)["margin"] for i,row in enumerate(maps)]
    require(all(codec.rational(m) < 0 for m in insufficient), "precision control unexpectedly passes")
    return {"schema": "oph-native-updates-verification-v1", "controls_sha256": codec.digest(artifact),
            "resources": resources, "read_bounds": read_bounds, "histories": audits,
            "chain_scaling": scaling_controls(), "schedule_control": schedule_control(edges),
            "insufficient_grid_margins": insufficient}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls", type=Path, default=codec.HERE/"controls.json")
    parser.add_argument("--receipt", type=Path, default=codec.HERE/"receipt.json")
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = verify(codec.load_artifact(args.controls))
    if args.write_receipt:
        args.receipt.write_bytes(codec.canonical(result))
    else:
        codec.equal(codec.load_artifact(args.receipt), result, "verification receipt mismatch")
    print("Verified native commit and 64 reads across 16 complete captured histories.")


if __name__ == "__main__":
    main()

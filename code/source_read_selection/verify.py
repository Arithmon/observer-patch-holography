"""Independent exact census, coefficient replay and probability verification."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib

if __package__:
    from . import codec,verify_constraints
else:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_read_selection import codec,verify_constraints


def require(ok, label):
    if not ok:
        raise ValueError(label)


def support():
    raw = codec.load(codec.ROOT/codec.SUPPORT)
    require(raw["schema"] == "oph.source_routing.w12_support.v1", "support schema")
    require(raw["carriers"] == 1280 and raw["level"] == 3, "support regulator")
    adjacency = {i:set() for i in range(15360)}
    inserted = 0
    def add(u,v):
        nonlocal inserted
        require(type(u) is int and type(v) is int and u != v, "invalid edge")
        require(u in adjacency and v in adjacency, "edge range")
        require(v not in adjacency[u], "duplicate primitive seam")
        adjacency[u].add(v)
        adjacency[v].add(u)
        inserted += 1
    for carrier in range(1280):
        for pair in raw["intra_carrier_seams"]:
            require(len(pair)==2 and all(type(p) is int and 0<=p<12 for p in pair), "local port")
            add(carrier*12+pair[0], carrier*12+pair[1])
    for c,u,d,v in raw["glued_pairs"]:
        require(all(type(p) is int and 0<=p<12 for p in (u,v)), "glued port")
        require(type(c) is int and type(d) is int and c != d, "glued carrier")
        add(c*12+u,d*12+v)
    alphabet = [[u,v] for u in adjacency for v in sorted(adjacency[u]) if u<v]
    require(len(alphabet)==inserted==46050, "captured alphabet census")
    cut = [[u,v] for u in range(12) for v in sorted(adjacency[u]) if v>=12]
    require(len(cut)==11, "captured cut census")
    return alphabet,cut


def probability_row(total, cut, h):
    # Integer powers and Euclidean division independently enclose the exact
    # avoidance mass. No logarithm, binary64 approximation or Monte Carlo.
    numerator,denominator = pow(total-cut,h),pow(total,h)
    units,remainder = divmod(numerator*(1<<32),denominator)
    require(units*denominator <= numerator*(1<<32) <= (units+(remainder!=0))*denominator,
            "invalid rational enclosure")
    return {"horizon":h,"avoid_mass":{"base":str(F(total-cut,total)),"exponent":h},
            "denominator":1<<32,"lower_units":units,"upper_units":units+(remainder!=0),
            "balanced_error_lower":str(F(units,1<<33))}


def toy(h):
    error_counts = [0,0]
    avoid_count = means = 0
    histories = hashlib.sha256()
    for number in range(3**h):
        remainder = number
        word = [0]*h
        for place in range(h-1,-1,-1):
            remainder,word[place] = divmod(remainder,3)
        safe = all(move != 0 for move in word)
        avoid_count += safe
        # Propagate the coefficient of the sole source variable, instead of
        # executing either input state as the producer does.
        coefficients = [F(1),F(0),F(0),F(0)]
        writers = [-1,-2,-3,-4]
        symbolic = []
        for event,u in enumerate(word):
            v = u+1
            average = coefficients[u]/2+coefficients[v]/2
            symbolic.append((u,v,coefficients[u],coefficients[v],writers[u],writers[v],average))
            coefficients[u] = coefficients[v] = average
            writers[u] = writers[v] = event
            if safe:
                require(all(c==0 for c in coefficients[1:]), "cut-prefix locality violation")
        require(not safe or coefficients[3]==0, "cut locality violation")
        for index,payload in enumerate((-1,1)):
            tape = [[u,v,str(a*payload),str(b*payload),wu,wv,str(out*payload)]
                    for u,v,a,b,wu,wv,out in symbolic]
            value = coefficients[3]*payload
            decoded = 0 if value==0 else (1 if value>0 else -1)
            error_counts[index] += decoded != payload
            means += h
            histories.update((codec.compact([word,payload,tape,[str(value)],decoded])+"\n").encode())
    require(avoid_count==2**h, "subalphabet word count")
    require(sum(error_counts)>=avoid_count, "two-input error bound")
    return {"horizon":h,"words":3**h,"avoid_words":avoid_count,"failure_counts":error_counts,
            "scalar_means":means,"scalar_reads":2*means,"scalar_writes":2*means,
            "histories_sha256":histories.hexdigest()}


def control(name, alphabet):
    if name=="routed":
        cells = [(0,9),(1,5),(14,40),(23,38),(45,39)]
        word = [[cells[i][rail],cells[i+1][rail]] for i in range(4) for rail in range(2)]
    else:
        require(name=="cut_avoiding", "control name")
        word = [[0,1],[14,23],[23,45],[38,39]]*4
    ports = [0,1,5,9,14,23,38,39,40,45]
    allowed = {tuple(e) for e in alphabet}
    require(all(tuple(sorted(e)) in allowed for e in word), "unsupported control")
    if name=="cut_avoiding":
        require(all((u<12)==(v<12) for u,v in word), "control crosses cut")
    coefficient = dict.fromkeys(ports,F(0))
    coefficient[0],coefficient[9] = F(1),F(-1)
    writers = {p:-p-1 for p in ports}
    symbolic = []
    for event,(u,v) in enumerate(word):
        mean = (coefficient[u]+coefficient[v])*F(1,2)
        symbolic.append((u,v,coefficient[u],coefficient[v],writers[u],writers[v],mean))
        coefficient[u] = coefficient[v] = mean
        writers[u] = writers[v] = event
        if name=="cut_avoiding":
            require(all(c==0 for p,c in coefficient.items() if p>=12), "exterior transcript changed")
    cases = []
    for payload in (-1,1):
        tape = [[u,v,str(2+a*payload),str(2+b*payload),wu,wv,str(2+out*payload)]
                for u,v,a,b,wu,wv,out in symbolic]
        local = [2+coefficient[r]*payload for r in (45,39)]
        decoded = 0 if local[0]==local[1] else (1 if local[0]>local[1] else -1)
        require(decoded==(payload if name=="routed" else 0), "control decoder")
        if name=="routed":
            require(local[0]-local[1]==F(payload,8), "routed gain")
        initial = [[p,str(F(2)+(payload if p==0 else -payload if p==9 else 0))] for p in ports]
        cases.append({"payload":payload,"initial":initial,"local":list(map(str,local)),
                      "decoded":decoded,"tape":tape,"tape_sha256":codec.digest(tape)})
    return {"name":name,"word":word,"cases":cases,"means_per_case":len(word),
            "reads_per_case":2*len(word),"writes_per_case":2*len(word)}


def verify(artifact):
    require(type(artifact) is dict, "artifact object")
    codec.equal(artifact.get("pins"),codec.pins(),"source pins")
    projection_receipt = verify_constraints.verify(artifact["projection_controls"])
    alphabet,cut = support()
    expected = {"schema":"oph-read-selection-v1","pins":codec.pins(),
                "scope":"conditional_classical_KL_constraints_cover_and_read_obstruction",
                "alphabet":{"identity":"distinct_undirected_scalar_seams","size":len(alphabet),
                            "sha256":codec.digest(alphabet)},
                "source_ports":list(range(12)),"receiver_ports":[45,39],"cut":cut,
                "probabilities":[probability_row(len(alphabet),len(cut),h) for h in (0,1,2,16,256,2279)],
                "toy":{"alphabet":[[0,1],[1,2],[2,3]],"source":0,"receiver":3,
                       "payloads":[-1,1],"decoder":"sign_with_zero_failure",
                       "rows":[toy(h) for h in range(7)]},
                "captured_controls":[control(name,alphabet) for name in ("routed","cut_avoiding")],
                "projection_controls":artifact["projection_controls"]}
    codec.equal(artifact,expected,"independently reconstructed evidence")
    means = sum(row["scalar_means"] for row in expected["toy"]["rows"])
    control_means = sum(2*c["means_per_case"] for c in expected["captured_controls"])
    return {"schema":"oph-read-selection-receipt-v1","controls_sha256":codec.digest(artifact),
            "seams":len(alphabet),"source_cut":len(cut),"noncrossing_seams":len(alphabet)-len(cut),
            "toy_histories":sum(2*row["words"] for row in expected["toy"]["rows"]),
            "toy_scalar_means":means,"captured_scalar_means":control_means,
            "toy_preparation_writes":4*sum(2*row["words"] for row in expected["toy"]["rows"]),
            "toy_receiver_samples":sum(2*row["words"] for row in expected["toy"]["rows"]),
            "captured_preparation_writes":40,"captured_receiver_samples":8,
            "constrained_projection":projection_receipt,
            "probabilities":expected["probabilities"],
            "status":"conditional_obstruction_not_full_axiom_countermodel"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls",type=Path,default=codec.HERE/"controls.json")
    parser.add_argument("--receipt",type=Path,default=codec.HERE/"receipt.json")
    parser.add_argument("--write-receipt",action="store_true")
    args = parser.parse_args()
    try:
        result = verify(codec.load_artifact(args.controls))
        if args.write_receipt:
            args.receipt.write_bytes(codec.canonical(result))
        else:
            codec.equal(codec.load_artifact(args.receipt),result,"receipt mismatch")
    except (ValueError,KeyError,TypeError,OSError) as error:
        parser.exit(1,f"verification failed: {error}\n")
    print("Verified constrained KL selection, cover controls and native-read obstruction.")


if __name__ == "__main__":
    main()

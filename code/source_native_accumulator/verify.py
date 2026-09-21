"""Independent positional replay, interval decoder and whole-history accounting."""
import argparse
from fractions import Fraction
import hashlib
from itertools import product
from pathlib import Path
import sys

if __package__:
    from . import codec
else:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_native_accumulator import codec

PORTS = (0,9,1,5,7,4,8,11,14,40,23,38,45,39)
GRID = 1 << 256


def require(condition, message):
    if not condition:
        raise ValueError(message)


def decode(plus, minus, exponent, count, limit):
    center = Fraction(plus-minus,2*(GRID-3))*(1 << exponent)
    allowance = Fraction(4+5*count,8*(GRID-3))*(1 << exponent)
    left, right = center-allowance, center+allowance
    low = max(-limit, -((-left.numerator)//left.denominator))
    high = min(limit, right.numerator//right.denominator)
    return low if low == high else None


def support_edges():
    capture = codec.load(codec.ROOT/codec.SUPPORT)
    require(capture["schema"] == "oph.source_routing.w12_support.v1" and
            capture["level"] == 3 and capture["carriers"] == 1280,"captured support identity")
    edges = {frozenset((12*c+a,12*c+b)) for c in range(capture["carriers"])
             for a,b in capture["intra_carrier_seams"]}
    edges.update(frozenset((12*c+a,12*d+b)) for c,a,d,b in capture["glued_pairs"])
    require(len(set(PORTS)) == 14,"overlapping physical records")
    return edges


def replay(sessions, payload, edges, disturbed=False):
    """No producer import. Closed-form loop counts and exact rational means.

    The optional deterministic signed disturbances exercise the separately
    proved error model; the checked-in receipt uses ordinary rounded means.
    """
    require(type(sessions) is list and sessions and all(type(ss) is list for ss in sessions),
            "nonempty episode inventory")
    require(all(type(s) is int and s in (2,3) for ss in sessions for s in ss),"request alphabet")
    require(type(payload) is int and payload in range(-2,3),"payload domain")
    ideal = [Fraction(2*GRID) for _ in PORTS]
    ideal[4], ideal[5] = 3*GRID-3, GRID+3
    ideal[6], ideal[7] = 2*GRID+payload*(GRID-3), 2*GRID-payload*(GRID-3)
    actual = list(ideal)
    if disturbed:
        actual = [v+Fraction(1 if i%2 else -1,4) for i,v in enumerate(actual)]
    writers = [-p-1 for p in PORTS]
    exponents = [0]*7
    tape, polls = [], []
    odd = 0
    version = position = 0
    expected = 0
    max_scale = 0

    def mean(a,b,phase):
        nonlocal odd
        u,v = PORTS[a],PORTS[b]
        require(frozenset((u,v)) in edges,"unsupported native mean")
        total = actual[a]+actual[b]
        rounded = round(Fraction(total,2))
        if not disturbed:
            odd += int(total)%2
        tape.append([version,position,phase,u,v,int(actual[a]) if not disturbed else actual[a],
                     int(actual[b]) if not disturbed else actual[b],writers[a],writers[b],rounded])
        actual[a] = actual[b] = rounded
        writers[a] = writers[b] = len(tape)-1
        midpoint = Fraction(ideal[a]+ideal[b],2)
        ideal[a] = ideal[b] = midpoint
        if disturbed:
            actual[:] = [x+Fraction(1 if (i+len(tape))%2 else -1,8)
                         for i,x in enumerate(actual)]
        # In grid units: preparation <=1/4 and each mean <=1/2+1/8.
        allowance = Fraction(1,4)+len(tape)*Fraction(5,8)
        require(all(abs(x-y)<=allowance for x,y in zip(actual,ideal)),"whole-word error exceeded")

    def pair(a,b,phase):
        mean(2*a,2*b,phase)
        mean(2*a+1,2*b+1,phase)

    def observations(cells):
        nonlocal max_scale
        for cell,bound,value in cells:
            shift = exponents[cell]
            max_scale = max(max_scale,shift)
            want = Fraction(value*(GRID-3),1 << shift)
            require(Fraction(ideal[2*cell]-ideal[2*cell+1],2) == want,
                    "native linear semantics or retained input changed")
            radius = Fraction((4+5*len(tape))*(1 << shift),8*(GRID-3))
            require(2*radius<1,"finite precision cannot identify this integer")
            plus,minus = actual[2*cell:2*cell+2]
            if disturbed:
                plus += Fraction(1,4)
                minus -= Fraction(1,4)
            result = decode(plus,minus,shift,len(tape),bound)
            require(result == value,"incorrect or missing local publication")
            polls.append([version,position,cell,len(tape),shift,bound,
                          int(plus) if not disturbed else plus,
                          int(minus) if not disturbed else minus,result])

    for version, requests in enumerate(sessions):
        before = len(tape)
        position = -1
        for a,b in ((0,2),(1,3),(0,3),(1,2)):
            mean(a,b,"reset")
        pair(2,0,"seed")
        exponents[0] = exponents[2]+1
        exponents[2] += 1
        expected = bound = 1
        observations(((0,bound,expected),(2,1,1),(3,2,payload)))
        for position, source in enumerate(requests):
            level = max(exponents[0],exponents[source]+1)
            for cell,repetitions,phase in (
                (0,level-exponents[0],"align_accumulator"),
                (source,level-exponents[source]-1,"align_input")):
                for _ in range(repetitions):
                    pair(cell,1,phase)
                    mean(2,3,phase)
            pair(source,1,"add")
            pair(0,1,"add")
            mean(2,3,"add")
            exponents[0], exponents[source] = level+1,level
            expected += 1 if source == 2 else payload
            bound += 1 if source == 2 else 2
            observations(((0,bound,expected),(2,1,1),(3,2,payload)))
        require(len(tape)-before >= 6+5*len(requests),"native work discounted")
    n = sum(map(len,sessions))
    s = len(sessions)
    require(max_scale <= s+2*n,"scale growth bound")
    require(len(tape) <= 6*s+n*(8+3*(s+2*n)),"polynomial native work bound")
    version,position = len(sessions),-1
    receiver_exponent = exponents[0]+4
    for a,b in ((0,1),(1,4),(4,5),(5,6)):
        pair(a,b,"remote")
    exponents[0] += 1
    exponents[6] = receiver_exponent
    observations(((0,bound,expected),(2,1,1),(3,2,payload),(6,bound,expected)))
    require(len(polls) == 3*(s+n)+4,"checkpoint inventory")
    return {"sessions":sessions,"payload":payload,"tape":tape,"polls":polls,
            "scales":exponents,"odd_sums":odd,
            "final":[[p,int(actual[i]) if not disturbed else actual[i],writers[i]]
                     for i,p in sorted(enumerate(PORTS),key=lambda pair:pair[1])]}


def summary(cases, label):
    checksum = hashlib.sha256()
    h=s=n=w=r=o=maximum=cross=0
    for case in cases:
        checksum.update((codec.compact(case)+"\n").encode("ascii"))
        h += 1
        s += len(case["sessions"])
        n += sum(len(x) for x in case["sessions"])
        w += len(case["tape"])
        r += 2*len(case["polls"])
        o += case["odd_sums"]
        maximum = max(maximum,*case["scales"])
        cross += sum(divmod(e[3],12)[0] != divmod(e[4],12)[0] for e in case["tape"])
    return {"label":label,"histories":h,"episodes":s,"requests":n,"scalar_means":w,
            "scalar_samples":r,"odd_sums":o,"max_scale":maximum,"cross_carrier_means":cross,
            "scalar_mean_reads":w*2,"scalar_mean_writes":w*2,"preparation_writes":h*14,
            "all_histories_sha256":checksum.hexdigest()}


def check_header(packet):
    require(type(packet) is dict,"artifact object required")
    codec.equal(sorted(packet),sorted(("schema","pins","m1_derived","source_selected",
                                      "source","exhaustive","reuse","examples")),"artifact fields")
    codec.equal(packet["schema"],"oph-native-accumulator-v1","schema")
    codec.equal(packet["pins"],codec.pins(),"source/proof custody")
    codec.equal(packet["m1_derived"],False,"unproved M1 promotion")
    codec.equal(packet["source_selected"],False,"unproved source selection")
    codec.equal(packet["source"],{
        "support":codec.SUPPORT,"rails":[list(PORTS[i:i+2]) for i in range(0,14,2)],
        "grid_Q":GRID,"amplitude_units":GRID-3,"baseline_units":GRID*2,
        "payloads":[-2,-1,0,1,2],"request_alphabet":[2,3],"one_time_unit_seed":True,
        "accumulator_retired_at_start":True,"logical_input_values_retained":True,
        "raw_input_amplitudes_immutable":False,"error_clock":"global_scalar_mean_count",
        "schedule":"supplied_requests_with_data_independent_scale_alignment",
        "remote":"one_final_read_on_prepared_bus"},"interface and retained premises")
    require(type(packet["exhaustive"]) is list and len(packet["exhaustive"]) == 8,
            "complete short request inventory")
    for h,row in enumerate(packet["exhaustive"]):
        codec.equal(row["histories"],5*2**h,"complete signed payload/request enumeration")
    codec.equal(packet["reuse"]["histories"],5,"complete reuse inventory")


def verify(packet):
    check_header(packet)
    pins = codec.pins()
    edges = support_edges()
    examples = [replay(ps,p,edges) for ps,p in (([[]],0),([[3,2,3]],-2),([[3],[],[2,3]],2))]
    codec.equal(packet["examples"],examples,"native tape, custody, scales or local publication")
    rows = []
    for h in range(8):
        row = summary((replay([list(word)],p,edges) for word in product((2,3),repeat=h)
                       for p in (-2,-1,0,1,2)),f"length_{h}")
        codec.equal(packet["exhaustive"][h],row,"complete histories and paid work")
        rows.append(row)
    reuse = [[2,3,3,2], [3,3,3,3,3,3,3], [], [2]*9, [3,2]*8]*3
    row = summary((replay(reuse,p,edges) for p in (-2,-1,0,1,2)),"fifteen_episodes")
    codec.equal(packet["reuse"],row,"repeated episodes, retained seed and global error clock")
    rows.append(row)
    codec.equal(codec.pins(),pins,"sources changed during verification")
    return {"schema":"oph-native-accumulator-receipt-v1","pins":pins,
            "controls_sha256":codec.digest(packet),"m1_derived":False,"source_selected":False,
            "histories_replayed":sum(r["histories"] for r in rows),
            "scalar_means_replayed":sum(r["scalar_means"] for r in rows),
            "scalar_samples_replayed":sum(r["scalar_samples"] for r in rows),
            "odd_sums_replayed":sum(r["odd_sums"] for r in rows),
            "maximum_scale":max(r["max_scale"] for r in rows),
            "preparation_writes_per_history":14,"dynamic_scalar_writes":"pair_means_only",
            "error_model":{"initial":"1/(4Q)","per_mean":"5/(8Q)","sample":"1/(4Q)"},
            "instrument":"supplied_local_integer_interval_comparator",
            "cost_exclusions":["controller computation and metadata", "publication instrument",
                               "physical isolation", "inactive support", "evidence storage"],
            "scope":"mutable_accumulator_retained_logical_inputs_single_final_remote_read"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls",type=Path,default=codec.HERE/"controls.json")
    parser.add_argument("--receipt",type=Path,default=codec.HERE/"receipt.json")
    parser.add_argument("--write-receipt",action="store_true")
    args = parser.parse_args()
    try:
        receipt = verify(codec.load_artifact(args.controls))
        if args.write_receipt:
            args.receipt.write_bytes(codec.canonical(receipt))
        else:
            codec.equal(codec.load_artifact(args.receipt),receipt,"receipt mismatch")
    except (ValueError,KeyError,TypeError,OSError,IndexError,ZeroDivisionError) as error:
        print(f"verification failed: {error}",file=sys.stderr)
        return 1
    print(f"Verified {receipt['histories_replayed']} histories and {receipt['scalar_means_replayed']} native means.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Independent replay, influence/first-hit counting and source-contract checks."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
from pathlib import Path
import sys

if __package__:
    from . import codec
else:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_read_acceptance import codec

PORTS = (0,9,1,5,14,40,23,38,45,39)
GRID = 16384
SIGNAL = F(GRID-3, GRID)


def require(condition, label):
    if not condition:
        raise ValueError(label)


def local_decision(positive, negative, time):
    # Opposite-sign models overlap on [-E,E]; closed interval endpoints abstain.
    twice_difference = 2*(positive-negative)
    limit = 2+5*time
    if abs(twice_difference)>4*(GRID-3)+limit:
        return 0
    return (twice_difference > limit)-(twice_difference < -limit)


def support_check():
    packet = codec.load(codec.ROOT/codec.SUPPORT)
    adjacency = {i:set() for i in range(12*packet["carriers"])}
    def insert(a,b):
        adjacency[a].add(b)
        adjacency[b].add(a)
    for c in range(packet["carriers"]):
        for a,b in packet["intra_carrier_seams"]:
            insert(12*c+a,12*c+b)
    for c,a,d,b in packet["glued_pairs"]:
        insert(12*c+a,12*d+b)
    require(len(set(PORTS))==10,"injective rail placement")
    edges = [(PORTS[2*i+r],PORTS[2*i+2+r]) for i in range(4) for r in range(2)]
    require(all(b in adjacency[a] for a,b in edges),"unsupported scalar mean")
    return {"scalar_seams":[list(e) for e in edges],"active_scalar_registers":10,
            "source_carrier":0,"receiver_carrier":3,"support_ports":len(adjacency),
            "support_undirected_seams":sum(map(len,adjacency.values()))//2}


def replay(word, payload):
    """No producer imports: positional registers, exact means, nearest-integer rounding."""
    state = [2*GRID]*10
    state[0] += payload*(GRID-3)
    state[1] -= payload*(GRID-3)
    writers = [-p-1 for p in PORTS]
    initial_decision = local_decision(state[8],state[9],0)
    require(initial_decision==0,"initial local observation cannot identify the payload")
    tape, polls = [], [[0,state[8],state[9],initial_decision]]
    frontier = 0
    decision = 0
    # Dyadic source influence is propagated independently of either payload.
    weights, denominator = [1,0,0,0,0], 1
    for t, move in enumerate(word, 1):
        require(type(move) is int and 0<=move<4,"word alphabet")
        for rail in (0,1):
            a,b = 2*move+rail,2*move+2+rail
            mean = F(state[a]+state[b],2)
            result = round(mean)
            require(abs(F(result)-mean)<=F(1,2),"rounding bound")
            tape.append([PORTS[a],PORTS[b],state[a],state[b],writers[a],writers[b],result])
            state[a] = state[b] = result
            writers[a] = writers[b] = len(tape)-1
        old = weights
        weights = [2*x for x in old]
        weights[move] = weights[move+1] = old[move]+old[move+1]
        denominator *= 2
        if move == frontier:
            frontier += 1
        require(sum(weights)==denominator,"ideal total-load conservation")
        require([x>0 for x in weights]==[i<=frontier for i in range(5)],"causal frontier")
        require(sum(state)==20*GRID,"paired fixed-point total-load conservation")
        error = F(2+5*t,4*GRID)
        for cell in range(5):
            ideal = payload*SIGNAL*F(weights[cell],denominator)
            for rail in (0,1):
                expected = 2+(-1)**rail*ideal
                require(abs(F(state[2*cell+rail],GRID)-expected)<=F(t,GRID),
                        "native rounding accumulation")
        ideal_contrast = payload*SIGNAL*F(weights[4],denominator)
        observed = F(state[8]-state[9],2*GRID)
        require(abs(observed-ideal_contrast)<=error,"receiver error enclosure")
        require(2*error<SIGNAL/F(2**t),"finite-horizon separation")
        decision = local_decision(state[8],state[9],t)
        expected_result = payload if frontier==4 else 0
        require(decision==expected_result,"local publication differs from first influence")
        # Both extreme local readout errors are allowed; the receiver still sees
        # no payload or route. Arbitrary per-operation errors are covered in Lean.
        for shift in (-F(1,4),F(0),F(1,4)):
            require(local_decision(F(state[8])+shift,F(state[9])-shift,t)==expected_result,
                    "readout-error endpoint changed acceptance")
        polls.append([t,state[8],state[9],decision])
        if decision:
            break
    return {"payload":payload,"word":list(word),"tape":tape,"polls":polls,
            "result":decision,"steps":len(polls)-1,
            "final":sorted([[p,x,w] for p,x,w in zip(PORTS,state,writers)])}


def arrival_law(edges,horizon):
    """Counts all words, including the absorbing state's unexecuted suffixes."""
    counts = [1]+[0]*edges
    first = [0]*(horizon+1)
    expected = F(0)
    for t in range(1,horizon+1):
        expected += F(sum(counts[:-1]),edges**(t-1))
        first[t] = counts[edges-1]*edges**(horizon-t)
        next_counts = [0]*(edges+1)
        for progress,number in enumerate(counts):
            if progress==edges:
                next_counts[progress] += edges*number
            else:
                next_counts[progress] += (edges-1)*number
                next_counts[progress+1] += number
        counts = next_counts
        require(sum(counts)==edges**t,"word mass dropped during stopping")
    abort = sum(counts[:-1])
    require(sum(first)+abort==edges**horizon,"first-hit partition")
    require(expected==F(sum(t*n for t,n in enumerate(first))+horizon*abort,edges**horizon),
            "unconditional tail/work identity")
    return first,abort,expected


def verify_row(row, horizon):
    first,aborts,expected = arrival_law(4,horizon)
    histogram = [0]*(horizon+1)
    lost = means = crossings = samples = wrong = 0
    histories = hashlib.sha256()
    for word in product(range(4),repeat=horizon):
        cases = [replay(word,payload) for payload in (-1,1)]
        require(cases[0]["steps"]==cases[1]["steps"],"payload-dependent stopping channel")
        require(cases[0]["result"]==-cases[1]["result"],"intervention asymmetry")
        if cases[0]["result"]:
            histogram[cases[0]["steps"]] += 1
        else:
            lost += 1
        for case in cases:
            histories.update((codec.compact(case)+"\n").encode("ascii"))
            wrong += case["result"] not in (0,case["payload"])
            means += 2*case["steps"]
            crossings += sum(e[0]//12 != e[1]//12 for e in case["tape"])
            samples += 2*(case["steps"]+1)
    require(histogram==first and lost==aborts,"execution/independent arrival law")
    require(means==4*expected*4**horizon,"failed-attempt work omitted")
    calculated = {
        "horizon":horizon,"words_per_payload":4**horizon,
        "first_publication_counts":histogram,"aborts_per_payload":lost,
        "wrong_publications":wrong,"paired_steps_per_payload":int(expected*4**horizon),
        "scalar_means_both_payloads":means,"cross_carrier_means_both_payloads":crossings,
        "scalar_reads_both_payloads":2*means,"scalar_writes_both_payloads":2*means,
        "receiver_samples_both_payloads":samples,"unconditional_expected_means":str(2*expected),
        "all_stopped_histories_sha256":histories.hexdigest(),
    }
    codec.equal(row,calculated,"exhaustive stopped histories, failures or paid work")
    return {"horizon":horizon,"coverage":str(1-F(aborts,4**horizon)),
            "abort_probability":str(F(aborts,4**horizon)),
            "expected_scalar_means":str(2*expected),"wrong_publications":wrong}


def check_count_claims(row,horizon):
    """Reject false denominators, omitted failures and discounted work before replay."""
    first,aborts,expected = arrival_law(4,horizon)
    values = {"horizon":horizon,"words_per_payload":4**horizon,
              "first_publication_counts":first,"aborts_per_payload":aborts,
              "wrong_publications":0,"paired_steps_per_payload":int(expected*4**horizon),
              "scalar_means_both_payloads":int(4*expected*4**horizon),
              "unconditional_expected_means":str(2*expected)}
    for key,value in values.items():
        codec.equal(row[key],value,"unconditional count/cost certificate")


def path_law(edges,horizon):
    _,abort,expected = arrival_law(edges,horizon)
    power = 1
    while power<=10*horizon+4:
        power *= 2
    grid = 4*2**horizon*power
    radius = F(2+5*horizon,4*grid)
    require(2*radius<F(1,2**horizon),"insufficient horizon-dependent precision")
    return {"edges":edges,"horizon":horizon,"abort_probability":str(F(abort,edges**horizon)),
            "expected_paired_steps":str(expected),"expected_scalar_means":str(2*expected),
            "sufficient_grid_Q":str(grid),"total_error_radius":str(radius),
            "minimum_nonzero_ideal_signal":str(F(1,2**horizon)),
            "ideal_unbounded_expected_scalar_means":2*edges*edges}


def verify(packet):
    require(type(packet) is dict,"artifact object required")
    codec.equal(sorted(packet),sorted(("schema","pins","classification","m1_derived",
        "full_axiom_source_contract","source","exhaustive","examples","path_laws")),"artifact fields")
    codec.equal(packet["schema"],"oph-local-read-acceptance-v1","schema")
    source_pins = codec.pins()
    codec.equal(packet["pins"],source_pins,"source/proof custody")
    codec.equal(packet["classification"],"exact_named_realization_with_conditional_schedule","classification")
    codec.equal(packet["m1_derived"],False,"unproved M1 promotion")
    codec.equal(packet["full_axiom_source_contract"],False,"unproved source contract")
    source = {"support":codec.SUPPORT,"rails":[list(PORTS[i:i+2]) for i in range(0,10,2)],
              "grid_Q":GRID,"amplitude_units":GRID-3,"baseline_units":2*GRID,
              "payloads":[-1,1],"macro_alphabet":[[i,i+1] for i in range(4)],
              "checkpoints":"after_each_complete_two_mean_copy",
              "reference":"uniform_full_word_counting","cover":"identity_on_complete_macro_histories",
              "feasible":"all_four_letter_words_including_aborts","source_selection":"not_established"}
    codec.equal(packet["source"],source,"source/interface/schedule specification")
    support = support_check()
    require(type(packet["exhaustive"]) is list and len(packet["exhaustive"])==9,"complete horizon inventory")
    for h,row in enumerate(packet["exhaustive"]):
        check_count_claims(row,h)
    laws = [path_law(d,h) for d in (1,2,4,8) for h in (0,1,8,32,64)]
    codec.equal(packet["path_laws"],laws,"all-horizon arrival/work/precision laws")
    words = [(0,1,2,3),(3,2,1,0,1,2,3),(0,)*8,(0,1,2,2,3,1,0,3)]
    examples = [replay(word,payload) for word in words for payload in (-1,1)]
    codec.equal(packet["examples"],examples,"native tape, custody or local publication")
    rows = [verify_row(row,h) for h,row in enumerate(packet["exhaustive"])]
    codec.equal(codec.pins(),source_pins,"sources changed during verification")
    return {"schema":"oph-local-read-acceptance-receipt-v1","pins":source_pins,
            "controls_sha256":codec.digest(packet),"support":support,"rows":rows,
            "m1_derived":False,"source_contract":"conditional_not_source_selected",
            "scalar_means_replayed":sum(r["scalar_means_both_payloads"] for r in packet["exhaustive"]),
            "histories_replayed":2*sum(4**h for h in range(9)),
            "preparation_writes_per_history":10,"comparator_calls_per_history":"paired_steps + 1",
            "publication_instrument":"supplied logical controller, not a scalar-mean construction",
            "exclusions":["physical isolation", "inactive support", "evidence storage", "immutable controller code"],
            "path_laws":laws}


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
    except (ValueError,KeyError,TypeError,OSError,IndexError,ZeroDivisionError) as error:
        print(f"verification failed: {error}",file=sys.stderr)
        return 1
    print(f"Verified {result['histories_replayed']} histories and {result['scalar_means_replayed']} native means; M1 remains conditional.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

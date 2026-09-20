"""Producer: exact scalar histories and a captured cut census."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib

if __package__:
    from . import codec,constraints
else:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_read_selection import codec,constraints


def edges():
    support = codec.load(codec.ROOT/codec.SUPPORT)
    result = {(12*c+u, 12*c+v) for c in range(support["carriers"])
              for u, v in support["intra_carrier_seams"]}
    result.update((12*c+u, 12*d+v) for c,u,d,v in support["glued_pairs"])
    return sorted({tuple(sorted(edge)) for edge in result})


def execute(initial, word, receiver):
    state = dict(initial)
    writers = {i: -i-1 for i in state}
    tape = []
    for event, (u,v) in enumerate(word):
        value = (state[u]+state[v])/2
        tape.append([u,v,str(state[u]),str(state[v]),writers[u],writers[v],str(value)])
        state[u] = state[v] = value
        writers[u] = writers[v] = event
    return [str(state[r]) for r in receiver], tape


def toy(horizon):
    alphabet = [(0,1),(1,2),(2,3)]
    failures = [0,0]
    avoiding = operations = 0
    histories = hashlib.sha256()
    for word_ids in product(range(3), repeat=horizon):
        word = [alphabet[i] for i in word_ids]
        avoiding += 0 not in word_ids
        for index, payload in enumerate((-1,1)):
            final,tape = execute({0:F(payload),1:F(0),2:F(0),3:F(0)},word,[3])
            value = F(final[0])
            decoded = (value > 0)-(value < 0)
            failures[index] += decoded != payload
            operations += len(tape)
            histories.update((codec.compact([list(word_ids),payload,tape,final,decoded])+"\n").encode())
    return {"horizon":horizon,"words":3**horizon,"avoid_words":avoiding,
            "failure_counts":failures,"scalar_means":operations,
            "scalar_reads":2*operations,"scalar_writes":2*operations,
            "histories_sha256":histories.hexdigest()}


def probability_row(total, cut, horizon):
    mass = F(total-cut,total)**horizon
    scaled = mass*2**32
    lower = scaled.numerator//scaled.denominator
    upper = lower + (scaled != lower)
    return {"horizon":horizon,"avoid_mass":{"base":str(F(total-cut,total)),"exponent":horizon},
            "denominator":2**32,"lower_units":lower,"upper_units":upper,
            "balanced_error_lower":str(F(lower,2**33))}


def captured_controls():
    route = [(0,1),(9,5),(1,14),(5,40),(14,23),(40,38),(23,45),(38,39)]
    avoiding = [(0,1),(14,23),(23,45),(38,39)]*4
    ports = sorted({p for word in (route,avoiding) for e in word for p in e} | {9})
    controls = []
    for name,word in (("routed",route),("cut_avoiding",avoiding)):
        cases = []
        for payload in (-1,1):
            initial = {p:F(2) for p in ports}
            initial[0] += payload
            initial[9] -= payload
            local,tape = execute(initial,word,[45,39])
            contrast = (F(local[0])-F(local[1]))/2
            decoded = (contrast > 0)-(contrast < 0)
            cases.append({"payload":payload,"initial":[[p,str(initial[p])] for p in ports],
                          "local":local,"decoded":decoded,"tape":tape,"tape_sha256":codec.digest(tape)})
        controls.append({"name":name,"word":[list(e) for e in word],"cases":cases,
                         "means_per_case":len(word),"reads_per_case":2*len(word),
                         "writes_per_case":2*len(word)})
    return controls


def build():
    alphabet = edges()
    cut = [list(e) for e in alphabet if (e[0]<12)!=(e[1]<12)]
    return {"schema":"oph-read-selection-v1","pins":codec.pins(),
            "scope":"conditional_classical_KL_constraints_cover_and_read_obstruction",
            "alphabet":{"identity":"distinct_undirected_scalar_seams","size":len(alphabet),
                        "sha256":codec.digest(alphabet)},
            "source_ports":list(range(12)),"receiver_ports":[45,39],"cut":cut,
            "probabilities":[probability_row(len(alphabet),len(cut),h) for h in (0,1,2,16,256,2279)],
            "toy":{"alphabet":[[0,1],[1,2],[2,3]],"source":0,"receiver":3,
                   "payloads":[-1,1],"decoder":"sign_with_zero_failure",
                   "rows":[toy(h) for h in range(7)]},
            "captured_controls":captured_controls(),"projection_controls":constraints.build()}


if __name__ == "__main__":
    (codec.HERE/"controls.json").write_bytes(codec.canonical(build()))

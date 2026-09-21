"""Producer: stopped native histories without conditioning on successful reads."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import hashlib

if __package__:
    from . import codec
else:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_read_acceptance import codec

RAILS = [(0,9), (1,5), (14,40), (23,38), (45,39)]
Q = 2**14
AMPLITUDE = Q-3
HORIZON = 8


def budget(t):
    return F(1,4*Q)+2*t*(F(1,2*Q)+F(1,8*Q))+F(1,4*Q)


def publish(positive_sample, negative_sample, checkpoint):
    """Only the receiver's two samples and public error budget enter here."""
    contrast = (positive_sample-negative_sample)/2
    radius = budget(checkpoint)
    if abs(contrast)>F(AMPLITUDE,Q)+radius:
        return 0
    return 1 if contrast > radius else -1 if contrast < -radius else 0


def execute(word, payload):
    ports = sorted(p for pair in RAILS for p in pair)
    state = {p:2*Q for p in ports}
    state[0] += payload*AMPLITUDE
    state[9] -= payload*AMPLITUDE
    writers = {p:-p-1 for p in ports}
    tape = []
    result = publish(F(state[45],Q),F(state[39],Q),0)
    polls = [[0, state[45], state[39], result]]
    for t, edge in enumerate(word, 1):
        for rail in (0,1):
            u,v = RAILS[edge][rail], RAILS[edge+1][rail]
            total = state[u]+state[v]
            value = total//2 + (total % 4 == 3)
            tape.append([u,v,state[u],state[v],writers[u],writers[v],value])
            state[u] = state[v] = value
            writers[u] = writers[v] = len(tape)-1
        result = publish(F(state[45],Q), F(state[39],Q), t)
        polls.append([t,state[45],state[39],result])
        if result:
            break
    return {"payload":payload, "word":list(word), "tape":tape, "polls":polls,
            "result":result, "steps":len(tape)//2,
            "final":[[p,state[p],writers[p]] for p in ports]}


def row(horizon):
    counts = [0]*(horizon+1)
    aborts = means = crossings = samples = wrong = 0
    streams = hashlib.sha256()
    for word in product(range(4), repeat=horizon):
        pair = [execute(word, payload) for payload in (-1,1)]
        first = pair[0]
        if first["result"]:
            counts[first["steps"]] += 1
        else:
            aborts += 1
        for case in pair:
            wrong += case["result"] not in (0,case["payload"])
            streams.update((codec.compact(case)+"\n").encode("ascii"))
            means += len(case["tape"])
            crossings += sum(u//12 != v//12 for u,v,*_ in case["tape"])
            samples += 2*len(case["polls"])
    paid = sum(i*n for i,n in enumerate(counts))+horizon*aborts
    return {"horizon":horizon, "words_per_payload":4**horizon,
            "first_publication_counts":counts, "aborts_per_payload":aborts,
            "wrong_publications":wrong, "paired_steps_per_payload":paid,
            "scalar_means_both_payloads":means, "cross_carrier_means_both_payloads":crossings,
            "scalar_reads_both_payloads":2*means, "scalar_writes_both_payloads":2*means,
            "receiver_samples_both_payloads":samples,
            "unconditional_expected_means":str(F(2*paid,4**horizon)),
            "all_stopped_histories_sha256":streams.hexdigest()}


def counting(d, h):
    failure = sum(comb(h,j)*(d-1)**(h-j) for j in range(min(d,h+1)))
    tail = F(failure,d**h)
    expected = sum(F(sum(comb(t,j)*(d-1)**(t-j)
                         for j in range(min(d,t+1))),d**t) for t in range(h))
    # A sufficient finite-horizon grid; the hardware-error allowance must scale too.
    q = 2**(h+(10*h+4).bit_length()+2)
    radius = F(1,2*q)+2*h*(F(1,2*q)+F(1,8*q))
    return {"edges":d,"horizon":h,"abort_probability":str(tail),
            "expected_paired_steps":str(expected),"expected_scalar_means":str(2*expected),
            "sufficient_grid_Q":str(q),"total_error_radius":str(radius),
            "minimum_nonzero_ideal_signal":str(F(1,2**h)),
            "ideal_unbounded_expected_scalar_means":2*d*d}


def build():
    examples = [(0,1,2,3), (3,2,1,0,1,2,3), (0,)*8, (0,1,2,2,3,1,0,3)]
    return {
        "schema":"oph-local-read-acceptance-v1", "pins":codec.pins(),
        "classification":"exact_named_realization_with_conditional_schedule",
        "m1_derived":False, "full_axiom_source_contract":False,
        "source":{"support":codec.SUPPORT,"rails":[list(x) for x in RAILS],
                  "grid_Q":Q,"amplitude_units":AMPLITUDE,"baseline_units":2*Q,
                  "payloads":[-1,1],"macro_alphabet":[[i,i+1] for i in range(4)],
                  "checkpoints":"after_each_complete_two_mean_copy",
                  "reference":"uniform_full_word_counting",
                  "cover":"identity_on_complete_macro_histories",
                  "feasible":"all_four_letter_words_including_aborts",
                  "source_selection":"not_established"},
        "exhaustive":[row(h) for h in range(HORIZON+1)],
        "examples":[execute(word,payload) for word in examples for payload in (-1,1)],
        "path_laws":[counting(d,h) for d in (1,2,4,8) for h in (0,1,8,32,64)],
    }


if __name__ == "__main__":
    packet = build()
    codec.equal(packet["pins"],codec.pins(),"sources changed during generation")
    (codec.HERE/"controls.json").write_bytes(codec.canonical(packet))

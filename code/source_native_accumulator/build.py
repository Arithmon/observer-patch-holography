"""Producer: prepare once, then write records only through scalar pair means."""
from fractions import Fraction as F
import hashlib
from itertools import product
from math import ceil, floor
from pathlib import Path

if __package__:
    from . import codec
else:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_native_accumulator import codec

RAILS = [(0,9),(1,5),(7,4),(8,11),(14,40),(23,38),(45,39)]
Q = 2**256
G = Q-3
STRESS = [[2,3,3,2], [3]*7, [], [2]*9, [3,2]*8]*3


def publish(p, m, scale, means, bound):
    """Local samples plus public metadata; no oracle value or source sample."""
    observation = F((p-m)*2**scale, 2*G)
    radius = F((4+5*means)*2**scale, 8*G)
    lower = max(-bound, ceil(observation-radius))
    upper = min(bound, floor(observation+radius))
    return lower if lower == upper else None


def execute(sessions, payload):
    if not sessions or any(type(s) is not int or s not in (2,3) for ss in sessions for s in ss):
        raise ValueError("nonempty list of episodes over inputs 2 and 3 required")
    if type(payload) is not int or payload not in range(-2,3):
        raise ValueError("integer payload outside declared preparation")
    ports = sorted(p for pair in RAILS for p in pair)
    state = {p:2*Q for p in ports}
    for cell, value in ((2,1),(3,payload)):
        state[RAILS[cell][0]] += value*G
        state[RAILS[cell][1]] -= value*G
    writers = {p:-p-1 for p in ports}
    scale = [0]*7
    tape, polls = [], []
    odd = 0
    session = step = 0

    def mean(u, v, phase):
        nonlocal odd
        total = state[u]+state[v]
        value = total//2 + (total % 4 == 3)
        odd += total % 2
        tape.append([session,step,phase,u,v,state[u],state[v],writers[u],writers[v],value])
        state[u] = state[v] = value
        writers[u] = writers[v] = len(tape)-1

    def copy(s, t, phase):
        for rail in (0,1):
            mean(RAILS[s][rail], RAILS[t][rail], phase)

    def clear(phase):
        mean(*RAILS[1], phase)

    def shrink(s, phase):
        copy(s, 1, phase)
        clear(phase)
        scale[s] += 1

    def poll(cell, bound):
        p, m = (state[i] for i in RAILS[cell])
        polls.append([session,step,cell,len(tape),scale[cell],bound,p,m,
                      publish(p,m,scale[cell],len(tape),bound)])

    for session, requests in enumerate(sessions):
        step = -1
        # Retire the previous accumulator; retain the original input records.
        copy(0,1,"reset")
        mean(RAILS[0][0],RAILS[1][1],"reset")
        mean(RAILS[0][1],RAILS[1][0],"reset")
        copy(2,0,"seed")
        scale[2] += 1
        scale[0] = scale[2]
        bound = 1
        for cell, limit in ((0,bound),(2,1),(3,2)):
            poll(cell,limit)
        for step, source in enumerate(requests):
            target = max(scale[0], scale[source]+1)
            while scale[0]<target:
                shrink(0,"align_accumulator")
            while scale[source]+1<target:
                shrink(source,"align_input")
            copy(source,1,"add")
            scale[source] += 1
            copy(0,1,"add")
            scale[0] = target+1
            clear("add")
            bound += 1 if source == 2 else 2
            for cell, limit in ((0,bound),(2,1),(3,2)):
                poll(cell,limit)
    session, step = len(sessions), -1
    receiver_scale = scale[0]+4
    for s,t in ((0,1),(1,4),(4,5),(5,6)):
        copy(s,t,"remote")
    scale[0] += 1
    scale[6] = receiver_scale
    for cell,limit in ((0,bound),(2,1),(3,2),(6,bound)):
        poll(cell,limit)
    return {"sessions":[list(ss) for ss in sessions],"payload":payload,
            "tape":tape,"polls":polls,"scales":scale,"odd_sums":odd,
            "final":[[p,state[p],writers[p]] for p in ports]}


def summarize(cases, label):
    stream = hashlib.sha256()
    totals = {"label":label,"histories":0,"episodes":0,"requests":0,
              "scalar_means":0,"scalar_samples":0,"odd_sums":0,"max_scale":0,
              "cross_carrier_means":0}
    for case in cases:
        stream.update((codec.compact(case)+"\n").encode("ascii"))
        totals["histories"] += 1
        totals["episodes"] += len(case["sessions"])
        totals["requests"] += sum(map(len,case["sessions"]))
        totals["scalar_means"] += len(case["tape"])
        totals["scalar_samples"] += 2*len(case["polls"])
        totals["odd_sums"] += case["odd_sums"]
        totals["max_scale"] = max(totals["max_scale"],max(case["scales"]))
        totals["cross_carrier_means"] += sum(e[3]//12 != e[4]//12 for e in case["tape"])
    totals["scalar_mean_reads"] = 2*totals["scalar_means"]
    totals["scalar_mean_writes"] = 2*totals["scalar_means"]
    totals["preparation_writes"] = 14*totals["histories"]
    totals["all_histories_sha256"] = stream.hexdigest()
    return totals


def build():
    return {
        "schema":"oph-native-accumulator-v1", "pins":codec.pins(),
        "m1_derived":False, "source_selected":False,
        "source":{"support":codec.SUPPORT,"rails":[list(p) for p in RAILS],
                  "grid_Q":Q,"amplitude_units":G,"baseline_units":2*Q,
                  "payloads":list(range(-2,3)),"request_alphabet":[2,3],
                  "one_time_unit_seed":True,"accumulator_retired_at_start":True,
                  "logical_input_values_retained":True,"raw_input_amplitudes_immutable":False,
                  "error_clock":"global_scalar_mean_count",
                  "schedule":"supplied_requests_with_data_independent_scale_alignment",
                  "remote":"one_final_read_on_prepared_bus"},
        "exhaustive":[summarize((execute([word],p) for word in product((2,3),repeat=h)
                                 for p in range(-2,3)), f"length_{h}") for h in range(8)],
        "reuse":summarize((execute(STRESS,p) for p in range(-2,3)),"fifteen_episodes"),
        "examples":[execute(sessions,p) for sessions,p in
                    (([[]],0), ([[3,2,3]],-2), ([[3],[],[2,3]],2))],
    }


if __name__ == "__main__":
    packet = build()
    codec.equal(packet["pins"],codec.pins(),"sources changed during generation")
    (codec.HERE/"controls.json").write_bytes(codec.canonical(packet))

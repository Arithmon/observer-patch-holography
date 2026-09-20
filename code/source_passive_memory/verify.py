"""Independently decode and replay retained controls, never importing build.

Both experiment identity and every executed read/write are checked. A valid
trace for a different initial state, port map or word is not this experiment.
"""
from collections import deque
from fractions import Fraction
if __package__:
    from . import codec
else:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_passive_memory import codec


def support_graph():
    support = codec.load(codec.ROOT/codec.SUPPORT)
    codec.equal(support["schema"], "oph.source_routing.w12_support.v1", "support schema")
    n = 12*support["carriers"]
    graph = [set() for _ in range(n)]

    def connect(u, v):
        if not (0 <= u < n and 0 <= v < n and u != v):
            raise ValueError("support edge")
        graph[u].add(v)
        graph[v].add(u)

    for c in range(support["carriers"]):
        for u,v in support["intra_carrier_seams"]:
            if not (0 <= u < 12 and 0 <= v < 12):
                raise ValueError("intra port")
            connect(12*c+u,12*c+v)
    for a,u,b,v in support["glued_pairs"]:
        if not (0 <= a < support["carriers"] and 0 <= b < support["carriers"]
                and 0 <= u < 12 and 0 <= v < 12):
            raise ValueError("glued port")
        connect(12*a+u,12*b+v)
    seen, sizes = set(), []
    for start in range(n):
        if start in seen:
            continue
        queue, count = deque([start]), 0
        seen.add(start)
        while queue:
            u = queue.popleft()
            count += 1
            for v in graph[u]:
                if v not in seen:
                    seen.add(v)
                    queue.append(v)
        sizes.append(count)
    census = {"vertices": n, "edges": sum(map(len,graph))//2,
              "component_sizes": sorted(sizes),
              "mapped_control_ports": [0,1,5,7,10,11]}
    for v in census["mapped_control_ports"][1:]:
        if v not in graph[0]:
            raise ValueError("captured control seam absent")
    return census


def recipes():
    # This independent specification is deliberately not supplied by build.py.
    for m in [0,1,4,12,32]:
        yield (f"cold_star_{m}", [Fraction(3,2)]+[Fraction(0)]*m,
               [(0,v) for v in range(1,m+1)], 0, None)
    yield ("reused_zero", [Fraction(3,2),Fraction(0)], [(0,1)]*8, 0, None)
    yield ("prepared_balanced_copy_and_cleanup", list(map(Fraction,[5,3,4,4,6,2])),
           [(2,4),(3,5),(2,3)], 4, None)
    yield ("captured_five_neighbours", [Fraction(3,2)]+[Fraction(0)]*5,
           [(0,v) for v in range(1,6)], 0, [0,1,5,7,10,11])


def vector(items):
    if not isinstance(items,list):
        raise ValueError("vector must be a list")
    return [codec.rational(v) for v in items]


def energy(state, baseline):
    return sum(((x-baseline)**2 for x in state), Fraction(0))


def replay(case, recipe):
    name, initial, word, baseline, ports = recipe
    expected_keys = {"name","kind","baseline","global_ports","allowed_edges","initial",
                     "events","final","final_writers","metrics"}
    codec.equal(sorted(case),sorted(expected_keys),"case fields")
    codec.equal(case["name"],name,"experiment name")
    codec.equal(case["kind"],"finite_pair_mean_execution","execution kind")
    codec.equal(case["baseline"],str(baseline),"baseline")
    codec.equal(case["global_ports"],ports,"port identity")
    codec.equal(case["initial"],list(map(str,initial)),"initial state")
    codec.equal(case["allowed_edges"],[list(e) for e in sorted(set(word))],"operation menu")
    codec.equal(len(case["events"]),len(word),"event count")
    state = vector(case["initial"])
    n, m = len(state), len(word)
    writers = [f"init:{i}" for i in range(n)]
    touches = [0]*n
    losses, prefix_states = [], [state.copy()]
    for k, ((u,v), row) in enumerate(zip(word,case["events"])):
        old = state.copy()
        result = (old[u]+old[v])/2
        defect = (old[u]-old[v])**2/2
        expected = {"index": k, "op": "mean", "edge": [u,v],
                    "inputs": [str(old[u]),str(old[v])],
                    "writers": [writers[u],writers[v]],
                    "output": str(result), "loss": str(defect)}
        codec.equal(row,expected,f"event {k} read/write")
        # Construct the whole post-state independently of the builder's updates.
        state = [result if i in (u,v) else old[i] for i in range(n)]
        writers = [f"event:{k}" if i in (u,v) else writers[i] for i in range(n)]
        touches = [touches[i]+int(i in (u,v)) for i in range(n)]
        if sum(state) != sum(old) or energy(state,baseline)+defect != energy(old,baseline):
            raise ValueError("scalar conservation or quadratic ledger")
        if min(state) < 0 or any(initial[i] > 2**touches[i]*state[i] for i in range(n)):
            raise ValueError("nonnegative native-zero lower bound")
        losses.append(defect)
        prefix_states.append(state.copy())
    metrics = {"registers": n, "initializations": n, "means": m,
               "reads": 2*m, "writes_including_initialization": n+2*m,
               "touches": touches, "total_load": str(sum(initial)),
               "initial_quadratic": str(energy(initial,baseline)),
               "final_quadratic": str(energy(state,baseline)), "loss": str(sum(losses))}
    codec.equal(case["final"],list(map(str,state)),"final state")
    codec.equal(case["final_writers"],writers,"final writer provenance")
    codec.equal(case["metrics"],metrics,"resource metrics")
    if name.startswith("cold_star_") or name == "captured_five_neighbours":
        exact = [initial[0]/2**m]+[initial[0]/2**j for j in range(1,m+1)]
        if state != exact:
            raise ValueError("cold construction did not attain the bound")
    if name == "reused_zero" and (state != [Fraction(3,4)]*2 or any(losses[1:])):
        raise ValueError("reused zero incorrectly treated as fresh")
    if name == "prepared_balanced_copy_and_cleanup":
        if prefix_states[2] != list(map(Fraction,[5,3,5,3,5,3])):
            raise ValueError("prepared copy target")
        if state != list(map(Fraction,[5,3,4,4,5,3])):
            raise ValueError("cleanup hid ancillary change")
        # The preloaded ancillas already encode the payload. This is an exact
        # budget witness, not a derivation of unknown-input copying.
        before, copied = prefix_states[0], prefix_states[2]
        if energy(copied[4:],4)+2+sum(losses[:2]) != energy(before[4:],4):
            raise ValueError("copy ancillary budget")
        hypothetical = list(map(Fraction,[5,3,5,3,6,2]))
        if sum(hypothetical) != sum(initial) or energy(hypothetical,4) != energy(initial,4)+2:
            raise ValueError("catalytic target control")
    return {"name":name,**metrics}


def verify_fibre(cert):
    codec.equal(sorted(cert), sorted(["kind","initial","concentrated","witnesses",
                "separate_prepared_inputs","means_across_separate_inputs"]), "fibre fields")
    codec.equal(cert["kind"],"equal_image_witnesses_not_a_trajectory","fibre scope")
    codec.equal(cert["initial"],["1","2","3","4"],"fibre start")
    codec.equal(cert["concentrated"],["10","0","0","0"],"fibre end")
    codec.equal(cert["separate_prepared_inputs"],12,"fibre preparations")
    codec.equal(cert["means_across_separate_inputs"],12,"fibre means")
    moves = [(3,2,4),(2,1,4),(1,0,4),(2,1,3),(1,0,3),(1,0,2)]
    codec.equal(len(cert["witnesses"]),len(moves),"fibre witness count")
    current = vector(cert["initial"])
    for (u,v,a), witness in zip(moves,cert["witnesses"]):
        codec.equal(sorted(witness),sorted(["edge","amount","left","right","common_image"]),
                    "fibre witness fields")
        codec.equal(witness["edge"],[u,v],"fibre seam")
        codec.equal(witness["amount"],str(a),"fibre addend")
        left, right, image = (vector(witness[k]) for k in ("left","right","common_image"))
        if left != current or len(right) != 4 or len(image) != 4:
            raise ValueError("fibre chain")
        background = [left[i]-(a if i == u else 0) for i in range(4)]
        expected = [background[i]+(a if i == v else 0) for i in range(4)]
        if min(background) < 0 or right != expected:
            raise ValueError("fibre requires a nonnegative shared background")
        for state in (left,right):
            post = [(state[u]+state[v])/2 if i in (u,v) else state[i] for i in range(4)]
            if image != post or min(state) < 0 or sum(state) != 10:
                raise ValueError("fibre images must agree under actual means")
        current = right
    if current != vector(cert["concentrated"]):
        raise ValueError("fibre did not concentrate")
    return {"equal_image_witnesses":len(moves),"separate_prepared_inputs":12}


def verify(data=None):
    if data is None:
        data = codec.load(codec.HERE/"controls.json")
    codec.equal(sorted(data),sorted(["schema","source_sha256","support","executions","fibre"]),
                "control fields")
    codec.equal(data["schema"],"oph.source_passive_memory.controls.v1","schema")
    codec.equal(data["source_sha256"],codec.source_pins(),"source pins")
    census = support_graph()
    codec.equal(data["support"],census,"captured support census")
    expected = list(recipes())
    codec.equal(len(data["executions"]),len(expected),"execution count")
    results = [replay(case,recipe) for case,recipe in zip(data["executions"],expected)]
    fibre = verify_fibre(data["fibre"])
    return {"schema":"oph.source_passive_memory.verified.v1",
            "support":census,"executions":results,"fibre":fibre,
            "executed_means":sum(r["means"] for r in results),
            "source_sha256":codec.source_pins()}


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true",
                        help="write the receipt after independently verifying retained controls")
    args = parser.parse_args()
    result = verify()
    if args.write_receipt:
        (codec.HERE/"receipt.json").write_bytes(codec.canonical(result))
    else:
        codec.equal(codec.load(codec.HERE/"receipt.json"), result,"retained verification receipt")
    print(f"Verified {len(result['executions'])} retained executions, "
          f"{result['executed_means']} means, and {result['fibre']['equal_image_witnesses']} fibre witnesses.")


if __name__ == "__main__":
    main()

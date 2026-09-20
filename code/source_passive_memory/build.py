"""Generate explicit finite controls. Does not certify their correctness.

Run from the repository root: python code/source_passive_memory/build.py
"""
from fractions import Fraction as F
if __package__:
    from . import codec
else:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_passive_memory import codec


def graph_census():
    data = codec.load(codec.ROOT/codec.SUPPORT)
    n = data["carriers"]*12
    edges = {(12*c+u, 12*c+v) for c in range(data["carriers"])
             for u, v in data["intra_carrier_seams"]}
    edges.update((12*a+u, 12*b+v) for a, u, b, v in data["glued_pairs"])
    edges = {tuple(sorted(e)) for e in edges}
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for u, v in sorted(edges):
        parent[find(u)] = find(v)
    counts = {}
    for i in range(n):
        root = find(i)
        counts[root] = counts.get(root, 0)+1
    return {"vertices": n, "edges": len(edges),
            "component_sizes": sorted(counts.values()),
            "mapped_control_ports": [0, 1, 5, 7, 10, 11]}


def trace(name, initial, word, edges, baseline=0, global_ports=None):
    x = list(map(F, initial))
    initial = list(x)
    writers = [f"init:{i}" for i in range(len(x))]
    touches = [0]*len(x)
    events = []
    total_loss = F(0)
    for k, (u, v) in enumerate(word):
        average = (x[u]+x[v])/2
        loss = (x[u]-x[v])**2/2
        events.append({"index": k, "op": "mean", "edge": [u, v],
                       "inputs": [str(x[u]), str(x[v])],
                       "writers": [writers[u], writers[v]],
                       "output": str(average), "loss": str(loss)})
        x[u] = x[v] = average
        writers[u] = writers[v] = f"event:{k}"
        touches[u] += 1
        touches[v] += 1
        total_loss += loss
    n, m = len(x), len(events)
    return {"name": name, "kind": "finite_pair_mean_execution",
            "baseline": str(F(baseline)), "global_ports": global_ports,
            "allowed_edges": [list(e) for e in sorted(set(word if edges is None else edges))],
            "initial": list(map(str, initial)), "events": events,
            "final": list(map(str, x)), "final_writers": writers,
            "metrics": {"registers": n, "initializations": n, "means": m,
                        "reads": 2*m, "writes_including_initialization": n+2*m,
                        "touches": touches, "total_load": str(sum(x)),
                        "initial_quadratic": str(sum((a-baseline)**2 for a in initial)),
                        "final_quadratic": str(sum((a-baseline)**2 for a in x)),
                        "loss": str(total_loss)}}


def fibre_certificate():
    x = list(map(F, [1, 2, 3, 4]))
    witnesses = []
    for u, v, amount in ((3,2,4), (2,1,4), (1,0,4), (2,1,3), (1,0,3), (1,0,2)):
        y = x.copy()
        y[u] -= amount
        y[v] += amount
        post = x.copy()
        post[u] = post[v] = (x[u]+x[v])/2
        witnesses.append({"edge": [u,v], "amount": str(amount),
                          "left": list(map(str,x)), "right": list(map(str,y)),
                          "common_image": list(map(str,post))})
        x = y
    return {"kind": "equal_image_witnesses_not_a_trajectory",
            "initial": ["1","2","3","4"], "concentrated": list(map(str,x)),
            "witnesses": witnesses, "separate_prepared_inputs": 12,
            "means_across_separate_inputs": 12}


def build():
    cases = []
    for m in (0, 1, 4, 12, 32):
        word = [(0,j) for j in range(1,m+1)]
        cases.append(trace(f"cold_star_{m}", [F(3,2)]+[0]*m, word, word))
    cases.append(trace("reused_zero", [F(3,2),0], [(0,1)]*8, [(0,1)]))
    cases.append(trace("prepared_balanced_copy_and_cleanup", [5,3,4,4,6,2],
                       [(2,4),(3,5),(2,3)], None, baseline=4))
    word = [(0,j) for j in range(1,6)]
    cases.append(trace("captured_five_neighbours", [F(3,2)]+[0]*5, word, word,
                       global_ports=[0,1,5,7,10,11]))
    return {"schema": "oph.source_passive_memory.controls.v1",
            "source_sha256": codec.source_pins(), "support": graph_census(),
            "executions": cases, "fibre": fibre_certificate()}


def main():
    destination = codec.HERE/"controls.json"
    destination.write_bytes(codec.canonical(build()))
    print(f"Wrote {destination.name}; run the independent verifier next.")


if __name__ == "__main__":
    main()

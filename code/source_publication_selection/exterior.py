"""Producer: positive exterior-power dynamics, independent of scalar replay."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations

from . import codec


def case(size, sources, horizon):
    indices = tuple(combinations(range(size), len(sources)))
    slot = {s: i for i, s in enumerate(indices)}
    fixed = slot[tuple(sources)]
    transfers = []
    for edge in range(size-1):
        move = []
        for ports in indices:
            present = (edge in ports)+(edge+1 in ports)
            if present == 2:
                move.append(None)
            elif present == 0:
                move.append((slot[ports], slot[ports]))
            else:
                replacement = tuple(sorted((set(ports)-{edge, edge+1}) |
                                          ({edge+1} if edge in ports else {edge})))
                move.append((slot[ports], slot[replacement]))
        transfers.append(move)
    levels = [{"words": 0, "rank_lost": 0, "minimum_fixed": None,
               "sum_fixed": F(0), "hash": sha256()} for _ in range(horizon+1)]

    def walk(word, minors):
        n = len(word)
        item = levels[n]
        values = [F(v, 2**n) for v in minors]
        item["words"] += 1
        item["rank_lost"] += not any(minors)
        item["minimum_fixed"] = values[fixed] if item["minimum_fixed"] is None else min(item["minimum_fixed"], values[fixed])
        item["sum_fixed"] += values[fixed]
        item["hash"].update(codec.canonical([list(word), list(map(str, values))]))
        if n == horizon:
            return
        for edge, move in enumerate(transfers):
            child = tuple(0 if t is None else minors[t[0]]+minors[t[1]] for t in move)
            walk((*word, edge), child)

    walk((), tuple(int(i == fixed) for i in range(len(indices))))
    return {"ports": size, "sources": list(sources), "horizon": horizon,
            "minor_indices": [list(i) for i in indices],
            "levels": [{"horizon": n, "words": item["words"], "rank_lost": item["rank_lost"],
                        "minimum_fixed": str(item["minimum_fixed"]), "sum_fixed": str(item["sum_fixed"]),
                        "census_sha256": item["hash"].hexdigest()} for n, item in enumerate(levels)]}

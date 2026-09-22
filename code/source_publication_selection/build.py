"""Exact rational response-tree producer; no stored expanded histories."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import isqrt

from . import codec


def rank(rows):
    nonzero = [r for r in rows if any(r)]
    if not nonzero:
        return 0
    a = nonzero[0]
    return 2 if any(a[0]*b[1] != a[1]*b[0] for b in nonzero) else 1


def native(weights, horizon):
    """Traverse the response tree once, accumulating each complete level."""
    if (not isinstance(weights, (tuple, list)) or len(weights) != 3 or
            any(type(w) is not int or w <= 0 for w in weights) or
            type(horizon) is not int or not 0 <= horizon <= 12):
        raise ValueError("native specification")
    p = tuple(F(w, sum(weights)) for w in weights)
    h0 = p[1]/(p[0]+p[1])
    levels = [{"words": 0, "published": 0, "viable": 0, "H": F(0),
               "selected_published": F(0), "cylinders": F(0), "first": [F(0)]*3,
               "nodes": []} for _ in range(horizon+1)]

    def walk(state, records, word, probability):
        n = len(word)
        pub = rank(records) == 2
        live = rank((*records, *state)) == 2
        first = next((a for a in word if a != 2), None)
        h = h0 if first is None else F(first == 1)
        if (h > 0) != live:
            raise ValueError("native invariant failed")
        q = probability*h/h0
        item = levels[n]
        item["words"] += 1
        item["published"] += pub
        item["viable"] += live
        item["H"] += probability*pub
        item["selected_published"] += q*pub
        item["cylinders"] += q
        if word and pub:
            item["first"][word[0]] += probability
        item["nodes"].append([list(word), int(pub), int(live), str(h), str(q)])
        if n == horizon:
            return
        for edge in range(3):
            child = list(state)
            row = tuple((state[edge][j]+state[edge+1][j])/2 for j in range(2))
            child[edge] = child[edge+1] = row
            walk(tuple(child), (*records, child[3]), (*word, edge), probability*p[edge])

    walk(((F(1), F(0)), (F(0), F(1)), (F(0), F(0)), (F(0), F(0))),
         ((F(0), F(0)),), (), F(1))
    result = []
    for n, item in enumerate(levels):
        result.append({"horizon": n, "words": item["words"], "published": item["published"],
                       "viable": item["viable"], "mass": str(item["H"]),
                       "selected_published": str(item["selected_published"]),
                       "cylinder_total": str(item["cylinders"]),
                       "conditional_first": [str(x/item["H"]) for x in item["first"]]
                       if item["H"] else None,
                       "census_sha256": codec.digest(item["nodes"])})
    return {"weights": list(weights), "potential": str(h0),
            "equal_prior_menu_posterior": [str(h0/(1+h0)), str(1/(1+h0))],
            "selected_first": ["0", str(p[0]+p[1]), str(p[2])],
            "local_guard_first": ["0", str(p[1]/(p[1]+p[2])), str(p[2]/(p[1]+p[2]))],
            "aggregate_potential": "1", "aggregate_first": list(map(str, p)), "levels": result}


def sign(a, b):
    """Sign of a+b*sqrt(5), by integer squaring with sign separation."""
    if b == 0:
        return (a > 0)-(a < 0)
    if a == 0:
        return (b > 0)-(b < 0)
    if a*b > 0:
        return (a > 0)-(a < 0)
    difference = a*a-5*b*b
    return ((a > 0)-(a < 0))*((difference > 0)-(difference < 0))


def counts(q, population):
    """Tensor cube convolution of exact one-dimensional distance types."""
    if type(q) is not int or q < 2 or population not in ("golden", "grid"):
        raise ValueError("population")
    hist = Counter()
    if population == "golden":
        # coordinate = (u+v*sqrt(5))/2; squared distances have denominator 4.
        points = [(b-2*((b+isqrt(5*b*b))//2), b) for b in range(q)]
        denominator = 4
        for x, y in product(points, repeat=2):
            u, v = x[0]-y[0], x[1]-y[1]
            hist[(u*u+5*v*v, 2*u*v)] += 1
    else:
        denominator = q*q
        for b, c in product(range(q), repeat=2):
            hist[((b-c)**2, 0)] += 1
    totals = [0, 0, 0]
    for (a, b), (c, d), (e, f) in product(hist, repeat=3):
        A, B = a+c+e, b+d+f
        multiplicity = hist[(a, b)]*hist[(c, d)]*hist[(e, f)]
        square = (A*A+5*B*B, 2*A*B)
        comparisons = ((q*square[0]-denominator**2, q*square[1]),
                       (q*A-denominator, q*B),
                       (q**3*square[0]-denominator**2, q**3*square[1]))
        for i, comparison in enumerate(comparisons):
            if sign(*comparison) <= 0:
                totals[i] += multiplicity
    return {"q": q, "population": population, "sites": q**3,
            "all_ordered_pairs": q**6, "ordered_reads_including_self": totals,
            "one_dimensional_distance_types": len(hist)}


def packet():
    from . import exterior, locality
    return codec.seal({"schema": 1, "scope": codec.SCOPE, "source_pins": codec.pins(),
                       "native": [native(w, codec.HORIZON) for w in codec.WEIGHTS],
                       "paths": [exterior.case(*spec) for spec in codec.PATH_CASES],
                       "geometry": [counts(q, p) for q in codec.QS for p in ("golden", "grid")],
                       "locality": {"intervals": [locality.interval(*s) for s in codec.INTERVAL_CASES],
                                    "stencils": [locality.stencil(*s) for s in codec.STENCIL_CASES]}})


def main():
    (codec.HERE/"controls.json").write_bytes(codec.canonical(packet()))


if __name__ == "__main__":
    main()

"""Independent integer scalar experiments and rational-enclosure geometry.

Does not import the producer or its response, radical-sign or tree routines.
Shared code is restricted to serialization, experiment constants and custody.
Finite checks certify finite controls, not infinite-path or full-axiom claims.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import isqrt, prod

from . import codec


def dimension(rows):
    pivot = next((x for x in rows if x != (0, 0)), None)
    if pivot is None:
        return 0
    for x, y in rows:
        if pivot[0]*y-pivot[1]*x:
            return 2
    return 1


def experiment(word):
    """Run two scalar basis preparations using integer-scaled means."""
    if not isinstance(word, (tuple, list)) or any(type(e) is not int or not 0 <= e < 3 for e in word):
        raise ValueError("native word")
    x, y = [1, 0, 0, 0], [0, 1, 0, 0]
    samples = [(0, 0)]
    for e in word:
        for values in (x, y):
            mean_numerator = values[e]+values[e+1]
            values[:] = [2*v for v in values]
            values[e] = values[e+1] = mean_numerator
        samples.append((x[3], y[3]))
    # Each retained row has its own nonzero dyadic denominator; row rank
    # and aggregate span are unchanged by those independent row scalings.
    joint = samples+list(zip(x, y))
    pub = dimension(samples) == 2
    live = dimension(joint) == 2
    if dimension(joint+[(1, 1)]) != dimension(joint):
        raise ValueError("lost conserved total")
    return pub, live


def native(weights, horizon):
    p = [F(w, sum(weights)) for w in weights]
    h0 = F(weights[1], weights[0]+weights[1])
    levels = []
    previous = None
    for n in range(horizon+1):
        nodes, by_prefix = [], Counter()
        H = Qpub = Qtotal = F(0)
        first = [F(0)]*3
        published = viable = 0
        for word in product(range(3), repeat=n):
            pub, live = experiment(word)
            neutral = all(a == 2 for a in word)
            h = h0 if neutral else F(live)
            # Check the claimed infinite-state invariant against native data.
            branch = next((a for a in word if a != 2), None)
            if live != (branch != 0):
                raise ValueError("native viability partition")
            r = prod(p[a] for a in word)
            q = r*h/h0
            published += pub
            viable += live
            H += r*pub
            Qpub += q*pub
            Qtotal += q
            if word:
                by_prefix[word[:-1]] += q
                if pub:
                    first[word[0]] += r
            nodes.append([list(word), int(pub), int(live), str(h), str(q)])
        if Qtotal != 1 or Qpub != H/h0:
            raise ValueError("publication cylinder normalization")
        if previous is not None and dict(by_prefix) != previous:
            raise ValueError("prefix inconsistency")
        previous = {tuple(row[0]): F(row[4]) for row in nodes}
        levels.append({"horizon": n, "words": 3**n, "published": published,
                       "viable": viable, "mass": str(H), "selected_published": str(Qpub),
                       "cylinder_total": str(Qtotal),
                       "conditional_first": [str(v/H) for v in first] if H else None,
                       "census_sha256": codec.digest(nodes)})
    return {"weights": list(weights), "potential": str(h0),
            "equal_prior_menu_posterior": [str(F(weights[1], weights[0]+2*weights[1])),
                                          str(F(weights[0]+weights[1], weights[0]+2*weights[1]))],
            "selected_first": ["0", str(p[0]+p[1]), str(p[2])],
            "local_guard_first": ["0", str(F(weights[1], weights[1]+weights[2])),
                                  str(F(weights[2], weights[1]+weights[2]))],
            "aggregate_potential": "1", "aggregate_first": list(map(str, p)), "levels": levels}


def geometry(q, population):
    # Rational enclosures replace the producer's exact sign-by-squaring test.
    scale = 10**60
    root5 = isqrt(5*scale*scale)
    hist = Counter()
    if population == "golden":
        floors = []
        for b in range(q):
            lo = b*(scale+root5)//(2*scale)
            hi = b*(scale+root5+1)//(2*scale)
            if lo != hi:
                raise ValueError("coordinate floor enclosure")
            floors.append(lo)
        denominator = 4
        for i in range(q):
            for j in range(q):
                d, k = i-j, floors[i]-floors[j]
                hist[(6*d*d-4*d*k+4*k*k, 2*d*d-4*d*k)] += 1
    elif population == "grid":
        denominator = q*q
        # Difference multiplicities, independently of enumerating point pairs.
        for d in range(-(q-1), q):
            hist[(d*d, 0)] += q-abs(d)
    else:
        raise ValueError("population")
    pair = Counter()
    for (a, b), m in hist.items():
        for (c, d), n in hist.items():
            pair[(a+c, b+d)] += m*n
    rootq = isqrt(q*scale*scale)
    thresholds = ((F(scale, rootq+1), F(scale, rootq)), (F(1, q), F(1, q)),
                  (F(scale, q*(rootq+1)), F(scale, q*rootq)))
    totals = [0, 0, 0]
    for (a, b), m in pair.items():
        for (c, d), n in hist.items():
            A, B = a+c, b+d
            ends = (F(A*scale+B*root5, denominator*scale),
                    F(A*scale+B*(root5+1), denominator*scale))
            lo, hi = min(ends), max(ends)
            for k, (tlo, thi) in enumerate(thresholds):
                if hi <= tlo:
                    totals[k] += m*n
                elif lo > thi:
                    pass
                elif k == 1 and B == 0 and q*A == denominator:
                    totals[k] += m*n  # inclusive exact rational boundary
                else:
                    raise ValueError("unresolved radius boundary")
    if sum(hist.values()) != q*q or not q**3 <= totals[2] <= totals[1] <= totals[0] <= q**6:
        raise ValueError("geometry census")
    return {"q": q, "population": population, "sites": q**3, "all_ordered_pairs": q**6,
            "ordered_reads_including_self": totals, "one_dimensional_distance_types": len(hist)}


def verify(packet):
    if type(packet) is not dict or set(packet) != {"schema", "scope", "source_pins", "native", "paths", "geometry", "locality", "sha256"}:
        raise ValueError("packet fields")
    codec.equal(packet["schema"], 1, "schema")
    codec.equal(packet["scope"], codec.SCOPE, "scope")
    codec.equal(packet["source_pins"], codec.pins(), "source pins")
    codec.equal(packet["sha256"], codec.digest({k: v for k, v in packet.items() if k != "sha256"}), "digest")
    from . import check_locality
    codec.equal(packet["locality"], {
        "intervals": [check_locality.interval(*s) for s in codec.INTERVAL_CASES],
        "stencils": [check_locality.stencil(*s) for s in codec.STENCIL_CASES]}, "locality semantics")
    if type(packet["native"]) is not list or len(packet["native"]) != len(codec.WEIGHTS):
        raise ValueError("native cases")
    if type(packet["geometry"]) is not list or len(packet["geometry"]) != 2*len(codec.QS):
        raise ValueError("geometry cases")
    if type(packet["paths"]) is not list or len(packet["paths"]) != len(codec.PATH_CASES):
        raise ValueError("path cases")
    for item, weights in zip(packet["native"], codec.WEIGHTS):
        codec.equal(item, native(weights, codec.HORIZON), "native semantics")
    for item, (q, pop) in zip(packet["geometry"], product(codec.QS, ("golden", "grid"))):
        codec.equal(item, geometry(q, pop), "geometry semantics")
    from . import check_exterior
    for item, spec in zip(packet["paths"], codec.PATH_CASES):
        codec.equal(item, check_exterior.case(*spec), "path minor semantics")
    return {"schema": 1, "verified_controls_sha256": packet["sha256"],
            "native_words_replayed": len(codec.WEIGHTS)*sum(3**n for n in range(codec.HORIZON+1)),
            "path_words_replayed": sum(sum((N-1)**n for n in range(T+1)) for N, _, T in codec.PATH_CASES),
            "geometry_ordered_pairs_counted": 2*sum(q**6 for q in codec.QS),
            "locality_graphs_reconstructed": len(codec.INTERVAL_CASES)+len(codec.STENCIL_CASES),
            "complete_prefix_consistency": True, "M1_derived": False}


def main():
    result = verify(codec.load(codec.HERE/"controls.json"))
    codec.equal(result, codec.load(codec.HERE/"receipt.json"), "receipt")
    print(codec.canonical(result).decode("ascii").strip())


if __name__ == "__main__":
    main()

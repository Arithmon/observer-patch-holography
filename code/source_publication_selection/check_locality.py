"""Independent finite DAG closure and metric grid breadth-first search."""
from fractions import Fraction as F
from itertools import product


def interval(N, K):
    if type(N) is not int or type(K) is not int or min(N, K) < 1:
        raise ValueError("positive integer interval specification")
    layers = [[(0, 0)]] + [[(j, i) for i in range(N)] for j in range(1, K+1)] + [[(K+1, 0)]]
    ancestors = {layers[0][0]: set()}
    for before, after in zip(layers, layers[1:]):
        for event in after:
            ancestors[event] = set(before)
            for predecessor in before:
                ancestors[event].update(ancestors[predecessor])
    events = sum(layers, [])
    assert all(layers[0][0] in ancestors[e] for e in events[1:])
    assert all(e in ancestors[layers[-1][0]] for e in events[:-1])
    comparable = sum(len(a) for a in ancestors.values())
    pairs = sum(1 for i, _ in enumerate(events) for _ in events[:i])
    fraction = F(comparable, pairs)
    assert 1-F(1, K) <= fraction <= 1
    return {"sites_per_layer": N, "interior_layers": K, "events": len(events),
            "comparable_pairs": comparable, "unordered_pairs": pairs,
            "ordering_fraction": str(fraction), "incomparable_pairs": pairs-comparable}


def stencil(q, exponent):
    if type(q) is not int or q < 4 or q % 4 or type(exponent) is not int or exponent not in (1, 2):
        raise ValueError("invalid grid control")
    radius = F(1, q**exponent)
    elapsed = F(3, 4)
    target = (q//2, q//2, 0)
    candidates = product(range(-1, 2), repeat=3)
    moves = [s for s in candidates if F(sum(v*v for v in s), q*q) <= radius*radius]
    seen = {(0, 0, 0): 0}
    frontier = [(0, 0, 0)]
    while frontier:
        next_frontier = []
        for x in frontier:
            for s in moves:
                y = tuple(a+b for a, b in zip(x, s))
                if all(0 <= a <= q for a in y) and y not in seen:
                    seen[y] = seen[x]+1
                    next_frontier.append(y)
        frontier = next_frontier
    distance = seen.get(target)
    steps = elapsed/radius
    assert steps.denominator == 1
    return {"q": q, "radius_exponent": exponent, "elapsed": str(elapsed),
            "steps": steps.numerator, "target": list(target),
            "euclidean_strict_timelike": F(sum(a*a for a in target), q*q) < elapsed*elapsed,
            "stencil_size_including_self": len(moves), "minimum_steps": distance,
            "target_reachable": distance is not None and distance <= steps}

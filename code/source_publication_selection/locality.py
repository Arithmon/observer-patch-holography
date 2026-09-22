"""Closed formulas for removing locality and retaining a fixed grid stencil."""
from fractions import Fraction as F


def interval(N, K):
    if type(N) is not int or type(K) is not int or min(N, K) < 1:
        raise ValueError("positive integer interval specification")
    n = K*N+2
    comparable = N*N*K*(K-1)//2+2*K*N+1
    return {"sites_per_layer": N, "interior_layers": K, "events": n,
            "comparable_pairs": comparable, "unordered_pairs": n*(n-1)//2,
            "ordering_fraction": str(F(2*comparable, n*(n-1))),
            "incomparable_pairs": K*N*(N-1)//2}


def stencil(q, exponent):
    if type(q) is not int or q < 4 or q % 4 or type(exponent) is not int or exponent not in (1, 2):
        raise ValueError("grid multiple of four and exponent one or two required")
    # Unit cube, h=1/q, radius=q^-exponent, duration=radius, c=1.
    # The diagonal target is (1/2,1/2,0). At t=3/4 it is strictly timelike.
    steps = 3*q**exponent//4
    return {"q": q, "radius_exponent": exponent, "elapsed": "3/4",
            "steps": steps, "target": [q//2, q//2, 0],
            "euclidean_strict_timelike": True,
            "stencil_size_including_self": 7 if exponent == 1 else 1,
            "minimum_steps": q if exponent == 1 else None,
            "target_reachable": False}

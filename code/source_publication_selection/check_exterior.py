"""Verifier: integer scalar means and fraction-free determinant elimination."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product

from . import codec


def determinant(matrix):
    """Bareiss elimination; zero pivot uses row exchange, with exact divisions."""
    a = [list(row) for row in matrix]
    n, sign, previous = len(a), 1, 1
    for k in range(n-1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*a[k][k]-a[i][k]*a[k][j]
                if numerator % previous:
                    raise ValueError("nonexact determinant division")
                a[i][j] = numerator//previous
        previous = a[k][k]
        for i in range(k+1, n):
            a[i][k] = 0
    return sign*a[-1][-1]


def case(size, sources, horizon):
    k = len(sources)
    indices = tuple(combinations(range(size), k))
    fixed = indices.index(tuple(sources))
    separated = all(b-a >= 2 for a, b in zip(sources, sources[1:]))
    levels = []
    for n in range(horizon+1):
        commitment, total, minimum, lost = sha256(), F(0), None, 0
        for word in product(range(size-1), repeat=n):
            # Independent scalar experiments, one per preparation coordinate.
            columns = [[int(i == port) for i in range(size)] for port in sources]
            for edge in word:
                for x in columns:
                    v = x[edge]+x[edge+1]
                    x[:] = [2*a for a in x]
                    x[edge] = x[edge+1] = v
            rows = list(zip(*columns))
            values = [F(determinant([rows[i] for i in ports]), 2**(n*k)) for ports in indices]
            if any(v < 0 for v in values):
                raise ValueError("ordered-minor positivity")
            if separated and values[fixed] < F(1, 2**n):
                raise ValueError("permanent separated minor")
            rank_lost = not any(values)
            if not separated:
                if k != 2 or sources[1] != sources[0]+1:
                    raise ValueError("unclassified control")
                middle = sources[0]
                first = next((e for e in word if e in (middle-1, middle, middle+1)), None)
                if rank_lost != (first == middle):
                    raise ValueError("adjacent-record competing-edge classification")
            lost += rank_lost
            minimum = values[fixed] if minimum is None else min(minimum, values[fixed])
            total += values[fixed]
            commitment.update(codec.canonical([list(word), list(map(str, values))]))
        levels.append({"horizon": n, "words": (size-1)**n, "rank_lost": lost,
                       "minimum_fixed": str(minimum), "sum_fixed": str(total),
                       "census_sha256": commitment.hexdigest()})
    return {"ports": size, "sources": list(sources), "horizon": horizon,
            "minor_indices": [list(i) for i in indices], "levels": levels}

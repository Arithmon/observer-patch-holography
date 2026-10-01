"""Construct labelled affine permutations and count length-12 walks exactly."""

import numpy as np
from .format import digest

CATALOG = (2, 3, 5, 8, 11, 16)


def permutations(m):
    if type(m) is not int or not 2 <= m <= 64:
        raise ValueError('finite control size must be an integer from 2 to 64')
    result = []
    for axis in (0, 1):
        for offset in (0, 1):
            for sign in (-1, 1):
                row = []
                for x in range(m):
                    for y in range(m):
                        a, b = x, y
                        if axis == 0:
                            a = (x+sign*(2*y+offset)) % m
                        else:
                            b = (y+sign*(2*x+offset)) % m
                        row.append(m*a+b)
                result.append(row)
    return result


def walks(m, power=12):
    if type(power) is not int or not 1 <= power <= 12:
        raise ValueError('bounded walk power')
    maps = permutations(m)
    a = np.zeros((m*m, m*m), dtype=np.int64)
    for row in maps:
        a[np.arange(m*m), row] += 1
    b = np.eye(m*m, dtype=np.int64)
    for _ in range(power):
        b = b @ a
    return a, b


def summary(m):
    a, b = walks(m)
    n, d = m*m, 8**12
    mixing = int(np.max(np.sum(np.abs(n*b-d), axis=1)))
    return dict(m=m, degree=d, permutation_sha256=digest(permutations(m)),
                adjacency_sha256=digest(a.tolist()), walks_sha256=digest(b.tolist()),
                trace=int(np.trace(b)), maximum=int(b.max()),
                mixing_numerator=mixing, mixing_denominator=n*d)


def damage_cases():
    rows = []
    for m in (8, 16):
        _, b = walks(m)
        n, d = m*m, 8**12
        for mode in ('consecutive', 'spaced', 'column'):
            initial = (list(range(n//16)) if mode == 'consecutive' else
                       list(range(0, n, 16)) if mode == 'spaced' else
                       [m*(j % m)+j//m for j in range(n//16)])
            fresh = [v for v in range(n) if v not in initial][:n//64]
            damaged = list(range(n-1, n-1-n//64, -1))
            union = sorted(set(initial+fresh))
            wrong = np.flatnonzero(2*b[:, union].sum(axis=1) >= d).tolist()
            output = sorted(set(wrong+damaged))
            rows.append(dict(m=m, mode=mode, initial=initial, shared_faults=fresh,
                             damaged_voters=damaged, potential_wrong=wrong,
                             output_wrong=output, maximum_live_wrong=len(union)))
    return rows


def candidate():
    return dict(graphs=[summary(m) for m in CATALOG], damage=damage_cases())

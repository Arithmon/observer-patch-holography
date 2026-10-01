"""Independent edge insertion, matrix powering and exact mixing certificates."""

from functools import lru_cache
import numpy as np
from .format import keys, need, integer, digest


@lru_cache(maxsize=12)
def reference(m, power=12):
    # Four positive affine maps plus their transposes. No producer maps.
    n = m*m
    forward = np.zeros((n, n), dtype=np.int64)
    for v in range(n):
        x, y = divmod(v, m)
        for axis, offset in ((0, 0), (0, 1), (1, 0), (1, 1)):
            target = (((x+2*y+offset) % m)*m+y if axis == 0 else
                      x*m+(y+2*x+offset) % m)
            forward[v, target] += 1
    a = forward+forward.T
    b = np.linalg.matrix_power(a, power)
    a.flags.writeable = b.flags.writeable = False
    return a, b


def expected_maps(m):
    maps = []
    for axis, shift in ((0, 0), (0, 1), (1, 0), (1, 1)):
        positive = []
        for v in range(m*m):
            x, y = divmod(v, m)
            positive.append(((x+2*y+shift) % m)*m+y if axis == 0 else
                            x*m+(y+2*x+shift) % m)
        inverse = [0]*(m*m)
        for v, w in enumerate(positive):
            inverse[w] = v
        maps.extend((inverse, positive))
    return maps


def check_maps(rows, m):
    need(type(rows) is list and len(rows) == 8, 'eight labelled maps')
    for row, expected in zip(rows, expected_maps(m)):
        need(type(row) is list and len(row) == m*m
             and all(type(v) is int for v in row) and row == expected,
             'exact affine permutation and multiplicity')


def verify_graphs(rows):
    need(type(rows) is list and len(rows) == 6, 'complete graph catalog')
    for row, m in zip(rows, (2, 3, 5, 8, 11, 16)):
        keys(row, 'm degree permutation_sha256 adjacency_sha256 walks_sha256 trace maximum '
                  'mixing_numerator mixing_denominator')
        integer(row['m'], m, m)
        integer(row['degree'], 8**12, 8**12)
        a, b = reference(m)
        n, d = m*m, 8**12
        need(np.array_equal(a, a.T) and np.all(a.sum(axis=1) == 8), 'regular symmetric multigraph')
        need(np.all(b.sum(axis=1) == d) and np.array_equal(b, b.T), 'complete powered walks')
        for name, expected in (('permutation_sha256', digest(expected_maps(m))),
                               ('adjacency_sha256', digest(a.tolist())),
                               ('walks_sha256', digest(b.tolist()))):
            need(row[name] == expected, 'independent '+name)
        mixing = max(sum(abs(n*int(x)-d) for x in line) for line in b)
        for name, expected in (('trace', sum(int(b[i, i]) for i in range(n))),
                               ('maximum', int(b.max())), ('mixing_numerator', mixing),
                               ('mixing_denominator', n*d)):
            integer(row[name], expected, expected)
        # Exact row-norm certificate for these finite controls, separate from
        # the imported all-size spectral theorem.
        need(4*mixing <= n*d, 'finite exact mixing certificate')


def verify_damage(rows):
    need(type(rows) is list and len(rows) == 6, 'complete damage catalog')
    for row, (m, mode) in zip(rows, ((m, s) for m in (8, 16)
                                   for s in ('consecutive', 'spaced', 'column'))):
        keys(row, 'm mode initial shared_faults damaged_voters potential_wrong output_wrong maximum_live_wrong')
        integer(row['m'], m, m)
        need(row['mode'] == mode, 'frozen damage mode')
        n = m*m
        for field in ('initial', 'shared_faults', 'damaged_voters', 'potential_wrong', 'output_wrong'):
            need(type(row[field]) is list and all(type(v) is int and 0 <= v < n for v in row[field])
                 and len(set(row[field])) == len(row[field]), 'distinct bounded damage indices')
        initial = (list(range(n//16)) if mode == 'consecutive' else
                   [16*i for i in range(n//16)] if mode == 'spaced' else
                   [(j % m)*m+j//m for j in range(n//16)])
        fresh = sorted(set(range(n))-set(initial))[:n//64]
        damaged = list(reversed(range(n-n//64, n)))
        need(row['initial'] == initial and row['shared_faults'] == fresh
             and row['damaged_voters'] == damaged, 'frozen adversarial catalog')
        e = set(initial) | set(fresh)
        _, b = reference(m)
        wrong = [v for v in range(n) if sum(int(b[v, w]) for w in e)*2 >= 8**12]
        output = sorted(set(wrong) | set(damaged))
        need(row['potential_wrong'] == wrong and row['output_wrong'] == output, 'full integer voter replay')
        integer(row['maximum_live_wrong'], len(e), len(e))
        need(9*len(wrong) <= 4*len(e) and len(output) <= n//16 and 2*len(e) < n,
             'contraction and live majority')


def verify(row):
    keys(row, 'graphs damage')
    verify_graphs(row['graphs'])
    verify_damage(row['damage'])

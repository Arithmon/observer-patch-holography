"""Independent endpoint/Choi review of the typed native observer interface.

Expected images use floor/ceiling endpoint arithmetic and independently
enumerated total-two configurations. No native-map constructor, record
recovery function or retained output supplies the expected channel. One
maximally entangled input fixes every matrix unit of each complete CP branch.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import json

from ..native_program import interpreter as program


def configurations(binary):
    choose = combinations if binary else combinations_with_replacement
    rows = []
    for a, b in choose(range(12), 2):
        state = [0]*12
        state[a] += 1
        state[b] += 1
        rows.append(tuple(state))
    return tuple(sorted(rows))


def review():
    blocks, matrix_units, endpoints = 0, 0, 0
    weights = (25, 25)+(8,)*10
    for shell in ('H78', 'H66'):
        states = configurations(shell == 'H66')
        index = {x: i for i, x in enumerate(states)}
        d = len(states)
        # Normalized Choi input on source x untouched reference.
        initial = {(i, i, j, j): F(1, d) for i in range(d) for j in range(d)}
        for a, b in program.local.SEAMS:
            for law in ('C', 'E'):
                admitted = program.program(shell=shell, commands=[{
                    'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [a, b]},
                    'law': law, 'repeat': 1}])
                actual = program.execute(admitted, initial)
                expected = defaultdict(dict)
                for coin in (0, 1):
                    weight = F(weights[a if coin else b], weights[a]+weights[b])
                    images = defaultdict(list)
                    for i, x in enumerate(states):
                        y = list(x)
                        total = x[a]+x[b]
                        y[a], y[b] = (total+coin)//2, (total+1-coin)//2
                        label = x[a]-y[a] if law == 'C' else x[a]-x[b]
                        images[label].append((i, index[tuple(y)]))
                    for label, pairs in images.items():
                        expected[coin, label] = {
                            (u, i, v, j): weight/F(d)
                            for i, u in pairs for j, v in pairs}
                        matrix_units += len(pairs)**2
                for branch in actual.branches:
                    r = branch.record
                    wanted = (expected[r.coin, r.label]
                              if r.kind == 'repair' and r.law == law else {})
                    if branch.density != wanted:
                        raise ArithmeticError(f'complete Choi block mismatch: {shell},{law},{a,b},{r}')
                    blocks += 1
                if sum(program.trace(x.density) for x in actual.branches) != 1:
                    raise ArithmeticError('independent complete channel lost trace')
                if program.reference_marginal(x.density for x in actual.branches) != {
                        (i, i): F(1, d) for i in range(d)}:
                    raise ArithmeticError('independent Choi reference was changed')
                endpoints += 1
    return {'verdict': 'PASS', 'complete_endpoint_instruments': endpoints,
            'raw_blocks_including_zeros': blocks, 'checked_branch_matrix_units': matrix_units,
            'expected_rule': 'independent floor/ceiling endpoint arithmetic',
            'reference': 'full normalized maximally entangled Choi input',
            'canonical_physical_premises_discharged': 0}


if __name__ == '__main__':
    print(json.dumps(review(), indent=2))

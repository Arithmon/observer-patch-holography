"""Bounded independent controls of the weighted-graph theorem and its premise."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product

import numpy as np
import pytest


GRAPHS = (
    (2, ((0, 1, F(2, 3)),)),
    (3, ((0, 1, F(1, 3)), (2, 1, F(7, 5)))),
    (4, ((0, 1, F(2)), (0, 2, F(1, 2)), (3, 0, F(4, 3)))),
    (5, ((0, 1, F(1, 2)), (1, 2, F(2)), (2, 3, F(3, 4)),
         (3, 4, F(5, 3)), (4, 0, F(1)), (0, 2, F(2, 7)))),
)


def _kraus_column(candidate, edges, x, y):
    # Derive outputs by endpoint floor/ceiling, without a current update
    # formula, signed-difference sectors, or any production module.
    result = defaultdict(F)
    total_weight = sum(w for _, _, w in edges)
    for a, b, weight in edges:
        for coin in (0, 1):
            outputs, records = [], []
            for state in (x, y):
                output = list(state)
                total = state[a] + state[b]
                output[a] = (total + coin) // 2
                output[b] = total - output[a]
                outputs.append(tuple(output))
                records.append(state[a] - output[a] if candidate == "C"
                               else state[a] - state[b])
            if records[0] == records[1]:
                result[tuple(outputs)] += weight / (2 * total_weight)
    return {pair: value for pair, value in result.items() if value}


def _sector_column(candidate, edges, x, y):
    z = tuple(a - b for a, b in zip(x, y, strict=True))
    total_weight = sum(w for _, _, w in edges)
    killing = sum(w * (abs(z[a] - z[b]) if candidate == "C"
                       else 2 * (z[a] != z[b])) for a, b, w in edges)
    result = defaultdict(F, {(x, y): 1 - killing / (2 * total_weight)})
    for a, b, weight in edges:
        if z[a] == z[b] == 0 and x[a] != x[b]:
            left, right = list(x), list(y)
            left[a], left[b] = left[b], left[a]
            right[a], right[b] = right[b], right[a]
            result[tuple(left), tuple(right)] += weight / (2 * total_weight)
            result[x, y] -= weight / (2 * total_weight)
    return {pair: value for pair, value in result.items() if value}


@pytest.mark.parametrize("candidate", ("C", "E"))
@pytest.mark.parametrize("n,edges", GRAPHS)
def test_weighted_graphs_all_nontrivial_loads_and_matrix_units(candidate, n, edges):
    graph_laplacian = np.zeros((n, n))
    for a, b, weight in edges:
        graph_laplacian[a, a] += float(weight)
        graph_laplacian[b, b] += float(weight)
        graph_laplacian[a, b] -= float(weight)
        graph_laplacian[b, a] -= float(weight)
    graph_gap = np.linalg.eigvalsh(graph_laplacian)[1]
    for k in range(1, n):
        states = [tuple(int(i in occupied) for i in range(n))
                  for occupied in combinations(range(n), k)]
        pairs = list(product(states, repeat=2))
        columns = [_kraus_column(candidate, edges, *pair) for pair in pairs]
        for pair, actual in zip(pairs, columns, strict=True):
            assert actual == _sector_column(candidate, edges, *pair)
            assert all(value >= 0 for value in actual.values())
            assert sum(value for (x, y), value in actual.items() if x == y) == (pair[0] == pair[1])
        exact = [[column.get(pair, F(0)) for column in columns] for pair in pairs]
        assert exact == list(map(list, zip(*exact, strict=True)))
        # Corroboration only: the proof imports the weighted exclusion theorem.
        spectrum = np.linalg.eigvalsh(np.array(exact, dtype=float))
        assert spectrum[0] >= -1e-14
        assert spectrum[-1] == pytest.approx(1, abs=1e-14, rel=0)
        assert 1 - spectrum[-2] == pytest.approx(
            graph_gap / (2 * float(sum(w for _, _, w in edges))), abs=1e-14, rel=0)


@pytest.mark.parametrize("candidate", ("C", "E"))
def test_disconnected_graph_has_an_exact_fixed_coherence(candidate):
    edges = ((0, 1, F(1, 3)), (2, 3, F(7, 2)))
    x, y = (1, 1, 0, 0), (0, 0, 1, 1)
    # Each component is full on one side and empty on the other. Neither
    # record distinguishes them; connectedness is essential to erasure.
    for left, right in product((x, y), repeat=2):
        expected = {(left, right): F(1)}
        assert _kraus_column(candidate, edges, left, right) == expected
        assert _sector_column(candidate, edges, left, right) == expected

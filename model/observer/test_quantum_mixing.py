"""Independent endpoint enumeration and exact coherent-sector controls."""
from collections import Counter
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement

import numpy as np
import pytest

from . import quantum_mixing as q


def _independent_events():
    states = []
    for a, b in combinations_with_replacement(range(12), 2):
        x = [0]*12
        x[a] += 1
        x[b] += 1
        states.append(tuple(x))
    states.sort()
    assert tuple(states) == q.STATES
    index = {x: i for i, x in enumerate(states)}
    events = []
    for a, b in q.seams():
        for coin in (0, 1):
            columns = []
            for x in states:
                y = list(x)
                low, high = (x[a]+x[b])//2, (x[a]+x[b]+1)//2
                y[a], y[b] = (low, high) if coin == 0 else (high, low)
                columns.append((index[tuple(y)], x[a]-y[a], x[a]-x[b]))
            events.append(tuple(columns))
    return tuple(events)


def test_every_full_matrix_unit_matches_independent_floor_ceiling_kraus_grammar():
    # No native update, sector reduction or predicted spectral data supplies
    # this expected action. Both candidates retain their own record equality.
    events = _independent_events()
    for candidate, record in (("C", 1), ("E", 2)):
        checked = 0
        for sector in q.sectors():
            block = q.native_block(candidate, sector.difference)
            for column, (left, right) in enumerate(sector.pairs):
                counts = Counter((event[left][0], event[right][0]) for event in events
                                 if event[left][record] == event[right][record])
                expected = {pair: F(count, 60) for pair, count in counts.items()}
                actual = {pair: block[row][column] for row, pair in enumerate(sector.pairs)
                          if block[row][column]}
                assert actual == expected
                checked += 1
        assert checked == 6084


def test_all_blocks_reduce_exactly_and_full_piles_are_only_transient():
    result = q.audit()
    assert result["full_block_sizes"] == {"1": 4422, "12": 132, "78": 1}
    assert result["binary_block_sizes"] == {"1": 2970, "10": 132, "66": 1}
    assert result["binary_matrix_units"] == 4356
    assert result["binary_diagonal_dimension"] == 66
    assert result["binary_centered_diagonal_dimension"] == 65
    assert result["binary_offdiagonal_dimension"] == 4290
    assert result["candidates"]["C"]["binary_offdiagonal_spectral_radius"] == "5/6"
    assert result["candidates"]["E"]["binary_offdiagonal_spectral_radius"] == "7/10"
    for candidate in ("C", "E"):
        row = result["candidates"][candidate]
        assert row["transient_matrix_units"] == 78**2-66**2
        assert F(row["transient_spectral_radius"]) == F(5, 6)
        assert F(row["binary_gershgorin_lower_bound"]) > 0


def test_every_transient_survival_and_transfer_bound_independently():
    events = _independent_events()
    binary = [all(n in (0, 1) for n in state) for state in q.STATES]
    for candidate, record, coherence_max, feed_count in (
        ("C", 1, F(23, 30), 1260), ("E", 2, F(7, 10), 60),
    ):
        survivals, offdiagonal, transfers = [], [], 0
        for left in range(78):
            for right in range(78):
                if binary[left] and binary[right]:
                    continue
                counts = Counter((event[left][0], event[right][0]) for event in events
                                 if event[left][record] == event[right][record])
                for (a, b), count in counts.items():
                    assert count > 0
                    assert (a, b) == (left, right) or (binary[a] and binary[b])
                survival = F(counts[left, right], 60)
                survivals.append(survival)
                if left != right:
                    offdiagonal.append(survival)
                transfers += sum(binary[a] and binary[b] for a, b in counts)
        assert len(survivals) == 1728
        assert max(survivals) == F(5, 6)
        assert max(offdiagonal) == coherence_max
        assert transfers == feed_count


def test_deleted_graph_blocks_have_the_claimed_exact_coherence_witness():
    edges = {frozenset(edge) for edge in q.seams()}
    for a, b in combinations(range(12), 2):
        delta = tuple((i == a)-(i == b) for i in range(12))
        sector = next(s for s in q.sectors(binary=True) if s.difference == delta)
        common = [next(k for k in range(12) if q.STATES[i][k] and q.STATES[j][k])
                  for i, j in sector.pairs]
        for candidate in ("C", "E"):
            rate = F(5, 6) if candidate == "C" else F(2, 3)+F(frozenset((a, b)) in edges, 30)
            expected = [[F(0)]*10 for _ in range(10)]
            for i, k in enumerate(common):
                neighbors = [j for j, l in enumerate(common) if frozenset((k, l)) in edges]
                expected[i][i] = rate-F(len(neighbors), 60)
                for j in neighbors:
                    expected[i][j] = F(1, 60)
            actual = q.native_block(candidate, delta, binary=True)
            assert actual == tuple(map(tuple, expected))
            # Sum_k |a,k><b,k| is an exact eigenoperator, independent of
            # the numerical diagonalization used for the separate control.
            assert all(sum(row, F(0)) == rate for row in actual)


def test_binary_diagonal_is_the_same_exclusion_chain_and_gap_is_corroborated():
    c = q.native_block("C", q.ZERO, binary=True)
    assert c == q.native_block("E", q.ZERO, binary=True)
    assert c == tuple(zip(*c, strict=True))
    assert all(sum(row, F(0)) == 1 for row in c)
    eigenvalues = np.linalg.eigvalsh(np.array(c, dtype=float))
    # Diagnostic only: the exact gap statement uses the graph/exclusion
    # theorem and the rational bounds, not this floating eigenproblem.
    assert abs(1-eigenvalues[-2]-(5-np.sqrt(5))/60) < 1e-14
    assert eigenvalues[0] > 0


def test_audit_rejects_a_corrupted_reduction(monkeypatch):
    original = q.binary_reduction

    def wrong(candidate, difference):
        matrix = [list(row) for row in original(candidate, difference)]
        matrix[0][0] += F(1, 60)
        return tuple(map(tuple, matrix))

    monkeypatch.setattr(q, "binary_reduction", wrong)
    with pytest.raises(ArithmeticError, match="differs"):
        q.audit()


@pytest.mark.parametrize("field", (0, 1))
def test_audit_rejects_a_corrupted_native_record(monkeypatch, field):
    original = q.native.native

    def wrong(state, seam, coin):
        output, record = original(state, seam, coin)
        record = list(record)
        record[field] += 1
        return output, tuple(record)

    monkeypatch.setattr(q.native, "native", wrong)
    q._images.cache_clear()
    try:
        # A constant record relabel preserves Kraus equality, but violates
        # the retained primitive's declared coin/current meanings.
        with pytest.raises(ArithmeticError, match="repair record"):
            q.audit()
    finally:
        q._images.cache_clear()


@pytest.mark.parametrize("candidate", (True, "c", "", None, [], 1))
def test_candidate_validation(candidate):
    for function in (q.native_block, q.binary_reduction, q.killing_rate):
        with pytest.raises(ValueError):
            function(candidate, q.ZERO)


@pytest.mark.parametrize("difference", (
    [], [0]*12, (0,)*11, (True,)+(0,)*11, (1,)+(0,)*11, (3, -3)+(0,)*10,
))
def test_difference_validation(difference):
    for function in (q.native_block, q.binary_reduction, q.killing_rate):
        with pytest.raises(ValueError):
            function("C", difference)


@pytest.mark.parametrize("binary", (None, 0, 1, "false", ()))
def test_binary_switch_validation(binary):
    q.sectors(binary=True)  # Validation must precede cached bool/int keys.
    with pytest.raises(ValueError):
        q.sectors(binary=binary)
    with pytest.raises(ValueError):
        q.native_block("C", q.ZERO, binary=binary)


def test_sector_and_matrix_outputs_cannot_mutate_cached_evidence():
    sector = q.sectors()[0]
    with pytest.raises(FrozenInstanceError):
        sector.difference = q.ZERO
    matrix = q.native_block("C", q.ZERO, binary=True)
    with pytest.raises(TypeError):
        matrix[0][0] = F(0)
    pile_difference = (2, -2)+(0,)*10
    assert len(q.native_block("C", pile_difference)) == 1
    with pytest.raises(ValueError):
        q.native_block("C", pile_difference, binary=True)

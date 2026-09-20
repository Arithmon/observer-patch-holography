"""Independent rational probes of the quantified protocol and its hypotheses.

These are test-time controls, not additions to the retained eight executions.
No producer, verifier or codec is imported.
"""
from fractions import Fraction as F
from itertools import product

import pytest


def pair_mean(state, edge):
    u, v = edge
    # Apply a separately constructed rational stochastic matrix to the state.
    n = len(state)
    matrix = [[F(i == j) for j in range(n)] for i in range(n)]
    for i in (u, v):
        matrix[i] = [F(int(j == u) + int(j == v), 2) for j in range(n)]
    return tuple(sum((c*x for c, x in zip(row, state)), F(0)) for row in matrix)


def code(amplitudes, baseline=F(2)):
    return tuple(x for a in amplitudes for x in (baseline+a, baseline-a))


def copy_edges(source, target):
    return ((2*source, 2*target), (2*source+1, 2*target+1))


def execute(state, word):
    for edge in word:
        state = pair_mean(state, edge)
    return state


@pytest.mark.parametrize("amplitude", [F(-5, 4), F(0), F(1, 3)])
def test_all_short_three_cell_walks_preserve_a_spectator(amplitude):
    # All words of at most five hops on three cells, including nontrivial
    # revisits; a fourth, unrelated occupied cell must be left unchanged.
    # This is a declared two-rail complete graph, not a captured W12 layout.
    count = 0
    for depth in range(6):
        for successors in product(range(3), repeat=depth):
            route = (0,)+successors
            if any(s == t for s, t in zip(route, route[1:])):
                continue
            state = code((amplitude, 0, 0, F(-2, 3)))
            current = amplitude
            for source, target in zip(route, route[1:]):
                for edge in copy_edges(source, target):
                    state = pair_mean(state, edge)
                    assert min(state) >= 0 and sum(state) == 16
                current /= 2
                copied = [F(0), F(0), F(0), F(-2, 3)]
                copied[source] = copied[target] = current
                assert state == code(copied)
                state = pair_mean(state, (2*source, 2*source+1))
                copied[source] = 0
                assert state == code(copied)
            expected = [F(0), F(0), F(0), F(-2, 3)]
            expected[route[-1]] = amplitude/2**depth
            assert state == code(expected)
            count += 1
    assert count == 63


@pytest.mark.parametrize("amplitude", [F(-3, 2), F(2, 3)])
def test_copy_edge_order_commutes_but_cleanup_cannot_interleave(amplitude):
    initial = code((amplitude, 0))
    plus, minus = copy_edges(0, 1)
    for clear in ((0, 1), (2, 3)):
        correct = execute(initial, (plus, minus, clear))
        assert execute(initial, (minus, plus, clear)) == correct
        # Cleanup between the two rail means leaves a residual at the cell
        # meant to be blank and invalidates the balanced-state invariant.
        for early, late in ((plus, minus), (minus, plus)):
            wrong = execute(initial, (early, clear, late))
            assert wrong != correct
            assert (wrong[clear[0]], wrong[clear[1]]) != (2, 2)


def test_same_cell_hop_and_unequal_blank_baselines_are_not_the_protocol():
    initial = code((F(3, 2), 0))
    assert execute(initial, (*copy_edges(0, 0), (0, 1))) == (2, 2, 2, 2)
    # Both destination rails agree, but their baseline differs. Copy and
    # cleanup produce blanks at a new common mode, not the declared b=2.
    wrong_blank = (F(7, 2), F(1, 2), F(3), F(3))
    final = execute(wrong_blank, (*copy_edges(0, 1), (0, 1)))
    assert final[:2] == (F(5, 2), F(5, 2))
    assert final != code((0, F(3, 4)))


def test_algebraic_distinct_walk_does_not_certify_a_physical_edge():
    state = code((F(3, 2), 0, 0))
    jump = (*copy_edges(0, 2), (0, 1))
    # The algebra works, but neither copy edge belongs to the nearest-neighbor
    # ladder. Lean's GoodWalk checks consecutive inequality, not adjacency.
    assert execute(state, jump) == code((0, 0, F(3, 4)))
    ladder = {(0, 1), (2, 3), (4, 5), (0, 2), (1, 3), (2, 4), (3, 5)}
    assert sum(tuple(sorted(e)) not in ladder for e in jump) == 2


@pytest.mark.parametrize("sign", [-1, 1])
def test_adversarial_readout_margin_is_strict(sign):
    for depth in (0, 1, 14, 15, 64):
        amplitude = F(3, 2**(depth+1))
        actual = code((sign*amplitude,))
        for error in (amplitude/2, amplitude, 3*amplitude/2):
            observed = (actual[0]-sign*error, actual[1]+sign*error)
            if error < amplitude:
                assert (observed[0] > observed[1]) == (sign > 0)
            elif error == amplitude:
                assert observed == (2, 2)
            else:
                assert (observed[0] > observed[1]) != (sign > 0)

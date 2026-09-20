"""Exact counterexamples to dropping the hypotheses of the M1 results.

These controls use rational matrices, independently of both tape engines.
They illustrate the boundaries; the universal positive results are Lean proofs.
"""
from fractions import Fraction as F
from itertools import product


def mean_matrix(n, u, v):
    matrix = [[F(i == j) for j in range(n)] for i in range(n)]
    for i in {u,v}:
        matrix[i] = [F(int(j == u)+int(j == v), 2) for j in range(n)]
    return matrix


def apply(matrix, state):
    return [sum((a*x for a,x in zip(row,state)), F(0)) for row in matrix]


def test_signed_loads_can_reset_to_native_zero():
    initial = [F(1),F(-1)]
    final = apply(mean_matrix(2,0,1), initial)
    assert final == [0,0]
    assert initial[0] > 2*final[0]  # The nonnegative hypothesis is indispensable.


def test_encoded_zero_can_be_reached_with_positive_raw_loads():
    final = apply(mean_matrix(2,0,1), [F(5),F(3)])
    assert final == [4,4]
    assert min(final) > 0
    assert final[0]-final[1] == 0  # Logical zero is not native scalar zero.


def test_connected_restricted_code_can_distinguish_equal_totals():
    # The restricted code is dyadic nonnegative triples of total 1, together
    # with the uniform triple. Means keep dyadics dyadic; 1/3 is not dyadic.
    # The exact-uniform detector is invariant on this closed code, not on the
    # whole positive cone. Exhaustive finite controls illustrate that boundary.
    uniform = [F(1,3)]*3
    seed = [F(1),F(0),F(0)]
    edges = [(0,1),(1,2),(0,2)]  # Connected, with every edge available.
    for length in range(5):
        for word in product(edges, repeat=length):
            state = seed
            for u,v in word:
                matrix = mean_matrix(3,u,v)
                state = apply(matrix,state)
                assert apply(matrix,uniform) == uniform
            assert sum(state) == sum(uniform) == 1
            assert state != uniform
            assert all(x.denominator & (x.denominator-1) == 0 for x in state)
    outside_code = [F(2,3),F(0),F(1,3)]
    assert outside_code != uniform
    assert apply(mean_matrix(3,0,1),outside_code) == uniform
    # This last equality violates global invariance of the uniform detector.


def test_disconnected_components_can_preserve_distinct_records():
    left, right = list(map(F,[1,1,2,2])), list(map(F,[2,2,1,1]))
    readout = lambda state: sum(state[:2])**2
    assert sum(left) == sum(right)
    assert readout(left) != readout(right)
    for u,v in [(0,1),(2,3)]:
        matrix = mean_matrix(4,u,v)
        # Each component total is fixed for any input: its row-functional
        # times the mean matrix is the same functional [1,1,0,0].
        assert [matrix[0][j]+matrix[1][j] for j in range(4)] == [1,1,0,0]
        for state in (left,right):
            assert readout(apply(matrix,state)) == readout(state)

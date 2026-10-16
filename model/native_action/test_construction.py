"""All-input action and locality controls, including invalid global feedback."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import construction as c
from ..federated_quantum_sector import construction as q


def test_exact_positive_local_action_and_nontrivial_interaction():
    result = c.local_action_certificate()
    assert result['full_four_port_matrix_units_per_candidate'] == 256
    # One endpoint population can change the other: this is not a sum of
    # independent on-site generators. The empty joint state remains fixed.
    for law in ('C', 'E'):
        assert c.action_column(law, 0, 0, ((0, 1),)) == {}
        assert c.action_column(law, 1, 1, ((0, 1),)) == {(1, 1): F(1, 2), (2, 2): F(-1, 2)}
        assert c.action_column(law, 3, 3, ((0, 1),)) == {}
    assert result['vacuum_full_coherence_is_fixed']


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_tensor_kraus_entries_match_native_outcome_semantics(law):
    groups = c.full_groups(law, (1, 2), ports=3)
    for left, right in product(range(8), repeat=2):
        for j, group in zip((-1, 0, 1), groups, strict=True):
            actual = defaultdict(F)
            for weight, matrix in group:
                for a, b in product(range(8), repeat=2):
                    value = weight*matrix[a][left]*matrix[b][right]
                    if value:
                        actual[a, b] += value
            assert dict(actual) == q.native_unit(law, left, right, (1, 2), j)


def _column_history(law, left, right, schedule, truncate):
    states = {((), left, right): F(1)}
    for edge in schedule:
        next_states = defaultdict(F)
        for (word, a, b), weight in states.items():
            if truncate and edge == 2:
                next_states[word, a, b] += weight
                continue
            for current in ((-1, 0, 1) if edge == 0 else (None,)):
                for (i, j), value in q.native_unit(law, a, b, (edge, edge+1), current).items():
                    newword = word+(current,) if edge == 0 else word
                    sign = (-1)**((i & 1)+(j & 1)) if edge == 0 and current == 0 else 1
                    next_states[newword, i, j] += weight*value*sign
        states = {key: value for key, value in next_states.items() if value}
    result = defaultdict(F)
    for (word, a, b), value in states.items():
        if a//4 == b//4:
            result[word, a % 4, b % 4] += value
    return {key: value for key, value in result.items() if value}


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_input_matrix_units_good_schedule_complete_history_equality(law):
    # The record-controlled phase is included, not just an unread channel.
    edges = ((0, 1), (1, 2), (2, 3))
    good = 0
    for schedule in product(range(3), repeat=3):
        if 3 in c.causal_cone(tuple(edges[e] for e in schedule), (0, 1)):
            continue
        good += 1
        for left, right in product(range(16), repeat=2):
            assert _column_history(law, left, right, schedule, False) == _column_history(
                law, left, right, schedule, True)
    assert 0 < good < 27


def test_complete_runtime_outputs_entangled_reference_and_checkpoint():
    rows = c.locality_execution()
    assert {row['candidate'] for row in rows} == {'C', 'E'}
    for row in rows:
        assert row['reference_dimension'] == 2
        assert row['largest_bad_schedule_trace_distance_diagnostic'] > 0
        assert row['largest_good_schedule_residual_diagnostic'] < 1e-12
        assert row['serialized_checkpoint_residual_diagnostic'] < 1e-12


def test_poissonization_preserves_actual_channels_not_the_deterministic_count():
    # A C local mode has retention1/2, so after n attempts it is2^-n.
    # Averaging the actual counts gives exp(-s/2), not exp(-n log2) at n=s.
    for s in (F(1, 8), F(1, 2), F(1)):
        partial = sum(math.exp(-float(s))*float(s)**n/math.factorial(n)*2**(-n) for n in range(20))
        assert abs(partial-math.exp(-float(s)/2)) < 1e-14
    assert abs(math.exp(-.5)-.5) > .1
    # E has killed modes: its discrete channel cannot equal any finite matrix
    # exponential, and the source map is not a unitary evolution.
    assert q.native_unit('E', 1, 0, (0, 1)) == {}
    assert q.native_unit('E', 1, 1, (0, 1)) == {(1, 1): F(1, 2), (2, 2): F(1, 2)}


def test_exact_bounds_zero_attempts_radius_and_size_independence():
    assert c.discrete_local_bound(1, 2, 3, 2, 2) == 0
    assert c.poisson_local_bound(2, 6, F(1, 4), 8) == F(243, 1146880)
    assert c.poisson_local_bound(2, 6, 0, 8) == 0
    # Conditional discrete count bound tends monotonically toward its common
    # fixed per-edge intensity bound, without inserting a physical clock.
    poisson = c.poisson_local_bound(2, 6, F(1, 4), 8)
    counts = [c.discrete_local_bound(2, 6, m, m//4, 8) for m in (400, 800, 1600)]
    assert counts[0] < counts[1] < counts[2] < poisson


def test_native_order_is_essential_and_remote_controller_defeats_locality():
    assert 3 in c.causal_cone(((2, 3), (1, 2)), (0, 1))
    assert 3 not in c.causal_cone(((1, 2), (2, 3)), (0, 1))
    result = c.feedback_counterexample()
    assert result['zero_native_attempt_local_trace_distance'] == '1'
    assert result['executed_remote_read_and_controlled_local_flip']


@pytest.mark.parametrize('bad', (True, -1, .25, '1/4', float('nan')))
def test_nonexact_or_negative_clock_parameter_rejected(bad):
    with pytest.raises(ValueError):
        c.poisson_local_bound(2, 6, bad, 4)


@pytest.mark.parametrize('args', ((0, 6, 30, 2, 1), (2, 31, 30, 2, 1),
    (2, 6, 0, 2, 1), (2, 6, 30, -1, 1), (2, 6, 30, True, 1), (2, 6, 30, 2, -1)))
def test_discrete_parameter_contract(args):
    with pytest.raises(ValueError):
        c.discrete_local_bound(*args)

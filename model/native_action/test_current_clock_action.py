"""Exact stopped-clock moments and complete native current/clock executions."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as sp

from . import current_clock_action as c


def test_palm_occupation_identity_and_uniform_waiting_tail():
    source = c.stopped_source()
    q0, q1 = source['q'][0][0], source['q'][1][0]
    one = sp.ones(66, 1)
    assert source['nu']*source['inverse'] == 6*source['mu']
    assert source['nu']*source['transition'] == source['nu']
    assert source['transition']*one == one
    assert max(q0**3*one) == sp.Rational(2963, 3200)
    assert (q0+q1)*one == one
    assert len(source['Palm_support']) == 51
    assert set(source['local_degree']) == {0, 1, 2, 4, 5}
    # Post-tick preparation is neither faithful nor a uniform arbitrary start.
    assert source['nu'] != sp.ones(1, 66)/66
    assert source['nu'] != source['mu']


def test_local_tick_occupation_flips_and_tick_chain_is_connected():
    source = c.stopped_source()
    states = source['kernel'].states
    signed_occupation = sp.Matrix([2*state[0]-1 for state in states])
    transition = source['transition']
    assert transition*signed_occupation == -signed_occupation
    support = source['Palm_support']
    # Every positive transition is bipartite; a two-step return exists. Thus
    # the irreducible chain has period two, not aperiodic pointwise mixing.
    for i in support:
        for j in support:
            if transition[i, j]:
                assert states[i][0] != states[j][0]
    assert any(transition[i, j]*transition[j, i] > 0 for i in support for j in support)
    for initial in (support[0], support[-1]):
        reached = {initial}
        for _ in support:
            reached |= {j for i in tuple(reached) for j in support if transition[i, j]}
        assert reached == set(support)


def test_complete_cycle_current_has_zero_conditional_drift_before_stopping():
    source = c.stopped_source()
    projection = sp.Matrix(c.circulation.exact_geometry()['cycle_projection'])
    degree = c.event.embedded()[0]
    drift = sp.zeros(66, 30)
    for i, _, edge, _, current, _ in c.circulation.source()[3]:
        drift[i, edge] += sp.Rational(current, degree[i])
    assert drift*projection == sp.zeros(66, 30)
    # Actual stopped occupation, not multiplication by a mean waiting time
    # under an arbitrary preparation, gives every edge second moment 1/5.
    occupation = source['nu']*source['inverse']
    second = [sp.Rational(0)]*30
    for i, _, edge, _, current, _ in c.circulation.source()[3]:
        second[edge] += occupation[i]*sp.Rational(current**2, degree[i])
    assert set(second) == {sp.Rational(1, 5)}
    assert projection*sp.diag(*second)*projection == projection/5


def test_stopped_moments_and_boundary_clock_parity():
    source = c.stopped_source()
    moments = c.stopped_moments(source)
    assert moments[0] == sp.ones(66, 1)
    assert moments[1] == sp.zeros(66, 1)
    assert (source['nu']*moments[2])[0] == sp.Rational(3, 5)
    assert (sp.ones(1, 66)*moments[2]/66)[0] != sp.Rational(3, 5)
    parity = sp.Matrix([2*state[0]-1 for state in source['kernel'].states])
    marked = c.stopped_moments(source, boundary=parity)
    for moment, endpoint_marked in zip(moments, marked):
        assert endpoint_marked == -sp.diag(*parity)*moment


def test_exact_finite_blocks_are_neither_gaussian_nor_independent():
    result = c.local_tick_law()
    assert result['one_local_tick_triangle_moments_0_through_4'][2] == '3/5'
    assert F(result['one_local_tick_triangle_fourth_cumulant']) > F(1, 4)
    assert F(result['adjacent_local_tick_squared_current_covariance']) > F(1, 50)
    assert result['Gaussian_limit_current_action_coefficient'] == '5/2'
    assert result['Gaussian_limit_scalar_triangle_action_coefficient'] == '5/6'
    assert not result['finite_block_increments_independent_or_Gaussian']


def test_independent_floor_ceiling_counts_match_movement_clock_source():
    source = c.stopped_source()
    counts = c.exact_joint_columns()
    degree = c.event.embedded()[0]
    q = [[sp.zeros(66) for _ in range(3)] for _ in range(2)]
    for (moving, tick, score), matrix in zip(c.MARKS, counts):
        if not moving:
            continue
        for i in range(66):
            for j, count in enumerate(matrix[i]):
                if count:
                    for order in range(3):
                        q[tick][order][i, j] += sp.Rational(count*score**order, degree[i])
    assert all(q[b][r] == source['q'][b][r] for b in range(2) for r in range(3))


@pytest.mark.parametrize('law', ('C', 'E'))
def test_joint_refinement_is_trace_preserving_on_full_native_sector(law):
    for controlled in (False, True):
        instrument = c.joint_instrument(law, feedback_seam=controlled)
        total = sum(k.conj().T@k for group in instrument.kraus for k in group)
        np.testing.assert_allclose(total, np.eye(66), atol=2e-15, rtol=0)


def test_runtime_keeps_herald_failure_fixed_window_and_feedback_dependence():
    rows = c.finite_execution()
    for row in rows:
        assert row['herald_and_two_attempt_complete_records'] == 98
        assert sum(F(value) for _, value in row['complete_native_record_probabilities']) == 1
        assert row['preparation_failure_mass'] == '193/198'
        assert row['conditional_no_local_tick_timeout_after_two_native_attempts'] == '2033/2250'
        assert row['feedback_two_attempt_complete_records'] == 49
        assert row['feedback_vs_unchanged_schedule_record_TV'] == '28/5625'
        assert row['feedback_generating_function_E_2_power_current_3_power_ticks'] == '153281/120000'
        assert row['complete_terminal_densities_and_serialized_continuation_checked']
    assert rows[0]['complete_native_record_probabilities'] == rows[1]['complete_native_record_probabilities']


@pytest.mark.parametrize('port', (True, -1, 12, 1.0, '0'))
def test_invalid_local_clock_fails(port):
    with pytest.raises(ValueError):
        c.stopped_source(port)


@pytest.mark.parametrize('maximum', (True, -1, 5, 2.0))
def test_invalid_stopped_moment_order_fails(maximum):
    with pytest.raises(ValueError):
        c.stopped_moments(c.stopped_source(), maximum)


def test_invalid_terminal_boundary_and_instrument_choices_fail():
    with pytest.raises(ValueError):
        c.stopped_moments(c.stopped_source(), boundary=[1])
    with pytest.raises(ValueError):
        c.joint_instrument('X')
    with pytest.raises(ValueError):
        c.joint_instrument(herald=1)

"""Clock conversion keeps remote transitions, complete outputs and controls."""
from itertools import product

import numpy as np
import pytest
import sympy as sp

from . import clock_change as cc
from . import construction as c


@pytest.mark.parametrize('law,feedback', product(('C', 'E'), (False, True)))
def test_complete_clock_elimination_is_associative_and_uniformly_terminating(law, feedback):
    result = cc.certificate(law, feedback)
    assert result['symbolic_remote_count_generating_instruments_equal_in_both_elimination_orders']
    assert result['uniform_no_A_probability_after2k_union_events_bound'] == '(2/3)^k'
    assert result['completed_A_instrument_trace_preserving']
    assert result['independent_identity_hold_activity_one_seventh_exactly_cancels']
    b = sum(cc.union_maps(law, feedback)[2:], sp.zeros(16))
    for cap in range(7):
        arrived = sum((sum(cc.coefficient_maps(law, k, feedback), sp.zeros(16))
                       for k in range(cap)), sp.zeros(16))
        assert cc.TRACE*(arrived+b**cap) == cc.TRACE
        assert max((cc.TRACE*b**(2*cap))[5*i] for i in range(4)) <= sp.Rational(2, 3)**cap


def test_mean_rate_does_not_fix_joint_clock_law_and_feedback_changes_the_prediction():
    passive = cc.relative_count_prediction()
    row = passive['event_boundary_predictions'][0]
    assert row['mean_B_ticks_before_next_A_tick'] == '2'
    assert row['variance_B_ticks_before_next_A_tick'] == '16'
    assert row['probabilities_B_counts0_through8'] == [
        '2/3', '0', '1/9', '0', '2/27', '0', '4/81', '0', '8/243']
    controlled = cc.relative_count_prediction(True)['event_boundary_predictions'][0]
    assert controlled['mean_B_ticks_before_next_A_tick'] == '3/4'
    assert controlled['variance_B_ticks_before_next_A_tick'] == '35/16'
    assert controlled['probabilities_B_counts0_through8'][1] == '1/6'
    assert cc.complete_maps('C') != cc.complete_maps('C', True)


def test_actual_reference_assemblages_survive_count_elimination_and_timeout():
    # Work maximally entangled with a four-dimensional inert reference. A
    # probability-only predictor would incorrectly leave each reference I/4.
    state = sp.zeros(16)
    for i, j in product(range(4), repeat=2):
        state[5*i, 5*j] = sp.Rational(1, 4)
    reference_sum = sp.zeros(4)
    for remote_count in range(4):
        for destination, operation in zip((1, 0), cc.coefficient_maps('C', remote_count, True), strict=True):
            output = c._apply_exact(operation, state, 4)
            reference = sp.Matrix(4, 4, lambda r, s: sum(output[4*i+r, 4*i+s] for i in range(4)))
            expected = sp.zeros(16)
            expected[4*destination:4*destination+4, 4*destination:4*destination+4] = reference
            assert output == expected
            if sp.trace(reference):
                assert reference != sp.trace(reference)*sp.eye(4)/4
            reference_sum += reference
    b = sum(cc.union_maps('C', True)[2:], sp.zeros(16))
    timeout = c._apply_exact(b**4, state, 4)
    reference_sum += sp.Matrix(4, 4, lambda r, s: sum(timeout[4*i+r, 4*i+s] for i in range(4)))
    assert reference_sum == sp.eye(4)/4


def test_erasing_remote_motion_changes_completed_clock_output_for_every_candidate():
    for law in ('C', 'E'):
        maps = cc.controlled_maps(law)
        # Hiding remote B events is legitimate only if their CP state updates
        # are retained. Replacing the whole no-A channel by a scalar is wrong.
        no_a = maps[0]+maps[3]+maps[4]
        assert no_a != sp.eye(16)*sp.trace(no_a)/16
        initial = sp.zeros(4); initial[2, 2] = 1
        exact = [sp.trace(c.unvector(m*c.vector(initial))) for m in cc.complete_maps(law)]
        assert exact == [0, 1]
        immediate = [sp.trace(c.unvector(m*c.vector(initial))) for m in maps[1:3]]
        assert immediate == [0, 0]


@pytest.mark.parametrize('law,feedback,plus,timeout', (
    ('C', False, '35/54', '85/108'), ('E', False, '14/27', '85/108'),
    ('C', True, '547/864', '679/864'), ('E', True, '437/864', '679/864')))
def test_primitive_stopped_runtime_keeps_every_state_feedback_and_timeout(law, feedback, plus, timeout):
    result = cc.execute(law, feedback=feedback)
    assert result['complete_physical_record_count'] == 106
    assert result['fixed_runtime_menu_record_count'] == 250
    assert result['A_timeout_probability'] == timeout
    assert result['unconditional_field_plus_probability'] == plus
    assert result['checkpoint_after_first_source_attempt_and_feedback_replayed']
    assert result['full_state_and_checkpoint_residual_diagnostic'] < 1e-12
    for row in result['complete_record_law']:
        native = row['native_attempt_marks']
        assert all(mark not in (1, 2) for mark in native[:-1])
        assert row['B_tick_count'] == sum(mark in (3, 4) for mark in native)
        assert row['A_tick_completed'] == (native[-1] in (1, 2))


@pytest.mark.parametrize('law', ('C', 'E'))
def test_reference_run_continues_actual_state_after_remote_feedback(law):
    result = cc.execute(law, feedback=True, reference=True, cap=2)
    assert result['reference_dimension'] == 2
    assert result['complete_physical_record_count'] == 34
    assert result['complete_probability'] == '1'
    assert result['full_state_and_checkpoint_residual_diagnostic'] < 1e-12


def test_identity_holds_change_native_censoring_but_not_count_generating_maps():
    maps = c.source_maps('C')
    a = sp.Rational(1, 7)
    initial = sp.zeros(4); initial[0, 0] = 1
    no_a = maps[0]+maps[3]+maps[4]
    held = (1-a)*sp.eye(16)+a*no_a
    assert (cc.TRACE*held**3*c.vector(initial))[0] != (cc.TRACE*no_a**3*c.vector(initial))[0]
    # Every remote-count coefficient is invariant after eliminating holds.
    no_union = (1-a)*sp.eye(16)+a*maps[0]
    resolvent = (sp.eye(16)-no_union).inv()
    held_union = tuple(a*m*resolvent for m in maps[1:])
    assert held_union == c.event_maps('C')


@pytest.mark.parametrize('bad', (0, 1, None, 'yes', np.bool_(True)))
def test_cached_feedback_api_does_not_accept_coerced_boolean(bad):
    # Warm valid keys first: bools/integers otherwise collide in tuple caches.
    cc.controlled_maps('C', False); cc.controlled_maps('C', True)
    cc.union_maps('C', False); cc.union_maps('C', True)
    cc.complete_maps('C', False); cc.complete_maps('C', True)
    for function in (cc.controlled_maps, cc.union_maps, cc.complete_maps, cc.certificate):
        with pytest.raises(ValueError, match='Boolean'):
            function('C', bad)


@pytest.mark.parametrize('bad', (True, -1, 65, 1.0, '1'))
def test_remote_count_requires_an_exact_bounded_integer(bad):
    with pytest.raises(ValueError, match='remote count'):
        cc.coefficient_maps('C', bad)


def test_stopped_runtime_cap_and_reference_are_validated_before_execution():
    for cap in (True, 0, 4, 9, 2.0):
        with pytest.raises(ValueError):
            cc.experiment('C', cap=cap)
    with pytest.raises(ValueError):
        cc.experiment('C', reference=1)

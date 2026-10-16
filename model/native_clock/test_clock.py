"""Independent finite clock, scheduler and complete observer controls."""
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest

from . import clock as c


@pytest.fixture(scope='module')
def source():
    kernel = c.transfer.native_counts()
    return kernel, c.lumping(kernel)


def test_every_native_row_certifies_three_state_marked_reduction(source):
    kernel, reduced = source
    assert reduced['sizes'] == (30, 30, 6)
    assert reduced['marked'][0] == ((F(13, 15), F(0), F(0)),
                                     (F(0), F(5, 6), F(0)), (F(0), F(0), F(5, 6)))
    assert reduced['marked'][1] == ((F(1, 15), F(1, 15), F(0)),
                                     (F(1, 15), F(1, 15), F(1, 30)),
                                     (F(0), F(1, 6), F(0)))
    assert reduced['stationary'] == (F(5, 11), F(5, 11), F(1, 11))
    assert reduced['palm'] == (F(2, 5), F(1, 2), F(1, 10))
    assert len(kernel.states) == 66


def test_bit_word_predictions_match_direct_native_path_expansion(source):
    kernel, reduced = source
    # Directly execute each primitive choice; no marked matrix is used here.
    for initial in (kernel.states[0], kernel.states[10], kernel.states[-1]):
        populations = {((), initial): F(1)}
        for _ in range(3):
            future = {}
            for (word, state), mass in populations.items():
                for seam in kernel.seams:
                    for coin in (0, 1):
                        output, record = c.native.native(state, seam, coin)
                        key = word+(abs(record[1]),), output
                        future[key] = future.get(key, F(0))+mass/60
            populations = future
        initial_class = reduced['classes'][kernel.states.index(initial)]
        prior = tuple(F(int(i == initial_class)) for i in range(3))
        for word in product((0, 1), repeat=3):
            direct = sum((mass for (history, _), mass in populations.items() if history == word), F(0))
            assert direct == sum(c.word_law(reduced, word, prior))


@pytest.mark.parametrize('n', [0, 1, 2, 6])
def test_complete_counts_equal_independently_grouped_bit_words(source, n):
    _, reduced = source
    observed = [F(0)]*(n+1)
    for word in product((0, 1), repeat=n):
        observed[sum(word)] += sum(c.word_law(reduced, word))
    probabilities, mean, variance = c.count_law(reduced, n)
    assert tuple(observed) == probabilities
    assert sum(probabilities) == 1
    assert mean == n*F(5, 33)
    assert variance >= 0
    if n > 1:
        assert variance > n*F(5, 33)*F(28, 33)


@pytest.mark.parametrize('start,weight', [('stationary', F(5, 11)), ('palm', F(2, 5))])
def test_waiting_law_is_exact_geometric_mixture_with_unconditioned_tail(source, start, weight):
    _, reduced = source
    probabilities, tail = c.waiting_law(reduced, 16, start=start)
    expected = tuple(weight*F(13, 15)**(n-1)*F(2, 15)
                     +(1-weight)*F(5, 6)**(n-1)*F(1, 6) for n in range(1, 17))
    assert probabilities == expected
    assert tail == weight*F(13, 15)**16+(1-weight)*F(5, 6)**16
    assert sum(probabilities)+tail == 1
    assert tail > 0
    assert c.waiting_law(reduced, 0, start=start) == ((), F(1))


def test_actual_event_palm_law_differs_from_stationary_time_sample(source):
    _, reduced = source
    result = c.interval_statistics(reduced)
    assert result['mean_interval'] == F(33, 5)
    assert result['stationary_start_wait_mean'] == F(147, 22)
    assert result['interval_variance'] == F(951, 25)
    assert result['adjacent_interval_covariance'] == F(9, 100)
    assert result['P_next_1_given_last_1'] == F(23, 150)
    assert result['P_next_1_given_last_two_11'] == F(53, 345)
    # The discrete stationary inspection identity is independent of the
    # three-class computation of the mean residual waiting time.
    mean = result['mean_interval']
    assert result['stationary_start_wait_mean'] == (result['interval_variance']+mean*mean+mean)/(2*mean)


def test_native_attempt_law_has_an_exact_seconds_scale_nonidentifiability_certificate(source):
    _, reduced = source
    audit = c.duration_scale_audit(reduced)
    assert audit['same_exact_record_law_for_every_positive_seconds_per_attempt']
    assert not audit['dimensionful_source_quantity_with_seconds_units']
    assert not audit['seconds_per_attempt_identified']
    assert not audit['matched_public_lab_dataset_identified']
    rows = audit['illustrative_rescalings']
    assert [row['illustrative_seconds_per_attempt'] for row in rows] == ['1/100', '1', '100']
    assert len({tuple(row['exact_native_one_attempt_record_law']) for row in rows}) == 1
    assert [row['implied_stationary_event_rate_per_second'] for row in rows] == [
        '500/33', '5/33', '1/660']
    assert 'one nonzero paired native-attempt difference' in audit['minimal_calibration_if_a_constant_attempt_clock_is_validated']


def test_equal_mean_independent_clock_misses_second_order_measurement(source):
    _, reduced = source
    p = F(5, 33)
    both = sum(c.word_law(reduced, (1, 1)))
    assert both == F(23, 990)
    assert both-p*p == F(1, 3630)
    count, _, _ = c.count_law(reduced, 2)
    assert count != ((1-p)**2, 2*p*(1-p), p*p)


def test_nonuniform_seams_preserve_mean_but_change_visible_correlations(source):
    kernel, _ = source
    control = c.scheduler_controls(kernel)
    assert control['all_30_single_seam_stationary_moving_rates'] == ['5/33']
    assert control['uniform_rate'] == control['double_first_seam_rate'] == '5/33'
    assert control['double_first_seam_mean_Palm_interval'] == '33/5'
    assert F(control['uniform_P_11']) != F(control['double_first_seam_P_11'])
    assert control['global_lazification_preserves_all_stopped_configuration_transitions']
    assert not control['changed_seam_weights_preserve_stopped_transition_law']
    assert control['half_active_rate'] == '5/66'
    assert F(control['half_active_P_11']) == F(control['uniform_P_11'])/4
    counterexample = control['single_seam_reducible_negative_control']
    assert counterexample['never_moving_configurations'] == 46
    assert counterexample['stationary_mean_move_rate'] == '5/33'
    assert counterexample['Palm_mean_interval_attempts'] == '2'
    assert counterexample['reciprocal_stationary_rate'] == '33/5'


def test_clock_packet_executes_full_native_records_and_checkpoints():
    result = c.run()
    execution = result['actual_observer_execution']
    assert result['all_pass']
    assert len(execution['complete_bit_records']) == 32
    assert sum(F(row['exact_probability']) for row in execution['complete_bit_records']) == 1
    assert execution['complete_terminal_states_checked']
    assert execution['all_no_current_records_retained']
    assert execution['terminal_quantum_dimension_per_record'] == 66
    assert execution['checkpoint_replay_cq_distance_diagnostic'] == 0
    for stopping in execution['stopping_laws']:
        mass = sum(map(F, stopping['first_move_at_attempt_1_to_4']))
        assert mass+F(stopping['no_move_through_four_attempts']) == 1
        assert stopping['runtime_residual_diagnostic'] < 3e-14
        assert not stopping['physical_halt_device_implemented']
    assert not result['scope']['physical_duration_or_proper_time_identified']
    assert not result['scope']['new_lab_prediction_or_empirical_validation']
    assert not result['scope']['state_dependent_or_adaptive_scheduler_covered']
    assert set(result['sources']) == set(c.SOURCE_PATHS)
    json.dumps(result, allow_nan=False)


def test_basis_preparation_has_same_future_bit_law_without_stationary_assumption(source):
    kernel, reduced = source
    visible, _ = c.transfer.instruments(kernel)
    state = np.zeros((66, 66), complex)
    state[0, 0] = 1
    program = c.runtime.Program('basis_clock_control', 66, state,
        tuple(c.runtime.Event(str(i), instruments=(visible,)) for i in range(3)))
    history = c.runtime.execute(program)
    prior = tuple(F(int(j == reduced['classes'][0])) for j in range(3))
    for branch in history.branches:
        assert abs(branch.probability-float(sum(c.word_law(reduced, branch.record, prior)))) < 3e-14


@pytest.mark.parametrize('value', [True, -1, 1., 257, None])
def test_invalid_attempt_counts_fail_closed(source, value):
    _, reduced = source
    with pytest.raises(ValueError):
        c.count_law(reduced, value)
    with pytest.raises(ValueError):
        c.waiting_law(reduced, value)


@pytest.mark.parametrize('word', [(True,), (2,), (-1,), [0, 1], (0,)*257])
def test_invalid_native_words_fail_closed(source, word):
    with pytest.raises(ValueError):
        c.word_law(source[1], word)


@pytest.mark.parametrize('value', [F(0), F(-1), F(2), True, .5])
def test_invalid_activity_is_not_a_clock_parameter(source, value):
    with pytest.raises(ValueError):
        c.weighted_native_marks(source[0], (F(1, 30),)*30, activity=value)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1., .01])
def test_failed_or_nonfinite_runtime_diagnostics_cannot_pass(value):
    with pytest.raises(ArithmeticError):
        c._bound(value, 'mutant')

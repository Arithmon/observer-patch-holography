"""Source-independent token controls and complete native observer comparisons."""
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import runtime, transfer as t


def test_native_kernel_equals_independent_lazy_token_exchange_operator():
    kernel = t.native_counts()
    index = {state: i for i, state in enumerate(kernel.states)}
    expected = [[0]*66 for _ in range(66)]
    for i, state in enumerate(kernel.states):
        expected[i][i] = 60
        for a, b in kernel.seams:
            if state[a] != state[b]:
                swapped = list(state)
                swapped[a], swapped[b] = swapped[b], swapped[a]
                expected[i][index[tuple(swapped)]] += 1
                expected[i][i] -= 1
    assert kernel.counts == tuple(map(tuple, expected))
    result = t.exact_audit(kernel)
    assert result['positive_operator_lower_bound'] == 'K >= (2/3)*I'
    assert result['exchange_degree_histogram'] == {'8': 30, '10': 36}
    assert len(result['connected_spanning_tree']) == 65
    assert result['primitive_branches_checked'] == 3960


def test_two_current_bits_equal_independent_path_enumeration_constants():
    # Independently enumerated all 66 initial states and both successive
    # 60-choice native events; neither the transfer matrix nor this producer
    # was used to obtain these four target fractions.
    law = t.exact_record_law(t.native_counts(), 0)
    assert [law[a, b][0] for a in (0, 1) for b in (0, 1)] == [
        F(713, 990), F(127, 990), F(127, 990), F(23, 990)]
    assert law[1, 0][0]+law[1, 1][0] == F(5, 33)
    assert law[1, 1][0]-F(5, 33)**2 == F(1, 3630)


def test_reconstructed_spectrum_matches_native_kernel_and_is_not_a_supplied_oscillator():
    kernel = t.native_counts()
    h, omega, diagnostics = t.spectral_reconstruction(kernel)
    assert diagnostics['gap_per_Euclidean_native_attempt'] > 0
    assert diagnostics['minimum_transfer_eigenvalue'] >= 2/3
    assert diagnostics['exp_minus_H_vs_native_K_residual'] < 2e-14
    assert diagnostics['H_vacuum_residual'] < 2e-15
    np.testing.assert_allclose(omega, np.ones(66)/np.sqrt(66), atol=0, rtol=0)
    assert h.shape == (66, 66)
    # Positive powers of a stochastic K are not a native unitary: contraction
    # of a nonconstant vector remains present in the original update.
    k = np.asarray(kernel.counts, float)/60
    assert np.linalg.norm(k.T@k-np.eye(66), 'fro') > 1


def test_complete_run_checks_current_records_full_terminal_states_and_checkpoint():
    result = t.run()
    assert result['all_pass']
    assert len(result['observer_current_records']) == 3
    for row in result['observer_current_records']:
        assert row['complete_records'] == [[0, 0], [0, 1], [1, 0], [1, 1]]
        assert sum(map(F, row['exact_probabilities'])) == 1
        assert row['single_event_moving_probability'] == '5/33'
        assert row['quantum_output_dimension_per_record'] == 66
        assert row['complete_terminal_density_residual_diagnostic'] < 2e-14
        assert row['runtime_probability_residual_diagnostic'] < 2e-14
        assert row['spectral_transfer_probability_residual_diagnostic'] < 2e-14
    assert result['checkpoint']['complete_cq_distance_diagnostic'] == 0
    assert result['scope']['retained_bit_is_coarse_absolute_current']
    assert not result['scope']['complete_raw_seam_coin_signed_current_dictionary_retained']
    assert not result['scope']['marked_insertions_are_claimed_POVM_effects']
    assert not result['scope']['classical_transfer_selects_complete_quantum_instrument']
    assert not result['scope']['real_time_unitary_identified_with_native_or_physical_dynamics']
    assert set(result['sources']) == set(t.SOURCE_PATHS)
    json.dumps(result, allow_nan=False)


def test_reflection_form_agrees_with_full_marked_history_and_exact_positive_bound():
    result = t._reflection_check(t.native_counts())
    direct = F(result['direct_native_marked_path_form'])
    assert direct == F(result['conditional_future_kernel_form'])
    assert direct >= F(result['exact_positive_lower_bound_for_fixture']) > 0
    assert result['non_Gaussian_stationary_history']


def test_negative_transfer_eigenvalue_can_break_bond_reflection_positivity():
    # A stochastic reversible chain alone is insufficient. The flip kernel
    # has a negative transfer eigenvalue and a negative reflection form.
    flip = np.array([[0., 1.], [1., 0.]])
    f = np.array([1., -1.])
    assert f@flip@f/2 == -1


def test_native_current_mutation_cannot_be_hidden_by_unchanged_configuration_outputs(monkeypatch):
    original = t.native.native
    def wrong_record(state, seam, coin):
        output, record = original(state, seam, coin)
        return output, (coin, -record[1])
    monkeypatch.setattr(t.native, 'native', wrong_record)
    with pytest.raises(ArithmeticError, match='current'):
        t.native_counts()


def test_incomplete_marked_law_fails_before_spectral_reconstruction():
    kernel = t.native_counts()
    corrupted = [list(row) for row in kernel.marked_counts[0]]
    corrupted[0][0] -= 1
    with pytest.raises(ArithmeticError, match='complete native kernel'):
        t.exact_audit(replace(kernel, marked_counts=(tuple(map(tuple, corrupted)), kernel.marked_counts[1])))


@pytest.mark.parametrize('gap', [True, -1, 1., 17, None])
def test_invalid_gap_is_rejected(gap):
    with pytest.raises(ValueError):
        t.exact_record_law(t.native_counts(), gap)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1., 1e-3])
def test_nonfinite_or_failed_numerical_reconstruction_cannot_pass(value):
    with pytest.raises(ArithmeticError):
        t._bounded(value, 'mutant')


def test_actual_runtime_coarse_current_instrument_is_complete():
    kernel = t.native_counts()
    visible, hidden = t.instruments(kernel)
    assert tuple(map(len, visible.kraus)) == (60, 60)
    assert len(hidden.kraus[0]) == 120
    history = runtime.execute(t.program(kernel, 0, (visible, hidden)))
    assert history.probability_mass == pytest.approx(1, abs=2e-14)
    assert all(dict(branch.resources) == {'native_attempts': 2, 'retained_current_bits': 2}
               for branch in history.branches)

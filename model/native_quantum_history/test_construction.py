"""Complete native quantum execution, independent controls and exact bounds."""
from fractions import Fraction as F

import numpy as np
import pytest

from ..observer.scalar_response.experiment import compare_receipt
from ..primitive_observer import runtime
from . import construction as q


@pytest.fixture(scope='module')
def result():
    return q.run()


def test_full_operator_log_uses_all_matrix_units_and_exact_positive_bounds(result):
    for law, lower in (('C', '2/3'), ('E', '11/30')):
        row = result['operator_transfer'][law]
        assert row['all_binary_matrix_units_checked'] == 4356
        assert row['block_dimensions'] == {'1': 2970, '10': 132, '66': 1}
        assert row['exact_lower_operator_bound'] == lower
        assert row['gap_vs_native_exclusion_formula_residual'] < 2e-14
        assert row['block_spectral_reconstruction_residual'] < 2e-14
    assert result['operator_transfer']['C']['exact_native_block_sha256'] != result['operator_transfer']['E']['exact_native_block_sha256']


def test_complete_phase_measurements_execute_exact_laws_and_optimal_readout(result):
    assert len(result['executed_phase_measurements']) == 20
    for row in result['executed_phase_measurements']:
        rate = F(5, 6) if row['candidate'] == 'C' else F(7, 10)
        expected = (1+row['preparation_sign']*rate**row['native_attempts'])/2
        assert list(map(F, row['exact_probabilities'])) == [expected, 1-expected]
        assert row['runtime_probabilities'] == pytest.approx([float(expected), float(1-expected)], abs=2e-14)
        assert row['complete_terminal_density_residual_diagnostic'] < 2e-14
        assert row['exp_minus_n_H_vs_runtime_density_residual_diagnostic'] < 2e-14
        assert row['expanded_tree_resources'].get('native_attempts', 0) == row['native_attempts']
        assert row['expanded_tree_resources']['phase_readouts'] == 1
    best = next(row for row in result['discrimination'] if row['native_attempts'] == 4)
    assert F(best['exact_C_minus_E_plus_probability']) == F(12259, 101250)
    assert F(best['exact_equal_prior_success_probability']) == F(113509, 202500)
    assert best['complete_post_readout_cq_distance'] == pytest.approx(12259/101250, abs=2e-14)
    assert result['optimal_count_for_declared_phase_pair_experiment'] == 4
    assert F(25, 21)**3 < F(9, 5) < F(25, 21)**4


def test_reference_audit_retains_states_even_when_marginal_pointer_laws_coincide(result):
    for row in result['reference_and_checkpoint']:
        assert row['complete_record_probabilities'] == pytest.approx([.5, .5], abs=2e-14)
        assert row['complete_record_probability_residual_diagnostic'] < 2e-14
        assert row['quantum_output_dimension_per_record'] == 132
        assert row['complete_terminal_reference_density_residual_diagnostic'] < 2e-14
        assert row['untouched_reference_marginal_residual_diagnostic'] < 2e-14
        assert row['checkpoint_independently_replayed']
        assert row['checkpoint_complete_cq_distance_diagnostic'] == 0


@pytest.mark.parametrize('law', ['C', 'E'])
def test_exact_blocks_cover_complex_nonpositive_operators_and_references(law):
    # A deliberately non-Hermitian matrix is not a density test. Matrix-unit
    # linearity must hold on the complete complex operator space, including
    # off-diagonal reference operators; no positivity projection is allowed.
    rng = np.random.default_rng(731)
    matrix = rng.normal(size=(132, 132))+1j*rng.normal(size=(132, 132))
    lifted = [np.kron(k, np.eye(2)) for k in q.native_instrument(law).kraus[0]]
    oracle = sum(k@matrix@k.conj().T for k in lifted)
    actual = q.apply_transfer(matrix, law, 1, reference_dimension=2)
    spectral = q.apply_transfer(matrix, law, 1, reference_dimension=2, spectral=True)
    np.testing.assert_allclose(actual, oracle, atol=2e-14, rtol=2e-14)
    np.testing.assert_allclose(spectral, oracle, atol=2e-14, rtol=2e-14)


def test_diagonal_preparation_erases_C_E_phase_discrimination():
    rho, _, read, _ = q.phase_apparatus()
    diagonal = np.diag(np.diag(rho))
    outputs = []
    for law in ('C', 'E'):
        p = q.program(law, 4)
        history = runtime.execute(p, initial=diagonal)
        assert [b.probability for b in history.branches] == pytest.approx([.5, .5], abs=2e-14)
        outputs.append(history)
    assert runtime.complete_distance(*outputs) < 2e-14
    assert np.linalg.norm(read.kraus[0][0].conj().T@read.kraus[0][0]
                          +read.kraus[1][0].conj().T@read.kraus[1][0]-np.eye(66)) < 2e-14


def test_half_step_is_exactly_nonpositive_on_a_positive_diagonal_input(result):
    row = result['fractional_channel_obstruction']
    assert F(row['exact_partial_sum']) == -F(809497755979, 4777574400000000)
    assert F(row['absolute_tail_bound']) == F(1, 13122)
    assert F(row['strictly_negative_upper_bound']) == -F(4008679803811, 42998169600000000)
    assert row['fractional_map_fails_positivity_already_on_a_diagonal_density']
    assert not row['integer_native_channels_affected']


def test_positive_quantum_channel_need_not_be_positive_in_operator_inner_product():
    # Conjugation by Z is CP and trace preserving but gives the negative
    # bond form <X, Z X Z>/2=-1. CP alone cannot justify a positive HS log.
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    assert np.trace(x.T@z@x@z)/2 == -1


def test_portable_comparison_detects_a_changed_measurement_prediction(result):
    import copy
    changed = copy.deepcopy(result)
    changed['executed_phase_measurements'][0]['runtime_probabilities'][0] -= .01
    with pytest.raises((ValueError, AssertionError)):
        compare_receipt(result, changed)


def test_native_record_corruption_fails_before_execution(monkeypatch):
    original = q.native.native
    def corrupt(state, seam, coin):
        output, record = original(state, seam, coin)
        return output, (1-coin, record[1])
    monkeypatch.setattr(q.native, 'native', corrupt)
    with pytest.raises(ArithmeticError, match='coin/current'):
        q.native_instrument('C')


@pytest.mark.parametrize('law', [None, 'c', 'X', 1, True])
def test_invalid_candidate_rejected(law):
    with pytest.raises(ValueError):
        q.native_instrument(law)


@pytest.mark.parametrize('count', [True, -1, 65, 1., None])
def test_invalid_count_rejected(count):
    with pytest.raises(ValueError):
        q.program('C', count)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1., 1e-3])
def test_numerical_failure_cannot_be_reported_as_pass(value):
    with pytest.raises(ArithmeticError):
        q.bounded(value, 'mutant')


def test_scope_does_not_promote_operator_log_to_physical_unitary_or_empirical_agreement(result):
    scope = result['scope']
    assert scope['native_channel_semantics_executed']
    assert scope['spectral_log_is_unique_on_full_binary_operator_space']
    for key in ('native_pulse_words_emitted', 'log_is_a_carrier_ket_Hamiltonian_or_Lindblad_generator',
                'fractional_transfer_powers_are_physical_channels',
                'all_quantum_history_reflection_axioms_proved',
                'physical_clock_or_continuum_QFT_derived',
                'empirical_dataset_used_or_empirical_agreement_claimed'):
        assert not scope[key]

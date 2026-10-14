"""Independent checks on native preparation, rotated instruments and records."""
from copy import deepcopy
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.observer.source_scalar import compiler as c
from model.scalar_observations.experimental_readout import executable_source as s


def _independent(dimension, alpha, thermal, u, v, theta):
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    displacement = expm(alpha * (a.T - a))
    ratio = thermal / (1 + thermal)
    weights = (1 - ratio) * ratio ** np.arange(dimension)
    weights[-1] += ratio ** dimension
    rho = (displacement * weights) @ displacement.T
    characteristic = np.trace(rho @ expm((u + 1j*v)*a.T - (u - 1j*v)*a))
    mean = float((np.exp(1j * theta) * characteristic).real)
    return np.array([(1 + mean) / 2, (1 - mean) / 2])


@pytest.mark.parametrize('dimension,alpha,thermal', [(2, .7, 0.), (4, -.4, .2), (8, .9, .1)])
def test_actual_source_probabilities_match_independent_full_finite_displacement(dimension, alpha, thermal):
    settings = [{'u': u, 'v': v, 'theta': theta}
                for u, v in [(0., 0.), (.2, -.7), (-.5, .3)]
                for theta in [0., -math.pi/2]]
    result = s.source_run(dimension, alpha, thermal, settings)
    for record, setting in zip(result['records'], settings, strict=True):
        expected = _independent(dimension, alpha, thermal, **setting)
        assert record['probabilities'] == pytest.approx(expected, abs=2e-12)
        assert record['probability_mass'] == pytest.approx(1., abs=2e-12)
        assert record['mean'] == pytest.approx(expected[0] - expected[1], abs=2e-12)
        assert record['both_detector_outcomes_retained'] is True
        assert record['postselection'] is False
    assert result['analytic_characteristic_used_for_prediction'] is False
    assert result['prediction_origin'] == 'norms_of_source_executed_pointer_branches'
    assert result['occupation_error_certified_here'] is False
    assert result['physical_service_error_certified'] is False


def test_thermal_cap_keeps_all_mass_and_has_claimed_initial_distance():
    dimension, thermal = 8, .3
    weights, metadata = s.thermal_weights(dimension, thermal)
    ratio = thermal/(1+thermal)
    infinite = (1-ratio)*ratio**np.arange(100)
    padded = np.pad(weights, (0, 100-dimension))
    assert weights.sum() == pytest.approx(1., abs=1e-15)
    assert abs(padded-infinite).sum()/2 == pytest.approx(ratio**dimension)
    assert weights[-1] == pytest.approx(ratio**(dimension-1))
    assert metadata['state_or_outcome_renormalization'] is False
    assert metadata['no_preparation_branch_discarded'] is True


def test_native_fock_preparation_retains_full_scalar_phase():
    prepared, records = s._fock_preparations(8)
    assert prepared == pytest.approx(np.eye(8), abs=1e-14)
    for n, record in enumerate(records):
        assert record['occupation'] == n
        assert record['resources']['native_pulses'] == 2*n.bit_count()
        assert record['resources']['identity_phase_pulses'] == n.bit_count()
        assert record['full_phase_target_residual'] < 1e-14


def test_branch_weights_are_retained_and_resources_include_every_preparation():
    result = s.source_run(4, .2, .5, [{'u': .1, 'v': .7, 'theta': 0.}])
    record = result['records'][0]
    weights = result['preparation']['thermal']['classical_branch_probabilities']
    by_branch = np.array(record['probabilities_by_initial_occupation'])
    assert np.array(weights) @ by_branch == pytest.approx(record['probabilities'])
    assert by_branch.sum(axis=1) == pytest.approx(np.ones(4), abs=2e-12)
    counts = record['native_pulses_by_initial_occupation']
    assert counts[1]-counts[0] == counts[2]-counts[0] == 2
    assert counts[3]-counts[0] == 4
    assert record['resources']['native_pulses']['maximum_per_shot'] == max(counts)
    assert record['resources']['native_pulses']['classical_weighted_mean_per_shot'] == pytest.approx(np.dot(weights, counts))
    assert record['resources']['proper_code_service_events']['maximum_per_shot'] == 2*max(counts)
    assert len(record['chronological_word_hashes_after_fock_preparation']) == 7
    assert all(len(digest) == 64 for digest in record['chronological_word_hashes_after_fock_preparation'])


def test_source_pulses_actually_control_predictions(monkeypatch):
    s._displacement.cache_clear()
    original = c.execute
    executed = []
    def instrumented(word, qubits, state=None):
        executed.append((len(word), qubits, state is not None))
        return original(word, qubits, state)
    monkeypatch.setattr(c, 'execute', instrumented)
    result = s.source_run(2, .11234, .07, [{'u': .1, 'v': .2, 'theta': 0.}])
    assert len(executed) >= 5
    assert any(count > 1 and bits == 1 for count, bits, _ in executed)
    assert result['preparation']['displacement']['matrix_origin'] == 'contracted_native_source_pulses'


def test_failed_contraction_cannot_be_replaced_by_target_matrix(monkeypatch):
    monkeypatch.setattr(c, 'execute', lambda word, qubits: np.zeros((2**qubits, 2**qubits), complex))
    with pytest.raises(ArithmeticError, match='native contraction differs|fixed rounding guard'):
        s._contract([c.pulse('X', .2, [0])], 1, np.eye(2), name='broken')


def test_returned_preparation_metadata_cannot_mutate_cached_operations():
    _, _, metadata = s.preparation(4, .123, .05)
    before = deepcopy(metadata)
    metadata['displacement']['resources']['native_pulses'] = -1
    metadata['fock_branches'][0]['occupation'] = -1
    _, _, again = s.preparation(4, .123, .05)
    assert again == before


@pytest.mark.parametrize('dimension', [True, 1, 3, 128, 4.0, -2])
def test_invalid_dimensions_rejected(dimension):
    with pytest.raises(ValueError):
        s.source_run(dimension, .1, .1, [{'u': 0, 'v': 0, 'theta': 0}])


@pytest.mark.parametrize('alpha,thermal', [(True, 0), (float('inf'), 0), (0, -.1), (0, True), (0, float('nan'))])
def test_invalid_preparation_parameters_rejected(alpha, thermal):
    with pytest.raises(ValueError):
        s.preparation(2, alpha, thermal)


@pytest.mark.parametrize('settings', [[], None, [{'u': 1, 'v': 0}],
    [{'u': 1, 'v': 0, 'theta': True}], [{'u': 1, 'v': float('inf'), 'theta': 0}],
    [{'u': 1, 'v': 0, 'theta': 0, 'observed_mean': .5}]])
def test_settings_accept_no_observed_values_or_invalid_coordinates(settings):
    with pytest.raises(ValueError):
        s.source_run(2, 0, 0, settings)


@pytest.mark.parametrize('dimension', [2, 4, 8, 16])
def test_canonical_gray_action_and_counts_match_existing_compiler_domain(dimension):
    from . import canonical_source
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    basis = expm(.23 * (ladder.T - ladder))
    canonical = canonical_source.displacement_word(dimension, .23)
    original = c.compile_real_orthogonal(basis)
    canonical_resources, original_resources = c.resources(canonical), c.resources(original)
    for key in s.COUNT_KEYS:
        assert canonical_resources[key] == original_resources[key]
    bits = dimension.bit_length() - 1
    actual = c.execute(canonical, bits)
    assert actual == pytest.approx(c.execute(original, bits), abs=2e-12)
    assert actual == pytest.approx(basis, abs=2e-12)


def test_canonical_64_level_displacement_passes_unchanged_source_guards():
    """Windows previously refused this actual million-pulse preparation."""
    from .data import load_observations
    from .executable_calibration import calibrate_state
    alpha = calibrate_state(load_observations()['calibration'], extent=4)['parameters']['alpha_real']
    s._displacement.cache_clear()
    actual, metadata = s._displacement(64, alpha)
    ladder = np.diag(np.sqrt(np.arange(1, 64)), 1)
    expected = expm(alpha*(ladder.T-ladder))
    assert np.linalg.norm(actual-expected, 2) < 2e-10
    assert metadata['resources']['native_pulses'] == 1162054
    assert metadata['numerical_stabilization']['maximum_raw_unitarity_defect'] == 1e-10
    assert metadata['numerical_stabilization']['maximum_projection_change'] == 5.1e-11
    ledger = metadata['blockwise_contraction']
    assert ledger['block_count'] == math.ceil(1162054/s.MAXIMUM_CONTRACTION_BLOCK_PULSES)
    assert sum(ledger['block_pulse_counts']) == 1162054
    assert ledger['maximum_block_raw_unitarity_defect'] <= 1e-10
    assert ledger['maximum_block_projection_change'] <= 5.1e-11
    assert ledger['maximum_block_restored_unitarity_defect'] <= 2e-13
    assert metadata['corrected_full_phase_target_residual'] < 2e-10


def test_long_word_retains_every_pulse_phase_resource_and_raw_comparison(monkeypatch):
    word = [event for _ in range(2050)
            for event in (c.pulse('X', .2, [0]), c.pulse('Z', -.11, [0]))]
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    target = np.linalg.matrix_power(expm(-.11j*z) @ expm(.2j*x), 2050)
    monolithic = c.execute(word, 1)
    original, executed = c.execute, []

    def traced(block, qubits):
        executed.extend(deepcopy(block))
        return original(block, qubits)

    monkeypatch.setattr(c, 'execute', traced)
    result, metadata = s._contract(word, 1, target, name='long_native_control')
    assert executed == word
    assert result == pytest.approx(monolithic, abs=3e-12)
    assert metadata['word_sha256'] == s._word_hash(word)
    assert metadata['resources'] == c.resources(word)
    assert metadata['full_phase_target_residual'] < 2e-10
    ledger = metadata['blockwise_contraction']
    assert ledger['block_pulse_counts'] == [4096, 4]
    assert not ledger['source_word_or_resource_counts_changed']
    assert not ledger['outward_roundoff_or_physical_error_certificate']
    assert ledger['complete_corrected_minus_raw_norm'] <= (
        ledger['algebraic_telescoping_expression_evaluated']+3e-14)


@pytest.mark.parametrize('fault', ['wrong_unitary', 'global_phase', 'nonunitary'])
def test_fixed_block_execution_retains_hostile_full_matrix_guards(monkeypatch, fault):
    word = [c.pulse('I', 0., [0])]*(s.MAXIMUM_CONTRACTION_BLOCK_PULSES+1)
    wrong = {'wrong_unitary': np.array([[0., 1.], [1., 0.]]),
             'global_phase': 1j*np.eye(2),
             'nonunitary': (1+1e-7)*np.eye(2)}[fault]
    calls = 0

    def corrupt(block, qubits):
        nonlocal calls
        calls += 1
        return wrong.copy() if calls == 1 else np.eye(2, dtype=complex)

    monkeypatch.setattr(c, 'execute', corrupt)
    with pytest.raises(ArithmeticError, match='independent target|fixed rounding guard'):
        s._contract(word, 1, np.eye(2), name='hostile')


def test_corrected_full_phase_gate_cannot_be_replaced_by_unitarity(monkeypatch):
    original = c.stabilize_executed_unitary

    def wrong_projection(raw):
        matrix, metadata = original(raw)
        return 1j*matrix, metadata

    monkeypatch.setattr(c, 'stabilize_executed_unitary', wrong_projection)
    with pytest.raises(ArithmeticError, match='corrected native contraction differs'):
        s._contract([c.pulse('I', 0., [0])], 1, np.eye(2), name='bad_projection')


def test_block_correction_account_uses_raw_prefix_norms(monkeypatch):
    epsilon = 2e-11
    block = np.diag([1+epsilon, 1-epsilon]).astype(complex)
    monkeypatch.setattr(c, 'execute', lambda word, qubits: block.copy())
    word = [c.pulse('I', 0., [0])]*(3*s.MAXIMUM_CONTRACTION_BLOCK_PULSES)
    matrix, metadata = s._contract(word, 1, np.eye(2), name='known_radial_drift')
    ledger = metadata['blockwise_contraction']
    expected = np.linalg.norm(np.linalg.matrix_power(block, 3)-np.eye(2), 2)
    assert matrix == pytest.approx(np.eye(2), abs=1e-15)
    assert ledger['complete_corrected_minus_raw_norm'] == pytest.approx(expected, abs=1e-15)
    assert ledger['algebraic_telescoping_expression_evaluated'] == pytest.approx(expected, abs=2e-15)
    assert ledger['sum_block_projection_changes'] > 5e-11
    assert not ledger['floating_matrix_product_roundoff_in_telescoping_expression']


def test_small_individual_corrections_cannot_evade_complete_raw_target_gate(monkeypatch):
    monkeypatch.setattr(c, 'execute', lambda word, qubits: (1+4e-11)*np.eye(2, dtype=complex))
    word = [c.pulse('I', 0., [0])]*(6*s.MAXIMUM_CONTRACTION_BLOCK_PULSES)
    with pytest.raises(ArithmeticError, match='native contraction differs from its independent target'):
        s._contract(word, 1, np.eye(2), name='accumulated_fault')

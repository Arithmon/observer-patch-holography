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
    with pytest.raises(ArithmeticError, match='native contraction differs'):
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
def test_extended_gray_algorithm_matches_existing_words_on_original_domain(dimension):
    rng = np.random.default_rng(184 + dimension)
    basis, _ = np.linalg.qr(rng.normal(size=(dimension, dimension)))
    extended = s._extended_gray_word(basis)
    original = c.compile_real_orthogonal(basis)
    assert extended == original
    # Independently contract smaller instances, including every scalar sign.
    if dimension <= 8:
        assert c.execute(extended, dimension.bit_length()-1) == pytest.approx(basis, abs=2e-12)

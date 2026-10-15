"""Independent circuit, hidden-label and complete-outcome source checks."""
from itertools import product
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.interacting_reconstruction import source
from model.relativistic_scalar.lattice import ScalarParameters


def oscillator(d, omega):
    # An independent assembly retains all intermediate ladder states before
    # compression; replacing q^4 by the fourth power of truncated q fails.
    a = np.diag(np.sqrt(np.arange(1, d+6, dtype=float)), 1)
    q = (a+a.T)/math.sqrt(2*omega)
    p = 1j*math.sqrt(omega/2)*(a.T-a)
    return {'q': q[:d, :d], 'q2': (q@q)[:d, :d],
            'p2': (p@p)[:d, :d], 'q4': np.linalg.matrix_power(q, 4)[:d, :d]}


def pointer_laws(state, initial, other):
    probabilities = []
    for sign in (-1, 1):
        vector = np.zeros(len(state), dtype=complex)
        vector[initial], vector[other] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
        yes = abs(np.vdot(vector, state))**2
        probabilities.append([np.vdot(state, state).real-yes, yes])
    return np.array(probabilities)


@pytest.mark.parametrize('d,occupation', [(2, 0), (4, 1), (4, 3), (8, 6)])
def test_source_preparation_is_complete_permutation_with_no_column_substitution(d, occupation):
    matrix, metadata = source.preparation(d, occupation)
    expected = np.eye(d)
    expected[:, [0, occupation]] = expected[:, [occupation, 0]]
    np.testing.assert_allclose(matrix, expected, atol=3e-12, rtol=0)
    assert metadata['resources']['native_pulses'] > 0
    assert metadata['full_phase_target_residual'] < 2e-10


@pytest.mark.parametrize('initial,other', [(0, 2), (1, 3), (0, 3), (3, 0)])
@pytest.mark.parametrize('sign', [-1, 1])
def test_pointer_matches_full_projector_and_keeps_disturbed_branches(initial, other, sign):
    kraus, metadata = source.instrument(4, initial, other, sign)
    vector = np.zeros(4, dtype=complex)
    vector[initial], vector[other] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
    projector = np.outer(vector, vector.conj())
    np.testing.assert_allclose(kraus[0], np.eye(4)-projector, atol=3e-12, rtol=0)
    np.testing.assert_allclose(kraus[1], projector, atol=3e-12, rtol=0)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(4), atol=3e-13, rtol=0)
    # Keep a reference system and an unnormalized input to disallow replacing
    # the branches with postselected probabilities or normalized state vectors.
    rng = np.random.default_rng(28)
    reference = .2*(rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3)))
    for actual, expected in zip(kraus, (np.eye(4)-projector, projector), strict=True):
        np.testing.assert_allclose(actual@reference, expected@reference, atol=3e-12, rtol=0)
    assert metadata['completeness_residual'] < 3e-13


@pytest.mark.parametrize('transition,duration,parameters', [
    (0, 1/256, ScalarParameters()),
    (1, .037, ScalarParameters()),
    (0, .031, ScalarParameters(mass=1.3, speed=.9, spacing=.8, coupling=.7, omega=1.7)),
    (1, .021, ScalarParameters(mass=.9, speed=1.1, coupling=0., omega=2.2)),
])
def test_onsite_record_matches_independent_compressed_quartic_circuit(transition, duration, parameters):
    op = oscillator(4, parameters.omega)
    h = op['p2']/2+(parameters.mass**2+3*parameters.kappa)*op['q2']/2+parameters.g*op['q4']
    state = expm(-1j*duration*h)[:, transition]
    actual, metadata = source.onsite_probe(5, transition, duration, parameters)
    np.testing.assert_allclose(actual, pointer_laws(state, transition, transition+2), atol=3e-12, rtol=0)
    np.testing.assert_allclose(actual.sum(axis=1), 1, atol=3e-13, rtol=0)
    assert metadata['fresh_trials'] == 2
    assert metadata['resources']['proper_code_service_events'] == 2*metadata['resources']['native_pulses']


def test_bonds_reveal_hidden_graph_by_detector_law_with_equal_resource_schedules():
    parameters = ScalarParameters(speed=.83, omega=1.6)
    positions = np.array(list(product((0, 1), repeat=3)))[list(source.PERMUTATION)]
    adjacency = np.sum(abs(positions[:, None]-positions[None, :]), axis=-1) == 1
    q = oscillator(2, parameters.omega)['q']
    counts = []
    for pair in source.PAIRS:
        h = -parameters.kappa*np.kron(q, q)*adjacency[pair]
        expected = pointer_laws(expm(-.13j*h)[:, 0], 0, 3)
        actual, metadata = source.bond_probe(pair, .13, parameters)
        np.testing.assert_allclose(actual, expected, atol=3e-12, rtol=0)
        counts.append(metadata['resources'])
    assert all(count == counts[0] for count in counts)


def test_calibration_retains_all_outcomes_and_does_not_export_action_or_clock():
    records = source.calibration_records(1/256)
    assert records['onsite_probabilities'].shape == (8, 2, 2, 2)
    assert records['bond_probabilities'].shape == (28, 2, 2)
    assert records['phase_signs'] == (-1, 1)
    assert records['pairs'] == source.PAIRS
    assert records['fresh_trials'] == 88
    assert not {'duration', 'parameters', 'hamiltonian', 'adjacency', 'mass', 'coupling'} & records.keys()
    for key in ('onsite_probabilities', 'bond_probabilities'):
        assert np.isfinite(records[key]).all()
        assert records[key].min() > 0
        np.testing.assert_allclose(records[key].sum(axis=-1), 1, atol=3e-13, rtol=0)
    assert records['maximum_local_target_residual'] < 2e-10


def test_unknown_label_permutation_changes_records_covariantly():
    original = source.calibration_records(.017)
    shuffle = (2, 5, 1, 7, 0, 6, 4, 3)
    permutation = tuple(source.PERMUTATION[index] for index in shuffle)
    changed = source.calibration_records(.017, permutation=permutation)
    table = {pair: index for index, pair in enumerate(source.PAIRS)}
    for index, pair in enumerate(source.PAIRS):
        old_pair = tuple(sorted(shuffle[x] for x in pair))
        np.testing.assert_allclose(changed['bond_probabilities'][index],
                                   original['bond_probabilities'][table[old_pair]], atol=1e-14, rtol=0)
    assert np.max(abs(changed['bond_probabilities']-original['bond_probabilities'])) > .001


def test_audit_accessor_is_not_a_probability_oracle(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('audit target leaked into source record production')

    monkeypatch.setattr(source, 'supplied_target', forbidden)
    records = source.calibration_records(.019)
    assert records['fresh_trials'] == 88


def test_corrupted_source_execution_fails_before_any_target_replacement(monkeypatch):
    original = source.c.execute

    def corrupted(*args, **kwargs):
        return np.exp(.001j)*original(*args, **kwargs)

    source._preparation.cache_clear()
    monkeypatch.setattr(source.c, 'execute', corrupted)
    try:
        with pytest.raises(ArithmeticError, match='full-phase'):
            source.preparation(4, 1)
    finally:
        source._preparation.cache_clear()


def test_preparation_and_instrument_caches_cannot_be_mutated_by_callers():
    matrix, metadata = source.preparation(4, 1)
    matrix[:] = 0
    metadata['resources']['native_pulses'] = -1
    repeated, repeated_metadata = source.preparation(4, 1)
    assert abs(repeated[1, 0]) == pytest.approx(1., abs=1e-13)
    assert repeated_metadata['resources']['native_pulses'] > 0
    kraus, _ = source.instrument(4, 1, 3, -1)
    kraus[1][:] = 0
    again, _ = source.instrument(4, 1, 3, -1)
    assert np.trace(again[1]).real == pytest.approx(1., abs=1e-13)


@pytest.mark.parametrize('call', [
    lambda: source.preparation(True),
    lambda: source.preparation(4, -1),
    lambda: source.instrument(4, 1, 1, 1),
    lambda: source.instrument(4, 0, 2, 0),
    lambda: source.instrument(4, 0, 2, True),
    lambda: source.onsite_probe(8, 0, .1),
    lambda: source.onsite_probe(0, 2, .1),
    lambda: source.onsite_probe(0, 0, -.1),
    lambda: source.bond_probe((0, 0), .1),
    lambda: source.bond_probe((0,), .1),
    lambda: source.calibration_records(float('nan')),
    lambda: source.calibration_records(True),
    lambda: source.calibration_records(.1, permutation=(0,)*8),
    lambda: source.calibration_records(.1, parameters={}),
])
def test_invalid_protocol_requests_are_rejected(call):
    with pytest.raises(ValueError):
        call()

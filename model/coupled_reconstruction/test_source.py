"""Independent pointer, tensor-trace and fully coupled source checks."""
from itertools import product
import math

import numpy as np
import pytest

from model.coupled_reconstruction import source
from model.interacting_reconstruction import source as isolated
from model.relativistic_scalar.lattice import ScalarLattice, ScalarParameters


def independent_density(state, dimensions, sites):
    spectators = tuple(i for i in range(len(dimensions)) if i not in sites)
    strides = [math.prod(dimensions[:i]) for i in range(len(dimensions))]
    local_configurations = list(product(*(range(dimensions[i]) for i in reversed(sites))))
    result = np.zeros((len(local_configurations),)*2, dtype=complex)
    for row, reverse_left in enumerate(local_configurations):
        left = tuple(reversed(reverse_left))
        for column, reverse_right in enumerate(local_configurations):
            right = tuple(reversed(reverse_right))
            for environment in product(*(range(dimensions[i]) for i in spectators)):
                offset = sum(n*strides[i] for i, n in zip(spectators, environment, strict=True))
                a = offset+sum(n*strides[i] for i, n in zip(sites, left, strict=True))
                b = offset+sum(n*strides[i] for i, n in zip(sites, right, strict=True))
                result[row, column] += state[a]*state[b].conjugate()
    return result


def independent_pointer_law(density, initial, other):
    result = []
    for sign in (-1, 1):
        vector = np.zeros(len(density), complex)
        vector[initial], vector[other] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
        yes = float(np.vdot(vector, density@vector).real)
        result.append([float(np.trace(density).real)-yes, yes])
    return np.array(result)


@pytest.mark.parametrize('sites', [(0,), (2,), (0, 2), (2, 0), (2, 1)])
def test_spectator_trace_preserves_order_coherences_and_unnormalized_input(sites):
    dimensions = (2, 3, 4)
    rng = np.random.default_rng(195)
    state = .13*(rng.normal(size=24)+1j*rng.normal(size=24))
    actual = source.reduced_density(state, dimensions, sites)
    expected = independent_density(state, dimensions, sites)
    np.testing.assert_allclose(actual, expected, atol=3e-15, rtol=0)
    assert np.trace(actual).real == pytest.approx(np.vdot(state, state).real, abs=3e-15)


@pytest.mark.parametrize('dimension,initial,other', [(4, 0, 2), (4, 1, 3), (16, 0, 5)])
@pytest.mark.parametrize('sign', [-1, 1])
def test_pointer_complete_columns_include_pair_spectator_occupations(dimension, initial, other, sign):
    kraus, metadata = source.instrument(dimension, initial, other, sign)
    vector = np.zeros(dimension, complex)
    vector[initial], vector[other] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
    yes = np.outer(vector, vector.conj())
    for actual, expected in zip(kraus, (np.eye(dimension)-yes, yes), strict=True):
        np.testing.assert_allclose(actual, expected, atol=4e-12, rtol=0)
    np.testing.assert_allclose(sum(k.conj().T@k for k in kraus), np.eye(dimension), atol=2e-14, rtol=0)
    mixed = np.eye(dimension)/dimension
    probabilities = [np.trace(k@mixed@k.conj().T).real for k in kraus]
    np.testing.assert_allclose(probabilities, [1-1/dimension, 1/dimension], atol=3e-14, rtol=0)
    assert metadata['full_phase_target_residual'] < 2e-10
    assert metadata['completeness_residual'] < 2e-14
    assert metadata['resources']['proper_code_service_events'] == 2*metadata['resources']['native_pulses']


@pytest.fixture(scope='module')
def small_time_records():
    return source.calibration_records(1/8192)


def test_full_field_record_agrees_with_independent_exact_hamiltonian_at_short_time(small_time_records):
    lattice = ScalarLattice((4,)*8)
    prepared = np.column_stack([lattice.vacuum_state()]+[
        lattice.basis_state(tuple(int(i == physical) for i in range(8))) for physical in source.PERMUTATION])
    states = lattice.evolve(prepared, 1/8192)
    for label, physical in enumerate(source.PERMUTATION):
        for transition in (0, 1):
            column = 0 if transition == 0 else label+1
            density = source.reduced_density(states[:, column], (4,)*8, [physical])
            expected = independent_pointer_law(density, transition, transition+2)
            np.testing.assert_allclose(small_time_records['onsite_probabilities'][label, transition],
                                       expected, atol=2e-11, rtol=0)
    for index, (left, right) in enumerate(source.PAIRS):
        density = source.reduced_density(states[:, 0], (4,)*8,
            [source.PERMUTATION[left], source.PERMUTATION[right]])
        expected = independent_pointer_law(density, 0, 5)
        np.testing.assert_allclose(small_time_records['bond_probabilities'][index], expected, atol=2e-11, rtol=0)


def test_record_protocol_and_fresh_trial_costs_include_whole_field(small_time_records):
    records = small_time_records
    assert records['onsite_probabilities'].shape == (8, 2, 2, 2)
    assert records['bond_probabilities'].shape == (28, 2, 2)
    assert records['fresh_trials'] == 88
    assert records['independent_preparation_columns_contracted'] == 9
    assert records['local_generators_per_trial'] == 28
    assert records['occupation_dimensions'] == (4,)*8
    assert records['fixed_symmetric_steps'] == 1
    assert records['complete_field_schedule_in_every_trial']
    assert not records['independently_isolated_action_services']
    assert not {'duration', 'parameters', 'adjacency', 'hamiltonian', 'mass'} & records.keys()
    assert len(records['pointer_word_sha256']) == 6
    instruments = [source.instrument(4, n, n+2, sign)[1] for n in (0, 1) for sign in (-1, 1)]
    pair_instruments = [source.instrument(16, 0, 5, sign)[1] for sign in (-1, 1)]
    preparations = [isolated.preparation(4, n)[1] for n in (0, 1)]
    assert preparations[0]['resources'] == preparations[1]['resources']
    for key in source.COUNT_KEYS:
        expected = (88*(records['evolution_resources_per_trial'][key]+8*preparations[0]['resources'][key])
                    +8*sum(info['resources'][key] for info in instruments)
                    +28*sum(info['resources'][key] for info in pair_instruments))
        assert records['resources'][key] == expected
    for key in ('onsite_probabilities', 'bond_probabilities'):
        laws = records[key]
        assert laws.min() > 0
        np.testing.assert_allclose(laws.sum(axis=-1), 1., atol=3e-13, rtol=0)


def test_finite_time_records_reveal_actual_coupling_beyond_isolated_generator():
    # This larger-duration diagnostic is numerical. It is not used for the
    # small-duration coefficient certificate or for fitting any parameters.
    duration = 1/8
    coupled = source.calibration_records(duration)
    for transition in (0, 1):
        old, _ = isolated.onsite_probe(0, transition, duration)
        assert np.max(abs(coupled['onsite_probabilities'][0, transition]-old)) > 1e-4
    adjacency = isolated.supplied_target()['adjacency']
    differences = []
    phase_asymmetries = []
    for index, pair in enumerate(source.PAIRS):
        if not adjacency[pair]:
            old, _ = isolated.bond_probe(pair, duration)
            actual = coupled['bond_probabilities'][index]
            differences.append(np.max(abs(actual-old)))
            phase_asymmetries.append(abs(actual[0, 1]-actual[1, 1]))
    assert max(differences) > 1e-4
    assert max(phase_asymmetries) > 1e-5


def test_label_permutation_changes_detector_table_covariantly(small_time_records):
    shuffle = (2, 5, 1, 7, 0, 6, 4, 3)
    changed = source.calibration_records(1/8192,
        permutation=tuple(source.PERMUTATION[i] for i in shuffle))
    np.testing.assert_allclose(changed['onsite_probabilities'],
        small_time_records['onsite_probabilities'][list(shuffle)], atol=3e-13, rtol=0)
    lookup = {pair: index for index, pair in enumerate(source.PAIRS)}
    for index, pair in enumerate(source.PAIRS):
        old_index = lookup[tuple(sorted(shuffle[i] for i in pair))]
        np.testing.assert_allclose(changed['bond_probabilities'][index],
            small_time_records['bond_probabilities'][old_index], atol=3e-13, rtol=0)
    assert np.max(abs(changed['bond_probabilities']-small_time_records['bond_probabilities'])) > 1e-5


def test_no_isolated_probe_or_audit_accessor_supplies_calibration(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('a hidden isolated probe or audit target was used')

    monkeypatch.setattr(isolated, 'onsite_probe', forbidden)
    monkeypatch.setattr(isolated, 'bond_probe', forbidden)
    monkeypatch.setattr(isolated, 'supplied_target', forbidden)
    records = source.calibration_records(1/16384, ScalarParameters(coupling=.7))
    assert records['fresh_trials'] == 88


def test_changed_executed_source_cannot_be_replaced_by_analytic_pointer(monkeypatch):
    original = source.c.execute

    def corrupted(*args, **kwargs):
        return np.exp(.001j)*original(*args, **kwargs)

    source._instrument.cache_clear()
    monkeypatch.setattr(source.c, 'execute', corrupted)
    try:
        with pytest.raises(ArithmeticError, match='full-phase'):
            source.instrument(16, 0, 5, 1)
    finally:
        source._instrument.cache_clear()


@pytest.mark.parametrize('call', [
    lambda: source.instrument(8, 0, 2, 1),
    lambda: source.instrument(16, 0, 0, 1),
    lambda: source.instrument(16, 0, 5, True),
    lambda: source.instrument(16, 0, 16, 1),
    lambda: source.calibration_records(-.1),
    lambda: source.calibration_records(np.nan),
    lambda: source.calibration_records(True),
    lambda: source.calibration_records(.1, steps=0),
    lambda: source.calibration_records(.1, steps=True),
    lambda: source.calibration_records(.1, parameters={}),
    lambda: source.calibration_records(.1, permutation=(0,)*8),
    lambda: source.reduced_density(np.ones(4), (2, 2), [0, 0]),
    lambda: source.reduced_density(np.ones(4), (2, 2), [2]),
    lambda: source.reduced_density(np.ones(3), (2, 2), [0]),
])
def test_invalid_protocols_rejected(call):
    with pytest.raises(ValueError):
        call()

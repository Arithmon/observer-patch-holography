from fractions import Fraction
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.anharmonic_observations.reference import readout
from . import intervals as iv
from .experiment import detector, run, thermal_state
from .solver import cell, evolve, generator


def test_interval_operations_contain_exact_rational_arithmetic():
    rng = np.random.default_rng(615)
    for x, y in rng.normal(size=(200, 2)):
        for answer, interval in ((Fraction(x)+Fraction(y), iv.add(iv.exact(x), iv.exact(y))),
                                 (Fraction(x)*Fraction(y), iv.mul(iv.exact(x), iv.exact(y))),
                                 (Fraction(x)-Fraction(y), iv.sub(iv.exact(x), iv.exact(y)))):
            assert Fraction(float(interval[0])) <= answer <= Fraction(float(interval[1]))
    tiny = np.nextafter(0., 1.)
    lo, hi = iv.mul(iv.exact(tiny), iv.exact(.5))
    assert Fraction(float(lo)) <= Fraction(tiny)/2 <= Fraction(float(hi))


def test_pairwise_sum_bounds_exact_long_sum_and_rejects_nonfinite():
    values = np.geomspace(1e-300, 1e100, 1001)
    assert Fraction(iv.sum_up(values)) >= sum(map(Fraction, values))
    with pytest.raises(ArithmeticError):
        iv.sum_up([float('inf')])
    with np.errstate(over='ignore'), pytest.raises(ArithmeticError):
        iv.sum_up([np.finfo(float).max, np.finfo(float).max])


@pytest.mark.parametrize('dimension,sectors', [(3, 2), (8, 2), (5, 4)])
def test_interval_sparse_generator_encloses_independent_dense_operator(dimension, sectors):
    rng = np.random.default_rng(dimension)
    raw = rng.normal(size=(dimension*sectors, 3))+1j*rng.normal(size=(dimension*sectors, 3))
    rho = raw@raw.conj().T
    rho = (rho+rho.conj().T)/2
    operation = generator(dimension, sectors=sectors, delta=1.3, force=-.7,
                          force_spin_flip=True, spin_z=.3, gamma=.1)
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    x = np.kron(np.eye(sectors//2), [[0., 1.], [1., 0.]])
    z = np.diag(np.tile([1., -1.], sectors//2))
    n = np.kron(np.eye(sectors), np.diag(np.arange(dimension)))
    h = 1.3*n-.7*np.kron(x, a+a.T)+.3*np.kron(z, np.eye(dimension))
    direct = -1j*(h@rho-rho@h)-.05*(n@n@rho+rho@n@n-2*n@rho@n)
    actual = operation.apply_float(rho)
    assert np.max(abs(direct-actual)) < 1e-12
    real, imag = operation.apply_interval(rho)
    # Dense BLAS itself is not interval arithmetic: allow its independently
    # measured rounding discrepancy when comparing the two implementations.
    assert np.all(direct.real >= real[0]-1e-12)
    assert np.all(direct.real <= real[1]+1e-12)
    assert np.all(direct.imag >= imag[0]-1e-12)
    assert np.all(direct.imag <= imag[1]+1e-12)


def test_full_defect_controls_actual_escape_beyond_cutoff():
    d, large = 3, 16
    rho = np.zeros((2*d, 2*d), complex); rho[d-1, d-1] = 1
    approximate, budget = evolve(rho, cells=4, force=.7, duration=.4)
    a = np.diag(np.sqrt(np.arange(1, large)), 1)
    u = expm(-.4j*.7*np.kron(np.eye(2), a+a.T))
    initial = np.zeros((2*large, 2*large), complex); initial[d-1, d-1] = 1
    reference = u@initial@u.conj().T
    embed = np.kron(np.eye(2), np.eye(large)[:, :d])
    error = np.linalg.norm(reference-embed@approximate@embed.T, 'nuc')/2
    assert error > .2
    assert error < iv.sum_up(list(budget.values()))
    assert budget['occupation_defect'] > budget['finite_defect']


def test_coarse_polynomial_is_charged_for_its_time_defect():
    rho = np.zeros((8, 8), complex); rho[0, 0] = 1
    approximate, budget = evolve(rho, cells=1, order=2, force=1., duration=.5)
    a = np.diag(np.sqrt(np.arange(1, 4)), 1)
    u = expm(-.5j*np.kron(np.eye(2), a+a.T))
    error = np.linalg.norm(u@rho@u.conj().T-approximate, 'nuc')/2
    assert error > .03
    assert error <= budget['finite_defect']+budget['endpoint_jump']


def test_preparation_tail_is_exact_capped_geometric_distance():
    d, thermal = 5, .2
    state, budget = thermal_state(d, thermal)
    ratio = Fraction(thermal)/(1+Fraction(thermal))
    assert Fraction(budget['thermal_tail']) >= ratio**d
    assert sum(budget.values()) < 2e-4
    assert np.trace(state).real == pytest.approx(1.)


def test_detector_includes_occupation_escape_and_has_correct_phase():
    state, _ = thermal_state(12, .13)
    state, _ = evolve(state, cells=8, force=-.9, force_imaginary=True)
    actual, budget = detector(state, .4)
    direct = readout(state, [{'u': 0., 'v': .4, 'theta': -math.pi/2}])['reads'][0]
    assert actual['selected_mean'] == pytest.approx(direct['selected_mean'], abs=3e-12)
    assert actual['selected_mean'] < -.5
    assert budget['occupation_defect'] > 0


def test_run_returns_a_nonvacuous_complete_certificate_without_native_or_data_access():
    result = run(dimension=32, case='b', steps=2)
    final = result['records'][-1]
    assert final['half_trace_error_upper'] < .02
    assert final['acceptance_positive_certified']
    assert final['drive_elapsed_seconds'] == pytest.approx(.000414)
    assert result['certificate']['preparation_and_detector_displacements_included']
    assert not result['certificate']['native_source_words_executed']
    assert not result['certificate']['empirical_outcomes_loaded']
    assert set(final['budget']) == {'finite_defect', 'occupation_defect', 'endpoint_jump',
                                   'thermal_tail', 'thermal_rounding', 'pointer_embedding_rounding'}


def test_nonhermitian_path_fails_closed():
    state = np.eye(4, dtype=complex); state[0, 1] = 1j
    with pytest.raises(ValueError):
        cell(state, generator(2))


@pytest.mark.parametrize('probe', [0., -1., float('inf'), float('nan'), True])
def test_detector_rejects_invalid_probe(probe):
    with pytest.raises(ValueError):
        detector(np.eye(8, dtype=complex)/8, probe)


def test_uncertified_acceptance_returns_no_conditional_radius(monkeypatch):
    from . import experiment
    monkeypatch.setattr(experiment, '_total', lambda budget: 1.)
    record = experiment.run(dimension=8, steps=0)['records'][0]
    assert not record['acceptance_positive_certified']
    assert record['selected_probability_error_upper'] is None
    assert record['selected_mean_error_upper'] is None
    assert record['estimator_error_upper'] is None


@pytest.mark.parametrize('consumer,representation', [
    (consumer, representation)
    for consumer in ('cell', 'evolve', 'detector')
    for representation in ('boolean_array', 'nested_boolean', 'masked_array', 'nested_masked')
    if consumer != 'detector' or representation in ('boolean_array', 'masked_array')
])
def test_original_state_missingness_and_booleans_are_not_erased(consumer, representation):
    state = np.zeros((4, 4), complex)
    state[0, 0] = 1
    if representation == 'boolean_array':
        supplied = state.astype(bool)
    elif representation == 'nested_boolean':
        supplied = state.tolist()
        supplied[0][0] = True
    elif representation == 'masked_array':
        supplied = np.ma.array(state, mask=np.eye(4, dtype=bool))
    else:
        supplied = [np.ma.array(state[0], mask=[True, False, False, False]),
                    *state[1:]]
    with pytest.raises(ValueError, match='Boolean|masked'):
        if consumer == 'cell':
            cell(supplied, generator(2))
        elif consumer == 'evolve':
            evolve(supplied, cells=1)
        else:
            detector(supplied, .4)


@pytest.mark.parametrize('consumer', ['cell', 'evolve', 'detector'])
def test_ordinary_complex_state_still_reaches_the_public_consumer(consumer):
    vector = np.array([.5, .5j, -.5, .5j])
    state = np.outer(vector, vector.conj())
    assert np.any(state.imag)
    if consumer == 'cell':
        actual, budget = cell(state, generator(2))
        np.testing.assert_array_equal(actual, state)
    elif consumer == 'evolve':
        actual, budget = evolve(state, cells=1)
        np.testing.assert_array_equal(actual, state)
    else:
        actual, budget = detector(state, .4)
        assert actual['acceptance'] > 0
        assert np.isfinite(actual['selected_mean'])
    assert all(np.isfinite(value) and value >= 0 for value in budget.values())

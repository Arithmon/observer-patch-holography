"""Independent local-channel, reference, exact-budget and custody controls."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import open_environment as e, protocol as p


def _program(d=2):
    return p.continuation_program(d, steps=2)


def _native_maps(d, site, channel):
    # Independent tensor construction, not the adapter's flattened indexing.
    return tuple(np.kron(np.kron(u, np.eye(d)), np.eye(4)) if site == 1
                 else np.kron(np.kron(np.eye(d), u), np.eye(4))
                 for u in channel.unitaries)


@pytest.mark.parametrize('d,site,reference', [(2, 0, 1), (2, 1, 3), (4, 0, 2), (4, 1, 2)])
def test_all_executed_kraus_branches_match_independent_tensor_reference_oracle(d, site, reference):
    program, noise = _program(d), e.PhaseNoise(2, site, F(1, 25), 8)
    random = np.random.default_rng(480+site+reference)
    vector = random.normal(size=(program.work_dimension*reference, 3))
    vector = vector+1j*random.normal(size=vector.shape)
    rho = vector @ vector.conj().T
    rho /= np.trace(rho)
    original = rho.copy()
    channel = e.source.channel(d, 1., float(noise.variance), noise.order)
    maps = _native_maps(d, site, channel)
    expected = sum(weight*np.kron(unitary, np.eye(reference)) @ rho
                   @ np.kron(unitary.conj().T, np.eye(reference))
                   for weight, unitary in zip(channel.weights, maps, strict=True))
    actual = e.apply(rho, program, noise)
    assert np.max(np.abs(actual-expected)) < 2e-15
    assert np.array_equal(rho, original)
    assert abs(np.trace(actual)-1) < 2e-15
    assert np.linalg.eigvalsh(actual).min() > -2e-15
    marginal = lambda matrix: np.einsum('aras->rs', matrix.reshape(
        program.work_dimension, reference, program.work_dimension, reference))
    assert np.max(np.abs(marginal(actual)-marginal(rho))) < 2e-15


@pytest.mark.parametrize('site', [0, 1])
def test_every_matrix_unit_has_correct_site_and_gauge_action(site):
    program, noise = _program(), e.PhaseNoise(0, site, F(1, 7), 4)
    size = program.work_dimension
    channel = e.source.channel(2, 1., float(noise.variance), noise.order)
    maps = _native_maps(2, site, channel)
    for a in range(size):
        for b in range(size):
            matrix = np.zeros((size, size), complex)
            matrix[a, b] = 1j
            expected = sum(weight*unitary @ matrix @ unitary.conj().T
                           for weight, unitary in zip(channel.weights, maps, strict=True))
            assert np.max(np.abs(e.apply(matrix, program, noise)-expected)) < 2e-15


def test_noise_preserves_fourth_occupation_moment_and_changes_coherence():
    d, program, noise = 4, _program(4), e.PhaseNoise(0, 1, F(1, 25), 8)
    numbers = np.diag(np.arange(d))
    total = np.kron(np.kron(numbers, np.eye(d))+np.kron(np.eye(d), numbers), np.eye(4))
    moment = np.linalg.matrix_power(total+2*np.eye(program.work_dimension), 4)
    channel = e.source.channel(d, 1., float(noise.variance), noise.order)
    adjoint = sum(weight*unitary.conj().T@moment@unitary
                  for weight, unitary in zip(channel.weights, _native_maps(d, 1, channel), strict=True))
    assert np.linalg.norm(adjoint-moment, 2) < 2e-12
    vector = np.zeros(program.work_dimension, complex)
    vector[0] = vector[4*d*(d-1)] = 1/math.sqrt(2)
    rho = np.outer(vector, vector.conj())
    actual = e.apply(rho, program, noise)
    assert np.trace(moment@actual) == pytest.approx(np.trace(moment@rho), abs=2e-12)
    assert abs(actual[0, 4*d*(d-1)]) < abs(rho[0, 4*d*(d-1)])


def test_fresh_draws_are_not_a_coherent_sum_or_one_reused_draw():
    program, noise = _program(4), e.PhaseNoise(0, 0, F(1, 25), 8)
    vector = np.zeros(program.work_dimension, complex)
    vector[0] = vector[12] = 1/math.sqrt(2)
    rho = np.outer(vector, vector.conj())
    channel = e.source.channel(4, 1., float(noise.variance), noise.order)
    maps = _native_maps(4, 0, channel)
    fresh = e.apply(e.apply(rho, program, noise), program, noise)
    reused = sum(weight*(unitary@unitary)@rho@(unitary.conj().T@unitary.conj().T)
                 for weight, unitary in zip(channel.weights, maps, strict=True))
    coherent = sum(math.sqrt(weight)*unitary for weight, unitary in zip(channel.weights, maps, strict=True))
    assert np.linalg.norm(fresh-reused, ord='nuc')/2 > .1
    assert abs(np.trace(coherent@rho@coherent.conj().T)-1) > 1
    assert abs(np.trace(fresh)-1) < 2e-15


@pytest.mark.parametrize('order', [1, 2, 8])
def test_gaussian_comparison_and_identity_radius_include_arbitrary_reference(order):
    program, noise = _program(4), e.PhaseNoise(0, 0, F(1, 25), order)
    vector = np.zeros(program.work_dimension*2, complex)
    vector[0] = vector[12*2+1] = 1/math.sqrt(2)
    rho = np.outer(vector, vector.conj())
    actual, gaussian = e.apply(rho, program, noise), e.apply(rho, program, noise, 'gaussian')
    expected = rho.copy()
    expected[0, 25] *= math.exp(-float(noise.variance)*9/2)
    expected[25, 0] = expected[0, 25].conjugate()
    assert np.max(np.abs(gaussian-expected)) < 2e-15
    metadata = e.event_metadata(program, noise)
    assert np.linalg.norm(actual-gaussian, ord='nuc')/2 <= float(F(
        metadata['ideal_finite_to_gaussian_half_diamond_upper']))+2e-14
    assert np.linalg.norm(actual-rho, ord='nuc')/2 <= float(F(
        metadata['ideal_finite_to_identity_half_diamond_upper']))+2e-14
    if order == 1:
        assert metadata['ideal_finite_to_identity_half_diamond_upper'] == '0'
        assert np.linalg.norm(actual-gaussian, ord='nuc')/2 > .08


def test_exact_fraction_conversion_and_underflow_are_charged():
    for variance in (F(1, 10), F(1, 10**400)):
        event = e.PhaseNoise(0, 0, variance, 4)
        result = e.event_metadata(_program(4), event)
        conversion = abs(variance-F(float(variance)))*F(9, 4)
        assert F(result['variance_conversion_half_diamond_upper']) == conversion > 0
        assert F(result['ideal_finite_to_gaussian_half_diamond_upper']) == min(F(1), conversion+F(
            result['ideal_quadrature_half_diamond_upper']))
        assert result['variance_exact'] == str(variance)
    assert F(result['ideal_finite_to_gaussian_half_diamond_upper']) > 0
    assert result['represented_variance_exact'] == '0'


def test_schedule_counts_each_application_and_detaches_every_field():
    program = _program(4)
    noise = e.PhaseNoise(0, 0, 0, 8)
    env = e.Environment((noise, e.PhaseNoise(2, 0, F(1, 25), 4),
                         e.PhaseNoise(2, 1, F(1, 25), 4)))
    checked = e.validated(program, env)
    object.__setattr__(noise, 'variance', F(1))
    assert env.events[0].variance == checked.events[0].variance == 0
    specification = env.specification()
    specification['events'][0]['variance'] = '99'
    assert env.specification()['events'][0]['variance'] == '0'
    result = e.schedule_budget(program, env)
    assert result['independent_environment_draws'] == 3
    assert result['resources_per_realization']['native_pulses'] == 9
    assert result['resources_for_all_branch_contractions']['native_pulses'] == 48
    for key in ('ideal_finite_to_gaussian_half_diamond_upper', 'ideal_finite_to_identity_half_diamond_upper'):
        assert F(result[key]) == min(F(1), sum(F(row[key]) for row in result['events']))
    result['events'][0]['source']['resources_per_realization']['native_pulses'] = 99
    assert e.event_metadata(program, env.events[0])['resources_per_realization']['native_pulses'] == 3
    assert not result['floating_execution_error_included']
    assert not result['physical_noise_parameters_or_service_errors_certified']
    assert e.schedule_budget(program, e.Environment())['independent_environment_draws'] == 0


@pytest.mark.parametrize('changes', [
    {'after': True}, {'after': -1}, {'after': .0}, {'site': True}, {'site': 2},
    {'variance': .1}, {'variance': True}, {'variance': -1}, {'variance': '1'},
    {'variance': F(10**400)}, {'order': True}, {'order': 0}, {'order': 65},
])
def test_malformed_noise_specification_refused(changes):
    arguments = dict(after=0, site=0, variance=F(1, 4), order=4)
    with pytest.raises(ValueError):
        e.PhaseNoise(**(arguments | changes))


def test_invalid_schedule_or_mutated_frozen_event_refused():
    a, b = e.PhaseNoise(0, 0, 1), e.PhaseNoise(0, 1, 1)
    for events in ([a], (a, a), (b, a), (object(),)):
        with pytest.raises(ValueError):
            e.Environment(events)
    with pytest.raises(ValueError, match='exceeds'):
        e.validated(_program(), e.Environment((e.PhaseNoise(99, 0, 1),)))
    env = e.Environment((a,))
    object.__setattr__(env.events[0], 'site', True)
    with pytest.raises(ValueError, match='site'):
        e.validated(_program(), env)


@pytest.mark.parametrize('mutation', ['weights', 'schur', 'unitary', 'all_identity', 'metadata', 'hash', 'cost', 'boolean_cost'])
def test_cached_source_mutation_cannot_silently_change_executed_channel(monkeypatch, mutation):
    program, noise = _program(4), e.PhaseNoise(0, 0, F(1, 25), 4)
    original = e.source.channel(4, 1., float(noise.variance), noise.order)
    weights, schur = original.weights.copy(), original.schur.copy()
    unitaries, metadata = tuple(u.copy() for u in original.unitaries), deepcopy(original.metadata)
    if mutation == 'weights':
        weights[0] += .001
        weights[1] -= .001
    elif mutation == 'schur':
        schur[0, 1] = 0
    elif mutation == 'unitary':
        unitaries[0][0, 0] *= -1
    elif mutation == 'all_identity':
        unitaries = tuple(np.eye(4) for _ in unitaries)
        schur[:] = 1
    elif mutation == 'metadata':
        metadata['duration_seconds'] *= 2
    elif mutation == 'hash':
        metadata['branches'][0]['source']['word_sha256'] = '0'*64
    elif mutation == 'cost':
        metadata['resources_per_realization']['native_pulses'] += 1
    else:
        metadata['resources_per_realization']['identity_phase_pulses'] = True
    poisoned = replace(original, weights=weights, schur=schur, unitaries=unitaries, metadata=metadata)
    monkeypatch.setattr(e.source, 'channel', lambda *args: poisoned)
    with pytest.raises(ValueError, match='source environment'):
        e.apply(np.eye(program.work_dimension)/program.work_dimension, program, noise)


@pytest.mark.parametrize('bad', [np.zeros((0, 0)), np.zeros((3, 3)), np.zeros(16),
                               np.full((16, 16), np.nan), np.eye(16, dtype=bool),
                               np.ma.array(np.eye(16), mask=False)])
def test_invalid_operators_refused_without_projection_or_renormalization(bad):
    with pytest.raises(ValueError):
        e.apply(bad, _program(), e.PhaseNoise(0, 0, 1))


@pytest.mark.parametrize('variance,order', [(F(1, 2**1074), 16), (F(10**300), 16)])
def test_oversized_uncapped_certificate_refused_before_source_execution(monkeypatch, variance, order):
    event = e.PhaseNoise(0, 0, variance, order)
    def forbidden(*args):
        raise AssertionError('source compilation must not precede certificate size validation')
    monkeypatch.setattr(e.source, 'channel', forbidden)
    with pytest.raises(ValueError, match='integer serialization budget'):
        e.event_metadata(_program(), event)
    with pytest.raises(ValueError, match='integer serialization budget'):
        e.apply(np.eye(16)/16, _program(), event)


def test_oversized_original_exact_variance_is_not_hidden_by_float_rounding():
    # This finite positive variance rounds to one, but its original exact
    # numerator and denominator cannot fit the declared receipt budget.
    variance = F(10**4100+1, 10**4100)
    assert float(variance) == 1.
    with pytest.raises(ValueError, match='integer serialization budget'):
        e.PhaseNoise(0, 0, variance, 1)


def test_serialization_refusal_preserves_both_valid_underflow_neighbors():
    program = _program()
    rho = np.ones((16, 16))/16
    for variance, order in ((F(1, 2**1074), 2), (F(1, 2**1075), 16)):
        event = e.PhaseNoise(0, 0, variance, order)
        metadata = e.event_metadata(program, event)
        actual = e.apply(rho, program, event)
        assert abs(np.trace(actual)-1) < 3e-15
        assert F(metadata['ideal_finite_to_gaussian_half_diamond_upper']) > 0
        assert metadata['maximum_exact_certificate_integer_digits'] <= 4000


def test_combined_radii_are_preflighted_before_any_source_execution(monkeypatch):
    def forbidden(*args):
        raise AssertionError('combined certificate validation must precede source compilation')
    monkeypatch.setattr(e.source, 'channel', forbidden)
    program = _program()
    # The original fraction and the dyadic quadrature radius separately fit;
    # adding the conversion charge gives a much larger common denominator.
    event = e.PhaseNoise(0, 0, F(10**3800, 10**3801+1), 64)
    with pytest.raises(ValueError, match='integer serialization budget'):
        e.event_metadata(program, event)
    # Individually admissible events can likewise exceed the budget only
    # after their distinct exact conversion denominators are combined.
    events = tuple(e.PhaseNoise(index//2, index % 2, F(1, 10**250+2*index+1), 2)
                   for index in range(20))
    for event in events:
        e.validated(program, e.Environment((event,)))
    with pytest.raises(ValueError, match='integer serialization budget'):
        e.schedule_budget(program, e.Environment(events))

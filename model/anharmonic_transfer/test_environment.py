"""Finite environments: native execution, reference systems and ideal bounds."""
from fractions import Fraction
import math

import numpy as np
import pytest

from model.anharmonic_transfer import environment as e
from model.scalar_observations.experimental_readout import canonical_source as canonical


@pytest.mark.parametrize('order', [1, 2, 4, 8, 16])
def test_normal_quadrature_has_positive_weights_and_correct_moments(order):
    context = canonical._context()
    nodes, weights = e.canonical_quadrature(order)
    nodes, weights = list(map(context.mpf, nodes)), list(map(context.mpf, weights))
    assert all(weight > 0 for weight in weights)
    for power in range(2*order):
        target = (context.zero if power % 2 else
                  context.factorial(power)/(2**(power//2)*context.factorial(power//2)))
        measured = context.fsum(weight*node**power for node, weight in zip(nodes, weights))
        assert abs(measured-target) < context.mpf('1e-55')*max(1, abs(target))


@pytest.mark.parametrize('dimension,order,rate,duration', [
    (2, 1, 2., .2), (4, 2, 3., .03), (8, 4, 2., .01),
    (16, 8, 4., .001), (64, 16, 18., 25e-6),
])
def test_high_precision_characteristic_errors_obey_ideal_bound(dimension, order, rate, duration):
    context = canonical._context()
    nodes, weights = e.canonical_quadrature(order)
    variance = canonical._parameter(context, rate, 'rate')*canonical._parameter(context, duration, 'duration')
    radius = Fraction(e.quadrature_bound(dimension, rate, duration, order)
                      ['exact_maximum_schur_entry_error_upper'])
    radius_mp = context.mpf(radius.numerator)/radius.denominator
    for frequency in range(dimension):
        finite = context.fsum(context.mpf(weight)*context.exp(
            -context.j*frequency*context.sqrt(variance)*context.mpf(node))
            for node, weight in zip(nodes, weights))
        target = context.exp(-variance*frequency**2/2)
        assert abs(finite-target) <= radius_mp+context.mpf('1e-65')


@pytest.mark.parametrize('dimension', [2, 4, 16, 64])
def test_every_branch_has_phase_exact_native_word_and_count(dimension):
    executed = e.channel(dimension, 18., 25e-6)
    for unitary, branch in zip(executed.unitaries, executed.metadata['branches']):
        target = np.diag(np.exp(-1j*branch['phase']*np.arange(dimension)))
        assert np.linalg.norm(unitary-target, 2) < 2e-12
        assert branch['source']['matrix_origin'] == 'contracted_native_source_pulses'
    expected = dimension.bit_length()
    assert executed.metadata['resources_per_realization']['native_pulses'] == expected
    assert executed.metadata['resources_for_all_branch_contractions']['native_pulses'] == 16*expected
    assert executed.metadata['all_phase_unitaries_source_executed']
    assert not executed.metadata['analytic_dephasing_schur_matrix_inserted']
    assert not executed.metadata['gaussian_noise_law_or_rate_derived_from_observer']


@pytest.mark.parametrize('reference', [1, 2, 3])
def test_schur_contraction_matches_all_kraus_branches_with_reference(reference):
    dimension = 4
    executed = e.channel(dimension, 2., .07, 8)
    random = np.random.default_rng(1729+reference)
    matrix = random.normal(size=(reference*dimension, 3))+1j*random.normal(size=(reference*dimension, 3))
    rho = matrix @ matrix.conj().T
    rho /= np.trace(rho)
    explicit = np.zeros_like(rho)
    completeness = np.zeros((dimension, dimension), complex)
    for weight, unitary in zip(executed.weights, executed.unitaries):
        full = np.kron(np.eye(reference), unitary)
        explicit += weight*full @ rho @ full.conj().T
        completeness += weight*unitary.conj().T @ unitary
    assert np.max(np.abs(executed.apply(rho)-explicit)) < 2e-15
    assert np.max(np.abs(completeness-np.eye(dimension))) < 3e-15
    assert abs(np.trace(explicit)-1) < 2e-15
    assert np.linalg.eigvalsh(explicit).min() > -2e-15


def test_phase_damping_preserves_entangled_reference_and_populations():
    dimension = 8
    executed = e.channel(dimension, 2., .01, 4)
    bell = np.zeros(2*dimension, complex)
    bell[0] = bell[-1] = 1/math.sqrt(2)
    rho = np.outer(bell, bell.conj())
    result = executed.apply(rho)
    target = rho.copy()
    target[0, -1] *= math.exp(-2*.01*(dimension-1)**2/2)
    target[-1, 0] = target[0, -1].conjugate()
    difference = np.linalg.norm(result-target, ord='nuc')/2
    assert 0 < difference <= executed.half_diamond_error_bound+2e-14
    assert np.max(np.abs(np.diag(result)-np.diag(rho))) < 2e-15
    assert abs(result[0, -1]) < abs(rho[0, -1])


@pytest.mark.parametrize('rate,duration', [(0., 1.), (18., 0.), (0., 0.)])
def test_zero_diffusion_preserves_input_and_declared_execution_schedule(rate, duration):
    executed = e.channel(4, rate, duration)
    rho = np.ones((8, 8), complex)/8
    assert np.max(np.abs(executed.apply(rho)-rho)) < 2e-15
    assert executed.half_diamond_error_bound == 0
    assert executed.metadata['resources_per_realization']['native_pulses'] == 3


def test_bound_is_exact_rational_outward_rounded_and_can_be_vacuous():
    result = e.quadrature_bound(64, 18., 25e-6)
    exact = Fraction(32*math.factorial(16), math.factorial(32))*(
        Fraction(18.)*Fraction(25e-6)*63**2)**16
    assert Fraction(result['exact_half_diamond_error_upper']) == exact
    assert Fraction(result['half_diamond_error_upper']) >= exact
    assert 2e-17 < result['half_diamond_error_upper'] < 3e-17
    assert e.quadrature_bound(64, 1., 1., 1)['half_diamond_error_upper'] == 1
    assert result['arbitrary_reference_system_included']
    assert not result['floating_execution_error_included']


@pytest.mark.parametrize('changes', [
    {'dimension': True}, {'dimension': 3}, {'dimension': 128},
    {'rate': True}, {'rate': -1}, {'rate': float('nan')},
    {'duration': True}, {'duration': -1}, {'duration': float('inf')},
    {'quadrature_order': True}, {'quadrature_order': 0},
    {'quadrature_order': 1.5}, {'quadrature_order': 65},
])
def test_invalid_environment_parameters_are_rejected_even_with_warm_cache(changes):
    e.channel(4, 1., 1., 1)
    parameters = dict(dimension=4, rate=1., duration=1., quadrature_order=1)
    parameters.update(changes)
    with pytest.raises(ValueError):
        e.channel(**parameters)
    with pytest.raises(ValueError):
        e.quadrature_bound(**parameters)


@pytest.mark.parametrize('rho', [np.zeros(4), np.zeros((3, 3)), np.zeros((0, 0)),
                                  np.zeros((4, 8)), np.full((4, 4), np.nan)])
def test_invalid_density_shape_and_nonfinite_entries_rejected(rho):
    with pytest.raises(ValueError):
        e.channel(4, 1., .01).apply(rho)


def test_no_host_random_sample_changes_reproducible_ensemble():
    executed = e.channel(4, 1., .01)
    rho = np.ones((4, 4), complex)/4
    assert np.array_equal(executed.apply(rho), executed.apply(rho))
    assert not executed.weights.flags.writeable
    assert not executed.schur.flags.writeable
    assert all(not unitary.flags.writeable for unitary in executed.unitaries)

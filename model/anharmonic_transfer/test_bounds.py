"""Independent channel, cutoff-boundary and conditioning regressions."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from .bounds import (channel_history_bound, conditional_readout_bound,
                     execution_error_budget, strang_bound)


def _operators(d, delta, omega):
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    n = np.kron(np.eye(2), np.diag(np.arange(d)))
    h = delta*n-omega/2*np.kron([[0, 1], [1, 0]], a+a.T)
    return n, h


def _trace_norm(matrix):
    return np.linalg.svd(matrix, compute_uv=False).sum()


@pytest.mark.parametrize('d,subdivisions', [(2, 1), (3, 2), (4, 3)])
def test_complete_split_vs_unsplit_channel_on_entangled_input(d, subdivisions):
    """A reference qubit is kept while actual finite Lindblad maps are compared."""
    delta, omega, gamma, tau, steps = 1.3, 2.1, .7, .07, 2
    n, h = _operators(d, delta, omega)
    size = 2*d
    a = -1j*(np.kron(h, np.eye(size))-np.kron(np.eye(size), h.T))
    numbers = np.tile(np.arange(d), 2)
    diffusion = np.diag(-gamma/2*(numbers[:, None]-numbers[None, :]).ravel()**2)
    subduration = tau/subdivisions
    exact = expm(subduration*(a+diffusion))
    half = expm(subduration*diffusion/2)
    split = half @ expm(subduration*a) @ half
    # Apply both maps to a state maximally entangled across two distinct
    # spin/occupation states and a spectator reference qubit.
    state = np.zeros((size, 2), complex)
    state[0, 0] = 1/math.sqrt(2)
    state[-1, 1] = 1j/math.sqrt(2)
    initial = np.outer(state.ravel(), state.ravel().conj()).reshape(size, 2, size, 2)
    outputs = []
    for channel in (exact, split):
        rho = initial.copy()
        for _ in range(4*steps*subdivisions):
            rho = np.einsum('abcd,crds->arbs', channel.reshape(size, size, size, size),
                            rho, optimize=True)
        outputs.append(rho.reshape(2*size, 2*size))
    error = _trace_norm(outputs[0]-outputs[1])/2
    certificate = strang_bound(d, steps, subdivisions, delta=delta, omega=omega,
                               pulse_duration=tau, dephasing_rate=gamma)
    assert error > 1e-8
    assert error <= certificate['half_diamond_error_upper']
    assert np.trace(outputs[0]).real == pytest.approx(1)
    assert np.trace(outputs[1]).real == pytest.approx(1)


def test_finite_CCR_boundary_cannot_be_discarded():
    d = 16
    n, h = _operators(d, 0, 1)
    comm = lambda x, y: x@y-y@x
    j = -1j*comm(h, n)
    ell = -1j*comm(h, j)
    cutoff = np.diag([.5]*(d-1)+[-(d-1)/2])
    assert np.allclose(ell, np.kron(np.eye(2), cutoff))
    a = lambda rho: -1j*comm(h, rho)
    diffusion = lambda rho: -.5*comm(n, comm(n, rho))
    # A trace-norm-one high/low coherence violates the purported bound that
    # incorrectly uses [a,a†]=I at this finite cutoff.
    rho = np.zeros((2*d, 2*d), complex)
    rho[d-1, 0] = 1
    actual = a(a(diffusion(rho)))-2*a(diffusion(a(rho)))+diffusion(a(a(rho)))
    assert _trace_norm(actual) > 4*(d-1)
    valid = strang_bound(d, 1, 1, delta=0, omega=1, dephasing_rate=1)
    assert _trace_norm(actual) < valid['hamiltonian_double_commutator_diamond_upper']


@pytest.mark.parametrize('zero', ['omega', 'dephasing_rate', 'steps'])
def test_commuting_or_empty_schedule_has_exactly_zero_split_bound(zero):
    parameters = dict(dimension=32, steps=78, subdivisions=128)
    parameters[zero] = 0
    result = strang_bound(**parameters)
    assert result['exact_half_diamond_error_upper'] == '0'


def test_exact_second_order_refinement_and_linear_horizon():
    first = strang_bound(32, 10, 128)
    second = strang_bound(32, 20, 256)
    assert F(second['exact_half_diamond_error_upper'])*2 == F(first['exact_half_diamond_error_upper'])
    assert second['drive_subpulses'] == 4*20*256


@pytest.mark.parametrize('d,subdivisions,low,high',
                         [(32, 128, .000529, .000531), (64, 256, .000478, .000480)])
def test_declared_long_horizon_certificate_is_informative(d, subdivisions, low, high):
    certificate = strang_bound(d, 78, subdivisions)
    assert low < certificate['half_diamond_error_upper'] < high
    for key, value in certificate.items():
        if key.startswith('exact_'):
            assert F(certificate[key[6:]]) >= F(value)
    selected = conditional_readout_bound(F(certificate['exact_half_diamond_error_upper']), F(1, 2))
    assert selected['conditional_binary_probability_error_upper'] < .0011


def test_channel_composition_preserves_exact_small_radii():
    result = channel_history_bound([(80_000, F(1, 10**50)), (3, F(2, 10**40))])
    assert F(result['exact_half_diamond_error_upper']) == F(80_000, 10**50)+F(6, 10**40)
    assert channel_history_bound([])['exact_half_diamond_error_upper'] == '0'
    assert channel_history_bound([(10, F(1, 2))])['exact_half_diamond_error_upper'] == '1'


def test_conditional_bound_is_sharp_for_retained_classical_records():
    q, epsilon = F(1, 5), F(1, 100)
    first = np.array([float(q), 0, float(1-q)])
    second = np.array([float(q-epsilon), float(epsilon), float(1-q)])
    joint = np.abs(first-second).sum()/2
    selected = np.abs(first[:2]/float(q)-second[:2]/float(q)).sum()/2
    result = conditional_readout_bound(epsilon, q)
    assert joint == pytest.approx(float(epsilon))
    assert selected == pytest.approx(result['conditional_detector_tv_upper'])
    assert F(result['exact_conditional_detector_tv_upper']) == epsilon/q
    assert result['reference_selection_guaranteed_positive']


def test_conditioning_handles_different_acceptance_masses():
    rng = np.random.default_rng(706)
    for _ in range(100):
        p, q = rng.dirichlet(np.ones(6), size=2)
        epsilon = np.abs(p-q).sum()/2
        mass = q[:3].sum()
        conditional = np.abs(p[:3]/p[:3].sum()-q[:3]/mass).sum()/2
        bound = conditional_readout_bound(epsilon, mass)
        assert conditional <= bound['conditional_detector_tv_upper']+1e-14


def test_zero_reference_acceptance_is_not_assumed_away():
    result = conditional_readout_bound(F(1, 5), F(1, 5))
    assert not result['reference_selection_guaranteed_positive']
    assert result['conditional_detector_tv_upper'] == 1
    assert result['reference_selected_probability_lower'] == 0


def test_combined_budget_counts_only_declared_environment_invocations():
    result = execution_error_budget(32, 10, 128, F(3, 4),
                                    environment_applications=[(10, F(1, 100_000))])
    exact = F(result['splitting']['exact_half_diamond_error_upper'])+F(1, 10_000)
    assert F(result['exact_complete_trace_distance_upper']) == exact
    assert F(result['selected_readout']['exact_conditional_detector_tv_upper']) == exact/F(3, 4)
    assert not result['occupation_removal_error_included']
    assert not result['floating_arithmetic_error_included']


@pytest.mark.parametrize('field,value', [('dimension', True), ('dimension', 1),
                                       ('steps', -1), ('subdivisions', 0),
                                       ('omega', -1), ('dephasing_rate', float('nan')),
                                       ('pulse_duration', 0), ('delta', float('inf'))])
def test_split_contract_validation(field, value):
    parameters = dict(dimension=32, steps=1, subdivisions=8)
    parameters[field] = value
    with pytest.raises(ValueError):
        strang_bound(**parameters)


@pytest.mark.parametrize('epsilon,acceptance', [(-1, .5), (.1, 0), (.1, 1.1),
                                              (1.1, .5), (.1, float('nan')),
                                              (True, .5)])
def test_conditional_contract_validation(epsilon, acceptance):
    with pytest.raises(ValueError):
        conditional_readout_bound(epsilon, acceptance)


@pytest.mark.parametrize('terms', [[(-1, F(1, 2))], [(1, F(2))], [(True, F(1, 2))]])
def test_history_contract_validation(terms):
    with pytest.raises(ValueError):
        channel_history_bound(terms)

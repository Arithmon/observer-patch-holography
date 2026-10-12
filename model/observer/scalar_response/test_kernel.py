"""Independent oscillator, full-mode and hostile-input retained-state checks."""
from copy import deepcopy
from itertools import product
import json

import numpy as np
import pytest
from scipy.linalg import expm

from .field import ScalarField
from .kernel import MAX_DIRECTIONS, MAX_TERMS, WeylKernel, WeylState


def joint_law(kernel):
    probabilities = {}
    states = {}
    for word in product((-1, 1), repeat=len(kernel.times)):
        state, weight = kernel.vacuum(), 1.
        for index, outcome in enumerate(word):
            probability, state = state.instrument(index, outcome)
            weight *= probability
        probabilities[word] = weight
        states[word] = state
    return probabilities, states


@pytest.fixture(scope="module")
def oscillator():
    """Schrodinger-picture finite oscillator, with a negligible occupation tail."""
    cutoff = 48
    omega, beta, alpha = 1.7, .73, .3+.2j
    times, couplings = (.11, .59, 1.03), (.8, -.6, 1.1)
    annihilation = np.diag(np.sqrt(np.arange(1, cutoff)), 1)
    vacuum = np.eye(cutoff)[:, 0].astype(complex)
    unitary = lambda time: np.diag(np.exp(-1j*omega*time*np.arange(cutoff)))
    operators = []
    for coupling in couplings:
        position = coupling*beta*(annihilation+annihilation.T)
        operators.append({y: (expm(-.5j*position)-1j*y*expm(.5j*position))/2
                          for y in (-1, 1)})
    return omega, beta, alpha, times, couplings, annihilation, vacuum, unitary, operators


@pytest.mark.parametrize("current", (-1, 0, 1))
@pytest.mark.parametrize("reads", (2, 3))
def test_every_joint_probability_conditional_vector_and_gram_against_fock(oscillator, current, reads):
    omega, beta, alpha, times, couplings, annihilation, vacuum, unitary, operators = oscillator
    times, couplings = times[:reads], couplings[:reads]
    kernel = WeylKernel.from_modes([omega], [beta], [alpha], times, current,
                                  couplings=couplings)
    prepared = expm(current*(alpha*annihilation.T-alpha.conjugate()*annihilation)) @ vacuum
    directions = 1j*beta*np.array(couplings)*np.exp(1j*omega*np.array(times))
    probabilities, states = joint_law(kernel)
    assert sum(probabilities.values()) == pytest.approx(1, abs=2e-14)
    basis_cache = {}
    for word, state in states.items():
        expected, previous = prepared.copy(), 0.
        for index, outcome in enumerate(word):
            expected = operators[index][outcome] @ unitary(times[index]-previous) @ expected
            previous = times[index]
        expected_probability = np.vdot(expected, expected).real
        assert probabilities[word] == pytest.approx(expected_probability, rel=0, abs=3e-14)
        # The kernel stores the interaction-picture state. Undo the final
        # free evolution of the independent Schrodinger-picture oracle.
        expected = unitary(-times[-1]) @ expected / np.sqrt(expected_probability)
        assert np.linalg.norm(expected[-8:]) < 2e-14
        basis = []
        for shifts in state.shifts:
            key = tuple(shifts)
            if key not in basis_cache:
                delta = np.dot(shifts, directions)/2
                displacement = expm(delta*annihilation.T-delta.conjugate()*annihilation)
                basis_cache[key] = displacement @ prepared
            basis.append(basis_cache[key])
        basis = np.column_stack(basis)
        np.testing.assert_allclose(state.gram(), basis.conj().T @ basis, rtol=0, atol=3e-14)
        # Comparing the complete vector, including phase, detects a wrong
        # CCR sign even when every two-read probability is unchanged.
        np.testing.assert_allclose(basis @ state.coefficients, expected, rtol=0, atol=4e-14)
    if reads == 2:
        assert kernel.two_read_law() == pytest.approx(probabilities, abs=3e-14)


def test_phase_and_real_only_negative_controls_are_resolved(oscillator):
    omega, beta, alpha, times, couplings, annihilation, vacuum, _, _ = oscillator
    kernel = WeylKernel.from_modes([omega], [beta], [alpha], times, couplings=couplings)
    actual, states = joint_law(kernel)
    commuting = WeylKernel(kernel.covariance.real, kernel.unit_means, times)
    real_only, _ = joint_law(commuting)
    assert max(abs(actual[word]-real_only[word]) for word in actual) > .005

    # A global conjugation of covariance may preserve the classical law;
    # the actual state and the Weyl composition phase still must be checked.
    prepared = expm(alpha*annihilation.T-alpha.conjugate()*annihilation) @ vacuum
    directions = 1j*beta*np.array(couplings)*np.exp(1j*omega*np.array(times))
    state = states[(1, -1, 1)]
    basis = np.column_stack([
        expm(delta*annihilation.T-delta.conjugate()*annihilation) @ prepared
        for delta in state.shifts @ directions/2
    ])
    assert np.max(np.abs(state.gram()-(basis.conj().T @ basis).conj())) > .1

    # Deliberately reverse only the composition phase, retaining the true
    # basis Gram. Three reads expose this broken state-update convention.
    bad = kernel.vacuum()
    bad_weight = 1.
    for index, outcome in enumerate((1, -1, 1)):
        shifts, coefficients = [], []
        for sign, factor in ((-1, .5), (1, -.5j*outcome)):
            moved = bad.shifts.copy()
            moved[:, index] += sign
            phase = np.exp(+.25j*sign*(bad.shifts @ kernel.covariance.imag[index]))
            shifts.append(moved)
            coefficients.append(bad.coefficients*factor*phase)
        bad = WeylState(kernel, np.concatenate(shifts), np.concatenate(coefficients), index)
        bad_weight *= bad.norm_squared()
        bad = bad.normalized()
    assert abs(bad_weight-actual[(1, -1, 1)]) > 1e-5


def test_retention_differs_from_reset_and_matched_donor(oscillator):
    omega, beta, alpha, times, couplings, *_ = oscillator
    kernel = WeylKernel.from_modes([omega], [beta], [alpha], times[:2], couplings=couplings[:2])
    probabilities, _ = joint_law(kernel)
    first = {y: sum(probabilities[y, z] for z in (-1, 1)) for y in (-1, 1)}
    second = {z: sum(probabilities[y, z] for y in (-1, 1)) for z in (-1, 1)}
    donor = {(y, z): first[y]*second[z] for y, z in probabilities}
    reset = {(y, z): first[y]*kernel.vacuum().instrument(1, z)[0] for y, z in probabilities}
    for replacement in (donor, reset):
        assert sum(replacement.values()) == pytest.approx(1, abs=2e-14)
        assert max(abs(probabilities[word]-replacement[word]) for word in probabilities) > .02
        for y in (-1, 1):
            assert sum(replacement[y, z] for z in (-1, 1)) == pytest.approx(first[y], abs=2e-14)
    for z in (-1, 1):
        assert sum(donor[y, z] for y in (-1, 1)) == pytest.approx(second[z], abs=2e-14)


@pytest.mark.parametrize("current", (-1, 0, 1))
def test_q13_tensor_kernel_and_full_mode_conditional_continuation_agree(current):
    field = ScalarField(13)
    times, couplings = (.2, .61, .98), (.6, -.9, 1.1)
    modes = field.mode_data()
    initial = field.prepared_state(0., current=current)
    unit_initial = field.prepared_state(0.).amplitudes[0]
    kernel = WeylKernel.from_field(field, times, current=current, couplings=couplings)
    explicit = WeylKernel.from_modes(modes["omega"], modes["field_vector"], unit_initial,
                                    times, current=current, couplings=couplings)
    assert kernel.mode_count == explicit.mode_count == 1728
    np.testing.assert_allclose(kernel.covariance, explicit.covariance, rtol=0, atol=3e-15)
    np.testing.assert_allclose(kernel.unit_means, explicit.unit_means, rtol=0, atol=3e-15)
    for word in product((-1, 1), repeat=3):
        coherent, retained, previous = initial, kernel.vacuum(), 0.
        for index, outcome in enumerate(word):
            coherent = coherent.evolve(modes["omega"], times[index]-previous)
            full_probability, coherent = coherent.centered_instrument(
                couplings[index]*modes["field_vector"], outcome)
            probability, retained = retained.instrument(index, outcome)
            assert probability == pytest.approx(full_probability, rel=0, abs=8e-15)
            previous = times[index]


def test_checkpoint_restores_future_and_rejects_every_changed_kernel_binding(oscillator):
    omega, beta, alpha, times, couplings, *_ = oscillator
    kernel = WeylKernel.from_modes([omega], [beta], [alpha], times, couplings=couplings)
    _, branch = kernel.vacuum().instrument(0, -1)
    payload = json.loads(json.dumps(branch.checkpoint(), allow_nan=False))
    restored = WeylState.restore(kernel, payload)
    np.testing.assert_array_equal(restored.shifts, branch.shifts)
    np.testing.assert_array_equal(restored.coefficients, branch.coefficients)
    for word in product((-1, 1), repeat=2):
        left, right = branch, restored
        for index, outcome in enumerate(word, start=1):
            lp, left = left.instrument(index, outcome)
            rp, right = right.instrument(index, outcome)
            assert lp == rp
            assert left.checkpoint() == right.checkpoint()
    other_kernels = (
        kernel.for_current(0),
        WeylKernel(kernel.covariance, kernel.unit_means+1e-8, times, mode_count=1),
        WeylKernel(1.01*kernel.covariance, kernel.unit_means, times, mode_count=1),
        WeylKernel(kernel.covariance, kernel.unit_means, (.1, .6, 1.1), mode_count=1),
        WeylKernel(kernel.covariance, kernel.unit_means, times, mode_count=2),
    )
    for other in other_kernels:
        with pytest.raises(ValueError, match="actual field kernel"):
            WeylState.restore(other, payload)


@pytest.mark.parametrize("mutation", (
    {"schema": "other"}, {"extra": 1}, {"last_index": True}, {"last_index": 2},
    {"half_shifts": [[-1., 0], [1., 0]]},
    {"half_shifts": [[True, 0], [1, 0]]},
    {"half_shifts": [[-1, 1], [1, 0]]},
    {"half_shifts": [[-8, 0], [1, 0]]},
    {"coefficients": [[float("nan"), 0], [0, 1]]},
    {"coefficients": [[0, 0], [0, 0]]},
    {"coefficients": [[1, 0, 0], [0, 1, 0]]},
))
def test_checkpoint_rejects_invalid_shape_phase_values_and_norm(mutation):
    kernel = WeylKernel.from_modes([1.], [.7], [.2j], (0., 1.))
    _, state = kernel.vacuum().instrument(0, 1)
    payload = deepcopy(state.checkpoint())
    payload.update(mutation)
    with pytest.raises((ValueError, ArithmeticError)):
        WeylState.restore(kernel, payload)


@pytest.mark.parametrize("times, couplings", (
    ((), None), ((-1.,), None), ((1., 0.), None), ((True,), None),
    ((float("nan"),), None), ((0.,)*8, None), ((0., 1.), (1.,)),
    ((0.,), (float("inf"),)), ((0.,), (False,)),
))
def test_rejects_invalid_schedule_and_couplings(times, couplings):
    with pytest.raises(ValueError):
        WeylKernel.from_modes([1.], [.7], [.2j], times, couplings=couplings)


@pytest.mark.parametrize("omega, beta, alpha", (
    ([0.], [.7], [.2j]), ([1j], [.7], [.2j]),
    ([1.], [.7j], [.2j]), ([1.], [.7], [.2j, 0.]),
    ([float("inf")], [.7], [.2j]), ([1.], [float("nan")], [.2j]),
    (np.ma.array([1.], mask=False), [.7], [.2j]),
))
def test_rejects_invalid_mode_arrays(omega, beta, alpha):
    with pytest.raises(ValueError):
        WeylKernel.from_modes(omega, beta, alpha, (0.,))


def test_rejects_indefinite_covariance_even_when_complex_modulus_overflows():
    covariance = [[1e308, 1.7e308+1.7e308j], [1.7e308-1.7e308j, 1e308]]
    with pytest.raises(ValueError, match="positive semidefinite"):
        WeylKernel(covariance, [0., 0.], (0., 1.))


def test_nonfinite_intermediates_cannot_become_finite_weyl_norms():
    huge = WeylKernel(np.diag([1e308, 1e308]), [0., 0.], (0., 1.))
    overflow = WeylState(huge, [[-7, -7], [7, 7]], [1., 1.], 1)
    with pytest.raises(ArithmeticError):
        overflow.gram()
    ordinary = WeylKernel([[1.]], [0.], (0.,))
    bad_coefficients = WeylState(ordinary, [[0]], [1e308])
    with pytest.raises(ArithmeticError):
        bad_coefficients.normalized()
    with pytest.raises(ArithmeticError):
        WeylKernel.from_modes([1e308], [.7], [.2j], (1e308,))


def test_declared_direction_mode_and_term_budgets_fail_explicitly():
    with pytest.raises(ValueError):
        WeylKernel.from_modes([1.], [.7], [.2j], (0.,)*(MAX_DIRECTIONS+1))
    with pytest.raises(ValueError):
        WeylKernel.from_modes(np.ones(32769), np.ones(32769), np.zeros(32769), (0.,))
    kernel = WeylKernel.from_modes([1.], [.7], [.2j], (0.,)*MAX_DIRECTIONS)
    with pytest.raises(ValueError):
        WeylState(kernel, np.zeros((MAX_TERMS+1, MAX_DIRECTIONS), int), np.ones(MAX_TERMS+1))
    coefficients = np.zeros(MAX_TERMS)
    coefficients[0] = 1
    full = WeylState(kernel, np.zeros((MAX_TERMS, MAX_DIRECTIONS), int), coefficients)
    with pytest.raises(ValueError, match="term budget"):
        full.instrument(0, 1)
    state = kernel.vacuum()
    for index in range(MAX_DIRECTIONS):
        _, state = state.instrument(index, 1)
    assert len(state.coefficients) == MAX_TERMS
    for index in (0, MAX_DIRECTIONS-1, MAX_DIRECTIONS, True, 1.):
        with pytest.raises(ValueError, match="read index"):
            state.instrument(index, 1)


def test_zero_coupling_retains_the_complete_input_ray_and_fair_outcomes():
    kernel = WeylKernel.from_modes([1.], [.7], [.2j], (0., .3), couplings=(.8, 0.))
    for y in (-1, 1):
        _, branch = kernel.vacuum().instrument(0, y)
        for z in (-1, 1):
            probability, final = branch.instrument(1, z)
            assert probability == pytest.approx(.5, abs=2e-15)
            # Zero second direction makes the two copies of every existing
            # basis vector identical; its Kraus scalar is a global phase.
            count = len(branch.coefficients)
            merged = final.coefficients[:count]+final.coefficients[count:]
            overlap = np.vdot(branch.coefficients, branch.gram() @ merged)
            assert abs(overlap)**2 == pytest.approx(1., abs=3e-15)

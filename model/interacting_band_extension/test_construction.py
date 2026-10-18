"""Independent Fourier, compression, matching and observation-bound controls."""
from fractions import Fraction as F
from math import pi

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import eigsh

from . import algebra as a
from .construction import (Lattice, Term, exact_apply, exact_norm_squared,
                           execute_term, history, leakage_witness, normalized,
                           high_parity_instrument, reduced_output_distance, tv)


@pytest.fixture(scope='module')
def extended():
    return Lattice()


def test_exact_fourier_integral_and_low_restriction():
    integral = a.quartic_coefficients()
    assert integral == a.declared_quartic()
    assert len(integral) == 49
    restricted = {powers[:7]: coefficient for powers, coefficient in integral.items()
                  if not any(powers[7:])}
    assert restricted == a.quartic_coefficients(False)
    assert integral[(0, 3, 0, 0, 0, 0, 0, 1, 0)] == 2
    assert integral[(0, 1, 2, 0, 0, 0, 0, 1, 0)] == -6


def test_quartic_grid_threshold_and_alias_counterexample():
    # The zero mode of (sqrt(2)cos(3x))**4 is3/2. At n=12 the
    # +/-12 Fourier components alias into zero, raising it to2.
    for n in (13, 17, 33):
        x = 2*pi*np.arange(n)/n
        assert np.mean((np.sqrt(2)*np.cos(3*x))**4) == pytest.approx(1.5, abs=1e-14)
    x = 2*pi*np.arange(12)/12
    assert np.mean((np.sqrt(2)*np.cos(3*x))**4) == pytest.approx(2.)
    with pytest.raises(ValueError):
        a.grid_history_bound(12, F(1, 8))


def test_pairing_preserves_target_translation_not_rectangular_code():
    values = np.array([.3, -.2, .4, .6, -.1, .2, -.3, .1, -.5])
    angle = .37
    changed = values.copy()
    for i, multiplier in ((1, 1), (7, 3)):
        c, s = np.cos(multiplier*angle), np.sin(multiplier*angle)
        changed[i:i+2] = np.array([[c, -s], [s, c]])@values[i:i+2]
    def polynomial(q):
        return sum(float(coefficient)*np.prod(q**np.array(powers))
                   for powers, coefficient in a.quartic_coefficients().items())
    assert polynomial(changed) == pytest.approx(polynomial(values), abs=2e-14)
    # Pair rotations include a_cos* a_sin: |3,1> has amplitude2 toward
    # |4,0>, which lies outside the rectangular four-by-four code.
    assert np.sqrt((3+1)*1) == 2
    assert 4 not in range(a.DIMENSIONS[1])


def test_polynomial_is_formed_before_occupation_compression():
    q = normalized(a.polynomial(4, 1))
    q2 = normalized(a.polynomial(4, 2))
    q3 = normalized(a.polynomial(4, 3))
    assert np.linalg.norm(q2-q@q) > .9
    assert np.linalg.norm(q3-q@q@q) > 1
    larger = normalized(a.polynomial(8, 1))
    assert np.allclose(q3, np.linalg.matrix_power(larger, 3)[:4, :4])


def test_projected_action_requires_mass_and_scalar_matching(extended):
    matched, bare = Lattice(False, True), Lattice(False)
    block = extended.hamiltonian[:matched.dimension, :matched.dimension]
    scalar = a.residual_bound()['matched_scalar_increment']
    residual = block-matched.hamiltonian-float(scalar)*sparse.eye(matched.dimension)
    assert max(abs(residual.data), default=0) < 1e-13
    difference = matched.hamiltonian-bare.hamiltonian
    expected = sum((bare.embed(Term((i,), (a.polynomial(d, 2),), F(1, 8), 'q2'))
                    for i, d in enumerate(a.LOW_DIMENSIONS)),
                   sparse.csr_matrix(difference.shape))
    assert np.linalg.norm((difference-expected).toarray()) < 1e-13


def test_residual_norm_bound_covers_actual_offdiagonal_block(extended):
    low = np.prod(a.LOW_DIMENSIONS)
    block = extended.hamiltonian[low:, :low]
    actual = eigsh(block.T@block, k=1, which='LA', return_eigenvectors=False)[0]
    assert 0 < actual < float(a.residual_bound()['squared_residual_norm_upper'])
    budget = a.history_bound(F(1, 8))
    assert 0 < budget['reduced_complete_output_trace_distance_upper'] < F(1, 8)
    assert budget['high_retained_complete_output_trace_distance_upper'] > budget['reduced_complete_output_trace_distance_upper']


def test_exact_resonant_leakage_is_interaction_specific(extended):
    first = exact_apply(extended, {(0,)*9: F(1)}, centered=True)
    odd = {key: value for key, value in first.items() if (key[7]+key[8]) % 2}
    assert exact_norm_squared(odd) == F(1, 768)
    free = Lattice(coupling=0)
    first_free = exact_apply(free, {(0,)*9: F(1)}, centered=True)
    assert not any((key[7]+key[8]) % 2 for key in first_free)
    # Free high even excitations are present and must not be called an
    # interaction witness.
    assert any(key[7]+key[8] for key in first_free)
    witness = leakage_witness(extended, execute=False)
    assert witness['ideal_source_odd_high_probability_lower'] > 0


def test_exact_sparse_derivative_matches_independent_hamiltonian(extended):
    vector = exact_apply(extended, {(0,)*9: F(1)})
    numeric = np.zeros(extended.dimension)
    strides = np.cumprod((1,)+a.DIMENSIONS[:-1])
    for key, value in vector.items():
        from math import factorial, prod
        numeric[int(np.dot(key, strides))] = float(value)*np.sqrt(prod(factorial(n) for n in key))
    assert np.max(abs(numeric-extended.hamiltonian[:, 0].toarray().ravel())) < 1e-13


def test_local_three_mode_word_and_full_phase_are_executed(extended):
    term = next(term for term in extended.terms if len(term.sites) == 3)
    duration = F(1, 2048)
    matrix, metadata = execute_term(term, duration)
    assert matrix.shape == (64, 64)
    assert np.linalg.norm(matrix-expm(-1j*float(duration)*term.matrix()), 2) < 2e-10
    assert metadata['resources']['native_pulses'] > 0
    assert metadata['local_register_qubits'] == 6


def test_resonance_readout_is_a_complete_executed_parity_instrument():
    kraus, metadata = high_parity_instrument()
    assert np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(16), 2) < 1e-12
    for index in range(16):
        bit = (index % 4+index//4) % 2
        assert np.linalg.norm(kraus[bit][:, index])**2 == pytest.approx(1.)
        assert np.linalg.norm(kraus[1-bit][:, index])**2 < 1e-24
    assert metadata['resources']['native_pulses'] > 0


def _partial_high(state, low=2):
    matrix = state.reshape((-1, low))
    return matrix.T@matrix.conj()


def _trace_distance(left, right):
    return np.sum(abs(np.linalg.eigvalsh(left-right)))/2


def test_quadratic_reduced_bound_and_linear_full_output_are_distinct():
    x = np.array([[0., 1.], [1., 0.]])
    angle = .14
    state = expm(-1j*angle*np.kron(x, x))@np.array([1., 0., 0., 0.])
    low = _partial_high(state)
    target = np.diag([1., 0.])
    distance = _trace_distance(low, target)
    full_distance = np.sqrt(1-abs(state[0])**2)
    assert distance == pytest.approx(np.sin(angle)**2)
    assert distance <= angle**2
    assert full_distance == pytest.approx(abs(np.sin(angle)))
    assert full_distance > angle**2  # Quadratic bound cannot cover high output.


def test_omitting_the_projected_action_can_fail_at_linear_order():
    z = np.diag([1., -1.])
    state = np.array([1., 1.])/np.sqrt(2)
    angle = .03
    evolved = expm(-1j*angle*z)@state
    distance = _trace_distance(np.outer(evolved, evolved.conj()), np.outer(state, state))
    assert distance == pytest.approx(np.sin(angle))
    assert distance > angle**2
    # H=Z_low tensor I_high has R=0: only the PHP target has zero error.
    assert np.linalg.norm(evolved-expm(-1j*angle*z)@state) == 0


def test_measurement_feedback_keeps_complete_quadratic_bound():
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    times = (.09, .07)
    # An untouched reference begins maximally entangled with the low qubit.
    # Tensor order is high, low, reference.
    h = np.kron(x, np.kron(x, np.eye(2)))
    kraus = (np.diag([np.sqrt(.8), np.sqrt(.3)]), np.diag([np.sqrt(.2), np.sqrt(.7)]))
    initial = np.array([1., 0., 0., 1., 0., 0., 0., 0.], complex)/np.sqrt(2)
    first = expm(-1j*times[0]*h)@initial
    difference = 0.
    mass = 0.
    for bit, k in enumerate(kraus):
        controlled = expm(-.23j*z) if bit else np.eye(2)
        low_operation = np.kron(controlled@k, np.eye(2))
        full = expm(-1j*times[1]*h)@np.kron(np.eye(2), low_operation)@first
        low = low_operation@initial[:4]
        rho = _partial_high(full, low=4)
        difference += np.sum(abs(np.linalg.eigvalsh(rho-np.outer(low, low.conj()))))/2
        mass += np.trace(rho).real
    assert mass == pytest.approx(1.)
    assert 0 < difference <= sum(times)**2


def test_same_low_protocol_retains_high_memory(extended):
    principal = history(extended, source=False)
    reset = history(extended, source=False, reset_high=True)
    assert sum(row['probability'] for row in principal) == pytest.approx(1.)
    assert sum(row['probability'] for row in reset) == pytest.approx(1.)
    assert tv(principal, reset) > 1e-6


def test_quantum_output_comparison_detects_identical_record_probabilities():
    low_size = np.prod(a.LOW_DIMENSIONS)
    full = np.zeros(np.prod(a.DIMENSIONS), complex)
    full[0] = 1
    same = np.zeros(low_size, complex)
    same[0] = 1
    different = np.zeros(low_size, complex)
    different[1] = 1
    assert reduced_output_distance([full], [same]) == pytest.approx(0.)
    assert reduced_output_distance([full], [different]) == pytest.approx(1.)
    # Both pairs have the same single-record probability one.
    assert np.linalg.norm(full) == np.linalg.norm(same) == np.linalg.norm(different)


def test_grid_bound_refines_quadratically_without_enlarging_band():
    coarse = a.grid_history_bound(65, F(1, 8))
    fine = a.grid_history_bound(130, F(1, 8))
    assert coarse['same_nine_mode_hamiltonian_norm_difference_upper'] == 4*fine['same_nine_mode_hamiltonian_norm_difference_upper']
    assert coarse['number_of_independent_modes'] == fine['number_of_independent_modes'] == 9
    assert not fine['grid_increases_bandwidth']
    for n in (13, 65):
        # Exact norm of the positive sum on separate oscillator factors,
        # evaluated independently from the rational coefficient estimate.
        actual = 0.
        for wave, d in zip(a.WAVES, a.DIMENSIONS, strict=True):
            k = np.linalg.norm(wave)
            kinetic_error = k*k-(n*np.sin(pi*k/n)/pi)**2
            actual += kinetic_error*np.linalg.eigvalsh(normalized(a.polynomial(d, 2)))[-1]/2
        certificate = a.grid_history_bound(n, F(1, 8))
        assert 0 < actual <= float(certificate['same_nine_mode_hamiltonian_norm_difference_upper'])


@pytest.mark.parametrize('value', [True, 1.5, -1])
def test_invalid_exact_inputs_are_rejected(value):
    with pytest.raises(ValueError):
        a.history_bound(value)
    with pytest.raises(ValueError):
        a.grid_history_bound(value, F(1, 8))

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
                           high_parity_instrument, reduced_output_distance, tv,
                           ProductBand, product_gibbs_preparation,
                           product_interaction_control, product_parameter_recovery,
                           product_incident_tapes, product_incident_clock_control,
                           product_incident_clock_bound)


def test_product_radial_integral_and_compression_against_padded_oscillators():
    from scipy.integrate import quad
    for left, right, expected in ((4, 0, 1.5), (2, 2, 1.), (3, 1, 0.), (1, 3, 0.)):
        integral = quad(lambda r: (np.sqrt(2)*np.sin(r))**left
                        *(np.sqrt(2)*np.sin(2*r))**right/pi, 0, pi)[0]
        assert integral == pytest.approx(expected, abs=1e-13)
    q = normalized(a.polynomial(8, 1))
    annihilation = np.diag(np.sqrt(np.arange(1, 8)), 1)
    p = 1j*(annihilation.T-annihilation)
    q2, q4, p2 = (matrix[:4, :4] for matrix in (q@q, np.linalg.matrix_power(q, 4), p@p))
    h1, h2 = (p2/2+k*k*q2/2+q4/16 for k in (1, 2))
    expected = np.kron(np.eye(4), h1)+np.kron(h2, np.eye(4))+np.kron(q2, q2)/4
    assert np.allclose(ProductBand().hamiltonian.toarray(), expected, atol=2e-14)
    assert np.linalg.eigvalsh(expected)[0] > 0
    assert np.linalg.norm(q4-np.linalg.matrix_power(q[:4, :4], 4)) > 1


def test_product_gibbs_source_prepares_the_interacting_state_and_keeps_reference():
    lattice = ProductBand()
    raw, target, metadata = product_gibbs_preparation(lattice)
    exact = expm(-lattice.hamiltonian.toarray()/4)
    exact /= np.trace(exact)
    assert raw.shape == target.shape == (16, 16)
    assert np.linalg.norm(raw-target) < 1e-9
    assert np.allclose(target@target.conj().T, exact, atol=1e-13)
    assert np.linalg.matrix_rank(raw) == 16
    assert metadata['source_resources']['native_pulses'] > 0
    assert not metadata['physical_preparation_error_certified']
    with pytest.raises(ArithmeticError, match='faithful Gibbs'):
        product_gibbs_preparation(lattice, 10**6)


def test_product_complete_records_interaction_clock_and_feedback_controls():
    result = product_interaction_control()
    assert result['complete_work_reference_record_trace_error_numerical'] < 1e-7
    assert result['ideal_splitting_complete_output_bound'] < F(1, 1000)
    assert result['interaction_removal_record_tv'] > 7e-4
    assert result['feedback_removal_record_tv'] > 1e-3
    assert result['ideal_Gibbs_modular_identity_residual'] < 1e-10
    assert result['raw_prepared_marginal_modular_identity_residual'] < 1e-9
    assert result['free_Gibbs_state_interacting_commutator_norm'] > .05
    assert result['native_history_count'] == 1056
    assert result['native_history_probability_sum'] == pytest.approx(1., abs=1e-10)
    assert not result['finite_band_is_invariant_under_untruncated_local_interaction']
    assert not result['observed_QFT_or_laboratory_clock_identified']


def test_product_nonlinearity_leaves_the_declared_radial_band():
    from scipy.integrate import quad
    # The cubic equation from even a single retained mode produces mode3.
    leaked = quad(lambda r: (np.sqrt(2)*np.sin(r))**3*np.sqrt(2)*np.sin(3*r)/pi, 0, pi)[0]
    assert leaked == pytest.approx(-.5, abs=1e-13)


def test_product_parameters_are_recovered_from_source_records_with_a_withheld_interval():
    raw, _, _ = product_gibbs_preparation(ProductBand())
    laws = [history(ProductBand(), sign, initial=raw) for sign in (-1, 1)]
    result = product_parameter_recovery(laws, raw)
    assert result['estimates_rounded_to_1e_minus_4'] == [1., .25, 1.]
    assert min(result['local_sensitivity_singular_values_rounded_to_1e_minus_6']) > 1e-5
    assert not result['physical_parameter_validation']
    assert not result['radial_topology_or_angular_tensor_identified']


def test_product_records_do_not_fix_the_absolute_energy_and_time_unit():
    h = ProductBand().hamiltonian.toarray()
    beta, time, scale = .25, .0625, 7.
    assert np.allclose(expm(-beta*h), expm(-(beta/scale)*(scale*h)))
    assert np.allclose(expm(-1j*time*h), expm(-1j*(time/scale)*(scale*h)))


def test_native_product_clock_keeps_every_endpoint_coin_and_zero_current_incident():
    from ..observer.quantum_mixing import seams
    tapes = product_incident_tapes()
    assert len(tapes) == len(set(tapes)) == 3600
    counts = [0, 0, 0]
    null_incidents = 0
    for first, second, final in tapes:
        loads = [0]*10+[1, 1]
        for edge, coin, current, incident in (first, second):
            left, right = seams()[edge]
            total = loads[left]+loads[right]
            # Independent balancing rule, expressed as the final endpoint
            # allocation instead of calling the native current implementation.
            expected_left = total//2 + coin*(total % 2)
            assert current == loads[left]-expected_left
            loads[left], loads[right] = expected_left, total-expected_left
            assert incident == int(10 in (left, right))
            null_incidents += int(incident and current == 0)
        assert tuple(loads) == final
        counts[first[3]+second[3]] += 1
    assert counts == [2500, 1000, 100]
    assert null_incidents > 0
    # This is the exact two-attempt binomial clock, not an inserted Poisson law.
    assert abs(counts[0]/3600-np.exp(-1/3)) > .01


def test_incident_controller_retains_timeout_work_and_reference_outputs():
    result = product_incident_clock_control()
    assert result['complete_native_tape_and_field_prefix_outcomes'] == 4900
    assert result['timeout_probability'] == F(35, 36)
    rows = result['complete_clock_and_field_record_law']
    assert [sum(row['incident_ticks'] == tick for row in rows) for tick in range(3)] == [1, 2, 4]
    assert sum(row['probability'] for row in rows) == pytest.approx(1, abs=1e-10)
    assert sum(row['probability'] for row in rows if row['timeout']) == pytest.approx(35/36, abs=1e-10)
    assert result['retained_work_addresses'] == result['retained_reference_addresses'] == 16
    assert result['complete_native_tape_work_reference_record_trace_error_numerical'] < 1e-7
    assert result['nonzero_current_tick_marginal_tv'] > .1
    assert result['feedback_removal_complete_output_distance'] > 1e-3
    assert result['source_pauses_during_field_services']
    assert not result['fixture_is_full_refined_screen_clock']
    assert not result['physical_service_or_lab_second_calibration_certified']


def test_clock_controller_oracle_rejects_counting_ticks_without_interacting_evolution(monkeypatch):
    monkeypatch.setattr(ProductBand, 'source_evolve', lambda self, state, duration, steps: state)
    with pytest.raises(ArithmeticError, match='joint incident-clock'):
        product_incident_clock_control()


def test_same_interacting_generator_clock_schedule_has_vanishing_complete_budget():
    rows = [product_incident_clock_bound(level) for level in (2, 5, 8)]
    errors = [row['nominal_work_law_with_retained_clock_tape_error_upper'] for row in rows]
    assert errors[0] > errors[1] > errors[2] > 0
    assert errors[-1] < F(1, 1000)
    for row in rows:
        rate = row['incident_rate_per_normalized_attempt_unit']
        assert rate == 10*4**row['level']
        assert row['model_time_per_zero_inclusive_incident'] == F(1, rate)
        assert row['modular_increment_per_incident_at_beta_one_quarter'] == -F(4, rate)
        assert row['total_source_instrument_error_upper'] == F(2, rate**2)
        assert row['same_occupation_H_at_every_clock_level']
        assert not row['source_services_and_storage_accuracy_certified']
    # Occupation and clock limits require a diagonal budget; a fixed clock
    # schedule cannot absorb an arbitrarily growing interacting generator.
    larger = product_incident_clock_bound(5, occupation=8)
    assert larger['centered_generator_norm_sum_upper'] > rows[1]['centered_generator_norm_sum_upper']
    assert larger['same_H_product_formula_error_upper'] > rows[1]['same_H_product_formula_error_upper']
    nonaligned = product_incident_clock_bound(3, interval=F(1, 17))
    assert nonaligned['nominal_read_time_rounding_error_upper'] > 0
    for kwargs in ({'degree': True}, {'degree': 7}, {'steps': 0}, {'tolerance': F(0)}):
        with pytest.raises(ValueError):
            product_incident_clock_bound(2, **kwargs)
    with pytest.raises(ValueError, match='occupation two, four or eight'):
        product_gibbs_preparation(ProductBand(occupation=6))


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


@pytest.mark.parametrize('coupling,beta', [(0, F(1, 4)), (1, F(1, 4)), (3, F(1))])
def test_product_occupation_harmonic_eigenvalue_and_exact_Gibbs_tail_bounds(coupling, beta):
    d = 12
    eigenvalues = np.linalg.eigvalsh(ProductBand(coupling, occupation=d).hamiltonian.toarray())
    harmonic = sorted(i+2*j+1.5 for i in range(d*d) for j in range(d*d))[:d*d]
    assert np.min(eigenvalues-np.array(harmonic)) > -1e-11
    weights = np.exp(-float(beta)*eigenvalues)
    partition = float(weights.sum())
    for k in (0, 4, 8, 16):
        bound = a.product_gibbs_tail(beta, coupling, k)
        assert float(bound['partition_function_lower']) <= partition
        assert partition <= float(bound['partition_function_upper'])
        rank = bound['retained_spectral_rank']
        assert weights[rank:].sum()/partition <= float(bound['normalized_spectral_tail_upper'])+1e-15
        assert not bound['occupation_projection_or_record_error_certified']
    assert a.product_gibbs_tail(beta, coupling, 32)['normalized_spectral_tail_upper'] < F(1, 1000)


def test_product_occupation_is_nested_full_polynomial_compression():
    small, large = ProductBand(occupation=4), ProductBand(occupation=8)
    addresses = [i+8*j for j in range(4) for i in range(4)]
    large_matrix = large.hamiltonian.toarray()
    assert np.allclose(large_matrix[np.ix_(addresses, addresses)], small.hamiltonian.toarray(), atol=1e-13)
    low = np.linalg.eigvalsh(small.hamiltonian.toarray())
    high = np.linalg.eigvalsh(large_matrix)
    assert np.all(high[:len(low)] <= low+1e-12)
    assert np.exp(-high/4).sum() > np.exp(-low/4).sum()
    # Taking powers after projection changes even the finite target action.
    q = normalized(a.polynomial(4, 1))
    assert np.linalg.norm(normalized(a.polynomial(4, 4))-q@q@q@q) > 1


def test_product_reference_alignment_is_needed_for_complete_state_convergence():
    energies, basis = np.linalg.eigh(ProductBand().hamiltonian.toarray())
    weights = np.exp(-energies/4); weights /= weights.sum()
    amplitude = basis*np.sqrt(weights)
    changed_basis = basis.copy(); changed_basis[:, ::2] *= -1
    changed = changed_basis*np.sqrt(weights)
    assert np.allclose(amplitude@amplitude.T, changed@changed.T, atol=1e-14)
    assert np.linalg.norm(amplitude-changed) > 1
    # Applying each conjugate eigenbasis to its retained reference gives
    # the same canonical square root, without tracing out any amplitude.
    assert np.allclose(amplitude@basis.T, changed@changed_basis.T, atol=1e-14)


def test_product_occupation_refinement_keeps_the_complete_ideal_history():
    result = a.product_occupation_refinement()
    assert [row['occupation_per_mode'] for row in result['dimensions']] == [4, 8, 12]
    for row in result['dimensions']:
        assert row['canonical_reference_norm_squared'] == pytest.approx(1., abs=1e-13)
        assert sum(x['probability'] for x in row['complete_two_read_law']) == pytest.approx(1., abs=1e-13)
        assert [x['bits'] for x in row['complete_two_read_law']] == [[0, 0], [0, 1], [1, 0], [1, 1]]
    distances = [row['complete_record_work_reference_trace_distance_numerical']
                 for row in result['successive_comparisons']]
    assert 0 < distances[1] < distances[0] < 1
    assert not result['numerical_comparisons_certify_infinite_occupation_error']
    assert not result['growing_occupation_source_execution']


def test_product_occupation_refinement_does_not_reuse_a_fixed_source_error_budget():
    from .construction import splitting_bound
    small = splitting_bound(ProductBand(), F(1, 16), 64)
    large = splitting_bound(ProductBand(occupation=12), F(1, 16), 64)
    assert large['centered_generator_norm_sum_upper'] > 5*small['centered_generator_norm_sum_upper']
    assert large['one_interval_operator_error_upper'] > small['one_interval_operator_error_upper']
    with pytest.raises(ValueError, match='two, four or eight'):
        product_gibbs_preparation(ProductBand(occupation=12))


@pytest.mark.parametrize('beta,coupling,shell', [(0, 1, 8), (True, 1, 8), (F(1, 4), -1, 8),
                                                (F(1, 4), 1, -1), (F(1, 4), 1, 2.5)])
def test_product_Gibbs_tail_requires_declared_exact_positive_temperature(beta, coupling, shell):
    with pytest.raises(ValueError):
        a.product_gibbs_tail(beta, coupling, shell)


def test_product_harmonic_shell_count_and_tail_are_an_exact_series_identity():
    q = F(4, 5)
    total = 1/((1-q)*(1-q*q))
    for k in (0, 1, 3, 8):
        labels = [(i, j) for i in range(2*k) for j in range(k) if i+2*j < 2*k]
        assert len(labels) == k*(k+1)
        head = sum((q**(i+2*j) for i, j in labels), F(0))
        bound = a.product_gibbs_tail(shell=k)
        assert bound['unnormalized_spectral_tail_upper'] == total-head

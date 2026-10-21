"""Full-register clock, intersector pulse and spectral detector certificates."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

from . import prepared_clock_source as m


@pytest.fixture(scope='module')
def source():
    return m.compiled_source()


def test_spectrum_comes_from_exact_native_seams_with_full_constant_and_padding():
    row = m.native_spectral_data()
    independent = sp.zeros(12)
    for a, b in row['seams']:
        vector = sp.zeros(12, 1)
        vector[a], vector[b] = 1, -1
        independent += vector*vector.T
    assert independent == row['laplacian']
    assert independent.eigenvals() == {0: 1, 5-sp.sqrt(5): 3, 6: 5, 5+sp.sqrt(5): 3}
    assert tuple(row['basis_exact'][:12, 0]) == (1/sp.sqrt(12),)*12
    for band, addresses in enumerate(row['band_addresses']):
        for address in addresses:
            vector = row['basis_exact'][:12, address]
            assert all(sp.simplify(x) == 0
                       for x in independent*vector-row['band_eigenvalues'][band]*vector)
    np.testing.assert_array_equal(row['basis'][12:, 12:], np.eye(4))
    np.testing.assert_array_equal(row['basis'][:12, 12:], np.zeros((12, 4)))


def test_reference_is_unnormalized_sampler_polynomial_before_global_normalization(source):
    row, spectrum = source['reference'], source['spectrum']
    lazy = sp.eye(12)-spectrum['laplacian']/120
    independent = sum((sp.Rational(w)*lazy**(2*n) for n, w in enumerate(row['weights'])),
                      sp.zeros(12))
    assert row['polynomial'] == independent
    assert row['density_exact'] == independent/sp.trace(independent)
    assert row['m0'] == sum(row['weights'])
    assert all(x > 0 for x in np.linalg.eigvalsh(np.array(row['density_exact'], float)))
    # Normalizing each sampled component before mixing changes the reference.
    wrong = sum((float(w)*np.array(lazy**(2*n), float)/float(sp.trace(lazy**(2*n)))
                 for n, w in enumerate(row['weights'])), np.zeros((12, 12)))
    wrong /= float(sum(row['weights']))
    assert np.linalg.norm(wrong-row['density'][:12, :12], 'fro') > 1e-5


def test_raw_complete_clocks_match_independent_matrix_functions_and_sign(source):
    laplacian = np.array(source['spectrum']['laplacian'], float)
    values, vectors = np.linalg.eigh(laplacian)
    values = np.maximum(values, 0)
    weights = source['reference']['weights']
    masses = sum(float(w)*(1-values/120)**(2*n) for n, w in enumerate(weights))
    hamiltonian = (vectors*np.log(float(sum(weights))/masses))@vectors.T
    square_root = (vectors*np.sqrt(values))@vectors.T
    for name, generator in (('clock1', hamiltonian), ('sqrt_clock1', square_root)):
        expected = np.eye(16, dtype=complex)
        expected[:12, :12] = expm(-1j*generator)
        np.testing.assert_allclose(source['raw'][name], expected,
                                   atol=float(source['guards'][name])+2e-8, rtol=0)
        assert np.linalg.norm(source['raw'][name]-expected.conj().T, 'fro') > .1
        np.testing.assert_allclose(source['raw'][name][12:, 12:], np.eye(4), atol=1e-9, rtol=0)
    assert np.linalg.norm(source['raw']['clock1']-source['raw']['sqrt_clock1'], 'fro') > 1
    state = source['reference']['density']
    actual = source['raw']['clock1']
    assert np.linalg.norm(actual@state@actual.conj().T-state, 'fro') < 1e-9
    # On the faithful support, sigma_1(A)=alpha_-1(A), not alpha_+1(A).
    density_eigenvalues, density_basis = np.linalg.eigh(state[:12, :12])
    modular = (density_basis*np.exp(1j*np.log(density_eigenvalues)))@density_basis.conj().T
    observable = np.zeros((12, 12)); observable[0, 1] = observable[1, 0] = 1
    lhs = modular@observable@modular.conj().T
    rhs = actual[:12, :12]@observable@actual[:12, :12].conj().T
    np.testing.assert_allclose(lhs, rhs, atol=1e-9, rtol=0)


@pytest.mark.parametrize('band', (1, 2, 3))
def test_same_register_kicks_and_complete_detectors_on_entangled_columns(source, band):
    basis = source['spectrum']['basis']
    address = source['pair_addresses'][band]
    x = np.zeros((16, 16), complex); x[0, address] = x[address, 0] = 1
    y = np.zeros((16, 16), complex); y[0, address], y[address, 0] = -1j, 1j
    for sign, name in ((1, f'kick_plus_{band}'), (-1, f'kick_minus_{band}')):
        expected = basis@expm(sign*1j*np.pi*x/4)@basis.T
        np.testing.assert_allclose(source['raw'][name], expected, atol=1e-9, rtol=0)
    zero_minus_pair = np.zeros((16, 16)); zero_minus_pair[0, 0] = 1
    zero_minus_pair[address, address] = -1
    for axis, observable in (('x', x), ('y', y)):
        detector = source['raw'][f'detector_{axis}_{band}']
        np.testing.assert_allclose(detector.conj().T@zero_minus_pair@detector,
                                   basis@observable@basis.T, atol=1e-9, rtol=0)
    # Every system/reference column continues linearly, including actual padding.
    rng = np.random.default_rng(180+band)
    columns = rng.normal(size=(16, 12))+1j*rng.normal(size=(16, 12))
    columns /= np.linalg.norm(columns)
    detector = source['raw'][f'detector_x_{band}']
    output = detector@source['raw']['clock1']@source['raw'][f'kick_plus_{band}']@columns
    expected = (source['targets'][f'detector_x_{band}']@source['targets']['clock1']
                @source['targets'][f'kick_plus_{band}']@columns)
    assert np.linalg.norm(output-expected, 'fro') < 1e-9
    probabilities = np.sum(abs(output)**2, axis=1)
    assert len(probabilities) == 16 and min(probabilities) > 0
    assert abs(sum(probabilities)-1) < 1e-9
    assert sum(probabilities[12:]) > .1


def test_wrong_phase_sign_and_padding_are_rejected_by_full_source_guards(source):
    name = 'modular_phase'
    wrong_sign = m.compiler.inverse(source['words'][name])
    with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
        m._execute(wrong_sign, source['targets'][name], source['guards'][name], 'wrong_sign')
    wrong_padding = source['raw']['clock1'].copy()
    wrong_padding[12, 12] *= -1
    with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
        m._check_matrix(wrong_padding, source['targets']['clock1'],
                        source['guards']['clock1'], 'wrong_padding')
    dropped_padding = source['raw']['detector_x_2'][:12, :12]
    with pytest.raises(ArithmeticError, match='complete finite register'):
        m._check_matrix(dropped_padding, source['targets']['detector_x_2'],
                        source['guards']['detector_x_2'], 'dropped_padding')


@pytest.mark.parametrize('target', (
    np.full((16, 16), np.nan, dtype=complex),
    np.full((16, 16), np.inf, dtype=complex),
    np.eye(12, dtype=complex),
    np.eye(16, dtype=bool),
    np.eye(16, dtype=object),
))
def test_matrix_guard_rejects_invalid_targets_before_comparing_norms(target):
    with pytest.raises(ArithmeticError, match='target .*complete finite register'):
        m._check_matrix(np.eye(16, dtype=complex), target, F(1, 10**9), 'bad_target')


@pytest.mark.parametrize('guard', (float('nan'), float('inf'), True, 1., F(-1), F(10**400)))
def test_matrix_guard_requires_a_finite_nonnegative_rational_allowance(guard):
    with pytest.raises(ArithmeticError, match='finite nonnegative rational guard'):
        m._check_matrix(2*np.eye(16, dtype=complex), np.eye(16), guard, 'bad_guard')


def test_matrix_guard_rejects_finite_input_overflow_without_runtime_warnings():
    import warnings
    raw = np.full((16, 16), 1e308, dtype=complex)
    assert np.isfinite(raw).all()
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
            m._check_matrix(raw, np.eye(16), F(1, 10**9), 'overflow')


def test_matrix_guard_accepts_an_exact_zero_allowance():
    result = m._check_matrix(np.eye(16, dtype=int), np.eye(16), F(0), 'exact_identity')
    assert result['complete_target_Frobenius_error'] == 0
    assert result['raw_unitarity_Frobenius_defect'] == 0


def test_matrix_guard_does_not_accept_unitarity_from_wrapped_integer_products():
    raw = np.eye(16, dtype=np.int64)*np.iinfo(np.int64).max
    np.testing.assert_array_equal(raw.T@raw, np.eye(16))
    with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
        m._check_matrix(raw, raw.copy(), F(0), 'integer_wrap')


def test_matrix_guard_does_not_erase_small_float32_errors_by_squaring():
    target = np.eye(16, dtype=np.float32)
    raw = target.copy()
    raw[0, 1] = 1e-25
    assert np.linalg.norm(raw-target, 'fro') == 0
    with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
        m._check_matrix(raw, target, F(1, 10**30), 'underflow')


@pytest.mark.parametrize('name', ('clock1', 'kick_plus_2', 'detector_x_2', 'detector_y_2'))
def test_literal_continuation_words_agree_with_raw_routines_without_polar_correction(
        source, monkeypatch, name):
    def forbidden(*_args, **_kwargs):
        raise AssertionError('polar correction must not replace executed source maps')
    monkeypatch.setattr(m.compiler, 'stabilize_executed_unitary', forbidden)
    literal = m.compiler.execute(source['words'][name], 4)
    np.testing.assert_allclose(literal, source['raw'][name], atol=1e-10, rtol=0)
    assert source['raw_source_routines_used_without_polar_correction']


def test_public_results_are_detached_and_api_rejects_bool_and_unselected_settings(source):
    changed = m.compiled_source()
    changed['raw']['clock1'][0, 0] = 17
    changed['words']['clock1'][0]['angle'] = 17
    fresh = m.compiled_source()
    np.testing.assert_array_equal(fresh['raw']['clock1'], source['raw']['clock1'])
    assert fresh['words']['clock1'][0] == source['words']['clock1'][0]
    for band in (True, 0, 4, 1.):
        with pytest.raises(ValueError):
            m.detector_word(band)
    with pytest.raises(ValueError):
        m.pair_word(1, True)
    with pytest.raises(ValueError):
        m.clock_word('fitted')
    with pytest.raises(ValueError):
        m.detector_word(1, 'Z')
    for value in (np.str_('X'), ['X'], 1, True):
        with pytest.raises(ValueError):
            m.detector_word(1, value)
        with pytest.raises(ValueError):
            m._pair_spectral_word(1, value, F(1, 4))
    with pytest.raises(ValueError):
        m.clock_word(np.str_('modular'))
    for wrong in (np.eye(16, dtype=bool), np.eye(16, dtype=object),
                  np.full((16, 16), '0')):
        with pytest.raises(ArithmeticError, match='complete finite register'):
            m._check_matrix(wrong, np.eye(16), F(1, 100), 'wrong_dtype')


def test_exact_basis_source_words_do_not_depend_on_platform_libm_angle_paths(source, monkeypatch):
    def forbidden(*_args):
        raise AssertionError('basis synthesis used platform libm angle arithmetic')
    for name in ('atan2', 'cos', 'sin', 'acos'):
        monkeypatch.setattr(m.math, name, forbidden)
    assert m._basis_word(source['spectrum']['basis_exact']) == source['words']['basis']
    certificate = source['spectrum']['exact_basis_certificate']
    assert certificate['symmetric_source_laplacian']
    assert certificate['all_twelve_eigenvector_equations_checked_exactly']
    assert certificate['within_each_band_orthogonality_and_positive_norms_checked_exactly']


@pytest.mark.parametrize('axis', ('X', 'Y'))
def test_controlled_ordered_pairs_are_actual_complete_pointer_readouts(source, axis):
    factors = ((1, 'X', F(1, 3)), (1, 'Y', F(1, 4)), (2, 'X', F(1, 5)))
    result = m.pair_pointer_source(factors, axis)
    spectral = np.eye(16, dtype=complex)
    for band, direction, angle in factors:
        address = source['pair_addresses'][band]
        generator = np.zeros((16, 16), complex)
        generator[0, address] = 1 if direction == 'X' else -1j
        generator[address, 0] = generator[0, address].conjugate()
        spectral = expm(1j*float(angle)*generator)@spectral
    basis = source['spectrum']['basis']
    ordered = basis@spectral@basis.T
    values, vectors = np.linalg.eigh(source['reference']['density'][:12, :12])
    initial = np.zeros((32, 12), complex)
    initial[:12] = vectors*np.sqrt(values)
    output = result['raw']@initial
    probabilities = np.sum(abs(output)**2, axis=1)
    characteristic = np.trace(source['reference']['density']@ordered)
    expected = characteristic.real if axis == 'X' else characteristic.imag
    assert abs(sum(probabilities[:16])-sum(probabilities[16:])-expected) < 1e-10
    assert abs(sum(probabilities)-1) < 1e-10
    assert len(probabilities) == result['complete_outcome_count'] == 32
    # Literal source execution checks the cached composition and pointer-one
    # inputs too, not just the prepared pointer-zero sector.
    literal = m.compiler.execute(result['word'], 5)
    np.testing.assert_allclose(literal, result['raw'], atol=1e-10, rtol=0)
    assert result['pulse_count'] == len(result['word'])
    assert m.compiler.resources(result['word'])['proper_code_service_events'] == 2*len(result['word'])
    assert result['guard'] > result['diagnostic']['complete_target_Frobenius_error']
    if axis == 'Y':
        reversed_order = m.pair_pointer_source(tuple(reversed(factors)), axis)
        flipped = reversed_order['raw']@initial
        contrast = np.sum(abs(flipped[:16])**2)-np.sum(abs(flipped[16:])**2)
        assert abs(contrast-expected) > 100*float(result['guard'])


@pytest.mark.parametrize('factors,axis', (
    ((), 'X'), (((True, 'X', F(1)),), 'X'), (((1, 'Z', F(1)),), 'X'),
    (((1, 'X', .25),), 'X'), (((1, 'X', F(65)),), 'X'),
    (((1, 'X', F(1)),), np.str_('Y')),
))
def test_pair_pointer_rejects_unaccounted_settings(factors, axis):
    with pytest.raises(ValueError):
        m.pair_pointer_source(factors, axis)


@pytest.mark.parametrize('jump_cap,count_cap', ((2, 2), (64, 8)))
def test_positive_log_modulus_bounds_the_actual_normalized_sampler_polynomial(jump_cap, count_cap):
    import mpmath as mp
    weights = m.inverse_clock.subordinated_filter_weights(jump_cap=jump_cap, count_cap=count_cap)
    result = m.prepared_log_limit_bound(jump_cap=jump_cap, count_cap=count_cap)
    assert result['sampler_coefficients_checked_by_existing_positive_evaluator']
    with mp.workdps(80):
        number = lambda x: mp.mpf(x.numerator)/x.denominator
        mass0 = sum(number(w) for w in weights)
        beta = mp.mpf('.5')/mp.sqrt(60)
        assert number(result['infinite_reference_mass_lower']) <= mp.exp(-beta*mp.sqrt(8))
        for value in range(17):
            lam = mp.mpf(value)/2
            mass = sum(number(w)*(1-lam/120)**(2*n) for n, w in enumerate(weights))
            generator = mp.log(mass0/mass)
            assert abs(generator-beta*mp.sqrt(lam)) <= number(
                result['ideal_generator_to_b_sqrt_D_upper'])
            if not value:
                assert generator == 0
                assert -mp.log(mass) > 0  # Omitting m0 changes even the zero clock.


def test_log_modulus_covers_noncommuting_source_filter_error():
    from scipy.linalg import logm
    epsilon = F(1, 1000)
    result = m.prepared_log_limit_bound(filter_operator_error=epsilon)
    weights = m.inverse_clock.subordinated_filter_weights()
    ideal = np.diag(1-np.array([0., 2., 6., 8.])/120)
    perturbation = np.roll(np.eye(4), 1, axis=0)
    actual = ideal+float(epsilon)*perturbation
    reference = sum(float(w)*np.linalg.matrix_power(ideal, 2*n) for n, w in enumerate(weights))
    raw = sum(float(w)*(np.linalg.matrix_power(actual, n)
                        @np.linalg.matrix_power(actual, n).conj().T)
              for n, w in enumerate(weights))
    assert np.linalg.norm(reference@raw-raw@reference, 2) > 1e-7
    assert np.linalg.norm(raw-reference, 2) <= float(result['raw_positive_mixture_operator_error_upper'])
    assert np.linalg.norm(logm(raw)-logm(reference), 2) <= float(
        result['raw_positive_mixture_log_error_upper'])
    unconstrained = m.prepared_log_limit_bound(filter_operator_error=F(1))
    assert unconstrained['raw_positive_mixture_log_error_upper'] is None
    assert unconstrained['raw_generator_to_b_sqrt_D_upper'] is None


def test_fixed_caps_collapse_instead_of_approaching_a_nonzero_square_root_clock():
    errors = []
    for rate in (4, 16, 32):
        weights = m.inverse_clock.subordinated_filter_weights(rate=F(rate))
        heat = F(1, 16*rate*rate)  # b=R sqrt(h)=1/4 stays fixed.
        lam = F(4)
        mass = sum(w*(1-heat*lam/2)**(2*n) for n, w in enumerate(weights))
        gap = float(sp.log(sp.Rational(sum(weights))/sp.Rational(mass)))
        result = m.positive_log_clock_bound(F(rate), heat, lam, 1-sum(weights), 4)
        assert gap <= float(result['fixed_cap_generator_upper'])
        errors.append(abs(gap-.5))
    assert errors[-1] > errors[0] > .45
    assert errors[-1] > .499


def test_refinement_caps_precision_and_small_mass_obstruction_are_explicit():
    coarse = m.positive_log_refinement_schedule(3)
    assert coarse['raw_generator_to_b_sqrt_D_upper'] is None
    rows = [m.positive_log_refinement_schedule(level, source_preserves_band=True)
            for level in (4, 6, 8)]
    errors = [float(row['raw_generator_to_b_sqrt_D_upper']) for row in rows]
    assert errors[0] > errors[1] > errors[2] > 0
    assert errors[2] < .0012
    for row in rows:
        assert row['dyadic_rounding_upper'] <= row['jump_tail_upper']
        assert row['inverse_temperature_interval'] == (F(1), F(1))
        assert not row['growing_sampler_evaluated']
        assert not row['growing_source_filter_or_phase_compiler_certified']
        assert not row['retained_request_algebra_has_this_modular_generator']
    first = rows[0]
    fixed_precision = m.positive_log_clock_bound(first['rate'], first['heat_step'], F(8),
        first['missing_sampler_mass_upper'], first['maximum_filter_requests'], F(1, 10**9))
    assert fixed_precision['raw_generator_to_b_sqrt_D_upper'] is None
    with pytest.raises(ValueError, match='bounded band'):
        m.positive_log_clock_bound(F(1), F(1), F(2), F(0), 1)


def test_a_band_floor_does_not_bound_log_errors_from_coupled_small_outside_modes():
    # In exact arithmetic M=diag(1,exp(-700)) and Mraw=R(theta) M R(theta)*.
    # Their operator distance is (1-exp(-700))*sin(theta), but the compressed
    # log difference on the first coordinate is 700*sin(theta)^2. Its large
    # value is visible without diagonalizing an ill-conditioned float matrix.
    import math
    theta = .01
    mixture_error = (1-math.exp(-700))*math.sin(theta)
    bad_band_only_log_bound = mixture_error/(1-mixture_error)
    actual_compressed_log_error = 700*math.sin(theta)**2
    assert actual_compressed_log_error > 6*bad_band_only_log_bound
    unverified = m.positive_log_refinement_schedule(6)
    assert unverified['ideal_generator_to_b_sqrt_D_upper'] is not None
    assert unverified['raw_positive_mixture_operator_error_upper'] is not None
    assert unverified['raw_positive_mixture_log_error_upper'] is None
    assert unverified['raw_generator_to_b_sqrt_D_upper'] is None
    assert not unverified['off_band_source_log_bound_claimed']
    with pytest.raises(ValueError, match='Boolean'):
        m.positive_log_refinement_schedule(6, source_preserves_band=1)


def test_coefficientwise_truncation_can_leave_an_unsupported_outside_mode():
    # w0=0,w1=1/8 lies below the R=1 compound law (its n=1 mass is1/(2e)).
    # D=diag(0,2/h) makes F=diag(1,0): the low-band logarithm exists while
    # the retained full polynomial is rank deficient. A small rotation of
    # its support cannot justify a full raw logarithm from a low-band floor.
    heat = F(1, 10)
    retained = np.diag([1/8, 0.])
    assert np.linalg.matrix_rank(retained) == 1
    result = m.positive_log_clock_bound(F(1), heat, F(1, 100), F(7, 8), 1, F(1, 1000))
    assert result['retained_band_mass_lower'] > 0
    assert result['ideal_generator_to_b_sqrt_D_upper'] is not None
    assert result['raw_generator_to_b_sqrt_D_upper'] is None


def test_native_success_herald_is_the_whole_twelve_dimensional_log_domain():
    _, data = m.inverse_clock._prepared_edge_gibbs_data()
    raw = data['raw_filter']
    assert raw.shape == (32768, 12)
    successful, failed = raw[:12], raw[12:]
    assert successful.shape == (12, 12)
    assert failed.shape == (32756, 12)  # Padding12:16 remains in this outcome.
    reference = m.prepared_reference()
    ideal = np.array(reference['lazy_filter'], float)
    assert np.linalg.norm(successful-ideal, 2) < 1e-9
    weights = reference['weights']
    mixed = sum(float(w)*(np.linalg.matrix_power(successful, n)
                          @np.linalg.matrix_power(successful, n).conj().T)
                for n, w in enumerate(weights))
    bound = m.prepared_log_limit_bound(filter_operator_error=F(1, 10**9))
    assert np.linalg.eigvalsh(mixed)[0] >= float(bound['retained_band_mass_lower']
        -bound['raw_positive_mixture_operator_error_upper'])
    assert bound['successful_filter_domain_selected_by_herald_projection']
    assert bound['successful_filter_domain_dimension'] == 12

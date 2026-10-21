"""Exact Fourier coefficients, compressed oscillator polynomials and bounds.

The new modes are cos(3x), sin(3x), not a complete momentum shell. All
occupation powers are formed before compression. Time is model time.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import factorial

from ..scalar_observations.bounds import sqrt_upper


LOW_DIMENSIONS = (2, 4, 4, 2, 2, 2, 2)
DIMENSIONS = LOW_DIMENSIONS + (4, 4)
MODES = ('zero', 'cos_x', 'sin_x', 'cos_y', 'sin_y', 'cos_z', 'sin_z',
         'cos_3x', 'sin_3x')
WAVES = ((0, 0, 0), (1, 0, 0), (1, 0, 0), (0, 1, 0), (0, 1, 0),
         (0, 0, 1), (0, 0, 1), (3, 0, 0), (3, 0, 0))


def integer(value, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError('an original integer in the admitted range is required')
    return value


def rational(value, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('an exact nonnegative rational is required')
    return F(value)


@lru_cache(maxsize=2)
def quartic_coefficients(extended=True):
    """Multinomial coefficients of V integral phi**4, in exact rationals.

    Nonconstant real modes have Fourier coefficients 1/sqrt(2) and
    +/-i/sqrt(2). Constant Fourier coefficients are selected explicitly.
    This is an independent integral calculation, not the model's formula.
    """
    if type(extended) is not bool:
        raise ValueError('boolean band selector required')
    count = 9 if extended else 7
    result = {}
    for indices in combinations_with_replacement(range(count), 4):
        factors = []
        nonconstant = 0
        for index in indices:
            if index == 0:
                factors.append([((0, 0, 0), 1)])
            else:
                nonconstant += 1
                wave = WAVES[index]
                sine = index % 2 == 0
                factors.append([(tuple(sign*x for x in wave), (-sign*1j if sine else 1))
                                for sign in (-1, 1)])
        total = 0j
        for choices in product(*factors):
            if all(sum(row[0][axis] for row in choices) == 0 for axis in range(3)):
                coefficient = 1
                for _, number in choices:
                    coefficient *= number
                total += coefficient
        if total == 0:
            continue
        if total.imag or nonconstant % 2:
            raise ArithmeticError('unexpected nonrational Fourier coefficient')
        powers = tuple(indices.count(i) for i in range(count))
        multiplicity = factorial(4)
        for power in powers:
            multiplicity //= factorial(power)
        result[powers] = F(int(total.real)*multiplicity, 2**(nonconstant//2))
    return result


def declared_quartic():
    """The explicit nine-mode action, independently checked by Fourier sums."""
    result = {}

    def term(coefficient, powers):
        key = tuple(powers.get(i, 0) for i in range(9))
        result[key] = result.get(key, F(0)) + F(coefficient)

    for i in range(9):
        term(1 if i == 0 else F(3, 2), {i: 4})
    for i in range(9):
        for j in range(i+1, 9):
            same_pair = (i, j) in ((1, 2), (3, 4), (5, 6), (7, 8))
            term(3 if same_pair else 6, {i: 2, j: 2})
    term(2, {1: 3, 7: 1})
    term(-6, {1: 1, 2: 2, 7: 1})
    term(6, {1: 2, 2: 1, 8: 1})
    term(-2, {2: 3, 8: 1})
    return result


@lru_cache(maxsize=32)
def polynomial(dimension, power, momentum=False):
    """Full Q**power (Q=(a+a*)/2) or P**2 in Bargmann coordinates."""
    d, power = integer(dimension), integer(power, 0)
    if type(momentum) is not bool or (momentum and power != 2):
        raise ValueError('only the squared momentum polynomial is admitted')
    matrix = [[F(0) for _ in range(d)] for _ in range(d)]
    for column in range(d):
        entries = {column: F(1)}
        for _ in range(power):
            following = {}
            for n, value in entries.items():
                following[n+1] = following.get(n+1, F(0)) + value*(1 if momentum else F(1, 2))
                if n:
                    following[n-1] = following.get(n-1, F(0)) + value*n*(-1 if momentum else F(1, 2))
            entries = following
        for row, value in entries.items():
            if row < d:
                matrix[row][column] = -value if momentum else value
    return tuple(tuple(row) for row in matrix)


def row_norm(matrix):
    return max(sum((abs(value)*sqrt_upper(F(factorial(i), factorial(j)))
                    for j, value in enumerate(row)), F(0))
               for i, row in enumerate(matrix))


def local_norm(d, power):
    return row_norm(polynomial(d, power))


def residual_bound(coupling=1):
    """R**2 upper bound for Q H J on the declared 9-mode occupation code.

    Only high labels (2,0),(0,2),(2,2),(1,0),(0,1) occur. Their ranges
    are orthogonal. The first two blocks include free high-blank mismatch;
    the odd blocks come exclusively from the resonant quartic terms.
    """
    g = rational(coupling)
    squares = sum((local_norm(d, 2) for d in LOW_DIMENSIONS), F(0))
    even_without_root_two = F(3, 4)+g/32+g*squares/16
    x, y = LOW_DIMENSIONS[1:3]
    odd_z = g*(local_norm(x, 3)+3*local_norm(x, 1)*local_norm(y, 2))/24
    odd_w = g*(3*local_norm(x, 2)*local_norm(y, 1)+local_norm(y, 3))/24
    squared = 4*even_without_root_two**2+(g/64)**2+odd_z**2+odd_w**2
    return {'squared_residual_norm_upper': squared,
            'residual_norm_upper': sqrt_upper(squared),
            'even_high_blocks_without_sqrt_two_bound': even_without_root_two,
            'odd_high_block_bounds': (odd_z, odd_w),
            'two_high_excitation_block_exact': g/64,
            'matched_mass_squared_increment': g/4,
            'matched_scalar_increment': F(7, 2)+g/32}


def history_bound(total_time, coupling=1):
    time = rational(total_time)
    result = residual_bound(coupling)
    squared = result['squared_residual_norm_upper']*time*time
    return {**result, 'total_model_time': time,
            'reduced_complete_output_trace_distance_upper': min(F(1), squared),
            'high_retained_complete_output_trace_distance_upper': min(F(1), sqrt_upper(squared)),
            'high_blank_survival_failure_probability_upper': min(F(1), squared),
            'postselection_used': False,
            'all_low_input_states_and_untouched_references': True,
            'arbitrary_identical_low_only_adaptive_instruments': True,
            'unrestricted_high_occupation_or_uv_limit_included': False}


def grid_history_bound(grid, total_time):
    """Same nine-mode code: finite differences versus exact gradients.

    Quartic Fourier quadrature is exact for grid>=13. For harmonic k,
    0<=k**2-kappa_grid(k)<=k**4*pi**2/(3*grid**2), using pi<22/7.
    """
    grid = integer(grid, 13)
    time = rational(total_time)
    weighted = sum((sum(k*k for k in wave)**2*local_norm(d, 2)
                    for wave, d in zip(WAVES, DIMENSIONS, strict=True)), F(0))
    generator = F(484, 294*grid*grid)*weighted
    return {'grid_per_axis': grid,
            'same_nine_mode_hamiltonian_norm_difference_upper': generator,
            'complete_history_trace_distance_upper': min(F(1), time*generator),
            'total_model_time': time,
            'quartic_quadrature_exact': True,
            'number_of_independent_modes': 9,
            'grid_increases_bandwidth': False}


def product_gibbs_tail(beta=F(1, 4), coupling=1, shell=8):
    """Exact spectral-rank tail for the two radial modes, at any occupation.

    H >= H_harm with energies n1+2*n2+3/2. The first k*(k+1)
    harmonic eigenvalues have n1+2*n2 < 2*k. Min-max transfers their
    ordered tail to every full-polynomial compression and the limit.
    This is not an occupation-projection error or a finite-cutoff rate.
    """
    beta, coupling, k = rational(beta, positive=True), rational(coupling), integer(shell, 0)
    # exp(-beta) <= 1/(1+beta); discard exp(-3*beta/2) in upper bounds.
    q = 1/(1+beta)
    upper = 1/((1-q)*(1-q*q))
    tail = q**(2*k)*((k+1)-k*q*q)*upper
    # Jensen on the frequency-two product blank gives Z >= exp(-beta E0).
    # (1-x/n)**n <= exp(-x), with n >= 2*x, is positive and rational.
    blank_energy = F(13, 8)+5*coupling/128
    x = beta*blank_energy
    n = max(1, (2*x.numerator+x.denominator-1)//x.denominator)
    lower = (1-x/n)**n
    return {'inverse_temperature': beta, 'coupling': coupling,
            'harmonic_shell_exclusive': 2*k,
            'retained_spectral_rank': k*(k+1),
            'product_blank_energy': blank_energy,
            'partition_function_lower': lower,
            'partition_function_upper': upper,
            'unnormalized_spectral_tail_upper': tail,
            'normalized_spectral_tail_upper': min(F(1), tail/lower),
            'all_full_polynomial_occupation_compressions': True,
            'unrestricted_two_mode_Gibbs_state': True,
            'occupation_projection_or_record_error_certified': False,
            'spatial_bandwidth_removed': False}


def product_occupation_refinement():
    """Independent finite-matrix diagnostics for the fixed two-mode limit.

    The proof is in README; these are not certified infinite-tail errors
    or new growing-code source executions. Canonical references are kept.
    """
    import numpy as np
    from .construction import ProductBand, history

    rows, comparisons, previous = [], [], None
    for d in (4, 8, 12):
        lattice = ProductBand(occupation=d)
        energies, basis = np.linalg.eigh(lattice.hamiltonian.toarray())
        unnormalized = np.exp(-energies/4)
        partition = float(unnormalized.sum())
        weights = unnormalized/partition
        canonical = (basis*np.sqrt(weights))@basis.T
        laws, branches = history(lattice, source=False, initial=canonical,
                                  return_states=True)
        rows.append({'occupation_per_mode': d,
                     'partition_function_numerical': partition,
                     'first_six_energy_levels_numerical': [float(x) for x in energies[:6]],
                     'complete_two_read_law': laws,
                     'canonical_reference_norm_squared': float(np.linalg.norm(canonical)**2)})
        if previous is not None:
            old_d, old_canonical, old_branches = previous
            addresses = [i+d*j for j in range(old_d) for i in range(old_d)]

            def embed(matrix):
                result = np.zeros_like(canonical, dtype=complex)
                result[np.ix_(addresses, addresses)] = matrix
                return result

            overlap = abs(np.vdot(embed(old_canonical), canonical))
            preparation_distance = np.sqrt(max(0., 1-overlap*overlap))
            complete_distance = 0.
            for old, new in zip(old_branches, branches, strict=True):
                old = embed(old)
                mass = np.linalg.norm(old)**2+np.linalg.norm(new)**2
                discriminant = mass*mass-4*abs(np.vdot(old, new))**2
                complete_distance += np.sqrt(max(0., discriminant))/2
            comparisons.append({'occupations': [old_d, d],
                'canonical_purification_trace_distance_numerical': float(preparation_distance),
                'complete_record_work_reference_trace_distance_numerical': float(complete_distance)})
        previous = d, canonical, branches
    return {'scope': 'Fixed two radial modes; ideal finite-matrix occupation checks with canonical untouched reference.',
            'preparation_sign': 1, 'model_time_per_interval': F(1, 16),
            'dimensions': rows, 'successive_comparisons': comparisons,
            'rational_spectral_tail_bounds': [product_gibbs_tail(shell=k) for k in (8, 16, 32)],
            'unrestricted_two_mode_limit_proved_in_README': True,
            'numerical_comparisons_certify_infinite_occupation_error': False,
            'growing_occupation_source_execution': False,
            'spatial_ultraviolet_limit_or_physical_identification': False}


def product_detector_kms_control(interaction=None):
    """Gibbs/KMS on the actual detector effect and actual kick generator.

    The cross weights are signed, not transition probabilities. The source
    comparison reuses the complete +/-Q/4 histories; no A-generated kick,
    equilibrium noise instrument, or extra physical experiment is compiled.
    """
    import numpy as np
    from scipy.linalg import expm
    from .construction import (ProductBand, normalized, product_gibbs_preparation,
                               product_interaction_control)

    lattice, beta, epsilon = ProductBand(), F(1, 4), F(1, 4)
    h = lattice.hamiltonian.toarray()
    energies, basis = np.linalg.eigh(h)
    weights = np.exp(-float(beta)*energies); weights /= weights.sum()
    rho = (basis*weights)@basis.T
    q = np.kron(normalized(polynomial(4, 1)), np.eye(4))
    values, vectors = np.linalg.eigh(q)
    effect = (vectors*(-np.sin(values)))@vectors.T
    ae, qe = basis.T@effect@basis, basis.T@q@basis
    # W_AQ(t)=sum_ij p_i A_ij Q_ji exp[-i(E_j-E_i)t].
    gaps = energies[None, :]-energies[:, None]
    cross = ae*qe.T
    forward, reverse = weights[:, None]*cross, weights[None, :]*cross
    kms = float(np.max(np.abs(reverse-np.exp(-float(beta)*gaps)*forward)))
    upper = np.triu_indices(16, 1)
    factor = np.exp(-float(beta)*gaps[upper])
    symmetric = (forward+reverse)/2
    antisymmetric = forward-reverse
    fdt = float(np.max(np.abs((1-factor)*symmetric[upper]
                              -(1+factor)*antisymmetric[upper]/2)))

    raw, target, preparation = product_gibbs_preparation(lattice, beta)
    raw_rho = raw@raw.conj().T  # No numerical renormalization.
    e = preparation['purification_vector_error_guard']
    trace_guard = e*(2+e)
    q_bound = local_norm(4, 1)
    raw_trace_error = float(np.linalg.norm(raw_rho-rho, 'nuc'))
    free_energies, free_basis = np.linalg.eigh(ProductBand(0).hamiltonian.toarray())
    free_weights = np.exp(-float(beta)*free_energies); free_weights /= free_weights.sum()
    free_rho = (free_basis*free_weights)@free_basis.T
    rows = []
    for time in (F(0), F(1, 16), F(1, 8)):
        u = expm(1j*float(time)*h)
        evolved = u@effect@u.conj().T
        ordered = np.trace(rho@evolved@q)
        reversed_order = np.trace(rho@q@evolved)
        response = float((1j*(ordered-reversed_order)).real)
        spectral = np.sum(forward*np.exp(-1j*float(time)*gaps))
        raw_ordered = np.trace(raw_rho@evolved@q)
        raw_response = float((1j*np.trace(raw_rho@(evolved@q-q@evolved))).real)
        wrap = expm(float(beta)*h)@evolved@expm(-float(beta)*h)
        free_defect = abs(np.trace(free_rho@wrap@q)-np.trace(free_rho@q@evolved))
        rows.append({'model_time': time,
            'ordered_detector_kick_correlation': [float(ordered.real), float(ordered.imag)],
            'symmetrized_detector_kick_correlation': float(((ordered+reversed_order)/2).real),
            'Q_kick_linear_response': response,
            'spectral_reconstruction_error': float(abs(ordered-spectral)),
            'raw_preparation_ordered_correlation_error': float(abs(raw_ordered-ordered)),
            'raw_preparation_response_error': float(abs(raw_response-response)),
            'free_Gibbs_under_interacting_clock_KMS_defect': float(free_defect)})

    if interaction is None:
        interaction = product_interaction_control()
    first_means = [sum((1-2*row['bits'][0])*row['probability'] for row in law)
                   for law in interaction['source_histories_by_native_sign']]
    source_slope = float((first_means[1]-first_means[0])/(2*float(epsilon)))
    time = F(1, 16)
    evolved = expm(1j*float(time)*h)@effect@expm(-1j*float(time)*h)
    kicks = [expm(sign*1j*float(epsilon)*q) for sign in (-1, 1)]
    exact_means = [float(np.trace(rho@kick.conj().T@evolved@kick).real) for kick in kicks]
    exact_slope = (exact_means[1]-exact_means[0])/(2*float(epsilon))
    # Central difference: ||ad_Q^3(A(t))|| <= 8 ||Q||^3 ||A||.
    finite_kick_guard = F(4, 3)*epsilon**2*q_bound**3
    history_guard = (interaction['ideal_splitting_complete_output_bound']
                     +interaction['complete_source_numerical_guard'])
    source_slope_guard = 2*history_guard/epsilon

    # The admitted menu preserves the first-mode parity. Central reweighting
    # changes its state but not the restricted modular flow or KMS relation.
    parity = np.kron(np.eye(4), np.diag([1., -1., 1., -1.]))
    even = (np.eye(16)+parity)/2
    changed = (np.eye(16)+even)@rho
    changed /= np.trace(changed)
    def logarithm(state):
        masses, axes = np.linalg.eigh(state)
        return (axes*np.log(masses))@axes.conj().T
    s = -.25
    modular = expm(1j*s*logarithm(rho))
    changed_modular = expm(1j*s*logarithm(changed))
    restricted_flow_error = float(np.linalg.norm(modular@effect@modular.conj().T
        -changed_modular@effect@changed_modular.conj().T, 2))
    changed_wrap = expm(float(beta)*h)@evolved@expm(-float(beta)*h)
    changed_kms = float(abs(np.trace(changed@changed_wrap@q)-np.trace(changed@q@evolved)))
    population_difference = float(abs(np.trace((changed-rho)@even)))
    cross_noise_difference = float(abs(np.trace((changed-rho)@effect@q)))
    if (max(kms, fdt, restricted_flow_error, changed_kms) > 1e-11
            or raw_trace_error > float(trace_guard)
            or max(row['raw_preparation_ordered_correlation_error'] for row in rows) > float(trace_guard*q_bound)
            or max(row['raw_preparation_response_error'] for row in rows) > float(2*trace_guard*q_bound)
            or abs(source_slope-exact_slope) > float(source_slope_guard)
            or abs(exact_slope-rows[1]['Q_kick_linear_response']) > float(finite_kick_guard)
            or min(population_difference, cross_noise_difference) < 1e-5
            or rows[1]['free_Gibbs_under_interacting_clock_KMS_defect'] < 1e-3):
        raise ArithmeticError('the same-state detector/kick KMS contract failed')
    return {'scope': 'Finite interacting Gibbs reference, actual detector effect and Q-generated kick; no extra noise-law fit.',
        'inverse_temperature': beta, 'detector_effect': '-sin(Q_mode_1)',
        'kick_generator': 'Q_mode_1', 'source_kick_amplitude': epsilon,
        'energy_levels_numerical': [float(x) for x in energies],
        'Gibbs_weights_numerical': [float(x) for x in weights],
        'cross_spectral_KMS_residual': kms, 'cross_fluctuation_response_identity_residual': fdt,
        'cross_spectral_weights_are_probabilities': False,
        'raw_prepared_trace_mass': float(np.trace(raw_rho).real),
        'raw_preparation_full_trace_error_numerical': raw_trace_error,
        'raw_preparation_full_trace_error_upper': trace_guard,
        'ordered_correlation_preparation_error_upper': trace_guard*q_bound,
        'response_preparation_error_upper': 2*trace_guard*q_bound,
        'time_correlations': rows,
        'source_first_read_means_for_minus_plus_Q_kicks': [float(x) for x in first_means],
        'source_central_kick_slope': source_slope,
        'ideal_central_kick_slope': exact_slope,
        'source_central_slope_error_upper': source_slope_guard,
        'central_slope_to_linear_response_error_upper': finite_kick_guard,
        'central_reweighting_control': {
            'even_sector_population_difference': population_difference,
            'detector_kick_symmetrized_noise_difference': cross_noise_difference,
            'restricted_modular_action_residual': restricted_flow_error,
            'restricted_KMS_residual': changed_kms},
        'all_source_readouts_and_references_retained': True,
        'equilibrium_cross_correlation_detector_executed': False,
        'noise_law_independently_fitted': False,
        'Gibbs_temperature_or_central_weights_selected_by_repair': False,
        'physical_clock_or_laboratory_dictionary_identified': False}


def serializable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    return value

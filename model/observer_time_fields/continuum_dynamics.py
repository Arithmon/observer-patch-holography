"""One source continuum operator for native decay, coherent waves and readouts.

The limit follows from the source compactness theorem, not finite numerical
extrapolation. Wave synthesis processes retained history blocks; its parameter
is distinct from passive native waiting and from a calibrated laboratory clock.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from . import spatial_dynamics as dynamics
from . import wave_reference
from ..physical_space_time import local_continuum as spatial
from tools import strict_json

ROOT = dynamics.ROOT
SOURCE_PATHS = tuple(sorted(set((*dynamics.SOURCE_PATHS, *wave_reference.SOURCE_PATHS,
    *spatial.SOURCE_PATHS, 'model/observer_time_fields/continuum_dynamics.py'))))


def _root_upper(value):
    value = F(value)
    if value < 0:
        raise ValueError('a nonnegative exact square is required')
    grid = 1 << 32
    root = math.isqrt(value.numerator*grid*grid//value.denominator)
    root += int(F(root*root, grid*grid) < value)
    return F(root, grid)


def propagation_bounds(level, band=16, horizon=4, energy=1):
    """Errors on a unit vector with full scaled native energy at most energy.

    These compare finite native evolutions with the same finite Laplacian.
    The subsequent spatial subsequence has no claimed convergence rate.
    """
    level = dynamics._level(level)
    band, horizon, energy = map(dynamics._positive, (band, horizon, energy))
    scale = 4**level
    q, n = 630*scale**2, 1260*scale**2
    if band >= n:
        raise ValueError('the retained band must lie below the native attempt normalization')
    tail = 2*_root_upper(energy/band)
    wave = dynamics.wave_bounds(level, band, horizon)
    wave_low = horizon*_root_upper(band)*wave['relative_wave_frequency_error']
    wave_rounding = _root_upper(band)/(25*scale)
    heat_low = horizon*band**2/(2*(n-band))+band/n
    return {'level': level, 'scaled_native_energy_upper': energy,
        'spectral_band_upper': band, 'parameter_horizon': horizon,
        'high_spectral_tail_amplitude_upper': tail,
        'integer_attempt_heat_error_on_band': heat_low,
        'integer_coherent_wave_error_on_band': wave_low+wave_rounding,
        'integer_coherent_parameter_rounding_upper': F(1, 25*scale),
        'heat_error_on_energy_domain_upper': min(F(2), heat_low+tail),
        'wave_error_on_energy_domain_upper': min(F(2), wave_low+wave_rounding+tail),
        'joint_cap_modular_cosine_error_upper': min(F(2), horizon*26/scale),
        'sufficient_native_return_cap': wave_reference.refinement_budget(level)['sufficient_native_length_cap'],
        'spatial_subsequence_convergence_rate_supplied': False}


def fixed_cap_collapse(level, cap=512, band=16):
    """A fixed length cap freezes every bounded scaled spatial band.

    P_h(x)/P_h(1)>=x^h gives Theta_h(x)<=-h log(x). This is an
    exact scalar inequality for the actual positive native transfer.
    """
    level = dynamics._level(level)
    if type(cap) is not int or cap < 0:
        raise ValueError('an original nonnegative length cap is required')
    band = dynamics._positive(band)
    scale, q = 4**level, 630*16**level
    if band >= 2*q:
        raise ValueError('the band must keep the native transfer positive')
    return {'level': level, 'fixed_cap': cap, 'band_upper': band,
        'scaled_modular_frequency_upper': F(26*cap*scale)*band/(2*q-band),
        'fixed_cap_low_band_generator_limit': 0,
        'uncapped_wave_generator_limit': 'sqrt(A)',
        'fixed_cap_interchanges_with_spatial_limit': False}


def wave_path_bounds(max_degree=10, horizon=4):
    """Uniform C([0,T],H^-2) tail and increment bounds for synthesized waves."""
    if type(max_degree) is not int or max_degree < 0:
        raise ValueError('an original nonnegative harmonic cutoff is required')
    horizon = dynamics._positive(horizon)
    return {'harmonic_cutoff': max_degree, 'parameter_horizon': horizon,
        'expected_squared_supremum_tail_upper': 4*(spatial.GREEN_BOUND+horizon**2)/(max_degree+1)**2,
        'expected_squared_increment_upper_per_dt_squared': F(2),
        'expected_fourth_increment_upper_per_dt_fourth': F(12),
        'initial_covariance': 'G=A^-1',
        'two_time_covariance': 'G cos((t-u) sqrt(A))',
        'velocity_covariance': 'identity on mean-zero L2(mu)',
        'preparation': 'two disjoint normalized blocks of one native-history functional CLT',
        'time_is_declared_spectral_readout_parameter': True}


def finite_speed_bounds(level, lipschitz=1, separation=1, heat_parameter=1):
    """Weighted native heat bounds and the resulting continuum support cone.

    Distances use the unit-circumradius flat comparison polyhedron. The
    support theorem acts on Abar=A plus its constant zero mode, not a
    nonlocal mean-zero projection. It is a continuum theorem only.
    """
    level = dynamics._level(level)
    lipschitz, separation, heat_parameter = map(dynamics._positive,
        (lipschitz, separation, heat_parameter))
    scale = 4**level
    ratio = lipschitz**2/(12*scale)
    if ratio >= 1:
        raise ValueError('refine until the squared edge weight is below twelve')
    # Successive cosh terms after x^2/2 have ratio at most x^2/12.
    native_exponent = heat_parameter*lipschitz**2/(2*(1-ratio))
    return {'level': level, 'weight_Lipschitz_constant': lipschitz,
        'heat_parameter': heat_parameter, 'comparison_support_separation': separation,
        'native_weighted_heat_log_norm_upper': native_exponent,
        'continuum_weighted_heat_log_norm_upper': heat_parameter*lipschitz**2/2,
        'continuum_separated_heat_log_norm_upper': -separation**2/(2*heat_parameter),
        'continuum_wave_speed_squared_upper': F(1, 2),
        'continuum_zero_wave_matrix_element_when_time_squared_below': 2*separation**2,
        'constant_zero_mode_restored_before_support_localization': True,
        'exact_finite_level_wave_support_cone_asserted': False,
        'comparison_units_are_metres_and_seconds': False}


def native_weighted_form_control():
    """Exact weighted full-source identity and two locality countercontrols."""
    graph, lap, _, _, _ = dynamics.source_green.matrices(0)
    matrix = sp.Matrix(lap.toarray())
    weights = [sp.Rational(2)**((p//12) % 3) for p in range(graph.ports)]
    weighted = sp.Matrix(graph.ports, graph.ports,
        lambda i, j: matrix[i, j]*(weights[i]/weights[j]+weights[j]/weights[i])/2)
    extra = sp.zeros(graph.ports)
    penalty = [sp.S.Zero]*graph.ports
    for a, b in graph.seams:
        eta = (weights[a]/weights[b]+weights[b]/weights[a])/2-1
        if a//12 == b//12 and eta:
            raise ArithmeticError('cell weights vary inside the native carrier')
        for p in (a, b):
            extra[p, p] += eta
            penalty[p] += eta
        extra[a, b] -= eta; extra[b, a] -= eta
    maximum = max((weights[a]/weights[b]+weights[b]/weights[a])/2-1
                  for a, b in graph.inter_seams)
    if (weighted != matrix+extra-sp.diag(*penalty)
            or min(penalty) < 0 or max(penalty) > maximum):
        raise ArithmeticError('the native weighted accretivity decomposition failed')
    # Without cell-repeated weights the internal edges invalidate the
    # cross-seam-only penalty: choose an internal port with no cross seam.
    endpoints = {p for edge in graph.inter_seams for p in edge}
    port = next(p for p in range(graph.ports) if p not in endpoints)
    degree = len(graph.neighbours[port])
    invalid_constant_form = -sp.Rational(degree, 2)
    if invalid_constant_form >= 0:
        raise ArithmeticError('the internal-weight countercontrol disappeared')
    # The finite continuous cosine has an immediate t^2 term across a
    # native cross seam. Exact finite-level support is not the limit theorem.
    a, b = graph.inter_seams[0]
    t_squared_coefficient = -matrix[a, b]/2
    if t_squared_coefficient != sp.Rational(1, 2):
        raise ArithmeticError('the finite wave-tail countercontrol changed')
    return {'source_level': 0, 'source_ports': graph.ports, 'source_cross_seams': len(graph.inter_seams),
        'exact_weighted_form_decomposition': True,
        'maximum_cross_seam_cosh_minus_one': str(maximum),
        'maximum_port_penalty': str(max(penalty)),
        'nonnegative_edge_and_diagonal_remainder': True,
        'internal_weight_mutation_cross_seam_penalty': '0',
        'internal_weight_mutation_constant_quadratic_form': str(invalid_constant_form),
        'finite_cross_seam_wave_t_squared_coefficient': str(t_squared_coefficient),
        'finite_continuous_wave_has_instantaneous_small_tail': True}


def exact_wave_synthesis():
    """Exact source20cell covariance and conserved velocity identities."""
    _, _, _, _, _, schur = dynamics._base_blocks()
    cells = schur.rows
    constant = sp.ones(cells)/cells
    p0 = sp.eye(cells)-constant
    green = (schur+constant).inv(method='DM')-constant
    correction = sp.Rational(1, 2520)
    covariance = green-correction*p0
    velocity = schur*covariance
    if (velocity != p0-correction*schur or covariance*schur != velocity
            or green*schur != p0 or covariance*sp.ones(cells, 1) != sp.zeros(cells, 1)):
        raise ArithmeticError('native wave-synthesis covariance identities changed')
    # S<=L_cell/12<=I/2 proves positivity of C and 0<=B C<=I.
    _, pivots = covariance[:cells-1, :cells-1].LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(cells-1)):
        raise ArithmeticError('the actual native block covariance is not positive')
    return {'source_cells': cells, 'source_ports': 12*cells,
        'native_attempt_covariance_sha256': dynamics._digest(covariance),
        'native_Green_sha256': dynamics._digest(green),
        'velocity_covariance_sha256': dynamics._digest(velocity),
        'attempt_covariance_correction': F(1, 2520),
        'velocity_covariance_deficit_norm_upper': F(1, 5040),
        'constant_mode_exactly_removed': True,
        'two_history_blocks_have_independent_Gaussian_limits': True,
        'synthesized_wave_is_unprocessed_native_occupation_history': False}


def native_propagator_control():
    """Actual240port native Chebyshev execution and an independent spectral solve."""
    graph, lap, _, _, _ = dynamics.source_green.matrices(0)
    full = lap.toarray().astype(float)
    eigen, vectors = np.linalg.eigh(full)
    if eigen[0] < -1e-12 or np.linalg.norm(full@vectors-vectors*eigen) > 1e-10:
        raise ArithmeticError('native spectral diagnostic has an excessive residual')
    eigen = np.maximum(eigen, 0.)
    # A prescribed mean-zero cell contrast lifted with the actual port isometry.
    field = np.repeat(np.array([1., -1.]+[0.]*18)/np.sqrt(24.), 12)
    coefficients = vectors.T@field
    kernel = np.eye(graph.ports)-full/1260
    theta = 2*np.arcsin(np.sqrt(eigen/2520))
    previous, current = field, kernel@field
    rows = []
    for count in range(51):
        if count == 0:
            result = field
        elif count == 1:
            result = current
        else:
            previous, current = current, 2*kernel@current-previous
            result = current
        if count not in (0, 1, 2, 25, 50):
            continue
        oracle = vectors@(np.cos(count*theta)*coefficients)
        residual = float(np.linalg.norm(result-oracle))
        parameter = count/np.sqrt(630.)
        target = vectors@(np.cos(parameter*np.sqrt(eigen))*coefficients)
        error = float(np.linalg.norm(result-target))
        allowance = dynamics.wave_bounds(0, band=12, parameter=max(F(1), F(count, 25)))['cosine_propagator_error_on_band']
        if residual > 3e-12 or error > float(allowance)+3e-12:
            raise ArithmeticError('native Chebyshev wave does not meet its source bound')
        amplitude = float(field@result)
        rows.append({'coherent_steps': count,
            'coherent_parameter_diagnostic': parameter,
            'Ramsey_plus_diagnostic': (1+amplitude)/2,
            'Ramsey_minus_diagnostic': (1-amplitude)/2,
            'spectral_oracle_residual_diagnostic': residual,
            'same_native_L_wave_error_diagnostic': error,
            'same_native_L_wave_error_upper': str(allowance)})
    return {'source_ports': graph.ports, 'full_native_seams': len(graph.seams),
        'fixed_preparation': '(cell0 minus cell1)/sqrt(24) on native ports',
        'readouts': rows, 'spatial_limit_inferred_from_finite_rows': False}


def _permutation_moments(vertices, particles, seams, adjacent):
    """Exact normalized traces, from fixed subsets of permutation cycles."""
    if (any(type(v) is not int for v in (vertices, particles, seams, adjacent))
            or vertices < 4 or not 0 < particles < vertices
            or not 0 < seams <= vertices*(vertices-1)//2
            or not 0 <= adjacent <= seams*(seams-1)):
        raise ValueError('a valid finite graph/sector census is required')
    holes = vertices-particles
    chi = F(particles*holes, vertices*(vertices-1))
    psi = F(particles*(particles-1)*holes*(holes-1),
            vertices*(vertices-1)*(vertices-2)*(vertices-3))
    transposition = 1-2*chi
    three_cycle = 1-3*chi
    two_pairs = 1-4*chi+4*psi
    disjoint = seams*(seams-1)-adjacent
    variance = (seams*(1-transposition**2)
                +adjacent*(three_cycle-transposition**2)
                +disjoint*(two_pairs-transposition**2))/(4*seams**2)
    if variance < 0:
        raise ValueError('the supplied census cannot have these permutation moments')
    return {'vertices': vertices, 'particles': particles, 'seams': seams,
        'ordered_adjacent_edge_pairs': adjacent, 'ordered_disjoint_edge_pairs': disjoint,
        'susceptibility': chi, 'transfer_trace_mean': 1-chi,
        'transfer_trace_variance': variance,
        'swap_trace': transposition, 'three_cycle_trace': three_cycle,
        'disjoint_swap_trace': two_pairs}


@lru_cache(maxsize=2, typed=True)
def _native_trace_counts(particles):
    if type(particles) is not int or particles not in (2, 4):
        raise ValueError('the retained native oracle uses total two or four')
    states = tuple(sorted(tuple(int(i in occupied) for i in range(12))
                          for occupied in combinations(range(12), particles)))
    index = {state: i for i, state in enumerate(states)}
    edges, _ = dynamics.source_green.passage.recovered_carrier()
    counts = np.zeros((len(states), len(states)), dtype=np.int64)
    for column, before in enumerate(states):
        for edge in edges:
            for coin in (0, 1):
                after, _ = dynamics.native.native(before, edge, coin)
                counts[index[after], column] += 1
    if not np.array_equal(counts, counts.T) or np.any(counts.sum(axis=0) != 60):
        raise ArithmeticError('the actual native full-configuration transfer changed')
    # Return immutable rows so the cached exact certificate cannot be mutated.
    return states, tuple(tuple(map(int, row)) for row in counts), tuple(edges)


def native_reference_state_control(particles=4):
    """Direct native branches check permutation moments and a nontracial signal."""
    states, rows, edges = _native_trace_counts(particles)
    dimension, denominator = len(states), 2*len(edges)
    degree = Counter(v for edge in edges for v in edge)
    adjacent = sum(d*(d-1) for d in degree.values())
    exact = _permutation_moments(12, particles, len(edges), adjacent)
    mean = F(sum(rows[i][i] for i in range(dimension)), dimension*denominator)
    second = F(sum(rows[i][j]*rows[j][i] for i in range(dimension)
                   for j in range(dimension)), dimension*denominator**2)
    if mean != exact['transfer_trace_mean'] or second-mean**2 != exact['transfer_trace_variance']:
        raise ArithmeticError('native branch trace moments differ from permutation counting')
    index = {state: i for i, state in enumerate(states)}
    commutator_numerator = 0
    for i, state in enumerate(states):
        if state[0] == state[1]:
            continue
        swapped = (state[1], state[0], *state[2:])
        commutator_numerator += rows[index[swapped]][i]
    # A=n0, B=i[A,swap01], rho_1=(I+K)/Tr(I+K).
    signal = F(commutator_numerator, dimension*denominator)/(1+mean)
    expected = exact['susceptibility']/(len(edges)*(2-exact['susceptibility']))
    if signal != expected or signal <= 0:
        raise ArithmeticError('the finite nontracial commutator control disappeared')
    counts = np.asarray(rows, dtype=float)/denominator
    eigenvalues = np.linalg.eigvalsh(counts)
    if eigenvalues[0] < -1e-12 or eigenvalues[-1] > 1+1e-12:
        raise ArithmeticError('the native transfer positivity control failed')
    weights = np.exp(-np.arccos(np.clip(eigenvalues, 0, 1)))
    distance = float(np.sum(np.abs(weights/weights.sum()-1/dimension))/2)
    bound = 5*_root_upper(exact['transfer_trace_variance'])/exact['susceptibility']
    if distance > float(bound)+1e-12:
        raise ArithmeticError('the finite reference exceeds its analytic trace-distance bound')
    # Independent two-configuration coherent-state control inside the same sector.
    occupation = sp.diag(1, 0)
    swap = sp.Matrix([[0, 1], [1, 0]])
    conjugate = sp.I*(occupation*swap-swap*occupation)
    pure = sp.ones(2)/2
    pure_signal = -sp.I*sp.trace(pure*(occupation*conjugate-conjugate*occupation))
    if pure_signal != 1:
        raise ArithmeticError('the supplied nontracial-state countercontrol changed')
    return {**exact, 'configuration_dimension': dimension,
        'native_endpoint_branches': dimension*denominator,
        'native_and_permutation_trace_moments_exactly_equal': True,
        'ideal_reference_half_trace_distance_diagnostic': distance,
        'ideal_reference_half_trace_distance_upper': min(F(1), bound),
        'cap_one_commutator_expectation_minus_i': signal,
        'finite_reference_exactly_tracial': False,
        'different_supplied_two_configuration_pure_state_commutator': F(pure_signal)}


def source_trace_census(level):
    """Check the analytic degree census on actual finite source seams."""
    if type(level) is not int or level not in range(3):
        raise ValueError('retained source census levels are zero through two')
    graph = dynamics.source_green.source.federation(level)
    degrees = Counter(Counter(v for edge in graph.seams for v in edge).values())
    scale = 4**level
    if degrees != {5: 180*scale, 6: 60*scale}:
        raise ArithmeticError('the actual source degree census changed')
    moments = _permutation_moments(graph.ports, graph.ports//3,
        len(graph.seams), sum(n*d*(d-1) for d, n in degrees.items()))
    if moments['transfer_trace_variance'] != reference_state_bounds(level)['transfer_trace_variance']:
        raise ArithmeticError('the source census and all-level trace formula differ')
    return {'level': level, 'degree_five_ports': degrees[5], 'degree_six_ports': degrees[6],
        'source_census_matches_all_level_formula': True}


def reference_state_bounds(level, modular_parameter=1, terms_per_vertex=1, overlaps=1, norm=1):
    """Full-configuration state bounds, not low-density-sector normalization.

    The growing cap is h_l=64*(4^l)^4. Its state-tail allowance is charged
    separately. All modular bounds keep the dimensionless parameter u fixed.
    """
    level = dynamics._level(level)
    u, a, b, norm = map(dynamics._positive,
                        (modular_parameter, terms_per_vertex, overlaps, norm))
    scale, vertices = 4**level, 240*4**level
    variance = F(2*(365*vertices**2-1631*vertices+1560),
                 11907*(vertices-2)*(vertices-1)**2)
    if not 0 < variance <= F(1, 16*vertices):
        raise ArithmeticError('the all-level trace-variance inequality failed')
    epsilon, tail = F(3, 8*2**level), F(1, 2*scale**2)
    # binom(V,V/3)>=3^(V/3)>=2^(V/3), and exp(pi/2)<5.
    probability_exponent = (5*(vertices-1)-1).bit_length()-vertices//3
    return {'level': level, 'ports': vertices, 'particles': vertices//3,
        'configuration_dimension': f'binom({vertices},{vertices//3})',
        'transfer_trace_mean': 1-F(2*vertices, 9*(vertices-1)),
        'transfer_trace_variance': variance, 'variance_upper': F(1, 16*vertices),
        'ideal_reference_half_trace_distance_to_normalized_trace_upper': epsilon,
        'growing_cap_reference_half_trace_distance_to_normalized_trace_upper': min(F(1), epsilon+tail),
        'fixed_modular_parameter': u,
        'modular_observable_operator_norm_upper': 1,
        'ideal_modular_automorphism_normalized_trace_L2_error_upper': min(F(2), 4*u/F(15*2**level)),
        'growing_cap_modular_automorphism_normalized_trace_L2_error_upper': min(F(2), 4*u/F(15*2**level)+2*u/scale**2),
        'density_mode_dimension': vertices-1,
        'density_mode_reference_probability_log2_upper': probability_exponent,
        'bounded_overlap_fluctuation_commutator_norm_upper': 2*a*b*norm**2,
        'growing_cap_fluctuation_commutator_expectation_upper': 4*a*b*norm**2*(epsilon+tail),
        'reference_half_trace_bound_holds_on_full_configuration_algebra': True,
        'modular_operator_norm_triviality_claimed': False,
        'refinement_rescaled_modular_time_triviality_claimed': False,
        'every_noncommutative_or_nontracial_limit_excluded': False}


def run():
    rows = {
        'schema': 'oph.observer_time_fields.continuum_dynamics.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'exact_native_wave_synthesis': exact_wave_synthesis(),
        'actual_native_propagators': native_propagator_control(),
        'all_level_propagation_bounds': [propagation_bounds(l) for l in (0, 1, 2, 4, 10)],
        'fixed_cap_order_of_limits_control': [fixed_cap_collapse(l) for l in (0, 2, 4, 10)],
        'continuum_path_bounds': [wave_path_bounds(j) for j in (0, 10, 100)],
        'native_weighted_form_control': native_weighted_form_control(),
        'continuum_finite_speed_bounds': [finite_speed_bounds(l) for l in (0, 1, 2, 4, 10)],
        'full_reference_state_limit': {
            'native_trace_and_nontracial_controls': [native_reference_state_control(q) for q in (2, 4)],
            'actual_source_degree_census': [source_trace_census(l) for l in range(3)],
            'all_level_reference_bounds': [reference_state_bounds(l) for l in (0, 1, 2, 4, 10, 100)],
            'bounded_algebra_readout_moments_and_products_must_converge': True,
            'fixed_parameter_bounded_algebra_limit_state_is_tracial': True,
            'density_wave_reference_weight_vanishes_without_sector_conditioning': True,
            'bounded_overlap_fluctuation_expected_commutator_tends_to_zero': True,
            'quantum_central_limit_existence_asserted_for_this_refinement': False},
        'common_limit': {'space': 'mean-zero L2(mu) on the source S2',
            'same_subsequence': 'operator-norm limit G of embedded native Schur inverses',
            'positive_self_adjoint_generator': 'A=G^-1 with compact resolvent',
            'constant_mode_extension_is_conservative_strongly_local_Dirichlet': True,
            'constant_mode_extension_has_finite_speed_in_comparison_distance': True,
            'source_heat_trace_spectral_dimension': 2,
            'native_attempt_parameter': 'tau=n/(1260*16^level)',
            'native_two_point_transfer': 'exp(-tau A)',
            'coherent_step_parameter': 't=n/(sqrt(630)*4^level)',
            'native_coherent_compression': 'cos(t sqrt(A))',
            'return_reference_modular_phase_magnitudes': 'sqrt(A)',
            'synthesized_wave_field_topology': 'C([0,T],H^-2(S2))',
            'uniform_propagator_convergence_on_compact_parameter_intervals': True,
            'propagator_convergence_topology': 'strong for each fixed input, uniform on compact parameter intervals',
            'same_source_covariance_resolvent_heat_and_wave_operator': True,
            'whole_refinement_sequence_limit_or_rate_proved': False},
        'scope': {'same_native_source_and_ancestry_dictionary': True,
            'stationary_sector_uniform_IID_seams_and_fair_repair_required': True,
            'full_initial_occupation_and_current_readback_required': True,
            'coherent_preparation_reflection_and_Ramsey_services_required': True,
            'Gaussian_wave_uses_declared_record_processing': True,
            'modular_reference_is_a_native_coherently_prepared_state': False,
            'faithful_continuum_modular_algebra_state_constructed': False,
            'full_modular_unitary_or_intervention_implementation_proved': False,
            'comparison_wave_support_cone_proved': True,
            'physical_signal_speed_or_Lorentzian_metric_identified': False,
            'observed_one_plus_three_spacetime_identified': False,
            'physical_interacting_QFT_CFT_or_vacuum_state_identified': False,
            'atomic_clock_or_public_lab_data_identification_proved': False}}
    # Retained diagnostics are explicitly numerical; theorem budgets are rational.
    def encode(value):
        if type(value) is F:
            return str(value)
        raise TypeError('unexpected receipt value')
    return json.loads(json.dumps(rows, default=encode))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        from model.observer.scalar_response.experiment import compare_receipt
        compare_receipt(strict_json.load(args.check), result)
        print('Common native continuum dynamics, modular cap and wave-path controls replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2, sort_keys=True); stream.write('\n')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

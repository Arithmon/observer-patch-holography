"""One source continuum operator for native decay, coherent waves and readouts.

The limit follows from the source compactness theorem, not finite numerical
extrapolation. Wave synthesis processes retained history blocks; its parameter
is distinct from passive native waiting and from a calibrated laboratory clock.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
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


def run():
    rows = {
        'schema': 'oph.observer_time_fields.continuum_dynamics.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'exact_native_wave_synthesis': exact_wave_synthesis(),
        'actual_native_propagators': native_propagator_control(),
        'all_level_propagation_bounds': [propagation_bounds(l) for l in (0, 1, 2, 4, 10)],
        'fixed_cap_order_of_limits_control': [fixed_cap_collapse(l) for l in (0, 2, 4, 10)],
        'continuum_path_bounds': [wave_path_bounds(j) for j in (0, 10, 100)],
        'common_limit': {'space': 'mean-zero L2(mu) on the source S2',
            'same_subsequence': 'operator-norm limit G of embedded native Schur inverses',
            'positive_self_adjoint_generator': 'A=G^-1 with compact resolvent',
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
            'finite_signal_speed_or_Lorentzian_metric_proved': False,
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

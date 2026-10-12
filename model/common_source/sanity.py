"""Internal physics diagnostics of existing, explicitly supplied model laws.

The scalar grid is the public scalar-response family's actual weighted
Dirichlet operator. It is not a spatial extension of the interacting two-site
common source. Its oscillator frequencies are read from the supplied action,
not independently measured. The clock check uses the common source's actual
post-sharp state and excludes only one proposed identification of flows.
Neither calculation compares the construction with experimental data.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.linalg import expm

from ..observer.scalar_response.field import ScalarField
from . import experiment, operational, scalar_coupling


GRIDS = (13, 34, 89)
MODE_LABELS = ((1, 1, 1), (2, 1, 1), (2, 2, 1), (3, 2, 1))


def scalar_mode_dispersion():
    """Compare the existing weighted scalar action with its declared KG target.

    The box is [0,1]^3 with zero Dirichlet boundaries. Length, mass, speed
    and hbar are set to one by the inherited dimensionless convention.
    Labels denote standing waves, not measured signed momentum eigenvalues.
    Only one-axis matrices are diagonalized; tensor separability gives the
    three-axis oscillator frequencies without constructing an alternative grid.
    """
    levels = []
    for q in GRIDS:
        field = ScalarField(q)
        diagonal = field.stiffness_diagonal / field.weights
        off = field.stiffness_off_diagonal / np.sqrt(field.weights[:-1] * field.weights[1:])
        laplacian = np.diag(diagonal) + np.diag(off, 1) + np.diag(off, -1)
        values, vectors = np.linalg.eigh(laplacian)
        operator_scale = float(np.linalg.norm(laplacian, 2))
        spectral_error = float(np.max(np.abs(values-field.laplacian_eigenvalues)))
        residual = float(np.linalg.norm(laplacian@vectors-vectors*values, 2) / operator_scale)
        rows = []
        for label in MODE_LABELS:
            eigenvalues = [float(values[n-1]) for n in label]
            lattice_k_squared = math.fsum(eigenvalues)
            k_squared = math.pi**2 * sum(n*n for n in label)
            frequency = math.sqrt(1+lattice_k_squared)
            target = math.sqrt(1+k_squared)
            rows.append({
                'standing_wave_label': list(label),
                'continuum_wavenumber_magnitude': math.sqrt(k_squared),
                'continuum_wavenumber_squared': k_squared,
                'actual_lattice_wavenumber_squared': lattice_k_squared,
                'actual_oscillator_frequency': frequency,
                'declared_continuum_KG_frequency': target,
                'absolute_frequency_error': abs(frequency-target),
                'relative_frequency_error': abs(frequency-target)/target,
                'frequency_squared_minus_continuum_wavenumber_squared': frequency**2-k_squared,
                'actual_frequency_over_wavenumber': frequency/math.sqrt(k_squared),
                'continuum_frequency_over_wavenumber': target/math.sqrt(k_squared),
                'finite_action_identity_residual': abs(frequency**2-1-lattice_k_squared),
            })
        levels.append({
            'q': q, 'one_axis_interior_sites': q-1,
            'tensor_oscillator_count': field.mode_count,
            'minimum_gap': float(np.min(field.gaps)),
            'maximum_gap': float(np.max(field.gaps)),
            'nonuniform_grid': bool(np.ptp(field.gaps) > 1e-12),
            'minimum_laplacian_eigenvalue': float(values[0]),
            'dense_vs_existing_tridiagonal_eigenvalue_error': spectral_error,
            'relative_eigenvector_residual': residual,
            'modes': rows,
            'maximum_selected_relative_frequency_error': max(row['relative_frequency_error'] for row in rows),
        })
    checks = {
        'actual_weighted_operator_is_positive': all(row['minimum_laplacian_eigenvalue'] > 0 for row in levels),
        'dense_solver_matches_existing_field_spectrum': all(row['dense_vs_existing_tridiagonal_eigenvalue_error'] < 2e-8 for row in levels),
        'eigenvector_residuals_small': all(row['relative_eigenvector_residual'] < 2e-13 for row in levels),
        'finite_action_identity': all(item['finite_action_identity_residual'] < 2e-12 for row in levels for item in row['modes']),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return {
        'source_operator': 'model/observer/scalar_response/field.py: ScalarField',
        'same_operator_used_by': 'model/observer/scalar_response/local_evolution.py: LocalEvolution.apply_action',
        'geometry': 'declared golden tensor grid in the unit three-dimensional box',
        'boundary_conditions': 'zero Dirichlet on all six faces',
        'declared_units': {'box_length': 1., 'mass': 1., 'propagation_speed': 1., 'hbar': 1.},
        'target_law': 'omega^2 = 1 + pi^2 (n_x^2+n_y^2+n_z^2)',
        'linear_dispersion_scope': 'omega approaches |k| only when |k| is large compared with mass and the mode remains resolved by the grid',
        'energy_readout': 'oscillator level spacing hbar*omega inferred from the supplied quadratic Hamiltonian',
        'levels': levels, 'checks': checks, 'all_pass': True,
        'scope': {
            'translation_momentum_operator_or_detector_constructed': False,
            'independent_detector_frequency_measurement': False,
            'Lorentz_invariance_derived': False,
            'common_interacting_source_spatial_extension': False,
            'comparison_with_experimental_data': False,
            'outward_numerical_error_certificate': False,
        },
    }


def clock_flow_comparison(d, steps=512):
    """Test the proposed equality on this actual complete post-sharp state.

    Every modular flow of a faithful state preserves that state, including
    after any scalar parameter rescaling. The Hamiltonian's nonzero local
    state derivative therefore obstructs that equality on all local
    expectations. This says nothing against other states, observer clocks,
    effective clock readings, or an independently calibrated time dictionary.
    """
    if type(d) is not int or d not in (2, 4):
        raise ValueError('executed scalar dimension 2 or 4 required')
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer split count required')
    history = experiment.history(d, steps)
    existing = operational.observation_receipt(history)
    vectors = operational.terminal_vectors(history)
    work = sum(np.outer(vector, vector.conj()) for vector in vectors.values())
    mass = float(np.trace(work).real)
    work = work / mass
    rho = operational.scalar_reduction(d, work)
    h = scalar_coupling.ideal_hamiltonian(d, history['coupling'])
    rate = operational.scalar_reduction(d, -1j*(h@work-work@h))
    rate = (rate+rate.conj().T)/2
    norm = float(np.linalg.norm(rate))
    spectral_norm = float(np.linalg.norm(rate, 2))
    if spectral_norm <= 1e-3:
        raise AssertionError('this retained state no longer witnesses the stated flow obstruction')
    witness = rate / spectral_norm
    expectation_rate = float(np.trace(rate@witness).real)
    modular_generator = operational.modular_hamiltonian(rho)
    modular_rate = -1j*(modular_generator@rho-rho@modular_generator)
    moved = operational.modular_flow(rho, witness, .7)
    modular_expectation_error = float(abs(np.trace(rho@moved)-np.trace(rho@witness)))
    delta = 1e-6
    forward, backward = expm(-1j*delta*h), expm(1j*delta*h)
    central_rate = operational.scalar_reduction(d,
        (forward@work@forward.conj().T-backward@work@backward.conj().T)/(2*delta))
    numerical_derivative_error = float(np.linalg.norm(central_rate-rate))
    checks = {
        'same_retained_history_derivative': abs(norm-existing['source_model_time_local_state_derivative_frobenius']/mass) < 2e-12,
        'finite_difference_matches_commutator': numerical_derivative_error < 2e-8,
        'modular_reference_state_is_stationary': float(np.linalg.norm(modular_rate)) < 2e-12,
        'modular_witness_expectation_is_stationary': modular_expectation_error < 2e-12,
        'actual_model_time_witness_changes': expectation_rate > 1e-3,
        'witness_is_norm_one_local_observable': abs(float(np.linalg.norm(witness, 2))-1) < 2e-12,
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return {
        'dimension': d, 'steps_per_segment': steps,
        'state': 'complete post-sharp history, forgetting records, restricted to scalar site zero',
        'pre_normalization_mass': mass,
        'state_eigenvalues': np.linalg.eigvalsh(rho).tolist(),
        'model_time_local_state_derivative_frobenius': norm,
        'modular_state_derivative_frobenius': float(np.linalg.norm(modular_rate)),
        'norm_one_witness': [[[float(z.real), float(z.imag)] for z in row] for row in witness],
        'model_time_witness_expectation_derivative': expectation_rate,
        'modular_witness_expectation_change_at_parameter_0_7': modular_expectation_error,
        'central_difference_model_time_increment': delta,
        'central_difference_derivative_error': numerical_derivative_error,
        'same_state_model_flow_equals_scalar_rescaled_modular_flow': False,
        'obstruction': 'the modular reference state is stationary while this local state changes under the actual common Hamiltonian',
        'checks': checks, 'all_pass': True,
        'scope': {
            'every_possible_observer_clock_rejected': False,
            'programmed_modular_clock_invalidated': False,
            'physical_clock_calibration_established': False,
            'rigorous_roundoff_enclosure_supplied': False,
            'OPH_empirically_falsified': False,
        },
    }


def report(steps=512):
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer split count required')
    dispersion = scalar_mode_dispersion()
    clocks = [clock_flow_comparison(d, steps) for d in (2, 4)]
    return {
        'schema': 'oph.common_source.internal_physics_sanity.v1',
        'scalar_mode_dispersion': dispersion, 'clock_flow_comparison': clocks,
        'all_pass': dispersion['all_pass'] and all(row['all_pass'] for row in clocks),
        'scientific_verdict': {
            'agreement_with_nature_established': False,
            'specified_model_time_modular_time_identification_rejected_on_tested_state': True,
            'massive_dispersion_is_a_check_of_supplied_action': True,
            'whole_theory_verdict': 'not supplied by these internal finite-model checks',
        },
    }


if __name__ == '__main__':
    print(json.dumps(report(), indent=2, sort_keys=True, allow_nan=False))

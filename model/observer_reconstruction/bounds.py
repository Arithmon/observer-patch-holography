"""Directed cutoff certificate for vacuum Weyl/rotation observer histories.

This concerns the ideal number-basis compression of each normal mode.  It
does not bound floating source execution, physical services, finite sampling
or error amplification in the inverse reconstruction.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import factorial

from ..scalar_observations.bounds import sqrt_upper


def coherent_history_cutoff(dimension, *, amplitude_squared_upper, modes=8,
                           digits=36):
    """Bound the full controlled-displacement history, including its field.

    Every mode starts in its vacuum.  Its history consists of displacements
    D(delta) and number rotations, with the sum of absolute displacement
    amplitudes in every pointer arm at most A.  ``amplitude_squared_upper``
    is an exact rational upper bound on A**2, common to all ``modes``.
    A common coherent pointer may control these displacements in one fixed
    basis, without intermediate rotations mixing its arms; both its arms
    and both eventual recorded outcomes are retained.

    With P retaining occupations 0,...,d-1, projected Duhamel gives
    ||finite - P infinite|| <= sqrt(d)*A**(d+1)/sqrt(d!).  Its only boundary
    residual is the lowering of the infinite n=d coefficient, whose
    modulus is at most A**d/sqrt(d!).  Rotations commute with P.  The
    final coherent tail is at most
    A**(2*d)/(d!*(1-A**2/(d+1))); low and high components are orthogonal.
    Thus L**2+T bounds the squared full single-mode vector error.  The
    mode isometries are diagonal in the common pointer's arms, so their
    error bounds remain valid with that pointer entangled with earlier
    modes.  Telescoping sums their norm errors.  Pure-state trace distance
    and the final recorded-output trace distance are no larger.

    The supplied amplitude bound must hold for the actual declared word;
    this function does not infer it from floating eigenvectors or angles.
    It neither normalizes a projected state nor discards a pointer arm.
    """
    if type(dimension) is not int or dimension < 2:
        raise ValueError("occupation dimension must be an integer of at least two")
    if type(amplitude_squared_upper) is not F or amplitude_squared_upper < 0:
        raise ValueError("squared amplitude bound must be an exact nonnegative Fraction")
    if amplitude_squared_upper >= dimension + 1:
        raise ValueError("squared amplitude bound must be smaller than dimension + 1")
    if type(modes) is not int or modes < 1:
        raise ValueError("number of modes must be a positive integer")
    if type(digits) is not int or digits < 1:
        raise ValueError("directed root precision must be a positive integer")

    a2 = amplitude_squared_upper
    boundary_weight = a2**dimension / factorial(dimension)
    projected_squared = dimension * a2 * boundary_weight
    tail = min(F(1), boundary_weight / (1-a2/F(dimension+1)))
    single_squared = projected_squared + tail
    single_vector = sqrt_upper(single_squared, digits=digits)
    joint_vector = modes * single_vector
    trace = min(F(1), joint_vector)
    return {
        "schema": "oph.observer_reconstruction.coherent_history_cutoff.v1",
        "occupation_dimension": dimension,
        "number_of_modes": modes,
        "squared_displacement_budget_upper": a2,
        "directed_root_decimal_digits": digits,
        "projected_path_error_squared_upper": projected_squared,
        "final_single_mode_tail_probability_upper": tail,
        "single_mode_vector_error_squared_upper": single_squared,
        "single_mode_vector_error_upper": single_vector,
        "joint_vector_error_upper": joint_vector,
        "complete_output_trace_distance_upper": trace,
        "pointer_probability_error_upper": trace,
        "pointer_mean_error_upper": 2*trace,
        "scope": {
            "vacuum_preparation": True,
            "all_controlled_displacements_in_budget": True,
            "pointer_control_basis_fixed_until_terminal_readout": True,
            "number_rotations_preserve_cutoff": True,
            "both_pointer_arms_and_terminal_field_included": True,
            "both_pointer_outcomes_retained": True,
            "rational_directed_bound_arithmetic": True,
            "arbitrary_initial_field_states_certified": False,
            "numerical_source_execution_error_included": False,
            "physical_service_error_included": False,
            "finite_sampling_error_included": False,
            "inverse_reconstruction_error_included": False,
            "spatial_continuum_limit_certified": False,
        },
    }

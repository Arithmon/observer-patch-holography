"""Finite rest-pair pole diagnosis for the retained quadratic candidate.

This is an independent four-state CAR calculation, not a replay or change
of the retained native evolution. It distinguishes the blank-pair quench
frequency from the charged linear-field pole, and tests reference-state
independence without selecting a physical vacuum.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from scipy.linalg import expm
import sympy as s

from ..chiral_candidate import action, history
from ..mass_response import model


def _zero(matrix):
    return all(s.simplify(value) == 0 for value in matrix)


def _ladders():
    # Complete two-low-mode Fock block of the existing CAR sign implementation.
    return tuple(np.asarray(a[:4, :4], complex) for a in model.annihilators()[:2])


def _exact_ladders():
    result = []
    for array in _ladders():
        if np.max(abs(array.imag)) != 0 or not np.array_equal(array.real, np.rint(array.real)):
            raise ArithmeticError("the retained CAR block lost its exact integer entries")
        result.append(s.ImmutableMatrix(np.rint(array.real).astype(int).tolist()))
    return tuple(result)


@lru_cache(maxsize=1)
def exact_identities():
    """Symbolic identities for arbitrary complex pairing and positive rest gap."""
    a, b = _exact_ladders()
    x, y = s.symbols("d_real d_imag", real=True)
    d, magnitude_squared = x+s.I*y, x*x+y*y
    h = d*a.H*b.H+s.conjugate(d)*b*a
    delta = s.Matrix([[0, d], [-d, 0]])
    nambu = s.BlockMatrix([[s.zeros(2), delta], [delta.H, s.zeros(2)]]).as_explicit()
    xi = (a, b, a.H, b.H)
    for i, left in enumerate(xi):
        for j, right in enumerate(xi):
            if not _zero(left*right.H+right.H*left-(s.eye(4) if i == j else s.zeros(4))):
                raise ArithmeticError("full equal-time Nambu CAR identity failed")
        if not _zero(left*h-h*left-sum((nambu[i, j]*value for j, value in enumerate(xi)), s.zeros(4))):
            raise ArithmeticError("exact Nambu Heisenberg equation failed")
    characteristic = h.charpoly()
    energy = characteristic.gen
    if s.simplify(characteristic.as_expr()-energy**2*(energy**2-magnitude_squared)) != 0:
        raise ArithmeticError("finite Fock rest-pair spectrum failed")
    if not _zero(nambu*nambu-magnitude_squared*s.eye(4)):
        raise ArithmeticError("quadratic Nambu squared-gap identity failed")
    m, z, time = s.symbols("m z time", positive=True)
    positive_h = m*(a.H*b.H+b*a)
    ground = s.Matrix([1, 0, 0, -1])/s.sqrt(2)
    odd = s.diag(0, 1, 1, 0)
    if not _zero(positive_h*ground+m*ground):
        raise ArithmeticError("illustrative pair ground state identity failed")
    residues = [(vector.H*odd*vector)[0] for vector in (a.H*ground, a*ground)]
    if residues != [s.Rational(1, 2), s.Rational(1, 2)]:
        raise ArithmeticError("exact charged ground CAR residues failed")
    if s.simplify(z/(z*z-m*m)-(1/(z-m)+1/(z+m))/2) != 0:
        raise ArithmeticError("ordinary CAR resolvent residue identity failed")
    if s.trigsimp(s.cos(m*time)**2-s.sin(m*time)**2-s.cos(2*m*time)) != 0:
        raise ArithmeticError("blank-pair occupation frequency identity failed")
    return {
        "Hamiltonian": "H=d*a^dagger*b^dagger+conjugate(d)*b*a",
        "Fock_basis": ["00", "10", "01", "11"],
        "characteristic_polynomial": "E^2*(E^2-|d|^2)",
        "Fock_energies": ["-|d|", "0", "0", "+|d|"],
        "ground_to_odd_CAR_excitation_cost": "|d|",
        "ground_to_even_pair_excitation_cost": "2*|d|",
        "Nambu_order": ["a", "b", "a^dagger", "b^dagger"],
        "equal_time_CAR_exact": True,
        "Heisenberg_commutator_identity_exact": "[xi,H]=B*xi",
        "Nambu_square_exact": "B^2=|d|^2*I4",
        "reference_independent_operator_identity": "{xi_i(t),xi_j^dagger}=[exp(-i*B*t)]_ij*I_Fock",
        "reference_independence_proof": "Finite linear Heisenberg evolution and equal-time CAR; the anticommutator is a scalar operator before taking any state expectation.",
        "retarded_resolvent": "G_R(z)=(z*I4+B)/(z^2-|d|^2), Im(z)>0",
        "normal_local_resolvent": "G_aa(z)=z/(z^2-|d|^2)=(1/2)/(z-|d|)+(1/2)/(z+|d|)",
        "normal_local_signed_poles": ["-|d|", "+|d|"],
        "normal_local_residues_exact": ["1/2", "1/2"],
        "ground_addition_and_removal_weights_exact": ["1/2", "1/2"],
        "complex_phase_scope": "An onsite phase changes arg(d) without changing |d| or either normal residue.",
        "blank_pair_occupation_signal": "P_empty(t)-P_occupied(t)=cos(2*|d|*t)",
        "blank_quench_frequency_is_twice_CAR_gap": True,
    }


def pair_matrices(coefficient):
    """Actual complete four-state CAR H and its independently constructed BdG B."""
    if isinstance(coefficient, (bool, np.bool_)) or isinstance(coefficient, (str, bytes)) or np.ndim(coefficient) != 0:
        raise ValueError("a finite nonzero complex pairing coefficient required")
    try:
        d = complex(coefficient)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("a finite nonzero complex pairing coefficient required") from exc
    if not math.isfinite(d.real) or not math.isfinite(d.imag) or not 0 < abs(d) < math.inf:
        raise ValueError("a finite nonzero complex pairing coefficient required")
    a, b = _ladders()
    h = d*a.conj().T@b.conj().T+d.conjugate()*b@a
    delta = np.array([[0, d], [-d, 0]])
    nambu = np.block([[np.zeros((2, 2)), delta], [delta.conj().T, np.zeros((2, 2))]])
    return h, nambu, (a, b, a.conj().T, b.conj().T)


def state_independence_check(coefficient, time=.37):
    """Compare the full operator anticommutator and four distinct references."""
    if isinstance(time, (bool, np.bool_)) or isinstance(time, (str, bytes)) or np.ndim(time) != 0:
        raise ValueError("a finite real model parameter required")
    try:
        time = float(time)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("a finite real model parameter required") from exc
    if not math.isfinite(time):
        raise ValueError("a finite real model parameter required")
    h, nambu, xi = pair_matrices(coefficient)
    unitary = expm(-1j*time*h)
    evolved = tuple(unitary.conj().T@value@unitary for value in xi)
    anticommutator = np.array([[left@right.conj().T+right.conj().T@left for right in xi]
                              for left in evolved])
    prediction = expm(-1j*time*nambu)
    operator_error = float(np.max(abs(anticommutator-prediction[:, :, None, None]*np.eye(4))))
    _, basis = np.linalg.eigh(h)
    ground = basis[:, 0]
    references = {
        "encoded_blank": np.diag([1., 0., 0., 0.]),
        "illustrative_pair_ground": np.outer(ground, ground.conj()),
        "maximally_mixed": np.eye(4)/4,
        "nonGaussian_even_pair_mixture": np.diag([.5, 0., 0., .5]),
    }
    expectation_errors, reference_stationarity, greater = {}, {}, {}
    for name, density in references.items():
        observed = np.einsum("ij,abji->ab", density, anticommutator)
        expectation_errors[name] = float(np.max(abs(observed-prediction)))
        reference_stationarity[name] = float(np.max(abs(h@density-density@h)))
        value = np.trace(density@evolved[0]@xi[0].conj().T)
        greater[name] = [float(value.real), float(value.imag)]
    return {"model_parameter": time,
            "operator_anticommutator_residual": operator_error,
            "reference_expectation_residuals": expectation_errors,
            "reference_stationarity_residuals": reference_stationarity,
            "normal_greater_correlations_do_depend_on_reference": greater,
            "all_references_are_physical_vacuum_candidates": False,
            "Gaussian_reference_required_for_operator_identity": False,
            "independent_Fock_matrix_exponential_evaluated": True,
            "native_source_history_replayed": False,
            "floating_residuals_outward_certified": False}


def pair_report(coefficient, *, time=.37):
    h, nambu, xi = pair_matrices(coefficient)
    energies, basis = np.linalg.eigh(h)
    ground, odd = basis[:, 0], np.diag([0., 1., 1., 0.])
    magnitude = abs(complex(coefficient))
    create, destroy = xi[0].conj().T@ground, xi[0]@ground
    addition = float(np.vdot(create, odd@create).real)
    removal = float(np.vdot(destroy, odd@destroy).real)
    unitary = expm(-1j*time*h)
    blank_output = unitary[:, 0]
    occupation = float(abs(blank_output[1])**2+abs(blank_output[3])**2)
    return {"pairing_coefficient": [complex(coefficient).real, complex(coefficient).imag],
            "charged_CAR_gap": magnitude, "Fock_energy_eigenvalues": energies.tolist(),
            "ground_to_odd_gap": float(-energies[0]),
            "ground_to_even_pair_gap": float(energies[-1]-energies[0]),
            "normal_signed_poles": [-magnitude, magnitude],
            "normal_pole_residues": [removal, addition],
            "maximum_normal_residue_error_from_one_half": max(abs(addition-.5), abs(removal-.5)),
            "blank_pair_occupation_cosine_residual": abs((1-2*occupation)-math.cos(2*magnitude*time)),
            "Nambu_squared_gap_residual": float(np.max(abs(nambu@nambu-magnitude**2*np.eye(4)))),
            "reference_independence": state_independence_check(coefficient, time)}


def report():
    """JSON-safe diagnosis using the actual three retained lepton pairings."""
    exact = exact_identities()
    _, layout = history.rest_layout()
    pairing = action.pairing_kernel()
    labels = action.mode_labels()
    rows = []
    for i, j in layout["probed_original_mode_pairs"]:
        coefficient = s.simplify(pairing[i, j])
        row = pair_report(complex(coefficient.evalf()))
        row.update({"native_family_axis": labels[i]["family"], "original_mode_pair": [i, j],
                    "pairing_coefficient_exact": str(coefficient),
                    "charged_CAR_gap_exact": str(s.simplify(s.sqrt(coefficient*s.conjugate(coefficient))))})
        rows.append(row)
    return {"schema": "oph.hierarchy_diagnosis.quadratic_response.v1",
            "exact_finite_CAR_identities": exact, "actual_lepton_pair_rows": rows,
            "interpretation": {
                "same_quadratic_generator_reference_change_moves_linear_retarded_poles": False,
                "constant_linear_probe_can_hide_poles_by_zero_overlap": True,
                "constant_linear_probe_moves_remaining_poles": False,
                "factor_two_changes_family_ratios": False,
                "normal_residue_anisotropy_explains_the_trial_hierarchy": False,
                "time_dependent_background_or_selfenergy_changes_covered": False,
                "nonlinear_composite_field_poles_classified": False,
                "physical_interacting_vacuum_selected": False,
                "observed_particles_or_masses_identified": False,
                "each_pair_is_one_rest_block_not_a_new_complete_Dirac_species": True,
            },
            "source_parameters_or_mass_maps_changed": False,
            "empirical_mass_inputs": [], "native_source_evolution_replayed": False}

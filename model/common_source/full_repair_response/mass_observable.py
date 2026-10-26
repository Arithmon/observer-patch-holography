"""Calibrate a charged mass observable on the actual declared chiral action.

The probe is a charged CAR field with nonzero spectral residue, observed at
several momenta against the existing source clock. A literal four-state rest
pair independently connects its charged pole, ground-state Euclidean decay
and empty-pair quench frequency. Neutral native CP relaxation is a different
operator response until an actual field/history intertwiner is supplied.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from scipy.linalg import expm
import sympy as s

from ..chiral_candidate import action, gaussian_source, history, kinetic


INTERVAL = F(1, 4)
MOMENTA = ((.17, 0., 0.), (0., .29, 0.), (0., 0., .43))


def _numeric(value):
    return np.asarray(value.evalf(17) if hasattr(value, "evalf") else value, complex)


def _complex(value):
    return [float(complex(value).real), float(complex(value).imag)]


def _clock_record(duration=INTERVAL):
    """Actual source preparation/evolution, with both clock purifiers retained."""
    initial, calibration, preparation_cost = history._clock_preparation()
    final = history._clock_advance(initial, duration, backend="source")
    rho = final@final.conj().T
    x, y, _ = kinetic.pauli_matrices()
    # Clock bit0 is the low bit in the native source's little-endian encoding.
    readings = [float(np.trace(rho@np.kron(np.eye(2), p)).real) for p in (x, y)]
    contrast = float(F(calibration["contrast"]))
    phase = math.atan2(readings[1]/contrast, readings[0]/contrast)
    if not 0 < phase < math.pi:
        raise ArithmeticError("the supplied source clock did not lie on its admitted phase branch")
    return {"X_expectation": readings[0], "Y_expectation": readings[1],
            "contrast": contrast, "decoded_phase": phase,
            "preparation_resources": preparation_cost,
            "purification_reference_columns_retained": initial.shape[1],
            "clock_is_actual_source_execution": True,
            "parameter": "programmed modular parameter; not laboratory seconds"}


def _charged_slots():
    operators = action.identification_operators()
    q, weak = operators["Q"], operators["weak_C2"]
    left = [i for i in range(45) if q[i, i] == -1 and weak[i, i] == s.Rational(3, 4)]
    singlet = [i for i in range(45) if q[i, i] == 1 and weak[i, i] == 0]
    if len(left) != 3 or len(singlet) != 3:
        raise ArithmeticError("the actual charge/Casimir projectors did not isolate three charged leptons")
    return tuple((a, b) for a, b in zip(left, singlet, strict=True))


def _ratio_from_response(car_response, clock_phase):
    cosine = float(complex(car_response).real)
    if not -1 < cosine < 1 or not 0 < clock_phase < math.pi:
        raise ArithmeticError("response/clock lies outside the supplied inverse-phase branch")
    return math.acos(cosine)/clock_phase


def charged_pole_diagnostics(clock=None):
    """CAR residues, clocked frequencies and invariant rest gaps at unequal p."""
    clock = _clock_record() if clock is None else clock
    omega = float(action.source_coefficients()["Omega"])
    kernel, charge = _numeric(action.mass_kernel()), kinetic.nambu_charge()
    rows = []
    for p in MOMENTA:
        symbol = kinetic.bdg_symbol(p)
        for family, (left, conjugate) in enumerate(_charged_slots()):
            indices = [spin*45+left for spin in range(2)] + [90+spin*45+conjugate for spin in range(2)]
            h = symbol[np.ix_(indices, indices)]
            q = charge[np.ix_(indices, indices)]
            if not np.allclose(q, -np.eye(4), atol=0, rtol=0):
                raise ArithmeticError("the purported charged pole is outside the charge-minus-one sector")
            values, vectors = np.linalg.eigh(h)
            plus = vectors[:, values>0]@vectors[:, values>0].conj().T
            minus = vectors[:, values<0]@vectors[:, values<0].conj().T
            # A fixed spin component of the actual charged left field. The
            # full CAR anticommutator is state independent for this fixed law.
            residue_plus, residue_minus = float(plus[0, 0].real), float(minus[0, 0].real)
            response = expm(-1j*float(INTERVAL)*h)[0, 0]
            ratio = _ratio_from_response(response, clock["decoded_phase"])
            positive_frequency = float(values[values>0].mean())
            rest_over_clock = float(abs(kernel[left, conjugate]))/omega
            invariant = ratio**2-.25*math.fsum(x*x for x in p)
            spectral_response = (residue_plus*np.exp(-1j*positive_frequency*float(INTERVAL))
                                 +residue_minus*np.exp(1j*positive_frequency*float(INTERVAL)))
            rows.append({"native_family_axis": family, "momentum": list(p),
                         "electric_charge_of_Nambu_branch": -1,
                         "charged_field_slot": left,
                         "positive_pole_residue": residue_plus,
                         "negative_pole_residue": residue_minus,
                         "CAR_anticommutator_response": _complex(response),
                         "clock_decoded_energy_over_Omega": ratio,
                         "pole_energy_over_Omega_diagnostic": positive_frequency/omega,
                         "invariant_E_over_Omega_squared_minus_p_squared_over4": invariant,
                         "independent_rest_gap_over_Omega": rest_over_clock,
                         "independent_rest_gap_over_Omega_squared": rest_over_clock**2,
                         "invariant_residual": abs(invariant-rest_over_clock**2),
                         "two_signed_pole_response_residual": float(abs(response-spectral_response))})
    return rows


def _annihilators(modes):
    result = []
    for mode in range(modes):
        a = np.zeros((1<<modes, 1<<modes), complex)
        for mask in range(1<<modes):
            if mask & (1<<mode):
                a[mask^(1<<mode), mask] = (-1)**((mask&((1<<mode)-1)).bit_count())
        result.append(a)
    return tuple(result)


def rest_pair_comparison():
    """Literal same-Hamiltonian comparison, including its charged ground state."""
    left, conjugate = _charged_slots()[0]
    coefficient = _numeric(action.mass_kernel())[left, conjugate]
    a, b = _annihilators(2)
    h = coefficient*a.conj().T@b.conj().T+coefficient.conjugate()*b@a
    values, vectors = np.linalg.eigh(h)
    ground, e0 = vectors[:, 0], float(values[0])
    mass = float(abs(coefficient))
    q = -a.conj().T@a+b.conj().T@b
    addition, removal = a.conj().T@ground, a@ground
    number = a.conj().T@a
    centered_number = number-float(np.vdot(ground, number@ground).real)*np.eye(4)
    neutral_vector = centered_number@ground
    euclidean = []
    for tau in (0., .7, 1.3, 2.1):
        contraction = expm(-tau*(h-e0*np.eye(4)))
        add = float(np.vdot(addition, contraction@addition).real)
        remove = float(np.vdot(removal, contraction@removal).real)
        neutral = float(np.vdot(neutral_vector, contraction@neutral_vector).real)
        euclidean.append({"Euclidean_parameter": tau, "addition_correlator": add,
                          "removal_correlator": remove,
                          "same_H_prediction": .5*math.exp(-mass*tau),
                          "neutral_connected_number_correlator": neutral,
                          "neutral_same_H_prediction": .25*math.exp(-2*mass*tau)})
    c1, c2 = euclidean[1], euclidean[-1]
    euclidean_rate = -math.log(c2["addition_correlator"]/c1["addition_correlator"])/(c2["Euclidean_parameter"]-c1["Euclidean_parameter"])
    neutral_rate = -math.log(c2["neutral_connected_number_correlator"]/c1["neutral_connected_number_correlator"])/(c2["Euclidean_parameter"]-c1["Euclidean_parameter"])
    time = float(INTERVAL)
    local = gaussian_source.local_gate("pairing", coefficient, time)
    evolved = local.fock_unitary.conj().T@a@local.fock_unitary
    anticommutator = evolved@a.conj().T+a.conj().T@evolved
    blank = np.array([1., 0., 0., 0.], complex)
    after = local.fock_unitary@blank
    occupied = float(np.vdot(after, a.conj().T@a@after).real)
    def charge(vector):
        return float(np.vdot(vector, q@vector).real/np.vdot(vector, vector).real)
    return {"native_family_axis": 0, "source_mass_parameter": mass,
            "literal_Fock_dimension": 4, "Fock_energies": values.tolist(),
            "ground_electric_charge": float(np.vdot(ground, q@ground).real),
            "addition_electric_charge": charge(addition), "removal_electric_charge": charge(removal),
            "positive_Lorentzian_pole": mass, "negative_Lorentzian_pole": -mass,
            "addition_residue": float(np.vdot(addition, addition).real),
            "removal_residue": float(np.vdot(removal, removal).real),
            "ground_Euclidean_correlations": euclidean,
            "Euclidean_rate_of_same_H_charged_correlator": euclidean_rate,
            "Euclidean_rate_equals_charged_gap_residual": abs(euclidean_rate-mass),
            "same_H_neutral_number_decay_rate": neutral_rate,
            "neutral_decay_equals_twice_charged_gap_residual": abs(neutral_rate-2*mass),
            "neutral_number_commutes_electric_charge_residual": float(np.max(abs(q@number-number@q))),
            "source_Lorentzian_anticommutator_cos_mt_residual": float(np.max(abs(anticommutator-math.cos(mass*time)*np.eye(4)))),
            "blank_quench_number_contrast": 1-2*occupied,
            "blank_quench_prediction_cos_2mt": math.cos(2*mass*time),
            "source_local_block": local.record(),
            "ground_reference_is_empty_Fock_state": False,
            "Euclidean_correlator_requires_ground_reference_and_same_H": True,
            "Euclidean_contraction_executed_as_native_repair": False}


def _nambu_from_majorana(orthogonal):
    modes = len(orthogonal)//2
    transform = np.zeros_like(orthogonal, dtype=complex)
    for j in range(modes):
        transform[2*j, j] = transform[2*j, modes+j] = 1
        transform[2*j+1, j], transform[2*j+1, modes+j] = -1j, 1j
    return transform.conj().T@orthogonal@transform/2


@lru_cache(maxsize=1)
def source_propagating_check():
    """Execute one ±p eight-CAR block through actual routed source words."""
    p = (.13, -.21, .34)
    left, conjugate = _charged_slots()[0]
    coefficient = _numeric(action.mass_kernel())[left, conjugate]
    internal = np.array([[0, coefficient], [coefficient, 0]], complex)
    h, delta = kinetic.momentum_pair_source(p, internal)
    target_generator = np.block([[h, delta], [delta.conj().T, -h.T]])
    target = expm(-1j*float(INTERVAL)*target_generator)
    clock = _clock_record()
    rows = []
    final = None
    for steps in (4, 8):
        source = gaussian_source.evolve(h, delta, float(INTERVAL), steps=steps)
        actual = _nambu_from_majorana(source.orthogonal)
        decoded = _ratio_from_response(actual[0, 0], clock["decoded_phase"])
        invariant = decoded**2-.25*math.fsum(x*x for x in p)
        target_mass_ratio = float(abs(coefficient))/float(action.source_coefficients()["Omega"])
        rows.append({"split_steps": steps,
                     "full_CAR_linear_map_max_residual": float(np.max(abs(actual-target))),
                     "source_CAR_anticommutator": _complex(actual[0, 0]),
                     "clock_decoded_energy_over_Omega": decoded,
                     "invariant_rest_squared": invariant,
                     "rest_squared_residual": abs(invariant-target_mass_ratio**2)})
        final = source
    return {"momentum": list(p), "work_CAR_modes_at_p_and_minus_p": len(h),
            "same_source_clock": clock,
            "native_source_convergence": rows,
            "final_source_Majorana_map": final.orthogonal.tolist(),
            "final_source_execution": final.record(),
            "target_BdG_not_executed_by_source": True,
            "charged_probe_instrument_for_this_new_momentum_menu_executed": False,
            "CAR_response_recovered_from_completed_source_Heisenberg_map": True,
            "interpretation": "additional finite source evolution/calibration diagnostic; no new laboratory or Standard Model identification"}


def report():
    clock = _clock_record()
    rows = charged_pole_diagnostics(clock)
    rest = rest_pair_comparison()
    return {
        "schema": "oph.full_repair_response.mass_observable.v1",
        "observable_definition": "charged CAR response pole with nonzero residue, invariant dispersion and calibrated observer-clock frequency; rest energy is the gap in those units",
        "source_clock": clock,
        "three_axis_charged_response": rows,
        "largest_invariant_residual": max(row["invariant_residual"] for row in rows),
        "smallest_charged_pole_residue": min(min(row["positive_pole_residue"], row["negative_pole_residue"]) for row in rows),
        "source_propagating_response": source_propagating_check(),
        "literal_rest_pair_and_Euclidean_comparison": rest,
        "native_repair_operator_boundary": {
            "native_population_resolvent": "neutral occupation mean response of a completion proposal, with its complement and timeout",
            "charged_response_here": "odd charged CAR-field anticommutator generated by the declared Weyl/Higgs Hamiltonian",
            "same_operator_or_same_history_intertwiner_established": False,
            "native_CP_relaxation_eigenvalue_equals_mass_proved": False,
            "positive_Euclidean_decay_can_equal_mass": "yes, for the same physical Hamiltonian, charged operator with nonzero ground-state spectral weight, rest momentum and matching clock; the literal rest-pair test demonstrates this conditional relation",
            "Laplace_variable_of_native_repair_equals_physical_energy": False},
        "scope": {
            "existing_candidate_mass_readout_method_calibrated": True,
            "new_native_interacting_mass_generation_derived": False,
            "field_and_clock_calibrations_are_declared_candidate_services": True,
            "phase_inverse_requires_admitted_no_alias_branch": True,
            "new_momentum_probe_statistics_measured_by_committed_detectors": False,
            "physical_clock_seconds_or_GeV_selected": False,
            "observed_Standard_Model_mass_spectrum_identified": False,
            "floating_diagnostics_are_outward_physical_error_bounds": False,
            "empirical_mass_or_coupling_inputs": []},
    }

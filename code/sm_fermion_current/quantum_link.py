"""Finite invariant charged-fermion sector of an infinite quantized U1 link.

The electric basis is the integer line, never a cyclic finite register. A
normalized pure fermion/link state is represented by a product of exact
Cayley factors and an injective Gauss-constrained basis embedding.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib

import current as parent

SCHEMA = "oph.sm_fermion_quantum_link.v1"
STEP = F(1, 2)
ACTIVE_EDGE = 0

SCOPE = {
    "finite_CAR_with_quantized_integer_electric_link": True,
    "complete_one_generation_representation_census": True,
    "normalized_pure_entangled_fermion_link_state": True,
    "operator_Gauss_on_every_supported_branch": True,
    "commuting_local_Cayley_factors_exactly_replayed": True,
    "same_Hamiltonian_current_and_electric_transfer": True,
    "exact_continuous_quantum_edge_solution": True,
    "finite_cyclic_electric_cutoff": False,
    "same_preparation_as_uniform_semiclassical_history": False,
    "all_144_links_dynamically_evolved": False,
    "nonabelian_gauge_or_Yukawa_dynamics": False,
    "nonabelian_quantum_Gauss": False,
    "electric_energy_or_plaquette_factors": False,
    "Cayley_to_continuous_error_bound": False,
    "Cayley_word_identified_with_continuous_exponential": False,
    "chiral_continuum_anomaly_measure_constructed": False,
    "native_source_action_state_clock_or_quantization_selected": False,
    "physical_current_or_laboratory_identification": False,
}

PREMISES = [
    "The prepared q5 source addresses, first oriented edge and Pauli frame are supplied.",
    "One registered generation and the finite canonical anticommutation relations are supplied.",
    "Each internal component occupies the same normalized left-end spin orbital; other vertices are empty.",
    "The active electric link is ell2(Z), with integer electric unit, bilateral unitary shift and hbar=1.",
    "Each channel has one fermion and two active spatial-spin modes; the orthogonal spin modes are empty and invariant.",
    "Hopping coefficients are 1 except 2 on the e_c channel; each Cayley factor has step 1/2.",
    "The initial electric state is |0>; all other link momenta remain zero and their hopping terms are omitted.",
    "Exact state amplitudes and the finite factor schedule are mathematical data, not sampled outcomes or source-selected operations.",
]

QUANTUM_LINK = {
    "gauge_constraint_group": "U1_hypercharge",
    "ambient_space": "ell2(Z)",
    "electric_action": "E|ell>=ell|ell>",
    "shift_action": "U|ell>=|ell+1>",
    "commutator": "[E,U]=U",
    "finite_wrap": False,
    "active_edge": ACTIVE_EDGE,
    "other_links": "sharp zero electric values; no other hopping factor acts",
    "Hamiltonian_channel": "kappa_i*(c_L_i^dagger*U^q_i*c_R_i+c_R_i^dagger*U^(-q_i)*c_L_i)",
    "current_channel": "i*kappa_i*q_i*(c_L_i^dagger*U^q_i*c_R_i-c_R_i^dagger*U^(-q_i)*c_L_i)",
    "Heisenberg_electric_rate": "i[H_i,E]=-J_i",
    "factor": "V_i=(I+i*step*H_i/2)^(-1)*(I-i*step*H_i/2)",
    "factor_current": "quadratic expectation in (state_before+state_after)/2; midpoint vector is not normalized",
    "discrete_current_integral": "step*sum_i J_mid_i = -(mean_E_final-mean_E_initial)",
    "instantaneous_current": "expectation of sum_i J_i in the normalized endpoint state",
    "factor_schedule": "channel ids 0 through 14, each once; no continuous-time exponential or physical clock identified",
}


def channels():
    answer = []
    for field in parent.multiplets():
        for component in range(field["multiplicity"]):
            index = len(answer)
            kappa = F(2 if field["name"] == "e_c" else 1)
            r = STEP*kappa/2
            a, beta = (1-r*r)/(1+r*r), 2*r/(1+r*r)
            answer.append({"id": index, "multiplet": field["name"],
                           "internal_component": component,
                           "integer_charge": field["integer_charge"],
                           "mode_indices": [2*index, 2*index+1],
                           "kappa": str(kappa), "step": str(STEP),
                           "left_amplitude": parent.encode(parent.C(a)),
                           "right_amplitude": parent.encode(parent.C(0, -beta))})
    return answer


def amplitudes(row):
    return [parent.C(F(z[0]), F(z[1]))
            for z in (row["left_amplitude"], row["right_amplitude"])]


def convolve(histogram, charge, right_probability):
    result = {}
    for flux, probability in histogram.items():
        for new_flux, weight in ((flux, 1-right_probability),
                                 (flux-charge, right_probability)):
            result[new_flux] = result.get(new_flux, F(0))+probability*weight
    return {flux: probability for flux, probability in result.items() if probability}


def histogram_readout(histogram):
    norm = sum(histogram.values(), F(0))
    if norm != 1:
        raise ArithmeticError("quantum link state normalization")
    mean = sum((flux*p for flux, p in histogram.items()), F(0))
    second = sum((flux*flux*p for flux, p in histogram.items()), F(0))
    purity = sum((p*p for p in histogram.values()), F(0))
    return {"norm_squared": str(norm), "mean_electric": str(mean),
            "electric_second_moment": str(second),
            "electric_variance": str(second-mean*mean),
            "mean_charge_left": str(mean), "mean_charge_right": str(-mean),
            "link_reduced_purity": str(purity),
            "flux_probabilities": [{"flux": flux, "probability": str(p)}
                                   for flux, p in sorted(histogram.items())]}


def branch_census(channel_rows):
    """Retain a digest of all amplitudes and check Gauss on every branch."""
    state = {0: parent.C(1)}
    for row in channel_rows:
        a, b = amplitudes(row)
        bit = 1 << row["id"]
        state = {new_mask: value*factor
                 for mask, value in state.items()
                 for new_mask, factor in ((mask, a), (mask | bit, b))}
    total_charge = sum(row["integer_charge"] for row in channel_rows)
    if total_charge:
        raise ArithmeticError("initial operator Gauss requires neutral generation")
    stream = hashlib.sha256()
    histogram = {}
    for mask, amplitude in sorted(state.items()):
        moved_charge = sum(row["integer_charge"] for row in channel_rows
                           if mask & (1 << row["id"]))
        flux = -moved_charge
        left_charge, right_charge = total_charge-moved_charge, moved_charge
        gauss = [flux-left_charge, -flux-right_charge]
        if gauss != [0, 0]:
            raise ArithmeticError("operator Gauss branch")
        record = {"mask": mask,
                  "occupied_modes": [2*row["id"]+int(bool(mask & (1 << row["id"])))
                                     for row in channel_rows],
                  "electric_flux": flux, "charge_left": left_charge,
                  "charge_right": right_charge,
                  "amplitude": parent.encode(amplitude), "gauss_residual": gauss}
        stream.update(parent.canonical(record))
        histogram[flux] = histogram.get(flux, F(0))+amplitude.norm()
    if set(histogram) != set(range(-18, 19)):
        raise ArithmeticError("complete integer flux support")
    return {"branch_count": len(state), "branch_sha256": stream.hexdigest(),
            "branch_order": "ascending integer mask; bit i selects right mode of channel i",
            "branch_record_fields": ["mask", "occupied_modes", "electric_flux", "charge_left",
                                     "charge_right", "amplitude", "gauss_residual"],
            "norm_squared": str(sum(histogram.values(), F(0))),
            "all_branch_Gauss_residuals": [0, 0],
            "other_vertex_Gauss_residuals": "zero, since their fermion modes and incident electric momenta are empty"}, histogram


def produce():
    geo, rows = parent.geometry(), channels()
    edge = geo["edges"][ACTIVE_EDGE]
    left = [parent.C(F(3, 5)), parent.C(F(4, 5))]
    right = parent.spin_action(edge["axis"], parent.C(1), left, adjoint=True)
    census, full_histogram = branch_census(rows)
    histogram = {0: F(1)}
    history = [{"completed_factors": 0, "factor": None,
                "readout": histogram_readout(histogram)}]
    sum_midpoint, instantaneous = F(0), F(0)
    for row in rows:
        q, kappa = row["integer_charge"], F(row["kappa"])
        a, b = amplitudes(row)
        beta = -b.imag
        midpoint_left, midpoint_right = (1+a.real)/2, beta/2
        current = 2*kappa*q*midpoint_left*midpoint_right
        integral = STEP*current
        old_mean = sum((flux*p for flux, p in histogram.items()), F(0))
        histogram = convolve(histogram, q, b.norm())
        result = histogram_readout(histogram)
        if F(result["mean_electric"])-old_mean != -integral:
            raise ArithmeticError("operator current/electric transfer")
        sum_midpoint += current
        instantaneous += 2*kappa*q*a.real*beta
        history.append({"completed_factors": row["id"]+1,
                        "factor": {"channel": row["id"],
                                   "midpoint_norm_squared": str(midpoint_left**2+midpoint_right**2),
                                   "midpoint_current": str(current),
                                   "discrete_integrated_current": str(integral),
                                   "mean_electric_change": str(-integral)},
                        "readout": result})
    if histogram != full_histogram:
        raise ArithmeticError("factorized versus branch probability census")
    final = histogram_readout(histogram)
    if F(final["link_reduced_purity"]) >= 1:
        raise ArithmeticError("fermion/link entanglement witness")
    if F(final["mean_electric"]) != -F(18144, 7225):
        raise ArithmeticError("nonzero quantum current witness")
    continuous_point_histogram = {0: F(1)}
    for row in rows:
        continuous_point_histogram = convolve(
            continuous_point_histogram, row["integer_charge"],
            F(3, 4) if row["multiplet"] == "e_c" else F(1, 4))
    continuous_point = histogram_readout(continuous_point_histogram)
    if F(continuous_point["mean_electric"]) != -3 or F(continuous_point["electric_variance"]) != F(45, 2):
        raise ArithmeticError("continuous exact quantum witness")
    return {"schema": SCHEMA, "scope": SCOPE, "premises": PREMISES,
            "carrier": geo, "quantum_link": QUANTUM_LINK,
            "channel_basis": {
                "active_edge": ACTIVE_EDGE, "ends": edge["ends"], "axis": edge["axis"],
                "left_spin": [parent.encode(z) for z in left],
                "right_spin": [parent.encode(z) for z in right],
                "spin_relation": "right_spin=(i*sigma_axis)^dagger*left_spin",
                "CAR_mode_order": "L0,R0,L1,R1,...,L14,R14",
                "spectators": "all orthogonal spin modes and all other source sites are empty and invariant under these factors",
                "branch_embedding": "|mask> maps to ordered occupied fermion modes tensor |electric=-sum moved integer charges>",
                "one_particle_modes_on_carrier": 1920, "active_fermion_modes": 30,
                "fermions": 15, "invariant_sector_dimension": 32768},
            "channels": rows,
            "initial_state": {"occupied_modes": list(range(0, 30, 2)), "electric_flux": 0,
                              "total_charge": 0, "norm_squared": "1",
                              "operator_Gauss": "zero at every vertex",
                              "preparation": "localized generation at the left endpoint, distinct from the uniform semiclassical Slater state"},
            "history": history, "final_state": {"branch_census": census, "readout": final},
            "continuous_solution": {
                "Hamiltonian": "sum_i H_i on the same infinite-link invariant sector",
                "commuting_channels": True,
                "channel_left_amplitude": "cos(kappa_i*t)",
                "channel_right_amplitude": "-i*sin(kappa_i*t)",
                "state_embedding": "the same ordered fermion/flux branch embedding; no cyclic cutoff",
                "mean_electric": "6*(sin(t)^2-sin(2*t)^2)",
                "instantaneous_current": "-6*sin(2*t)+12*sin(4*t)",
                "electric_variance": "84*sin(t)^2*cos(t)^2+36*sin(2*t)^2*cos(2*t)^2",
                "electric_derivative": "d(mean_electric)/dt=-instantaneous_current",
                "operator_Gauss": "zero for every real t",
                "exact_point": {"model_time": "pi/6", "readout": continuous_point,
                                "instantaneous_current_Qsqrt3": ["0", "3"]},
                "physical_time_calibrated": False,
                "Cayley_word_is_this_exact_time": False},
            "observables": {"sum_factor_midpoint_currents": str(sum_midpoint),
                            "discrete_integrated_current": str(STEP*sum_midpoint),
                            "final_instantaneous_current": str(instantaneous),
                            "initial_instantaneous_current": "0",
                            "mean_electric_change": final["mean_electric"],
            "operator_Gauss_norm_squared": "0",
                            "entanglement_witness": "pure joint state with mixed electric reduction; exact reduced purity below 1",
                            "finite_support": [-18, 18], "ambient_electric_cutoff": None}}

"""A clocked neutral-pair readout of the full trial's three lepton rates.

At rest the trial has42 disjoint charged pair blocks and6 neutral zero
modes. A fermionic wire ordering puts the three probed blocks first. We
retain the full90-mode Gaussian history and independently execute its
six-mode factor and two quantum clocks on actual source Fock columns.
No physical vacuum preparation, dynamical Higgs selection or mass fit is
inferred from this supplied quadratic background experiment.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from . import action, gaussian_source as g
from ..mass_response import history as instruments, model
from ..mass_response import experiment as old
from ...observer.source_scalar import compiler


INTERVAL = F(1, 4)


def rest_layout():
    delta = np.array(action.pairing_kernel().evalf(), dtype=complex)
    labels = action.mode_labels()
    selected = []
    for label in labels:
        if label["spin"] == 0 and label["internal_type"] == "L" and F(label["electric_Q"]) == -1:
            i = label["index"]
            partners = np.flatnonzero(delta[i] != 0)
            if len(partners) != 1:
                raise ArithmeticError("the selected rest charged mode is not a disjoint pair")
            selected.append((i, int(partners[0])))
    if len(selected) != 3:
        raise ArithmeticError("three native family lepton pairs required")
    pairs, used = list(selected), {i for pair in selected for i in pair}
    for i in range(90):
        if i in used:
            continue
        partners = np.flatnonzero(delta[i] != 0)
        if len(partners) == 1:
            j = int(partners[0])
            if j in used:
                raise ArithmeticError("overlapping quadratic rest blocks")
            pairs.append((i, j)); used.update((i, j))
        elif len(partners) != 0:
            raise ArithmeticError("the source trial is not in its factorized rest domain")
    order = [i for pair in pairs for i in pair]+[i for i in range(90) if i not in used]
    if len(pairs) != 42 or len(set(order)) != 90:
        raise ArithmeticError("missing charged blocks or neutrino modes")
    return delta[np.ix_(order, order)], {
        "initial_fermionic_wire_order": order, "probed_original_mode_pairs": selected,
        "layout_is_initial_encoding_choice": True, "modes_discarded": 0,
        "charged_pair_blocks": 42, "neutral_zero_modes_retained": 6,
        "selected_pair_coefficients": [[float(delta[i,j].real), float(delta[i,j].imag)] for i,j in selected]}


def _clock_preparation():
    calibration = model.ideal_clock_calibration()
    amplitudes = np.sqrt(list(map(lambda x: float(F(x)), calibration["populations"])))
    word = []
    for clock, purifier in ((0, 2), (1, 3)):
        word += compiler.remap(instruments.preparation_word(amplitudes), [clock])
        word += compiler.cnot(clock, purifier)
        word += [compiler.pulse("H", math.pi/2, [clock])]
    initial = np.zeros((16, 1), complex); initial[0,0] = 1
    output, cost = instruments.execute_word(word, 4, initial)
    return output[:,0].reshape(4,4).T, calibration, cost


def _clock_advance(state, duration, *, backend="source", offset=0, qubits=2):
    angle = -model.parameters().clock*float(duration)/2
    if backend == "source":
        return compiler.execute([compiler.pulse("Z", angle, [offset+i]) for i in (0,1)], qubits, state)
    for i in (0,1):
        state = model.apply_subsystem(state, np.diag(np.exp(1j*angle*np.array([1.,-1.]))), [offset+i], qubits)
    return state


def _dense_probe(coefficients, *, backend):
    clock, calibration, clock_cost = _clock_preparation()
    blank = np.zeros(64); blank[0] = 1
    initial = np.stack([np.kron(clock[:,j], blank) for j in range(4)], axis=1)
    source = [g.local_gate("pairing", value, float(INTERVAL)) for value in coefficients]
    c0, c1 = [a[:4,:4] for a in model.annihilators()[:2]]
    targets = [expm(-1j*float(INTERVAL)*(value*c0.conj().T@c1.conj().T+
                                       value.conjugate()*c1@c0)) for value in coefficients]
    def advance(state):
        for j, gate in enumerate(source):
            unitary = gate.fock_unitary if backend == "source" else targets[j]
            state = model.apply_subsystem(state, unitary, [2*j,2*j+1], 8)
        return _clock_advance(state, INTERVAL, backend=backend, offset=6, qubits=8)
    number_cost,clock_read_costs = None,{}
    def number(state, bit):
        nonlocal number_cost
        if backend == "source":
            read = instruments.diagonal_read(state, [1.,0.], targets=[bit])
            number_cost = read.resources
            return read.outputs
        mask = ((np.arange(256) >> bit)&1).astype(bool)
        empty, occupied = state.copy(), state.copy()
        empty[mask], occupied[~mask] = 0, 0
        return empty, occupied
    branches = {(): advance(initial)}
    for bit in (0,2,4):
        branches = {(*record,value): output for record,state in branches.items()
                    for value,output in enumerate(number(state,bit))}
    prefix = {key:old.norm(state) for key,state in branches.items()}
    for record,state in branches.items():
        if record[0]:
            # A declared neutral controller phase; no additional force is
            # inserted into the autonomous interval Hamiltonian.
            if backend == "source":
                state = compiler.execute([compiler.pulse("I", math.pi/4,[2]),
                                          compiler.pulse("Z", -math.pi/4,[2])],8,state)
            else:
                state = model.apply_subsystem(state,np.diag([1.,1j]),[2],8)
        branches[record] = advance(state)
    branches = {(*record,value):output for record,state in branches.items()
                for value,output in enumerate(number(state,2))}
    for bit,axis in ((6,"X"),(7,"Y")):
        outputs = {}
        for record,state in branches.items():
            if backend == "source":
                read = instruments.pauli_read(state,bit,axis=axis)
                states = read.outputs
                clock_read_costs[axis] = read.resources
            else:
                pauli = np.array([[0,1],[1,0]]) if axis == "X" else np.array([[0,-1j],[1j,0]])
                changed = model.apply_subsystem(state,pauli,[bit],8)
                states = ((state+changed)/2,(state-changed)/2)
            for value,output in enumerate(states):
                outputs[(*record,value)] = output
        branches = outputs
    ledger = {}
    if backend == "source":
        base_pulses = (clock_cost["native_pulses"]+2*sum(gate.resources["native_pulses"] for gate in source)
                       +4+4*number_cost["native_pulses"]
                       +sum(cost["native_pulses"] for cost in clock_read_costs.values()))
        ledger = {"base_native_pulses_per_realized_history":base_pulses,
                  "additional_native_pulses_if_first_record1":2,
                  "base_proper_code_service_events":2*base_pulses,
                  "additional_proper_code_service_events_if_first_record1":4,
                  "fresh_pointer_blanks":6,"committed_classical_bits":6,
                  "initial_field_clock_and_purifier_blank_qubits":10,
                  "maximum_including_fresh_pointer_and_untouched_reference":11,
                  "scope":"literal six-mode Fock factor and clocks; not the full90-mode implementation",
                  "floating_and_physical_services_certified":False}
    return branches, prefix, calibration, clock_cost, source, ledger


def infer_ratios(probabilities, calibration, phase_contract=None):
    """Decode three pair frequencies using only records and clock calibration."""
    if phase_contract is None or phase_contract.get("pair_phases_principal_branch_proved_for_declared_trial") is not True:
        raise ValueError("a prior admitted pair-frequency window is required")
    def mean(position):
        return math.fsum((1 if key[position] == 0 else -1)*value for key,value in probabilities.items())
    x,y = mean(4)/calibration["contrast"],mean(5)/calibration["contrast"]
    if abs(math.hypot(x,y)-1) > 1e-7:
        raise ValueError("clock differs from the declared calibrated circle")
    angle = math.atan2(y,x)
    if not 0 < angle < math.pi:
        raise ValueError("ambiguous clock window")
    phases = []
    for position in range(3):
        cosine = mean(position)
        if not -1 <= cosine <= 1:
            raise ValueError("unresolved pair phase")
        phases.append(math.acos(cosine))
    return {"native_family_order_gap_over_clock": [phase/angle for phase in phases],
            "measured_pair_phases": phases, "measured_clock_phase": angle,
            "decoder_uses_Hamiltonian_coefficients_or_simulator_time": False,
            "finite_shot_or_physical_service_errors_certified": False}


def execute():
    delta, layout = rest_layout()
    coefficients = [complex(*value) for value in layout["selected_pair_coefficients"]]
    # Exact elementary inequalities give J<1, Omega<2 and max(F)<1/4.
    # Hence m_l=J*Omega*F<1/2 and 2*m_l*tau<1/4<pi. They concern
    # this declared trial; cosine data alone cannot rule out frequency aliases.
    w0,w1,w2 = model.parameters().weights
    if not (1 < w1*w1/(w0*w2) < 2 and 1 < w2/w0 < F(19,3)
            and F(41,225)+F(7,750)*F(9,4) < F(1,4)
            and all(0 < 2*abs(value)*float(INTERVAL) < .25 for value in coefficients)):
        raise ArithmeticError("the admitted pair-frequency window changed")
    phase_contract = {"pair_phases_principal_branch_proved_for_declared_trial":True,
                      "pair_phase_exact_upper":"1/4",
                      "proof":"J<1, Omega<2, sqrt(5)<9/4, max(F)<1/4; tau=1/4; pi>3",
                      "other_laws_or_unbounded_frequencies_admitted":False}
    print("Executing the full90-mode rest Gaussian law and the literal clocked probe", flush=True)
    first = g.evolve(np.zeros_like(delta),delta,float(INTERVAL),covariance=g.vacuum_covariance(90))
    gaussian = {(): (1.,first.covariance)}
    pointer_cost = None
    for bit in (0,2,4):
        outputs = {}
        for record,(weight,state) in gaussian.items():
            read = g.occupation_read(state,bit,weight=weight)
            pointer_cost = read.resources
            for branch in read.branches:
                outputs[(*record,branch.occupation)] = (branch.weight,branch.covariance)
        gaussian = outputs
    prefix = {key:weight for key,(weight,_) in gaussian.items()}
    outputs = {}
    kick = g.local_gate("onsite",-math.pi/2,1.)
    rotation = first.orthogonal
    for record,(weight,state) in gaussian.items():
        if record[0]:
            transform = np.eye(180); transform[4:6,4:6] = kick.orthogonal[:2,:2]
            state = transform@state@transform.T
        state = rotation@state@rotation.T
        read = g.occupation_read(state,2,weight=weight)
        for branch in read.branches:
            outputs[(*record,branch.occupation)] = (branch.weight,branch.covariance)
    actual,old_prefix,calibration,clock_cost,local,literal_ledger = _dense_probe(coefficients,backend="source")
    target,*_ = _dense_probe(coefficients,backend="continuous")
    probabilities = {key:old.norm(state) for key,state in actual.items()}
    gaussian_error = max(abs(weight-math.fsum(p for key,p in probabilities.items() if key[:4] == record))
                         for record,(weight,_) in outputs.items())
    prefix_error = max(abs(weight-old_prefix[key]) for key,weight in prefix.items())
    old_record_error = max(abs(weight-math.fsum(p for key,p in probabilities.items() if key[:3] == record))
                           for record,weight in old_prefix.items())
    source_cq = old.Result(actual,{}, {},{}, {})
    target_cq = old.Result(target,{}, {},{}, {})
    distance = old.quantum_distance(source_cq,target_cq)
    if max(gaussian_error,prefix_error,old_record_error,distance) > 2e-9:
        raise ArithmeticError("source Gaussian law and literal continuing observer probe disagree")
    decoded = infer_ratios(probabilities,calibration,phase_contract)
    expected = [abs(value)/model.parameters().clock for value in coefficients]
    if max(abs(a-b) for a,b in zip(expected,decoded["native_family_order_gap_over_clock"],strict=True)) > 1e-9:
        raise ArithmeticError("the record decoder does not recover the trial's charged gaps")
    records,field_records = sorted(actual),sorted(outputs)
    arrays = {"records":np.array(["".join(map(str,key)) for key in records]),
              "probe_outputs":np.stack([actual[key] for key in records]),
              "field_records":np.array(["".join(map(str,key)) for key in field_records]),
              "gaussian_weights":np.array([outputs[key][0] for key in field_records]),
              "gaussian_covariances":np.stack([outputs[key][1] for key in field_records])}
    return {
        "schema":"oph.chiral_candidate.clocked_history.v1", "layout":layout,
        "record_order":["native_family0_empty","native_family1_empty","native_family2_empty",
                        "continued_family1_empty","clock_X_plus","clock_Y_plus"],
        "field_record_convention":"0 is empty,1 occupied",
        "clock_calibration":calibration,"clock_final_modular_parameter":str(2*INTERVAL),
        "source_record_probabilities":{"".join(map(str,key)):value for key,value in probabilities.items()},
        "total_record_probability":math.fsum(probabilities.values()),
        "normal_records":len(records),"gaussian_final_field_records":len(field_records),
        "work_preparations":1,"work_resets":0,"all_outcomes_retained":True,
        "old_record_marginal_residual":old_record_error,
        "full90_mode_and_probe_prefix_residual":prefix_error,
        "full90_mode_and_probe_final_residual":gaussian_error,
        "complete_probe_CQ_distance_to_independent_continuous":distance,
        "observer_inverse":decoded,"independent_analytic_rest_gaps_over_clock":expected,
        "phase_window_contract":phase_contract,
        "full_rest_source_evolution":first.record(),
        "clock_preparation_source_resources":clock_cost,
        "full_field_evolution_words_per_nonabort_history":2,
        "full_field_source_pulses_for_both_intervals":2*first.resources["native_pulses"],
        "conditional_phase_source_pulses_if_first_record1":kick.resources["native_pulses"],
        "gaussian_field_pointer_source_pulses_for_four_reads":4*pointer_cost["native_pulses"],
        "literal_probe_source_resources":literal_ledger,
        "equivalent_source_realizations":"full90 Gaussian uses two half pair blocks per interval; literal factor uses one full block; equal in exact arithmetic, floating residuals and costs are separate",
        "initial_field_blank_modes":90,"quantum_clock_qubits":2,"untouched_clock_reference_qubits":2,
        "source_factorization_domain":"quadratic fixed-neutral-Higgs rest law in its declared initial fermion wire order",
        "no_2_to_90_Fock_vector_allocated":True,"floating_errors_are_diagnostics":True,
        "physical_units_or_observed_mass_identification":False,
    },arrays

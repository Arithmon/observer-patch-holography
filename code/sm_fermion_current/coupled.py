"""Exact rational coupled CAR/rotor transport on the prepared golden carrier.

The convex electric kinetic law is declared and nonquadratic. Its Cayley
phase drift is exact over the declared duration, while the fermionic factor
uses implicit midpoint. This is a supplied semiclassical action, not a
quantized gauge field or an identification with physical Maxwell dynamics.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json

import current as base

ROOT = base.ROOT
SCHEMA = "oph.sm_fermion_coupled.v1"
PARENT = "code/sm_fermion_current/current_receipt.json"
KINETIC_LAMBDA = F(1, 2)
ELECTRIC_MASS = F(1)
ELECTRIC_DURATION = F(1, 2)
MAGNETIC_DURATION = F(1, 2)
MAGNETIC_WEIGHT = F(1)
RETURN_EDGE = 14
MAX_RUN_BYTES = 4_500_000
MAX_VALUE_BITS = 12_000

SCOPE = {
    "same_finite_hybrid_action_for_CAR_electric_and_Wilson_factors": True,
    "all_144_electric_links_and_108_plaquettes_executed": True,
    "actual_electric_to_link_to_matter_read_after_write": True,
    "exact_discrete_continuity_and_expectation_Gauss": True,
    "full_one_generation_internal_multiplicities_retained": True,
    "number_conserving_Slater_and_classical_field_interpretation": True,
    "nonquadratic_electric_kinetic_law_explicitly_declared": True,
    "quadratic_Maxwell_electric_kinetic_law": False,
    "exact_full_Hamiltonian_exponential": False,
    "full_word_energy_conservation_claimed": False,
    "quantum_operator_Gauss_or_quantized_gauge_history": False,
    "three_generation_execution": False,
    "global_Z6_quotient_selected": False,
    "nonzero_Higgs_or_nonabelian_field_execution": False,
    "continuous_trajectory_error_bound": False,
    "native_action_state_clock_or_physical_identification": False,
}

LAW = {
    "gauge_group": "SU3xSU2xU1_direct_product_cover",
    "Hamiltonian": "sum_e T_lambda(E_e)+sum_f beta*(1-Re(W_f))+sum_e <H_CAR,e>",
    "symplectic_form": "sum_e dE_e wedge dtheta_e + i sum_s multiplicity_s dconj(psi_s) wedge dpsi_s",
    "kinetic_lambda": str(KINETIC_LAMBDA), "electric_mass": str(ELECTRIC_MASS),
    "electric_duration": str(ELECTRIC_DURATION),
    "magnetic_duration": str(MAGNETIC_DURATION), "magnetic_weight": str(MAGNETIC_WEIGHT),
    "matter_duration": str(base.STEP), "matter_hopping": str(base.HOPPING),
    "electric_kinetic": "T_lambda(E)=2*E/lambda*atan(lambda*E/(2*M))-2*M/lambda^2*log(1+(lambda*E/(2*M))^2)",
    "kinetic_derivative": "2/lambda*atan(lambda*E/(2*M))",
    "kinetic_second_derivative": "1/(M*(1+(lambda*E/(2*M))^2))",
    "kinetic_parameter_boundary": "lambda is a fixed action parameter, not an adaptive step; this fixture separately chooses electric duration=lambda",
    "electric_drift": "U_new=(1+i*lambda*E/(2*M))/(1-i*lambda*E/(2*M))*U_old; E unchanged",
    "Wilson_loop": "W_f=product_boundary U_e^orientation",
    "magnetic_kick": "E_e_new=E_e_old-magnetic_duration*beta*orientation*Im(W_f); links unchanged",
    "matter_factor": "the same CAR midpoint factor and variational expectation current as the pinned prefix",
    "Gauss": "outward_divergence(E)-sum_s multiplicity_s*q_s*orbital_site_norm_squared",
    "schedule": "pinned full matter sweep; all electric drifts; all Wilson kicks; one matter return at edge14; public readout",
    "return_edge": RETURN_EDGE,
    "control": "drift_disabled retains each link unchanged with a link-only read; it is a counterfactual schedule, not the full declared action word",
    "energy_boundary": "symplectic splitting and exact Gauss do not assert conservation of the summed Hamiltonian over the finite word",
}

REDUCTION = {
    "families": 3,
    "occupied_family": 0,
    "vacuum_families": [1, 2],
    "fermion_normal_covariance": "P_s=I_internal_multiplicity tensor |psi_s><psi_s|",
    "anomalous_covariance": "<c_i c_j>=0 in the fixed-number Slater state",
    "nonabelian_links": "identity",
    "nonabelian_electric_momenta": "zero",
    "Higgs_and_conjugate_momentum": "zero",
    "Yukawa_coefficients": "all 27 arbitrary complex coefficients of three declared families",
    "vacuum_family_invariance": "Higgs zero removes Yukawa pair terms; family-diagonal number-conserving hopping preserves the other two Fock vacua",
    "omitted_mean_force_boundary": "traceless internal generators have zero current expectation; Higgs Yukawa force uses the zero anomalous covariance; this is not zero operator fluctuation",
    "physical_boundary": "supplied finite mean-field restriction, not a chiral regulator or a gauge-constrained quantum state",
}


def carrier():
    pilot = base.load_parent("code/sm_abelian_reduction/pilot.py", "coupled_golden_geometry")
    geometry = pilot.geometry()
    return {"base": base.geometry(),
            "plaquettes": [{"boundary": row["boundary"]} for row in geometry["plaquettes"]]}


def decode(value):
    return base.C(F(value[0]), F(value[1])) if isinstance(value, list) else F(value)


def loop_value(values, boundary):
    result = base.C(1)
    for edge, orientation in boundary:
        phase = values[f"link/{edge}"]
        result = result*(phase if orientation == 1 else phase.conj())
    return result


def readout(values, geometry, fields):
    result = base.observables(values, geometry["base"], fields)
    loops = [loop_value(values, row["boundary"]) for row in geometry["plaquettes"]]
    matter_energy_terms = []
    for edge, row in enumerate(geometry["base"]["edges"]):
        u, v = row["ends"]
        edge_energy = F(0)
        for field in fields:
            name, q = field["name"], field["integer_charge"]
            product = base.inner(base.spinor(values, name, u), base.spin_action(
                row["axis"], values[f"link/{edge}"]**q, base.spinor(values, name, v)))
            edge_energy += 2*base.HOPPING*field["multiplicity"]*product.real
        matter_energy_terms.append(edge_energy)
    magnetic_energy_terms = [MAGNETIC_WEIGHT*(1-z.real) for z in loops]
    def enclose(terms):
        denominator = 2**40
        lower_numerator = sum((term.numerator*denominator)//term.denominator for term in terms)
        return [str(F(lower_numerator, denominator)), str(F(lower_numerator+len(terms), denominator))]
    result.update({"Wilson_loops": [base.encode(z) for z in loops],
                   "magnetic_energy_interval": enclose(magnetic_energy_terms),
                   "magnetic_energy_terms": [str(value) for value in magnetic_energy_terms],
                   "matter_energy_interval": enclose(matter_energy_terms),
                   "matter_energy_terms": [str(value) for value in matter_energy_terms],
                   "energy_interval_grid_denominator": str(2**40),
                   "electric_kinetic_arguments": [str(values[f"electric/{edge}"])
                                                   for edge in range(144)],
                   "total_energy_expression": "sum_e T_lambda(electric_kinetic_arguments[e])+sum magnetic_energy_terms+sum matter_energy_terms"})
    if any(value != "0" for value in result["gauss_residual"]):
        raise ArithmeticError("coupled expectation Gauss")
    if any(value != "1" for value in result["orbital_norms"].values()):
        raise ArithmeticError("coupled occupied orbital norm")
    if any(z.norm() != 1 for z in loops):
        raise ArithmeticError("Wilson loop unit phase")
    return result


def require_small(value):
    numbers = (value.real, value.imag) if isinstance(value, base.C) else (value,)
    if any(max(x.numerator.bit_length(), x.denominator.bit_length()) > MAX_VALUE_BITS for x in numbers):
        raise ValueError("declared exact-arithmetic bit budget exceeded")


def execute(prefix, geometry, fields, cohort):
    values, versions, writers = {}, {}, {}
    for event in prefix["events"]:
        for row in event["writes"]:
            name = row["port"]
            values[name], versions[name], writers[name] = decode(row["value"]), row["version"], event["id"]
    offset, previous = len(prefix["events"]), prefix["final_event_hash"]
    events, checkpoints = [], []
    byte_count = len(base.canonical(prefix))

    def emit(operation, read_names, writes, metadata):
        nonlocal previous, byte_count
        names = sorted(set(read_names))
        for value in writes.values():
            require_small(value)
        event = {"id": offset+len(events), "operation": operation, "metadata": metadata,
                 "reads": [{"port": name, "version": versions[name], "writer": writers[name],
                            "value": base.encode(values[name])} for name in names],
                 "parents": sorted({writers[name] for name in names}),
                 "writes": [], "previous_hash": previous}
        for name, value in sorted(writes.items()):
            values[name], versions[name], writers[name] = value, versions[name]+1, event["id"]
            event["writes"].append({"port": name, "version": versions[name], "value": base.encode(value)})
        event["event_hash"] = base.digest(event)
        byte_count += len(base.canonical(event))
        if byte_count > MAX_RUN_BYTES:
            raise ValueError("declared bounded run byte budget exceeded")
        previous = event["event_hash"]
        events.append(event)

    def checkpoint(label):
        result = readout(values, geometry, fields)
        reads = [name for name in values if not name.startswith("factor_current/")]
        emit("coupled_readout", reads, {}, {"checkpoint": label})
        checkpoints.append({"event": events[-1]["id"], "readout": result})

    initial_electric = values[f"electric/{RETURN_EDGE}"]
    for edge in range(144):
        link, electric = f"link/{edge}", f"electric/{edge}"
        if cohort == "drift_disabled":
            emit("electric_drift_disabled", [link], {link: values[link]}, {"edge": edge})
        else:
            x = KINETIC_LAMBDA*values[electric]/(2*ELECTRIC_MASS)
            rotation = (1+base.I*x)/(1-base.I*x)
            if rotation.norm() != 1:
                raise ArithmeticError("electric drift unit phase")
            emit("electric_drift", [electric, link], {link: rotation*values[link]}, {"edge": edge})
    after_drift_link = values[f"link/{RETURN_EDGE}"]
    for face, row in enumerate(geometry["plaquettes"]):
        boundary = row["boundary"]
        phase = loop_value(values, boundary)
        names = [f"{kind}/{edge}" for edge, _ in boundary for kind in ("link", "electric")]
        writes = {f"electric/{edge}": values[f"electric/{edge}"]-
                  MAGNETIC_DURATION*MAGNETIC_WEIGHT*orientation*phase.imag
                  for edge, orientation in boundary}
        emit("Wilson_kick", names, writes, {"plaquette": face})
    checkpoint("after_gauge_factors")

    edge, row = RETURN_EDGE, geometry["base"]["edges"][RETURN_EDGE]
    u, v = row["ends"]
    names, updates, current = [f"link/{edge}", f"electric/{edge}"], {}, F(0)
    for field in fields:
        name, q = field["name"], field["integer_charge"]
        left, right = base.spinor(values, name, u), base.spinor(values, name, v)
        phase = values[f"link/{edge}"]**q
        new_left, new_right = base.cayley(left, right, row["axis"], phase)
        middle_left = [(x+y)/2 for x, y in zip(left, new_left, strict=True)]
        middle_right = [(x+y)/2 for x, y in zip(right, new_right, strict=True)]
        j = -2*base.HOPPING*field["multiplicity"]*q*base.inner(
            middle_left, base.spin_action(row["axis"], phase, middle_right)).imag
        old_p = sum((z.norm() for z in left), F(0))
        new_p = sum((z.norm() for z in new_left), F(0))
        if field["multiplicity"]*q*(new_p-old_p) != -base.STEP*j:
            raise ArithmeticError("coupled local current/charge transfer")
        current += j
        for site, output in ((u, new_left), (v, new_right)):
            for spin, value in enumerate(output):
                port = base.orbital_port(name, site, spin)
                names.append(port)
                updates[port] = value
    updates[f"electric/{edge}"] = values[f"electric/{edge}"]-base.STEP*current
    updates[f"factor_current/{edge}"] = current
    emit("coupled_matter_return", names, updates, {"edge": edge})
    return_event = events[-1]["id"]
    checkpoint("final")
    return {"cohort": cohort, "prefix": prefix, "events": events, "checkpoints": checkpoints,
            "feedback_witness": {"edge": edge, "electric_before_drift": str(initial_electric),
                                 "link_after_drift": base.encode(after_drift_link),
                                 "electric_drift_event": offset+edge,
                                 "matter_return_event": return_event,
                                 "return_midpoint_current": str(current)},
            "final_event_hash": previous}


def produce():
    raw = (ROOT/PARENT).read_bytes()
    parent_packet = json.loads(raw)
    prefixes = {run["cohort"]: run for run in parent_packet["runs"]}
    geometry, fields = carrier(), base.multiplets()
    runs = [execute(prefixes["phase_intervention" if cohort == "drift_disabled" else cohort],
                    geometry, fields, cohort)
            for cohort in ("baseline", "phase_intervention", "gauge_copy", "drift_disabled")]
    final = [run["checkpoints"][-1]["readout"] for run in runs]
    if final[1] != final[2]:
        raise ArithmeticError("coupled gauge-copy public mismatch")
    return_difference = F(runs[1]["feedback_witness"]["return_midpoint_current"])-F(
        runs[3]["feedback_witness"]["return_midpoint_current"])
    if not return_difference:
        raise ArithmeticError("vacuous electric-to-matter feedback")
    differences = {key: [str(F(b)-F(a)) for a, b in zip(final[0][key], final[1][key], strict=True)]
                   for key in ("charge", "electric", "instantaneous_current")}
    return {"schema": SCHEMA, "scope": SCOPE, "law": LAW, "mean_field_restriction": REDUCTION,
            "parent": {"path": PARENT, "sha256": hashlib.sha256(raw).hexdigest()},
            "carrier": geometry, "multiplets": fields, "runs": runs,
            "resource_contract": {"maximum_run_bytes": MAX_RUN_BYTES, "maximum_register_integer_bits": MAX_VALUE_BITS,
                                  "maximum_receipt_bytes": 20_000_000},
            "comparisons": {"gauge_invariant_public_records_equal": True,
                            "intervention_minus_baseline_final": differences,
                            "return_current_minus_disabled_drift": str(return_difference),
                            "mean_field_only": True, "operator_Gauss_claimed": False}}

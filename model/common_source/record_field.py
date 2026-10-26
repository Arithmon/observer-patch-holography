"""A specified CAR walk joined to the record-corrected native code.

Six copies of record_code's six-dimensional code each factor as a qubit and
a three-dimensional spectator: logical_label=2*spectator+occupation.
Native repair and recovery act identically on all six directions. The
source-compiled fermionic swaps are an additional,
fixed operational choice; the native classical transcript cannot select them.

The native tensor instrument is evaluated by its exact all-column local
certificate, not by materializing 512**6 work amplitudes. Its complete record
tree is retained in product-factor form, including raw zero outcomes. The
64-dimensional field word and reference columns are actually executed.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from ..observer.source_scalar import compiler
from ..observer.source_discrimination import circuit
from . import record_code


CELLS = 3
MODES = 2*CELLS
DIMENSION = 1 << MODES
CAPS = ((0, 1, 5, 7), (0, 5, 7, 10))


def mode(cell, component):
    if type(cell) is not int or type(component) is not int or component not in (0, 1):
        raise ValueError("integer cell and component zero or one required")
    return 2*(cell % CELLS)+component


def annihilator(port):
    if type(port) is not int or not 0 <= port < MODES:
        raise ValueError("one of the six logical modes required")
    result = np.zeros((DIMENSION, DIMENSION), dtype=np.int64)
    for mask in range(DIMENSION):
        if mask >> port & 1:
            result[mask ^ (1 << port), mask] = (-1)**((mask & ((1 << port)-1)).bit_count())
    return result


def second_quantization(permutation):
    """Independent signed exterior lift, including every occupied sector."""
    if (type(permutation) not in (list, tuple) or len(permutation) != MODES
            or any(type(j) is not int for j in permutation)
            or sorted(permutation) != list(range(MODES))):
        raise ValueError("a permutation of all six modes required")
    result = np.zeros((DIMENSION, DIMENSION), dtype=np.int64)
    for mask in range(DIMENSION):
        occupied = [permutation[j] for j in range(MODES) if mask >> j & 1]
        inversions = sum(a > b for i, a in enumerate(occupied) for b in occupied[i+1:])
        result[sum(1 << j for j in occupied), mask] = (-1)**inversions
    return result


def swap_word(a, b):
    """Gamma(transposition), with the spectator parity even at the ring edge.

    exp(-i*pi*(c_a^dag c_b+c_b^dag c_a)/2) exp(i*pi*(n_a+n_b)/2).
    The Jordan--Wigner strings are an explicit source addressing cost.
    """
    if type(a) is not int or type(b) is not int or not 0 <= a < MODES or not 0 <= b < MODES or a == b:
        raise ValueError("two different field modes required")
    a, b = sorted((a, b))
    word = []
    for endpoint in ("X", "Y"):
        label = ["I"]*MODES
        label[a] = label[b] = endpoint
        for j in range(a+1, b):
            label[j] = "Z"
        word += compiler.pauli_rotation(label, -math.pi/4)
    word += [compiler.pulse("I", math.pi/2, [a]),
             compiler.pulse("Z", -math.pi/4, [a]),
             compiler.pulse("Z", -math.pi/4, [b])]
    return word


def walk_word():
    first = [(mode(x, 0), mode(x, 1)) for x in range(CELLS)]
    second = [(mode(x, 1), mode(x+1, 0)) for x in range(CELLS)]
    return [event for a, b in first+second for event in swap_word(a, b)]


def decoded_register_walk_word():
    """Literal wiring after the six actual nine-bit native code decoders.

    Bit9*j is the field occupation; the other bits include both spectator
    bits and the retained decoder syndrome. No syndrome is selected/reset.
    Full54-bit tensor execution is replaced by the exact repeated-block
    identity and the source-executed64-column common field block.
    """
    return compiler.remap(walk_word(), tuple(9*j for j in range(MODES)))


def permutation():
    return tuple(mode(x+(1 if component == 0 else -1), component)
                 for x in range(CELLS) for component in (0, 1))


@lru_cache(maxsize=1)
def source_walk():
    actual = compiler.execute(walk_word(), MODES)
    target = second_quantization(permutation())
    error = float(np.max(abs(actual-target)))
    if not np.isfinite(actual).all() or not math.isfinite(error) or error > 2e-11:
        raise ArithmeticError("source walk differs from signed exterior transport")
    return np.frombuffer(actual.tobytes(), dtype=actual.dtype).reshape(actual.shape)


def algebra_and_continuity():
    local_field = np.zeros((6, 6), dtype=np.int64)
    for spectator in range(3):
        local_field[2*spectator, 2*spectator+1] = 1
    if not np.array_equal(local_field@local_field.T+local_field.T@local_field, np.eye(6)):
        raise ArithmeticError("the full six-code spectator factor lost its CAR identity")
    fields = [annihilator(j) for j in range(MODES)]
    identity = np.eye(DIMENSION, dtype=np.int64)
    for i in range(MODES):
        for j in range(MODES):
            if (np.any(fields[i]@fields[j]+fields[j]@fields[i])
                    or not np.array_equal(fields[i]@fields[j].T+fields[j].T@fields[i], identity*(i == j))):
                raise ArithmeticError("the complete six-mode CAR algebra failed")
    number = [a.T@a for a in fields]
    u = second_quantization(permutation())
    first_current = [number[mode(x, 0)]-number[mode(x, 1)] for x in range(CELLS)]
    second_current = [number[mode(x, 0)]-number[mode(x+1, 1)] for x in range(CELLS)]
    for x in range(CELLS):
        da = u.T@number[mode(x, 0)]@u-number[mode(x, 0)]
        db = u.T@number[mode(x, 1)]@u-number[mode(x, 1)]
        if (not np.array_equal(da, -first_current[x]+second_current[(x-1) % CELLS])
                or not np.array_equal(db, first_current[x]-second_current[x])):
            raise ArithmeticError("the two local swap currents fail continuity")
    if not np.array_equal(u.T@sum(number)@u, sum(number)):
        raise ArithmeticError("the walk failed logical charge conservation")
    return {"full_CAR_dimension": DIMENSION, "CAR_modes": MODES,
            "all_CAR_and_continuity_identities_exact": True,
            "charge": "Q=sum_x(n_Ax+n_Bx), a declared logical charge",
            "first_layer_current_A_to_B": "n_Ax-n_Bx",
            "second_layer_current_Bx_to_Axplus1": "n_Ax-n_B(x+1)",
            "cell_continuity": "Delta Q_x=J2_(x-1)-J2_x",
            "field_carrier": "each code6=C^3_spectator tensor C^2_occupation, label=2*s+n",
            "total_corrected_code_dimension": 6**MODES,
            "CAR_representation_multiplicity": 3**MODES,
            "maximum_exact_CAR_modes_on_entire_code_by_dimension_divisibility": MODES,
            "dimension_obstruction": "a representation of CAR_m=M_(2**m) with unit I_(6**6) requires 2**m to divide 6**6; a larger selected subcode is a different encoding",
            "code_directions_discarded_per_carrier": 0,
            "spectator_sectors_identical_by_exact_tensor_factorization": True,
            "native_total_occupation_changes_under_logical_creation": False,
            "native_electric_charge_identification": False,
            "spatial_locality": "two nearest-cell fermionic swap layers on a three-cell ring",
            "source_addressing": "Jordan-Wigner parity strings are compiled, including the ring wrap; physical local hardware routing is not derived"}


def dispersion():
    """All lattice momenta and both orientations, in the same cycle clock."""
    u = source_walk()
    one = u[np.ix_([1 << j for j in range(MODES)], [1 << j for j in range(MODES)])]
    rows = []
    for integer in (-1, 0, 1):
        k = 2*math.pi*integer/CELLS
        for component in (0, 1):
            velocity = 1 if component == 0 else -1
            vector = np.zeros(MODES, complex)
            for x in range(CELLS):
                vector[mode(x, component)] = np.exp(1j*k*x)/math.sqrt(CELLS)
            multiplier = np.exp(-1j*velocity*k)
            residual = float(np.max(abs(one@vector-multiplier*vector)))
            if not math.isfinite(residual) or residual > 2e-11:
                raise ArithmeticError("source dispersion differs from the fixed translation")
            rows.append({"momentum_integer": integer, "k": k, "orientation": velocity,
                         "quasienergy_per_cycle": velocity*k,
                         "multiplier": [float(multiplier.real), float(multiplier.imag)],
                         "source_residual": residual})
    # The anticommutator response is state-independent for this number-preserving
    # quadratic unitary. Compare actual all-sector conjugation, not just one-body data.
    response = np.empty((MODES, MODES), complex)
    fields = [annihilator(j) for j in range(MODES)]
    for i in range(MODES):
        evolved = u.conj().T@fields[i]@u
        for j in range(MODES):
            response[i, j] = np.trace(evolved@fields[j].T+fields[j].T@evolved)/DIMENSION
    if not np.isfinite(response).all() or np.max(abs(response-one)) > 3e-11:
        raise ArithmeticError("charged response and one-particle walk disagree")
    return {"modes": rows, "clock": "one complete corrected-repair plus two-layer walk cycle",
            "dispersion": "epsilon_A(k)=k, epsilon_B(k)=-k modulo 2*pi",
            "charged_retarded_anticommutator_response_residual": float(np.max(abs(response-one))),
            "zero_rest_quasienergy_exact": True, "nontrivial_spatial_dispersion": True,
            "physical_chirality_or_3plus1_spacetime_claim": False,
            "mass_parameters_fitted": False}


def record_factors(law):
    """Complete raw local record alphabet; each probability includes its coin.

    A common five-bit selector is uniform on 0..31. Accepted addresses select
    all thirty source seams; addresses30,31 halt without another work cycle.
    Conditional on a seam, the six local record factors multiply. This is an
    exact compressed representation of every full native tensor branch.
    """
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("C or E required")
    rows = []
    for address, seam in enumerate(record_code.SEAMS):
        branches = []
        for coin in (0, 1):
            for raw in range(1 << circuit.RECORD_BITS[law]):
                weight = record_code.branch_weight(law, seam, coin, raw)
                images = record_code.branch_images(law, seam, coin, raw)
                if len(images) != 6 or any(len(row) != len(images[0]) for row in images):
                    raise ArithmeticError("native code branch ranks depend on logical input")
                expected = F(len(images[0]), 12)
                if weight != expected:
                    raise ArithmeticError("native code factor and exact branch images disagree")
                labels = circuit.record_values(law, coin)
                branches.append({"coin": coin, "raw": raw,
                                 "value": labels.get(raw),
                                 "abort": raw not in labels,
                                 "weight_exact": str(weight), "zero_outcome_retained": weight == 0})
        if sum((F(b["weight_exact"]) for b in branches), F()) != 1:
            raise ArithmeticError("a local repair factor lost an outcome")
        rows.append({"selector_address": address, "seam": list(seam),
                     "selector_weight_exact": "1/32", "local_record_factor": branches})
    if len(rows) != 30:
        raise ArithmeticError("the declared finite selector requires all thirty native seams")
    return {"law": law, "accepted_seams": rows,
            "terminal_selector_addresses": [{"address": 30, "weight_exact": "1/32"},
                                             {"address": 31, "weight_exact": "1/32"}],
            "native_carrier_factors_per_accepted_seam": MODES,
            "joint_branch_weight": "(1/32)*product_site(code_branch_weight)",
            "complete_native_record_addresses_per_accepted_cycle": 30*(2 << circuit.RECORD_BITS[law])**MODES,
            "raw_abort_policy": "complete the six-site native record batch, retain all work and records, then halt before field decode/walk if any raw record is an abort; these branches have exact zero weight on the ideal code",
            "valid_zero_weight_record_policy": "retain its work output with identity recovery; no conditioning or deletion",
            "corrected_branch_action": "sqrt(joint_branch_weight)*I on all six logical factors, including arbitrary references",
            "representation": "exact factorized complete instrument; not materialized 512**6 native tensor states",
            "source_code_recovery_certificate": "model.common_source.record_code; local all-six-column identities compose tensorially"}


def record_summary(law):
    """Deduplicate a common factor only after checking all source seams."""
    full = record_factors(law)
    common = full["accepted_seams"][0]["local_record_factor"]
    if any(row["local_record_factor"] != common for row in full["accepted_seams"]):
        raise ArithmeticError("native record distributions depend on the seam")
    return {**{key: value for key, value in full.items() if key != "accepted_seams"},
            "accepted_selector_seam_menu": [row["seam"] for row in full["accepted_seams"]],
            "each_accepted_selector_address_weight_exact": "1/32",
            "common_local_record_factor": common,
            "factor_equality_checked_on_all_thirty_source_seams": True}


def clock_histories(cap, deadline=2):
    """Every finite selector-address path through first incident or timeout."""
    if type(deadline) is not int or not 1 <= deadline <= 3:
        raise ValueError("deadline from one through three required")
    if (type(cap) not in (tuple, list) or not cap
            or any(type(x) is not int or not 0 <= x < 12 for x in cap)
            or len(set(cap)) != len(cap)):
        raise ValueError("a nonempty native port access region required")
    result = []
    def visit(prefix, elapsed, probability):
        for address in range(32):
            path, weight = prefix+(address,), probability/32
            if address >= 30:
                flag, cycles = "selector_timeout", elapsed
            else:
                cycles = elapsed+1
                incident = bool(set(cap).intersection(record_code.SEAMS[address]))
                flag = "success" if incident else "deadline_timeout" if cycles == deadline else None
            if flag is None:
                visit(path, cycles, weight)
            else:
                result.append({"selector_addresses": list(path), "elapsed_cycles": cycles,
                               "incident_ticks": int(flag == "success"), "flag": flag,
                               "weight_exact": str(weight)})
    visit((), 0, F(1))
    if sum((F(row["weight_exact"]) for row in result), F()) != 1:
        raise ArithmeticError("finite clock paths were omitted")
    return result


def clock_report():
    rows = []
    k = 2*math.pi/CELLS
    for cap in CAPS:
        histories = clock_histories(cap)
        weights = {flag: sum((F(r["weight_exact"]) for r in histories if r["flag"] == flag), F())
                   for flag in ("success", "selector_timeout", "deadline_timeout")}
        response = sum(float(F(r["weight_exact"]))*np.exp(-1j*k*r["elapsed_cycles"])
                       for r in histories if r["flag"] == "success")
        success_cycles = sorted({r["elapsed_cycles"] for r in histories if r["flag"] == "success"})
        groups = {}
        for history in histories:
            key = history["flag"], history["elapsed_cycles"], history["incident_ticks"]
            count, weight = groups.get(key, (0, F()))
            groups[key] = count+1, weight+F(history["weight_exact"])
        aggregate = [{"flag": flag, "elapsed_cycles": cycles, "incident_ticks": ticks,
                      "path_count": count, "weight_exact": str(weight)}
                     for (flag, cycles, ticks), (count, weight) in sorted(groups.items())]
        rows.append({"cap": list(cap), "complete_selector_history_count": len(histories),
                     "terminal_history_groups": aggregate,
                     "complete_local_record_factors_attached_to_each_accepted_address": True,
                     "flag_probabilities_exact": {key: str(value) for key, value in weights.items()},
                     "success_elapsed_cycle_values": success_cycles,
                     "unconditional_success_mode_response": [float(response.real), float(response.imag)],
                     "same_one_incident_tick_has_different_elapsed_phases": len(success_cycles) > 1})
    joint = joint_clock_histories()
    groups = {}
    for history in joint["histories"]:
        key = history["global_terminal_flag"], history["elapsed_cycles"], tuple(history["first_incident_cycle_by_observer"])
        count, weight = groups.get(key, (0, F()))
        groups[key] = count+1, weight+F(history["weight_exact"])
    aggregate = [{"global_terminal_flag": flag, "elapsed_cycles": cycles,
                  "first_incident_cycle_by_observer": list(hits), "path_count": count,
                  "weight_exact": str(weight)}
                 for (flag, cycles, hits), (count, weight) in sorted(groups.items(), key=lambda item: repr(item[0]))]
    joint_summary = {**{key: value for key, value in joint.items() if key != "histories"},
                     "complete_selector_history_count": len(joint["histories"]),
                     "terminal_history_groups": aggregate,
                     "source_executed_stopped_channel": stopped_clock_execution(aggregate)}
    return {"deadline_cycles": 2, "observers": rows,
            "joint_clock_histories": joint_summary,
            "enumeration_contract": "clock_histories and joint_clock_histories enumerate every raw five-bit selector path; report groups identical terminal clock labels only after exact summation, with path counts retained",
            "per_history_response": "G_h(i,j)=p_h*(u**elapsed_cycles)_ij, including terminal flags",
            "clock_compatibility": "both observers use the recorded elapsed cycle count; their incident counters alone do not determine phase",
            "timeouts_restarted_or_success_renormalized": False,
            "physical_duration_of_a_cycle_derived": False}


def stopped_clock_execution(groups):
    """Source-execute every clock-conditioned work block, including failures.

    Grouping is classical: sqrt(group_probability)*U**T is a compressed
    representation of identical conditional work maps. Raw selector and native
    outcomes remain in their separately enumerated factors.
    """
    unitary = source_walk()
    target = second_quantization(permutation())
    actual_powers = [np.eye(DIMENSION), unitary, unitary@unitary]
    target_powers = [np.eye(DIMENSION, dtype=np.int64), target, target@target]
    vector = np.zeros(DIMENSION, complex)
    for x in range(CELLS):
        vector[1 << mode(x, 0)] = np.exp(2j*math.pi*x/CELLS)/math.sqrt(CELLS)
    differences, gram, rows = [], np.zeros((DIMENSION, DIMENSION), complex), []
    total = F()
    for group in groups:
        weight, elapsed = F(group["weight_exact"]), group["elapsed_cycles"]
        if weight < 0 or type(elapsed) is not int or elapsed not in (0, 1, 2):
            raise ValueError("positive clock probabilities and at most two elapsed cycles required")
        total += weight
        root = math.sqrt(float(weight))
        actual = root*actual_powers[elapsed]
        expected = root*target_powers[elapsed]
        differences.append(actual-expected)
        gram += actual.conj().T@actual
        phase = np.vdot(vector, actual_powers[elapsed]@vector)
        first_phases = []
        for first in group["first_incident_cycle_by_observer"]:
            if first is None:
                first_phases.append(None)
            else:
                value = np.vdot(vector, actual_powers[first]@vector)
                first_phases.append([float(value.real), float(value.imag)])
        rows.append({"elapsed_cycles": elapsed, "weight_exact": str(weight),
                     "global_terminal_flag": group["global_terminal_flag"],
                     "actual_mode_characteristic": [float(phase.real), float(phase.imag)],
                     "actual_first_incident_characteristic_by_observer": first_phases})
    if total != 1:
        raise ArithmeticError("stopped source channel does not include every clock outcome")
    error = float(np.linalg.norm(np.concatenate(differences), 2))
    gram_error = float(np.linalg.norm(gram-np.eye(DIMENSION), 2))
    if any(not math.isfinite(value) for value in (error, gram_error)) or error > 5e-11 or gram_error > 1e-10:
        raise ArithmeticError("stopped source clock channel differs from the independent charged walk")
    return {"all_work_input_columns_per_clock_block": DIMENSION,
            "complete_clock_block_count": len(groups),
            "complete_retained_stopped_operator_error": error,
            "complete_stopped_gram_error": gram_error,
            "terminal_source_characteristics": rows,
            "zero_cycle_selector_timeouts_and_two_cycle_deadlines_retained": True,
            "per_block_action": "sqrt(group_probability)*source_walk**elapsed_cycles",
            "classical_grouping_coherently_recombines_records": False}


def joint_clock_histories():
    """Two passive first-hit counters on one actual finite selector tape.

    Both clocks observe the same source advance. A first hit stores its cycle
    timestamp but does not stop the other observer or disturb the field. The
    common run ends at a selector timeout or the two-cycle deadline.
    """
    result = []
    def visit(prefix, elapsed, first_hits, probability):
        for address in range(32):
            path, weight = prefix+(address,), probability/32
            hits = list(first_hits)
            if address >= 30:
                cycles, flag = elapsed, "selector_timeout"
            else:
                cycles = elapsed+1
                seam = record_code.SEAMS[address]
                for observer, cap in enumerate(CAPS):
                    if hits[observer] is None and set(cap).intersection(seam):
                        hits[observer] = cycles
                flag = "deadline" if cycles == 2 else None
            if flag is None:
                visit(path, cycles, tuple(hits), weight)
            else:
                result.append({"selector_addresses": list(path), "elapsed_cycles": cycles,
                               "first_incident_cycle_by_observer": hits,
                               "incident_tick_by_observer": [int(t is not None) for t in hits],
                               "global_terminal_flag": flag, "weight_exact": str(weight)})
    visit((), 0, (None, None), F(1))
    if sum((F(row["weight_exact"]) for row in result), F()) != 1:
        raise ArithmeticError("the shared two-observer tape lost a history")
    return {"histories": result, "one_source_tape_for_both_observers": True,
            "actual_source_control": "accepted selector address executes one corrected repair and one source walk; addresses30/31 halt before either operation",
            "clock_read_service": "first incident seam stores the current completed-cycle counter; clocks do not measure the field",
            "conditional_clock_response": "u**first_incident_cycle at that observer's tagged event",
            "clock_error": "zero in the declared finite classical register semantics; no physical duration calibration"}


def detector_experiment():
    """Two local number reads, complete update, feedback and reference memory.

    Initial (|A0>R0+|B0>R1)/sqrt2. After one walk read n_A1; on
    outcome1 swap A1 and B1. After another corrected repair and walk read n_B0.
    Native record factors multiply every work branch; they are never summed
    coherently. Selector timeouts at either cycle retain their current state.
    """
    u = source_walk()
    initial = np.zeros((DIMENSION, 2), complex)
    initial[1 << mode(0, 0), 0] = initial[1 << mode(0, 1), 1] = 1/math.sqrt(2)
    nfirst = annihilator(mode(1, 0)).T@annihilator(mode(1, 0))
    nsecond = annihilator(mode(0, 1)).T@annihilator(mode(0, 1))
    identity = np.eye(DIMENSION)
    first_effects, second_effects = (identity-nfirst, nfirst), (identity-nsecond, nsecond)
    feedback = compiler.execute(swap_word(mode(1, 0), mode(1, 1)), MODES)
    accepted = F(15, 16)
    branches, arrays = [], {"initial": initial}
    for first in (0, 1):
        after_read = first_effects[first]@u@initial
        continued = feedback@after_read if first else after_read
        arrays[f"after_first_{first}"] = continued
        for raw in (30, 31):
            arrays[f"selector_timeout_second_{first}_{raw}"] = continued*math.sqrt(float(accepted/F(32)))
        for second in (0, 1):
            state = second_effects[second]@u@continued
            arrays[f"terminal_{first}_{second}"] = state
            arrays[f"weighted_terminal_{first}_{second}"] = state*float(accepted)
            probability = float(np.vdot(state, state).real)
            exact_conditional = F(1, 2) if first == second else F()
            if not np.isfinite(state).all() or not math.isfinite(probability) or abs(probability-float(exact_conditional)) > 3e-11:
                raise ArithmeticError("source read-feedback-read joint law changed")
            branches.append({"first_read": first, "second_read": second,
                             "unconditional_probability_exact": str(accepted**2*exact_conditional),
                             "source_probability_before_selector_weights": probability,
                             "zero_probability_output_retained": exact_conditional == 0})
    timeout_first = F(1, 16)
    for raw in (30, 31):
        arrays[f"selector_timeout_first_{raw}"] = initial/math.sqrt(32)
    timeout_second_each = accepted*F(1, 2)*F(1, 16)
    if timeout_first+2*timeout_second_each+sum((F(row["unconditional_probability_exact"]) for row in branches), F()) != 1:
        raise ArithmeticError("the two-read experiment lost a terminal outcome")
    no_feedback = []
    identity_control = []
    for first in (0, 1):
        for second in (0, 1):
            state = second_effects[second]@u@first_effects[first]@u@initial
            no_feedback.append(float(np.vdot(state, state).real))
            unchanged = first_effects[first]@initial
            if first:
                unchanged = feedback@unchanged
            unchanged = second_effects[second]@unchanged
            identity_control.append(float(np.vdot(unchanged, unchanged).real))
    weighted_arrays = [array for name, array in arrays.items()
                       if name.startswith("selector_timeout_") or name.startswith("weighted_terminal_")]
    total_source_probability = sum(float(np.vdot(array, array).real) for array in weighted_arrays)
    if not math.isfinite(total_source_probability) or abs(total_source_probability-1) > 4e-11:
        raise ArithmeticError("source-weighted terminal arrays omit an outcome")
    # Duplicate reads on the same event are compatible because the effects are
    # orthogonal projectors. This is checked as an operator identity, including
    # both outcomes; it does not assert repeatability after further transport.
    for effects in (first_effects, second_effects):
        for i in (0, 1):
            for j in (0, 1):
                if not np.array_equal(effects[i]@effects[j], effects[i] if i == j else np.zeros_like(identity)):
                    raise ArithmeticError("local observer readout effects disagree")
    return {"initial_state": "(|A0>R0+|B0>R1)/sqrt2", "reference_columns_retained": 2,
            "first_local_read": "n_A1", "record_feedback": "if first_read=1 apply FSWAP(A1,B1)",
            "second_local_read": "n_B0 after another corrected native repair and walk",
            "complete_detector_outputs": branches,
            "first_selector_timeout_probability_exact": str(timeout_first),
            "second_selector_timeout_probability_each_first_read_exact": str(timeout_second_each),
            "both_selector_timeout_raw_addresses_retained_in_record_factors": True,
            "spectator_code_directions_source_padding_and_export_failure_outputs_retained": True,
            "ideal_source_padding_and_export_failure_probability_exact": "0",
            "no_feedback_detector_probabilities_before_selector_weights": no_feedback,
            "identity_walk_control_detector_probabilities_before_selector_weights": identity_control,
            "same_site_nondemolition_duplicate_read_verified_exact": True,
            "weighted_terminal_source_probability_sum": total_source_probability,
            "terminal_array_scope": "both raw timeout addresses and all read outcomes include their complete coarse probability weights; earlier native and selector records remain separate product factors, never coherently recombined",
            "conditioning": "each unnormalized branch is retained; a nonzero observed branch uses its own normalized post-read state",
            "detector_source_scope": "logical number projector and committed record are declared read services; walk and feedback words are source executed",
            "native_tensor_execution_scope": "exact local record_code branch identity composed on six carriers; native 512**6 state is not materialized"}, arrays


def declared_services(law, pulse_resources):
    """Count supplied service invocations, without inventing their precision.

    Each epsilon denotes a bound on that service's half diamond error. The
    conditional sum is valid if such physical bounds are supplied. Storage
    epsilon is for the complete experiment window, whose duration is not
    calibrated here. Native record commits include the corresponding write.
    """
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("C or E required")
    raw_bits = circuit.RECORD_BITS[law]
    archive = 10+12*(1+raw_bits)+2
    coefficients = {
        "logical_reference_input_preparation": 1,
        "six_qubit_code_blank_register_preparation": 6,
        "fair_classical_selector_bit": 10,
        "native_coin_blank_register": 12,
        "native_raw_record_blank_register": 12,
        "native_coin_or_raw_record_commit": 24,
        "selector_five_bit_commit": 2,
        "local_field_number_read_and_commit": 2,
        "source_advance_or_halt_control": 2,
        "native_seam_routing_control": 12,
        "native_record_recovery_control": 12,
        "first_field_read_feedback_control": 1,
        "clock_cap_incidence_test": 4,
        "clock_counter_or_first_hit_write": 4,
        "classical_storage_bit_for_whole_window": archive+6,
        "quantum_work_and_reference_storage_for_whole_window": 1,
        "source_analogue_pulse": pulse_resources["native_pulses"],
        "proper_code_coherent_load": pulse_resources["proper_code_load_events"],
        "proper_code_complete_export": pulse_resources["proper_code_export_events"],
    }
    return {"maximum_selector_proposals": 2, "selector_bits_per_proposal": 5,
            "native_attempts_and_coin_record_pairs": 12,
            "native_raw_record_bits_per_attempt": raw_bits,
            "initial_blank_code_registers": 6, "blank_qubits_per_code_register": 6,
            "field_detector_invocations": 2,
            "maximum_archive_bits": archive,
            "derived_clock_storage_bits": 6,
            "maximum_archive_and_clock_bits": archive+6,
            "native_work_qubits": 54, "untouched_reference_dimension_in_experiment": 2,
            "physical_half_diamond_error_sum_coefficients": coefficients,
            "conditional_physical_budget_formula": "sum_service invocation_count[service]*epsilon[service]; every epsilon is an independently supplied whole-service half-diamond bound",
            "physical_epsilon_values_available": False,
            "calibrated_physical_duration_energy_or_total_error_available": False,
            "cost_scope": "source pulses and proper-code transfers are counted; preparation, recording, control and memory counts are supplied-service obligations, without calibrated time or energy costs",
            "zero_probability_raw_and_service_failure_outcomes_discarded": False}


def source_composition_report(native_budget=None):
    """All-branch numerical budget for two complete corrected field cycles.

    Local all30 schedule bounds are composed before taking a reference or
    instrument marginal. This is a floating calculation, not outward-rounded
    analysis or a physical apparatus error estimate. Projective read/blank/
    selector/archive services use their explicitly declared ideal contracts.
    """
    native_budget = record_code.source_composition_report() if native_budget is None else native_budget
    encoder_log, decoder_log = (native_budget[name] for name in ("encoder_log_error", "decoder_log_error"))
    walk_error = float(np.linalg.norm(source_walk()-second_quantization(permutation()), 2))
    fb_permutation = list(range(MODES))
    fb_permutation[2], fb_permutation[3] = fb_permutation[3], fb_permutation[2]
    feedback_word = swap_word(2, 3)
    feedback_error = float(np.linalg.norm(compiler.execute(feedback_word, MODES)
                                         -second_quantization(fb_permutation), 2))
    numbers = [encoder_log, decoder_log, walk_error, feedback_error]
    numbers += [native_budget["laws"][law]["maximum_corrected_attempt_log_error"] for law in ("C", "E")]
    if any(not math.isfinite(value) or value < 0 for value in numbers):
        raise ArithmeticError("nonfinite or negative source composition input")
    encode_resources, decode_resources = record_code.encoding_resources(), record_code.decoding_resources()
    walk_resources, feedback_resources = compiler.resources(walk_word()), compiler.resources(feedback_word)
    keys = ("native_pulses", "single_code_pulses", "pair_code_pulses", "identity_phase_pulses",
            "proper_code_load_events", "proper_code_export_events", "proper_code_service_events")
    rows = {}
    for law in ("C", "E"):
        local = native_budget["laws"][law]
        log_error = (6*encoder_log+12*(local["maximum_corrected_attempt_log_error"]+decoder_log+encoder_log)
                     +2*math.log1p(walk_error)+math.log1p(feedback_error))
        resources = {key: 6*encode_resources[key]
                     +12*(local["maximum_corrected_attempt_resources"][key]+decode_resources[key]+encode_resources[key])
                     +2*walk_resources[key]+feedback_resources[key] for key in keys}
        epsilon = math.expm1(log_error)
        rows[law] = {"maximum_complete_two_cycle_operator_error_diagnostic": epsilon,
                     "maximum_complete_cq_diamond_error_diagnostic": epsilon*(2+epsilon),
                     "maximum_complete_cq_half_trace_error_diagnostic": epsilon*(1+epsilon/2),
                     "maximum_two_cycle_resources": resources,
                     "declared_service_inventory": declared_services(law, resources),
                     "all_native_seams_raw_records_and_timeout_prefixes_covered": True,
                     "reference_extension_changes_bound": False}
    return {"laws": rows, "walk_operator_error": walk_error, "feedback_operator_error": feedback_error,
            "source_composition": "six initial encoders; each of two cycles has six corrected native attempts, six actual decoders, one walk, optional first-read feedback and six encoders",
            "literal_decoded_field_register_targets": [9*j for j in range(MODES)],
            "total_native_work_qubits": 54,
            "actual_field_common_block_input_columns": 64,
            "spectator_and_syndrome_directions_are_retained": True,
            "field_read_service": "committed local low-bit number read; all outcomes retained",
            "selector_service": "five declared fair classical bits per proposal; codes30/31 halt before source advance",
            "failure_control": "raw native abort halts after its complete record batch; a proper-code export failure retains its failure flag and halts; ideal valid-code probability is zero, while numerical or physical error budgets are applied to the complete flagged channel",
            "native_conditional_source_words_selected_by_archived_records": True,
            "complete_instrument_error_argument": "coin preparation and controlled native dilation precede an exact orthogonal record split; record-controlled recoveries form one direct-sum operator whose error is the maximum block error, including identity on zero/abort blocks; measurement and classical copying contract trace norm",
            "channel_error_conversion": "||A-B||<=epsilon and ||B||=1 imply ||A(.)A*-B(.)B*||_diamond<=epsilon*(2+epsilon), preserved by final classical dephasing",
            "numerical_bound_is_outward_rounded_or_physical": False,
            "logical_input_preparation_selector_archive_and_detector_physical_costs_included": False}


def report(*, native_source_budget=None):
    detector, _ = detector_experiment()
    target = second_quantization(permutation())
    return {"schema": "oph.common_source.record_field.v1",
            "algebra": algebra_and_continuity(), "dispersion": dispersion(),
            "source": {"all_field_input_columns_executed": DIMENSION,
                       "source_walk_max_entry_residual": float(np.max(abs(source_walk()-target))),
                       "walk_resources": compiler.resources(walk_word()),
                       "feedback_resources": compiler.resources(swap_word(mode(1, 0), mode(1, 1))),
                       "numerical_residuals_are_physical_error_certificates": False},
            "complete_source_composition": source_composition_report(native_source_budget),
            "native_records": {law: record_summary(law) for law in ("C", "E")},
            "clocks": clock_report(), "two_read_experiment": detector,
            "nonselection": {"identity_and_walk_have_identical_complete_native_record_laws": True,
                             "reason": "every corrected native record is a scalar identity on the code; intervening logical unitaries and normalized interventions cannot change its distribution",
                             "charged_field_interventions_distinguish_them": True,
                             "executed_detector_countercontrol": detector["identity_walk_control_detector_probabilities_before_selector_weights"],
                             "walk_selected_by_native_current_records": False},
            "scope": {"finite_encoded_CAR_transport_with_complete_records": True,
                      "logical_charge_is_native_occupation_charge": False,
                      "physical_spacetime_clock_or_Standard_Model_identification": False,
                      "empirical_mass_inputs": []}}

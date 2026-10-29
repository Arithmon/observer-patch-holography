"""Move the sharp-q=1 code to its recorded core, keeping every stop.

The 60 live core labels index supplied active encoders, not coherent spatial
position states. Positive endpoint biases change flag probabilities but not
the scalar flagged-code identities. All probabilities below are exact.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math

import numpy as np
import sympy as sp

from ..common_source import native_routing as routing, record_code
from ..exchange_sector import sector
from ..local_repair_family import local_bias as local, mixing
from ..observer import record_recovery as records
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from ..particle_screening import screening
from ..screened_recovery import protocol as fixed

CORES = tuple(combinations(range(12), 2))
LIVE = tuple(c for c in CORES if screening.sharp_code(c) is not None)
STOP = tuple(c for c in CORES if c not in LIVE)
SEAMS = tuple(records.seams())
SIZE = 512
INITIAL = (0, 1)


def parameters():
    """The fixed port decoration used by the record-derived clock branch."""
    return replace(mixing.canonical_parameters(), c_weight=F(1, 2))


def _core(center):
    center = sector.core(center)
    if center not in LIVE:
        raise ValueError("a nonopposite core admitting an exact q=1 qubit required")
    return center


def supports(center):
    return screening.sharp_code(_core(center))


def _label(center, law, seam, coin, raw):
    center = _core(center)
    law, seam, coin, raw = fixed._label(law, seam, coin, raw)
    return center, law, seam, coin, raw


def branch(center, law, seam, coin, raw):
    """Exact branch data, including valid zero slots and raw abort slots."""
    return _branch(*_label(center, law, seam, coin, raw))


@lru_cache(None)
def _branch(center, law, seam, coin, raw):
    columns = supports(center)
    label = circuit.record_values(law, coin).get(raw)
    images = []
    moved = set()
    for column in columns:
        row = []
        for pair in column:
            x = tuple(int(i in pair) for i in range(12))
            y, record = records.recorded_step(law, x, seam, coin)
            if label == record[-1]:
                row.append(circuit.encode(y))
                moved.add(sector.moved_core(center, law, record))
        images.append(tuple(row))
    if len(images[0]) != len(images[1]) or len(set(sum(images, ()))) != sum(map(len, images)):
        raise ArithmeticError("a native flag lost its scalar injective code restriction")
    if len(moved) > 1 or bool(moved) != bool(images[0]):
        raise ArithmeticError("the complete flag does not determine the transported core")
    target = next(iter(moved), None)
    return tuple(images), target, label is None


def branch_probability(center, law, seam, coin, raw, value=None):
    p = parameters() if value is None else value
    mixing._parameters(p)
    images, _, _ = branch(center, law, seam, coin, raw)
    return local.coin_probability(p, seam, coin)*F(len(images[0]), len(supports(center)[0]))


def encoder_permutation(center):
    return _encoder_permutation(_core(center))


@lru_cache(None)
def _encoder_permutation(center):
    return routing.complete_permutation(tuple(
        (logical+8*k, circuit.encode(tuple(int(i in pair) for i in range(12))))
        for logical, row in enumerate(supports(center)) for k, pair in enumerate(row)))


def branch_permutation(center, law, seam, coin, raw):
    images, _, _ = branch(center, law, seam, coin, raw)
    return routing.complete_permutation(tuple(
        (logical+8*k, address) for logical, row in enumerate(images)
        for k, address in enumerate(row)))


def ideal_code(center):
    result = np.zeros((SIZE, 2), complex)
    for logical, row in enumerate(supports(center)):
        for pair in row:
            result[circuit.encode(tuple(int(i in pair) for i in range(12))), logical] = 1/math.sqrt(len(row))
    return result


def _encode_register(center, value):
    mixed = compiler.execute(fixed._mixing_word(len(supports(center)[0])), 9, value)
    return record_code._apply_permutation(encoder_permutation(center), mixed)


def encode(center, logical):
    center = _core(center)
    value = np.asarray(logical)
    if (value.ndim not in (1, 2) or value.shape[0] != 2
            or (value.ndim == 2 and value.shape[1] == 0)
            or value.dtype.kind not in "iufc" or not np.isfinite(value).all()):
        raise ValueError("finite logical/reference amplitudes required")
    padded = np.zeros((SIZE,)+value.shape[1:], complex)
    padded[:2] = value
    return _encode_register(center, padded)


def recover(center, law, seam, coin, raw, state):
    """Full-register emitted source unitary, or identity on stops/zero slots.

    On a forbidden transported core the native branch state is kept and the
    controller halts. There is no postselection, re-anchoring or silent retry.
    """
    value = routing._state(state)
    images, target, _ = branch(center, law, seam, coin, raw)
    if not images[0] or target not in LIVE:
        return value
    value = record_code._apply_permutation(branch_permutation(center, law, seam, coin, raw), value, True)
    value = compiler.execute(fixed._mixing_word(len(images[0]), True), 9, value)
    return _encode_register(target, value)


def recovery_resources(center, law, seam, coin, raw):
    images, target, _ = branch(center, law, seam, coin, raw)
    zero = dict.fromkeys(routing.COUNT_KEYS, 0)
    if not images[0] or target not in LIVE:
        return zero
    return fixed._counts(
        record_code._permutation_resources(branch_permutation(center, law, seam, coin, raw)),
        fixed._mixing_resources(len(images[0]), True),
        fixed._mixing_resources(len(supports(target)[0])),
        record_code._permutation_resources(encoder_permutation(target)))


def source_branches(center, law, seam, state, value=None):
    """Replay native source dilations and all raw recovery/stop outputs.

    The source backend compiles its biased placement coin; the recovery
    controller never reads those weights. This is the source law's declared
    apparatus input, not information supplied to the record-clock estimator.
    """
    center, law, seam, _, _ = _label(center, law, seam, 0, 0)
    p = parameters() if value is None else value
    mixing._parameters(p)
    value = routing._state(state)
    route = record_code._route_permutation(seam)
    routed = record_code._apply_permutation(route, value)
    columns = routed.reshape(SIZE, -1)
    count = columns.shape[1]
    local_columns = columns.reshape(64, 8, count).transpose(1, 0, 2).reshape(8, 64*count)
    slots = 1 << circuit.RECORD_BITS[law]
    blank = np.zeros((8*slots, 64*count), complex)
    blank[:8] = local_columns
    angle = -math.asin(math.sqrt(float(local.coin_probability(p, seam, 1))))
    coin_amplitudes = compiler.execute([compiler.pulse("Y", angle, [0])], 1, [1., 0.])
    result = []
    for coin in (0, 1):
        executed = (recorded._repair(law, coin)@blank).reshape(slots, 8, 64, count)
        for raw in range(slots):
            amplitude = executed[raw].transpose(1, 0, 2).reshape(SIZE, count)*coin_amplitudes[coin]
            amplitude = record_code._apply_permutation(route, amplitude, True)
            if value.ndim == 1:
                amplitude = amplitude[:, 0]
            images, target, abort = branch(center, law, seam, coin, raw)
            status = "raw_abort" if abort else "zero" if not images[0] else "opposite_core_stop" if target in STOP else "continue"
            result.append({"coin": coin, "raw": raw, "core_before": center,
                           "core_after": target, "status": status,
                           "probability": branch_probability(center, law, seam, coin, raw, p),
                           "amplitude": recover(center, law, seam, coin, raw, amplitude)})
    return result


def kernel(value=None, law="C"):
    """60 live rows x66 target cores. Column order is CORES; stops stay in."""
    p = parameters() if value is None else value
    mixing._parameters(p)
    law = records.candidate(law)
    rows = []
    total = sum(p.seam_weights)
    for center in LIVE:
        row = Counter({center: 1-p.repair_probability})
        for seam, seam_weight in zip(SEAMS, p.seam_weights, strict=True):
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    probability = branch_probability(center, law, seam, coin, raw, p)
                    if probability:
                        _, target, _ = branch(center, law, seam, coin, raw)
                        row[target] += p.repair_probability*F(seam_weight, total)*probability
        if sum(row.values()) != 1:
            raise ArithmeticError("the complete moving-code row loses a stop")
        rows.append(tuple(row[c] for c in CORES))
    return tuple(rows)


def stopping_certificate(value=None):
    p = parameters() if value is None else value
    full = kernel(p)
    if full != kernel(p, "E"):
        raise ArithmeticError("C/E moving-core marginals differ")
    q = sp.Matrix([[sp.Rational(row[CORES.index(c)]) for c in LIVE] for row in full])
    hazard = tuple(sum(row[CORES.index(c)] for c in STOP) for row in full)
    two_step_stop = tuple(hazard[i]+sum(F(q[i,j])*hazard[j] for j in range(60)) for i in range(60))
    epsilon = min(two_step_stop)
    if epsilon <= 0:
        raise ArithmeticError("the two-step absorption witness is not positive")
    # Positive two-step hazard makes every putative closed live class leak.
    # Thus (I-Q)^-1 exists, and its exact solution is the finite expected exit.
    t = (sp.eye(60)-q).inv()*sp.ones(60, 1)
    if (sp.eye(60)-q)*t != sp.ones(60, 1) or any(v < 1 for v in t):
        raise ArithmeticError("expected stop count failed its original equations")
    start = LIVE.index(INITIAL)
    survival = []
    distribution = [F(i == start) for i in range(60)]
    wanted = {0, 1, 2, 10, 30, 100}
    for count in range(101):
        if count in wanted:
            survival.append({"attempts": count, "survival_probability": str(sum(distribution))})
        if count < 100:
            distribution = [sum(distribution[i]*F(q[i,j]) for i in range(60)) for j in range(60)]
    return {
        "parameters": p.description(), "live_cores": [list(c) for c in LIVE],
        "stopping_cores": [list(c) for c in STOP],
        "complete_kernel_target_order": [list(c) for c in CORES],
        "complete_kernel": [[str(v) for v in row] for row in full],
        "law_C_E_core_kernel_equal": True,
        "one_attempt_stop_probabilities": [str(v) for v in hazard],
        "minimum_two_attempt_stop_probability": str(epsilon),
        "survival_bound_any_initial_core": "Pr(T>2n) <= (1-epsilon)^n",
        "closed_nonstopping_classes": 0,
        "eventual_stop_probability": "1",
        "expected_stopping_attempts": [str(F(v)) for v in t],
        "initial_core": list(INITIAL), "initial_expected_stopping_attempts": str(F(t[start])),
        "initial_movement_probability": str(1-full[start][CORES.index(INITIAL)]),
        "initial_survival": survival,
        "stopping_parameter": "native repair attempt count, including declared holds",
        "calibrated_physical_lifetime": False,
    }


def exact_service_certificate():
    census = Counter()
    transitions = set()
    examples = {}
    for center in LIVE:
        cols = supports(center)
        for law in ("C", "E"):
            for seam in SEAMS:
                total = F()
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        images, target, abort = branch(center, law, seam, coin, raw)
                        census["all_raw_slots"] += 1
                        probability = branch_probability(center, law, seam, coin, raw)
                        total += probability
                        if not images[0]:
                            census["raw_abort_slots" if abort else "zero_valid_slots"] += 1
                            continue
                        status = "stop" if target in STOP else "continue"
                        census[f"{law}_{status}_positive_flags"] += 1
                        if status == "continue":
                            # Exact flagged Knill-Laflamme and target map:
                            # decode two orthogonal normalized uniform columns,
                            # then encode V_target; no logical normalization.
                            if any(sector.charge(tuple(int(i in pair) for i in range(12)), target) != 1
                                   for row in supports(target) for pair in row):
                                raise ArithmeticError("transported encoder leaves sharp q=1")
                            transitions.add((center, target))
                        else:
                            examples.setdefault(law, {"core": list(center), "seam": list(seam),
                                "coin": coin, "raw": raw, "target": list(target),
                                "branch_probability_given_seam": str(probability),
                                "normalized_native_image_columns": [list(row) for row in images]})
                if total != 1:
                    raise ArithmeticError("scalar code probabilities do not sum to one")
    if len(LIVE) != 60 or len(STOP) != 6 or census["all_raw_slots"] != 43200:
        raise ArithmeticError("the complete moving code domain changed")
    return {"counts": dict(sorted(census.items())),
            "directed_successful_core_transitions_including_holds": len(transitions),
            "branch_identity": "R_(c,h) K_h V_c = sqrt(p_(c,h)) V_c' on every continuing flag",
            "failure_identity": "opposite-core flags retain K_h V_c and halt; no renormalization",
            "arbitrary_logical_reference_preserved_on_each_known_flag": True,
            "probabilities_independent_of_unknown_logical_state": True,
            "identity_holds_and_history_controlled_C_E_seams_supported": True,
            "code_recovery_does_not_read_port_weights": True,
            "successful_core_label_reset": False,
            "coherent_superposition_of_core_labels_constructed": False,
            "opposite_stop_examples": examples}


def exterior_certificate():
    """Actual biased reference, not the old transported sector weights.

    The sharp-q=1 exterior number is one; measuring it proves the exact
    universal lower bound 1-pi(q=1), attained by the conditioned reference.
    The code's worst one-port distance is a separate, weaker diagnostic.
    """
    p = parameters()
    pi = local.reference(p)
    states = local.basis()
    rows = []
    for center in LIVE:
        outside = tuple(i for i in range(12) if i not in center)
        mass = sum(w for x,w in zip(states, pi, strict=True) if sector.charge(x, center) == 1)
        marginal = tuple(sum(w*x[i] for x,w in zip(states, pi, strict=True)) for i in range(12))
        occupation = tuple(F(sum(i in pair for pair in supports(center)[0]), len(supports(center)[0])) for i in range(12))
        if any(occupation[i] != F(sum(i in pair for pair in supports(center)[1]), len(supports(center)[1])) for i in range(12)):
            raise ArithmeticError("port occupations reveal logical input")
        rows.append({"core": list(center), "reference_q1_mass": str(mass),
                     "sharp_q1_full_exterior_distance_minimum": str(1-mass),
                     "code_worst_single_exterior_port_distance": str(max(abs(occupation[i]-marginal[i]) for i in outside)),
                     "safe_exterior_ports": list(screening.safe_exterior(center))})
    return {"reference": "fixed PR1151 product-pair pi(theta), theta=(25,25,8,...,8)",
            "port_weights_are_not_transported_with_core": True,
            "universal_bound_domain": "every quantum state sharp q=1 at the specified classical core",
            "bound_witness": "full exterior number measurement; conditioned pi attains bound",
            "rows": rows, "exterior_reference_equality_obtained": False}


def clock_compatibility():
    """Code invariance under the *same fixed* modular diagonal group."""
    theta = parameters().port_weights
    rows = []
    for center in LIVE:
        distributions = [Counter(theta[a]*theta[b] for a,b in row) for row in supports(center)]
        # Disjoint supports mean code invariance iff each logical column has
        # a single pair weight. Equal distributions make the compressed clock
        # scalar; they do not imply invariance when multiple phases occur.
        invariant = all(len(d) == 1 for d in distributions)
        fraction = F(distributions[0].get(200, 0), len(supports(center)[0]))
        if distributions[0] != distributions[1] or any(set(d)-{64,200} for d in distributions):
            raise ArithmeticError("the two-weight scalar clock compression changed")
        rows.append({"core": list(center), "pair_weight_counts": [{str(k):v for k,v in sorted(d.items())} for d in distributions],
                     "modular_code_leakage_sine_squared_coefficient": str(4*fraction*(1-fraction)),
                     "code_invariant_under_all_modular_parameters": invariant,
                     "compressed_modular_flow_scalar": distributions[0] == distributions[1]})
    return {"rows": rows,
            "invariant_core_count": sum(r["code_invariant_under_all_modular_parameters"] for r in rows),
            "time_parameter": "dimensionless reference modular parameter",
            "code_leakage_at_parameter_s": "4r(1-r) sin^2(s log(25/8)/2); r=count(pair weight200)/column size",
            "all_moving_codes_share_a_preserved_modular_subspace": all(r["code_invariant_under_all_modular_parameters"] for r in rows)}


def compiled_layout_certificate():
    """Check actual emitted Gray schedules, including all512 padding rows."""
    checked = set()
    inverse_checked = set()

    def verify(permutation, inverse=False):
        cache = inverse_checked if inverse else checked
        if permutation in cache:
            return
        actual = record_code._integer_schedule_action(routing.permutation_edges(permutation, inverse=inverse))
        expected = tuple(permutation.index(i) for i in range(SIZE)) if inverse else permutation
        if actual != expected:
            raise ArithmeticError("an emitted full-register Gray schedule differs from its map")
        cache.add(permutation)

    for center in LIVE:
        verify(encoder_permutation(center))
        verify(encoder_permutation(center), True)
        for law in ("C", "E"):
            for seam in SEAMS:
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        images, target, _ = branch(center, law, seam, coin, raw)
                        if not images[0] or target in STOP:
                            continue
                        verify(branch_permutation(center, law, seam, coin, raw), True)
    return {"encoders": 60, "all_source_addresses_per_schedule": SIZE,
            "distinct_forward_schedules_checked": len(checked),
            "distinct_inverse_schedules_checked": len(inverse_checked),
            "covers_every_positive_continuing_flag": True,
            "padding_and_zero_amplitude_rows_retained": True}


def resource_certificate():
    """Counts per core/law, including stops, active recovery and native words."""
    p = parameters()
    rows = []
    maximum = 0
    initial_encoding = fixed._counts(
        record_code._permutation_resources(encoder_permutation(INITIAL)),
        fixed._mixing_resources(len(supports(INITIAL)[0])))
    for center in LIVE:
        laws = {}
        for law in ("C", "E"):
            mean, peak = F(), 0
            for seam in SEAMS:
                native = fixed.repair_resources(law, seam)["native_pulses"]
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        probability = branch_probability(center, law, seam, coin, raw, p)
                        if not probability:
                            continue
                        cost = native+recovery_resources(center, law, seam, coin, raw)["native_pulses"]
                        mean += probability*cost/30
                        peak = max(peak, cost)
            maximum = max(maximum, peak)
            laws[law] = {"uniform_seam_expected_source_pulses_per_attempt": str(mean),
                         "maximum_source_pulses_one_positive_flag": peak,
                         "fixed_format_archive_bits_per_attempt": 1+5+1+circuit.RECORD_BITS[law]+7+7+1}
        rows.append({"core": list(center), "laws": laws})
    full = kernel(p)
    q = sp.Matrix([[sp.Rational(row[CORES.index(c)]) for c in LIVE] for row in full])
    mean_cost = sp.Matrix([sum(sp.Rational(row["laws"][law]["uniform_seam_expected_source_pulses_per_attempt"])
                                 for law in ("C", "E"))/2 for row in rows])
    expected = (sp.eye(60)-q).inv()*mean_cost
    if (sp.eye(60)-q)*expected != mean_cost:
        raise ArithmeticError("expected complete stopping-work equations failed")
    start = LIVE.index(INITIAL)
    return {"initial_encoder": initial_encoding, "rows": rows,
            "equal_C_E_expected_source_pulses_until_and_including_stop": [str(F(v)) for v in expected],
            "initial_expected_source_pulses_including_encoder_and_stopping_attempt": str(F(expected[start])+initial_encoding["native_pulses"]),
            "maximum_source_pulses_per_attempt_all_live_cores_and_laws": maximum,
            "archive_format": "law1, seam5, coin1, raw2or3, old_core7, new_core7, stop1",
            "initial_core_label_bits": 7, "live_source_register_qubits": 9,
            "native_record_register_peak_qubits": 3,
            "additional_six_level_source_processor_supplied": True,
            "physical_cost_of_load_export_storage_and_classical_control_established": False,
            "cache_reuse_is_charged_per_invocation": True,
            "units": "source pulses and classical bits; no energy or duration conversion"}


def witness_path():
    """A C then E path with changing live core, followed by an explicit stop."""
    for seam in SEAMS:
        for coin in (0, 1):
            for raw in range(1 << circuit.RECORD_BITS["C"]):
                images, target, _ = branch(INITIAL, "C", seam, coin, raw)
                if not images[0] or target not in LIVE or len(supports(target)[0]) != 2:
                    continue
                for second_seam in SEAMS:
                    for second_coin in (0, 1):
                        for second_raw in range(1 << circuit.RECORD_BITS["E"]):
                            second, stop, _ = branch(target, "E", second_seam, second_coin, second_raw)
                            if second[0] and stop in STOP:
                                return ((INITIAL, "C", seam, coin, raw),
                                        (target, "E", second_seam, second_coin, second_raw))
    raise ArithmeticError("no two-event live-motion/stop witness")


def source_replay_diagnostics():
    """Replay all60 encoders and every raw flag on a changing-core C/E path.

    This finite numerical gate is separate from the exact all43,200-slot
    identities and all-address schedule check. No full43,200-case floating
    replay or physical processor-failure certificate is asserted.
    """
    worst_encoder = 0.
    for center in LIVE:
        error = float(np.linalg.norm(encode(center, np.eye(2))-ideal_code(center), 2))
        if not math.isfinite(error):
            raise ArithmeticError("nonfinite source encoder")
        worst_encoder = max(worst_encoder, error)
    # A unitary logical basis makes these comparisons operator norms on the
    # entire logical space. Dividing its vectorization by sqrt(2) is a
    # complex maximally entangled state with an untouched reference.
    logical = np.array([[1, 1j], [1j, 1]], complex)/math.sqrt(2)
    state = encode(INITIAL, logical)
    weight = F(1)
    rows = []
    worst_branch, worst_mass = 0., 0.
    for center, law, seam, selected_coin, selected_raw in witness_path():
        outputs = source_branches(center, law, seam, state)
        total = sum(float(np.linalg.norm(row["amplitude"])**2) for row in outputs)
        mass_error = abs(total-float(np.linalg.norm(state)**2))
        if not math.isfinite(mass_error):
            raise ArithmeticError("nonfinite complete raw output mass")
        worst_mass = max(worst_mass, mass_error)
        for row in outputs:
            images, target, _ = branch(center, law, seam, row["coin"], row["raw"])
            if not images[0]:
                expected = np.zeros_like(state)
            elif target in LIVE:
                expected = math.sqrt(float(weight*row["probability"]))*(ideal_code(target)@logical)
            else:
                expected = np.zeros_like(state)
                for j, addresses in enumerate(images):
                    for address in addresses:
                        expected[address] = math.sqrt(float(weight*local.coin_probability(parameters(), seam, row["coin"]))/len(supports(center)[0]))*logical[j]
            error = float(np.linalg.norm(row["amplitude"]-expected, 2))
            if not math.isfinite(error):
                raise ArithmeticError("nonfinite source branch comparison")
            worst_branch = max(worst_branch, error)
        chosen = next(row for row in outputs if (row["coin"],row["raw"]) == (selected_coin,selected_raw))
        state = chosen["amplitude"]
        weight *= chosen["probability"]
        rows.append({"core_before": list(center), "law": law, "seam": list(seam),
                     "coin": selected_coin, "raw": selected_raw,
                     "core_after": list(chosen["core_after"]), "status": chosen["status"],
                     "retained_raw_slots": len(outputs), "cumulative_flag_weight_given_selected_seams": str(weight)})
    return {"all_sixty_encoder_max_operator_error": worst_encoder,
            "complete_path_max_branch_operator_error": worst_branch,
            "complete_path_max_total_mass_error": worst_mass,
            "changing_core_C_E_path": rows}

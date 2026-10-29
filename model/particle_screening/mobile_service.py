"""Source-word transport of the sharp-q=1 mobile code family.

The native attempted repair is followed by an explicitly supplied coherent
recovery.  Legal core moves are retained.  A forbidden antipodal core is
reanchored to the addressed seam, with both labels reconstructible from the
complete native tape and initial core.  Logical
transport is the identity in the declared equivariant frames.  This is an
actively serviced quantum memory, not autonomous motion or particle statistics.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
import math

import numpy as np

from . import mobile
from ..common_source import native_routing as routing, record_code
from ..observer import record_recovery as records
from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_discrimination import experiment as recorded
from ..observer.source_scalar import compiler
from ..screened_recovery import protocol as fixed


SIZE = 512
DIMENSION = 2
SEAMS = tuple(records.seams())
PROTOTYPES = ((0, 1), (0, 2))


def _check(value, name, ceiling=None):
    if not math.isfinite(value) or value < 0 or (ceiling is not None and value > ceiling):
        raise ArithmeticError(f"invalid mobile source diagnostic: {name}")
    return value


def _count(core):
    return len(mobile.supports(core)[0])


def encoder_permutation(core):
    core = mobile._legal(core)
    _validate_address_access()
    return _encoder_permutation(core)


@lru_cache(maxsize=60)
def _encoder_permutation(core):
    return routing.complete_permutation(tuple(
        (logical+8*term, circuit.encode(tuple(int(p in pair) for p in range(12))))
        for logical, pairs in enumerate(mobile.supports(core))
        for term, pair in enumerate(pairs)))


def branch_permutation(core, law, seam, coin, raw):
    args = mobile._label(core, law, seam, coin, raw)
    _validate_address_access()
    return _branch_permutation(*args)


@lru_cache(maxsize=None)
def _branch_permutation(core, law, seam, coin, raw):
    return routing.complete_permutation(tuple(
        (logical+8*term, address)
        for logical, row in enumerate(mobile.branch_images(core, law, seam, coin, raw))
        for term, address in enumerate(row)))


def _encode(core, value):
    return record_code._apply_permutation(encoder_permutation(core),
        compiler.execute(fixed._mixing_word(_count(core)), 9, value))


def code_matrix(core):
    """Actually execute the source words on both unknown logical columns."""
    core = mobile._legal(core)
    _validate_address_access()
    # Read-only array data does not freeze ndarray shape/stride metadata.
    # Detach the public result, including its bytes, so neither it nor its
    # .base chain exposes the source-amplitude cache to caller mutation.
    return routing._immutable(_code_matrix(core))


@lru_cache(maxsize=60)
def _code_matrix(core):
    return routing._immutable(_encode(core, np.eye(SIZE, DIMENSION, dtype=complex)))


def ideal_code_matrix(core):
    result = np.zeros((SIZE, DIMENSION), complex)
    for logical, pairs in enumerate(mobile.supports(core)):
        for pair in pairs:
            result[circuit.encode(tuple(int(p in pair) for p in range(12))), logical] = (
                1/math.sqrt(len(pairs)))
    return result


def encoding(core, logical):
    routing._original_scalars(logical)
    value = np.asarray(logical)
    if (value.ndim not in (1, 2) or value.shape[0] != DIMENSION
            or (value.ndim == 2 and value.shape[1] < 1)
            or value.dtype.kind not in "iufc" or not np.isfinite(value).all()):
        raise ValueError("two finite logical amplitudes or arbitrary reference columns required")
    return routing._immutable(code_matrix(core)@np.asarray(value, complex))


def decoding(core, state):
    value = record_code._apply_permutation(encoder_permutation(core), routing._state(state), True)
    return routing._immutable(compiler.execute(fixed._mixing_word(_count(core), True), 9, value))


def apply_recovery(core, law, seam, coin, raw, state):
    """All512 unitary extension; no conditional renormalization or projection."""
    status = mobile.transition(core, law, seam, coin, raw)
    value = routing._state(state)
    if not status["admitted"]:
        raise ValueError("raw abort flags terminate instead of selecting recovery")
    rows = mobile.branch_images(core, law, seam, coin, raw)
    count = len(rows[0])
    if not count:
        return routing._immutable(value)
    value = record_code._apply_permutation(branch_permutation(core, law, seam, coin, raw), value, True)
    value = compiler.execute(fixed._mixing_word(count, True), 9, value)
    return routing._immutable(_encode(tuple(status["target_core"]), value))


def repair_branches(core, law, seam, state):
    """The full source instrument with every zero, abort and recovery output."""
    mobile.transition(core, law, seam, 0, 0)  # Validate the complete address.
    _validate_address_access()
    route = record_code._route_permutation(seam)
    value = record_code._apply_permutation(route, routing._state(state))
    branches = recorded.repair_branches(law, value)
    _validate_native_grid(law, branches)
    for row in branches:
        coin, raw = row["coin"], row["record_code"]
        status = mobile.transition(core, law, seam, coin, raw)
        value = record_code._apply_permutation(route, row["amplitude"], True)
        native = routing._immutable(value)
        if status["admitted"]:
            value = apply_recovery(core, law, seam, coin, raw, value)
        row.update(amplitude=routing._immutable(value), native_amplitude=native,
                   core_before=core, **status,
                   branch_weight=mobile.branch_weight(core, law, seam, coin, raw),
                   recovery_resources=recovery_resources(core, law, seam, coin, raw))
    return branches


def _validate_native_grid(law, branches):
    """Zero weight is never permission to omit a raw source output."""
    expected = {(coin, raw) for coin in (0, 1)
                for raw in range(1 << circuit.RECORD_BITS[law])}
    if not isinstance(branches, (tuple, list)) or len(branches) != len(expected):
        raise ArithmeticError("incomplete native raw output grid")
    seen = set()
    for row in branches:
        if not isinstance(row, dict) or not {
                "coin", "record_code", "current", "contrast", "abort", "amplitude"} <= row.keys():
            raise ArithmeticError("invalid native raw output metadata")
        coin, raw = row["coin"], row["record_code"]
        if (type(coin) is not int or type(raw) is not int
                or (coin, raw) not in expected or (coin, raw) in seen):
            raise ArithmeticError("incomplete or duplicate native raw output grid")
        seen.add((coin, raw))
        value = circuit.record_values(law, coin).get(raw)
        current = value if law == "C" or value is None else (
            value-(2*coin-1)*(abs(value) % 2))//2
        decoded = {"abort": value is None, "current": current,
                   "contrast": value if law == "E" else None}
        if any(type(row[key]) is not type(item) or row[key] != item
               for key, item in decoded.items()):
            raise ArithmeticError("native raw output metadata disagrees with decoded record")
    if seen != expected:
        raise ArithmeticError("incomplete native raw output grid")


def _sum_counts(*rows):
    return {key: sum(row[key] for row in rows) for key in routing.COUNT_KEYS}


def encoding_resources(core):
    return _sum_counts(record_code._permutation_resources(encoder_permutation(core)),
                       fixed._mixing_resources(_count(core)))


def recovery_resources(core, law, seam, coin, raw):
    status = mobile.transition(core, law, seam, coin, raw)
    count = len(mobile.branch_images(core, law, seam, coin, raw)[0])
    if not count or not status["admitted"]:
        return dict.fromkeys(routing.COUNT_KEYS, 0)
    return _sum_counts(
        record_code._permutation_resources(branch_permutation(core, law, seam, coin, raw)),
        fixed._mixing_resources(count, True), encoding_resources(tuple(status["target_core"])))


def _addresses():
    for core in mobile.legal_cores():
        for law in ("C", "E"):
            for seam in SEAMS:
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        yield core, law, seam, coin, raw


def _integer_permutation(images, inverse=False):
    """Execute the actual emitted Gray-flip word on every integer address."""
    word = routing.permutation_edges(images, inverse=inverse)
    return _checked_integer_word(images, inverse, word)


@lru_cache(maxsize=None)
def _checked_integer_word(images, inverse, word):
    # Cache by the emitted word, so a changed compiler schedule cannot reuse
    # an admission obtained for an earlier word with the same target map.
    answer = record_code._integer_schedule_action(word)
    expected = tuple(images.index(i) for i in range(SIZE)) if inverse else images
    if tuple(answer) != expected:
        raise ArithmeticError("mobile emitted full512 permutation differs from its declared map")
    return tuple(answer)


def compiled_layout_certificate():
    address_contract = primitive_address_certificate()
    inherited = record_code._compiled_layout_certificate()
    if inherited["native_routes_checked_on_all78_states_and_padding"] != 30:
        raise ArithmeticError("native full-carrier routing contract changed")
    encoder_slots = image_slots = positive = 0
    for core in mobile.legal_cores():
        forward = _integer_permutation(encoder_permutation(core))
        inverse = _integer_permutation(encoder_permutation(core), True)
        if any(inverse[forward[i]] != i for i in range(SIZE)):
            raise ArithmeticError("mobile encoder inverse failed")
        for logical, pairs in enumerate(mobile.supports(core)):
            for term, pair in enumerate(pairs):
                if forward[logical+8*term] != circuit.encode(tuple(int(p in pair) for p in range(12))):
                    raise ArithmeticError("mobile encoder altered a logical support")
                encoder_slots += 1
    for args in _addresses():
        rows = mobile.branch_images(*args)
        if not rows[0]:
            continue
        decoded = _integer_permutation(branch_permutation(*args), True)
        for logical, row in enumerate(rows):
            for term, address in enumerate(row):
                if decoded[address] != logical+8*term:
                    raise ArithmeticError("mobile source recovery changed a branch logical coordinate")
                image_slots += 1
        positive += 1
    return {"primitive_address_contract": address_contract,
            "core_encoders": len(mobile.legal_cores()), "encoder_slots": encoder_slots,
            "positive_branch_inverse_schedules": positive, "native_image_slots": image_slots,
            "source_register_dimension": SIZE, "native_routes_full78_and_padding": 30,
            "each_distinct_emitted_permutation_checked_on_all512_addresses": True,
            "exact_identity": "R[c,h] K[h] V[c] = sqrt(p[c,h]) V[target(c,h)]",
            "logical_holonomy_on_every_closed_corrected_core_path": "I2",
            "arbitrary_reference_preserved_under_ideal_services": True,
            "quantum_exchange_phase_or_particle_identity": False}


def primitive_address_certificate():
    """Independent total-two census and literal source register contract.

    This refuses a jointly mutated geometry/encoder interpretation. The
    native carrier is rebuilt from unordered port multisets, and source
    addresses from the documented endpoint and spectator register layout.
    These targets never enter numerical source-word execution.
    """
    states = tuple(sorted(tuple(pair.count(p) for p in range(12))
                          for pair in combinations_with_replacement(range(12), 2)))
    endpoints = ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0))
    if tuple(sorted(native.states())) != states or circuit.SEAM_STATES != endpoints:
        raise ArithmeticError("primitive total-two carrier/address contract changed")
    raw_values = {("C", 0): {0: 0, 1: 1, 2: -1},
                  ("C", 1): {0: 0, 1: 1, 2: -1},
                  ("E", 0): {0: 0, 1: 2, 2: -2, 4: -1, 5: 1},
                  ("E", 1): {0: 0, 1: 2, 2: -2, 4: 1, 6: -1}}
    if (circuit.RECORD_BITS != {"C": 2, "E": 3}
            or any(circuit.record_values(*key) != value for key, value in raw_values.items())):
        raise ArithmeticError("primitive complete native raw-record contract changed")
    shells = []
    for total in range(3):
        entries = sorted({state[2:] for state in states if sum(state[2:]) == total})
        if total:
            first = ((1,)+(0,)*9) if total == 1 else ((1, 1)+(0,)*8)
            entries.remove(first)
            entries.insert(0, first)
        shells.append(entries)
    for state in states:
        address = endpoints.index(state[:2])+8*shells[sum(state[2:])].index(state[2:])
        if circuit.encode(state) != address or circuit.decode(address) != state:
            raise ArithmeticError("primitive source address differs from independent register layout")
    return {"independently_enumerated_total_two_states": 78,
            "binary_states": sum(max(state) == 1 for state in states),
            "each_source_encode_decode_checked_against_independent_layout": True}


def _validate_address_access():
    # A warm source-amplitude cache may never bypass changed public geometry
    # or source encoding. The token includes all globals consumed by encode /
    # decode, their current implementations, and a fresh carrier value.
    token = (circuit.encode, getattr(circuit.encode, "__code__", None),
             circuit.decode, getattr(circuit.decode, "__code__", None),
             circuit.record_values, getattr(circuit.record_values, "__code__", None),
             tuple(sorted(circuit.RECORD_BITS.items())),
             circuit.SEAM_STATES, circuit.SPECTATOR_SHELLS, tuple(native.states()))
    _validated_address_token(token)


@lru_cache(maxsize=8)
def _validated_address_token(token):
    primitive_address_certificate()


def resource_report():
    """Charge every emitted word even if its numerical amplitudes were cached."""
    laws = {}
    for law in ("C", "E"):
        stationary_counts = dict.fromkeys(routing.COUNT_KEYS, F())
        stationary_recovery_counts = dict.fromkeys(routing.COUNT_KEYS, F())
        maximum_counts = dict.fromkeys(routing.COUNT_KEYS, 0)
        maximum_recovery_counts = dict.fromkeys(routing.COUNT_KEYS, 0)
        per_core = []
        for core in mobile.legal_cores():
            total = dict.fromkeys(routing.COUNT_KEYS, F())
            recovery = dict.fromkeys(routing.COUNT_KEYS, F())
            for seam in SEAMS:
                native = record_code.repair_resources(law, seam)
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        weight = mobile.branch_weight(core, law, seam, coin, raw)/len(SEAMS)
                        if not weight:
                            continue
                        cost = recovery_resources(core, law, seam, coin, raw)
                        for key in routing.COUNT_KEYS:
                            total[key] += weight*(native[key]+cost[key])
                            recovery[key] += weight*cost[key]
                            maximum_counts[key] = max(maximum_counts[key], native[key]+cost[key])
                            maximum_recovery_counts[key] = max(maximum_recovery_counts[key], cost[key])
            # The exact native seam-reset controller has class masses3/5,2/5;
            # each orbit has30 cores. This is a supplied uniform scheduler.
            stationary = F(1, 50) if core in SEAMS else F(1, 75)
            for key in routing.COUNT_KEYS:
                stationary_counts[key] += stationary*total[key]
                stationary_recovery_counts[key] += stationary*recovery[key]
            per_core.append({"core": list(core),
                             "mean_native_plus_recovery_pulses": str(total["native_pulses"]),
                             "mean_native_plus_recovery_resources": {key: str(value) for key, value in total.items()}})
        laws[law] = {
            "stationary_mean_native_plus_recovery_pulses": str(stationary_counts["native_pulses"]),
            "stationary_mean_recovery_pulses": str(stationary_recovery_counts["native_pulses"]),
            "stationary_mean_native_plus_recovery_resources": {key: str(value) for key, value in stationary_counts.items()},
            "stationary_mean_recovery_resources": {key: str(value) for key, value in stationary_recovery_counts.items()},
            "maximum_one_attempt_native_plus_recovery_pulses": maximum_counts["native_pulses"],
            "maximum_recovery_pulses": maximum_recovery_counts["native_pulses"],
            "maximum_by_resource_key_over_all_branches": maximum_counts,
            "maximum_recovery_by_resource_key_over_all_branches": maximum_recovery_counts,
            "per_core_uniform_seam_means": per_core,
            "native_tape_bits_per_attempt": 5+1+circuit.RECORD_BITS[law],
            "optional_debug_format_bits_per_attempt": 5+1+circuit.RECORD_BITS[law]+15,
            "native_blank_preparations_per_attempt": {"coin_qubits": 1,
                                                       "record_qubits": circuit.RECORD_BITS[law]},
            "seam_selection_service_per_attempt": 1,
            "native_coin_and_record_readout_bits_per_attempt": 1+circuit.RECORD_BITS[law],
            "native_ancilla_reuse_requires_blank_preparation_each_attempt": True,
            "stationary_core_reset_rate": "1/150"}
    return {"encoders": [{"core": list(c), **encoding_resources(c)} for c in mobile.legal_cores()],
            "initial_logical_qubits": 1, "initial_blank_code_qubits": 8,
            "additional_recovery_scratch_qubits": 0, "initial_classical_core_bits": 7,
            "archive_format": "seam5, coin1, raw2/3; both core labels and reset status reconstructed",
            "optional_debug_format": "add reset-status1, native-core7, target-core7; not required by the controller",
            "classical_reset_service_has_no_assigned_physical_work": True,
            "native_blank_reset_readout_storage_services_have_no_assigned_physical_work": True,
            "counts_are_not_energy_or_time": True, "laws": laws}


def source_composition_report():
    """Fresh floating product bounds covering every emitted mobile word.

    The prototype packet supplies independently executed full512 macro and
    local dilation diagnostics. We count the actual new mobile schedules,
    all mixing sizes and all retained orthogonal classical control blocks.
    These are numerical bounds, not outward rounding or hardware guarantees.
    """
    inherited = fixed.source_composition_report()
    macro = {False: inherited["forward_full512_macro_operator_error_plus_polar_change"],
             True: inherited["inverse_full512_macro_operator_error_plus_polar_change"]}
    xerror = inherited["source_X_conjugation_operator_error"]
    mixing = {}
    for count in (1, 2, 3, 4):
        synthesis = inherited["mixing_angle_synthesis_diagnostics"][str(count)]
        for inverse in (False, True):
            mixing[count, inverse] = (record_code._word_log_error(fixed._mixing_word(count, inverse))
                +math.log1p(synthesis["uniform_target_unitary_completion_error_upper"]))
    def permutation_log(images, inverse=False):
        edges = routing.permutation_edges(images, inverse=inverse)
        return (len(edges)*math.log1p(macro[inverse])
                +2*sum(pattern.bit_count() for _, pattern in edges)*math.log1p(xerror))
    encoders = {c: permutation_log(encoder_permutation(c))+mixing[_count(c), False]
                for c in mobile.legal_cores()}
    decoders = {c: permutation_log(encoder_permutation(c), True)+mixing[_count(c), True]
                for c in mobile.legal_cores()}
    laws = {}
    for law in ("C", "E"):
        native = inherited["laws"][law]["complete_local_native_dilation_error"]
        prefix = {seam: (permutation_log(record_code._route_permutation(seam))
                    +permutation_log(record_code._route_permutation(seam), True)
                    +math.log1p(inherited["fair_coin_preparation_operator_error"])
                    +math.log1p(native)) for seam in SEAMS}
        worst = 0.
        for args in _addresses():
            core, branch_law, seam, coin, raw = args
            if branch_law != law:
                continue
            count = len(mobile.branch_images(*args)[0])
            recovery = 0.
            if count:
                target = tuple(mobile.transition(*args)["target_core"])
                recovery = (permutation_log(branch_permutation(*args), True)
                            +mixing[count, True]+encoders[target])
            worst = max(worst, _check(prefix[seam]+recovery, "attempt log"))
        bounds = {str(t): _check(math.expm1(max(encoders.values())+t*worst),
                                f"{law} {t} attempts", 1e-5) for t in (1, 2, 10, 100)}
        laws[law] = {"maximum_attempt_log_error": worst,
                     "encoded_t_attempt_operator_error_upper": bounds,
                     "all_sixty_cores_thirty_seams_and_raw_slots_counted": True}
    return {"maximum_encoder_operator_error": _check(math.expm1(max(encoders.values())), "encode", 1e-8),
            "maximum_decoder_operator_error": _check(math.expm1(max(decoders.values())), "decode", 1e-8),
            "method": "product(1+factor error)-1; maximum over orthogonal classical flag blocks",
            "full512_macro_inverse_mixing_and_native_dilation_diagnostics_reused": True,
            "arbitrary_reference_has_same_operator_bound": True,
            "numerical_only_not_interval_hardware_or_classical_reset_certificate": True,
            "laws": laws}


def _positive_for_target(core, law, target, *, reset=False):
    for seam in SEAMS:
        for coin in (0, 1):
            for raw in range(1 << circuit.RECORD_BITS[law]):
                status = mobile.transition(core, law, seam, coin, raw)
                if (status["positive"] and status["reset_performed"] == reset
                        and (target is None or tuple(status["native_core"]) == target)):
                    return law, seam, coin, raw
    raise ArithmeticError("the requested positive mobile path is unavailable")


def witness_schedules():
    """Source incidence selects one open/reset path and one closed core loop."""
    start = PROTOTYPES[0]
    distance_two = next(c for c in mobile.legal_cores() if c not in SEAMS
        and any(mobile.transition(start, "C", seam, coin, raw)["positive"]
                and tuple(mobile.transition(start, "C", seam, coin, raw)["native_core"]) == c
                for seam in SEAMS for coin in (0, 1) for raw in range(4)))
    moving = _positive_for_target(start, "C", distance_two)
    reset = _positive_for_target(distance_two, "E", None, reset=True)
    # Three cores in a native triangular face exchange the two core labels.
    a, b = start
    third = next(v for v in range(12) if tuple(sorted((a, v))) in SEAMS
                 and tuple(sorted((b, v))) in SEAMS)
    loop = (start, tuple(sorted((b, third))), tuple(sorted((a, third))), start)
    word = tuple(_positive_for_target(c, law, d)
                 for c, d, law in zip(loop[:-1], loop[1:], ("C", "E", "C"), strict=True))
    return {"motion_then_forbidden_reset": (moving, reset), "closed_core_exchange": word}


def execute_stopped_history(schedule, start=PROTOTYPES[0]):
    """Complete stopped instrument: every nonselected outcome terminates.

    Each prefix executes all raw source outputs and recovers every admitted
    output. One declared flag continues; all complements, zero slots and
    aborts are terminal unnormalized quantum outputs. Nothing is restarted.
    The two logical columns prove the same linear identity with any reference.
    """
    if (type(schedule) is not tuple or not schedule
            or any(type(entry) is not tuple or len(entry) != 4 for entry in schedule)):
        raise ValueError("a nonempty tuple of complete native continuation addresses required")
    value, core, weight = code_matrix(start), start, F(1)
    terminal, steps = [], []
    for index, (law, seam, selected_coin, selected_raw) in enumerate(schedule):
        mobile._label(core, law, seam, selected_coin, selected_raw)
        branches = repair_branches(core, law, seam, value)
        incoming = weight
        total = sum(float(np.linalg.norm(row["amplitude"])**2)/2 for row in branches)
        _check(abs(total-float(incoming)), "complete prefix output mass", 1e-8)
        for row in branches:
            branch_weight = incoming*row["branch_weight"]
            is_selected = (row["coin"], row["record_code"]) == (selected_coin, selected_raw)
            target = tuple(row["target_core"]) if row["admitted"] else core
            expected = math.sqrt(float(branch_weight))*ideal_code_matrix(target)
            error = _check(float(np.linalg.norm(row["amplitude"]-expected, 2)), "complete source branch", 1e-8)
            row.update(cumulative_weight=branch_weight, prefix_length=index+1,
                       selected_for_continuation=is_selected and index+1 < len(schedule),
                       source_operator_error=error)
            if not row["selected_for_continuation"]:
                terminal.append(row)
        chosen = next(row for row in branches if (row["coin"], row["record_code"]) == (selected_coin, selected_raw))
        if not chosen["positive"] or not chosen["admitted"]:
            raise ValueError("a continuing prefix must select a positive admitted flag")
        steps.append({"core_before": list(core), "native_core": list(chosen["native_core"]),
                      "target_core": list(chosen["target_core"]), "law": law,
                      "seam": list(seam), "coin": selected_coin, "raw": selected_raw,
                      "reset_performed": chosen["reset_performed"],
                      "retained_raw_slots": 2*(1 << circuit.RECORD_BITS[law]),
                      "selected_cumulative_probability": str(chosen["cumulative_weight"])})
        value, core, weight = chosen["amplitude"], tuple(chosen["target_core"]), chosen["cumulative_weight"]
    exact_total = sum((row["cumulative_weight"] for row in terminal), F())
    expected_terminal_slots = sum(2*(1 << circuit.RECORD_BITS[law])
                                  for law, _, _, _ in schedule)-len(schedule)+1
    if len(terminal) != expected_terminal_slots:
        raise ArithmeticError("incomplete stopped native raw output grid")
    if exact_total != 1:
        raise ArithmeticError("stopped history silently discarded a complementary output")
    decoded = decoding(core, value)
    target = np.zeros((SIZE, DIMENSION), complex)
    target[:2] = math.sqrt(float(weight))*np.eye(DIMENSION)
    error = _check(float(np.linalg.norm(decoded-target, 2)), "terminal decoder", 1e-8)
    return {"steps": steps, "terminal": terminal, "selected_output": value,
            "selected_decoded_output": decoded, "terminal_decoder_error": error,
            "terminal_exact_probability_sum": str(exact_total)}


def source_replay_report():
    primitive_address_certificate()
    prototype = []
    for core in PROTOTYPES:
        error = _check(float(np.linalg.norm(code_matrix(core)-ideal_code_matrix(core), 2)),
                       "prototype encoder", 1e-8)
        prototype.append({"core": list(core), "source_encoder_error": error})
    histories = {}
    for name, schedule in witness_schedules().items():
        result = execute_stopped_history(schedule)
        histories[name] = {"steps": result["steps"], "terminal_outputs": len(result["terminal"]),
            "terminal_exact_probability_sum": result["terminal_exact_probability_sum"],
            "maximum_terminal_source_operator_error": max(row["source_operator_error"] for row in result["terminal"]),
            "selected_terminal_decoder_error": result["terminal_decoder_error"],
            "all512_source_rows_and_untouched_reference_columns_retained": True,
            "complementary_outputs_terminate_without_restart": True}
    return {"source_executed_prototypes": prototype, "complete_stopped_histories": histories,
            "unexecuted_full_allcore_tree_replaced_by_exact_branch_composition": True,
            "ordinary_native_history_rank_one_exchange_theorem_not_applied_across_recovery": True,
            "closed_corrected_core_exchange_logical_holonomy": "I2",
            "physical_spin_statistics_or_particle_identification": False}


def report():
    """Fresh result for the parent packet; custody belongs to its one receipt."""
    return {"compiled_layout": compiled_layout_certificate(),
            "resources": resource_report(),
            "source_composition": source_composition_report(),
            "source_replay": source_replay_report()}

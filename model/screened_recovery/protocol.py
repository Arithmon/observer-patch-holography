"""Complete flagged C/E service for the screened sharp-q=1 qubit.

The source encoder and conditional recovery are emitted through the existing
nine-qubit source compiler.  Each branch acts on all 512 source rows.  The
classical core label is reset after recovery, with both labels archived.  The
finite exact statement is about actively serviced memory at native attempt
count, not an autonomous particle or a calibrated lifetime.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from ..common_source import native_routing as routing, record_code
from ..observer import record_recovery as records
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from ..particle_screening import screening
from ..exchange_sector import sector


CORE = screening.CORE
SUPPORTS = screening.sharp_code(CORE)
SEAMS = tuple(records.seams())
SIZE = 512
DIMENSION = 2
assert SUPPORTS is not None and all(len(row) == 4 for row in SUPPORTS)


def _label(law, seam, coin, raw):
    law = routing._candidate(law)
    if type(seam) is not tuple or seam not in SEAMS:
        raise ValueError("an ordered seam from the complete native menu required")
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("a native binary coin required")
    if type(raw) is not int or not 0 <= raw < 1 << circuit.RECORD_BITS[law]:
        raise ValueError("a raw record slot from the complete menu required")
    return law, seam, coin, raw


@lru_cache(maxsize=None)
def branch_images(law, seam, coin, raw):
    """Independently enumerate native endpoint maps and their raw record."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    values = circuit.record_values(law, coin)
    rows = []
    for pairs in SUPPORTS:
        row = []
        for pair in pairs:
            state = tuple(int(port in pair) for port in range(12))
            output, record = records.recorded_step(law, state, seam, coin)
            if values.get(raw) == record[-1]:
                row.append(circuit.encode(output))
        rows.append(tuple(row))
    if len(rows[0]) != len(rows[1]) or len(set(rows[0]+rows[1])) != sum(map(len, rows)):
        raise ArithmeticError("screened code failed a scalar, injective raw branch")
    return tuple(rows)


def branch_weight(law, seam, coin, raw):
    """Exact probability of a raw branch, including its fair coin."""
    return F(len(branch_images(law, seam, coin, raw)[0]), 8)


@lru_cache(maxsize=1)
def _encoder_permutation():
    return routing.complete_permutation(tuple(
        (logical+8*term, circuit.encode(tuple(int(port in pair) for port in range(12))))
        for logical, pairs in enumerate(SUPPORTS) for term, pair in enumerate(pairs)))


@lru_cache(maxsize=None)
def _branch_permutation(law, seam, coin, raw):
    return routing.complete_permutation(tuple(
        (logical+8*term, address)
        for logical, row in enumerate(branch_images(law, seam, coin, raw))
        for term, address in enumerate(row)))


@lru_cache(maxsize=8)
def _mixing_word(count, inverse=False):
    if type(count) is not int or count not in (1, 2, 3, 4) or type(inverse) is not bool:
        raise ValueError("a support size 1..4 and explicit inverse flag required")
    if count == 1:
        word = []
    elif count == 2:
        word = [compiler.pulse("Y", -math.pi/4, [3])]
    elif count == 4:
        word = [compiler.pulse("Y", -math.pi/4, [3]),
                compiler.pulse("Y", -math.pi/4, [4])]
    else:
        word = compiler.remap(compiler.compile_real_state(
            [1/math.sqrt(3)]*3+[0.]*5), (3, 4, 5))
    return compiler.inverse(word) if inverse else word


def _encode(value):
    return record_code._apply_permutation(
        _encoder_permutation(), compiler.execute(_mixing_word(4), 9, value))


@lru_cache(maxsize=1)
def code_matrix():
    """Source-execute both logical columns, retaining all padding rows."""
    value = _encode(np.eye(SIZE, DIMENSION, dtype=complex))
    return routing._immutable(value)


def ideal_code_matrix():
    """Independent geometric target for comparison; never used to encode."""
    result = np.zeros((SIZE, DIMENSION), complex)
    for logical, pairs in enumerate(SUPPORTS):
        for pair in pairs:
            state = tuple(int(port in pair) for port in range(12))
            result[circuit.encode(state), logical] = .5
    return routing._immutable(result)


def encoding(state):
    """Encode arbitrary two-component logical/reference columns from blanks."""
    value = np.asarray(state)
    if (value.ndim not in (1, 2) or value.shape[0] != 2
            or (value.ndim == 2 and value.shape[1] < 1)
            or value.dtype.kind not in "iufc" or not np.isfinite(value).all()):
        raise ValueError("finite two-component logical amplitudes/reference columns required")
    return routing._immutable(code_matrix()@np.asarray(value, complex))


def decoding(state):
    value = routing._state(state)
    value = record_code._apply_permutation(_encoder_permutation(), value, True)
    return routing._immutable(compiler.execute(_mixing_word(4, True), 9, value))


def apply_recovery(law, seam, coin, raw, state):
    """Execute a record-selected unitary on all 512 source rows."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    value = routing._state(state)
    if raw not in circuit.record_values(law, coin):
        raise ValueError("raw abort record cannot select an active recovery")
    count = len(branch_images(law, seam, coin, raw)[0])
    if count == 0:
        return routing._immutable(value)
    value = record_code._apply_permutation(_branch_permutation(law, seam, coin, raw), value, True)
    value = compiler.execute(_mixing_word(count, True), 9, value)
    return routing._immutable(_encode(value))


def repair_branches(law, seam, state, *, recover=True):
    """Source-execute the native instrument, route, and all raw recoveries."""
    law = routing._candidate(law)
    _label(law, seam, 0, 0)
    if type(recover) is not bool:
        raise ValueError("recover must be an explicit Boolean")
    value = routing._state(state)
    route = record_code._route_permutation(seam)
    local = record_code._apply_permutation(route, value)
    branches = recorded.repair_branches(law, local)
    for branch in branches:
        coin, raw = branch["coin"], branch["record_code"]
        amplitude = record_code._apply_permutation(route, branch["amplitude"], True)
        applied = recover and not branch["abort"]
        if applied:
            amplitude = apply_recovery(law, seam, coin, raw, amplitude)
        branch.update(amplitude=routing._immutable(amplitude), seam=seam,
                      recovery_applied=applied,
                      code_branch_weight=branch_weight(law, seam, coin, raw),
                      recovery_resources=recovery_resources(law, seam, coin, raw)
                      if applied else dict.fromkeys(routing.COUNT_KEYS, 0))
    return branches


def _counts(*rows):
    return {key: sum(row[key] for row in rows) for key in routing.COUNT_KEYS}


def _mixing_resources(count, inverse=False):
    row = compiler.resources(_mixing_word(count, inverse))
    return {key: row[key] for key in routing.COUNT_KEYS}


def encoding_resources():
    row = _counts(record_code._permutation_resources(_encoder_permutation()),
                  _mixing_resources(4))
    return {**row, "supplied_logical_levels": 2, "supplied_register_qubits": 1,
            "blank_code_qubits": 8, "full_register_dimension": SIZE,
            "new_classical_records": 0}


def recovery_resources(law, seam, coin, raw):
    law, seam, coin, raw = _label(law, seam, coin, raw)
    count = len(branch_images(law, seam, coin, raw)[0])
    admitted = raw in circuit.record_values(law, coin)
    row = dict.fromkeys(routing.COUNT_KEYS, 0)
    if count and admitted:
        row = _counts(encoding_resources(),
                      record_code._permutation_resources(_branch_permutation(law, seam, coin, raw)),
                      _mixing_resources(count, True))
    return {**row, "code_branch_weight_exact": str(F(count, 8)),
            "raw_abort": not admitted, "identity_on_valid_zero_weight_branch": admitted and not count,
            "input_dependent_normalization": False, "archive_record_rewrites": 0}


def repair_resources(law, seam):
    return record_code.repair_resources(law, seam)


def core_after_branch(law, seam, coin, raw):
    """Classical co-moving core on an admitted positive branch if unique."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    if raw not in circuit.record_values(law, coin):
        raise ValueError("raw abort has no transported core")
    outputs = []
    for pairs in SUPPORTS:
        for pair in pairs:
            state = tuple(int(port in pair) for port in range(12))
            final, record = records.recorded_step(law, state, seam, coin)
            if circuit.record_values(law, coin).get(raw) == record[-1]:
                outputs.append((final, record))
    if not outputs:
        return None
    moved = {sector.moved_core(CORE, law, record) for _, record in outputs}
    if len(moved) != 1:
        raise ArithmeticError("transported core is not determined by the complete flag")
    return moved.pop()


def archive_entry(law, seam, coin, raw):
    """Fixed-format reset record; readout and reset are supplied services."""
    moved = core_after_branch(law, seam, coin, raw)
    if moved is None:
        return None
    return {"core_before": list(CORE), "core_after_native": list(moved),
            "core_after_recovery_and_reset": list(CORE),
            "reset_performed": moved != CORE,
            "archive_bits": 15,
            "format": "status + seven-bit old core label + seven-bit new core label"}


def exact_service_certificate():
    """Check all 720 raw slots and return exact classical-composition data."""
    counts, positive, weights, resets, unreset_q1, nonedge_core = {}, {}, {}, {}, {}, {}
    for law in ("C", "E"):
        slot_hist, active = Counter(), 0
        seam_weights, reset_probability, q1_probability, nonedge_probability = [], F(), F(), F()
        for seam in SEAMS:
            total = F()
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    rows = branch_images(law, seam, coin, raw)
                    count = len(rows[0])
                    slot_hist[count] += 1
                    p = branch_weight(law, seam, coin, raw)
                    total += p
                    if p:
                        active += 1
                        moved = core_after_branch(law, seam, coin, raw)
                        reset_probability += p*int(moved != CORE)/len(SEAMS)
                        if moved not in screening.SEAMS:
                            nonedge_probability += p/len(SEAMS)
                        q1_columns = tuple(F(sum(sector.charge(
                            tuple(int(port in pair) for port in range(12)), moved) == 1
                            for pair in pairs), 4) for pairs in SUPPORTS)
                        if q1_columns[0] != q1_columns[1]:
                            raise ArithmeticError("unreset q=1 probability depends on logical input")
                        q1_probability += p*q1_columns[0]/len(SEAMS)
            if total != 1:
                raise ArithmeticError("raw branch weights do not sum to one")
            seam_weights.append(str(total))
        counts[law], positive[law] = dict(sorted(slot_hist.items())), active
        weights[law] = seam_weights
        resets[law] = str(reset_probability)
        unreset_q1[law] = str(q1_probability)
        nonedge_core[law] = str(nonedge_probability)
    if counts != {"C": {0: 150, 1: 20, 2: 20, 3: 20, 4: 30},
                  "E": {0: 362, 1: 40, 2: 46, 3: 20, 4: 12}}:
        raise ArithmeticError("complete screened branch census changed")
    if resets != {"C": "1/15", "E": "1/15"}:
        raise ArithmeticError("core reset frequency differs from prior screening")
    if unreset_q1 != {"C": "29/30", "E": "29/30"}:
        raise ArithmeticError("unreset fixed-code charge differs from prior screening")
    if nonedge_core != {"C": "1/30", "E": "1/30"}:
        raise ArithmeticError("core nonedge transport differs from prior screening")
    return {
        "code_supports": [[list(pair) for pair in row] for row in SUPPORTS],
        "raw_slots_checked": 720,
        "positive_slots_by_law": positive,
        "support_count_census_by_law": {law: {str(k): v for k, v in row.items()}
                                        for law, row in counts.items()},
        "raw_probability_sum_per_seam": weights,
        "uniform_seam_core_reset_probability_by_law": resets,
        "uniform_seam_core_nonedge_probability_by_law": nonedge_core,
        "uniform_seam_q1_probability_after_recovery_without_label_reset": unreset_q1,
        "q1_probability_after_label_reset": "1",
        "exact_branch_identity": "R_h K_h V = sqrt(p_h) V; zero and abort slots have p_h=0",
        "arbitrary_unknown_reference_preserved": True,
        "record_adaptive_composition": "For any t native attempts with a history-controlled C/E seam choice, provided each positive branch receives its selected unitary recovery and archived core reset, the corrected carrier/reference is unchanged; flag probabilities are logical-state-independent.",
        "active_service_required_after_each_attempt": True,
        "native_attempt_count_horizon_under_ideal_services": "all finite t",
        "physical_time_or_fault_tolerance_bound": False,
    }


def compiled_layout_certificate():
    """Execute emitted integer Gray schedules on every source address.

    This exact check is independent of floating gate execution.  It checks
    the actual controlled-flip schedule, not only the requested permutation.
    The inherited all-78 route certificate supplies the native endpoint and
    padding check for each of the thirty physical seam routes.
    """
    inherited = record_code._compiled_layout_certificate()
    if inherited["native_routes_checked_on_all78_states_and_padding"] != 30:
        raise ArithmeticError("the inherited native route certificate changed")
    permutation = _encoder_permutation()
    forward = record_code._integer_schedule_action(routing.permutation_edges(permutation))
    inverse = record_code._integer_schedule_action(
        routing.permutation_edges(permutation, inverse=True))
    if any(inverse[forward[i]] != i for i in range(SIZE)):
        raise ArithmeticError("emitted encoder inverse is not exact on all 512 addresses")
    slots = 0
    for logical, pairs in enumerate(SUPPORTS):
        for term, pair in enumerate(pairs):
            state = tuple(int(port in pair) for port in range(12))
            if forward[logical+8*term] != circuit.encode(state):
                raise ArithmeticError("emitted encoder altered a declared code slot")
            slots += 1
    positive, images = 0, 0
    for law in ("C", "E"):
        for seam in SEAMS:
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    rows = branch_images(law, seam, coin, raw)
                    if not rows[0]:
                        continue
                    permutation = _branch_permutation(law, seam, coin, raw)
                    decoded = record_code._integer_schedule_action(
                        routing.permutation_edges(permutation, inverse=True))
                    for logical, row in enumerate(rows):
                        for term, address in enumerate(row):
                            if decoded[address] != logical+8*term:
                                raise ArithmeticError("emitted recovery changed a native image slot")
                            images += 1
                    positive += 1
    if positive != 208 or slots != 8:
        raise ArithmeticError("a compiled screened branch or encoder slot is missing")
    return {
        "encoder_schedule_and_inverse_checked_on_all512_addresses": True,
        "declared_encoder_slots_checked": slots,
        "native_routes_reused_from_exact_all78_and_padding_certificate": 30,
        "positive_branch_inverse_schedules_checked": positive,
        "native_branch_image_slots_checked": images,
        "source_register_dimension": SIZE,
        "no_target_permutation_substituted_for_execution": True,
    }


def resource_report():
    """Exact source pulse and classical-tape counts for the complete menu."""
    rows = {}
    for law in ("C", "E"):
        pulse_sum, recovery_sum = F(), F()
        maximum_total, maximum_recovery = 0, 0
        zero = dict.fromkeys(routing.COUNT_KEYS, 0)
        for seam in SEAMS:
            native = repair_resources(law, seam)
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    p = branch_weight(law, seam, coin, raw)
                    if p == 0:
                        continue
                    correction = recovery_resources(law, seam, coin, raw)
                    recovery_sum += p*correction["native_pulses"]/len(SEAMS)
                    pulse_sum += p*(native["native_pulses"]+correction["native_pulses"])/len(SEAMS)
                    maximum_recovery = max(maximum_recovery, correction["native_pulses"])
                    maximum_total = max(maximum_total,
                                        native["native_pulses"]+correction["native_pulses"])
        rows[law] = {
            "uniform_seam_expected_native_plus_recovery_pulses": str(pulse_sum),
            "uniform_seam_expected_recovery_pulses": str(recovery_sum),
            "maximum_one_branch_native_plus_recovery_pulses": maximum_total,
            "maximum_one_branch_recovery_pulses": maximum_recovery,
            "retained_coin_and_raw_bits_per_attempt": 1+circuit.RECORD_BITS[law],
            "retained_uniform_seam_selector_bits_per_attempt": 5,
            "retained_core_reset_archive_bits_per_attempt": 15,
            "fixed_format_total_classical_bits_per_uniform_attempt":
                1+circuit.RECORD_BITS[law]+5+15,
            "active_core_reset_probability_uniform_seam": "1/15",
        }
    return {
        "encoder": encoding_resources(),
        "decoder": {**encoding_resources(), "output_rows_retained": SIZE,
                    "syndrome_postselection_applied": False},
        "recovery_zero_weight_branch_resources": zero,
        "source_pulse_reuse_of_numerical_cache_is_not_free": True,
        "initial_core_label_bits": 7,
        "classical_core_reset_is_a_supplied_service_without_source_pulse_cost": True,
        "active_recovery_has_no_new_quantum_scratch_register": True,
        "laws": rows,
    }


def source_composition_report():
    """Numerical full-register product bound for every emitted source word.

    This is a floating diagnostic.  It is neither a directed-rounding bound
    nor a hardware error rate.  It does not include classical reset faults.
    """
    exact_flip = np.eye(SIZE, dtype=complex)
    exact_flip[[0, 256]] = exact_flip[[256, 0]]
    macro_errors = {}
    for inverse in (False, True):
        actual, guard = routing._executed_macro(inverse)
        macro_errors[inverse] = (float(np.linalg.norm(actual-exact_flip, 2))
                                 +guard["projection_change"])
    exact_ix = np.asarray([[0, 1j], [1j, 0]])
    xerror = max(float(np.linalg.norm(
        compiler.native_pulse("X", sign*math.pi/2)[:2, :2]-sign*exact_ix, 2))
                 for sign in (1, -1))
    mixing = {}
    for count in (1, 2, 3, 4):
        literal = tuple(tuple((event["kind"], event["angle"], tuple(event["targets"]))
                              for event in compiler._validate(_mixing_word(count, inverse), 9))
                        for inverse in (False, True))
        if literal[1] != tuple((kind, -angle, targets)
                               for kind, angle, targets in reversed(literal[0])):
            raise ArithmeticError("mixing inverse is not the emitted reversed word")
        synthesis = record_code._mixing_synthesis_report(count, *literal)
        mixing[count] = synthesis
    mixing_logs = {(count, inverse): (
        record_code._word_log_error(_mixing_word(count, inverse))
        +math.log1p(mixing[count]["uniform_target_unitary_completion_error_upper"]))
        for count in (1, 2, 3, 4) for inverse in (False, True)}

    def permutation_log(images, inverse=False):
        edges = routing.permutation_edges(images, inverse=inverse)
        xcount = 2*sum(pattern.bit_count() for _, pattern in edges)
        return len(edges)*math.log1p(macro_errors[inverse])+xcount*math.log1p(xerror)

    encoder_log = permutation_log(_encoder_permutation())+mixing_logs[4, False]
    decoder_log = permutation_log(_encoder_permutation(), True)+mixing_logs[4, True]
    coin = np.asarray([[1., -1.], [1., 1.]])/math.sqrt(2)
    coin_error = float(np.linalg.norm(compiler.native_pulse("Y", -math.pi/4)[:2, :2]-coin, 2))
    rows = {}
    for law in ("C", "E"):
        native_error = max(record_code._native_dilation_error(law, coin_bit)
                           for coin_bit in (0, 1))
        worst_log = 0.
        for seam in SEAMS:
            route = record_code._route_permutation(seam)
            prefix = (permutation_log(route)+permutation_log(route, True)
                      +math.log1p(coin_error)+math.log1p(native_error))
            for coin_bit in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    count = len(branch_images(law, seam, coin_bit, raw)[0])
                    recovery_log = 0.
                    if count:
                        recovery_log = (
                            permutation_log(_branch_permutation(law, seam, coin_bit, raw), True)
                            +mixing_logs[count, True]+encoder_log)
                    worst_log = max(worst_log, prefix+recovery_log)
        rows[law] = {
            "maximum_one_corrected_attempt_log_error": worst_log,
            "maximum_one_corrected_attempt_operator_error": math.expm1(worst_log),
            "encoded_t_attempt_operator_error_upper": {
                str(t): math.expm1(encoder_log+t*worst_log) for t in (1, 2, 10, 100)},
            "all_thirty_seams_and_raw_control_blocks_counted": True,
            "complete_local_native_dilation_error": native_error,
        }
    if any(row["encoded_t_attempt_operator_error_upper"]["100"] > 1e-4
           for row in rows.values()):
        raise ArithmeticError("numerical source composition diagnostic exceeds its gate")
    return {
        "method": "product(1+factor_operator_error)-1; max over retained classical control blocks",
        "forward_full512_macro_operator_error_plus_polar_change": macro_errors[False],
        "inverse_full512_macro_operator_error_plus_polar_change": macro_errors[True],
        "source_X_conjugation_operator_error": xerror,
        "fair_coin_preparation_operator_error": coin_error,
        "mixing_angle_synthesis_diagnostics": {str(k): v for k, v in mixing.items()},
        "encoder_operator_error_bound": math.expm1(encoder_log),
        "decoder_operator_error_bound": math.expm1(decoder_log),
        "numerical_source_composition_not_interval_hardware_or_reset_certificate": True,
        "arbitrary_reference_extension_has_same_operator_bound": True,
        "laws": rows,
    }


def _selected_source_step(law, seam, coin, raw, state):
    """Execute full local raw instrument; then select and correct one flag."""
    _label(law, seam, coin, raw)
    route = record_code._route_permutation(seam)
    local = record_code._apply_permutation(route, routing._state(state))
    branches = recorded.repair_branches(law, local)
    selected = next(b for b in branches if (b["coin"], b["record_code"]) == (coin, raw))
    if selected["abort"] or branch_weight(law, seam, coin, raw) == 0:
        raise ValueError("source replay needs an admitted positive record")
    native_output = record_code._apply_permutation(route, selected["amplitude"], True)
    corrected = apply_recovery(law, seam, coin, raw, native_output)
    return corrected, branches, native_output


def source_replay_report():
    """Full512 numerical replay of all correctors and a C/E continuation.

    Every positive branch's emitted inverse permutation, inverse mixing word
    and encoder acts on its two native output columns. The two-step sample
    additionally executes the native C and E source instruments with routing.
    """
    code = code_matrix()
    target = ideal_code_matrix()
    encoder_error = float(np.linalg.norm(code-target, 2))
    samples = []
    # The E branch moves the co-moving core from (0,1) to (1,5); the source
    # recovery returns the fixed code and the classical core is then reset.
    schedule = (("C", (0, 1), 0, 0), ("E", (0, 5), 0, 5))
    state, probability = code, F(1)
    for law, seam, coin, raw in schedule:
        state, all_raw, native_output = _selected_source_step(law, seam, coin, raw, state)
        p = branch_weight(law, seam, coin, raw)
        probability *= p
        expected = math.sqrt(float(probability))*target
        raw_norm = sum(float(np.linalg.norm(b["amplitude"])**2)/2 for b in all_raw)
        samples.append({
            "law": law, "seam": list(seam), "coin": coin, "raw": raw,
            "branch_weight_exact": str(p), "cumulative_branch_weight_exact": str(probability),
            "retained_raw_slots": len(all_raw),
            "raw_instrument_total_probability_for_two_input_columns": raw_norm,
            "post_native_column_gram_error": float(np.linalg.norm(
                native_output.conj().T@native_output-float(probability)*np.eye(2), 2)),
            "corrected_against_exact_code_operator_error": float(np.linalg.norm(state-expected, 2)),
            "core_archive": archive_entry(law, seam, coin, raw),
        })
    all_branch = {}
    for law in ("C", "E"):
        count_hist, maximum, checked = Counter(), 0., 0
        for seam in SEAMS:
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    rows = branch_images(law, seam, coin, raw)
                    count = len(rows[0])
                    if not count:
                        continue
                    image = np.zeros((SIZE, DIMENSION), complex)
                    for logical, row in enumerate(rows):
                        image[list(row), logical] = 1/math.sqrt(8)
                    output = apply_recovery(law, seam, coin, raw, image)
                    error = float(np.linalg.norm(
                        output-math.sqrt(float(branch_weight(law, seam, coin, raw)))*target, 2))
                    maximum = max(maximum, error)
                    count_hist[count] += 1
                    checked += 1
        all_branch[law] = {
            "positive_source_recoveries_executed": checked,
            "support_count_histogram": {str(k): v for k, v in sorted(count_hist.items())},
            "maximum_full512_source_recovery_error": maximum,
        }
    if sum(row["positive_source_recoveries_executed"] for row in all_branch.values()) != 208:
        raise ArithmeticError("a positive source branch was omitted")
    errors = [encoder_error]
    errors.extend(row["corrected_against_exact_code_operator_error"] for row in samples)
    errors.extend(row["maximum_full512_source_recovery_error"] for row in all_branch.values())
    if max(errors) > 1e-8:
        raise ArithmeticError("screened source replay exceeds its numerical gate")
    return {
        "encoder_operator_error": encoder_error,
        "corrected_two_attempt_CE_continuation": samples,
        "all_positive_raw_branch_source_recoveries": all_branch,
        "maximum_source_replay_error": max(errors),
        "maximum_allowed_numerical_execution_error": 1e-8,
        "all_512_output_rows_retained": True,
        "actual_native_raw_instrument_executed_on_continuation": True,
        "all_thirty_seam_recovery_source_replay": True,
        "all_thirty_seam_native_instrument_source_replay": False,
        "physical_service_or_roundoff_certificate": False,
    }

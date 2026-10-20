"""Continue the source-prepared edge reference with a clock and feedback.

The complete evidence keeps the sampler request, every first herald failure,
both native seam/coin/current records, all sixteen detector outcomes and a
twelve-dimensional untouched reference. Modular identification concerns the
successful work marginal, not the larger classical-record/reference algebra.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer import native, quantum_mixing
from ..observer.source_scalar import compiler
from . import inverse_clock as preparation
from . import prepared_clock_source as source

ROOT = Path(__file__).resolve().parents[2]
INITIAL = (0,)*10+(1, 1)
CLOCK_PORT = 10
BAND = 2  # One vector of the lambda=6 sector, not its multiplicity.
FILTER_GUARD = F(1, 10**9)
ARITHMETIC_ALLOWANCE = F(1, 10**10)


def _array(value, shape):
    if (type(value) is not np.ndarray or value.shape != shape
            or value.dtype != np.dtype("complex128") or not np.isfinite(value).all()):
        raise ValueError(f"original finite complex128 array of shape {shape} required")
    return np.frombuffer(value.tobytes(), dtype=np.complex128).reshape(shape)


def _choice(axis, intervention):
    if type(axis) is not str or axis not in ("X", "Y"):
        raise ValueError("the complete X or Y detector setting is required")
    if type(intervention) is not int or intervention not in (0, 1):
        raise ValueError("an original binary intervention setting is required")


def _native_record(state, edge, coin):
    seam = quantum_mixing.seams()[edge]
    after, (_, current) = native.native(state, seam, coin)
    return after, (edge, coin, current, int(CLOCK_PORT in seam))


@lru_cache(maxsize=1)
def native_histories():
    """All 3600 native tapes, grouped only by identical work operations.

    Sampling the request first uses 64 fair placement coins on the empty
    seam (0,1). It leaves INITIAL fixed for every word, so it introduces no
    unrecorded dependence between that request and the later native carrier.
    """
    if quantum_mixing.seams()[0] != (0, 1):
        raise ArithmeticError("the fixed empty sampler seam changed")
    for coin in (0, 1):
        after, record = native.native(INITIAL, (0, 1), coin)
        if after != INITIAL or record != (coin, 0):
            raise ArithmeticError("sampling coins disturb the subsequent native carrier")
    groups = defaultdict(list)
    for e1 in range(30):
        for c1 in (0, 1):
            middle, first = _native_record(INITIAL, e1, c1)
            # The first committed local current is the entire feedback input.
            left, right = quantum_mixing.seams()[e1]
            local_current = first[2]*(int(right == CLOCK_PORT)-int(left == CLOCK_PORT))
            feedback = int(local_current < 0)
            for e2 in range(30):
                for c2 in (0, 1):
                    final, second = _native_record(middle, e2, c2)
                    control = first[3], feedback, second[3]
                    groups[control].append((first, second, final))
    return tuple((key, tuple(rows)) for key, rows in sorted(groups.items()))


def _ideal_filter():
    """Independent endpoint dilation in the declared event chart.

    Only the event basis completion is reused from the existing preparation.
    The shift is evaluated directly from seam endpoints, without its circuit,
    permutation executor, filter compression or native transfer matrix.
    """
    _, event_basis = preparation._lazy_filter_program(F(1, 2))
    output = np.zeros((1 << 15, 12), dtype=complex)
    for event in range(60):
        seam = quantum_mixing.seams()[event//2]
        for site in range(12):
            destination = seam[0]+seam[1]-site if event & 1 and site in seam else site
            for final_event in range(64):
                output[destination+16*final_event, site] += (
                    event_basis[event, final_event]*event_basis[event, 0])
    return output


@dataclass(frozen=True)
class Preparation:
    weights: tuple[F, ...]
    rejection: F
    successes: tuple[np.ndarray, ...]
    failures: tuple[np.ndarray, ...]


def prepare(*, ideal=False):
    """Unnormalized amplitudes on all work/environment/reference outputs."""
    if type(ideal) is not bool:
        raise ValueError("an explicit ideal comparison flag is required")
    _, data = preparation._prepared_edge_gibbs_data()
    raw = _ideal_filter() if ideal else data["raw_filter"]
    raw = _array(raw, (1 << 15, 12))
    weights = preparation.subordinated_filter_weights()
    successful = np.eye(12, dtype=complex)/np.sqrt(12.)
    successes, failures = [], []
    for n in range(len(weights)):
        embedded = np.zeros((16, 12), dtype=complex)
        embedded[:12] = successful
        successes.append(_array(embedded, (16, 12)))
        if n+1 < len(weights):
            full = raw@successful
            failures.append(_array(full[12:].copy(), ((1 << 15)-12, 12)))
            successful = full[:12]
    return Preparation(weights, 1-sum(weights), tuple(successes), tuple(failures))


def _maps(ideal):
    compiled = source.compiled_source()
    return compiled["targets" if ideal else "raw"]


def _prefix(amplitude, control, intervention, maps):
    tick, feedback = control[:2]
    result = maps["clock1"]@amplitude if tick else amplitude
    if feedback and intervention:
        result = maps[f"kick_plus_{BAND}"]@result
    return result


def live_prefix(*, intervention=1, ideal=False):
    """Prepare, execute the first native tick and its committed feedback."""
    _choice("X", intervention)
    prepared, maps = prepare(ideal=ideal), _maps(ideal)
    controls = sorted({key[:2] for key, _ in native_histories()})
    return {(n, *control): _array(_prefix(amplitude, control, intervention, maps), (16, 12))
            for n, amplitude in enumerate(prepared.successes) for control in controls}


def checkpoint(*, intervention=1):
    """Serialize a replay-bound live prefix without normalizing its arrays.

    Failed/rejected outputs are the immutable completed preparation log; the
    checkpoint names that log by its raw binary digest. Replay checks its
    custody before continuing the serialized live arrays.
    """
    _choice("X", intervention)
    prefix, prepared = live_prefix(intervention=intervention), prepare()
    rows = []
    for key, value in sorted(prefix.items()):
        rows.append([list(key), value.astype("<c16", copy=False).tobytes().hex()])
    payload = {
        "schema": "oph.prepared_clock.prefix.v1", "intervention": intervention,
        "source_sha256": bindings(),
        "completed_preparation_sha256": _completed_digest(prepared),
        "live_arrays_complex128_little_endian": rows,
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _completed_digest(prepared):
    digest = hashlib.sha256()
    digest.update(str((prepared.weights, prepared.rejection)).encode())
    # Rejection leaves the supplied joint work/reference input untouched.
    digest.update(prepared.successes[0].tobytes())
    for value in prepared.failures:
        digest.update(value.tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class History:
    axis: str
    intervention: int
    preparation: Preparation
    # Each row retains every native tape with the same CP work operation.
    # The terminal address read projects each of the sixteen rows separately.
    successful: dict
    resumed: bool = False


def execute(*, axis="Y", intervention=1, ideal=False, saved=None):
    _choice(axis, intervention)
    if type(ideal) is not bool:
        raise ValueError("an explicit ideal comparison flag is required")
    prepared, maps = prepare(ideal=ideal), _maps(ideal)
    prefix = live_prefix(intervention=intervention, ideal=ideal)
    if saved is not None:
        if ideal or type(saved) is not str or saved != checkpoint(intervention=intervention):
            raise ValueError("checkpoint differs from the complete original source-prefix replay")
        payload = json.loads(saved)
        prefix = {tuple(key): _array(np.frombuffer(bytes.fromhex(raw), dtype="<c16")
                                    .reshape(16, 12).copy(), (16, 12))
                  for key, raw in payload["live_arrays_complex128_little_endian"]}
    result = {}
    for n in range(len(prepared.weights)):
        for control, _ in native_histories():
            amplitude = prefix[n, *control[:2]]
            if control[2]:
                amplitude = maps["clock1"]@amplitude
            amplitude = maps[f"detector_{axis.lower()}_{BAND}"]@amplitude
            result[n, *control] = _array(amplitude, (16, 12))
    return History(axis, intervention, prepared, result, saved is not None)


def _pure_distance(a, b):
    """Half trace norm of |a><a|-|b><b|, using a stable upper bound.

    Vector phases are irrelevant to a classical quantum-output block.
    The independent algebraic rank-two norm formula is checked in tests.
    """
    if np.array_equal(a, b):
        return 0.
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    overlap = np.vdot(b, a)
    phase = overlap/abs(overlap) if overlap else 1.
    return float((na+nb)*np.linalg.norm(a-phase*b)/2)


def _validate(actual):
    if type(actual) is not History:
        raise ValueError("a complete typed clock history is required")
    if type(actual.resumed) is not bool:
        raise ValueError("the continuation marker must be an original Boolean")
    _choice(actual.axis, actual.intervention)
    p = actual.preparation
    weights = preparation.subordinated_filter_weights()
    if (type(p) is not Preparation or type(p.weights) is not tuple
            or any(type(w) is not F or type(w.numerator) is not int
                   or type(w.denominator) is not int or w.denominator <= 0
                   or math.gcd(w.numerator, w.denominator) != 1 for w in p.weights)
            or p.weights != weights
            or type(p.rejection) is not F or p.rejection != 1-sum(weights)
            or type(p.successes) is not tuple or len(p.successes) != len(weights)
            or type(p.failures) is not tuple or len(p.failures) != len(weights)-1):
        raise ValueError("the complete unnormalized preparation menu is required")
    for state in p.successes:
        _array(state, (16, 12))
    for state in p.failures:
        _array(state, ((1 << 15)-12, 12))
    expected = {(n, *control) for n in range(len(weights)) for control, _ in native_histories()}
    if (type(actual.successful) is not dict
            or any(type(key) is not tuple or any(type(x) is not int for x in key)
                   for key in actual.successful)
            or set(actual.successful) != expected):
        raise ValueError("every request and native-control outcome must be retained")
    for state in actual.successful.values():
        _array(state, (16, 12))


def complete_distance(actual, expected):
    """Compare every failure and every terminal quantum-reference block."""
    _validate(actual); _validate(expected)
    if (actual.axis, actual.intervention) != (expected.axis, expected.intervention):
        raise ValueError("identical observation menus are required for comparison")
    p, q = actual.preparation, expected.preparation
    result = float(p.rejection)*_pure_distance(p.successes[0], q.successes[0])
    for n, weight in enumerate(p.weights):
        for k in range(n):
            result += float(weight)*_pure_distance(p.failures[k], q.failures[k])
        for control, tapes in native_histories():
            a, b = actual.successful[n, *control], expected.successful[n, *control]
            result += float(weight*F(len(tapes), 3600))*sum(
                _pure_distance(a[j], b[j]) for j in range(16))
    return result


def probabilities(actual):
    """Complete observer law with native tapes factored by identical outputs."""
    _validate(actual)
    result = {("sampler_reject",): float(actual.preparation.rejection)*float(
        np.linalg.norm(actual.preparation.successes[0])**2)}
    for n, weight in enumerate(actual.preparation.weights):
        for k in range(n):
            result["herald_failure", n, k+1] = float(weight)*float(
                np.linalg.norm(actual.preparation.failures[k])**2)
        for control, tapes in native_histories():
            output = actual.successful[n, *control]
            for address in range(16):
                result["read", n, *control, address] = float(weight*F(len(tapes), 3600))*float(
                    np.vdot(output[address], output[address]).real)
    return result


def bindings():
    paths = {*preparation.SOURCE_PATHS,
             "model/predictions/prepared_clock.py", "model/predictions/prepared_clock_source.py",
             "tools/strict_json.py"}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def error_budget():
    """Complete CP composition bound under the checked numerical guards.

    These are finite numerical checks, not directed-rounding certificates.
    All first failures are outcomes of each filter. Thus n filters have gain
    at most g**n and half-diamond error (g**n-1)/2, g=(1+delta)**2.
    A triangle through ideal preparation charges continuation error only on
    its successful mass; no successful branch is renormalized.
    """
    _, data = preparation._prepared_edge_gibbs_data()
    residual = float(np.linalg.norm(data["raw_filter"]-_ideal_filter(), "fro"))
    if not math.isfinite(residual) or residual > float(FILTER_GUARD):
        raise ArithmeticError("the complete raw filter exceeds its numerical guard")
    compiled = source.compiled_source()
    guards = compiled["guards"]
    # Both axes and every branch use at most two clocks, one kick and one read.
    delta = ((1+guards["clock1"])**2*(1+guards[f"kick_plus_{BAND}"])
             *(1+max(guards[f"detector_{a}_{BAND}"] for a in ("x", "y")))-1)
    gain = (1+delta)**2
    continuation_error = (gain-1)/2
    weights = preparation.subordinated_filter_weights()
    g = (1+FILTER_GUARD)**2
    prep_error = sum((w*(g**n-1)/2 for n, w in enumerate(weights)), F())
    success = F(source.prepared_reference()["normalization"])/12
    combined = gain*prep_error+success*continuation_error+ARITHMETIC_ALLOWANCE
    # Outward rational reporting keeps the receipt reviewable without changing
    # the exact product calculation that establishes each upper allowance.
    upper = lambda value: F(math.ceil(value*10**12), 10**12)
    return {"filter_operator_guard": FILTER_GUARD,
        "filter_complete_Frobenius_residual": residual,
        "stopped_preparation_half_diamond_upper_under_guards": upper(prep_error),
        "continuation_half_diamond_upper_under_guards": upper(continuation_error),
        "continuation_channel_norm_upper_under_guards": upper(gain),
        "complete_half_trace_upper_under_guards": upper(combined),
        "complete_event_weight_error_upper_under_guards": upper(2*combined),
        "numerical_arithmetic_allowance": ARITHMETIC_ALLOWANCE,
        "ideal_success_probability": success,
        "conditioned_success_half_trace_upper_under_guards": min(F(1), upper(2*combined/success)),
        "directed_rounding_or_physical_fidelity_certificate": False}


def audit(actual, target):
    """The acceptance gate checks complete quantum outputs, not just masses."""
    error = complete_distance(actual, target)
    budget = error_budget()
    allowance = float(budget["complete_half_trace_upper_under_guards"])
    probabilities_actual = probabilities(actual)
    mass = sum(probabilities_actual.values())
    if (not math.isfinite(error) or not math.isfinite(mass)
            or error > allowance or abs(mass-1) > 2*allowance
            or any(not math.isfinite(p) or p < 0 for p in probabilities_actual.values())):
        raise ArithmeticError("complete source clock/history quantum outputs exceed the budget")
    # Returned preparation states are part of the exposed continuation data.
    for a, b in zip(actual.preparation.successes, target.preparation.successes, strict=True):
        if _pure_distance(a, b) > 4*float(FILTER_GUARD)+float(ARITHMETIC_ALLOWANCE):
            raise ArithmeticError("a source preparation checkpoint exceeds its allowance")
    return {"axis": actual.axis, "intervention": actual.intervention,
        "complete_quantum_half_trace_upper_diagnostic": error,
        "complete_raw_mass_diagnostic": mass,
        "retained_coarse_records": len(probabilities_actual),
        "native_tapes_per_successful_request": 3600,
        "complete_uncompressed_records": 1+10+len(actual.preparation.weights)*3600*16,
        "all_terminal_reference_entries_compared": True,
        "all_first_failure_environment_rows_compared": True,
        "resumed_from_serialized_live_prefix": actual.resumed}


def resources():
    edge, _ = preparation._prepared_edge_gibbs_data()
    compiled = source.compiled_source()
    words = compiled["words"]
    continuation_word = (words["clock1"]+words[f"kick_plus_{BAND}"]
                         +words["clock1"]+words[f"detector_y_{BAND}"])
    continuation_pulses = max(
        2*len(words["clock1"])+len(words[f"kick_plus_{BAND}"])+len(words[f"detector_{a}_{BAND}"])
        for a in ("x", "y"))
    weights = preparation.subordinated_filter_weights()
    reference = source.prepared_reference()
    lazy = reference["lazy_filter"]
    expected_filters = sum((w*sum((F((lazy**(2*k)).trace())/12 for k in range(n)), F())
                            for n, w in enumerate(weights)), F())
    success = F(reference["normalization"])/12
    expected_continuation = success*(
        F(1, 3)*len(words["clock1"])+F(1, 15)*len(words[f"kick_plus_{BAND}"])
        +max(len(words[f"detector_{a}_{BAND}"]) for a in ("x", "y")))
    filter_pulses = edge["one_filter_resources"]["native_pulses"]
    maximum = 4*filter_pulses+continuation_pulses
    info = compiler.resources(continuation_word)
    filter_program, _ = preparation._lazy_filter_program(F(1, 2))
    program_bytes = len(json.dumps(
        {"filter": filter_program, "continuation_words": words},
        sort_keys=True, separators=(",", ":"), allow_nan=False).encode())
    return {"maximum_filter_requests": 4,
        "one_filter_source_pulses": filter_pulses,
        "maximum_successful_continuation_source_pulses": continuation_pulses,
        "maximum_total_source_pulses": maximum,
        "maximum_proper_code_load_events": maximum,
        "maximum_proper_code_export_events": maximum,
        "mean_executed_filters_with_first_failure_stopping": expected_filters,
        "ideal_mean_total_source_pulses_upper": expected_filters*filter_pulses+expected_continuation,
        "fair_native_sampler_coin_events": 64,
        "subsequent_native_attempts_maximum": 2,
        "ideal_subsequent_native_attempts_mean": 2*success,
        "clock_ticks_per_successful_trial_mean": F(1, 3),
        "feedback_kicks_per_successful_trial_mean": F(1, 15),
        "native_carriers": 1, "native_ports": 12,
        "protected_work_qubits": 4, "reference_storage_qubits": 4,
        "reference_populated_dimension": 12,
        "event_and_scratch_qubits_during_filter": 11,
        "peak_encoded_qubits_including_reference": 19,
        "source_processor_dimension": 6,
        "classical_record_storage_bits": 64+2*8+4+3+3+1,
        "record_encoding": "64 sampler bits; two seam(5)/coin(1)/current(2) records; "
                           "address(4), request(3), failure index(3), status(1); schema distinguishes branches",
        "fixed_attempt_horizon_not_wait_for_two_ticks": True,
        "zero_tick_trials_retained": True,
        "mode_basis_and_detector_costs_included": True,
        "declared_source_program_utf8_bytes": program_bytes,
        "program_cursor_bits_upper": maximum.bit_length(),
        "filter_request_counter_bits": 3,
        "attempt_counter_bits": 2,
        "herald_measurements_maximum": 4,
        "herald_measurements_ideal_mean": expected_filters,
        "terminal_address_measurements_maximum": 1,
        "terminal_address_measurements_ideal_mean": success,
        "initial_native_state_work_I12_and_purification_supplied": True,
        "herald_and_terminal_readout_services_supplied": True,
        "source_filter_routine_reuse_removes_no_service_invocations": True,
        "example_continuation_load_export_events": info["proper_code_service_events"],
        "physical_duration_or_error_of_supplied_services_bounded": False}


def sampler_contract():
    weights = preparation.subordinated_filter_weights()
    rows, boundary = [], 0
    for label, weight in enumerate((*weights, 1-sum(weights))):
        width = weight*(1 << 64)
        if width.denominator != 1:
            raise ArithmeticError("the request sampler is not on its declared coin grid")
        width = int(width)
        if not width:
            continue
        selected = label if label < len(weights) else None
        for integer in (boundary, boundary+width-1):
            state, records = INITIAL, []
            for bit in map(int, f"{integer:064b}"):
                state, record = native.native(state, (0, 1), bit)
                records.append(record)
            decoded = preparation.select_subordinated_request(
                tuple(coin for coin, _ in records), weights)
            if decoded != selected or state != INITIAL:
                raise ArithmeticError("request sampler does not continue the same native carrier")
        rows.append({"request": selected, "first_word": boundary,
                     "past_last_word": boundary+width, "probability": weight})
        boundary += width
    if boundary != 1 << 64:
        raise ArithmeticError("request and rejection intervals omit native coin words")
    return rows


def run():
    budget = error_budget()
    compiled = source.compiled_source()
    rows, executions = [], {}
    for axis in ("X", "Y"):
        for intervention in (0, 1):
            saved = checkpoint(intervention=1) if (axis, intervention) == ("Y", 1) else None
            actual = execute(axis=axis, intervention=intervention, saved=saved)
            target = execute(axis=axis, intervention=intervention, ideal=True)
            rows.append(audit(actual, target))
            executions[axis, intervention] = actual
    # A source-specific nonzero phase signal on one *unconditional* native
    # record event: first tick loses the port load, second attempt also ticks.
    def contrast(history):
        law = probabilities(history)
        address = compiled["pair_addresses"][BAND]
        return sum(law["read", n, 1, 1, 1, 0]-law["read", n, 1, 1, 1, address]
                   for n in range(len(history.preparation.weights)))
    reference = compiled["reference"]
    mode = compiled["pair_addresses"][BAND]
    m0, mj = float(reference["m0"]), float(reference["masses"][mode])
    gap = reference["modular_gaps"][mode]
    analytic = (m0-mj)*math.sin(gap)/(12*90)
    signal = contrast(executions["X", 1])-contrast(executions["X", 0])
    allowance = 4*float(budget["complete_half_trace_upper_under_guards"])
    if analytic <= allowance or abs(signal-analytic) > allowance:
        raise ArithmeticError("the retained source clock/feedback signal is not resolved")
    # Execute the clock-off control through the same source preparation and
    # noncommuting feedback; do not infer its absence from a clock marginal.
    null_histories = []
    for ideal in (False, True):
        prepared, maps = prepare(ideal=ideal), _maps(ideal)
        outputs = {}
        for n, state in enumerate(prepared.successes):
            for control, _ in native_histories():
                kicked = maps[f"kick_plus_{BAND}"]@state if control[1] else state
                outputs[n, *control] = _array(maps[f"detector_x_{BAND}"]@kicked, (16, 12))
        null_histories.append(History("X", 1, prepared, outputs))
    audit(null_histories[0], null_histories[1])
    null_signal = contrast(null_histories[0])-contrast(executions["X", 0])
    if abs(null_signal) > allowance:
        raise ArithmeticError("source clock-off control falsely retains a phase signal")
    # Common H is not modular time of the retained request/work algebra:
    # request zero has a tracial success state, while H is nonconstant.
    request_zero_gap = 0.
    if gap <= .01 or request_zero_gap == gap:
        raise ArithmeticError("the retained-request modular countercontrol disappeared")
    groups = native_histories()
    native_json = json.dumps(groups, separators=(",", ":"))
    return preparation._json({
        "schema": "oph.prepared_clock.v1", "source_sha256": bindings(),
        "complete_sampler_intervals": sampler_contract(),
        "native_record_groups": [{"control": list(key), "tape_count": len(tapes)}
                                  for key, tapes in groups],
        "complete_native_tape_sha256": hashlib.sha256(native_json.encode()).hexdigest(),
        "executions": rows, "error_budget": budget, "resources": resources(),
        "source_words": {"sha256": compiled["word_sha256"],
                         "pulse_counts": compiled["pulse_counts"],
                         "full_register_guards": compiled["guards"]},
        "clock_field_response": {
            "event": "first incident attempt loses port10 load, second attempt is incident",
            "unconditional_native_event_probability": F(1, 90),
            "observable": "X detector address0 minus selected lambda6 address, feedback on minus off",
            "analytic_signal": analytic, "executed_signal": signal,
            "resolved_signal_lower_under_guards": analytic-allowance,
            "executed_clock_off_signal": null_signal,
            "zero_signal_when_clock_or_noncommuting_feedback_is_absent": True,
            "formula": "(m0-m6)*sin(log(m0/m6))/(12*90)"},
        "state_and_clock": {
            "prepared_state": "M/Tr(M), M=sum_n w_n(I-L/120)^(2n)",
            "same_register_clock": "exp(-iH), H=-log(M/m0)",
            "modular_identification": "sigma_s=alpha_(-s) on successful work B(C^12)",
            "reference_and_sampler_request_are_traced_for_this_algebra": True,
            "sampler_request_retained_in_complete_evidence": True,
            "request_zero_modular_gap": request_zero_gap,
            "work_marginal_lambda6_modular_gap": gap,
            "one_clock_is_modular_for_every_retained_request_block": False,
            "faithful_modular_state_on_padding_or_reference_algebra": False,
            "sqrt_L_comparison_clock_compiled": True,
            "finite_prepared_reference_is_exact_sqrt_L_Gibbs_state": False},
        "scope": {
            "one_source_prepared_register_through_clock_feedback_and_readout": True,
            "complete_failed_preparation_and_native_histories_retained": True,
            "same_native_carrier_continues_after_sampler": True,
            "raw_arrays_normalized_or_projected_to_ideal_unitaries": False,
            "source_preparation_filter_and_clock_are_supplied_observer_services": True,
            "growing_product_field_encoder_or_original_scalar_benchmark_completed": False,
            "physical_time_or_empirical_field_dictionary_identified": False}})


def check_receipt(recorded, actual):
    """Strict replay additionally enforces the producer's scientific domains."""
    from ..observer.scalar_response.experiment import compare_receipt
    compare_receipt(recorded, actual)
    budget = actual["error_budget"]
    allowance = float(F(budget["complete_half_trace_upper_under_guards"]))
    filter_guard = float(F(budget["filter_operator_guard"]))
    residual = recorded["error_budget"]["filter_complete_Frobenius_residual"]
    if not 0 <= residual <= filter_guard:
        raise ValueError("stored filter norm violates its scientific guard")
    for row in recorded["executions"]:
        if (not 0 <= row["complete_quantum_half_trace_upper_diagnostic"] <= allowance
                or abs(row["complete_raw_mass_diagnostic"]-1) > 2*allowance):
            raise ValueError("stored complete quantum output violates its scientific guard")
    signal = recorded["clock_field_response"]
    if (signal["resolved_signal_lower_under_guards"] <= 0
            or abs(signal["executed_signal"]-signal["analytic_signal"]) > 4*allowance
            or abs(signal["executed_clock_off_signal"]) > 4*allowance):
        raise ValueError("stored clock signal violates its scientific guard")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path); group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        check_receipt(strict_json.load(args.check), result)
        print("Prepared reference, native clock, feedback and complete outputs replayed.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

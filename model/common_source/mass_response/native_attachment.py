"""Diagonal native repair drives and receives the complete quantum history.

The first committed native current selects whether the later site kick is
enabled. The receiver's committed occupation bit selects the second seam.
Integer loads are carried forward unchanged between those two native events.
Quantum outputs are represented by their existing branch label and an exact
classical weight; no quantum state is prepared, normalized or discarded here.

This attachment covers a supplied uniform DIAGONAL binary native prior.
It does not implement or discriminate the coherent C/E instruments. Native
integer loads and encoded fermion occupation are distinct registers.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math

import numpy as np

from model.observer import native
from . import experiment


RECORDS = tuple(product((0, 1), repeat=5))
RECORD_ORDER = ["weak_charge", "mass_return", "site1_occupied", "clock_X", "clock_Y"]
TOLERANCE = 1e-9


def case_from_current(current):
    """Only the committed current, never hidden loads, reaches this selector."""
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError("a committed native current -1,0,1 required")
    return "kick" if abs(current) > 0 else "no_kick"


def seam_from_receiver(receiver_bit, feedback=True):
    """Receiver outcome zero denotes occupied in the actual source read."""
    if type(receiver_bit) is not int or receiver_bit not in (0, 1) or type(feedback) is not bool:
        raise ValueError("a committed receiver bit and Boolean feedback required")
    return native.F if feedback and receiver_bit == 0 else native.E


def _norm(state):
    if np.ma.isMaskedArray(state):
        raise ValueError("masked quantum branches are inadmissible")
    raw = np.asarray(state)
    if (raw.ndim != 2 or not all(raw.shape) or raw.dtype.kind not in "iufc"
            or not np.isfinite(raw).all()):
        raise ValueError("finite work/reference columns required in every branch")
    with np.errstate(over="ignore", invalid="ignore"):
        probability = float(np.vdot(raw, raw).real)
    if not math.isfinite(probability):
        raise ValueError("finite branch probability required")
    return probability


def _validated(kick, no_kick):
    histories = {"kick": kick, "no_kick": no_kick}
    signature = None
    shape = None
    probabilities = {}
    for name, result in histories.items():
        if type(result) is not experiment.Result or type(result.branches) is not dict or type(result.aborts) is not dict:
            raise ValueError("two actual complete quantum Result objects required")
        if set(result.branches) != set(RECORDS) or any(
                type(record) is not tuple or any(type(bit) is not int for bit in record)
                for record in result.branches):
            raise ValueError("all32 ordered binary quantum records must be retained")
        if result.specification.get("kick_enabled") is not (name == "kick"):
            raise ValueError("quantum histories have mismatched kick settings")
        if result.specification.get("record_order") != RECORD_ORDER:
            raise ValueError("the receiver and clock record meanings must match the experiment")
        spec = {key: value for key, value in result.specification.items() if key != "kick_enabled"}
        encoded = json.dumps(spec, sort_keys=True, allow_nan=False)
        if signature is not None and encoded != signature:
            raise ValueError("the histories must share one preparation, law and clock protocol")
        signature = encoded
        table = {}
        for aborted, branches in ((False, result.branches), (True, result.aborts)):
            for record, state in branches.items():
                if type(record) is not tuple:
                    raise ValueError("branch records must be tuples")
                probability = _norm(state)
                branch_shape = np.asarray(state).shape
                if shape is not None and branch_shape != shape:
                    raise ValueError("all histories must preserve the same work/reference registers")
                shape = branch_shape
                table[aborted, record] = probability
        if abs(math.fsum(table.values()) - 1) > TOLERANCE:
            raise ValueError("quantum histories must be complete without renormalization")
        probabilities[name] = table
    if set(kick.aborts) != set(no_kick.aborts):
        raise ValueError("the two histories must retain the same scratch-failure menu")
    if set(kick.prefix_probabilities) != set(product((0, 1), repeat=2)):
        raise ValueError("all four pre-kick prefix probabilities required")
    if set(kick.prefix_probabilities) != set(no_kick.prefix_probabilities):
        raise ValueError("pre-kick prefix menus differ")
    for prefix in kick.prefix_probabilities:
        values = [result.prefix_probabilities[prefix] for result in histories.values()]
        if any(type(value) not in (int, float) or not math.isfinite(value) or value < 0 for value in values):
            raise ValueError("finite nonnegative prefix probabilities required")
        if abs(values[0] - values[1]) > TOLERANCE:
            raise ValueError("the matched histories changed a pre-kick prefix")
        for name in histories:
            terminal = math.fsum(probabilities[name][False, record] for record in RECORDS if record[:2] == prefix)
            if abs(terminal - histories[name].prefix_probabilities[prefix]) > TOLERANCE:
                raise ValueError("the continuation lost an old prefix probability")
    return histories, probabilities


@dataclass(frozen=True)
class JointPath:
    initial_index: int
    initial_loads: tuple[int, ...]
    after_first_loads: tuple[int, ...]
    first_record: tuple[int, int]
    quantum_case: str
    quantum_record: tuple
    quantum_scale_squared: F
    quantum_probability: float
    aborted: bool
    second_seam: tuple[int, int] | None
    second_record: tuple[int, int] | None
    final_loads: tuple[int, ...]

    @property
    def probability(self):
        return float(self.quantum_scale_squared) * self.quantum_probability


def _paths(histories, probabilities, feedback):
    configurations = native.states(binary=True)
    prior = F(1, len(configurations))
    for initial_index, loads in enumerate(configurations):
        for coin in (0, 1):
            retained, first_record = native.native(loads, native.E, coin)
            name = case_from_current(first_record[1])
            result = histories[name]
            for record in RECORDS:
                seam = seam_from_receiver(record[2], feedback)
                for final_coin in (0, 1):
                    final, second_record = native.native(retained, seam, final_coin)
                    yield JointPath(initial_index, loads, retained, first_record, name, record,
                                    prior / 4, probabilities[name][False, record], False,
                                    seam, second_record, final)
            for record in sorted(result.aborts, key=lambda key: json.dumps(key)):
                # An abort has no invented final event or unused last coin.
                yield JointPath(initial_index, loads, retained, first_record, name, record,
                                prior / 2, probabilities[name][True, record], True,
                                None, None, retained)


def iter_paths(kick, no_kick, feedback=True):
    """Enumerate every joint path, including all zero quantum outcomes."""
    if type(feedback) is not bool:
        raise ValueError("Boolean native feedback required")
    histories, probabilities = _validated(kick, no_kick)
    yield from _paths(histories, probabilities, feedback)


def quantum_output(path, kick, no_kick):
    """Recover a path's unchanged full quantum output with its exact weight.

    The classical prior and fair coins multiply the squared amplitude by
    1/264 on normal paths and1/132 on aborts. No division by a branch
    probability occurs; untouched reference columns remain intact.
    """
    if type(path) is not JointPath:
        raise ValueError("an enumerated JointPath required")
    result = {"kick": kick, "no_kick": no_kick}[path.quantum_case]
    branch = (result.aborts if path.aborted else result.branches)[path.quantum_record]
    return math.sqrt(float(path.quantum_scale_squared)) * branch


def attach(kick, no_kick, feedback=True):
    """Summarize the literal native/work/native execution and fixed-E control."""
    if type(feedback) is not bool:
        raise ValueError("Boolean native feedback required")
    histories, probabilities = _validated(kick, no_kick)
    configurations = native.states(binary=True)
    index = {loads: i for i, loads in enumerate(configurations)}
    first_law = defaultdict(F)
    movement = F(0)
    for loads in configurations:
        for coin in (0, 1):
            _, record = native.native(loads, native.E, coin)
            first_law[record] += F(1, 132)
            movement += F(abs(record[1]) > 0, 132)
    if movement != F(5, 33):
        raise ArithmeticError("the actual binary native movement prior changed")

    def census(enabled):
        final_distribution = defaultdict(list)
        observed_first = defaultdict(list)
        public_transcript = defaultdict(list)
        normal, abort = [], []
        records = defaultdict(list)
        count_normal = count_abort = 0
        digest = hashlib.sha256()
        for path in _paths(histories, probabilities, enabled):
            if (sum(path.initial_loads) != 2 or sum(path.after_first_loads) != 2
                    or sum(path.final_loads) != 2 or min(path.final_loads) < 0):
                raise ArithmeticError("the actual integer native state lost total load")
            weight = path.probability
            final_distribution[index[path.final_loads]].append(weight)
            observed_first[path.first_record].append(weight)
            if path.aborted:
                count_abort += 1
                abort.append(weight)
                observed = [path.first_record, ["terminal_abort", path.quantum_record]]
            else:
                count_normal += 1
                normal.append(weight)
                records[path.quantum_record].append(weight)
                observed = [path.first_record, path.quantum_record[2], path.second_record]
            # The policy's seam label and the hidden integer loads are not
            # observations in this distinguishability diagnostic.
            public_transcript[json.dumps(observed, separators=(",", ":"))].append(weight)
            digest.update(json.dumps([path.initial_index, path.first_record, path.quantum_case,
                path.quantum_record, str(path.quantum_scale_squared), path.aborted,
                path.second_seam, path.second_record, path.final_loads],
                separators=(",", ":")).encode())
        final = {key: math.fsum(values) for key, values in final_distribution.items()}
        mass = math.fsum([*normal, *abort])
        if abs(mass - 1) > TOLERANCE:
            raise ArithmeticError("the complete joint native/quantum law lost mass")
        return {"normal_paths": count_normal, "abort_paths": count_abort,
                "normal_probability": math.fsum(normal), "abort_probability": math.fsum(abort),
                "total_probability": mass, "enumeration_sha256": digest.hexdigest(),
                "first_record_marginal_residual": max(abs(math.fsum(observed_first[key]) - float(value))
                                                      for key, value in first_law.items()),
                "port0_mean": math.fsum(value * configurations[key][0] for key, value in final.items()),
                "final_native_distribution": {str(key): value for key, value in sorted(final.items())},
                "public_native_receiver_transcript": {key: math.fsum(values)
                    for key, values in sorted(public_transcript.items())},
                "quantum_record_marginal": {record: math.fsum(values) for record, values in records.items()}}

    actual, fixed = census(feedback), census(False)
    mixed = actual.pop("quantum_record_marginal")
    fixed.pop("quantum_record_marginal")
    case_weights = {"kick": movement, "no_kick": 1 - movement}
    mixture_error = max(abs(mixed[record] - math.fsum(float(case_weights[name]) * probabilities[name][False, record]
                                                    for name in histories)) for record in RECORDS)
    first = actual["final_native_distribution"]
    second = fixed["final_native_distribution"]
    total_variation = math.fsum(abs(first.get(key, 0.) - second.get(key, 0.))
                              for key in first.keys() | second.keys()) / 2
    observed, control = (case["public_native_receiver_transcript"] for case in (actual, fixed))
    observed_total_variation = math.fsum(abs(observed.get(key, 0.) - control.get(key, 0.))
                                        for key in observed.keys() | control.keys()) / 2
    selector_invariance = {}
    for position, label in ((1, "mass_return"), (3, "clock_X"), (4, "clock_Y")):
        source = {name: math.fsum(probabilities[name][False, record] for record in RECORDS
                                 if record[position] == 0) for name in histories}
        native_mean = math.fsum(value for record, value in mixed.items() if record[position] == 0)
        selector_invariance[label] = {"kick": source["kick"], "no_kick": source["no_kick"],
                                     "native_mixture": native_mean,
                                     "between_settings_residual": abs(source["kick"] - source["no_kick"])}
    return {
        "schema": "oph.common_source.mass_response.native_attachment.v1",
        "native_diagonal_prior": {"binary_configurations": 66, "probability_each_exact": "1/66",
                                  "initial_load": 2, "placement_coins_exact": "1/2"},
        "kick_probability_exact": str(movement), "no_kick_probability_exact": str(1 - movement),
        "normal_path_factor_exact": "1/264", "abort_path_factor_exact": "1/132",
        "first_controller_input": "abs(committed native current)>0 enables the later kick",
        "final_controller_input": "committed receiver record[2]==0 selects seam F; otherwise E",
        "feedback_enabled": feedback, "actual": actual, "fixed_E_control": fixed,
        "final_native_classical_total_variation": total_variation,
        "public_native_receiver_transcript_total_variation": observed_total_variation,
        "public_transcript_record_order": {
            "normal": ["first native coin/current", "receiver bit", "last native coin/current"],
            "abort": ["first native coin/current", "terminal abort and retained quantum record"],
        },
        "port0_mean_feedback_difference": actual["port0_mean"] - fixed["port0_mean"],
        "quantum_mixture_residual": mixture_error, "pre_kick_mass_and_clock_invariance": selector_invariance,
        "native_events_per_nonabort_path": 2, "native_events_per_abort_path": 1,
        "additional_work_preparations": 0, "quantum_outputs_retained_by_label_and_exact_factor": True,
        "native_initial_coherence_claimed": False, "C_E_coherence_discrimination_claimed": False,
        "native_loads_identified_with_fermion_occupations": False,
        "work_law_and_readout_architecture_are_declared": True,
        "physical_particle_or_clock_identification": False,
    }

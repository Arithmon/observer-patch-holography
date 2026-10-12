"""Execute preparation, recorded repair and a nondemolition binary readout.

Every coherent operation below comes from a source pulse word. Resets and
classical record/pointer measurements use the declared code services. The exact prediction
package is consulted only after execution, for independent comparisons and
the unchanged count decision rule. No quantum branch is renormalized.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction
import hashlib
import math
from pathlib import Path

import numpy as np

from model.observer.source_scalar import compiler as source
from model.observer.scalar_response.experiment import canonical, compare_receipt
from model.predictions import observer_instrument as target, protocol
from . import circuit

ROOT = Path(__file__).resolve().parents[3]
SEED = 2026101203
TOLERANCE = 2e-10


def _state(value):
    if np.ma.isMaskedArray(value):
        raise ValueError("a finite unmasked source state or reference-column array required")
    result = np.asarray(value)
    if (result.ndim not in (1, 2) or result.shape[0] != 512
            or result.dtype.kind not in "iufc" or not np.all(np.isfinite(result))):
        raise ValueError("a finite source state on the nine-qubit code required")
    return np.array(result, dtype=complex, copy=True)


@lru_cache(maxsize=1)
def _coin():
    result = source.execute([source.pulse("Y", -math.pi/4, [0])], 1, [1., 0.])
    result.flags.writeable = False
    return result


@lru_cache(maxsize=4)
def _repair(candidate, coin):
    result = source.execute(circuit.repair_word(candidate, coin), 3+circuit.RECORD_BITS[candidate])
    result.flags.writeable = False
    return result


def repair_branches(candidate, state):
    """Apply emitted local dilations, retaining all record codes and references.

    The spectator register is a literal untouched tensor factor. Record
    padding has explicit abort labels, even when its ideal probability is zero.
    """
    if candidate not in circuit.RECORD_BITS:
        raise ValueError("candidate must be C or E")
    state = _state(state)
    vector = state.ndim == 1
    columns = state.reshape(512, -1)
    count = columns.shape[1]
    local = columns.reshape(64, 8, count).transpose(1, 0, 2).reshape(8, 64*count)
    records = 1 << circuit.RECORD_BITS[candidate]
    blank = np.zeros((8*records, 64*count), complex)
    blank[:8] = local
    branches = []
    for coin in (0, 1):
        inverse = circuit.record_values(candidate, coin)
        executed = (_repair(candidate, coin)@blank).reshape(records, 8, 64, count)
        for code in range(records):
            amplitude = executed[code].transpose(1, 0, 2).reshape(512, count)*_coin()[coin]
            value = inverse.get(code)
            current = value if candidate == "C" or value is None else (
                value-(2*coin-1)*(abs(value) % 2))//2
            branches.append({"coin": coin, "record_code": code, "current": current,
                             "contrast": value if candidate == "E" else None,
                             "abort": value is None,
                             "amplitude": amplitude[:, 0] if vector else amplitude})
    return branches


@lru_cache(maxsize=1)
def _readout_local_columns():
    # The fixed preparations and local repairs occupy only spectator label 0.
    # Cache ACTUAL word execution on its eight source basis columns. This is
    # linear reuse of computed amplitudes, not insertion of the target projector.
    blank = np.zeros((1024, 8), complex)
    blank[:8] = np.eye(8)
    result = source.execute(circuit.readout_word(), 10, blank)
    result.flags.writeable = False
    return result


def readout_branches(state):
    """Read only the emitted pointer; keep the entire minus subspace coherent."""
    state = _state(state)
    vector = state.ndim == 1
    columns = state.reshape(512, -1)
    if np.count_nonzero(columns[8:]) == 0:
        output = _readout_local_columns()@columns[:8]
    else:
        blank = np.concatenate((columns, np.zeros_like(columns)))
        output = source.execute(circuit.readout_word(), 10, blank)
    return {outcome: block[:, 0] if vector else block for outcome, block in
            ((-1, output[:512]), (1, output[512:]))}


def events(candidate, setting):
    """All unnormalized terminal branches, including preparation and old records."""
    if candidate not in circuit.RECORD_BITS or setting not in circuit.SETTINGS:
        raise ValueError("a fixed candidate and preparation setting required")
    components = [(None, setting, 1.)]
    if setting == "diagonal_control":
        components = [(0, "phase_0", _coin()[0]), (1, "phase_pi", _coin()[1])]
    result = []
    for preparation_record, pure_setting, weight in components:
        blank = np.zeros(512)
        blank[0] = 1
        prepared = weight*source.execute(circuit.preparation_word(pure_setting), 9, blank)
        for branch in repair_branches(candidate, prepared):
            for outcome, amplitude in readout_branches(branch["amplitude"]).items():
                result.append({**branch, "preparation_record": preparation_record,
                               "outcome": outcome, "amplitude": amplitude})
    return result


def probability(event):
    return float(np.vdot(event["amplitude"], event["amplitude"]).real)


def _check_events(candidate, setting, executed):
    total = sum(map(probability, executed))
    if not math.isfinite(total) or abs(total-1) > TOLERANCE:
        raise ArithmeticError("complete source instrument is not normalized")
    expected = target.event_outputs(candidate, setting)
    maximum = 0.
    for event in expected:
        record = event["record"]
        selected = [row for row in executed if not row["abort"]
                    and row["coin"] == record["coin"] and row["current"] == record["current"]
                    and (candidate == "C" or row["contrast"] == record["contrast"])
                    and row["outcome"] == event["outcome"]]
        actual = sum((np.outer(row["amplitude"], row["amplitude"].conj())
                      for row in selected), np.zeros((512, 512), complex))
        ideal = np.zeros_like(actual)
        for item in event["source_density"]:
            ideal[circuit.SOURCE_CODES[item["row"]], circuit.SOURCE_CODES[item["column"]]] = complex(
                float(Fraction(item["real"])), float(Fraction(item["imaginary"])))
        maximum = max(maximum, float(np.linalg.norm(actual-ideal, "fro")))
    abort = sum(probability(row) for row in executed if row["abort"])
    if maximum > TOLERANCE or abort > TOLERANCE:
        raise ArithmeticError("source-executed branch densities differ from the exact instrument")
    return {"raw_total_probability": total, "abort_probability": abort,
            "maximum_branch_density_frobenius_error": maximum,
            "plus_probability": sum(probability(row) for row in executed if row["outcome"] == 1 and not row["abort"])}


@lru_cache(maxsize=8)
def resources(candidate, setting):
    if candidate not in circuit.RECORD_BITS or setting not in circuit.SETTINGS:
        raise ValueError("a fixed candidate and setting required")
    # Both repair words have equal resource costs. The pi randomizer branch
    # maximizes integrated pulse norm; zero-angle requests remain charged.
    pure_setting = "phase_pi" if setting == "diagonal_control" else setting
    word = ([source.pulse("Y", -math.pi/4, [0])]
            + circuit.preparation_word(pure_setting)+circuit.repair_word(candidate, 0)+circuit.readout_word())
    if setting == "diagonal_control":
        word += [source.pulse("Y", -math.pi/4, [0])]
    result = source.resources(word)
    result.update({"source_code_qubits": 9,
                   "maximum_work_qubits": 9+circuit.RECORD_BITS[candidate],
                   "maximum_live_code_qubits": 9+circuit.RECORD_BITS[candidate],
                   "retained_classical_bits": 2+circuit.RECORD_BITS[candidate]+int(setting == "diagonal_control"),
                   "blank_resets": 4+int(setting == "diagonal_control"),
                   "record_readout_channels": 3+int(setting == "diagonal_control"),
                   "processor_dimension": 6,
                   "resource_branch": "phase_pi for diagonal control; either fair-coin branch; local addresses are reused between stages",
                   "physical_loading_and_service_errors_certified": False})
    return result


def run(seed=SEED):
    """Reproduce source states, then feed source-generated counts to the fixed test."""
    if type(seed) is not int:
        raise ValueError("an integer seed required")
    from .certificate import (channel_certificate, readout_certificate, parity_refinement_certificate,
                              apparatus_error_contract, ERROR_TERMS)
    plan = protocol.checked_protocol()
    rng = np.random.default_rng(seed)
    records = {}
    for candidate in circuit.RECORD_BITS:
        data = protocol.data_template(plan)
        data["data_kind"] = "synthetic"
        data["provenance"] = f"Draws from emitted source-word amplitudes; seed {seed}; no laboratory data."
        rows = []
        for count_row in data["settings"]:
            setting = count_row["setting"]
            executed = events(candidate, setting)
            diagnostics = _check_events(candidate, setting, executed)
            probabilities = np.array([probability(event) for event in executed])
            total = float(probabilities.sum())
            # Only the classical sampler needs a normalized floating vector.
            # Its correction is exposed; raw quantum states above are untouched.
            draws = rng.multinomial(protocol.TRIALS, probabilities/total)
            counts = dict.fromkeys(protocol.BINS, 0)
            aborts = 0
            branches = []
            for event, chance, count in zip(executed, probabilities, draws, strict=True):
                if event["abort"]:
                    aborts += int(count)
                else:
                    counts[event["coin"], event["current"], event["outcome"]] += int(count)
                branches.append({key: value for key, value in event.items() if key != "amplitude"}
                                | {"probability": float(chance), "count": int(count)})
            count_row["aborts"] = aborts
            for event in count_row["events"]:
                event["count"] = counts[event["coin"], event["current"], event["outcome"]]
            rows.append({"setting": setting, **diagnostics, "resources": resources(candidate, setting),
                         "classical_sampling_mass_correction": abs(total-1), "branches": branches})
        verdict = protocol.evaluate(data, plan)
        for name, decision in verdict["candidate_verdicts"].items():
            expected = "not_rejected" if name == candidate else "rejected_under_contract"
            if decision["verdict"] != expected:
                raise ArithmeticError("source-sampled discrimination control failed")
        records[candidate] = {"channel_certificate": channel_certificate(candidate), "executions": rows,
                              "counts": data, "verdict": verdict}
    pins = ["model/observer/native.py", "model/observer/source_scalar/compiler.py",
            "model/predictions/observer_instrument.py", "model/predictions/protocol.json"]
    worst = resources("E", "diagonal_control")
    allocation = Fraction(1, 64*(3*worst["native_pulses"]+len(ERROR_TERMS)-3))
    required_bounds = dict.fromkeys(ERROR_TERMS, allocation)
    return {"schema": "oph.source_native_discrimination.v1", "seed": seed,
            "source": source.source_summary(), "sources": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in pins},
            "protocol_sha256": protocol.digest(plan), "candidates": records,
            "readout_certificate": readout_certificate(),
            "parity_refinement_certificate": parity_refinement_certificate(),
            "illustrative_required_error_allocation": {
                "each_bound": str(allocation), "independently_established": False,
                "contract": apparatus_error_contract(worst, required_bounds)},
            "all_pass": True, "scope": {"construction": "encoded realization of both complete native record instruments and the fixed preparation/readout test",
                "encoding": "nine-qubit embedding of 78 occupations; physical occupation-to-code loading is not asserted",
                "quantum_execution": "actual full six-level source pulses; complete old-record and pointer branches, without quantum-state renormalization",
                "physical_services": "proper-code load/export, blanks, classical recording/control, storage and computational readout supplied",
                "numerics": "floating execution tested against independent exact channels; no certified roundoff or laboratory calibration",
                "comparison": "synthetic counts only; unchanged prospective protocol; no empirical confirmation or OPH-wide verdict"}}

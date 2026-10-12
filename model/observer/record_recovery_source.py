"""Restore native branch translations with the existing source inverse words.

The original coin and repair record are already classical and remain on the
archive.  CNOTs copy the repair bits into blank scratch; the inverse repair
word restores the source and clears that scratch.  Coherences erased by the
earlier record commitment are not restored.  Only seam (0,1) is executed;
other seams require a separately implemented port relabelling.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from . import native
from .source_discrimination import circuit
from .source_scalar import compiler as source


TOLERANCE = F(1, 10**10)


def _candidate(candidate):
    if type(candidate) is not str or candidate not in ("C", "E"):
        raise ValueError("candidate C or E required")
    return candidate


def _coin(coin):
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("an integer placement coin zero or one required")
    return coin


def _record(candidate, coin, record):
    candidate, coin = _candidate(candidate), _coin(coin)
    if type(record) is not int or record not in circuit.record_values(candidate, coin):
        raise ValueError("an admitted raw repair record required; padding is an abort")
    return record


def inverse_word(candidate, coin):
    """Fresh chronological inverse of the existing phase-exact repair word."""
    return source.inverse(circuit.repair_word(_candidate(candidate), _coin(coin)))


def copy_word(candidate):
    """Copy an already classical archive into blank scratch, preserving it.

    In this local word scratch occupies bits 0..b-1 and the repair archive
    bits b..2b-1.  The externally retained coin selects the inverse word.
    This operation does not assert copying of an unknown quantum record.
    """
    bits = circuit.RECORD_BITS[_candidate(candidate)]
    return [pulse for bit in range(bits) for pulse in source.cnot(bits + bit, bit)]


def _forward_index(candidate, coin, address):
    for left, right in circuit.REPAIR_SWAPS["C", coin]:
        if address % 32 == left:
            address += right-left
        elif address % 32 == right:
            address += left-right
    if candidate == "E" and (address >> coin) & 1:
        address ^= 32
    return address


def _inverse_index(candidate, coin, address):
    if candidate == "E" and (address >> coin) & 1:
        address ^= 32
    for left, right in reversed(circuit.REPAIR_SWAPS["C", coin]):
        if address % 32 == left:
            address += right-left
        elif address % 32 == right:
            address += left-right
    return address


def exact_certificate(candidate):
    """Every recorded matrix unit restores precisely its input-history block.

    The coefficient 1/2 is the squared fair-coin amplitude.  Different raw
    records remain separate blocks throughout; no coherent sum across them
    is used.  Matrix-unit equality also proves arbitrary-reference equality.
    """
    candidate = _candidate(candidate)
    states = native.states()
    if len(states) != 78 or len(set(states)) != 78:
        raise ArithmeticError("the full 78-state native carrier is required")
    reports = []
    for coin in (0, 1):
        columns = []
        for state in states:
            encoded = circuit.encode(state)
            local = _forward_index(candidate, coin, encoded % 8)
            record, output = local // 8, local % 8
            _record(candidate, coin, record)
            repaired = circuit.decode((encoded // 8) * 8 + output)
            total = state[0] + state[1]
            expected_endpoints = ((total // 2, (total + 1) // 2) if coin == 0
                                  else ((total + 1) // 2, total // 2))
            expected_label = (state[0] - expected_endpoints[0] if candidate == "C"
                              else state[0] - state[1])
            if (repaired != expected_endpoints + state[2:]
                    or circuit.record_values(candidate, coin)[record] != expected_label):
                raise ArithmeticError("source permutation differs from independent endpoint repair")
            # The archive is a separate classical register: copying it only
            # initializes scratch.  The inverse never acts on the archive.
            restored = _inverse_index(candidate, coin, output + 8 * record)
            if restored != encoded % 8:
                raise ArithmeticError("inverse failed to restore source and blank scratch")
            columns.append((record, circuit.decode((encoded // 8) * 8 + restored)))
        if len(set(columns)) != len(states):
            raise ArithmeticError("record/source restoration must be injective")
        checked, surviving = 0, 0
        for left, x in enumerate(states):
            for right, y in enumerate(states):
                r, restored_x = columns[left]
                s, restored_y = columns[right]
                actual = {(r, restored_x, restored_y): F(1, 2)} if r == s else {}
                # Independently derive the input projector's label.  The
                # target is a pinched input, never the unmeasured input.
                labels = []
                for state in (x, y):
                    total = state[0] + state[1]
                    endpoint = (total + coin) // 2
                    labels.append(state[0] - endpoint if candidate == "C"
                                  else state[0] - state[1])
                expected = {(r, x, y): F(1, 2)} if labels[0] == labels[1] else {}
                if actual != expected:
                    raise ArithmeticError("recovered recorded channel differs on a matrix unit")
                checked += 1
                surviving += bool(actual)
        reports.append({"coin": coin, "native_inputs": len(states),
                        "matrix_units_checked": checked,
                        "surviving_matrix_units": surviving,
                        "branch_amplitude_squared": "1/2"})
    return {"candidate": candidate, "seam": [0, 1], "coins": reports,
            "restored_instrument": "(1/2) direct_sum_(c,r) P_(c,r) rho P_(c,r)",
            "archive_unchanged": True, "scratch_returns_blank": True,
            "arbitrary_reference_extension": True,
            "cross_record_coherence_restored": False}


@lru_cache(maxsize=4)
def _executed(candidate, coin):
    """Private immutable cache of actual local source-word executions."""
    bits = circuit.RECORD_BITS[candidate]
    forward = source.execute(circuit.repair_word(candidate, coin), 3 + bits)
    inverse = source.execute(inverse_word(candidate, coin), 3 + bits)
    copied = source.execute(copy_word(candidate), 2 * bits)
    for matrix in (forward, inverse, copied):
        matrix.flags.writeable = False
    return forward, inverse, copied


def _recover_local(inverse, copied, record, endpoints):
    """Contract executed copy and inverse, retaining all numerical leakage."""
    records = inverse.shape[0] // 8
    after_copy = np.kron(copied[:, records * record], endpoints)
    # Flattened order is endpoint, scratch, archive (little endian).
    return (inverse @ after_copy.reshape(records, 8 * records).T).T.reshape(-1)


def recover_local(candidate, coin, record, endpoints):
    """Executed local recovery with the supplied committed classical label.

    Return endpoint+scratch+archive amplitudes.  This helper rejects invalid
    labels; a complete instrument must retain their raw branches as aborts,
    as checked in compiled_certificate.  It never normalizes a branch.
    """
    candidate, coin = _candidate(candidate), _coin(coin)
    record = _record(candidate, coin, record)
    if np.ma.isMaskedArray(endpoints):
        raise ValueError("an unmasked finite endpoint amplitude vector required")
    value = np.asarray(endpoints)
    if value.shape != (8,) or value.dtype.kind not in "iufc" or not np.all(np.isfinite(value)):
        raise ValueError("a finite eight-component endpoint amplitude vector required")
    _, inverse, copied = _executed(candidate, coin)
    return _recover_local(inverse, copied, record, np.array(value, dtype=complex, copy=True))


def compiled_certificate(candidate):
    """Numerical source replay against the exact full local permutations.

    Absolute matrix comparisons guard relative and global phase conventions.
    They are numerical diagnostics, not interval or physical error bounds.
    Record-padding branches are included in the error check as retained aborts.
    """
    candidate = _candidate(candidate)
    bits = circuit.RECORD_BITS[candidate]
    records, local_size = 1 << bits, 8 << bits
    reports = []
    for coin in (0, 1):
        forward, inverse, copied = _executed(candidate, coin)
        expected_forward = np.zeros_like(forward)
        expected_inverse = np.zeros_like(inverse)
        for address in range(local_size):
            expected_forward[_forward_index(candidate, coin, address), address] = 1
            expected_inverse[_inverse_index(candidate, coin, address), address] = 1
        expected_copy = np.zeros_like(copied)
        for address in range(records**2):
            scratch, archive = address % records, address // records
            expected_copy[(scratch ^ archive) + records * archive, address] = 1
        for actual, expected in ((forward, expected_forward), (inverse, expected_inverse),
                                 (copied, expected_copy)):
            if (not np.all(np.isfinite(actual))
                    or np.linalg.norm(actual - expected, "fro") > float(TOLERANCE)):
                raise ArithmeticError("compiled source word violates the full phase-exact permutation")

        # A complete column of the recorded channel includes all raw labels,
        # including padding.  Comparisons never normalize or delete branches.
        columns_checked = branches_checked = aborts_checked = 0
        for state in native.states():
            encoded = circuit.encode(state) % 8
            ideal_forward = _forward_index(candidate, coin, encoded)
            ideal_record = ideal_forward // 8
            for record in range(records):
                endpoint = forward[8 * record:8 * (record + 1), encoded]
                if record not in circuit.record_values(candidate, coin):
                    if np.linalg.norm(endpoint) > float(TOLERANCE):
                        raise ArithmeticError("record-padding abort amplitude exceeds tolerance")
                    aborts_checked += 1
                else:
                    actual = _recover_local(inverse, copied, record, endpoint)
                    expected = np.zeros(8 * records**2, complex)
                    if record == ideal_record:
                        expected[encoded + 8 * records * record] = 1
                    if (not np.all(np.isfinite(actual))
                            or np.linalg.norm(actual - expected) > float(TOLERANCE)):
                        raise ArithmeticError("compiled complete branch failed input-projector recovery")
                branches_checked += 1
            columns_checked += 1
        inverse_resources = source.resources(inverse_word(candidate, coin))
        copy_resources = source.resources(copy_word(candidate))
        counts = {key: inverse_resources[key] + copy_resources[key]
                  for key in ("native_pulses", "proper_code_load_events", "proper_code_export_events")}
        reports.append({"coin": coin, "native_input_columns": columns_checked,
                        "all_raw_branches_checked": branches_checked,
                        "retained_padding_abort_branches_checked": aborts_checked,
                        "full_local_unitary_columns": local_size,
                        "inverse_native_pulses": inverse_resources["native_pulses"],
                        "copy_native_pulses": copy_resources["native_pulses"],
                        "recovery_word_resources": counts})
    return {"candidate": candidate, "coins": reports,
            "numerical_residual_threshold": str(TOLERANCE),
            "full_complex_phase_checked": True,
            "source_work_qubits": 9, "scratch_qubits": bits,
            "retained_classical_bits": 1 + bits,
            "classical_record_copy_bits": bits,
            "copy_cnot_words": bits,
            "blank_scratch_initializations": bits,
            "scratch_reusable_after_ideal_recovery": True,
            "archive_record_rewrites": 0,
            "new_coin_or_measurement_required_for_recovery": False,
            "initialization_storage_and_original_recording_are_supplied_services": True,
            "certified_roundoff_or_physical_error_bound": False}


def audit():
    """Replay the existing source controls; no new hardware claim or seam tape."""
    return {"schema": "oph.observer.record_recovery_source.v1",
            "candidates": {candidate: {"exact": exact_certificate(candidate),
                                       "compiled": compiled_certificate(candidate)}
                           for candidate in ("C", "E")},
            "executed_seam": [0, 1],
            "other_seams": "abstract port relabelling only; no executed addressing permutation",
            "record_requirement": "coin and repair labels were already committed classically",
            "new_protected_quantum_memory": False,
            "new_hardware_execution": False,
            "scope": "branchwise reversal of the native translation, retaining history measurement disturbance"}

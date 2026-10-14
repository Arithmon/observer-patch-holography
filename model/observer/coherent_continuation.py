"""Exact coherent continuations of the binary native repair calibration.

This is a finite law-selection experiment, not the committed N1 observer.
Arc coordinates are (binary load state, seam, placement coin). Real rational
arrays store amplitudes multiplied by sqrt(60); real and imaginary parts and
spectator-reference columns extend by linearity. No source pulses or physical
clock are claimed. The two continuation phases have the SAME first coherent
native event but different later Ramsey probabilities.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from . import native, quantum_mixing

STATES = native.states(binary=True)
INDEX = {state: i for i, state in enumerate(STATES)}
EVENTS = 60
DIMENSION = len(STATES) * EVENTS
ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    "model/observer/coherent_continuation.py", "model/observer/native.py",
    "model/observer/quantum_mixing.py", "model/common_source/slow_geometry.py",
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py",
)


def _vector(values, dimension):
    values = tuple(values)
    if len(values) != dimension or any(type(x) not in (int, F) for x in values):
        raise ValueError(f"{dimension} exact real rational amplitudes required")
    return tuple(map(F, values))


def _index(value, dimension):
    if type(value) is not int or not 0 <= value < dimension:
        raise ValueError("index outside the declared finite carrier")
    return value


@lru_cache(maxsize=1)
def branches():
    """Actual N1 branches: (output index, seam index, coin, signed current)."""
    result = []
    for state in STATES:
        for edge, seam in enumerate(quantum_mixing.seams()):
            for coin in (0, 1):
                output, record = native.native(state, seam, coin)
                if (type(record) is not tuple or len(record) != 2
                        or any(type(v) is not int for v in record)
                        or record != (coin, state[seam[0]] - output[seam[0]])):
                    raise ArithmeticError("native record does not encode its actual current")
                if output not in INDEX:
                    raise ArithmeticError("native repair leaves the binary shell")
                result.append((INDEX[output], edge, *record))
    if len(result) != DIMENSION or len(set(result)) != DIMENSION:
        raise ArithmeticError("native coherent columns do not have disjoint event support")
    return tuple(result)


def encode_record(record):
    """Re-encode a valid outgoing N1 record as its reversed arc, coherently.

    The current is not erased by measurement: it remains recoverable from the
    arc. This re-encoding is an added workspace service, not record commitment.
    """
    if (type(record) is not tuple or len(record) != 4
            or any(type(v) is not int for v in record)):
        raise ValueError("an exact outgoing (state, seam, coin, current) tuple required")
    output, edge, coin, current = record
    _index(output, len(STATES))
    _index(edge, len(quantum_mixing.seams()))
    _index(coin, 2)
    if current not in (-1, 0, 1):
        raise ValueError("invalid binary current")
    reverse_coin = coin if current == 0 else 1 - coin
    position = output * EVENTS + edge * 2 + reverse_coin
    back, back_edge, back_coin, back_current = branches()[position]
    recovered_coin = back_coin if back_current == 0 else 1 - back_coin
    if (back_edge, recovered_coin, -back_current) != (edge, coin, current):
        raise ValueError("record is not an attained binary native branch")
    # The inverse output is the input; no extra stored load configuration.
    if branches()[back * EVENTS + edge * 2 + coin] != record:
        raise ValueError("record has no native branch inverse")
    return position


def decode_arc(position):
    position = _index(position, DIMENSION)
    _, edge, reverse_coin, reverse_current = branches()[position]
    coin = reverse_coin if reverse_current == 0 else 1 - reverse_coin
    return position // EVENTS, edge, coin, -reverse_current


@lru_cache(maxsize=1)
def reversal():
    result = tuple(encode_record(record) for record in branches())
    if any(result[result[i]] != i for i in range(DIMENSION)):
        raise ArithmeticError("native branch reversal is not an involution")
    if any(decode_arc(result[i]) != record for i, record in enumerate(branches())):
        raise ArithmeticError("re-encoding loses a native record")
    return result


def embed(values):
    """Uniform coherent event preparation A, in sqrt(60)-scaled coordinates."""
    return tuple(value for value in _vector(values, len(STATES)) for _ in range(EVENTS))


def compress(values):
    values = _vector(values, DIMENSION)
    return tuple(sum(values[i:i + EVENTS], F(0)) / EVENTS
                 for i in range(0, DIMENSION, EVENTS))


def advance(values, phase=-1):
    """W_phase = S [Pi + phase (I-Pi)], for phase in {-1,+1}."""
    if type(phase) is not int or phase not in (-1, 1):
        raise ValueError("a real unitary continuation phase -1 or +1 required")
    values = _vector(values, DIMENSION)
    means = compress(values)
    reflected = tuple(means[i // EVENTS] + phase * (value - means[i // EVENTS])
                      for i, value in enumerate(values))
    return tuple(reflected[i] for i in reversal())


def norm_squared(values):
    values = _vector(values, DIMENSION)
    return sum((v * v for v in values), F(0)) / EVENTS


def ramsey_outputs(values, steps, phase=-1):
    """Both unnormalized work outputs of an added controlled-walk Ramsey read.

    The reference arm stores A|f>; the other applies the stated continuation.
    Summing both outcomes preserves input norm. No outcome is postselected.
    """
    if type(steps) is not int or steps < 0:
        raise ValueError("a nonnegative integer number of continuation steps required")
    if type(phase) is not int or phase not in (-1, 1):
        raise ValueError("a real unitary continuation phase -1 or +1 required")
    initial = current = embed(values)
    for _ in range(steps):
        current = advance(current, phase)
    return tuple(tuple((a + sign * b) / 2 for a, b in zip(initial, current, strict=True))
                 for sign in (1, -1))


def kernel_counts(record_phase=False):
    """60 A_eta^* S A_eta; eta is a sign on one committed record label.

    eta=-1 only for seam (0,1) and current +1. A phase constant on each
    classical Kraus label leaves that complete committed instrument unchanged.
    """
    if type(record_phase) is not bool:
        raise ValueError("record-phase switch must be Boolean")
    return _kernel_counts(record_phase)


@lru_cache(maxsize=2)
def _kernel_counts(record_phase):
    seams = quantum_mixing.seams()
    if (0, 1) not in seams:
        raise ArithmeticError("declared record-phase control seam is absent")
    signs = tuple(-1 if record_phase and seams[e] == (0, 1) and j == 1 else 1
                  for _, e, _, j in branches())
    result = [[0] * len(STATES) for _ in STATES]
    for arc, back in enumerate(reversal()):
        result[arc // EVENTS][back // EVENTS] += signs[arc] * signs[back]
    return tuple(map(tuple, result))


def kernel_apply(values):
    values = _vector(values, len(STATES))
    return tuple(sum((count * value for count, value in zip(row, values, strict=True)), F(0)) / EVENTS
                 for row in kernel_counts())


def _trace_cube(matrix):
    return sum(matrix[i][j] * matrix[j][k] * matrix[k][i]
               for i in range(len(matrix)) for j in range(len(matrix))
               for k in range(len(matrix)) if matrix[i][j] and matrix[j][k])


def audit():
    records, shift, counts = branches(), reversal(), kernel_counts()
    direct = [[0] * len(STATES) for _ in STATES]
    groups = {}
    for arc, (output, edge, coin, current) in enumerate(records):
        direct[arc // EVENTS][output] += 1
        sign = -1 if quantum_mixing.seams()[edge] == (0, 1) and current == 1 else 1
        groups.setdefault((edge, coin, current), set()).add(sign)
    if (tuple(map(tuple, direct)) != counts or counts != tuple(zip(*counts, strict=True))
            or any(sum(row) != EVENTS for row in counts)):
        raise ArithmeticError("coherent discriminant differs from the symmetric native kernel")
    if any(len(signs) != 1 for signs in groups.values()):
        raise ArithmeticError("record phase changes a committed Kraus branch")

    lap = [[0] * 12 for _ in range(12)]
    for a, b in quantum_mixing.seams():
        lap[a][a] += 1
        lap[b][b] += 1
        lap[a][b] -= 1
        lap[b][a] -= 1
    loads = tuple(tuple(6 * x[p] - 1 for x in STATES) for p in range(12))
    for p, column in enumerate(loads):
        wanted = tuple(F(column[i]) - sum(F(lap[p][q] * loads[q][i], EVENTS)
                                         for q in range(12)) for i in range(len(STATES)))
        if kernel_apply(column) != wanted:
            raise ArithmeticError("native load modes do not intertwine with I-L/60")

    # Execute rational arc dynamics independently of the compressed recurrence.
    initial = loads[0]
    current, previous, target = embed(initial), None, tuple(map(F, initial))
    recurrences = []
    for step in range(9):
        actual = compress(current)
        if actual != target or norm_squared(current) != sum(F(x*x) for x in initial):
            raise ArithmeticError("coherent compression recurrence or norm failed")
        recurrences.append({"step": step, "port_zero_pairing": str(sum(a*b for a, b in zip(initial, actual, strict=True)))})
        following = kernel_apply(target)
        if previous is not None:
            following = tuple(2*a-b for a, b in zip(following, previous, strict=True))
        previous, target, current = target, following, advance(current)

    basis = tuple(int(i == 0) for i in range(len(STATES)))
    prepared = embed(basis)
    if advance(prepared, 1) != advance(prepared, -1):
        raise ArithmeticError("the two laws do not have the same first coherent event")
    two_plus, two_minus = advance(advance(prepared, 1), 1), advance(advance(prepared, -1), -1)
    k2 = kernel_apply(kernel_apply(basis))
    if two_plus != prepared or compress(two_minus) != tuple(2*a-b for a,b in zip(k2,basis,strict=True)):
        raise ArithmeticError("two-step continuation phase control failed")
    probabilities = []
    for phase in (1, -1):
        for steps in (0, 1, 2):
            outputs = ramsey_outputs(basis, steps, phase)
            pair = tuple(norm_squared(v) for v in outputs)
            if min(pair) < 0 or sum(pair) != 1:
                raise ArithmeticError("complete Ramsey output is not normalized")
            probabilities.append({"phase": phase, "steps": steps, "plus": str(pair[0]), "minus": str(pair[1])})
    if probabilities[1] | {"phase": -1} != probabilities[4] or probabilities[2]["plus"] == probabilities[5]["plus"]:
        raise ArithmeticError("Ramsey read does not distinguish only the continuation")

    twisted = kernel_counts(True)
    trace_difference = F(_trace_cube(twisted) - _trace_cube(counts), EVENTS**3)
    if trace_difference != -F(1, 900):
        raise ArithmeticError("same-instrument phase control lost its spectral distinction")

    full_input = (2,) + (0,) * 11
    full_output, full_record = native.native(full_input, (0, 1), 0)
    if any(native.native(full_output, seam, coin)[0] == full_input
           for seam in quantum_mixing.seams() for coin in (0, 1)):
        raise ArithmeticError("the full-carrier one-way transition unexpectedly reverses")
    return {
        "branch_space": {"binary_states": len(STATES), "arcs": DIMENSION,
                         "self_arcs": sum(i == j for i, j in enumerate(shift)),
                         "reversed_pairs": sum(i < j for i, j in enumerate(shift)),
                         "complete_record_reencoding_bijective": True,
                         "reversal_rule": "coin unchanged for j=0, flipped otherwise; current becomes -j"},
        "discriminant": {"equal_to_actual_native_kernel": True,
                         "all_load_columns": 12, "load_intertwiner": "K G = G (I-L/60), G[x,p]=6*x[p]-1",
                         "compression": "A* W_minus^n A = T_n(K)",
                         "load_second_difference": "G[n+1]-2*G[n]+G[n-1] = -G[n]*L/30, G[n]=A* W_minus^n A G",
                         "amplitude_not_mean_or_probability": True,
                         "exact_arc_recurrence_samples": recurrences},
        "continuation_control": {"same_first_coherent_native_event": True,
                                 "two_step_compressions": {"phase_plus": "I", "phase_minus": "2*K^2-I"},
                                 "basis_zero_two_step_return_amplitude": str(compress(two_minus)[0]),
                                 "basis_zero_two_step_return_probability": str(compress(two_minus)[0]**2),
                                 "complete_ramsey_probabilities": probabilities},
        "committed_instrument_phase_control": {"phase_rule": "-1 on seam (0,1), current +1; +1 otherwise",
                                               "each_CP_record_branch_unchanged": True,
                                               "changed_undirected_kernel_entries": sum(counts[i][j] != twisted[i][j] for i in range(len(STATES)) for j in range(i+1,len(STATES))),
                                               "trace_D_cubed_minus_trace_K_cubed": str(trace_difference)},
        "full_carrier_counterexample": {"input": list(full_input), "output": list(full_output),
                                         "record": list(full_record), "any_native_one_event_reverse": False},
        "all_pass": True,
    }


def selection_audit():
    """Structural criteria on a fixed coherent experiment do not select a law.

    Full event-permutation covariance is a conditional stronger restriction,
    not a source premise. Physical locality, A1--A3 admission, a selected
    vacuum and a calibrated time reversal are not claimed by these checks.
    """
    from ..common_source.slow_geometry import source

    records, shift, counts = branches(), reversal(), kernel_counts()
    seams = quantum_mixing.seams()
    edge_index = {edge: i for i, edge in enumerate(seams)}
    signs = tuple(-1 if current else 1 for _, _, _, current in records)
    signed = [[0] * len(STATES) for _ in STATES]
    groups = ({}, {})
    for arc, (output, edge, coin, current) in enumerate(records):
        x, y = STATES[arc // EVENTS], STATES[output]
        a, b = seams[edge]
        if sum(x) != sum(y) or any(x[p] != y[p] for p in range(12) if p not in (a, b)):
            raise ArithmeticError("continuation support is not a conservative native seam move")
        signed[arc // EVENTS][output] += signs[arc]
        for table, label in zip(groups, (current, x[a]-x[b]), strict=True):
            table.setdefault((edge, coin, label), set()).add(signs[arc])
    if (any(signs[i] != signs[shift[i]] for i in range(DIMENSION))
            or any(len(v) != 1 for table in groups for v in table.values())):
        raise ArithmeticError("signed reversal fails involution or changes a committed CP branch")
    if any(signed[i][j] != (counts[i][j] if i == j else -counts[i][j])
           for i in range(len(STATES)) for j in range(len(STATES))):
        raise ArithmeticError("signed discriminant differs from 2*diag(K)-K")

    rotations = source().rotations()
    if len(rotations) != 60 or len(set(rotations)) != 60:
        raise ArithmeticError("the pinned source lacks its sixty distinct proper rotations")
    for permutation in rotations:
        state_map = []
        for x in STATES:
            y = [0] * 12
            for old, new in enumerate(permutation):
                y[new] = x[old]
            state_map.append(INDEX[tuple(y)])
        events = tuple(2*edge_index[tuple(sorted((permutation[a], permutation[b])))]
                       + (coin ^ int(permutation[a] > permutation[b]))
                       for a, b in seams for coin in (0, 1))
        transport = tuple(EVENTS*y+event for y in state_map for event in events)
        if (len(set(transport)) != DIMENSION
                or any(transport[shift[i]] != shift[transport[i]] or signs[transport[i]] != signs[i]
                       for i in range(DIMENSION))):
            raise ArithmeticError("native reversal or moving-record phase breaks proper covariance")

    reverse = lambda values: tuple(values[i] for i in shift)
    probe = tuple(F(i % 7-3, 5) for i in range(DIMENSION))
    means = compress(probe)
    disagreement = sum((v-means[i//EVENTS])**2 for i, v in enumerate(probe))
    uniform = (F(1),) * DIMENSION
    coins = []
    for phase in (1, -1):
        diagonal, off_diagonal = F(1+phase*(EVENTS-1), EVENTS), F(1-phase, EVENTS)
        if (diagonal+(EVENTS-1)*off_diagonal != 1
                or diagonal**2+(EVENTS-1)*off_diagonal**2 != 1
                or 2*diagonal*off_diagonal+(EVENTS-2)*off_diagonal**2 != 0):
            raise ArithmeticError("permutation-covariant coin is not orthogonal and uniform-fixing")
        coin_output = reverse(advance(probe, phase))
        expected = tuple((diagonal-off_diagonal)*v+EVENTS*off_diagonal*means[i//EVENTS]
                         for i, v in enumerate(probe))
        if (coin_output != expected or compress(coin_output) != means
                or sum((v-means[i//EVENTS])**2 for i, v in enumerate(coin_output)) != disagreement
                or advance(reverse(advance(probe, phase)), phase) != reverse(probe)
                or advance(uniform, phase) != uniform):
            raise ArithmeticError("continuation fails the shared structural selection criteria")
        coins.append({"phase": phase, "diagonal": str(diagonal), "off_diagonal": str(off_diagonal)})

    trace_difference = F(_trace_cube(signed)-_trace_cube(counts), EVENTS**3)
    amplitude = F(sum(signs), DIMENSION)
    remainder = sum((F(sign)-amplitude)**2 for sign in signs)/DIMENSION
    ramsey = tuple(norm_squared(tuple((1+sign*value)/F(2) for value in signs))/len(STATES)
                   for sign in (1, -1))
    if (signs.count(-1) != 600 or trace_difference != -F(1, 90)
            or amplitude != F(23, 33) or amplitude**2+remainder != 1
            or ramsey != (F(28, 33), F(5, 33))):
        raise ArithmeticError("proper-covariant signed source control lost its complete signal")
    return {
        "fixed_coherent_interface": "same 3960-dimensional arc space, embedding A and complete Ramsey read",
        "phase_endpoint_criteria": {"coins": coins, "same_first_coherent_event": True,
                                    "no_continuous_parameter": True, "real_orthogonal": True,
                                    "uniform_coherent_state_fixed": True,
                                    "support_locality": "coin fixes loads; shift makes one conservative native seam move",
                                    "common_antiunitary": "Theta=S conjugation; Theta^2=I; Theta W Theta^-1=W^*",
                                    "disagreement_norm": "||(I-Pi) C v||=||(I-Pi) v||; no strict repair"},
        "conditional_coin_classification": {"premises": "real orthogonal; fixes uniform vector; commutes with every S60 event permutation",
                                            "equations": "C=aI+bJ; a+60b=1; a^2=1",
                                            "solutions": ["I", "2*P_uniform-I"],
                                            "full_event_symmetry_is_a_source_premise": False},
        "proper_covariance": {"rotations": len(rotations), "arc_checks_per_reversal": len(rotations)*DIMENSION,
                              "unsigned_and_signed_reversals": True, "orientation_coin_flip_included": True},
        "moving_record_sign_control": {"negative_arcs": signs.count(-1), "phase_rule": "-1 iff native current is nonzero, on every seam",
                                       "signed_reversal_involutive": True, "each_committed_C_and_E_CP_branch_unchanged": True,
                                       "same_classical_stationary_distribution": True, "uniform_coherent_state_fixed": False,
                                       "same_first_coherent_event": False, "discriminant": "D=2*diag(K)-K",
                                       "trace_D_cubed_minus_trace_K_cubed": str(trace_difference),
                                       "uniform_return_amplitude": str(amplitude),
                                       "complete_uniform_return_probabilities": [str(amplitude**2), str(remainder)],
                                       "complete_same_Ramsey_probabilities": list(map(str, ramsey)),
                                       "unsigned_uniform_read_probabilities": ["1", "0"]},
        "committed_interface_scope": "scalar phases cancel in every CP branch, hence every finite adaptive committed-record history including reference systems; coherent record access is a different interface",
        "A1_A2_A3_source_admission_proved": False,
    }


def run():
    return {"schema": "oph.observer.coherent_continuation.v1",
            "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS},
            "audit": audit(),
            "selection": selection_audit(),
            "contract": {"native_rule": "N1 current-record instrument on its invariant 66-state binary shell",
                         "seams": "the existing thirty seams, with coherent uniform event amplitudes",
                         "workspace": "uncommitted event amplitudes, coherent record re-encoding and reflection",
                         "readout": "added controlled-continuation Ramsey arm, both outcomes and quantum outputs retained",
                         "parameter": "coherent continuation iteration; neither committed native attempts nor physical time",
                         "nearest_theorem": "Szegedy quantum walk / unitary-involution spectral mapping",
                         "spectral_mapping_source": "https://arxiv.org/abs/1506.06457"},
            "scope": {"same_as_existing_committed_observer": False,
                      "coherent_extension_selected_by_native_instrument": False,
                      "no_new_supplied_field_Hamiltonian": True,
                      "supplied_coherent_control_law": True,
                      "source_pulses_emitted": False,
                      "physical_clock_or_spatial_refinement": False,
                      "unknown_quantum_input_restored_after_commitment": False,
                      "physical_predictions_or_empirical_validation": False}}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def read_receipt(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate receipt key: {key}")
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f"nonfinite receipt constant: {value}")
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique, parse_constant=invalid)


def compare_receipt(recorded, actual):
    if canonical(recorded) != canonical(actual):
        raise ValueError("exact coherent continuation receipt mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    encoded = canonical(result)
    if args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    if args.check:
        compare_receipt(read_receipt(args.check), result)
    print(json.dumps({"exact_replay": True, "arcs": DIMENSION,
                      "coherent_extension_selected": False,
                      "result_sha256": hashlib.sha256(encoded.encode()).hexdigest()}, sort_keys=True))


if __name__ == "__main__":
    main()

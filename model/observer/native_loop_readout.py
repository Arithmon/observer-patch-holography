"""Exact loop test for the coin/current information emitted by native repair.

The fixed observer sees two adaptive seams. An additional, explicitly ideal
protocol measures both coin-conditioned current laws at every labelled seam
on fresh identical total-two preparations. Those laws select integer edge
contrasts. This gate asks whether those selected contrasts contain a
nontrivial closed-loop transport, without supplying a link phase.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from model.common_source import slow_geometry
from model.observer import native
from model.observer.source_discrimination import circuit
from model.observer.source_scalar import compiler


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE / "results" / "native_loop_readout.json"


def carrier() -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int, int], ...]]:
    source = slow_geometry.source()
    seams = tuple(source.seams())
    faces = tuple(source.oriented_faces())
    if (len(seams) != 30 or len(faces) != 20 or len(set(seams)) != 30
            or any(a >= b or not 0 <= a < b < 12 for a, b in seams)
            or any(len(set(face)) != 3 or any(not 0 <= x < 12 for x in face)
                   for face in faces)):
        raise ArithmeticError("native carrier incidence changed")
    edge_set = set(seams)
    if any(tuple(sorted((u, v))) not in edge_set
           for face in faces for u, v in zip(face, face[1:] + face[:1])):
        raise ArithmeticError("native face has an absent seam")
    return seams, faces


def cycle_basis(seams: tuple[tuple[int, int], ...]) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, int], ...]]:
    """A deterministic spanning-tree basis of all independent graph loops."""
    neighbors: dict[int, list[int]] = {i: [] for i in range(12)}
    for a, b in seams:
        neighbors[a].append(b)
        neighbors[b].append(a)
    parent = {0: None}
    queue = [0]
    for a in queue:
        for b in sorted(neighbors[a]):
            if b not in parent:
                parent[b] = a
                queue.append(b)
    if len(parent) != 12:
        raise ArithmeticError("native seam graph disconnected")
    tree = {tuple(sorted((a, b))) for a, b in parent.items() if b is not None}
    chords = tuple(edge for edge in seams if edge not in tree)
    cycles = []
    for a, b in chords:
        ancestors_a = []
        node = a
        while node is not None:
            ancestors_a.append(node)
            node = parent[node]
        ancestors_b = []
        node = b
        while node is not None:
            ancestors_b.append(node)
            node = parent[node]
        meet = next(node for node in ancestors_a if node in ancestors_b)
        path_a = ancestors_a[:ancestors_a.index(meet) + 1]
        path_b = ancestors_b[:ancestors_b.index(meet)]
        cycles.append(tuple(path_a + path_b[::-1] + [a]))
    if len(tree) != 11 or len(chords) != 19 or len(cycles) != 19:
        raise ArithmeticError("native cycle rank changed")
    return tuple(cycles), chords


def read_contrasts(loads: tuple[int, ...], seams: tuple[tuple[int, int], ...]) -> dict[tuple[int, int], int]:
    """Derive the integer edge field from both actual native current records."""
    if loads not in native.states():
        raise ValueError("native total-two load required")
    contrasts = {}
    for seam in seams:
        records = [native.native(loads, seam, coin)[1] for coin in (0, 1)]
        if tuple(record[0] for record in records) != (0, 1):
            raise ArithmeticError("native coin labels changed")
        contrasts[seam] = sum(record[1] for record in records)
    return contrasts


def integer_holonomy(contrasts: dict[tuple[int, int], int], closed: tuple[int, ...]) -> int:
    """Oriented additive loop value; each character maps zero to unit phase."""
    if len(closed) < 4 or closed[0] != closed[-1]:
        raise ValueError("a closed loop with at least three seams is required")
    total = 0
    for u, v in zip(closed, closed[1:]):
        edge = tuple(sorted((u, v)))
        if edge not in contrasts:
            raise ValueError("loop uses a seam without a measured contrast")
        total += contrasts[edge] if u < v else -contrasts[edge]
    return total


def binary_transfer(seams: tuple[tuple[int, int], ...]) -> tuple[list[list[int]], tuple[int, int, int]]:
    """Exact 60-branch transfer counts for the separately uniform-seam protocol."""
    states = native.states(binary=True)
    index = {loads: i for i, loads in enumerate(states)}
    counts = [[0] * len(states) for _ in states]
    for i, loads in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                counts[i][index[native.native(loads, seam, coin)[0]]] += 1
    if any(sum(row) != 60 for row in counts):
        raise ArithmeticError("native transfer row lost mass")
    if any(counts[i][j] != counts[j][i]
           for i in range(len(states)) for j in range(len(states))):
        raise ArithmeticError("uniform-seam native transfer is not reversible")
    triangle = next(((i, j, k) for i in range(len(states))
                     for j in range(i + 1, len(states)) if counts[i][j]
                     for k in range(j + 1, len(states))
                     if counts[j][k] and counts[k][i]), None)
    if triangle is None:
        raise ArithmeticError("native transfer triangle witness disappeared")
    return counts, triangle


def same_endpoint_phase_readout() -> dict:
    """Source-execute a declared eraser on the fixed observer's two-event tape.

    This is a phase-sensitive *path-record* control, not a carrier face-loop
    holonomy. The first coin is committed in the fixed observer; the eraser
    requires keeping that register coherent until its final source read.
    """
    loads = tuple(int(i in (1, 11)) for i in range(12))
    paths = [(first_weight * second_weight, last)
             for first_weight, first in native.tick(native.State(loads))
             for second_weight, last in native.tick(first)]
    if len(paths) != 4 or any(weight != Q(1, 4) for weight, _ in paths):
        raise ArithmeticError("fixed native path probabilities changed")

    # Little-endian tape: c1, two C current bits, c2, two C current bits.
    # The current-code dictionary is the one used by the source-executed C
    # dilation. The final source state is retained as a separate fibre.
    def code(word: tuple[tuple[int, int], ...]) -> int:
        return sum((coin << (3 * event))
                   + (circuit.CURRENT_CODES[current] << (3 * event + 1))
                   for event, (coin, current) in enumerate(word))

    endpoints: dict[tuple[int, ...], np.ndarray] = {}
    for weight, last in paths:
        assert weight == Q(1, 4)
        vector = endpoints.setdefault(last.loads, np.zeros(64, dtype=complex))
        record = code(last.prefix)
        if vector[record] != 0:
            raise ArithmeticError("native record histories collided before readout")
        vector[record] = .5  # the native isometry's positive square-root weight
    returning = [(code(last.prefix), last.prefix)
                 for _, last in paths if last.loads == loads]
    if returning != [
        (0, ((0, 0), (0, 0))),
        (21, ((1, -1), (0, 1))),
    ]:
        raise ArithmeticError("same-endpoint native record words changed")

    # Recode only the second returning word to |1>, then use the ordinary
    # source-compiled real mixer to read |0> +/- |1>. The scalar first-coin
    # phase is another executable source word. Neither service is selected
    # by the committed native observer.
    phase = [compiler.pulse("Z", math.pi / 2, [0]),
             compiler.pulse("I", -math.pi / 2, [0])]
    eraser = (circuit.basis_swap_word(returning[1][0], 1, 6)
              + [compiler.pulse("Y", math.pi / 4, [0])])
    source_rows = {}
    for label, word in (("base", eraser), ("first_coin_phase", phase + eraser)):
        outputs = {endpoint: compiler.execute(word, 6, vector)
                   for endpoint, vector in endpoints.items()}
        mass = sum(float(np.vdot(v, v).real) for v in outputs.values())
        target = (Q(1, 2), Q(0)) if label == "base" else (Q(0), Q(1, 2))
        returned = tuple(float(abs(outputs[loads][bit]) ** 2) for bit in (0, 1))
        remainder = sum(float(abs(value) ** 2)
                        for endpoint, vector in outputs.items()
                        for bit, value in enumerate(vector)
                        if endpoint != loads or bit not in (0, 1))
        if (abs(mass - 1) > 2e-10 or
                max(abs(returned[bit] - float(target[bit])) for bit in (0, 1)) > 2e-10 or
                abs(remainder - .5) > 2e-10):
            raise ArithmeticError("source pulse word lost a path or the phase contrast")
        source_rows[label] = {
            "exact_returning_plus_minus_probabilities": list(map(str, target)),
            "other_endpoint_probability": "1/2",
            "all_outcome_probability_mass_within_2e_10": abs(mass - 1) <= 2e-10,
            "source_readout_probability_residual_below_2e_10": all(
                abs(returned[bit] - float(target[bit])) <= 2e-10 for bit in (0, 1)),
        }
    # First-coin Z multiplies each complete committed branch by a scalar.
    # This identity holds for arbitrary source/reference input and survives
    # any later adaptive classical control on the first recorded coin.
    committed = {last.prefix: weight for weight, last in paths}
    if len(committed) != 4 or sum(committed.values(), Q(0)) != 1:
        raise ArithmeticError("the complete native current law changed")
    phase_branch_counts = {}
    for law in ("C", "E"):
        signs: dict[tuple[int, ...], set[int]] = {}
        branches = 0
        for initial in native.states():
            for _, first in native.tick(native.State(initial)):
                first_seam = native.E
                first_label = (first.prefix[0][1] if law == "C" else
                               initial[first_seam[0]] - initial[first_seam[1]])
                second_seam = native.selected_seam(first)
                for _, last in native.tick(first):
                    second_label = (last.prefix[1][1] if law == "C" else
                                    first.loads[second_seam[0]] - first.loads[second_seam[1]])
                    key = (first.prefix[0][0], first_label, *second_seam,
                           last.prefix[1][0], second_label)
                    signs.setdefault(key, set()).add(1 - 2 * first.prefix[0][0])
                    branches += 1
        if branches != 312 or any(len(values) != 1 for values in signs.values()):
            raise ArithmeticError("first-coin phase is not scalar on a complete committed branch")
        phase_branch_counts[law] = branches
    return {
        "first_and_second_seams_on_both_returning_paths": [list(native.E), list(native.E)],
        "returning_tape_words": [[list(record) for record in word] for _, word in returning],
        "returning_tape_codes": [value for value, _ in returning],
        "all_four_committed_word_probabilities": {
            str(code(word)): str(weight) for word, weight in sorted(committed.items())},
        "first_coin_scalar_phase_leaves_every_committed_CP_branch_unchanged": True,
        "complete_C_E_adaptive_branches_checked": phase_branch_counts,
        "source_pulse_counts": {"first_coin_phase": len(phase), "record_eraser": len(eraser)},
        "source_executed_eraser": source_rows,
        "read_code_is_a_declared_coherent_tape_analyzer_not_a_native_current": True,
        "phase_selected_by_committed_native_records": False,
        "carrier_face_loop_or_nontrivial_link_holonomy": False,
        "reversed_seam_pair_has_unit_wilson_product_for_any_U1_link_assignment": True,
        "physical_defect_or_mass_identification": False,
    }


def run() -> dict:
    seams, faces = carrier()
    cycles, chords = cycle_basis(seams)
    binding = slow_geometry.registered_binding()
    if not binding["all_pass"] or binding["registered_gram_entries_checked_exactly"] != 144:
        raise ArithmeticError("native port chart lost its registered carrier binding")
    all_states = native.states()
    for loads in all_states:
        contrasts = read_contrasts(loads, seams)
        if any(contrasts[a, b] != loads[a] - loads[b] for a, b in seams):
            raise ArithmeticError("native branch records are not a port-potential gradient")
        if any(integer_holonomy(contrasts, loop) != 0 for loop in cycles):
            raise ArithmeticError("native edge contrast has a nonzero fundamental loop")
        if any(integer_holonomy(contrasts, face + (face[0],)) != 0 for face in faces):
            raise ArithmeticError("native edge contrast has a nonzero face loop")

    witness_loads = tuple(int(i in (0, 1)) for i in range(12))
    witness_contrasts = read_contrasts(witness_loads, seams)
    corrupted = dict(witness_contrasts)
    corrupted[chords[0]] += 1
    changed_cycles = [i for i, loop in enumerate(cycles)
                      if integer_holonomy(corrupted, loop) != 0]
    if len(changed_cycles) != 1 or changed_cycles[0] != 0:
        raise ArithmeticError("single-chord flux control was not detected")

    counts, triangle = binary_transfer(seams)
    i, j, k = triangle
    forward = counts[i][j] * counts[j][k] * counts[k][i]
    reverse = counts[i][k] * counts[k][j] * counts[j][i]
    if not forward or forward != reverse:
        raise ArithmeticError("native transition-loop balance changed")
    # This stronger binary witness returns to its own source configuration
    # along distinct native tape paths, one with two nonzero currents.
    collision_loads = tuple(int(i in (1, 11)) for i in range(12))
    collision_paths = {
        last.prefix: first_weight * second_weight
        for first_weight, first in native.tick(native.State(collision_loads))
        for second_weight, last in native.tick(first)
        if last.loads == collision_loads
    }
    expected_paths = (((0, 0), (0, 0)), ((1, -1), (0, 1)))
    if tuple(sorted(collision_paths)) != expected_paths or any(
        weight != Q(1, 4) for weight in collision_paths.values()
    ):
        raise ArithmeticError("same-endpoint record-labelled path witness changed")
    return {
        "schema": "oph.native_loop_readout.v1",
        "status": "NATIVE_CURRENT_CONTRAST_FLAT__COHERENT_HOLONOMY_OPEN",
        "source": {
            "native_rule": "model/observer/native.py: native",
            "carrier_incidence": "model/common_source/slow_geometry.py: source().seams(), oriented_faces()",
            "registered_carrier_binding": binding,
            "mass_or_particle_target_data_consumed": False,
        },
        "ideal_all_seam_readout": {
            "protocol": "both exact coin-conditioned current laws on each labelled seam, with fresh identical preparations",
            "native_states_exhausted": len(all_states),
            "branch_records_checked": len(all_states) * len(seams) * 2,
            "edge_contrast_identity": "j_0(a,b)+j_1(a,b)=n_a-n_b",
            "independent_fundamental_loops": len(cycles),
            "oriented_face_loops": len(faces),
            "all_integer_loop_sums": "0",
            "every_character_holonomy": "1 for every U(1) character of the measured additive integer contrast",
            "formal_connection": (
                "The measured gradient is a pure gauge under every character; "
                "Lean/Screen/SeamU1HolonomyClassification.lean proves the abstract "
                "coboundary_range_eq_chordHolonomyHom_ker criterion on the "
                "registered seam carrier. The Python record-to-gradient check is "
                "not itself a Lean theorem."
            ),
            "single_chord_plus_one_control": {
                "chord": list(chords[0]),
                "nonzero_fundamental_loop_indices": changed_cycles,
                "native_branch_packet": False,
            },
        },
        "uniform_seam_binary_transfer": {
            "protocol": "independent uniform choice of 30 seams and fair native coin",
            "state_count": len(counts),
            "nonzero_unordered_transition_pairs": sum(bool(counts[a][b]) for a in range(len(counts))
                                                      for b in range(a + 1, len(counts))),
            "all_transition_counts_symmetric": True,
            "example_closed_transition_triangle": list(triangle),
            "forward_reverse_count_products": [forward, reverse],
            "all_positive_transition_cycle_affinities": "0 because P(x,y)=P(y,x)",
        },
        "fixed_observer": {
            "first_seam": list(native.E),
            "possible_second_seams": [list(native.E), list(native.F)],
            "distinct_seam_count": 2,
            "enclosed_carrier_face_count": 0,
            "same_endpoint_path_control": {
                "initial_and_final_occupied_ports": [1, 11],
                "separate_retained_record_words": [[list(record) for record in path]
                                                   for path in expected_paths],
                "probability_each": "1/4",
                "relative_phase_readout_in_committed_tape": False,
            },
            "coherent_path_record_control": same_endpoint_phase_readout(),
        },
        "scope": (
            "Exact native current records select a flat additive port contrast, not a "
            "U(1) connection with nontrivial flux. The two-event observer does not "
            "traverse a carrier face. The ideal all-seam protocol and uniform-seam "
            "transfer are separate declared protocols; their flat/balanced outputs "
            "do not rule out a different coherent source-derived holonomy. A "
            "record-accessible relative phase in the fixed committed interface, "
            "geometric path composition, protected flux sector and continuum "
            "spectral pole remain unconstructed. A declared coherent tape "
            "eraser distinguishes scalar branch phases on the two returning "
            "E/E paths; this is not a carrier-loop holonomy."
        ),
        "all_pass": True,
    }


def canonical(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", type=Path, default=DEFAULT_RECEIPT)
    group.add_argument("--write", type=Path)
    args = parser.parse_args()
    encoded = canonical(run())
    if args.write is not None:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    elif args.check.read_text(encoding="utf-8") != encoded:
        raise SystemExit("Native loop-readout receipt mismatch")
    print(json.dumps({"all_pass": True, "sha256": hashlib.sha256(encoded.encode()).hexdigest()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

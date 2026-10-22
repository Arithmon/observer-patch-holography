"""Exact native-link selection gate on the registered twelve-port carrier.

The C/E repair instrument is a recorded CP map.  Multiplying an entire
coin/current (or coin/contrast) Kraus branch by a scalar unit phase leaves
that instrument unchanged on every input, including inputs entangled with a
reference.  This module uses an exact +/-1 witness: one completion assigns
+1 to each moving one-token branch and another assigns -1.  Both have the
same committed repair law, whereas their candidate triangular link products
are +1 and -1.  These link products are *declared coherent completions*, not
measurements of a native holonomy.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

from model.common_source import slow_geometry
from model.observer import native


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE / "results" / "native_link_selection.json"


def _branches(loads, seam, law):
    if law not in ("C", "E"):
        raise ValueError("C or E record law required")
    images = []
    for coin in (0, 1):
        final, (_, current) = native.native(loads, seam, coin)
        label = (coin, current) if law == "C" else (coin, loads[seam[0]] - loads[seam[1]])
        images.append((label, final, current))
    return tuple(images)


def _matrix_unit(columns_left, columns_right, *, negative_moving_phase):
    """Full recorded image; a common sign on each branch cancels exactly."""
    result = {}
    for (label_l, final_l, current_l), (label_r, final_r, current_r) in zip(
            columns_left, columns_right, strict=True):
        if label_l != label_r:
            continue
        if current_l != current_r:
            raise ArithmeticError("equal record labels conceal distinct currents")
        sign_l = -1 if negative_moving_phase and current_l else 1
        sign_r = -1 if negative_moving_phase and current_r else 1
        result[label_l, final_l, final_r] = F(sign_l * sign_r, 2)
    return result


def _one_token_at(port):
    return tuple(int(index == port) for index in range(12))


def moving_branch(seam, source_port, law):
    """Return the unique current-resolved one-token transfer branch."""
    a, b = seam
    if source_port not in seam:
        raise ValueError("one-token source must be an endpoint of the seam")
    target_port = b if source_port == a else a
    coin = 0 if source_port == a else 1
    label, final, current = _branches(_one_token_at(source_port), seam, law)[coin]
    if final != _one_token_at(target_port) or current != (1 if source_port == a else -1):
        raise ArithmeticError("native one-token transfer branch changed")
    competing_label, _, _ = _branches(_one_token_at(target_port), seam, law)[coin]
    if label == competing_label:
        raise ArithmeticError("moving branch is not source-selective")
    return label


def normalized_moving_matrix(seam, source_port, law):
    """The selected branch's one-token Kraus matrix after removing sqrt(1/2)."""
    label = moving_branch(seam, source_port, law)
    coin = label[0]
    matrix = [[0, 0], [0, 0]]
    for column, port in enumerate(seam):
        candidate, final, _ = _branches(_one_token_at(port), seam, law)[coin]
        if candidate == label:
            for row, output_port in enumerate(seam):
                matrix[row][column] = int(final == _one_token_at(output_port))
    return tuple(tuple(row) for row in matrix)


def _multiply_two(left, right):
    return tuple(tuple(sum(left[i][k] * right[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def _link_sign(source_port, target_port, law, negative_moving_phase):
    seam = tuple(sorted((source_port, target_port)))
    label = moving_branch(seam, source_port, law)
    current = label[1] if law == "C" else native.native(
        _one_token_at(source_port), seam, label[0])[1][1]
    return -1 if negative_moving_phase and current else 1


def _loop_sign(loop, law, negative_moving_phase):
    sign = 1
    for source_port, target_port in zip(loop[:-1], loop[1:], strict=True):
        sign *= _link_sign(source_port, target_port, law, negative_moving_phase)
    return sign


def run():
    seams = tuple(slow_geometry.source().seams())
    faces = tuple(slow_geometry.source().oriented_faces())
    states = native.states()
    if len(seams) != 30 or len(faces) != 20 or len(states) != 78:
        raise ArithmeticError("registered carrier or native state count changed")

    # The full 78^2 matrix-unit comparison on one seam includes every total-two
    # input and arbitrary-reference extension.  The branch rule is identical
    # on each labelled seam; verify its columns and both moving directions on
    # all thirty seams as a separate transport check.
    matrix_units = 0
    moving_columns = 0
    for law in ("C", "E"):
        columns = [_branches(loads, seams[0], law) for loads in states]
        for left, right in product(columns, repeat=2):
            if _matrix_unit(left, right, negative_moving_phase=False) != _matrix_unit(
                    left, right, negative_moving_phase=True):
                raise ArithmeticError("a branch scalar phase changed the committed instrument")
            matrix_units += 1
        for seam in seams:
            for port in seam:
                moving_branch(seam, port, law)
                forward = normalized_moving_matrix(seam, port, law)
                reverse = normalized_moving_matrix(
                    seam, seam[1] if port == seam[0] else seam[0], law)
                source_projection = ((1, 0), (0, 0)) if port == seam[0] else ((0, 0), (0, 1))
                if _multiply_two(reverse, forward) != source_projection:
                    raise ArithmeticError("selected moving branch became a two-port unitary")
                moving_columns += 1

    # Both completions obey inverse-link, reversal, and backtrack identities.
    # The all-minus completion is port-label symmetric but gives -1 around
    # every native triangular face, so even those constraints cannot select
    # the source phase from committed C/E records.
    for law, negative in product(("C", "E"), (False, True)):
        for a, b in seams:
            if (_link_sign(a, b, law, negative) * _link_sign(b, a, law, negative) != 1
                    or _loop_sign((a, b, a), law, negative) != 1):
                raise ArithmeticError("inverse-link or backtrack identity failed")
        for face in faces:
            loop = face + (face[0],)
            expected = -1 if negative else 1
            if (_loop_sign(loop, law, negative) != expected
                    or _loop_sign(tuple(reversed(loop)), law, negative) != expected):
                raise ArithmeticError("candidate face phase or reversal changed")
            gauge = {port: (-1 if port % 2 else 1) for port in face}
            transformed = 1
            for a, b in zip(loop[:-1], loop[1:], strict=True):
                transformed *= gauge[a] * _link_sign(a, b, law, negative) * gauge[b]
            if transformed != expected:
                raise ArithmeticError("candidate face phase failed a port rephasing")

        # The equal-length route comparison 0-1-5 versus 0-11-5 closes
        # around two adjacent triangles.  It misses the uniform pi phase on
        # each individual triangle because the four link signs multiply to +1.
        double_face_loop = (0, 1, 5, 11, 0)
        if any(tuple(sorted((a, b))) not in seams for a, b in zip(
                double_face_loop[:-1], double_face_loop[1:], strict=True)):
            raise ArithmeticError("declared two-face loop left the native seam graph")
        if _loop_sign(double_face_loop, law, negative) != 1:
            raise ArithmeticError("uniform moving-branch phase altered the even loop")

    return {
        "schema": "oph.native_link_selection.v1",
        "status": "COMMITTED_NATIVE_LINK_PHASE_UNSELECTED",
        "native_states": len(states),
        "native_seams": len(seams),
        "oriented_triangular_faces": len(faces),
        "complete_committed_matrix_units_checked": matrix_units,
        "one_token_moving_columns_checked": moving_columns,
        "moving_current_branch_amplitude_squared": "1/2",
        "moving_branch_normalized_operator": "|target><source|; rank one, not a two-port unitary",
        "inverse_branch_normalized_product": "|source><source|, not identity on the one-token two-port space",
        "same_complete_C_E_committed_instruments_under_both_branch_phase_completions": True,
        "all_plus_completion_face_products": "+1 on all 20 faces",
        "all_minus_moving_branch_completion_face_products": "-1 on all 20 faces",
        "two_face_0_1_5_11_0_loop_product_in_both_completions": "+1",
        "two_face_even_loop_detects_uniform_pi_per_triangle": False,
        "both_completions_inverse_reversal_backtrack_and_port_rephasing_checks": True,
        "all_minus_phase_rule_depends_on_port_labels": False,
        "source_selected_or_observed_native_link_holonomy": False,
        "scope": "A coherent branch-phase interface can distinguish these completions, "
                 "but committed native C/E repair and its current records do not select one. "
                 "The full 78-state matrix-unit test and the one-token number-sector "
                 "extension are distinct domains.",
        "all_pass": True,
    }


def canonical(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main():
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
        raise SystemExit("Native link-selection receipt mismatch")
    print(json.dumps({"all_pass": True, "sha256": hashlib.sha256(encoded.encode()).hexdigest()},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

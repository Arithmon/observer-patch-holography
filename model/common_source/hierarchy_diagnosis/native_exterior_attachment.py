"""Exact same-carrier exterior-matter attachment test on native two-load repair.

The native 66-state binary source is the exterior square of its twelve port
coefficient space as a vector space.  This tests the source's rank-one and
rank-five spatial bands and its existing protected core as a candidate
carrier for the independently decoded fifteen-dimensional matter module.
Neither an exterior action on port coefficients nor its charge split is
assumed to be the physical gauge current.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path

import sympy as s

from .. import slow_geometry
from ...native_equilibration import cut_completion
from ...native_modular_reference import common_caps, local_preparation
from ...observer import native
from ...repair_scheduler import coins


ROOT = Path(__file__).resolve().parents[3]
MATTER = ROOT / "model/observer_particles/results/representation_identification.json"
STATES = native.states(binary=True)
PAIRS = tuple(tuple(i for i, value in enumerate(state) if value) for state in STATES)
INDEX = {state: i for i, state in enumerate(STATES)}


def _rational(value):
    a, b = value
    if b != 0:
        raise ArithmeticError("the selected native band is not rational")
    return s.Rational(a)


def _band(index):
    c = slow_geometry.source()
    rows = c.band_projector_exact(index)
    return s.ImmutableMatrix(12, 12, lambda i, j: _rational(rows[i][j]))


def _exterior(projector):
    """Second exterior power in the source's actual binary-pair ordering."""
    return s.ImmutableMatrix(66, 66, lambda a, b:
        projector[PAIRS[a][0], PAIRS[b][0]] * projector[PAIRS[a][1], PAIRS[b][1]]
        - projector[PAIRS[a][0], PAIRS[b][1]] * projector[PAIRS[a][1], PAIRS[b][0]])


def _candidate():
    p1, p5 = _band(0), _band(2)
    if (p1*p1 != p1 or p5*p5 != p5 or p1*p5 != s.zeros(12)
            or s.trace(p1) != 1 or s.trace(p5) != 5):
        raise ArithmeticError("native constant/five-dimensional bands changed")
    core = tuple(local_preparation.CORE)
    if len(core) != 2:
        raise ArithmeticError("a two-port protected core is required")
    columns = p5[:, list(core)]
    gram = columns.T*columns
    if gram.det() != s.Rational(1, 6):
        raise ArithmeticError("the projected core no longer spans a two-plane")
    p2 = s.ImmutableMatrix(columns*gram.inv()*columns.T)
    p3 = s.ImmutableMatrix(p5-p2)
    if (p2*p2 != p2 or p3*p3 != p3 or p2*p3 != s.zeros(12)
            or s.trace(p2) != 2 or s.trace(p3) != 3):
        raise ArithmeticError("the core-defined three/two split failed")
    q15, q10 = _exterior(p1+p5), _exterior(p5)
    q5 = s.ImmutableMatrix(q15-q10)
    if (q15*q15 != q15 or q10*q10 != q10 or q5*q5 != q5
            or q5*q10 != s.zeros(66)
            or (s.trace(q15), s.trace(q10), s.trace(q5)) != (15, 10, 5)):
        raise ArithmeticError("native exterior 10+5 carrier failed")
    return p1, p2, p3, p5, q15, q10, q5, gram


def _accepted_transfer(weights):
    matrix = s.zeros(66)
    for column, state in enumerate(STATES):
        for seam, weight in weights.items():
            for coin in (0, 1):
                final, _ = native.native(state, seam, coin)
                matrix[INDEX[final], column] += s.Rational(weight)/2
    if matrix != matrix.T or len(STATES) != 66:
        raise ArithmeticError("native two-load transfer lost its symmetric source law")
    return matrix


def _one_body_exterior_action(operator):
    """Declared signed CAR lift of a port operator, as a comparison control."""
    return s.ImmutableMatrix(66, 66, lambda a, b:
        operator[PAIRS[a][0], PAIRS[b][0]] * int(PAIRS[a][1] == PAIRS[b][1])
        - operator[PAIRS[a][1], PAIRS[b][0]] * int(PAIRS[a][0] == PAIRS[b][1])
        + operator[PAIRS[a][1], PAIRS[b][1]] * int(PAIRS[a][0] == PAIRS[b][0])
        - operator[PAIRS[a][0], PAIRS[b][1]] * int(PAIRS[a][1] == PAIRS[b][0]))


def _commutator_certificate(transfer, projector):
    commutator = transfer*projector-projector*transfer
    return {"rank_exact": commutator.rank(),
            "frobenius_norm_squared_exact": str(s.trace(commutator.T*commutator)),
            "entry_0_1_exact": str(commutator[0, 1])}


def _actual_quantum_leakage(law, projector, weights):
    """Unconditional accepted P-outside probability from rho=P/rank(P)."""
    result = s.Rational(0)
    rank = s.trace(projector)
    for seam, weight in weights.items():
        for _, _, mapping in coins.coin_maps(law, seam):
            inputs = tuple(source for source, _ in mapping)
            outputs = tuple(target for _, target in mapping)
            remaining = sum((projector[i, i] for i in inputs), s.Rational(0))
            overlap = sum((projector[outputs[a], outputs[b]]
                           *projector[inputs[a], inputs[b]]
                           for a in range(len(mapping)) for b in range(len(mapping))),
                          s.Rational(0))
            result += s.Rational(weight)*(remaining-overlap)/(2*rank)
    return result


def invariant_rank_certificate():
    """Apply the existing connected-cap fixed algebra to the finite selector."""
    core, ideal, schedule = common_caps.cap_schedule()
    if core != tuple(local_preparation.CORE) or schedule["gap_lower_per_attempt"] != "1/600":
        raise ArithmeticError("the common native cap theorem changed")
    selected = local_preparation.selectors(1)
    weights = selected["completion_accepted_seam_weights"]
    timeout = selected["completion_timeout_probability"]
    if (set(weights) != set(ideal) or min(weights.values()) <= 0
            or sum(weights.values(), F())+timeout != 1):
        raise ArithmeticError("finite selector lacks the connected-cap support")
    minorization = min(F(weights[edge])/ideal[edge] for edge in ideal)
    if minorization != F(3, 4) or any(weights[edge] < minorization*ideal[edge] for edge in ideal):
        raise ArithmeticError("finite selector no longer minorizes the ideal cap law")
    # Each branch K=(1/sqrt2)M is a partial injection. Summing all current
    # labels and both coins counts every input and output twice, for both
    # native laws. The accepted mixture plus timeout identity is bistochastic.
    for law in ("C", "E"):
        for seam in weights:
            inputs, outputs = [0]*66, [0]*66
            for _, _, mapping in coins.coin_maps(law, seam):
                for source, target in mapping:
                    inputs[source] += 1
                    outputs[target] += 1
            if inputs != [2]*66 or outputs != [2]*66:
                raise ArithmeticError("native source Kraus sum is not bistochastic")
    sectors = cut_completion.cut_sectors(core)
    sizes = tuple((len(sector.cap_states), len(sector.complement_states)) for sector in sectors)
    if sizes != ((1, 45), (2, 10), (1, 1)):
        raise ArithmeticError("protected-core fixed algebra blocks changed")
    possible = {0}
    for cap_dimension, complement_dimension in sizes:
        possible = {rank+complement_dimension*block_rank
                    for rank in possible for block_rank in range(cap_dimension+1)}
    ranks = sorted(possible)
    if 15 in possible or ranks != [0, 1, 10, 11, 20, 21, 45, 46, 55, 56, 65, 66]:
        raise ArithmeticError("a rank-fifteen fixed projector appeared")
    return {
        "finite_selector_same_positive_twenty_seam_support": True,
        "selector_timeout_weight_exact": str(timeout),
        "ideal_common_cap_gap_lower_exact": schedule["gap_lower_per_attempt"],
        "finite_selector_ideal_mix_minorization_exact": str(minorization),
        "finite_selector_gap_lower_exact": str(minorization*F(1, 600)),
        "C_and_E_complete_work_channels_bistochastic_exact": True,
        "existing_fixed_algebra_blocks": "C*I45 direct_sum (M2 tensor I10) direct_sum C*I1",
        "all_fixed_projector_ranks_exact": ranks,
        "rank_fifteen_noiseless_invariant_work_subspace_possible": False,
        "proof_scope": "The unread finite C/E mixed-cap channel on the fixed binary 66-state work carrier. Its positive weighted seams have the same fixed algebra as the established connected ideal cap law; timeout adds identity. A projector invariant under every work Kraus branch is fixed because the complete channel is unital and trace preserving. No conclusion about record-dependent recovery, ancillas, noiseless subsystems, operator responses, or another coherent law.",
    }


def report():
    """Same-carrier possibility and two exact failures of its direct reading."""
    p1, p2, p3, p5, q15, q10, q5, gram = _candidate()
    n3, n2 = int(s.trace(p3)), int(s.trace(p2))
    actual = Counter({-4: n3*(n3-1)//2, 1: n3*n2, 6: n2*(n2-1)//2,
                      -2: n3, 3: n2})
    decoded = json.loads(MATTER.read_text())["decoded_even_internal_multiplets"]
    expected = Counter({int(row["q6Y"]): int(row["dimension"]) for row in decoded})
    if expected != Counter({-4: 3, 1: 6, 6: 1, 2: 3, -3: 2}):
        raise ArithmeticError("decoded fifteen-component comparison changed")
    line_shifts = (2-(-2), -3-3)
    if actual == expected or line_shifts[0] == line_shifts[1]:
        raise ArithmeticError("natural one-body exterior charges unexpectedly match")
    selector = local_preparation.selectors(1)
    selected = _accepted_transfer(selector["completion_accepted_seam_weights"])
    timeout = s.Rational(selector["completion_timeout_probability"])
    complete = selected+timeout*s.eye(66)
    uniform = _accepted_transfer({seam: F(1, 30) for seam in slow_geometry.source().seams()})
    chosen = _commutator_certificate(selected, q15)
    if timeout != s.Rational(5, 32) or _commutator_certificate(complete, q15) != chosen:
        raise ArithmeticError("complete selector timeout changed the exterior leakage claim")
    quantum_leakage = {law: _actual_quantum_leakage(
        law, q15, selector["completion_accepted_seam_weights"]) for law in ("C", "E")}
    if quantum_leakage != {"C": s.Rational(379, 2560), "E": s.Rational(703, 2560)}:
        raise ArithmeticError("actual native quantum branch leakage certificate changed")
    control = _commutator_certificate(uniform, q15)
    signed = _one_body_exterior_action(s.Matrix(slow_geometry.source().laplacian()))
    if (signed*q5 != 6*q5 or signed*q10 != 12*q10
            or signed*q15 != q15*signed):
        raise ArithmeticError("signed exterior comparison action failed its band law")
    if (sum(selected[:, 0]) != s.Rational(27, 32)
            or sum(uniform[:, 0]) != 1
            or chosen != {"rank_exact": 30, "frobenius_norm_squared_exact": "961/8192",
                          "entry_0_1_exact": "1/64"}
            or control != {"rank_exact": 30, "frobenius_norm_squared_exact": "29/600",
                           "entry_0_1_exact": "1/120"}):
        raise ArithmeticError("direct native exterior dynamics certificate changed")
    return {
        "native_binary_two_load_dimension": 66,
        "source_bands": {"constant_rank": int(s.trace(p1)),
                         "antipodal_even_rank": int(s.trace(p5)),
                         "projected_protected_core_gram_determinant_exact": str(gram.det()),
                         "core_plane_rank": int(s.trace(p2)),
                         "orthogonal_complement_rank": int(s.trace(p3))},
        "same_carrier_exterior_candidate": {"projector_rank": int(s.trace(q15)),
                                             "wedge_five_rank": int(s.trace(q10)),
                                             "constant_wedge_five_rank": int(s.trace(q5)),
                                             "all_projector_identities_exact": True},
        "conditional_additive_one_body_charge": {
            "q6Y_on_core_defined_three_plane": -2,
            "q6Y_on_core_defined_two_plane": 3,
            "native_exterior_charge_multiplicities": {str(k): v for k, v in sorted(actual.items())},
            "decoded_M15_charge_multiplicities": {str(k): v for k, v in sorted(expected.items())},
            "ten_component_wedge_charges_match": True,
            "five_component_wedge_charges_match": False,
            "single_constant_line_charge_shift_needed_for_three_and_two": list(line_shifts),
            "one_additive_one_body_charge_action_recovers_M15": False},
        "direct_native_transfer_attachment": {
            "mixed_cap_one_round_accepted_weight_exact": "27/32",
            "selector_timeout_weight_exact": str(timeout),
            "complete_map_with_timeout_has_same_projector_commutator": True,
            "mixed_cap_projector_commutator": chosen,
            "complete_C_E_accepted_quantum_leakage_from_Q15_over_15_exact": {
                law: str(value) for law, value in quantum_leakage.items()},
            "uniform_all_seam_control_projector_commutator": control,
            "fifteen_dimensional_candidate_invariant_under_either_transfer": False},
        "fixed_work_carrier_invariant_rank_test": invariant_rank_certificate(),
        "declared_signed_exterior_control": {
            "one_body_operator": "the same native thirty-seam port Laplacian",
            "constant_wedge_five_eigenvalue_exact": 6,
            "wedge_five_eigenvalue_exact": 12,
            "fifteen_dimensional_carrier_invariant_exact": True,
            "this_CAR_lift_selected_by_recorded_native_repair": False},
        "interpretation": "A specified exterior vector-space map on the actual 66-state source gives a 15-dimensional 10+5 carrier from its rank-one/five bands and protected core. Exterior orientation signs are an additional representation choice; the native unordered-pair repair does not inherit them. The additive one-body charge gives 10+5 rather than the decoded 10+bar5. The complete finite C/E channel has no invariant rank-15 work subspace. This excludes the stated direct work-carrier attachment, not encoded, record-corrected, interacting or continuum field maps.",
    }

#!/usr/bin/env python3
"""Coordinate-corrected constructive SL(2,F5) -> PORT-SPIN-LIFT bridge.

The #314 current fixture and the committed Lean port table use isomorphic but
not identical vertex labellings.  This revision derives an exact
orientation-preserving carrier relabelling before comparing the two sixty-row
actions.  No row equality is assumed across coordinate systems.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import sl2f5_port_spin_bridge_certificate as base

c784 = base.c784
m314 = base.m314
MODULE_DIR = base.MODULE_DIR
SCHEMA = base.SCHEMA
Matrix = base.Matrix
SpinMatrix = base.SpinMatrix
PortRow = base.PortRow
require = base.require


def inverse_permutation(p: Sequence[int]) -> PortRow:
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def relabel_port_row(row: PortRow, phi: Sequence[int]) -> PortRow:
    """Conjugate a permutation from carrier labels to committed labels.

    ``phi`` maps old carrier vertex labels to committed port labels, so the
    transported row is ``phi ∘ row ∘ phi^{-1}``.
    """

    inv = inverse_permutation(phi)
    return tuple(phi[row[inv[j]]] for j in range(len(row)))


def current_carrier_relabel(
    algebra: m314.CurrentAlgebra,
    fixture: c784.Fixture,
) -> tuple[PortRow, dict[str, int]]:
    """Derive an orientation-preserving relabelling of the #314 carrier.

    The current fixture's carrier has its own p00..p11 coordinate choice.
    Compare only incidence and orientation, then transport every current-row
    permutation into the committed Lean coordinates.
    """

    isos = c784.graph_isomorphisms(algebra.carrier.adjacency, fixture.adjacency)
    require(bool(isos), "SPIN_CARRIER_RELABEL", "#314 carrier is not graph-isomorphic to the committed carrier")

    committed = c784.relabel_faces(fixture.oriented_faces, tuple(range(12)))
    matching = [
        tuple(phi)
        for phi in isos
        if c784.relabel_faces(algebra.carrier.faces, phi) == committed
    ]
    require(
        bool(matching),
        "SPIN_CARRIER_ORIENTATION",
        "no carrier relabelling preserves the committed face orientation",
    )

    phi = min(matching)
    transported = {relabel_port_row(tuple(row), phi) for row in algebra.plus}
    require(
        transported == set(fixture.port_action),
        "SPIN_PORT_ACTION",
        "relabelled PORT-SPIN-LIFT rows differ from the committed port action",
    )
    return phi, {
        "graph_isomorphisms": len(isos),
        "orientation_matching_relabellings": len(matching),
    }


def build_certificate() -> dict[str, Any]:
    # --- Canonical SL(2,F5) source and its exact committed-port quotient -----
    coset_manifest = c784.load_json(MODULE_DIR / "manifests" / "coset_carrier_reference.json")
    source = c784.build_source(coset_manifest)
    quotient, quotient_summary = c784.central_quotient(source)
    fixture = c784.load_fixture()

    solutions = base.presentation_solutions(source)
    require(len(solutions) == 120, "PRESENTATION", f"expected 120 presentation solutions, got {len(solutions)}")
    source_generators = min(solutions)
    placement = base.placement_of_solution(quotient, source_generators)
    incidence = c784.build_incidence(quotient, placement)
    require(c784.classify(incidence) == "CARRIER", "PRESENTATION", "chosen presentation solution is not a carrier placement")
    phi, relabel_summary = c784.relabel(quotient, placement, incidence, fixture)
    source_to_port = base.port_action_for_placement(quotient, placement, phi)

    source_port_rows = {source_to_port(g) for g in source.elements}
    require(source_port_rows == set(fixture.port_action), "SOURCE_PORT_ACTION", "SL2F5 quotient action differs from the committed port rows")
    source_kernel = [g for g in source.elements if source_to_port(g) == tuple(range(12))]
    require(
        set(source_kernel) == {source.identity, source.minus_identity},
        "SOURCE_PORT_KERNEL",
        "source-to-port kernel is not exactly {+I,-I}",
    )

    # --- Recompute PORT-SPIN-LIFT, then transport its carrier coordinates ----
    matter_manifest = m314.load_json(MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json")
    upstream = m314.load_upstream(matter_manifest, MODULE_DIR)
    algebra = m314.CurrentAlgebra(upstream["current_manifest"], MODULE_DIR)
    spin = m314.spin_lift_certificate(algebra)
    raw_lifts: Mapping[PortRow, SpinMatrix] = spin["lifts"]
    require(
        set(raw_lifts) == {tuple(row) for row in algebra.plus},
        "SPIN_INTERNAL_ACTION",
        "PORT-SPIN-LIFT keys differ from the current fixture's own rotation rows",
    )

    current_phi, current_relabel_summary = current_carrier_relabel(algebra, fixture)
    lifts: dict[PortRow, SpinMatrix] = {}
    for raw_row, lift in raw_lifts.items():
        committed_row = relabel_port_row(tuple(raw_row), current_phi)
        require(
            committed_row not in lifts,
            "SPIN_CARRIER_RELABEL",
            "carrier relabelling identifies two distinct rotation rows",
        )
        lifts[committed_row] = lift
    require(
        set(lifts) == set(fixture.port_action),
        "SPIN_PORT_ACTION",
        "transported PORT-SPIN-LIFT rows differ from the committed port action",
    )
    spin_elements, spin_to_port = base.build_spin_group_index(lifts)

    # --- Match one exact presentation triple through the eight lift signs ----
    a, b, c = source_generators
    generator_rows = (source_to_port(a), source_to_port(b), source_to_port(c))
    require(
        tuple(c784.element_order(source.mul, source.identity, x) for x in source_generators) == (4, 6, 10),
        "PRESENTATION_ORDERS",
        "chosen SL2F5 presentation generators do not have orders (4,6,10)",
    )
    base_lifts = tuple(lifts[row] for row in generator_rows)
    choices = base.signed_generator_choices(base_lifts)
    require(
        len(choices) == 1,
        "GENERATOR_SIGN_CHOICE",
        f"expected one exact sign choice for the spin presentation, got {len(choices)}",
    )
    signs, spin_generators = choices[0]

    # --- Extend by words, then verify all 120^2 products ----------------------
    source_to_spin = base.extend_generator_map(source, source_generators, spin_generators)
    image_keys = {m314.matrix_key(matrix) for matrix in source_to_spin.values()}
    require(
        image_keys == set(spin_elements),
        "ISOMORPHISM_SURJECTIVE",
        "generator map does not hit exactly the 120 PORT-SPIN-LIFT elements",
    )

    product_checks = 0
    port_homomorphism_checks = 0
    for g in source.elements:
        for h in source.elements:
            gh = source.mul(g, h)
            require(
                base.spin_equal(source_to_spin[gh], m314.cmul(source_to_spin[g], source_to_spin[h])),
                "ISOMORPHISM_HOMOMORPHISM",
                "SL2F5-to-spin map fails a multiplication check",
            )
            product_checks += 1
            require(
                source_to_port(gh)
                == base.permutation_compose(source_to_port(g), source_to_port(h)),
                "SOURCE_PORT_HOMOMORPHISM",
                "SL2F5-to-port action fails a multiplication check",
            )
            port_homomorphism_checks += 1

    identity2 = m314.cidentity(2)
    minus_identity2 = base.spin_neg(identity2)
    require(base.spin_equal(source_to_spin[source.identity], identity2), "CENTER_MAP", "+I does not map to +I2")
    require(base.spin_equal(source_to_spin[source.minus_identity], minus_identity2), "CENTER_MAP", "-I does not map to -I2")

    # --- Verify the cover square element by element ---------------------------
    square_checks = 0
    for g in source.elements:
        spin_key = m314.matrix_key(source_to_spin[g])
        require(spin_key in spin_to_port, "SPIN_PORT_MAP", "isomorphism image is not a registered spin lift")
        require(
            spin_to_port[spin_key] == source_to_port(g),
            "COVER_SQUARE",
            "the SL2F5-to-spin isomorphism does not commute with the port quotient",
        )
        square_checks += 1

    spin_kernel_keys = {key for key, row in spin_to_port.items() if row == tuple(range(12))}
    require(
        spin_kernel_keys == {m314.matrix_key(identity2), m314.matrix_key(minus_identity2)},
        "SPIN_PORT_KERNEL",
        "spin-to-port kernel is not exactly {+I2,-I2}",
    )

    return {
        "schema": SCHEMA,
        "construction": "explicit_generator_correspondence_with_derived_carrier_relabelling_and_exhaustive_group_check",
        "source": {
            "group": "SL(2,F5)",
            "order": len(source.elements),
            "center": ["+I", "-I"],
            "central_quotient_order": quotient_summary["order"],
            "presentation": "a^2 = b^3 = c^5 = abc = -I",
            "presentation_solution_count": len(solutions),
            "chosen_generators": [base.matrix_as_list(x) for x in source_generators],
            "chosen_generator_orders": [4, 6, 10],
        },
        "carrier_bridge": {
            "committed_port_rows": len(source_port_rows),
            "source_to_port_kernel": ["+I", "-I"],
            "coset_to_committed_relabel": relabel_summary,
            "current_to_committed_relabel": {
                **current_relabel_summary,
                "map": list(current_phi),
            },
            "chosen_generator_port_rows": [base.row_as_list(row) for row in generator_rows],
            "source_port_homomorphism_checks": port_homomorphism_checks,
        },
        "port_spin_lift": {
            "order": len(spin_elements),
            "center": ["+I2", "-I2"],
            "base_generator_lift_signs": list(signs),
            "order_profile": spin["order_profile"],
            "exact_field": "Q(sqrt(5), i)",
            "faithful_two_dimensional_complex_representation": True,
        },
        "isomorphism": {
            "source_elements_mapped": len(source_to_spin),
            "distinct_spin_images": len(image_keys),
            "multiplication_checks": product_checks,
            "center_map": "+I -> +I2; -I -> -I2",
            "bijective": True,
        },
        "commuting_cover_square": {
            "checks": square_checks,
            "source_kernel": ["+I", "-I"],
            "spin_kernel": ["+I2", "-I2"],
            "commutes_on_all_source_elements": True,
        },
        "claim_boundary": {
            "proves": (
                "an explicit exact group isomorphism from the canonical SL(2,F5) source "
                "to the executable PORT-SPIN-LIFT, after deriving the coordinate change "
                "between the #314 carrier labels and the committed Lean port labels"
            ),
            "method_not_used": [
                "same-order inference",
                "element-order-profile classification",
                "uniqueness of the non-split central extension",
            ],
            "does_not_prove": [
                "McKay or E8 correspondence",
                "selection or derivation of phi",
                "mass relation",
                "physical identification of the finite spin action",
                "source selection of the #566 charged-double-triplet current fixture",
            ],
        },
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

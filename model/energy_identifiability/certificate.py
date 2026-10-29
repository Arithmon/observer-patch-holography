"""Dynamical energy-offset comparisons require the declared observable access.

This package supplies a finite pre-fit diagnostic, not an energy calibration.
Native instruments and field Hamiltonians retain their distinct carriers.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as F
import json


def comparison_graph(vertices, edges):
    """Exact rank/nullity of the declared matrix-unit comparison graph.

An edge i--j means access to the off-diagonal matrix units E_ij and E_ji,
including their generator response, not merely a population transition.
This function validates a DECLARATION; it cannot certify source access.
"""
    if type(vertices) is not int or vertices < 1:
        raise ValueError("a positive integer vertex count required")
    adjacency = [set() for _ in range(vertices)]
    for edge in edges:
        if (not isinstance(edge, (tuple, list)) or len(edge) != 2
                or any(type(i) is not int or not 0 <= i < vertices for i in edge)
                or edge[0] == edge[1]):
            raise ValueError("two distinct admitted vertex indices required")
        i, j = edge
        adjacency[i].add(j)
        adjacency[j].add(i)
    components, seen = [], set()
    for root in range(vertices):
        if root in seen:
            continue
        queue, component = deque((root,)), []
        seen.add(root)
        while queue:
            i = queue.popleft()
            component.append(i)
            for j in sorted(adjacency[i]):
                if j not in seen:
                    seen.add(j)
                    queue.append(j)
        components.append(sorted(component))
    return {"vertices": vertices, "edges": sum(map(len, adjacency))//2,
            "components": components, "incidence_rank": vertices-len(components),
            "offset_gauge_dimension": len(components),
            "relative_offsets_unidentified": len(components)-1}


def reconstruct_offsets(vertices, measurements):
    """Recover a_j-a_i from exact oriented differences, or reject a bad cycle.

One zero is chosen per connected component. Cross-component differences
remain unidentified; the returned zeros must not be interpreted as equal
physical energy origins. Rational inputs make cycle checks exact.
"""
    rows = tuple(measurements)
    if any(not isinstance(row, (tuple, list)) or len(row) != 3
           or type(row[2]) not in (int, F) for row in rows):
        raise ValueError("exact rational oriented differences required")
    graph = comparison_graph(vertices, [(i, j) for i, j, _ in rows])
    adjacency = [[] for _ in range(vertices)]
    for i, j, value in rows:
        adjacency[i].append((j, F(value)))
        adjacency[j].append((i, -F(value)))
    offset = [None]*vertices
    for component in graph["components"]:
        root = component[0]
        offset[root] = F()
        queue = deque((root,))
        while queue:
            i = queue.popleft()
            for j, difference in adjacency[i]:
                proposed = offset[i]+difference
                if offset[j] is None:
                    offset[j] = proposed
                    queue.append(j)
                elif offset[j] != proposed:
                    raise ValueError("inconsistent energy-offset cycle")
    return {**graph, "offsets_with_one_zero_per_component": tuple(offset)}


def native_sector_certificate():
    """Check QK=KQ on every native raw source/core matrix entry.

The core is classical record-driven metadata. Q counts occupations in
that transported core. This does not assign a new intrinsic source charge.
"""
    from ..completion import observer
    from ..exchange_sector import sector

    binary = set(sector.BASIS)
    checks, slots = Counter(), Counter()
    for law in ("C", "E"):
        for seam in observer.SEAMS:
            for _, _, current, pairs in observer._raw_maps(law, seam):
                slots[law] += 1
                for source, target in pairs:
                    x, y = observer.STATES[source], observer.STATES[target]
                    if x not in binary:
                        continue
                    if y not in binary:
                        raise ArithmeticError("native branch left the binary shell")
                    for core in sector.CORES:
                        a, b = seam
                        moved = (tuple(sorted(b if p == a else a if p == b else p
                                              for p in core)) if current else core)
                        if sum(x[p] for p in core) != sum(y[p] for p in moved):
                            raise ArithmeticError("native branch mixes charge sectors")
                        checks[law] += 1
    return {"source_core_dimension": len(binary)*len(sector.CORES),
            "raw_record_slots": dict(slots), "matrix_entry_checks": dict(checks),
            "all_raw_branch_maps_intertwine_charge_projectors": True,
            "arbitrary_scalar_sector_offsets_invisible_under_declared_protocol": True,
            "block_diagonal_initial_state_required_for_equal_full_CQ_states": True,
            "same_initial_preparation_compared_for_both_generators": True,
            "every_intervention_must_preserve_charge_projectors": True,
            "separate_clock_factor_removes_offset_gauge": False,
            "all_source_hardware_interventions_are_sector_preserving": False}


def native_archive_certificate():
    """Induction certificate for zero energy excess on the entire native tape.

The reference at length n is Phi_n(pi), with the SAME complete instrument
and retained records as Phi_n(rho). It is not a stationary growing archive.
Every partial injection preserves equality of diagonals; -log Phi_n(pi)
is diagonal on its support, so the modular-energy difference is exactly 0.
"""
    from ..local_repair_family import local_bias as local
    from ..native_localization import localization
    from ..observer import adaptive_record_recovery as adaptive

    rho = localization.initial_density()
    pi = local.reference(localization.parameters())
    if any(rho[i][i] != value or value <= 0 for i, value in enumerate(pi)):
        raise ArithmeticError("localized coherence changed the reference diagonal")
    slots, nonzero, entries = 0, 0, 0
    for branch in adaptive.branch_maps(binary=True):
        slots += 1
        inputs, outputs = branch.inputs, branch.outputs
        if (len(inputs) != len(outputs) or len(set(inputs)) != len(inputs)
                or len(set(outputs)) != len(outputs)
                or any(type(i) is not int or not 0 <= i < len(pi)
                       for i in (*inputs, *outputs))):
            raise ArithmeticError("complete raw branch is not a partial injection")
        nonzero += bool(inputs)
        entries += len(inputs)
        # A target diagonal receives exactly one input diagonal. Positive,
        # source-independent coin weights multiply both states identically.
        if sum((rho[i][i]-pi[i] for i in inputs), F()) != 0:
            raise ArithmeticError("a raw record probability changed")
    return {"complete_raw_slots": slots, "nonzero_slots": nonzero,
            "partial_injection_entries": entries,
            "initial_off_diagonal_difference_nonzero": any(
                rho[i][j] for i in range(len(pi)) for j in range(i)),
            "all_lengths_source_and_complete_archive_modular_energy_excess": "0",
            "reference": "same-history pushforward Phi_n(pi), logarithm on its support",
            "equal_record_probabilities_and_branch_diagonals": True,
            "reference_on_growing_archive_is_stationary": False,
            "zero_energy_implies_zero_relative_entropy": False,
            "additional_coherent_analyzers_are_included": False}


def clock_history_control():
    """Exact rational three-step CQ control with a nontrivial clock factor.

The q-preserving rotations may depend on q. Hence the clock can already
carry sector information, but a scalar q-energy shift is still invisible.
The Ramsey control deliberately permits a cross-sector preparation/read.
"""
    import sympy as sp

    f = sp.Rational
    rotation = sp.Matrix([[f(3, 5), -f(4, 5)], [f(4, 5), f(3, 5)]])
    base = sp.diag(rotation, rotation.T)
    shift = sp.diag(1, 1, sp.I, sp.I)
    rho = sp.diag(sp.Matrix([[f(1, 4), f(1, 12)], [f(1, 12), f(1, 4)]]),
                  sp.Matrix([[f(1, 3), f(1, 12)], [f(1, 12), f(1, 6)]]))
    operators = (sp.diag(f(3, 5), f(4, 5), f(3, 5), f(4, 5)),
                 sp.diag(f(4, 5), f(3, 5), f(4, 5), f(3, 5)))

    def histories(unitary):
        states = [rho]
        for _ in range(3):
            states = [k*unitary*s*unitary.H*k.H for s in states for k in operators]
        return states

    first, second = histories(base), histories(shift*base)
    if any(a != b for a, b in zip(first, second, strict=True)):
        raise ArithmeticError("clock factor exposed a conserved-sector scalar offset")
    if sum(sp.trace(a) for a in first) != 1:
        raise ArithmeticError("clock history control dropped an outcome")
    plus = sp.Matrix([[f(1, 2), f(1, 2)], [f(1, 2), f(1, 2)]])
    relative_phase = sp.diag(1, sp.I)
    before = sp.trace(plus*plus)
    after = sp.simplify(sp.trace(plus*relative_phase*plus*relative_phase.H))
    if (before, after) != (1, f(1, 2)):
        raise ArithmeticError("cross-sector coherence countercontrol failed")
    return {"retained_histories": len(first), "complete_CQ_states_exactly_equal": True,
            "clock_rotations_depend_on_sector": True,
            "Ramsey_plus_probability_without_offset": str(before),
            "Ramsey_plus_probability_with_pi_over_2_phase": str(after),
            "Ramsey_preparation_and_read_mix_sectors": True,
            "physical_time_or_mass_calibrated": False}


def field_reference_control():
    """All 64 actual field fibers; exact graph rank, floating spectral control.

The 63 star edges are mathematical matrix units. This function does not
construct their preparation/readout using the local source gate grammar.
"""
    import numpy as np
    from ..band_clock.clock import reference_spectrum

    reference = reference_spectrum()
    offsets = reference["fiber_offsets"]
    count = len(offsets)-1
    graph = comparison_graph(count, [(0, i) for i in range(1, count)])
    log_z = np.log(reference["fiber_partitions"])
    energy = reference["energy"]
    beta = reference["beta"]
    independent = np.concatenate([
        beta*energy[a:b]+log_z[i]
        for i, (a, b) in enumerate(zip(offsets, offsets[1:]))])
    global_k = -reference["log_probability"]
    within = max(float(np.max(abs((independent[a:b]-global_k[a:b])
                                  -(independent[a]-global_k[a]))))
                 for a, b in zip(offsets, offsets[1:]))
    # Edge (0,i) determines log Z_i-log Z_0 after the common beta H
    # response is subtracted. All relative fiber constants are then fixed.
    observed = np.array([independent[a]-global_k[a] for a in offsets[:-1]])
    tree_differences = observed[1:]-observed[0]
    residual = float(np.max(abs(tree_differences-(log_z[1:]-log_z[0]))))
    if within > 1e-12 or residual > 1e-12 or graph["incidence_rank"] != 63:
        raise ArithmeticError("actual field offset comparison failed")
    spread = float(np.ptp(log_z))
    if spread < .01:
        raise ArithmeticError("independent fiber normalization countercontrol vanished")
    return {"field_dimension": len(energy), "momentum_fibers": count,
            "fiber_dimensions": [b-a for a, b in zip(offsets, offsets[1:])],
            "intrafiber_response_residual_float": within,
            "partition_log_spread_float": spread,
            "tree_reconstruction_residual_float": residual,
            "additional_comparison_rank_for_one_global_offset": graph["incidence_rank"],
            "relative_fiber_offsets_unidentified_without_cross_comparisons": count-1,
            "tree_matrix_units_source_accessibility_certified": False,
            "spectral_values_are_float_diagnostics": True,
            "fixed_particle_number": 2,
            "vacuum_sector_is_in_this_carrier": False,
            "additive_c_times_N_invisible_on_fixed_N2": True,
            "common_reference_or_physical_energy_selected": False}


def report(include_field=False):
    result = {"scope": "finite pre-fit energy-offset identifiability",
              "native_sector_gauge": native_sector_certificate(),
              "native_complete_archive_energy": native_archive_certificate(),
              "clock_and_Ramsey_controls": clock_history_control()}
    if include_field:
        result["actual_field_reference"] = field_reference_control()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", action="store_true", help="also diagonalize all 64 field fibers")
    args = parser.parse_args()
    print(json.dumps(report(args.field), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

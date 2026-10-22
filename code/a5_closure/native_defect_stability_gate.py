#!/usr/bin/env python3
"""Exact native-record and declared Z6 defect gate on the icosahedral carrier.

The coefficient group, puncture pair, and link updates below are diagnostic
choices. Only the port loads, seam currents, and carrier incidence come from
the executable observer. No link field or puncture is added to native repair.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import flux_defect_criterion_certificate as flux  # noqa: E402
from model.observer import native, native_loop_readout  # noqa: E402


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    """Gaussian elimination over F2 or F3, without floating point."""
    if prime not in (2, 3) or not matrix or not matrix[0]:
        raise ValueError("a nonempty matrix over F2 or F3 is required")
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix rows have unequal lengths")
    rows = [[value % prime for value in row] for row in matrix]
    pivot = 0
    for column in range(width):
        target = next((r for r in range(pivot, len(rows)) if rows[r][column]), None)
        if target is None:
            continue
        rows[pivot], rows[target] = rows[target], rows[pivot]
        inverse = pow(rows[pivot][column], -1, prime)
        rows[pivot] = [(inverse * x) % prime for x in rows[pivot]]
        for r in range(len(rows)):
            if r != pivot and rows[r][column]:
                factor = rows[r][column]
                rows[r] = [(x - factor * y) % prime
                           for x, y in zip(rows[r], rows[pivot], strict=True)]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def support_and_punctures() -> tuple[dict, tuple[int, int]]:
    seams, faces = native_loop_readout.carrier()
    support = flux.build_complex(faces, 12)
    if support["edges"] != list(seams) or len(faces) != 20:
        raise ArithmeticError("native and finite-flux carrier incidence differ")
    if support["vertex_count"] - len(seams) + len(faces) != 2:
        raise ArithmeticError("native carrier Euler characteristic changed")
    for vertex in range(12):
        link: dict[int, set[int]] = {}
        for face in faces:
            if vertex in face:
                a, b = (other for other in face if other != vertex)
                link.setdefault(a, set()).add(b)
                link.setdefault(b, set()).add(a)
        if len(link) != 5 or any(len(neighbors) != 2 for neighbors in link.values()):
            raise ArithmeticError("a native vertex link is not a five-cycle")
        reached = {next(iter(link))}
        frontier = list(reached)
        for neighbor in frontier:
            for other in link[neighbor]:
                if other not in reached:
                    reached.add(other)
                    frontier.append(other)
        if reached != set(link):
            raise ArithmeticError("a native vertex link is disconnected")
    distances = flux.graph_distances(support["adjacency"])
    antipodes = []
    for vertex in range(12):
        distant = [other for other in range(12) if distances[vertex][other] == 3]
        if len(distant) != 1:
            raise ArithmeticError("native carrier has no unique vertex antipode")
        antipodes.append(distant[0])
    opposite = frozenset(antipodes[v] for v in faces[0])
    candidates = [i for i, face in enumerate(faces) if frozenset(face) == opposite]
    if len(candidates) != 1:
        raise ArithmeticError("native face has no unique antipodal face")
    return support, (0, candidates[0])


def homology_ranks(support: dict, punctures: tuple[int, int]) -> dict[str, dict]:
    """Compute H1 of the full and two-face-punctured cell complexes over F2/F3."""
    edges = support["edges"]
    face_count = len(support["face_rows"])
    vertex_count = support["vertex_count"]
    if vertex_count - len(edges) + face_count != 2 or len(set(punctures)) != 2:
        raise ValueError("a closed sphere complex and two distinct faces are required")
    boundary_one = [[int(vertex == v) - int(vertex == u) for u, v in edges]
                    for vertex in range(vertex_count)]
    boundary_two = support["boundary"]
    if any(sum(boundary_one[v][e] * boundary_two[e][f]
               for e in range(len(edges))) != 0
           for v in range(vertex_count)
           for f in range(face_count)):
        raise ArithmeticError("the carrier boundaries do not compose to zero")
    retained = [f for f in range(face_count) if f not in punctures]
    punctured_two = [[row[f] for f in retained] for row in boundary_two]
    singly_retained = [f for f in range(face_count) if f != punctures[0]]
    singly_punctured_two = [[row[f] for f in singly_retained]
                            for row in boundary_two]
    result = {}
    for prime in (2, 3):
        rank_one = rank_mod(boundary_one, prime)
        rank_two = rank_mod(boundary_two, prime)
        rank_single = rank_mod(singly_punctured_two, prime)
        rank_punctured = rank_mod(punctured_two, prime)
        full_h1 = len(edges) - rank_one - rank_two
        single_h1 = len(edges) - rank_one - rank_single
        punctured_h1 = len(edges) - rank_one - rank_punctured
        if (rank_one, rank_two, rank_single, rank_punctured,
                full_h1, single_h1, punctured_h1) != (
            vertex_count - 1, face_count - 1, face_count - 1,
            face_count - 2, 0, 0, 1
        ):
            raise ArithmeticError("the full/once/twice-punctured sphere homology changed")
        result[f"F{prime}"] = {
            "boundary_one_rank": rank_one,
            "full_boundary_two_rank": rank_two,
            "one_puncture_boundary_two_rank": rank_single,
            "punctured_boundary_two_rank": rank_punctured,
            "full_H1_dimension": full_h1,
            "one_puncture_H1_dimension": single_h1,
            "two_puncture_H1_dimension": punctured_h1,
        }
    return result


def flat_gauge_potential(support: dict, chain: dict[int, int]) -> list[int]:
    """Recover a vertex gauge for a flat Z6 chain on the closed carrier."""
    faces = len(support["face_rows"])
    if flux.face_holonomies(chain, support["boundary"], faces) != [0] * faces:
        raise ValueError("the difference chain has nonzero face curvature")
    potential: dict[int, int] = {0: 0}
    queue = [0]
    while queue:
        u = queue.pop(0)
        for edge, (a, b) in enumerate(support["edges"]):
            if u == a and b not in potential:
                potential[b] = (potential[a] + chain.get(edge, 0)) % 6
                queue.append(b)
            elif u == b and a not in potential:
                potential[a] = (potential[b] - chain.get(edge, 0)) % 6
                queue.append(a)
    if len(potential) != support["vertex_count"]:
        raise ArithmeticError("native seam graph disconnected")
    gauge = [potential[v] for v in range(support["vertex_count"])]
    if any(flux.coboundary_chain(gauge, support["edges"]).get(e, 0)
           != chain.get(e, 0) % 6 for e in range(len(support["edges"]))):
        raise ArithmeticError("a flat chain has nontrivial full-sphere holonomy")
    return gauge


def verify_local_link_path(
    support: dict, punctures: tuple[int, int],
    path: list[tuple[int, int]], charge: int,
) -> dict[int, int]:
    """Build and then erase a declared flux pair by single-seam link changes.

    This is a countermodel to protection by gauge invariance alone. It is not
    a native repair history: native repair has no Z6 seam-link register.
    """
    if charge not in range(1, 6):
        raise ValueError("a nonzero Z6 charge is required")
    start, end = punctures
    boundary = support["boundary"]
    face_count = len(support["face_rows"])
    chain: dict[int, int] = {}
    cursor = start
    prefix = []
    for following, edge in path:
        if (not 0 <= following < face_count or not 0 <= edge < len(support["edges"])
                or boundary[edge][cursor] not in (-1, 1)
                or boundary[edge][following] != -boundary[edge][cursor]
                or edge in chain):
            raise ArithmeticError("a declared dual step is not one adjacent seam update")
        chain[edge] = (boundary[edge][cursor] * charge) % 6
        cursor = following
        prefix.append((cursor, edge))
        expected = [0] * face_count
        expected[start] = charge
        expected[cursor] = (-charge) % 6
        if flux.face_holonomies(chain, boundary, face_count) != expected:
            raise ArithmeticError("a local seam update failed to move the flux endpoint")
    if cursor != end:
        raise ArithmeticError("the declared dual path misses its second puncture")
    complete = dict(chain)
    for index in range(len(prefix) - 1, -1, -1):
        _, edge = prefix[index]
        del chain[edge]
        endpoint = start if index == 0 else prefix[index - 1][0]
        expected = [0] * face_count
        if endpoint != start:
            expected[start] = charge
            expected[endpoint] = (-charge) % 6
        if flux.face_holonomies(chain, boundary, face_count) != expected:
            raise ArithmeticError("one local inverse update failed to erase the pair")
    if chain:
        raise ArithmeticError("the declared pair did not annihilate")
    return complete


def run() -> dict:
    support, punctures = support_and_punctures()
    ranks = homology_ranks(support, punctures)
    refined, child_faces = flux.edge_midpoint_refinement(support)
    refined_punctures = (child_faces[punctures[0]][0],
                         child_faces[punctures[1]][0])
    refined_ranks = homology_ranks(refined, refined_punctures)
    if (refined["vertex_count"], len(refined["edges"]),
            len(refined["face_rows"])) != (42, 120, 80):
        raise ArithmeticError("the exact midpoint refinement counts changed")
    child_unit_chain = flux.flux_tube_chain(refined, *refined_punctures, 1)
    for charge in range(6):
        child_chain = {e: (charge * value) % 6
                       for e, value in child_unit_chain.items()}
        child_holonomies = flux.face_holonomies(child_chain,
                                                refined["boundary"], 80)
        expected = [0] * 80
        expected[refined_punctures[0]] = charge
        expected[refined_punctures[1]] = (-charge) % 6
        if child_holonomies != expected:
            raise ArithmeticError("a declared child flux sector is absent")
    if tuple(native.State.__dataclass_fields__) != ("loads", "prefix"):
        raise ArithmeticError("the native state has an unexamined register")
    seams, _ = native_loop_readout.carrier()
    checked = 0
    for loads in native.states():
        contrast = native_loop_readout.read_contrasts(loads, seams)
        chain = {e: contrast[edge] % 6 for e, edge in enumerate(seams)}
        if flux.face_holonomies(chain, support["boundary"], 20) != [0] * 20:
            raise ArithmeticError("native current generated a face flux")
        flat_gauge_potential(support, chain)
        checked += len(seams) * 2
    path = flux.dual_path(support, *punctures)
    alternative = flux.dual_path(support, *punctures,
                                 avoid=frozenset({path[0][1]}))
    if path == alternative:
        raise ArithmeticError("a second declared dual route was not found")
    for charge in range(1, 6):
        first = verify_local_link_path(support, punctures, path, charge)
        second = verify_local_link_path(support, punctures, alternative, charge)
        difference = {e: (first.get(e, 0) - second.get(e, 0)) % 6
                      for e in range(len(seams))}
        flat_gauge_potential(support, difference)
        if any(flux.face_holonomies(first, support["boundary"], 20)[f]
               for f in range(20) if f not in punctures):
            raise ArithmeticError("declared puncture connection has bulk face curvature")
    return {
        "schema": "oph.native_defect_stability_gate.v1",
        "carrier": {"vertices": 12, "seams": 30, "faces": 20,
                    "euler_characteristic": 2,
                    "all_vertex_links_five_cycles": True,
                    "declared_antipodal_puncture_faces": list(punctures)},
        "homology": ranks,
        "one_step_topology_refinement": {
            "child_vertices": 42, "child_seams": 120, "child_faces": 80,
            "declared_child_puncture_faces": list(refined_punctures),
            "declared_child_flux_menu": list(range(6)),
            "homology": refined_ranks,
            "native_repair_on_child_complex": False,
            "coarse_link_connection_transport_checked": False,
        },
        "native_record": {"load_preparations": len(native.states()),
                          "coin_current_branches_checked": checked,
                          "native_state_registers": ["loads", "prefix"],
                          "all_face_fluxes_mod_6": 0,
                          "current_selected_nontrivial_Z6_defect": False},
        "declared_Z6_link_control": {
            "two_routes_checked_per_nonzero_charge": 2,
            "charges_checked": [1, 2, 3, 4, 5],
            "route_lengths": [len(path), len(alternative)],
            "flat_on_all_unpunctured_faces": True,
            "route_difference_is_vertex_regauge": True,
            "single_seam_updates_move_and_annihilate_pair": True,
            "annihilation_updates_include_puncture_boundary_seams": True,
            "single_seam_link_updates_are_native_repair": False,
            "all_depth_refinement_or_dynamical_protection_proved": False,
        },
        "scope": (
            "The full closed carrier has no H1 loop sector; a declared two-face "
            "puncture gives one independent Z6-valued loop generator. Native "
            "current contrasts remain "
            "pure gauge. The constructed nonzero flux pair is a separate link "
            "model and is removable if arbitrary local seam-link updates are "
            "admitted. Gauge invariance alone does not protect a defect."
        ),
        "all_pass": True,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))

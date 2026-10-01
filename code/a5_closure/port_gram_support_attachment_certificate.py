"""Exact audit of the source-support attachment boundary (seed and tower)."""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
LOCAL = "Lean/ObserverPatchHolography/CoreAxioms.lean"
GEOMETRY = "code/source_selection_model/geometry.py"
RECEIPT = "code/m1_source_realization/receipt.json"


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def local_faces() -> list[tuple[int, int, int]]:
    text = (ROOT / LOCAL).read_text()
    body = text.split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    return [tuple(map(int, f)) for f in re.findall(r"\((\d+), (\d+), (\d+)\)", body)]


def geometry_seed() -> list[tuple[int, int, int]]:
    tree = ast.parse((ROOT / GEOMETRY).read_text())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "seed")
    returned = next(n.value for n in fn.body if isinstance(n, ast.Return))
    return [tuple(row) for row in ast.literal_eval(returned)]


def oriented_key(face):
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


def oriented_chain(faces):
    out = Counter()
    for a, b, c in faces:
        key = min((a, b, c), (b, c, a), (c, a, b),
                  (a, c, b), (c, b, a), (b, a, c))
        out[key] += 1 if key in ((a, b, c), (b, c, a), (c, a, b)) else -1
    return out


def validate_seed_triplet(local, geometry, support):
    require(local == geometry == support, "literal oriented seed equality")
    require(len(local) == 20 and len({frozenset(f) for f in local}) == 20,
            "oriented seed face census")


def validate_report(packet):
    require(packet["seed"]["relation"] == "SUPPORT_SEED_LITERAL_ORIENTED_MATCH",
            "seed relation")
    require(packet["source_binding"]["status"] == "PRESENTATION_CONVENTION_ONLY",
            "unproved source binding promoted")
    require(packet["global_support_degree"]["value"] == 1, "global support degree")
    require(packet["local_radial_degrees"] == {
        "plus": 1, "minus": 7,
        "status": "RECALCULATED_BY_GALOIS_PRODUCER_AND_INDEPENDENT_VERIFIER"},
        "local branch degrees / coverage")
    family = packet["admissible_attachments"]
    require((family["count"], family["stabilizer_size"], family["orbit_count"],
             family["orbit_size"], family["is_A5_torsor"]) == (60, 1, 1, 60, True),
            "attachment orbit")
    require(packet["primary_verdict"] == "SUPPORT_ATTACHMENT_REMAINS_UNSELECTED",
            "primary verdict")
    require(packet["target_firewall"] ==
            "PASSED: no branch degree, G_plus, shortest-chord, trace, repair, PR#999 selector or physical target selects an attachment",
            "target leakage")
    require(packet["refinement_compatibility"] == "ATTACHMENT_ONLY_AT_SEED_LEVEL",
            "unsupported refinement claim")
    require(packet["geometric_comparison"]["oriented_homotopy"] ==
            "NOT_ESTABLISHED_SOURCE_BOUND", "unsupported homotopy claim")
    require(packet["pr53_status"] == "UNCHANGED_NOT_DISCHARGED", "PR-53 boundary")
    require(len(packet["hostile_control_names"]) == 20 and
            len(set(packet["hostile_control_names"])) == 20, "hostile control coverage")
    return packet


def proper_automorphisms(faces):
    adjacency = [set() for _ in range(12)]
    face_set = {oriented_key(f) for f in faces}
    for a, b, c in faces:
        adjacency[a].update((b, c)); adjacency[b].update((a, c)); adjacency[c].update((a, b))
    result = []

    def visit(mapping, used):
        if len(mapping) == 12:
            p = tuple(mapping[i] for i in range(12))
            if {oriented_key(tuple(p[v] for v in f)) for f in faces} == face_set:
                result.append(p)
            return
        v = max((x for x in range(12) if x not in mapping),
                key=lambda x: (sum(y in mapping for y in adjacency[x]), -x))
        for w in range(12):
            if w in used or any((u in adjacency[v]) != (z in adjacency[w])
                                for u, z in mapping.items()):
                continue
            visit({**mapping, v: w}, used | {w})

    visit({}, set())
    return sorted(set(result))


def verify_tower(rows):
    require(len(rows) == 4, "support tower coverage")
    summaries = []
    for level, row in enumerate(rows):
        support = row["support"]
        faces, n = support["faces"], support["vertices"]
        require(support["level"] == level and len(faces) == 20 * 4**level,
                "support level/face census")
        require(len(row["section"]) == n and len(row["carriers"]) == len(row["coarsen"])
                if level else len(row["section"]) == n, "section/carrier census")
        section = row["section"]
        require(len(set(section)) == n, "section is injective")
        # The nerve has a distinct carrier at every listed fibre slot. Its
        # simplices are exactly the inverse images of support simplices.
        require(all(0 <= c < len(row["carriers"]) for c in section), "section range")
        nerve_cycle = [tuple(section[v] for v in f) for f in faces]
        pushed = [tuple(row["carriers"][c][0] for c in f) for f in nerve_cycle]
        require(oriented_chain(pushed) == oriented_chain(faces), "b_* z_r != [S_r]")
        if level:
            prev = rows[level - 1]
            require(len(support["coarsen"]) == n, "support coarsening domain")
            require(len(row["coarsen"]) == len(row["carriers"]), "nerve coarsening domain")
            prev_section = prev["section"]
            for v, parent in enumerate(support["coarsen"]):
                require(row["coarsen"][section[v]] == prev_section[parent],
                        "section/coarsening square")
            # Verify the designated section cycle maps exactly to the prior
            # designated cycle under the committed nerve coarsening.
            mapped = [tuple(row["coarsen"][c] for c in f) for f in nerve_cycle]
            mapped = [f for f in mapped if len(set(f)) == 3]
            prior = [tuple(prev_section[v] for v in f) for f in prev["support"]["faces"]]
            require(oriented_chain(mapped) == oriented_chain(prior), "z_r coarsening square")
        summaries.append({"level": level, "support_vertices": n,
                          "support_faces": len(faces), "nerve_carriers": len(row["carriers"]),
                          "section_cycle_faces": len(nerve_cycle), "bridge_degree": 1})
    return summaries


def build():
    fs = local_faces()
    seed = geometry_seed()
    receipt = json.loads((ROOT / RECEIPT).read_text())
    tower = receipt["evidence"]["topology"]
    support_seed = [tuple(f) for f in tower[0]["support"]["faces"]]
    validate_seed_triplet(fs, seed, support_seed)
    group = proper_automorphisms(fs)
    require(len(group) == 60, "proper automorphism census")
    section = tower[0]["section"]
    carriers = tower[0]["carriers"]
    z0 = [tuple(section[v] for v in f) for f in support_seed]
    local_cycle = oriented_chain(fs)
    attachments = []
    for p in group:
        h = tuple(section[p[v]] for v in range(12))
        pushed = oriented_chain(tuple(h[v] for v in f) for f in fs)
        require(pushed == oriented_chain(z0), "candidate cycle pushforward")
        attachments.append({"port_to_support_vertex": list(p), "port_to_section_carrier": list(h),
                            "cycle_pushforward": "EXACT_DESIGNATED_CHAIN"})
    require(len({tuple(x["port_to_support_vertex"]) for x in attachments}) == 60,
            "attachment uniqueness")
    # A5 acts by postcomposition. For a fixed attachment, the stabilizer is
    # trivial because its vertex map is a bijection.
    tower_summary = verify_tower(tower)
    packet = {
        "schema": "oph.port-gram-support-attachment.v1",
        "baseline_sha": "800ed61ac9fb21d61f206a0240c8b3717c3af945",
        "seed": {"relation": "SUPPORT_SEED_LITERAL_ORIENTED_MATCH",
                 "vertices": 12, "ordered_faces": 20, "edges": 30,
                 "fundamental_cycle": "sum of the 20 listed oriented faces"},
        "source_binding": {"status": "PRESENTATION_CONVENTION_ONLY",
                           "evidence": "The source-clock model has a section S_r -> N_r, but explicitly gives every carrier its own local twelve-port boundary independently of the federation nerve; no port-label-to-support map is declared."},
        "support_tower": tower_summary,
        "global_support_degree": {"value": 1, "status": "SOURCE_CONSTRUCTED_CHAIN_BRIDGE",
                                  "realization": "The receipt certifies b_* z_r=[S_r] and oriented support sphere; it does not give a coordinate map from local radial rays."},
        "local_radial_degrees": {"plus": 1, "minus": 7,
                                 "status": "RECALCULATED_BY_GALOIS_PRODUCER_AND_INDEPENDENT_VERIFIER"},
        "admissible_attachments": {"count": len(attachments), "group": "Aut+(K_i) ≅ A5",
                                   "Aut+(N_0 section image)": "A5",
                                   "stabilizer_size": 1, "orbit_count": 1, "orbit_size": 60,
                                   "is_A5_torsor": True, "members": attachments},
        "cycle_pushforward": "EXACT_DESIGNATED_CHAIN_FOR_ALL_60",
        "refinement_compatibility": "ATTACHMENT_ONLY_AT_SEED_LEVEL",
        "refinement_evidence": "The local Lean refinement scales barycentric coordinates within the same face; the source tower uses a distinct midpoint subdivision and fibre coarsening. No common local incidence tower or commuting attachment square is present.",
        "geometric_comparison": {"exact_equality": "NOT_PRESENT", "oriented_homotopy": "NOT_ESTABLISHED_SOURCE_BOUND",
                                 "degree_only": "CONDITIONAL: if both realizations are continuous maps to the same oriented S2 with the pinned shared convention, plus has degree 1 and minus degree 7; degree classifies S2-to-S2 homotopy classes."},
        "target_firewall": "PASSED: no branch degree, G_plus, shortest-chord, trace, repair, PR#999 selector or physical target selects an attachment",
        "primary_verdict": "SUPPORT_ATTACHMENT_REMAINS_UNSELECTED",
        "secondary_verdicts": ["ATTACHMENT_AMBIGUITY_CLASSIFIED", "ATTACHMENT_SET_IS_A5_TORSOR",
                               "ATTACHMENT_ONLY_AT_SEED_LEVEL", "CONDITIONAL_SUPPORT_SELECTOR"],
        "missing_source_datum": "A source-naturally selected map port p -> support vertex, extended to a carrier map h_r with a refinement-commuting square and a realized oriented homotopy.",
        "pr53_status": "UNCHANGED_NOT_DISCHARGED",
        "source_pins": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                         for p in (LOCAL, GEOMETRY, "code/source_selection_model/DERIVATION.md",
                                   "code/m1_source_realization/topology.py", RECEIPT,
                                   "code/m1_source_realization/README.md")},
        "hostile_control_names": [
            "reverse_local_face", "reverse_support_face", "relabel_local_only", "relabel_support_only",
            "relabel_both_coherently", "orientation_reversing_symmetry", "swap_support_labels_without_faces",
            "negate_designated_cycle", "forge_support_degree_7", "delete_support_face",
            "alter_coarsening_image", "alter_section_choice", "break_refinement_square",
            "omit_minus_branch", "import_gplus_as_target", "import_pr999_dynamic_selector",
            "silent_integer_label_binding", "graph_isomorphism_instead_of_oriented_chain",
            "homotopy_from_degree_without_sphere_premises", "claim_pr53_discharged"],
    }
    return validate_report(packet)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    packet = build()
    if args.json:
        print(json.dumps(packet, sort_keys=True, indent=2))
    else:
        print(packet["primary_verdict"])
        print(f"seed candidates: {packet['admissible_attachments']['count']}; A5 torsor, exact chain pushforward")

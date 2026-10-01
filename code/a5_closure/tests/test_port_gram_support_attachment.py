"""Hostile and metamorphic controls for the support-attachment audit."""

import copy
import json
from pathlib import Path
import sys

import pytest

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import port_gram_support_attachment_certificate as producer
import verify_port_gram_support_attachment_independent as independent


ROOT = Path(__file__).resolve().parents[3]


def packet():
    return producer.build()


def require_rejects(call, *args):
    with pytest.raises((ValueError, KeyError, IndexError)):
        call(*args)


def test_literal_oriented_seed_equality_and_census():
    result = packet()
    assert result["seed"]["relation"] == "SUPPORT_SEED_LITERAL_ORIENTED_MATCH"
    assert (result["seed"]["vertices"], result["seed"]["ordered_faces"],
            result["seed"]["edges"]) == (12, 20, 30)


def test_exact_attachment_torsor_and_cycle_pushforward():
    result = packet()["admissible_attachments"]
    assert result["count"] == result["orbit_size"] == 60
    assert result["stabilizer_size"] == 1 and result["orbit_count"] == 1
    assert result["is_A5_torsor"]
    assert len({tuple(x["port_to_support_vertex"]) for x in result["members"]}) == 60
    assert all(x["cycle_pushforward"] == "EXACT_DESIGNATED_CHAIN" for x in result["members"])
    assert result["fixed_full_support_row_aut_size"] == 5
    assert result["raw_maps_mod_fixed_mark_aut_orbits"] == 12
    assert result["host_mark_orbit_size_under_A5"] == 12


def test_independent_reconstruction_does_not_consume_candidate_list():
    result = independent.verify()
    assert result["producer_candidate_list_consumed"] is False
    assert result["candidate_count"] == result["exact_chain_pushforwards"] == 60


def test_source_tower_degree_and_refinement_square():
    result = packet()
    assert [x["bridge_degree"] for x in result["support_tower"]] == [1, 1, 1, 1]
    assert result["refinement_compatibility"] == "ATTACHMENT_ONLY_AT_SEED_LEVEL_IN_PINNED_AUDITED_CORPUS"


def test_reverse_one_local_face_fails_literal_seed_gate():
    local = producer.local_faces()
    support = producer.geometry_seed()
    local[0] = (local[0][0], local[0][2], local[0][1])
    section = list(range(12))
    assert producer.classify_seed_attachments(local, support, section) == []


def test_reverse_one_support_face_fails_literal_seed_gate():
    local = producer.local_faces()
    support = producer.geometry_seed()
    support[0] = (support[0][0], support[0][2], support[0][1])
    section = list(range(12))
    assert producer.classify_seed_attachments(local, support, section) == []


@pytest.mark.parametrize("which", ["local", "support", "both"])
def test_proper_a5_relabelings_recompute_exact_candidate_pipeline(which):
    faces = producer.local_faces()
    group = producer.proper_automorphisms(faces)
    permutation = group[1]
    local, support = list(faces), list(faces)
    if which in ("local", "both"):
        local = [tuple(permutation[v] for v in f) for f in local]
    if which in ("support", "both"):
        support = [tuple(permutation[v] for v in f) for f in support]
    row = json.loads((ROOT / "code/m1_source_realization/receipt.json").read_text())["evidence"]["topology"][0]
    attachments = producer.classify_seed_attachments(local, support, row["section"])
    assert len(attachments) == 60
    expected = producer.oriented_chain(
        tuple(row["section"][v] for v in f) for f in support)
    for item in attachments:
        h = item["port_to_section_carrier"]
        pushed = producer.oriented_chain(tuple(h[v] for v in f) for f in local)
        assert pushed == expected


def test_orientation_reversing_icosahedral_map_is_not_admissible():
    faces = producer.local_faces()
    adjacency = [set() for _ in range(12)]
    unoriented = {frozenset(f) for f in faces}
    for a, b, c in faces:
        adjacency[a].update((b, c)); adjacency[b].update((a, c)); adjacency[c].update((a, b))
    found = []
    def visit(mapping, used):
        if len(mapping) == 12:
            p = tuple(mapping[i] for i in range(12))
            candidate = [tuple(p[v] for v in f) for f in faces]
            if ({frozenset(f) for f in candidate} == unoriented and
                    producer.oriented_chain(candidate) != producer.oriented_chain(faces)):
                found.append(p)
            return
        v = next(x for x in range(12) if x not in mapping)
        for w in range(12):
            if w in used or len(adjacency[v]) != len(adjacency[w]):
                continue
            if all((u in adjacency[v]) == (z in adjacency[w]) for u, z in mapping.items()):
                visit({**mapping, v: w}, used | {w})
    visit({}, set())
    assert found
    reversing_map = found[0]
    with pytest.raises(ValueError):
        producer.classify_seed_attachments(faces, faces, list(range(12)),
                                           proposed_map=reversing_map)


def test_unaccompanied_support_label_swap_fails_seed_gate():
    faces = producer.local_faces()
    swap = [1, 0] + list(range(2, 12))
    with pytest.raises(ValueError):
        producer.classify_seed_attachments(faces, faces, list(range(12)), proposed_map=swap)


def test_negated_designated_cycle_fails_integer_chain_equality():
    row = json.loads((ROOT / "code/m1_source_realization/receipt.json").read_text())["evidence"]["topology"][0]
    z = [tuple(row["section"][v] for v in f) for f in row["support"]["faces"]]
    negated = [tuple(reversed(face)) for face in z]
    with pytest.raises(ValueError):
        producer.classify_seed_attachments(producer.local_faces(), producer.local_faces(),
                                           row["section"], designated_cycle=negated)


def test_support_degree_seven_forgery_is_rejected_by_pinned_value():
    result = copy.deepcopy(packet())
    forged = copy.deepcopy(result["global_support_degree"])
    forged["value"] = 7
    result["global_support_degree"] = forged
    require_rejects(producer.validate_report, result)


def test_deleted_support_face_fails_seed_census():
    faces = producer.geometry_seed()[:-1]
    assert producer.classify_seed_attachments(producer.local_faces(), faces,
                                              list(range(12))) == []


@pytest.mark.parametrize("mutation", ["coarsening", "section", "square"])
def test_support_tower_mutations_fail_chain_or_commuting_square(mutation):
    rows = json.loads((ROOT / "code/m1_source_realization/receipt.json").read_text())["evidence"]["topology"]
    rows = copy.deepcopy(rows)
    if mutation == "coarsening":
        rows[1]["support"]["coarsen"][0] = 1
    elif mutation == "section":
        rows[0]["section"][0] = 1
    else:
        rows[1]["coarsen"][0] = 1
    require_rejects(producer.verify_tower, rows)


def test_both_galois_branches_are_mandatory():
    result = copy.deepcopy(packet())
    del result["local_radial_degrees"]["minus"]
    require_rejects(producer.validate_report, result)


def test_branch_degree_and_dynamic_selector_are_not_attachment_inputs():
    result = copy.deepcopy(packet())
    assert result["source_binding"]["status"] == "NO_ATTACHMENT_FOUND_IN_PINNED_AUDITED_CORPUS"
    assert "PR#999 selector" in result["target_firewall"]
    assert "G_plus" in result["target_firewall"]
    result["target_firewall"] = "PR#999 selects the attachment"
    require_rejects(producer.validate_report, result)


def test_integer_label_coincidence_is_not_source_binding():
    result = copy.deepcopy(packet())
    result["source_binding"]["status"] = "EXPLICIT_SOURCE_MAP"
    require_rejects(producer.validate_report, result)


def test_expected_source_hashes_are_enforced_fail_closed():
    hashes = producer.verify_source_pins()
    assert hashes == producer.PINNED_SOURCE_SHA256
    path = "code/source_selection_model/DERIVATION.md"
    require_rejects(producer.verify_source_pins, ROOT, {path: b"invented port-to-support theorem"})


def test_galois_degrees_are_consumed_from_pinned_reference():
    refpath = ROOT / "code/a5_closure/manifests/galois_port_frame_reference.json"
    reference = json.loads(refpath.read_text())
    assert producer.galois_degrees(reference, producer.local_faces()) == {"plus": 1, "minus": 7}
    report = packet()["local_radial_degrees"]
    assert report["status"] == "CONSUMED_FROM_PINNED_GALOIS_REFERENCE_AND_INDEPENDENTLY_REPLAYED"
    reference["family"]["plus"]["degree_controls"][0]["degree"] = 7
    forged = json.dumps(reference).encode()
    require_rejects(producer.verify_source_pins, ROOT,
                    {"code/a5_closure/manifests/galois_port_frame_reference.json": forged})


def test_attachment_absence_is_scoped_to_byte_pinned_corpus():
    source = packet()["source_binding"]
    assert source["status"] == "NO_ATTACHMENT_FOUND_IN_PINNED_AUDITED_CORPUS"
    assert tuple(source["corpus"]) == producer.ATTACHMENT_AUDIT_CORPUS


def test_unoriented_graph_match_does_not_replace_chain_equality():
    local = producer.local_faces()
    reversed_one = list(local)
    a, b, c = reversed_one[0]
    reversed_one[0] = (a, c, b)
    assert {frozenset(f) for f in local} == {frozenset(f) for f in reversed_one}
    assert producer.oriented_chain(local) != producer.oriented_chain(reversed_one)


def test_degree_comparison_does_not_claim_homotopy_without_shared_map():
    result = copy.deepcopy(packet())
    result["geometric_comparison"]["oriented_homotopy"] = "PROVED"
    require_rejects(producer.validate_report, result)


def test_pr53_remains_unchanged():
    result = copy.deepcopy(packet())
    result["pr53_status"] = "DISCHARGED"
    require_rejects(producer.validate_report, result)


def test_selector_and_target_leak_mutations_fail_closed():
    result = copy.deepcopy(packet())
    result["target_firewall"] = "choose h because G_plus has degree one"
    result["primary_verdict"] = "SOURCE_SUPPORT_ATTACHMENT_PRESENT"
    require_rejects(producer.validate_report, result)


def test_claimed_refinement_extension_fails_closed():
    result = copy.deepcopy(packet())
    result["refinement_compatibility"] = "ATTACHMENT_REFINES_NATURALLY"
    require_rejects(producer.validate_report, result)

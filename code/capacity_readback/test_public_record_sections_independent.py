"""Exhaustive gluing controls and producer-independent completeness replay."""
from copy import deepcopy
from itertools import product
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from correctable_public_record_capacity import evaluate_terminal, section_id
from public_record_csp import decode_section_id, encode_section_id, public_global_sections_csp
from verify_public_record_sections import ReplayBudgetExceeded, verify_public_record_sections


MAPS = tuple(dict(zip(("0", "1"), values)) for values in product(("x", "y"), repeat=2))
READOUT_PAIRS = tuple(product(MAPS, repeat=2))


def seam(left, right, pair):
    a, b = READOUT_PAIRS[pair]
    return {"left_observer": left, "right_observer": right,
            "left_readout": a, "right_readout": b, "interface_atoms": ["x", "y"]}


def brute(observers, interfaces):
    ids = sorted(observers)
    expected = set()
    for values in product(*(observers[x] for x in ids)):
        candidate = dict(zip(ids, values))
        if all(s["left_readout"][candidate[s["left_observer"]]] ==
               s["right_readout"][candidate[s["right_observer"]]] for s in interfaces):
            expected.add(values)
    return expected


def checked(observers, interfaces):
    sections = public_global_sections_csp(observers, interfaces)
    ids = sorted(observers)
    actual = {tuple(s[x] for x in ids) for s in sections}
    assert len(actual) == len(sections)
    assert actual == brute(observers, interfaces)
    assert verify_public_record_sections(observers, interfaces, sections)
    return sections


@pytest.mark.parametrize("pairs", tuple(product(range(16), repeat=3)))
def test_every_binary_readout_triangle(pairs):
    observers = {x: ["0", "1"] for x in "abc"}
    interfaces = [seam(a, b, k) for (a, b), k in zip((("a", "b"), ("b", "c"), ("c", "a")), pairs)]
    checked(observers, interfaces)


@pytest.mark.parametrize("pairs", tuple(product(range(16), repeat=2)))
@pytest.mark.parametrize("self_interfaces", [False, True])
def test_every_pair_of_parallel_or_self_interfaces(pairs, self_interfaces):
    observers = {x: ["0", "1"] for x in ("a" if self_interfaces else "ab")}
    target = "a" if self_interfaces else "b"
    checked(observers, [seam("a", target, k) for k in pairs])


def test_pairwise_feasible_odd_cycle_has_no_global_record():
    observers = {x: ["0", "1"] for x in "abc"}
    interfaces = [seam(a, b, 6) for a, b in [("a", "b"), ("b", "c"), ("c", "a")]]
    # Pair 6 is identity versus flip, so each edge separately has two sections.
    for s in interfaces:
        local = {x: observers[x] for x in [s["left_observer"], s["right_observer"]]}
        assert len(public_global_sections_csp(local, [s])) == 2
    assert checked(observers, interfaces) == []
    assert verify_public_record_sections(observers, interfaces, [], max_candidates=2)


def test_orientation_and_opaque_relabeling_do_not_change_sections():
    observers = {"a": ["0", "1"], "b": ["0", "1"]}
    original = [seam("a", "b", 6)]
    expected = checked(observers, original)
    names = {"a": "observer|a", "b": "observer=b"}
    atoms = {"0": "value%7C", "1": "value|b=x"}
    relabeled = {names[x]: [atoms[a] for a in values] for x, values in observers.items()}
    s = original[0]
    reverse = [{"left_observer": names["b"], "right_observer": names["a"],
                "left_readout": {atoms[a]: y for a, y in s["right_readout"].items()},
                "right_readout": {atoms[a]: y for a, y in s["left_readout"].items()}}]
    actual = checked(relabeled, reverse)
    assert {encode_section_id(x) for x in actual} == {
        encode_section_id({names[o]: atoms[a] for o, a in x.items()}) for x in expected}


@pytest.mark.parametrize("labels", tuple(product(["x", "%", "|", "=", "%7C", "λ", "x|b=y"], repeat=3)))
def test_identifier_has_a_unique_canonical_round_trip(labels):
    observer, a, b = labels
    section = {observer: a, observer+"second": b}
    identifier = section_id(section)
    assert identifier == encode_section_id(section)
    assert decode_section_id(identifier) == section
    assert encode_section_id(decode_section_id(identifier)) == identifier


@pytest.mark.parametrize("bad", ["", "a=", "=x", "a=x=a", "a=x|a=y", "b=x|a=y", "a=%7c", "a=%", "a=%20"])
def test_identifier_rejects_ambiguous_or_noncanonical_spellings(bad):
    with pytest.raises(ValueError):
        decode_section_id(bad)


def test_delimiter_labels_preserve_four_records_through_capacity():
    observers = {"a": ["x|b=y", "x"], "b": ["z", "y|b=z"]}
    sections = checked(observers, [])
    ids = [section_id(s) for s in sections]
    assert len(set(ids)) == 4
    packet = {"observers": observers, "interfaces": [], "capacity_dimension": 4,
              "reachability_witnesses": {x: ["write"] for x in ids},
              "publicness_policy": [["a", "b"]], "continuation_manifest_complete": True,
              "local_marginal_consistency_passed": True,
              "global_checkpoint_kernels": [{"authorized_observers": ["a", "b"],
                  "continuation_id": "identity", "rows": {x: {x: 1} for x in ids}}],
              "projection_supports": {x: [i] for i, x in enumerate(ids)}}
    result = evaluate_terminal(packet)
    assert result["status"] == "PASS"
    assert result["public_global_section_count"] == result["exact_zero_error_capacity"] == 4


def test_long_connected_diagram_propagates_with_three_search_nodes():
    observers = {f"o{i}": ["0", "1"] for i in range(1200)}
    interfaces = [seam(f"o{i}", f"o{i+1}", 5) for i in range(1199)]
    sections = public_global_sections_csp(observers, interfaces, max_search_nodes=3)
    assert len(sections) == 2
    assert all(len(set(s.values())) == 1 for s in sections)
    assert verify_public_record_sections(observers, interfaces, sections, max_candidates=2)


@pytest.mark.parametrize("limit", [{"max_search_nodes": 1}, {"max_sections": 1}])
def test_exhaustion_does_not_return_a_partial_record_set(limit):
    with pytest.raises(ValueError, match="no complete record set certified"):
        public_global_sections_csp({"a": ["0", "1"]}, [], **limit)


def test_independent_replay_refuses_unsupported_size():
    observers = {"a": ["0", "1"], "b": ["0", "1"]}
    interfaces = [seam("a", "b", 0)]  # constant maps: four sections, not a rooted bijection
    with pytest.raises(ReplayBudgetExceeded):
        verify_public_record_sections(observers, interfaces, [], max_candidates=2)


def test_spanning_replay_checks_self_interfaces():
    observers = {"a": ["0", "1"], "b": ["0", "1"]}
    interfaces = [seam("a", "b", 5), seam("a", "a", 6)]
    assert verify_public_record_sections(observers, interfaces, [], max_candidates=2)


def test_inconsistent_prefix_does_not_hide_a_malformed_later_interface():
    observers = {"a": ["0", "1"]}
    interfaces = [seam("a", "a", 6), {"left_observer": "a"}]
    with pytest.raises(ValueError):
        public_global_sections_csp(observers, interfaces)
    assert not verify_public_record_sections(observers, interfaces, [])


@pytest.mark.parametrize("mutation", ["constant_pass", "empty", "omit", "duplicate", "extra",
                                     "missing_observer", "unknown_atom", "bad_codomain", "bad_readout"])
def test_completeness_replay_rejects_forged_or_incomplete_sections(mutation):
    observers = {"a": ["0", "1"], "b": ["0", "1"]}
    interfaces = [deepcopy(seam("a", "b", 5))]
    sections = [{"a": "0", "b": "0"}, {"a": "1", "b": "1"}]
    assert verify_public_record_sections(observers, interfaces, sections)
    if mutation == "constant_pass": sections = {"status": "PASS"}
    elif mutation == "empty": sections = []
    elif mutation == "omit": sections.pop()
    elif mutation == "duplicate": sections.append(sections[0])
    elif mutation == "extra": sections.append({"a": "0", "b": "1"})
    elif mutation == "missing_observer": sections[0].pop("a")
    elif mutation == "unknown_atom": sections[0]["a"] = "unknown"
    elif mutation == "bad_codomain": interfaces[0]["interface_atoms"] = ["z"]
    elif mutation == "bad_readout": interfaces[0]["left_readout"].pop("0")
    assert not verify_public_record_sections(observers, interfaces, sections)


def test_full_source_diagram_replays_without_producer_under_optimized_python(tmp_path):
    from source_derived_public_checkpoint_packet import build_source_derived_packet
    packet = build_source_derived_packet()
    observers, interfaces = packet["observers"], packet["interfaces"]
    sections = public_global_sections_csp(observers, interfaces, max_search_nodes=25)
    assert len(sections) == 24
    assert verify_public_record_sections(observers, interfaces, sections, max_candidates=24)
    shutil.copy2(Path(__file__).with_name("verify_public_record_sections.py"), tmp_path)
    (tmp_path/"diagram.json").write_text(json.dumps([observers, interfaces, sections]), encoding="utf-8")
    script = """
import json, sys
from verify_public_record_sections import verify_public_record_sections
with open('diagram.json', encoding='utf-8') as stream:
    observers, interfaces, sections = json.load(stream)
if not verify_public_record_sections(observers, interfaces, sections, max_candidates=24):
    raise SystemExit('valid complete record set rejected')
if verify_public_record_sections(observers, interfaces, sections[:-1], max_candidates=24):
    raise SystemExit('incomplete record set accepted')
if 'public_record_csp' in sys.modules or 'source_derived_public_checkpoint_packet' in sys.modules:
    raise SystemExit('producer imported')
"""
    result = subprocess.run([sys.executable, "-O", "-c", script], cwd=tmp_path,
                            capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout+result.stderr

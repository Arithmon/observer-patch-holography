"""Adversarial verifier and independent full-prefix law checks."""
import ast
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
import subprocess
import sys

import pytest

from . import build, codec, verify


@pytest.fixture(scope="module")
def packet():
    return codec.load(codec.HERE/"controls.json")


def test_complete_committed_receipt(packet):
    codec.equal(verify.verify(packet), codec.load(codec.HERE/"receipt.json"), "receipt")


def test_producer_independence():
    tree = ast.parse((codec.HERE/"verify.py").read_text())
    modules = [n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
    assert not any("build" in name or "tomography" in name for name in modules)


@pytest.mark.parametrize("weights", ([1, 1, 1], [2, 3, 5]))
@pytest.mark.parametrize("targets", ([[1, 0], [0, 1]], [[1, 1]], [[0, 0]]))
def test_every_prefix_mass_against_full_words(weights, targets):
    spec = {"ports": [0, 1, 2, 3], "weights": weights, "targets": targets, "horizon": 5}
    planner = build.Planner(spec)
    for horizon in range(6):
        words = list(product(range(3), repeat=horizon))
        good = [w for w in words if verify.classify(4, w, targets) is not None]
        for k in range(horizon+1):
            for prefix in product(range(3), repeat=k):
                state, observed = planner.initial, build.basis((planner.initial[-1],))
                for a in prefix:
                    state, observed = planner.advance(state, observed, a)
                expected = 0
                for word in good:
                    if word[:k] == prefix:
                        term = 1
                        for a in word[k:]:
                            term *= weights[a]
                        expected += term
                assert planner.mass(horizon-k, state, observed) == expected
                if expected and horizon > k:
                    masses = [weight*planner.mass(horizon-k-1, *planner.advance(state, observed, a))
                              for a, weight in enumerate(weights)]
                    assert sum(masses) == expected


def test_deadline_selection_is_not_projective(packet):
    rows = {c["spec"]["horizon"]: c for c in packet["cases"] if c["spec"]["name"] == "chain"}
    assert rows[4]["accepted_mass"] == 0
    assert rows[5]["first_edge_masses"] == [0, 2, 0]
    assert rows[6]["first_edge_masses"] == [0, 19, 2]
    # A shorter-deadline law cannot be obtained by truncating this longer one.
    assert F(rows[5]["first_edge_masses"][2], rows[5]["accepted_mass"]) == 0
    assert F(rows[6]["first_edge_masses"][2], rows[6]["accepted_mass"]) == F(2, 21)


def test_aggregate_constraint_does_not_publish_individual_records(packet):
    row = next(c for c in packet["cases"] if c["spec"]["name"] == "aggregate" and c["spec"]["horizon"] == 4)
    assert row["accepted_mass"] == 1
    word = row["histories"][0]["word"]
    assert verify.classify(5, word, [[1, 1]]) == 4
    assert verify.classify(5, word, [[1, 0], [0, 1]]) is None


@pytest.mark.parametrize("mutation", [
    "accepted_mass", "reference_mass", "full_word_count", "first_edge", "first_time", "expected_time",
    "mean_count", "reads", "writes", "samples", "preparation", "omit_records", "deadline", "duplicate_move",
    "ports", "reference", "omit_history", "omit_prefix", "future_sample", "response", "decoder",
    "zero_witness", "gain", "decision_mass", "decision_ticket", "decision_edge", "word", "bool_word",
    "sample", "writer", "payload", "deleted_operation", "reorder_operations", "extra_field",
])
def test_resealed_semantic_mutations_rejected(packet, mutation):
    original = packet["cases"][1]  # feasible chain, five executed means
    item = deepcopy(original)
    history = item["histories"][0]
    if mutation in ("accepted_mass", "reference_mass", "full_word_count"):
        item[mutation] += 1
    elif mutation == "first_edge": item["first_edge_masses"][0] += 1
    elif mutation == "first_time": item["first_publication_masses"][0] += 1
    elif mutation == "expected_time": item["first_publication_expectation"] = "0"
    elif mutation == "mean_count": item["selected_mean_count"] -= 1
    elif mutation in ("reads", "writes"): item["selected_mean_"+mutation] -= 1
    elif mutation == "samples": item["selected_sample_count"] -= 1
    elif mutation == "preparation": item["spec"]["preparation_writes"] -= 1
    elif mutation == "omit_records": item["spec"]["targets"].pop()
    elif mutation == "deadline": item["spec"]["horizon"] -= 1
    elif mutation == "duplicate_move": item["spec"]["weights"].append(1)
    elif mutation == "ports": item["spec"]["ports"][0] = 999
    elif mutation == "reference": item["spec"]["weights"][0] = 2
    elif mutation == "omit_history": item["histories"].pop()
    elif mutation == "omit_prefix": history["certificates"].pop(0)
    elif mutation == "future_sample": history["certificates"][0] = deepcopy(history["certificates"][-1])
    elif mutation == "response": history["responses"][-1][0] = "0"
    elif mutation == "decoder": history["certificates"][-1][0]["coefficients"][0] = "1"
    elif mutation == "zero_witness": history["certificates"][0][0]["perturbation"] = ["0", "0"]
    elif mutation == "gain": history["certificates"][-1][0]["sample_error_gain"] = "0"
    elif mutation == "decision_mass": history["decisions"][0]["masses"][0] = 1
    elif mutation == "decision_ticket": history["decisions"][0]["ticket"] = -1
    elif mutation == "decision_edge": history["decisions"][0]["edge"] = 0
    elif mutation == "word": history["word"][0] = 0
    elif mutation == "bool_word": history["word"][0] = True
    elif mutation == "sample": history["payloads"][0]["samples"][-1] = "0"
    elif mutation == "writer": history["payloads"][0]["events"][-1][2] = -999
    elif mutation == "payload": history["payloads"][0]["payload"][0] = 42
    elif mutation == "deleted_operation": history["payloads"][0]["events"].pop()
    elif mutation == "reorder_operations": history["payloads"][0]["events"].reverse()
    elif mutation == "extra_field": item["unverified_override"] = True
    # A new outer checksum is deliberately not a way to bypass semantic checks.
    assert codec.digest(item) != codec.digest(original)
    with pytest.raises(ValueError):
        verify.check_case(item, original["spec"])


@pytest.mark.parametrize("mutation", ("scope", "pins", "schema", "missing_case", "extra"))
def test_packet_mutations(packet, mutation):
    forged = deepcopy(packet)
    if mutation == "scope": forged["scope"]["canonical_axiom_selection"] = True
    elif mutation == "pins": forged["pins"][codec.SUPPORT] = "0"*64
    elif mutation == "schema": forged["schema"] = True
    elif mutation == "missing_case": forged["cases"].pop()
    else: forged["unchecked"] = 1
    with pytest.raises(ValueError):
        verify.verify(forged)


@pytest.mark.parametrize("text", ('{"x":1,"x":2}', '{"x":1.0}', '{"x":NaN}', '{"x":Infinity}'))
def test_noncanonical_json_rejected(tmp_path, text):
    target = tmp_path/"trash.json"
    target.write_text(text)
    with pytest.raises(ValueError):
        codec.load(target)


def test_cli_rejects_resealed_forgery(packet, tmp_path):
    forged = deepcopy(packet)
    forged["cases"][0]["accepted_mass"] = 1
    path = tmp_path/"forged.json"
    path.write_bytes(codec.canonical(forged))
    result = subprocess.run([sys.executable, "-m", "source_checkpoint_selection.verify", "--controls", str(path)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "full-word census" in result.stderr


def test_false_success_at_infeasible_deadline(packet):
    item = deepcopy(packet["cases"][0])
    item["histories"] = deepcopy(packet["cases"][1]["histories"])
    with pytest.raises(ValueError):
        verify.check_case(item, packet["cases"][0]["spec"])


@pytest.mark.parametrize("field,value", (("means", 0), ("successful_words", 2),
    ("samples", 0), ("first_response", ["1", "0"]), ("final_response", ["0", "1"]),
    ("first_record_noise_gain", "1"), ("second_record_noise_gain", "0"),
    ("word_response_sha256", "0"*64)))
def test_pipeline_semantic_mutations(packet, field, value):
    items = deepcopy(packet["pipelines"])
    items[0][field] = value
    with pytest.raises(ValueError):
        verify.check_pipelines(items)


def test_minimal_deadline_reference_independence(packet):
    rows = {c["spec"]["name"]: c for c in packet["cases"]
            if c["spec"]["horizon"] == 7 and c["spec"]["name"] in ("captured", "tilted")}
    assert rows["captured"]["accepted_mass"] == 5
    assert rows["tilted"]["accepted_mass"] == 5*(1*2**2*3**2*5**2)
    assert [h["word"] for h in rows["captured"]["histories"]] == [h["word"] for h in rows["tilted"]["histories"]]

"""Exact receipt replay, complete polynomial exports and CLI boundary checks."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import pytest

from tools.strict_json import load, loads
from . import halfline, infinite_history as history


@pytest.fixture(scope="module")
def expanded_replay():
    # One actual full recomputation supplies both compact and expanded
    # receipts. CLI boundary tests reuse it without substituting evidence.
    return history.run(retain_polynomials=True)


@pytest.fixture(scope="module")
def compact_replay(expanded_replay):
    result = deepcopy(expanded_replay)
    for configuration in result["fine"]["configurations"]:
        for branch in configuration["retained_first_records"]:
            branch.pop("pre_readout_polynomial")
    return result


def test_retained_exact_receipt_matches_full_recomputation(compact_replay):
    retained = load(Path(__file__).with_name("infinite_results.json"))
    history.compare_receipt(retained, compact_replay)
    # Object key order has no scientific content; retained list order does.
    history.compare_receipt(dict(reversed(list(retained.items()))), compact_replay)


@pytest.mark.parametrize("mutation", (
    "empty", "all_pass", "bool_as_integer", "certified", "initial_norm",
    "branch_norm", "state_hash", "coefficient_count", "first_record_bool",
    "Taylor", "discarded_support", "coefficient_rounding", "event_error",
    "full_output_error", "tiny_interval_edit", "finite_readout_substitution",
    "physical_status", "data_status", "parameters", "coarse_control",
    "disconnected_control", "occupation_table", "missing_branch", "extra_field",
    "nan", "infinity",
))
def test_receipt_rejects_changed_exact_evidence(compact_replay, mutation):
    changed = deepcopy(compact_replay)
    fine = changed["fine"]
    configuration = fine["configurations"][0]
    branch = configuration["retained_first_records"][0]
    if mutation == "empty":
        changed = {}
    elif mutation == "all_pass":
        changed["all_pass"] = False
    elif mutation == "bool_as_integer":
        changed["all_pass"] = 1
    elif mutation == "certified":
        fine["nonzero_response_certified"] = False
    elif mutation == "initial_norm":
        fine["initial_polynomial_norm_squared"] = "1"
    elif mutation == "branch_norm":
        branch["approximate_pre_readout_norm_squared"] = str(
            F(branch["approximate_pre_readout_norm_squared"])+F(1, 10**100))
    elif mutation == "state_hash":
        branch["pre_readout_polynomial_sha256"] = "0"*64
    elif mutation == "coefficient_count":
        branch["pre_readout_nonzero_coefficients"] = float(branch["pre_readout_nonzero_coefficients"])
    elif mutation == "first_record_bool":
        branch["first_record"] = False
    elif mutation in ("Taylor", "discarded_support", "coefficient_rounding"):
        fine["first_segment_vector_errors"][mutation] = "0"
    elif mutation == "event_error":
        configuration["event_probability_error_upper"] = "0"
    elif mutation == "full_output_error":
        configuration["complete_output_half_trace_norm_error_upper"] = "0"
    elif mutation == "tiny_interval_edit":
        fine["infinite_feedback_contrast_interval"][0] = str(
            F(fine["infinite_feedback_contrast_interval"][0])+F(1, 10**100))
    elif mutation == "finite_readout_substitution":
        fine["scope"]["readout"] = "sign of the compressed finite coordinate matrix"
    elif mutation == "physical_status":
        fine["scope"]["spatial_limit_or_physical_services"] = True
    elif mutation == "data_status":
        changed["physical_clock_calibration_or_independent_data"] = True
    elif mutation == "parameters":
        fine["parameters"]["retained_total_occupation"] -= 1
    elif mutation == "coarse_control":
        changed["coarse_control"]["infinite_feedback_contrast_interval"] = ["0", "0"]
    elif mutation == "disconnected_control":
        changed["disconnected_control"]["nonzero_response_certified"] = True
    elif mutation == "occupation_table":
        changed["rectangular_cutoff_comparisons"].pop("64")
    elif mutation == "missing_branch":
        configuration["retained_first_records"].pop()
    elif mutation == "extra_field":
        branch["discarded_failure_outcome"] = True
    elif mutation == "nan":
        configuration["event_probability_error_upper"] = float("nan")
    else:
        configuration["event_probability_error_upper"] = float("inf")
    with pytest.raises(ValueError):
        history.compare_receipt(changed, compact_replay)


@pytest.mark.parametrize("left,right", (
    ({"values": (1, 2)}, {"values": [1, 2]}),
    ({1: "value"}, {"1": "value"}),
    ({"flag": True}, {"flag": 1}),
    ({"counter": 2}, {"counter": 2.}),
    ({"diagnostic": .125}, {"diagnostic": .125}),
    ({"rational": F(1, 2)}, {"rational": "1/2"}),
))
def test_public_replay_does_not_coerce_non_json_or_distinct_json_types(left, right):
    with pytest.raises(ValueError):
        history.compare_receipt(left, right)


def test_public_replay_rejects_json_container_subclasses():
    class Dictionary(dict):
        pass

    with pytest.raises(ValueError):
        history.compare_receipt(Dictionary(all_pass=True), {"all_pass": True})


@pytest.mark.parametrize("text", (
    '{"all_pass":false,"all_pass":true}',
    '{"fine":{"scope":{"independent_empirical_validation":true,"independent_empirical_validation":false}}}',
))
def test_duplicate_keys_are_rejected_at_every_receipt_depth(text):
    with pytest.raises(ValueError, match="duplicate"):
        loads(text)


def test_exported_polynomials_reconstruct_all_raw_branch_weights_and_hashes(expanded_replay):
    fine = expanded_replay["fine"]
    normalization = F(fine["initial_polynomial_norm_squared"])
    cutoff = fine["parameters"]["retained_total_occupation"]
    for configuration in fine["configurations"]:
        assert [row["first_record"] for row in configuration["retained_first_records"]] == [0, 1]
        for branch in configuration["retained_first_records"]:
            rows = branch["pre_readout_polynomial"]
            keys = [tuple(row[:3]) for row in rows]
            assert keys == sorted(set(keys))
            assert len(keys) == branch["pre_readout_nonzero_coefficients"]
            assert all(n0+n1 <= cutoff and 0 <= gauge < 4 for n0, n1, gauge in keys)
            digest = hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()
            assert digest == branch["pre_readout_polynomial_sha256"]
            vector = {(n0, n1, gauge): (F(real), F(imag))
                      for n0, n1, gauge, real, imag in rows}
            mass = halfline.norm_squared(vector)/normalization
            assert mass == F(branch["approximate_pre_readout_norm_squared"])
            positive = halfline.positive_weight_interval(vector, normalization)
            negative = halfline.weight_interval(vector, 0, normalization)
            for name, exact in (("positive", positive), ("negative", negative)):
                low, high = map(F, branch[f"approximate_{name}_weight"])
                assert low <= exact[0] <= exact[1] <= high
            assert positive[0]+negative[1] == mass
            assert positive[1]+negative[0] == mass


def test_polynomial_export_changes_are_binding_even_if_hash_is_replaced(expanded_replay):
    changed = deepcopy(expanded_replay)
    branch = changed["fine"]["configurations"][0]["retained_first_records"][0]
    row = branch["pre_readout_polynomial"][0]
    row[3] = str(F(row[3])+F(1, 10**100))
    # An attacker can update a local digest; recomputation must still bind
    # the actual propagated coefficients rather than trust that digest.
    branch["pre_readout_polynomial_sha256"] = hashlib.sha256(json.dumps(
        branch["pre_readout_polynomial"], separators=(",", ":")).encode()).hexdigest()
    with pytest.raises(ValueError):
        history.compare_receipt(changed, expanded_replay)


def _cli(monkeypatch, args, replay):
    calls = []

    def recompute(*, retain_polynomials=False):
        calls.append(retain_polynomials)
        return deepcopy(replay)

    monkeypatch.setattr(history, "run", recompute)
    monkeypatch.setattr(sys, "argv", ["infinite_history", *map(str, args)])
    history.main()
    return calls


def test_cli_creates_checks_and_refuses_overwrite(tmp_path, monkeypatch, capsys, compact_replay):
    path = tmp_path/"receipt.json"
    assert _cli(monkeypatch, ["--write", path], compact_replay) == [False]
    original = path.read_bytes()
    assert load(path) == compact_replay
    assert json.loads(capsys.readouterr().out)["all_pass"] is True
    assert _cli(monkeypatch, ["--check", path], compact_replay) == [False]
    summary = json.loads(capsys.readouterr().out)
    assert summary["contrast_interval"] == compact_replay["fine"]["infinite_feedback_contrast_interval"]
    assert summary["physical_clock_calibration_or_independent_data"] is False
    with pytest.raises(FileExistsError):
        _cli(monkeypatch, ["--write", path], compact_replay)
    assert path.read_bytes() == original
    assert capsys.readouterr().out == ""


def test_cli_polynomial_export_requires_matching_check_mode(
        tmp_path, monkeypatch, expanded_replay, compact_replay, capsys):
    path = tmp_path/"expanded.json"
    assert _cli(monkeypatch, ["--polynomials", "--write", path], expanded_replay) == [True]
    capsys.readouterr()
    assert _cli(monkeypatch, ["--polynomials", "--check", path], expanded_replay) == [True]
    capsys.readouterr()
    with pytest.raises(ValueError):
        _cli(monkeypatch, ["--check", path], compact_replay)
    assert capsys.readouterr().out == ""


def test_cli_duplicate_key_fails_before_expensive_recomputation(tmp_path, monkeypatch):
    path = tmp_path/"duplicate.json"
    path.write_text('{"all_pass":false,"all_pass":true}', encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["infinite_history", "--check", str(path)])

    def forbidden(**kwargs):
        raise AssertionError("malformed JSON must fail before recomputation")

    monkeypatch.setattr(history, "run", forbidden)
    with pytest.raises(ValueError, match="duplicate"):
        history.main()


@pytest.mark.parametrize("constant", ("NaN", "Infinity", "-Infinity"))
def test_cli_never_accepts_nonfinite_receipt_fields(
        tmp_path, monkeypatch, compact_replay, capsys, constant):
    path = tmp_path/"nonfinite.json"
    text = history.canonical(compact_replay).replace('"all_pass": true', f'"all_pass": {constant}', 1)
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        _cli(monkeypatch, ["--check", path], compact_replay)
    assert capsys.readouterr().out == ""


def test_cli_write_and_check_are_mutually_exclusive(tmp_path, monkeypatch):
    path = tmp_path/"receipt.json"
    monkeypatch.setattr(sys, "argv", ["infinite_history", "--write", str(path), "--check", str(path)])
    with pytest.raises(SystemExit) as error:
        history.main()
    assert error.value.code == 2
    assert not path.exists()

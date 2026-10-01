"""False-green controls for finite menu accounting and probability claims."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import sys

import pytest

from selection_accounting import DEFAULT_OUT, PARENTS, ROOT, SCORECARD, TRUNK, build, main
from verify_selection_accounting import verify


def test_canonical_receipt_reproduces_and_independent_verifier_accepts():
    stored = json.loads(DEFAULT_OUT.read_text(encoding="utf-8"))
    assert build() == stored
    verify(stored)


def test_accounting_is_independent_of_windows_ansi_text_default(monkeypatch, tmp_path):
    original_read = Path.read_text

    def ansi_default(path, encoding=None, errors=None):
        # Intercept before pathlib/io resolve a missing encoding using the
        # host's UTF-8 mode, which need not match a Windows ANSI environment.
        return original_read(path, encoding=encoding or "cp1252", errors=errors)

    monkeypatch.setattr(Path, "read_text", ansi_default)
    # Reproduce the actual failure mode: this UTF-8 scientific text contains
    # bytes that cannot be decoded using the Windows cp1252 default.
    with pytest.raises(UnicodeDecodeError):
        (ROOT / SCORECARD).read_text()
    expected = json.loads(DEFAULT_OUT.read_text(encoding="utf-8"))
    assert build() == expected
    verify(expected)
    output = tmp_path / "selection.json"
    monkeypatch.setattr(sys, "argv", ["selection_accounting.py", "--output", str(output)])
    main()
    assert output.read_bytes() == DEFAULT_OUT.read_bytes()


@pytest.mark.parametrize("mutation", [
    "missing_exponent", "wrong_minimizer", "trunk_certified", "hidden_pair_defect",
    "chance_probability", "probability_from_count", "alias_overcount",
    "invented_hit", "omitted_tail", "same_probability_countermodels", "bad_probability_measure",
])
def test_false_green_rewrites_are_rejected(mutation):
    data = deepcopy(json.loads(DEFAULT_OUT.read_text(encoding="utf-8")))
    if mutation == "missing_exponent":
        data["stage5_search"]["electron_candidates"].pop()
    elif mutation == "wrong_minimizer":
        data["stage5_search"]["selected_vector"][0] = 8
    elif mutation == "trunk_certified":
        data["historical_outputs"][2]["support"] = "interval_certified_fixed_point"
    elif mutation == "hidden_pair_defect":
        data["historical_outputs"][2]["relative_printed_pair_defect"] = "0"
    elif mutation == "chance_probability":
        data["interpretation"]["global_chance_probability"] = "7e-5"
    elif mutation == "probability_from_count":
        data["interpretation"]["counts_are_probabilities"] = True
    elif mutation == "alias_overcount":
        data["substitution_accounting"]["distinct_numeric_pair_count"] = 48
    elif mutation == "invented_hit":
        data["substitution_accounting"]["certified_alternative_hits"] = 1
    elif mutation == "omitted_tail":
        data["historical_outputs"][0]["tail_bounds_included"] = False
    elif mutation == "same_probability_countermodels":
        data["probability_nonidentification"]["probability_models"][1] = deepcopy(data["probability_nonidentification"]["probability_models"][0])
    else:
        data["probability_nonidentification"]["probability_models"][1]["weights"] = ["1/2", "3/4"]
    with pytest.raises(ValueError):
        verify(data)


@pytest.mark.parametrize("mutation", ["phi", "sqrt_pi", "zero_residual", "physical_prediction"])
def test_rehashed_parent_corruption_is_rejected(tmp_path, mutation):
    for relative in PARENTS:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    path = tmp_path / TRUNK
    trunk = json.loads(path.read_text(encoding="utf-8"))
    if mutation in ("phi", "sqrt_pi"):
        trunk["closed_form_candidate"][mutation] = "1.11111111111111111111111111111"
    elif mutation == "zero_residual":
        trunk["fixed_point_candidate"]["alpha_fixed_point_residual"] = "0"
    else:
        trunk["consumer_policy"]["may_feed_live_particle_predictions"] = True
    path.write_text(json.dumps(trunk), encoding="utf-8", newline="\n")
    # Rebuilding refreshes custody hashes, so rejection requires semantics.
    with pytest.raises(ValueError):
        verify(build(tmp_path), tmp_path)

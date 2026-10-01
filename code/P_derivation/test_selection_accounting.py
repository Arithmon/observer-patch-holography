"""False-green controls for finite menu accounting and probability claims."""
from copy import deepcopy
import json
import shutil

import pytest

from selection_accounting import DEFAULT_OUT, PARENTS, ROOT, TRUNK, build
from verify_selection_accounting import verify


def test_canonical_receipt_reproduces_and_independent_verifier_accepts():
    stored = json.loads(DEFAULT_OUT.read_text())
    assert build() == stored
    verify(stored)


@pytest.mark.parametrize("mutation", [
    "missing_exponent", "wrong_minimizer", "trunk_certified", "hidden_pair_defect",
    "chance_probability", "probability_from_count", "alias_overcount",
    "invented_hit", "omitted_tail", "same_probability_countermodels", "bad_probability_measure",
])
def test_false_green_rewrites_are_rejected(mutation):
    data = deepcopy(json.loads(DEFAULT_OUT.read_text()))
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
    trunk = json.loads(path.read_text())
    if mutation in ("phi", "sqrt_pi"):
        trunk["closed_form_candidate"][mutation] = "1.11111111111111111111111111111"
    elif mutation == "zero_residual":
        trunk["fixed_point_candidate"]["alpha_fixed_point_residual"] = "0"
    else:
        trunk["consumer_policy"]["may_feed_live_particle_predictions"] = True
    path.write_text(json.dumps(trunk))
    # Rebuilding refreshes custody hashes, so rejection requires semantics.
    with pytest.raises(ValueError):
        verify(build(tmp_path), tmp_path)

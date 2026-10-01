"""Reject plausible false-positive reconciliation rewrites."""
import copy
import json
import runpy
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
CHECK = runpy.run_path(str(HERE / "verify.py"))["verify"]


@pytest.fixture
def receipt():
    return json.loads((HERE / "receipt.json").read_text())


def test_registered_reconciliation(receipt):
    assert CHECK(receipt, HERE / "refinement.csv") == {"regions": 96, "trends": 352, "clock_outcomes": 120}


@pytest.mark.parametrize("mutation", [
    "drop_region", "duplicate_region", "drop_control", "erase_marginal_error", "invent_mean_error",
    "flip_deviation", "drop_reference", "change_cdf", "promote_trend", "discard_failed_clock",
    "invent_covariance", "promote_continuum", "recover_missing_driver", "rewrite_input_pin",
    "promote_clock_scope", "add_scope_claim", "add_interpretation_claim", "change_future_contract",
])
def test_corruption_fails(receipt, mutation):
    bad = copy.deepcopy(receipt)
    region = bad["regions"][0]
    if mutation == "drop_region":
        bad["regions"].pop()
    elif mutation == "duplicate_region":
        bad["regions"][1] = copy.deepcopy(region)
    elif mutation == "drop_control":
        bad["regions"] = [r for r in bad["regions"] if r["spatial_dimension"] != 1]
    elif mutation == "erase_marginal_error":
        region["statistics"]["C4_coefficient"]["marginal_standard_error"] = 0.0
    elif mutation == "invent_mean_error":
        region["statistics"]["interval_mean"]["marginal_standard_error"] = 0.0001
    elif mutation == "flip_deviation":
        region["statistics"]["interval_mean"]["signed_deviation"] *= -1
    elif mutation == "drop_reference":
        del region["cdf_signed_deviations_by_reference_dimension"]["6"]
    elif mutation == "change_cdf":
        region["cdf_signed_deviations_by_reference_dimension"]["4"][2] += 0.01
    elif mutation == "promote_trend":
        next(t for t in bad["trends"] if not t["strictly_shrinks_at_both_steps"])["strictly_shrinks_at_both_steps"] = True
    elif mutation == "discard_failed_clock":
        bad["clock_outcomes"] = [r for r in bad["clock_outcomes"] if r["status"] == "timelike_and_graph_related"]
    elif mutation == "invent_covariance":
        bad["interpretation"]["joint_chain_covariance"] = "reconstructed"
    elif mutation == "promote_continuum":
        bad["interpretation"]["continuum_manifoldlikeness_established"] = True
    elif mutation == "recover_missing_driver":
        bad["interpretation"]["q144_chunk_parallel_driver"] = "retained"
    elif mutation == "promote_clock_scope":
        bad["clock_boundary"] = "These are independently executed q144 samples."
    elif mutation == "add_scope_claim":
        bad["continuum_proved"] = True
    elif mutation == "add_interpretation_claim":
        bad["interpretation"]["covariance_recovered"] = True
    elif mutation == "change_future_contract":
        bad["interpretation"]["future_uncertainty_contract"]["weighted_spectrum"] = "The existing receipt suffices for joint significance."
    else:
        bad["inputs"][next(iter(bad["inputs"]))] = "0" * 64
    with pytest.raises((ValueError, KeyError)):
        CHECK(bad)


def test_table_cannot_silently_omit_rows(receipt, tmp_path):
    path = tmp_path / "truncated.csv"
    path.write_text("\n".join((HERE / "refinement.csv").read_text().splitlines()[:-1]) + "\n")
    with pytest.raises(ValueError, match="complete CSV"):
        CHECK(receipt, path)

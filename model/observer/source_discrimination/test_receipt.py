"""Retained receipts must reject changes to evidence and scientific status."""
from copy import deepcopy
from pathlib import Path

import pytest

from tools.strict_json import load, loads
from . import experiment


@pytest.fixture(scope="module")
def replay():
    return experiment.run()


def test_retained_source_receipt(replay):
    experiment.compare_receipt(load(Path(__file__).with_name("results.json")), replay)


@pytest.mark.parametrize("mutation", ("empty", "status", "nan", "count", "resource", "state_error", "physical"))
def test_receipt_rejects_changed_evidence(replay, mutation):
    changed = deepcopy(replay)
    if mutation == "empty":
        changed = {}
    elif mutation == "status":
        changed["all_pass"] = False
    elif mutation == "nan":
        changed["candidates"]["C"]["executions"][0]["plus_probability"] = float("nan")
    elif mutation == "count":
        changed["candidates"]["E"]["counts"]["settings"][0]["events"][0]["count"] += 1
    elif mutation == "resource":
        changed["candidates"]["C"]["executions"][0]["resources"]["proper_code_load_events"] -= 1
    elif mutation == "state_error":
        changed["candidates"]["C"]["executions"][0]["maximum_branch_density_frobenius_error"] = .001
    else:
        changed["illustrative_required_error_allocation"]["independently_established"] = True
    with pytest.raises(ValueError):
        experiment.compare_receipt(changed, replay)


def test_receipt_input_rejects_duplicate_keys():
    with pytest.raises(ValueError, match="duplicate"):
        loads('{"all_pass": false, "all_pass": true}')

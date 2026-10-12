"""Retained receipts must reject changes to evidence and scientific status."""
from copy import deepcopy
from pathlib import Path
import random

import numpy as np
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


def test_uniform_sampler_is_stable_across_the_binomial_half_boundary():
    # The former multinomial sampler maps a one-ulp sign change around 1/2
    # into a complemented draw (4019 versus 4173 at this retained seed).
    epsilon = np.spacing(.5)
    below = np.array([.5 - epsilon, .5 + epsilon])
    above = below[::-1]
    old_below = np.random.default_rng(experiment.SEED).multinomial(8192, below)
    old_above = np.random.default_rng(experiment.SEED).multinomial(8192, above)
    assert not np.array_equal(old_below, old_above)
    left = experiment.sample_counts(below, 8192, random.Random(experiment.SEED))
    right = experiment.sample_counts(above, 8192, random.Random(experiment.SEED))
    assert np.array_equal(left, right)
    assert sum(left) == 8192


def test_retained_source_sampler_tolerates_small_execution_roundoff(replay):
    baseline_rng = random.Random(experiment.SEED)
    perturbed_rng = random.Random(experiment.SEED)
    jitter = np.random.default_rng(34)
    for candidate in ("C", "E"):
        for execution in replay["candidates"][candidate]["executions"]:
            probabilities = np.array([row["probability"] for row in execution["branches"]])
            perturbed = probabilities * (1 + 1e-13 * jitter.standard_normal(len(probabilities)))
            expected = np.array([row["count"] for row in execution["branches"]])
            actual = experiment.sample_counts(probabilities, 8192, baseline_rng)
            changed = experiment.sample_counts(perturbed, 8192, perturbed_rng)
            assert np.array_equal(actual, expected)
            assert np.array_equal(changed, expected)


def test_sampler_uses_supplied_execution_probabilities_and_retains_zero_bins():
    assert experiment.sample_counts([0., 2., 0.], 128, random.Random(7)).tolist() == [0, 128, 0]
    assert experiment.sample_counts([2., 0., 0.], 128, random.Random(7)).tolist() == [128, 0, 0]

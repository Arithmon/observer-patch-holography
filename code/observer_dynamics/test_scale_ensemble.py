"""Independent controls for the retained scale and angular estimator population.

These checks start from the pinned integer mass or original small vectors.
They do not claim to reconstruct the unavailable large terminal arrays.
"""
from fractions import Fraction
from collections import Counter
from itertools import combinations
import json
from math import comb
from statistics import mean

import numpy as np
import pytest

from .archive_adapter import SIM
from . import retained_measurements as retained


@pytest.fixture
def population():
    return json.loads((SIM / "data/codex_audit_20260925/scale_ensemble.json").read_text())


def reject(population):
    with pytest.raises(ValueError):
        retained.scale_ensemble(population)


def angular_first(population):
    return next(iter(population["levels"][0]["angular"].values()))["rows"][0]


def test_complete_frozen_population_and_independent_parent_nulls(population):
    result = retained.scale_ensemble(population)
    assert result["terminal_digest_crosspins"] == 48
    assert result["scale_rows"] == 51
    assert result["raw_terminal_arrays_rehashed"] is False
    checked = 0
    for level in population["levels"]:
        parent = json.loads((SIM / (
            "data/codex_audit_20260925/inputs/data/receipts/"
            f"L{level['level']}/receipt.json"
        )).read_text())
        ports = parent["ports"]
        raised = parent["expected"]["component_total"][0] % ports
        for scale in level["scales"]:
            if not scale["name"].startswith("ancestor_"):
                continue
            group_ports = 12 * int(scale["name"].split("_")[1])
            # Scalar hypergeometric variance of an exhaustive group's mean.
            # The population-wide mean is fixed by the original integer mass.
            expected = Fraction(
                raised * (ports - raised) * (ports - group_ports),
                ports * ports * group_ports * (ports - 1),
            )
            assert scale["equilibrium_variance"] == pytest.approx(float(expected), rel=1e-13, abs=0)
            checked += 1
    assert checked == 27


def test_coherent_null_rescaling_cannot_reverse_enhancement(population):
    scale = population["levels"][0]["scales"][0]
    assert mean(row["variance_over_exchangeable_equilibrium"] for row in scale["per_schedule"]) > 1
    scale["equilibrium_variance"] *= 2
    for row in scale["per_schedule"]:
        row["variance_over_exchangeable_equilibrium"] /= 2
    assert mean(row["variance_over_exchangeable_equilibrium"] for row in scale["per_schedule"]) < 1
    reject(population)


def test_differing_responses_require_positive_schedule_noise(population):
    scale = population["levels"][0]["scales"][0]
    responses = [Fraction.from_float(row["response"]) for row in scale["per_schedule"]]
    center = sum(responses) / len(responses)
    initial = Fraction.from_float(scale["per_schedule"][0]["initial_variance"])
    # Projection onto the fixed initial vector cannot increase variance.
    minimum = initial * sum((value - center) ** 2 for value in responses) / (len(responses) - 1)
    assert float(minimum) == pytest.approx(7.028510514184177e-6, rel=1e-13)
    scale["conditional_schedule_variance_unbiased"] = 0.0
    scale["cross_schedule_common_power_unbiased"] = scale["total_terminal_variance"]
    scale["schedule_noise_fraction"] = 0.0
    reject(population)


@pytest.mark.parametrize("field", ["initial_variance_over_iid_null", "groups", "cells_min", "cells_max"])
def test_ancestral_geometry_and_iid_normalization_are_bound_to_declared_groups(population, field):
    population["levels"][0]["scales"][0][field] *= 2
    reject(population)


def test_reported_coherence_must_agree_with_response_and_power(population):
    population["levels"][0]["scales"][0]["per_schedule"][0]["coherence_squared"] *= 0.5
    reject(population)


def test_schedules_share_one_initial_field(population):
    row = population["levels"][0]["scales"][0]["per_schedule"][0]
    row["initial_variance"] *= 1.1
    row["power_ratio"] = row["terminal_variance"] / row["initial_variance"]
    row["coherence_squared"] = row["response"] ** 2 / row["power_ratio"]
    row["residual_power_ratio"] = row["power_ratio"] - row["response"] ** 2
    reject(population)


@pytest.mark.parametrize("field", ["ports", "carriers", "balanced_fraction", "initial_load_seed"])
def test_original_parent_population_metadata_cannot_change(population, field):
    population["levels"][0][field] *= 2
    reject(population)


def test_geometry_digest_matches_parent_receipt(population):
    population["levels"][0]["geometry_file_sha256"] = "0" * 64
    reject(population)


@pytest.mark.parametrize("change", ["remove", "duplicate", "empty"])
def test_specification_requires_the_complete_scale_population(population, change):
    scales = population["levels"][0]["scales"]
    if change == "remove":
        scales.pop()
    elif change == "duplicate":
        scales[-1] = scales[0]
    else:
        scales.clear()
    reject(population)


def test_changed_embedded_specification_cannot_excuse_missing_scale(population):
    population["spec"]["primary_ancestral_group_cells"].remove(65536)
    for level in population["levels"]:
        level["scales"] = [scale for scale in level["scales"] if scale["name"] != "ancestor_65536"]
    reject(population)


@pytest.mark.parametrize("change", ["grid", "row", "duplicate_ell", "empty"])
def test_specification_requires_all_angular_grids_and_multipoles(population, change):
    angular = population["levels"][0]["angular"]
    first = next(iter(angular))
    if change == "grid":
        del angular[first]
    elif change == "row":
        angular[first]["rows"].pop()
    elif change == "duplicate_ell":
        angular[first]["rows"][1]["ell"] = 1
    else:
        angular.clear()
    reject(population)


@pytest.mark.parametrize("field", [
    "C_terminal_mean", "cross_response_mean", "cross_response_schedule_sd",
    "total_power_ratio_mean", "coherence_squared_mean", "orthogonal_residual_power_ratio_mean",
])
def test_angular_observables_agree_with_retained_schedule_population(population, field):
    angular_first(population)[field] *= 2
    reject(population)


def test_angular_differing_responses_cannot_have_zero_schedule_noise(population):
    row = angular_first(population)
    assert len(set(row["per_schedule_response"])) > 1
    row["conditional_schedule_noise_power"] = 0.0
    row["cross_schedule_common_power"] = row["C_terminal_mean"]
    reject(population)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), True])
def test_nonphysical_coherence_scalars_are_rejected(population, bad):
    population["levels"][0]["scales"][0]["per_schedule"][0]["coherence_squared"] = bad
    reject(population)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), True])
def test_nonphysical_angular_scalars_are_rejected(population, bad):
    angular_first(population)["coherence_squared_mean"] = bad
    reject(population)


def test_reordered_schedule_entries_cannot_break_digest_alignment(population):
    rows = population["levels"][0]["scales"][0]["per_schedule"]
    rows[0], rows[1] = rows[1], rows[0]
    reject(population)


def test_reordered_terminal_inputs_cannot_break_parent_alignment(population):
    rows = population["levels"][0]["inputs"]
    rows[0], rows[1] = rows[1], rows[0]
    reject(population)



"""Tests for the compare-only criticality distance surface."""

from __future__ import annotations

import copy
import functools
import json

import pytest

import derive_d11_criticality_comparison as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_row_class_is_compare_only():
    artifact = _artifact()
    assert artifact["row_class"] == "compare_only_never_a_prediction_ancestor"
    assert artifact["promotion_allowed"] is False
    assert artifact["physical_fit_evaluable"] is False


def test_pole_comparison_does_not_consume_direct_monte_carlo_average():
    scan = json.loads(lane.SCAN.read_text(encoding="utf-8"))
    references = json.loads(lane.REFERENCE_JSON.read_text(encoding="utf-8"))
    expected = lane.build_artifact(scan, references)
    assert expected["references"]["mt_source_id"] == "Q007TP4"
    assert expected["references"]["mt_gev"] == 172.1
    changed = copy.deepcopy(references)
    changed["entries"]["top_quark_direct_aux"]["value_gev"] = 190.0
    assert lane.build_artifact(scan, changed) == expected
    changed["entries"]["top_quark"]["source"]["summary_id"] = "Q007TP"
    with pytest.raises(ValueError, match="cross-section pole-mass"):
        lane.build_artifact(scan, changed)


def test_archived_boundary_distances_match_known_deficits():
    archived = _artifact()["one_loop_named_distances"]["mu_U_gauge_unification"]
    assert abs(archived["mt_relative"] + 0.049) < 0.005
    assert abs(archived["mh_relative"] + 0.080) < 0.005


def test_single_anchor_boundary_improves_both_channels():
    named = _artifact()["one_loop_named_distances"]
    archived = named["mu_U_gauge_unification"]
    single = named["E_cell_pixel_energy"]
    assert abs(single["mt_relative"]) < abs(archived["mt_relative"])
    assert abs(single["mh_relative"]) < abs(archived["mh_relative"])
    assert abs(single["mt_relative"]) < 0.02
    assert abs(single["mh_relative"]) < 0.03


def test_curve_nearest_point_carries_its_boundary_scale():
    best = _artifact()["one_loop_curve_nearest_point"]
    assert 1.0e16 < best["boundary_scale_gev"] < 1.4e19
    assert abs(best["mt_relative"]) < 0.03
    assert abs(best["mh_relative"]) < 0.03

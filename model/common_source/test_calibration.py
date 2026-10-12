"""Synthetic sampling, exact dimensionless distances and energy-cutoff checks."""
import json
import math

import numpy as np
import pytest

from . import calibration as c


@pytest.fixture(scope="module")
def receipt():
    return c.report(steps=16)


def test_all_clock_settings_are_sampled_and_interval_is_simultaneous(receipt):
    clock = receipt["clock_sampling"]
    assert len(clock["settings"]) == 5
    assert clock["shots_per_setting"] == 100000
    expected_radius = math.sqrt(math.log(2 * 5 / .05) / (2 * 100000))
    for row in clock["settings"]:
        expected_probability = (1 + clock["state_derived_normalized_contrast"] * math.sin(row["programmed_modular_phase"])) / 2
        assert row["normalized_source_probability"] == pytest.approx(expected_probability, abs=2e-13)
        assert row["hoeffding_probability_radius"] == pytest.approx(expected_radius)
        assert row["observed_frequency"] == row["number_one_count"] / 100000
        assert row["shots"] == 100000
        assert row["programmed_parameter_inside_interval"]
    assert clock["all_programmed_parameters_inside_reported_intervals"]
    assert clock["sampling_is_synthetic"]
    assert clock["fresh_complete_history_per_shot_and_setting"]
    assert not clock["sequential_use_of_one_system"]
    assert not clock["physical_seconds_calibrated"]
    assert not clock["confidence_covers_hardware_or_state_estimation_error"]


def test_seeded_clock_counts_are_reproducible_without_refitting_the_clock(receipt):
    first = receipt["clock_sampling"]
    second = c.clock_benchmark(16)
    assert [row["number_one_count"] for row in first["settings"]] == [row["number_one_count"] for row in second["settings"]]
    assert first["state_derived_modular_frequency"] == second["state_derived_modular_frequency"]


def test_sampling_keeps_source_roundoff_without_flipping_the_half_probability_count(monkeypatch):
    clock = {"contrast": .25, "mass": 1., "modular_frequency": .5,
             "monotone_parameter_window": [0., math.pi]}
    emitted = {"clock_parameters": clock, "all_pass": True, "settings": []}
    monkeypatch.setattr(c.experiment, "history", lambda *args: None)
    monkeypatch.setattr(c.modular_clock, "receipt", lambda _: emitted)
    probabilities = (math.nextafter(.5, 0.), math.nextafter(.5, 1.))
    results = []
    for probability in probabilities:
        emitted["settings"] = [{"clock_bit_marginal": [1-probability, probability],
                                "modular_parameter": 0., "phase": 0.,
                                "analytic_probability_error": abs(probability-.5)}]
        results.append(c.clock_benchmark()["settings"][0])
    assert [row["normalized_source_probability"] for row in results] == list(probabilities)
    assert results[0]["number_one_count"] == results[1]["number_one_count"]
    assert min(row["minimum_sampling_threshold_distance"] for row in results) > 1e-10
    # Reconstruct uniforms directly from the independent raw integer stream.
    raw = np.random.PCG64(20261012).random_raw(100000)
    expected = sum(int(value) < 2**63 for value in raw)
    assert results[0]["number_one_count"] == expected


def test_sampling_counts_actual_threshold_crossings_without_rounding_probabilities():
    class FixedUniforms:
        def random(self, shots):
            assert shots == 3
            return np.array([.25, .5, .75])

    below = math.nextafter(.5, 0.)
    above = math.nextafter(.5, 1.)
    count_below, margin_below = c._sample_clock_count(FixedUniforms(), 3, below)
    count_above, margin_above = c._sample_clock_count(FixedUniforms(), 3, above)
    assert (count_below, count_above) == (1, 2)
    assert margin_below == .5-below
    assert margin_above == above-.5


def test_interval_inversion_handles_boundaries_and_empty_model_intersections():
    clock = {"contrast": .25, "mass": 1., "modular_frequency": .5,
             "monotone_parameter_window": [0., math.pi]}
    lower = c.binomial_parameter_interval(49500, 100000, clock)
    assert lower["point_estimate_projected_to_declared_window"]
    assert lower["inferred_modular_parameter"] == 0
    assert lower["modular_parameter_interval"][0] == 0
    upper = c.binomial_parameter_interval(62800, 100000, clock)
    assert upper["point_estimate_projected_to_declared_window"]
    assert upper["modular_parameter_interval"][1] == pytest.approx(math.pi)
    invalid = c.binomial_parameter_interval(0, 100000, clock)
    assert invalid["modular_parameter_interval"] is None
    assert not invalid["interval_intersects_declared_model_range"]


def test_probability_interval_inversion_contains_interior_values_with_exact_frequencies():
    clock = {"contrast": .25, "mass": 1., "modular_frequency": .5,
             "monotone_parameter_window": [0., math.pi]}
    for count in (50000, 53000, 55000, 57000, 60000, 62500):
        row = c.binomial_parameter_interval(count, 100000, clock)
        point = math.asin((2 * count / 100000 - 1) / .25) / .5
        assert row["modular_parameter_interval"][0] <= point <= row["modular_parameter_interval"][1]


@pytest.mark.parametrize("window", ([0., 2 * math.pi / 3], [.25, 2 * math.pi / 3]))
def test_narrower_declared_window_bounds_every_inference_including_high_counts(window):
    clock = {"contrast": .25, "mass": 1., "modular_frequency": .5,
             "monotone_parameter_window": window}
    for count in (0, 49500, 50000, 53000, 60000, 60800, 61000, 62500, 62800, 100000):
        row = c.binomial_parameter_interval(count, 100000, clock)
        assert window[0] <= row["inferred_modular_parameter"] <= window[1]
        interval = row["modular_parameter_interval"]
        if interval is not None:
            assert window[0] <= interval[0] <= interval[1] <= window[1]
    top = c.binomial_parameter_interval(61000, 100000, clock)
    assert top["modular_parameter_interval"][1] == window[1]
    assert top["point_estimate_projected_to_declared_window"]
    assert c.binomial_parameter_interval(62500, 100000, clock)["modular_parameter_interval"] is None


@pytest.mark.parametrize("window", (None, [0., 4.], [-1., 1.], [1., 1.], [0., float("inf")]))
def test_missing_or_nonmonotone_parameter_window_is_rejected(window):
    clock = {"contrast": .25, "mass": 1., "modular_frequency": .5,
             "monotone_parameter_window": window}
    with pytest.raises(ValueError):
        c.binomial_parameter_interval(55000, 100000, clock)


def test_exact_response_distances_match_icosahedron_and_native_current_displacements(receipt):
    geometry = receipt["response_distances"]
    assert geometry["exact_pair_comparisons"] == 66
    assert geometry["native_edge_count"] == 30
    assert geometry["native_current_displacements_checked_exactly"] == 90
    assert geometry["maximum_numerical_coordinate_distance_error"] < 2e-15
    expected_edge = math.sqrt(2 - 2 / math.sqrt(5))
    assert geometry["expected_unit_radius_edge_length"] == pytest.approx(expected_edge, abs=2e-15)
    edges = [row for row in geometry["pairs"] if row["native_edge"]]
    assert len(edges) == 30
    assert all(abs(row["response_distance"] - expected_edge) < 2e-15 for row in edges)
    distances = [row["response_distance"] for row in geometry["pairs"]]
    assert sum(abs(value - 2) < 2e-15 for value in distances) == 6
    assert sum(abs(value - math.sqrt(2 + 2 / math.sqrt(5))) < 2e-15 for value in distances) == 30
    assert not geometry["physical_length_calibrated"]
    assert not geometry["three_dimensional_spatial_continuum_established"]


def test_smooth_vacuum_energy_cutoffs_converge_without_labeling_matrix_target_as_execution(receipt):
    energy = receipt["smooth_energy_cutoff"]
    rows = energy["rows"]
    assert [row["dimension"] for row in rows] == [2, 4, 8, 16, 32]
    assert rows[0]["unconditional_energy_increment"] == pytest.approx(math.sin(1 / (2 * math.sqrt(2)))**2, abs=2e-15)
    assert rows[0]["absolute_error_from_infinite_increment"] > rows[1]["absolute_error_from_infinite_increment"] > rows[2]["absolute_error_from_infinite_increment"]
    for row in rows[2:]:
        assert row["absolute_error_from_infinite_increment"] < 2e-12
    for row in rows:
        assert row["outcome_probabilities"] == pytest.approx([.5, .5], abs=2e-14)
        assert row["all_outcomes_retained"]
        assert not row["source_compiler_execution"]
    assert not energy["common_interacting_vacuum_used_in_this_component_benchmark"]
    assert not energy["physical_energy_calibrated"]


def test_report_does_not_promote_sampling_coverage_to_a_guaranteed_pass(receipt):
    json.dumps(receipt, allow_nan=False)
    assert receipt["all_pass"]
    assert not any("inside" in key or "coverage" in key for key in receipt["checks"])


@pytest.mark.parametrize("count,shots", ((True, 100), (-1, 100), (101, 100), (1, 0), (1, 100.)))
def test_invalid_binomial_inputs_rejected(count, shots):
    with pytest.raises(ValueError):
        c.binomial_parameter_interval(count, shots, {"contrast": .2, "mass": 1., "modular_frequency": .5})

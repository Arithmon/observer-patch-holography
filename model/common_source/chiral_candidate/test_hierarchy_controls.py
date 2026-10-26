from fractions import Fraction as F

from . import hierarchy_controls as h


def test_fixed_source_textures_are_certified_and_fit_no_parameters():
    prediction = h.predictions()
    assert prediction["empirical_mass_inputs"] == []
    assert F(prediction["linear_F_ratios"][-1]["upper_exact"]) < F(3, 2)
    for row in prediction["nonlinear_controls"]:
        assert not row["parameters_fitted_to_masses"]
        gaps = row["holon_thresholds_over_clock"]
        assert all(F(item["lower_exact"]) > 0 for item in gaps)
        assert all(F(a["upper_exact"]) < F(b["lower_exact"]) for a, b in zip(gaps, gaps[1:]))


def test_two_independent_ratios_are_required_and_every_fixed_control_is_reported():
    comparison = h.lepton_comparison(h.predictions())
    assert len(comparison["comparisons"]) == 4
    assert not any(row["both_ratios_pass"] for row in comparison["comparisons"])
    assert not comparison["an_overall_mass_scale_can_change_these_ratios"]
    assert not comparison["other_nonlinear_or_interacting_actions_excluded"]

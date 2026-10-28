"""Fixed linear and nonlinear source textures, without mass-target fitting.

The Hubbard function is a declared NUMERICAL texture prescription here.
Its one-dimensional theorem is not an interaction-generation theorem for
this three-dimensional chiral candidate. No observed mass enters these
prediction functions. Comparisons belong in a separate benchmark report.
"""
from fractions import Fraction as F
from functools import lru_cache

from ..mass_response import gap_certificate as g
from ..particle_identification import family_attachment


@lru_cache(maxsize=1)
def predictions():
    _, _, _, _, eigenvalues, *_ = family_attachment._data()
    five = g.sqrt_bounds(g.RationalInterval.enclose(5))
    features = [g.RationalInterval.enclose(a)+b*five for a, b in eigenvalues]
    source = g.native_reference_data()
    w0, w1, w2 = source["weights"]
    scales = {"a": g.log_bounds(w1/w0), "b": g.log_bounds(w2/w1),
              "J": g.log_bounds(w1*w1/(w0*w2))}
    omega = g.log_bounds(w2/w0)
    rows = []
    for name, scale in scales.items():
        gaps = [g.hubbard_charge_gap_bounds(omega/2, scale*feature)/2 for feature in features]
        if not all(value.lower > 0 for value in gaps):
            raise ArithmeticError("source threshold too small for the declared interval precision")
        if not all(a.upper < b.lower for a, b in zip(gaps, gaps[1:])):
            raise ArithmeticError("the three source thresholds are not resolved in order")
        rows.append({"repulsion_coefficient": name, "U_assignment": name+"*eigenvalue(F)",
                     "holon_thresholds_over_clock": [(gap/omega).as_json() for gap in gaps],
                     "ascending_ratios": [(gap/gaps[0]).as_json() for gap in gaps],
                     "same_function_for_all_three_families": True,
                     "parameters_fitted_to_masses": False})
    return {"schema": "oph.chiral_candidate.hierarchy_controls.v1",
            "linear_F_ratios": [(feature/features[0]).as_json() for feature in features],
            "nonlinear_controls": rows,
            "empirical_mass_inputs": [],
            "scope": {"three_fixed_source_prescriptions_tested": True,
                      "continuous_parameter_search_performed": False,
                      "Hubbard_function_is_declared_texture_prescription": True,
                      "three_spatial_dimensional_dynamical_gap_theorem": False,
                      "these_tiny_thresholds_measured_in_finite_history": False}}


def lepton_comparison(prediction):
    """Diagnose hierarchy AFTER prediction against broad pole-ratio windows.

    These dimensionless windows are empirical pole-mass benchmark inputs.
    Neither selects a source prescription or changes a predicted value.
    No matching map identifies a source coefficient or declared Hubbard
    texture threshold with a physical lepton pole mass. Window overlap is
    therefore a diagnostic only, not an admissible physical mass comparison.
    """
    windows = ((F(206), F(207)), (F(3470), F(3485)))
    comparisons = []
    prescriptions = [("linear_F", prediction["linear_F_ratios"])]
    prescriptions += [(row["repulsion_coefficient"]+"_Hubbard_function", row["ascending_ratios"])
                      for row in prediction["nonlinear_controls"]]
    for name, ratios in prescriptions:
        overlaps = [not (F(ratio["upper_exact"]) < low or F(ratio["lower_exact"]) > high)
                    for ratio, (low, high) in zip(ratios[1:], windows, strict=True)]
        comparisons.append({"prescription": name, "ratio_intersects_comparison_windows": overlaps,
                            "both_ratios_pass": all(overlaps)})
    return {"benchmark": "charged lepton hierarchy, ascending e/mu/tau masses",
            "source": "https://pdg.lbl.gov/2025/download/db2024.pdf",
            "conservative_ratio_windows": [[str(low), str(high)] for low, high in windows],
            "comparison_contract": {
                "reference_quantity": "charged_lepton_pole_mass_ratios",
                "prediction_quantity": "dimensionless_source_coefficients_or_declared_Hubbard_texture_threshold_ratios",
                "prediction_renormalization_scheme_and_scale": None,
                "source_to_pole_matching_supplied": False,
                "physical_comparison_admissible": False,
                "status": "unmatched_hierarchy_diagnostic_only",
                "required_for_physical_comparison": "derive_mass_and_kinetic_response_together_then_supply_scheme_scale_and_pole_matching",
                "passing_windows_would_establish_physical_mass_agreement": False,
            },
            "empirical_data_used_only_for_comparison": True, "comparisons": comparisons,
            "an_overall_mass_scale_can_change_these_ratios": False,
            "other_nonlinear_or_interacting_actions_excluded": False}

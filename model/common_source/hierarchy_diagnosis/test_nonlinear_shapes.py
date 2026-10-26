"""Exact log enclosures and frozen-prescription shape rejection."""
from copy import deepcopy
from fractions import Fraction as F
import json

import pytest

from ..chiral_candidate import hierarchy_controls
from ..mass_response import gap_certificate as g
from . import nonlinear_shapes as n


def test_range_reduced_log_encloses_independent_high_precision_values():
    import mpmath as mp
    with mp.workdps(100):
        for value in (F(1), F(1, 3), F(3, 2), F(2), F(10**70), F(1, 10**60)):
            interval = n.log_rational_bounds(value)
            target = mp.log(mp.mpf(value.numerator)/value.denominator)
            lower = mp.mpf(interval.lower.numerator)/interval.lower.denominator
            upper = mp.mpf(interval.upper.numerator)/interval.upper.denominator
            assert lower <= target <= upper
            assert interval.upper-interval.lower < F(1, 10**40)


def test_exact_dyadic_log_range_reduction_preserves_powers():
    base = n.log_rational_bounds(2)
    for exponent in (-40, -1, 0, 1, 40):
        value = F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))
        result = n.log_rational_bounds(value)
        expected = exponent*base
        assert result.lower <= expected.lower <= expected.upper <= result.upper


def test_retained_controls_have_disjoint_certified_shapes_without_mutation():
    prediction = hierarchy_controls.predictions()
    windows = [["206", "207"], ["3470", "3485"]]
    before = deepcopy(prediction)
    report = n.report(prediction, windows)
    assert prediction == before
    json.dumps(report, allow_nan=False)
    assert report["all_fixed_controls_have_disjoint_benchmark_shape_exact"]
    observed = report["benchmark_Q_window_enclosure"]
    assert F(observed["lower_exact"]) > F(528, 1000)
    assert F(observed["upper_exact"]) < F(532, 1000)
    for row in report["fixed_source_controls"]:
        shape = row["adjacent_log_spacing_enclosure"]
        assert F(shape["lower_exact"]) > F(108, 100)
        assert F(shape["lower_exact"]) > F(observed["upper_exact"])
        assert not row["three_dimensional_dynamical_gap_derived"]
    correction = report["fixed_source_controls"][0]["necessary_individual_mass_correction_factor"]
    assert F(correction["C_must_be_at_least_lower_endpoint_exact"]) > F(434, 100)
    assert not correction["bound_on_coupling_perturbations"]
    assert not report["new_parameters_fitted_or_scanned"]
    assert not report["source_predictions_changed"]


def test_full_receipt_input_uses_original_exact_endpoints_not_float_displays():
    prediction = deepcopy(hierarchy_controls.predictions())
    for row in prediction["nonlinear_controls"]:
        for interval in row["ascending_ratios"]:
            interval["midpoint"] = -12345.0
    retained = {"source_hierarchy_predictions": prediction,
                "separate_empirical_comparison": hierarchy_controls.lepton_comparison(prediction)}
    report = n.report(retained)
    assert report["all_fixed_controls_have_disjoint_benchmark_shape_exact"]
    assert report["fixed_source_controls"][0]["ascending_ratios"][1] > 200


def test_log_and_shape_input_domains_are_explicit():
    for invalid in (0, -1):
        with pytest.raises(ValueError):
            n.log_rational_bounds(invalid)
    with pytest.raises(ValueError):
        n.logarithm_interval(g.RationalInterval.enclose(0))
    with pytest.raises(ValueError):
        n.original_interval({"midpoint": 2.})
    one, two = (g.RationalInterval.enclose(x).as_json() for x in (1, 2))
    with pytest.raises(ValueError):
        n.adjacent_log_shape(one, two)
    with pytest.raises(ValueError):
        n.benchmark_shape([[1, 2], [3, 4]])

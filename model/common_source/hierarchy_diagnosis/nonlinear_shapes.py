"""Certified log-spacing diagnosis of retained fixed nonlinear prescriptions.

No source prediction is changed or selected. Original rational ratio endpoints
enter range-reduced rational log bounds. Empirical windows are consulted only
for comparison and necessary mass-level correction sizes.
"""
from __future__ import annotations

from fractions import Fraction as F

from ..mass_response import gap_certificate as g


def log_rational_bounds(value):
    """Outward log enclosure, with dyadic range reduction to [1,2).

    The existing positive atanh tail converges rapidly on this range. Applying
    it directly to a huge mass ratio would still bound log, but very loosely.
    """
    value = F(value)
    if value <= 0:
        raise ValueError("a strictly positive exact logarithm argument required")
    exponent = value.numerator.bit_length()-value.denominator.bit_length()
    scale = F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))
    reduced = value/scale
    if reduced < 1:
        exponent -= 1
        reduced *= 2
    if not 1 <= reduced < 2:
        raise ArithmeticError("dyadic logarithm range reduction failed")
    return exponent*g.log_bounds(2)+g.log_bounds(reduced)


def original_interval(record):
    """Read authoritative exact endpoints; displayed floats are never inputs."""
    try:
        lower, upper = F(record["lower_exact"]), F(record["upper_exact"])
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("original exact rational interval endpoints required") from exc
    return g.RationalInterval(lower, upper)


def logarithm_interval(interval):
    if interval.lower <= 0:
        raise ValueError("the logarithm interval must be strictly positive")
    lower = log_rational_bounds(interval.lower)
    upper = log_rational_bounds(interval.upper)
    return g.RationalInterval(lower.lower, upper.upper)


def adjacent_log_shape(ratio_middle, ratio_heaviest):
    """Q=log(m3/m2)/log(m2/m1)=log(r3)/log(r2)-1."""
    r2, r3 = original_interval(ratio_middle), original_interval(ratio_heaviest)
    if r2.lower <= 1 or r3.lower <= r2.upper:
        raise ValueError("strictly resolved ascending ratios 1<r2<r3 required")
    return logarithm_interval(r3)/logarithm_interval(r2)-1


def _windows(windows):
    try:
        (mu_lo, mu_hi), (tau_lo, tau_hi) = tuple(tuple(F(x) for x in row) for row in windows)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("two positive ordered exact ratio windows required") from exc
    if not 1 < mu_lo <= mu_hi < tau_lo <= tau_hi:
        raise ValueError("two positive ordered exact ratio windows required")
    return (mu_lo, mu_hi), (tau_lo, tau_hi)


def benchmark_shape(windows):
    """Monotonic endpoint enclosure using only the declared benchmark windows."""
    (mu_lo, mu_hi), (tau_lo, tau_hi) = _windows(windows)
    log_mu = g.RationalInterval(log_rational_bounds(mu_lo).lower,
                                log_rational_bounds(mu_hi).upper)
    log_tau = g.RationalInterval(log_rational_bounds(tau_lo).lower,
                                 log_rational_bounds(tau_hi).upper)
    return log_tau/log_mu-1


def necessary_mass_correction(heavy_ratio, windows):
    """Necessary C>=1 if each mass correction lies in [1/C,C].

    A free common scale cancels from the ratio. Such corrections can change
    r3 by at most C^2. This constrains masses, not the size of couplings:
    exponentially generated thresholds can be sensitive to small couplings.
    """
    ratio = original_interval(heavy_ratio)
    if ratio.lower <= 0:
        raise ValueError("positive source mass ratios required")
    _, (tau_lo, tau_hi) = _windows(windows)
    square_floor = max(F(1), ratio.lower/tau_hi, tau_lo/ratio.upper)
    threshold = g.sqrt_bounds(g.RationalInterval.enclose(square_floor))
    return {"necessary_threshold_enclosure": threshold.as_json(),
            "C_must_be_at_least_lower_endpoint_exact": str(threshold.lower),
            "scope": "necessary only; each individual mass changes by a factor in[1/C,C], with any common scale allowed",
            "bound_on_coupling_perturbations": False,
            "sufficiency_or_successful_corrected_action_claimed": False}


def report(retained, benchmark_windows=None):
    """Use a retained candidate receipt, or its prediction dict plus windows."""
    if "source_hierarchy_predictions" in retained:
        prediction = retained["source_hierarchy_predictions"]
        comparison = retained["separate_empirical_comparison"]
        windows = comparison["conservative_ratio_windows"] if benchmark_windows is None else benchmark_windows
        old_comparisons = comparison.get("comparisons", [])
        benchmark_source = comparison.get("source")
    else:
        prediction = retained
        if benchmark_windows is None:
            raise ValueError("separate empirical ratio windows required")
        windows, old_comparisons, benchmark_source = benchmark_windows, [], None
    benchmark = benchmark_shape(windows)
    rows = []
    for row in prediction["nonlinear_controls"]:
        ratios = row["ascending_ratios"]
        if len(ratios) != 3:
            raise ValueError("three retained ascending mass ratios required")
        shape = adjacent_log_shape(ratios[1], ratios[2])
        intervals = list(map(original_interval, ratios))
        disjoint = shape.upper < benchmark.lower or benchmark.upper < shape.lower
        rows.append({
            "prescription": row["repulsion_coefficient"]+"_Hubbard_function",
            "original_ratio_intervals": ratios,
            "ascending_ratios": [interval.as_json()["midpoint"] for interval in intervals],
            "next_family_ratio_enclosure": (intervals[2]/intervals[1]).as_json(),
            "next_family_ratio": (intervals[2]/intervals[1]).as_json()["midpoint"],
            "adjacent_log_spacing_enclosure": shape.as_json(),
            "adjacent_log_spacing_diagnostic": shape.as_json()["midpoint"],
            "benchmark_shape_disjoint_exact": disjoint,
            "necessary_individual_mass_correction_factor": necessary_mass_correction(ratios[2], windows),
            "three_dimensional_dynamical_gap_derived": False,
        })
    linear = list(map(original_interval, prediction["linear_F_ratios"]))
    reciprocal_spacing = (1/linear[1]-1/linear[2])/(1/linear[0]-1/linear[1])
    return {
        "schema": "oph.hierarchy_diagnosis.nonlinear_shapes.v1",
        "definition": "Q=log(m3/m2)/log(m2/m1)",
        "fixed_source_controls": rows,
        "benchmark_Q_window_enclosure": benchmark.as_json(),
        "benchmark_Q_window_diagnostic": [float(benchmark.lower), float(benchmark.upper)],
        "benchmark_ratio_windows": [[str(x) for x in row] for row in _windows(windows)],
        "benchmark_source": benchmark_source,
        "all_fixed_controls_have_disjoint_benchmark_shape_exact": all(row["benchmark_shape_disjoint_exact"] for row in rows),
        "log_method": "exact dyadic range reduction, finite positive atanh sum and rigorous geometric remainder; outward rational interval arithmetic",
        "floating_logarithms_are_not_directed_interval_certificates": False,
        "float_displays_authoritative": False,
        "weak_coupling_explanation": {
            "threshold_asymptotic": "Delta/2 ~ (4/pi)*sqrt(t*s*lambda)*exp(-2*pi*t/(s*lambda)), U=s*lambda",
            "log_threshold_asymptotic": "constant + (1/2)log(lambda) - (pi*Omega/s)/lambda, t=Omega/2",
            "reciprocal_spacing_limit_definition": "(1/lambda2-1/lambda3)/(1/lambda1-1/lambda2)",
            "reciprocal_spacing_limit_enclosure": reciprocal_spacing.as_json(),
            "primary_source": "https://arxiv.org/html/cond-mat/9911281",
            "asymptotic_used_for_exact_exclusion": False,
            "exclusion_for_every_possible_repulsion_scale_claimed": False},
        "original_exact_interval_comparisons": old_comparisons,
        "new_parameters_fitted_or_scanned": False,
        "source_predictions_changed": False,
        "empirical_data_used_only_for_comparison": True,
    }

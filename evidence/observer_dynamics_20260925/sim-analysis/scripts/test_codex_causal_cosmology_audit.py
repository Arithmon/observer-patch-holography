"""Independent identities, uncertainty calibration, and receipt audit regressions."""
import importlib.util
import itertools
import math
from fractions import Fraction
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location("causal_audit", Path(__file__).with_name("causal_cosmology_audit.py"))
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_exact_layer_sum_identity():
    for k in range(1, 251):
        direct = Fraction(32 * sum(min(j, k - j) ** 3 for j in range(k + 1)), k ** 4)
        assert AUDIT.parity_factor(k) == direct
    with pytest.raises(ValueError):
        AUDIT.parity_factor(0)


def polynomial_multiply(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def integrate_polynomial(coeffs, lo, hi):
    return sum(c * (hi ** (n + 1) - lo ** (n + 1)) / (n + 1) for n, c in enumerate(coeffs))


def test_radiation_counterexample_by_exact_piecewise_polynomial_integration():
    scale_fourth = [Fraction(math.comb(4, j)) for j in range(5)]
    for t in [Fraction(1, 3), Fraction(1, 2), Fraction(1), Fraction(2)]:
        past = polynomial_multiply(scale_fourth, [Fraction(0)] * 3 + [Fraction(1)])
        future = polynomial_multiply(scale_fourth, [(-1) ** j * math.comb(3, j) * t ** (3 - j) for j in range(4)])
        exact = integrate_polynomial(past, Fraction(0), t / 2) + integrate_polynomial(future, t / 2, t)
        assert exact == AUDIT.radiation_mass_exact(t)
        numeric = AUDIT.continuum_mass(lambda x: 1 + x, float(t))
        assert math.isclose(numeric, float(exact), rel_tol=1e-13)
    ratio4 = AUDIT.radiation_mass_exact(Fraction(1)) / AUDIT.radiation_mass_exact(Fraction(1, 2))
    assert ratio4 == Fraction(1516288, 44451)
    assert ratio4 != Fraction(12, 5) ** 4


def test_stratified_design_variance_against_exhaustive_sampling():
    # Exhaust every without-replacement sample, independently calibrating the
    # estimator's mean and reported variance against its actual design distribution.
    population = [0, 1, 5, 10]
    estimates, variances = [], []
    for chosen in itertools.combinations(population, 2):
        total, var = AUDIT.stratified_moments([{"population": 4, "sample": 2,
            "sum": sum(chosen), "sum_of_squares": sum(x * x for x in chosen)}])
        estimates.append(total)
        variances.append(var)
    truth = sum(population)
    assert sum(estimates) / len(estimates) == truth
    true_design_variance = sum((x - truth) ** 2 for x in estimates) / len(estimates)
    assert sum(variances) / len(variances) == true_design_variance
    with pytest.raises(ValueError):
        AUDIT.stratified_moments([{"population": 4, "sample": 1, "sum": 5, "sum_of_squares": 25}])


def test_pair_normalization_and_uncertainty_mutations():
    import json
    receipt = json.loads((AUDIT.ROOT / AUDIT.INPUTS[1]).read_text())
    row = receipt["levels"][0]["families"][0]["vertical_intervals"][-1]
    result = AUDIT.pair_audit(row)
    assert result["sampling_starts"] == 2004
    assert 0.00023 < result["standard_error"] < 0.00024
    assert result["sampling_only_approximate_95pct_dimension_interval"][0] > 4
    mutated = {**row, "ordering_fraction": str(Fraction(row["ordering_fraction"]) / 2)}
    with pytest.raises(ValueError, match="normalization"):
        AUDIT.pair_audit(mutated)
    mutated = {**row, "ordering_fraction_standard_error": row["ordering_fraction_standard_error"] / 2}
    with pytest.raises(ValueError, match="standard error"):
        AUDIT.pair_audit(mutated)


@pytest.fixture(scope="module")
def result():
    return AUDIT.build()


def test_receipt_clock_baselines_and_correction(result):
    rows = {r["q"]: r for r in result["clock_rows"]}
    assert rows[21]["model_time_ratio"] == 2.5
    assert rows[21]["raw_relative_error"] < -0.19
    assert rows[5]["parity_corrected_clock"] is None  # no fictitious J=1 correction
    assert rows[144]["parity_corrected_relative_error"] == pytest.approx(-0.0010070369, abs=1e-10)
    assert not any(r["flat_count_clock_enclosure_hypotheses_satisfied"] for r in rows.values())
    assert not any(r["outermost_continuum_diamond_inside_cube"] for r in rows.values())
    assert rows[144]["outermost_fraction_of_continuum_volume_clipped"] < 1e-6


def test_continuum_clock_is_not_proper_time_and_methods_agree(result):
    profiles = result["continuum_FLRW_with_supplied_doubling_profiles"]
    assert profiles["constant"]["clock_over_proper_time"] == pytest.approx(1, abs=1e-13)
    assert profiles["de_sitter"]["clock_over_proper_time"] == pytest.approx(0.984046792234, abs=1e-11)
    assert profiles["radiation"]["clock_over_proper_time"] > 1.006
    assert profiles["matter"]["clock_over_proper_time"] > 1.002
    assert all(p["quadrature_max_absolute_disagreement"] < 1e-13 for p in profiles.values())
    assert profiles["de_sitter"]["window_reading_over_midpoint_scale_ratio"] > 1.003

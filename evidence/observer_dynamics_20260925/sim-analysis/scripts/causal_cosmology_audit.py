"""Receipt-only audit of layered clocks, pair sampling, and supplied FLRW profiles.

No large simulation, source-graph reconstruction, or parameter fitting is performed.
The parity correction uses continuum spatial balls at discrete times, not a theorem
that the finite graph's spatial error vanishes at the observed rate.

Run: python3 sim-analysis/scripts/causal_cosmology_audit.py --write
Check: python3 sim-analysis/scripts/causal_cosmology_audit.py --check
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
INPUTS = (
    "data/receipts/source_net_causal_limit_receipt_2026-09-24.json",
    "data/receipts/source_net_q144_sampled_receipt.json",
    "data/receipts/flrw_record_density_readout.json",
)
OUTPUT = ROOT / "data/codex_audit_20260925/causal_cosmology_audit.json"


def parity_factor(k: int) -> Fraction:
    """Exact 32 sum(min(j,k-j)^3)/k^4 for continuum spatial slices."""
    if k < 1:
        raise ValueError("a diamond needs at least one time step")
    return 1 + Fraction(4, k * k) if k % 2 == 0 else (1 - Fraction(1, k * k)) ** 2


def mm_fraction(d: float) -> float:
    return math.exp(math.lgamma(d + 1) + math.lgamma(d / 2)
                    - math.log(2) - math.lgamma(3 * d / 2))


def mm_dimension(fraction: float) -> float:
    lo, hi = 1.0, 20.0
    if not mm_fraction(hi) <= fraction <= mm_fraction(lo):
        raise ValueError("ordering fraction outside inversion range")
    for _ in range(90):
        mid = (lo + hi) / 2
        if mm_fraction(mid) > fraction:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def stratified_moments(rows: list[dict]) -> tuple[Fraction, Fraction]:
    """Rebuild total and design variance with exact rational arithmetic."""
    total = Fraction(0)
    variance = Fraction(0)
    for row in rows:
        n, m, s, ss = (row[k] for k in ("population", "sample", "sum", "sum_of_squares"))
        if not (1 <= m <= n and (m > 1 or m == n)):
            raise ValueError("sample too small for a design variance")
        total += Fraction(n * s, m)
        if m > 1:
            sample_var = Fraction(m * ss - s * s, m * (m - 1))
            if sample_var < 0:
                raise ValueError("inconsistent sum of squares")
            variance += Fraction(n * (n - m), m) * sample_var
    return total, variance


def pair_audit(row: dict) -> dict:
    n = row["inclusive_event_count"]
    estimate = row.get("strict_pair_count_estimate")
    if estimate:
        c, var = stratified_moments(estimate["strata"])
        if c != Fraction(estimate["value"]):
            raise ValueError("stratified total mismatch")
        fraction_se = 2 * math.sqrt(float(var)) / (n * (n - 1))
        if not math.isclose(fraction_se, row["ordering_fraction_standard_error"], rel_tol=1e-10, abs_tol=1e-14):
            raise ValueError("stratified standard error mismatch")
    else:
        c, fraction_se = Fraction(row["strict_pair_count"]), 0.0
    f = 2 * c / (n * (n - 1))
    if f != Fraction(row["ordering_fraction"]):
        raise ValueError("strict-pair normalization mismatch")
    if not math.isclose(mm_dimension(float(f)), row["myrheim_meyer_dimension"], abs_tol=1e-10):
        raise ValueError("Myrheim-Meyer inversion mismatch")
    ci = [max(mm_fraction(20), float(f) - 1.96 * fraction_se), min(1, float(f) + 1.96 * fraction_se)]
    return {"layers": row["layers"], "event_count": n, "ordering_fraction": float(f),
            "standard_error": fraction_se, "sampling_starts": None if not estimate else estimate["sample_size"],
            "sampling_only_approximate_95pct_fraction_interval": ci,
            "mm_dimension": mm_dimension(float(f)),
            "sampling_only_approximate_95pct_dimension_interval": [mm_dimension(ci[1]), mm_dimension(ci[0])],
            "standard_errors_below_0.1": (0.1 - float(f)) / fraction_se if fraction_se else None,
            "continuum_diamond_inside_cube": row.get("continuum_diamond_inside_cube"),
            "interpretation": "normal-approximation sampling interval only; no spatial, layer, boundary or physical uncertainty"}


def profile_functions() -> dict:
    return {"constant": lambda t: 1.0, "de_sitter": lambda t: 1 / (1 - t / 2),
            "radiation": lambda t: 1 + t, "matter": lambda t: (1 + (math.sqrt(2) - 1) * t) ** 2}


def continuum_mass(sigma, duration: float, offset: float = 0, method: str = "quad") -> float:
    """FLRW volume divided by 4pi/3; c=1, comoving vertical diamond."""
    if method == "quad":
        return quad(lambda t: sigma(t + offset) ** 4 * min(t, duration - t) ** 3,
                    0, duration, points=[duration / 2], epsabs=1e-13, epsrel=1e-13)[0]
    if method != "gauss":
        raise ValueError("unknown integration method")
    # Independent fixed-order Gaussian quadrature, split at the cone waist.
    nodes, weights = np.polynomial.legendre.leggauss(48)
    total = 0.0
    for a, b in [(0, duration / 2), (duration / 2, duration)]:
        for node, weight in zip(nodes, weights):
            t = (b + a) / 2 + (b - a) * node / 2
            total += (b - a) * weight / 2 * sigma(t + offset) ** 4 * min(t, duration - t) ** 3
    return float(total)


def continuum_profile(sigma, reference_duration: float = 0.5) -> dict:
    j = reference_duration
    masses = [continuum_mass(sigma, 1), continuum_mass(sigma, j), continuum_mass(sigma, j, 1 - j)]
    independent = [continuum_mass(sigma, 1, method="gauss"), continuum_mass(sigma, j, method="gauss"),
                   continuum_mass(sigma, j, 1 - j, method="gauss")]
    mismatch = max(abs(a - b) for a, b in zip(masses, independent))
    if mismatch > 2e-12:
        raise ValueError("independent quadratures disagree")
    clock = (masses[0] / masses[1]) ** 0.25
    proper_ratio = quad(sigma, 0, 1, epsabs=1e-13)[0] / quad(sigma, 0, j, epsabs=1e-13)[0]
    reading = (masses[2] / masses[1]) ** 0.25
    midpoint = sigma(1 - j / 2) / sigma(j / 2)
    return {"reference_duration": j, "continuum_count_clock": clock, "proper_time_ratio": proper_ratio,
            "clock_over_proper_time": clock / proper_ratio,
            "finite_window_redshift_reading": reading, "midpoint_scale_ratio": midpoint,
            "window_reading_over_midpoint_scale_ratio": reading / midpoint,
            "quadrature_max_absolute_disagreement": mismatch}


def radiation_mass_exact(t: Fraction) -> Fraction:
    """Exact integral for a(t)=1+t, verified independently by polynomial tests."""
    return t ** 4 * (99 * t ** 4 + 672 * t ** 3 + 1792 * t ** 2 + 2240 * t + 1120) / 35840


def build(input_root: Path = ROOT) -> dict:
    inputs = [json.loads((input_root / p).read_text()) for p in INPUTS]
    clocks, pair_rows, even_height = [], [], []
    for receipt in inputs[:2]:
        for level in receipt["levels"]:
            family = next(f for f in level["families"] if f["dimension"] == 3)
            clock = family["count_clock"]
            k, j = clock["interval_layers"], clock["reference_layers"]
            raw = (clock["interval_count"] / clock["reference_count"]) ** 0.25
            b_k, b_j = parity_factor(k), parity_factor(j)
            correction = float(b_k / b_j) ** 0.25 if b_j else None
            top = family["vertical_intervals"][-1]
            clipping = 1 - top["clipped_diamond_volume"] / top["continuum_diamond_volume"]
            clocks.append({"q": family["q"], "K": k, "J": j, "model_time_ratio": k / j,
                           "count_clock": raw, "raw_relative_error": raw / (k / j) - 1,
                           "layer_quadrature_B_K": str(b_k), "layer_quadrature_B_J": str(b_j),
                           "predicted_layer_clock": (k / j) * correction if correction else None,
                           "parity_corrected_clock": raw / correction if correction else None,
                           "parity_corrected_relative_error": raw / (correction * k / j) - 1 if correction else None,
                           "correction_unavailable_reason": None if correction else "J=1 continuum-spatial discrete-time diamond has zero volume",
                           "outermost_continuum_diamond_inside_cube": top["continuum_diamond_inside_cube"],
                           "outermost_fraction_of_continuum_volume_clipped": clipping,
                           "flat_count_clock_enclosure_hypotheses_satisfied": clock["enclosure_hypotheses_satisfied"]})
            for row in family["vertical_intervals"]:
                if row["layers"] >= 2:
                    pair_rows.append({"q": family["q"], **pair_audit(row)})
                if family["q"] == 144 and row["layers"] >= 4 and row["layers"] % 2 == 0:
                    deficit = 0.1 - row["ordering_fraction_float"]
                    even_height.append({"q": 144, "K": row["layers"], "fraction": row["ordering_fraction_float"],
                                        "deficit": deficit, "K_squared_times_deficit": row["layers"] ** 2 * deficit})
    continuum = {name: continuum_profile(sigma) for name, sigma in profile_functions().items()}
    finite_comparison = []
    for level in inputs[2]["levels"]:
        j = level["reference_layers"] / level["interval_layers"]
        for name, sigma in profile_functions().items():
            reference = continuum_profile(sigma, j)
            measured = level["profiles"][name]
            finite_comparison.append({"q": level["q"], "profile": name, "reference_duration": j,
                                     "finite_count_clock": measured["expanding_count_clock"]["physical_clock"],
                                     "continuum_count_clock": reference["continuum_count_clock"],
                                     "finite_over_continuum_clock": measured["expanding_count_clock"]["physical_clock"] / reference["continuum_count_clock"],
                                     "continuum_clock_over_proper_time": reference["clock_over_proper_time"]})
    ratio4 = radiation_mass_exact(Fraction(1)) / radiation_mass_exact(Fraction(1, 2))
    return {"schema": "oph.sim-analysis.causal-cosmology-audit.v1",
            "input_pins": [{"path": p, "sha256": hashlib.sha256((input_root / p).read_bytes()).hexdigest()} for p in INPUTS],
            "method": "receipt arithmetic; continuum-spatial layered-time quadrature; adaptive and Gaussian continuum FLRW integration",
            "clock_rows": clocks, "pair_rows": pair_rows,
            "q144_even_height_diagnostics": even_height,
            "even_height_boundary": "different-height nested diamonds at fixed q144, not independent refinement runs; no covariance available for a joint fitted uncertainty",
            "continuum_FLRW_with_supplied_doubling_profiles": continuum,
            "finite_FLRW_against_matched_continuum": finite_comparison,
            "exact_radiation_counterexample": {"clock_ratio_to_fourth_power": str(ratio4),
                                               "proper_time_ratio_to_fourth_power": str(Fraction(12, 5) ** 4),
                                               "difference": str(ratio4 - Fraction(12, 5) ** 4)},
            "nonclaims": ["no physical clock or scale selected", "no cosmological prediction or fitted CMB observable",
                          "parity correction is for ideal spatial balls, not an exact finite graph law",
                          "a persistent continuum volume-clock bias is not discretization error",
                          "stratified sampling uncertainties exclude systematic geometry and discretization errors",
                          "no Lean files changed; exact parity proof is algebraic, checked by rational tests"]}


def canonical(value: dict) -> str:
    # Stable scientific precision across supported numerical-library patch releases.
    def rounded(x):
        if isinstance(x, float):
            return float(f"{x:.12g}")
        if isinstance(x, list):
            return [rounded(v) for v in x]
        if isinstance(x, dict):
            return {k: rounded(v) for k, v in x.items()}
        return x
    return json.dumps(rounded(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--input-root", type=Path, default=ROOT)
    args = parser.parse_args()
    content = canonical(build(args.input_root))
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(content)
    if args.check and (not OUTPUT.exists() or OUTPUT.read_text() != content):
        raise SystemExit("causal cosmology audit receipt is stale")
    print(f"causal cosmology audit: {len(content)} bytes; sha256={hashlib.sha256(content.encode()).hexdigest()}")


if __name__ == "__main__":
    main()

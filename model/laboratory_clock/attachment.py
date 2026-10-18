"""Exact interval calibration and held-out tests of a declared affine clock.

The candidate clock law is an additional, falsifiable apparatus hypothesis.
Neither timestamp provenance nor physical correspondence is authenticated here.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
import re


@dataclass(frozen=True)
class Endpoint:
    name: str
    count_lo: int
    count_hi: int
    seconds_lo: F
    seconds_hi: F

    def __post_init__(self):
        if (type(self.name) is not str or not self.name
                or type(self.count_lo) is not int or type(self.count_hi) is not int
                or not 0 <= self.count_lo <= self.count_hi
                or type(self.seconds_lo) is not F or type(self.seconds_hi) is not F
                or self.seconds_lo > self.seconds_hi):
            raise ValueError("ordered exact counter and timestamp bounds required")


@dataclass(frozen=True)
class RateFit:
    lower: F
    upper: F | None
    impossible: bool = False

    @property
    def feasible(self):
        return (not self.impossible and (self.upper is None
                or (self.upper > 0 and self.lower <= self.upper)))

    @property
    def identified(self):
        return self.feasible and self.lower > 0 and self.upper is not None

    def report(self):
        return {"feasible": self.feasible, "bounded_away_from_zero_and_infinity": self.identified,
                "seconds_per_attempt_lower": str(self.lower) if self.feasible else None,
                "seconds_per_attempt_upper": (str(self.upper)
                    if self.feasible and self.upper is not None else None),
                "zero_rate_excluded": True}


def fit_rate(points):
    """Eliminate the common offset exactly from rectangular endpoint bounds.

    For rate tau>0, each endpoint permits offsets in
    [t_lo-tau*n_hi, t_hi-tau*n_lo]. Pairwise intersection inequalities
    characterize a common offset. Uncertain integer counts are relaxed to
    their real intervals: infeasibility is still a sound rejection, whereas
    feasibility need not give a compatible choice of integer counts.
    """
    points = tuple(points)
    if not points or any(type(p) is not Endpoint for p in points):
        raise ValueError("at least one validated endpoint required")
    lower, upper, impossible = F(0), None, False
    for a in points:
        for b in points:
            delta_count = a.count_hi-b.count_lo
            delta_seconds = a.seconds_lo-b.seconds_hi
            if delta_count > 0:
                lower = max(lower, delta_seconds/delta_count)
            elif delta_count < 0:
                bound = delta_seconds/delta_count
                upper = bound if upper is None else min(upper, bound)
            elif delta_seconds > 0:
                impossible = True
    return RateFit(lower, upper, impossible)


def offset_bounds(points, rate):
    """Allowed common origins at a specified positive seconds-per-attempt."""
    if type(rate) is not F or rate <= 0:
        raise ValueError("a positive exact rate required")
    return (max(p.seconds_lo-rate*p.count_hi for p in points),
            min(p.seconds_hi-rate*p.count_lo for p in points))


def predict_seconds(calibration, count_lo, count_hi):
    """Sharp outer timestamp interval, retaining origin/rate correlation."""
    if (type(count_lo) is not int or type(count_hi) is not int
            or not 0 <= count_lo <= count_hi):
        raise ValueError("ordered integer counter bounds required")
    calibration = tuple(calibration)
    rate = fit_rate(calibration)
    if not rate.identified:
        raise ValueError("calibration must bound its positive rate")
    # Piecewise affine extrema occur at interval ends or line intersections.
    lower_lines = [(p.seconds_lo, count_lo-p.count_hi) for p in calibration]
    upper_lines = [(p.seconds_hi, count_hi-p.count_lo) for p in calibration]

    def candidates(lines):
        rates = {rate.lower, rate.upper}
        for (a, b), (c, d) in combinations(lines, 2):
            if b != d:
                value = (c-a)/(b-d)
                if rate.lower <= value <= rate.upper:
                    rates.add(value)
        return rates

    lower = min(max(a+b*r for a, b in lower_lines) for r in candidates(lower_lines))
    upper = max(min(a+b*r for a, b in upper_lines) for r in candidates(upper_lines))
    return lower, upper


def validate(calibration, held_out):
    """Test untouched endpoints against the original calibration set."""
    calibration, held_out = tuple(calibration), tuple(held_out)
    if len(calibration) < 2 or not held_out:
        raise ValueError("at least two calibration and one held-out endpoint required")
    if len({p.name for p in calibration+held_out}) != len(calibration+held_out):
        raise ValueError("endpoint identifiers must be unique and splits disjoint")
    relaxed = any(p.count_lo != p.count_hi for p in calibration+held_out)
    fitted = fit_rate(calibration)
    result = {"calibration": fitted.report(), "calibration_ids": [p.name for p in calibration],
              "held_out_ids": [p.name for p in held_out], "held_out_used_in_calibration": False,
              "count_bounds_relaxed_to_real_intervals": relaxed,
              "integer_count_realization_certified": False,
              "true_count_order_enforced_with_uncertain_counts": False,
              "compatibility_scope": ("outer feasibility only; integer counts and their joint order are not certified"
                  if relaxed else "one common affine clock at the exact counter endpoints"),
              "per_endpoint": [], "joint_validation": "unidentified_calibration"}
    if not fitted.feasible:
        result["joint_validation"] = "incompatible_calibration"
        return result
    if not fitted.identified:
        return result
    for point in held_out:
        lo, hi = predict_seconds(calibration, point.count_lo, point.count_hi)
        result["per_endpoint"].append({"id": point.name,
            "prediction_seconds": [str(lo), str(hi)],
            "observed_seconds": [str(point.seconds_lo), str(point.seconds_hi)],
            "compatible_with_calibration": lo <= point.seconds_hi and point.seconds_lo <= hi})
    joint = fit_rate(calibration+held_out)
    result["joint_validation"] = ("compatible_with_some_calibration" if joint.feasible
                                   else "incompatible_with_every_calibration")
    result["integer_count_realization_certified"] = not relaxed and joint.feasible
    # Deliberately do not export a new rate estimated using validation data.
    return result


def _keys(value, expected, label):
    if type(value) is not dict or set(value) != set(expected):
        raise ValueError(f"{label} has missing or unrecognized fields")


def _rational(value, label):
    if type(value) is not str or re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value) is None:
        raise ValueError(f"{label} must be an exact rational string")
    return F(value)


def analyse(data):
    """Read a bounded prospective endpoint contract; unresolved bounds stay null."""
    _keys(data, ("schema", "data_kind", "counter", "provenance", "endpoints"), "clock input")
    if data["schema"] != "oph.affine_clock_input.v1":
        raise ValueError("unrecognized clock input schema")
    if data["data_kind"] not in ("synthetic_fixture", "host_measurement", "laboratory_measurement"):
        raise ValueError("declare synthetic, host or laboratory input")
    if data["counter"] != "native_attempts_including_zero_currents":
        raise ValueError("this contract requires the native attempt counter, including zero currents")
    provenance_names = ("event_mapping", "preparation_and_scheduler", "detection_bounds",
                        "timestamp_bounds", "traceability", "split_precommitment")
    _keys(data["provenance"], provenance_names, "provenance")
    for value in data["provenance"].values():
        if value is not None and (type(value) is not str or not value.strip()):
            raise ValueError("provenance must be a nonempty reference or null")
    rows = data["endpoints"]
    if type(rows) is not list or not 3 <= len(rows) <= 64:
        raise ValueError("supply between three and 64 endpoints")
    calibration, held_out, missing, ids = [], [], [], set()
    previous_count = -1
    for row in rows:
        _keys(row, ("id", "role", "attempt_count", "count_error", "seconds", "seconds_error"), "endpoint")
        name, role, n = row["id"], row["role"], row["attempt_count"]
        if type(name) is not str or not name or name in ids:
            raise ValueError("unique nonempty endpoint ids required")
        ids.add(name)
        if role not in ("calibration", "validation"):
            raise ValueError("each endpoint needs a fixed calibration or validation role")
        if type(n) is not int or n <= previous_count:
            raise ValueError("observed native attempt counts must increase strictly from zero or above")
        previous_count = n
        seconds = _rational(row["seconds"], "seconds")
        count_error, seconds_error = row["count_error"], row["seconds_error"]
        if count_error is not None and (type(count_error) is not int or count_error < 0):
            raise ValueError("count errors must be nonnegative integers or null")
        if seconds_error is not None:
            seconds_error = _rational(seconds_error, "seconds error")
            if seconds_error < 0:
                raise ValueError("timestamp errors must be nonnegative or null")
        if count_error is None or seconds_error is None:
            missing.append(name)
            continue
        endpoint = Endpoint(name, max(0, n-count_error), n+count_error,
                            seconds-seconds_error, seconds+seconds_error)
        (calibration if role == "calibration" else held_out).append(endpoint)
    if sum(row["role"] == "calibration" for row in rows) < 2 or not any(
            row["role"] == "validation" for row in rows):
        raise ValueError("declare at least two calibration and one validation endpoints")
    return {"data_kind": data["data_kind"], "counter": data["counter"],
            "provenance": data["provenance"], "provenance_authenticated": False,
            "unresolved_provenance": [key for key, value in data["provenance"].items() if value is None],
            "unresolved_error_endpoints": missing,
            "numerical_test": (None if missing else validate(calibration, held_out)),
            "laboratory_correspondence_established": False,
            "interpretation": "conditional compatibility test, not an identification of physical time"}

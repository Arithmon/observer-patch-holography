#!/usr/bin/env python3
"""Validate finite Z2 receipt bindings and portable numerical replay.

This verifier imports no numerical producer. It checks the complete declared
grid, source and payload digests, each reported component, and the analytic free
and Kogut--Susskind identities. Comparing two interacting results is portable
replay, not an independent interacting eigensolve or a certified interval bound.
Independent original-input numerical controls belong in the scientific tests.
Self-digests bind local bytes; they do not authenticate their origin.

Substantive numbers have relative-only replay tolerance. Absolute error is
allowed only for diagnostics whose exact value is zero: the free residual,
influence and mass outside the single-flip support. Their budgets
scale with binary64 epsilon, matrix dimension and fresh rates, including the
number of entries in an accumulated off-diagonal mass. These are conservative
portable roundoff allowances for the resolved default grid, not error proofs.
The KS construction preserves its structural zero entries exactly, so its mass
outside the single-flip support must be exactly zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

REPLAY_RTOL = 1e-7
DEFAULT_GRID = {
    "L": [2, 3],
    "wilson_diagonal_beta_s_eq_beta_t": [0.1, 0.3, 0.5, 0.7, 1.0],
    "kogut_susskind_lambda": [0.5, 1.0, 2.0],
    "universal_no_go": False,
}
METADATA = {
    "schema": "oph.yang_mills.z2_finite_transfer_receipt.v2",
    "scope": "finite_gauge_diagnostic",
    "physical_clay_receipt": False,
    "system": "Z2 lattice gauge theory, L x L periodic spatial torus, gauge-invariant sector",
    "receipt_under_test": (
        "finite ground-state-transform and cross-fiber receipt: "
        "U_r H_r U_r^{-1} = sum_C c_C (I - E_C) with c_C independent of the repaired value"
    ),
    "ground_state_transform": "Doob transform by the Perron vector, pi = Omega^2",
    "collars": "one collar per spatial link, fiber {o, X_l o}, pi-preserving heat bath",
}
FLOOR_TEXT = {
    "identity": "c_l(o) = lambda * (r_l(o) + 1/r_l(o))",
    "analytic_lower_bound": "c_l(o) >= 2 * lambda by AM-GM",
    "scope": (
        "finite Kogut-Susskind Doob transform; quotient-space "
        "approximate tensorization and continuum transfer not proved"
    ),
}


class ReceiptValidationError(ValueError):
    """Malformed, stale, inconsistent or numerically unreproduced receipt."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ReceiptValidationError(message)


def _object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result = {}
    for key, value in pairs:
        _require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _nonfinite_constant(value: str) -> None:
    raise ReceiptValidationError(f"nonfinite JSON constant: {value}")


def _json_float(value: str) -> float:
    number = float(value)
    _require(math.isfinite(number), f"nonfinite JSON number: {value}")
    return number


def parse_receipt(serialized: str | bytes) -> dict[str, Any]:
    """Parse strict JSON, including when the caller binds the same raw bytes."""
    value = json.loads(
        serialized,
        object_pairs_hook=_object_pairs,
        parse_constant=_nonfinite_constant,
        parse_float=_json_float,
    )
    _require(type(value) is dict, "receipt must be an object")
    return value


def load_receipt(path: Path) -> dict[str, Any]:
    """Read strict JSON, refusing duplicate keys and nonfinite numbers."""
    return parse_receipt(Path(path).read_bytes())


def _keys(value: Any, expected: set[str], path: str) -> None:
    _require(type(value) is dict, f"{path}: expected object")
    _require(set(value) == expected, f"{path}: missing or unexpected fields")


def _number(value: Any, path: str, *, positive: bool = False,
            nonnegative: bool = False) -> float:
    _require(type(value) in (int, float), f"{path}: expected real number, not bool")
    try:
        result = float(value)
    except OverflowError as error:
        raise ReceiptValidationError(f"{path}: number is outside finite range") from error
    _require(math.isfinite(result), f"{path}: expected finite number")
    if positive:
        _require(result > 0, f"{path}: expected positive number")
    if nonnegative:
        _require(result >= 0, f"{path}: expected nonnegative number")
    return result


def _bool(value: Any, path: str) -> None:
    _require(type(value) is bool, f"{path}: expected bool")


def _same(value: Any, expected: Any, path: str) -> None:
    _require(type(value) is type(expected) and value == expected,
             f"{path}: unexpected value")


def _relative(actual: float, expected: float, path: str,
              tolerance: float = REPLAY_RTOL) -> None:
    # Dividing avoids an underflowing absolute tolerance for subnormal values.
    actual, expected = float(actual), float(expected)
    if expected == 0:
        _require(actual == 0, f"{path}: expected exact zero")
    else:
        _require(abs(actual / expected - 1.0) <= tolerance,
                 f"{path}: relative numerical mismatch")


def _grid(value: Any, expected: dict[str, Any]) -> list[tuple[int, str, dict, str | None]]:
    _keys(value, set(DEFAULT_GRID), "grid_scope")
    _same(value["universal_no_go"], False, "grid_scope.universal_no_go")
    for key in ("L", "wilson_diagonal_beta_s_eq_beta_t", "kogut_susskind_lambda"):
        values = value[key]
        _require(type(values) is list and bool(values), f"grid_scope.{key}: empty/non-list grid")
        if key == "L":
            _require(all(type(x) is int and x in (2, 3) for x in values),
                     "grid_scope.L: verifier supports the resolved L=2,3 tori")
        else:
            for item in values:
                _number(item, f"grid_scope.{key}", positive=True)
        _require(len(set(values)) == len(values), f"grid_scope.{key}: duplicate entry")
    _require(value == expected, "grid_scope: differs from requested replay grid")
    rows = []
    for size in value["L"]:
        rows.append((size, "wilson", {"beta_s": 0.0, "beta_t": 0.5}, "beta_s_zero"))
        rows.extend((size, "wilson", {"beta_s": b, "beta_t": b}, None)
                    for b in value["wilson_diagonal_beta_s_eq_beta_t"])
        rows.extend((size, "kogut_susskind", {"lam": lam}, None)
                    for lam in value["kogut_susskind_lambda"])
    return rows


def _zero_budget(run: dict[str, Any], kind: str, *, rate_scale: float | None = None) -> float:
    dimension = run["n_orbits"]
    if kind == "outside":
        if rate_scale is None:
            rate_scale = max(abs(x) for x in run["constant_rate_fit"]["rates"])
        return 64 * sys.float_info.epsilon * dimension**2 * run["n_links"] * rate_scale
    return 128 * sys.float_info.epsilon * dimension


def _validate_run(run: Any, specification: tuple, index: int) -> None:
    size, transfer, parameters, control = specification
    path = f"runs[{index}]"
    fields = {
        "L", "transfer", "parameters", "n_orbits", "n_links",
        "doob_generator_rows_sum_zero", "doob_generator_offdiagonal_nonpositive",
        "constant_rate_fit", "fiber_dependent_rates", "dobrushin", "spectral",
    }
    fields.add("lambda_max" if transfer == "wilson" else "ground_energy")
    if transfer == "kogut_susskind":
        fields.add("variable_rate_floor")
    if control is not None:
        fields.add("control")
    _keys(run, fields, path)
    for key, expected in (("L", size), ("transfer", transfer),
                          ("n_links", 2 * size**2), ("n_orbits", 2 ** (size**2 + 1))):
        _same(run[key], expected, f"{path}.{key}")
    if control is not None:
        _same(run["control"], control, f"{path}.control")
    _keys(run["parameters"], set(parameters), f"{path}.parameters")
    for key, expected in parameters.items():
        value = _number(run["parameters"][key], f"{path}.parameters.{key}")
        _require(value == expected, f"{path}.parameters.{key}: wrong grid entry")
    for key in ("doob_generator_rows_sum_zero", "doob_generator_offdiagonal_nonpositive"):
        _bool(run[key], f"{path}.{key}")
    _same(run["doob_generator_rows_sum_zero"], True, f"{path}.doob_generator_rows_sum_zero")
    if control is not None or transfer == "kogut_susskind":
        _same(run["doob_generator_offdiagonal_nonpositive"], True,
              f"{path}.doob_generator_offdiagonal_nonpositive")

    fit = run["constant_rate_fit"]
    _keys(fit, {"rates", "rate_min", "rate_max", "relative_frobenius_residual"}, path + ".constant_rate_fit")
    _require(type(fit["rates"]) is list and len(fit["rates"]) == run["n_links"],
             path + ".constant_rate_fit.rates: wrong length")
    for position, value in enumerate(fit["rates"]):
        _number(value, f"{path}.constant_rate_fit.rates[{position}]", positive=True)
    for key in ("rate_min", "rate_max"):
        _number(fit[key], f"{path}.constant_rate_fit.{key}", positive=True)
    _require(fit["rate_min"] == min(fit["rates"]) and fit["rate_max"] == max(fit["rates"]),
             path + ".constant_rate_fit: inconsistent rate extrema")
    _number(fit["relative_frobenius_residual"], path + ".constant_rate_fit.relative_frobenius_residual", nonnegative=True)

    fiber = run["fiber_dependent_rates"]
    _keys(fiber, {"rate_min", "rate_max", "spread_max_over_min", "all_rates_positive",
                  "offdiagonal_mass_outside_single_flip", "offdiagonal_mass_single_flip"},
          path + ".fiber_dependent_rates")
    for key in ("rate_min", "rate_max", "spread_max_over_min", "offdiagonal_mass_single_flip"):
        _number(fiber[key], f"{path}.fiber_dependent_rates.{key}", positive=True)
    _number(fiber["offdiagonal_mass_outside_single_flip"], path + ".fiber_dependent_rates.offdiagonal_mass_outside_single_flip", nonnegative=True)
    _same(fiber["all_rates_positive"], True, path + ".fiber_dependent_rates.all_rates_positive")
    _require(fiber["rate_min"] <= fiber["rate_max"], path + ".fiber_dependent_rates: inverted extrema")
    _relative(fiber["spread_max_over_min"], fiber["rate_max"] / fiber["rate_min"],
              path + ".fiber_dependent_rates.spread_max_over_min", 8 * sys.float_info.epsilon)

    dob = run["dobrushin"]
    _keys(dob, {"eta_star", "dobrushin_condition_holds", "unit_rate_floor_c_star_times_1_minus_eta"}, path + ".dobrushin")
    eta = _number(dob["eta_star"], path + ".dobrushin.eta_star", nonnegative=True)
    _bool(dob["dobrushin_condition_holds"], path + ".dobrushin.dobrushin_condition_holds")
    _require(dob["dobrushin_condition_holds"] == (eta < 1), path + ".dobrushin: inconsistent condition")
    floor = _number(dob["unit_rate_floor_c_star_times_1_minus_eta"], path + ".dobrushin.unit_rate_floor_c_star_times_1_minus_eta", nonnegative=True)
    _require(floor == max(0.0, 1 - eta), path + ".dobrushin: inconsistent unit-rate floor")

    spectral = run["spectral"]
    _keys(spectral, {"gap_H", "gap_unit_rate_heat_bath", "pi_min", "pi_max"}, path + ".spectral")
    for key, value in spectral.items():
        _number(value, f"{path}.spectral.{key}", positive=True)
    _require(spectral["pi_min"] <= 1 / run["n_orbits"] <= spectral["pi_max"] <= 1,
             path + ".spectral: impossible probability extrema")
    if transfer == "wilson":
        _number(run["lambda_max"], path + ".lambda_max", positive=True)
    else:
        _number(run["ground_energy"], path + ".ground_energy")
        # Off-support matrix entries are zero before diagonal Doob scaling;
        # this quantity is not a cancellation-prone eigensolver residual.
        _require(fiber["offdiagonal_mass_outside_single_flip"] == 0,
                 path + ".fiber_dependent_rates: KS outside-support mass must be exactly zero")
        block = run["variable_rate_floor"]
        _keys(block, set(FLOOR_TEXT) | {"lower_bound_value", "numerical_min_respects_bound"}, path + ".variable_rate_floor")
        for key, expected in FLOOR_TEXT.items():
            _same(block[key], expected, f"{path}.variable_rate_floor.{key}")
        bound = _number(block["lower_bound_value"], path + ".variable_rate_floor.lower_bound_value", positive=True)
        _require(bound == 2 * parameters["lam"], path + ".variable_rate_floor: wrong analytic bound")
        _same(block["numerical_min_respects_bound"], True, path + ".variable_rate_floor.numerical_min_respects_bound")
        _require(fiber["rate_min"] / bound >= 1 - REPLAY_RTOL,
                 path + ".variable_rate_floor: minimum contradicts analytic bound")

    if control is not None:
        # The shortest nontrivial electric cycle has length L on these two tori.
        rate = math.log1p(2.0 / math.expm1(2.0 * parameters["beta_t"]))
        dimensionless_budget = _zero_budget(run, "dimensionless")
        for key, value in (
            ("constant_rate_fit.relative_frobenius_residual", fit["relative_frobenius_residual"]),
            ("dobrushin.eta_star", eta),
        ):
            _require(value <= dimensionless_budget,
                     f"{path}.{key}: free zero diagnostic exceeds roundoff budget")
        # Bind the allowance to the analytic input rate, not a supplied fit.
        _require(fiber["offdiagonal_mass_outside_single_flip"] <=
                 _zero_budget(run, "outside", rate_scale=rate),
                 path + ".fiber_dependent_rates.offdiagonal_mass_outside_single_flip: free zero diagnostic exceeds roundoff budget")
        for position, value in enumerate(fit["rates"]):
            _relative(value, rate, f"{path}.constant_rate_fit.rates[{position}]: free identity")
        for key in ("rate_min", "rate_max"):
            _relative(fiber[key], rate, f"{path}.fiber_dependent_rates.{key}: free identity")
        for key, expected in (("gap_H", size * rate), ("gap_unit_rate_heat_bath", size),
                              ("pi_min", 1 / run["n_orbits"]), ("pi_max", 1 / run["n_orbits"])):
            _relative(spectral[key], expected, f"{path}.spectral.{key}: free identity")
        _relative(run["lambda_max"], (2 * math.cosh(parameters["beta_t"])) ** run["n_links"],
                  path + ".lambda_max: free identity")


def verify_bindings(receipt: dict[str, Any], producer_path: Path, *,
                    expected_grid: dict[str, Any] | None = None) -> None:
    """Check complete structure, local byte bindings and scalar consistency."""
    _keys(receipt, set(METADATA) | {"grid_scope", "producer_sha256", "sha256_of_runs", "runs"}, "receipt")
    for key, expected in METADATA.items():
        _same(receipt[key], expected, key)
    expected_digest = hashlib.sha256(Path(producer_path).read_bytes()).hexdigest()
    _same(receipt["producer_sha256"], expected_digest, "producer_sha256")
    specifications = _grid(receipt["grid_scope"], DEFAULT_GRID if expected_grid is None else expected_grid)
    _require(type(receipt["runs"]) is list and len(receipt["runs"]) == len(specifications),
             "runs: incomplete or extra grid rows")
    for index, (run, specification) in enumerate(zip(receipt["runs"], specifications)):
        _validate_run(run, specification, index)
    payload = json.dumps(receipt["runs"], sort_keys=True, allow_nan=False).encode("utf-8")
    _same(receipt["sha256_of_runs"], hashlib.sha256(payload).hexdigest(), "sha256_of_runs")


def _compare_run(retained: dict, fresh: dict, index: int) -> None:
    free = fresh.get("control") == "beta_s_zero"
    absolute_paths = {}
    if free:
        absolute_paths[("constant_rate_fit", "relative_frobenius_residual")] = _zero_budget(fresh, "dimensionless")
        absolute_paths[("dobrushin", "eta_star")] = _zero_budget(fresh, "dimensionless")
        absolute_paths[("fiber_dependent_rates", "offdiagonal_mass_outside_single_flip")] = _zero_budget(fresh, "outside")

    def compare(actual: Any, expected: Any, path: tuple) -> None:
        label = f"runs[{index}]." + ".".join(str(x) for x in path)
        if type(expected) is dict:
            for key in expected:
                compare(actual[key], expected[key], path + (key,))
        elif type(expected) is list:
            for position, value in enumerate(expected):
                compare(actual[position], value, path + (position,))
        elif type(expected) in (int, float) and path not in (("L",), ("n_links",), ("n_orbits",)) and path[:1] != ("parameters",):
            if path in absolute_paths:
                # Bound each operand, rather than allowing a shared false zero claim.
                budget = absolute_paths[path]
                _require(actual <= budget and expected <= budget, label + ": zero diagnostic exceeds roundoff budget")
            else:
                _relative(actual, expected, label)
        else:
            _require(actual == expected, label + ": metadata or classification mismatch")

    compare(retained, fresh, ())


def verify_receipt(receipt: dict[str, Any], replay: dict[str, Any], producer_path: Path, *,
                   expected_grid: dict[str, Any] | None = None) -> None:
    """Validate bindings, then compare every component with a supplied fresh run."""
    verify_bindings(receipt, producer_path, expected_grid=expected_grid)
    verify_bindings(replay, producer_path, expected_grid=expected_grid)
    for index, (retained, fresh) in enumerate(zip(receipt["runs"], replay["runs"])):
        _compare_run(retained, fresh, index)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("receipt", type=Path)
    parser.add_argument("replay", type=Path, help="fresh output from the canonical producer")
    parser.add_argument("--producer", type=Path, default=Path(__file__).with_name("z2_finite_transfer_receipt.py"))
    args = parser.parse_args()
    verify_receipt(load_receipt(args.receipt), load_receipt(args.replay), args.producer)
    print("finite Z2 receipt bindings, analytic identities and portable replay verified")


if __name__ == "__main__":
    main()

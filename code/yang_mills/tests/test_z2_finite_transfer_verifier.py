"""Adversarial tests for componentwise finite-transfer receipt replay."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import verify_z2_finite_transfer_receipt as verifier  # noqa: E402


def _rehash(receipt: dict) -> None:
    receipt["sha256_of_runs"] = hashlib.sha256(
        json.dumps(receipt["runs"], sort_keys=True).encode("utf-8")
    ).hexdigest()


@pytest.fixture
def specimen(tmp_path):
    # No eigensolve is needed to challenge the serialized verifier. These
    # retained L=2 values are rebound to explicit fixture bytes, not presented
    # as a new production receipt or used as an independent numerical oracle.
    receipt = verifier.load_receipt(HERE / "receipts/z2_finite_transfer_receipt.json")
    receipt["grid_scope"]["L"] = [2]
    receipt["runs"] = [run for run in receipt["runs"] if run["L"] == 2]
    producer = tmp_path / "fixture_producer.py"
    producer.write_bytes(b"bounded serialized-verifier fixture\n")
    receipt["producer_sha256"] = hashlib.sha256(producer.read_bytes()).hexdigest()
    _rehash(receipt)
    return receipt, producer, copy.deepcopy(receipt["grid_scope"])


def _verify(specimen, receipt=None, replay=None):
    original, producer, grid = specimen
    verifier.verify_receipt(
        original if receipt is None else receipt,
        original if replay is None else replay,
        producer,
        expected_grid=grid,
    )


def test_valid_receipt_and_portable_roundoff(specimen):
    original, _, _ = specimen
    changed = copy.deepcopy(original)
    run = changed["runs"][1]
    for key in run["spectral"]:
        run["spectral"][key] *= 1 + 1e-10
    run["lambda_max"] *= 1 + 1e-10
    _rehash(changed)
    assert changed != original
    _verify(specimen, changed)


@pytest.mark.parametrize("key", ["gap_H", "gap_unit_rate_heat_bath", "pi_min", "pi_max"])
def test_each_spectral_component_is_checked(specimen, key):
    changed = copy.deepcopy(specimen[0])
    changed["runs"][1]["spectral"][key] *= 1.01
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match=key):
        _verify(specimen, changed)


@pytest.mark.parametrize("scale", [1e-200, 1e-310])
@pytest.mark.parametrize("value", [0.0, 1.01])
@pytest.mark.parametrize("field", ["gap_H", "pi_min"])
def test_tiny_nonzero_components_have_no_absolute_floor(specimen, scale, value, field):
    reference = copy.deepcopy(specimen[0])
    reference["runs"][1]["spectral"][field] = scale
    _rehash(reference)
    changed = copy.deepcopy(reference)
    changed["runs"][1]["spectral"][field] = scale * value
    _rehash(changed)
    _verify(specimen, reference, reference)
    with pytest.raises(verifier.ReceiptValidationError, match=field):
        _verify(specimen, changed, reference)


def test_correlated_component_errors_cannot_hide_in_rate_spread(specimen):
    changed = copy.deepcopy(specimen[0])
    fiber = changed["runs"][1]["fiber_dependent_rates"]
    fiber["rate_min"] *= 1.01
    fiber["rate_max"] *= 1.01
    fiber["spread_max_over_min"] = fiber["rate_max"] / fiber["rate_min"]
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="rate_min|rate_max"):
        _verify(specimen, changed)


def test_every_substantive_reported_component_is_checked(specimen):
    # One interacting Wilson row and one KS row exercise both output variants.
    # Each mutation changes one reported scalar while retaining the payload's
    # self-hash; internal consistency checks may detect it before replay does.
    original = specimen[0]
    positions = [1, next(i for i, row in enumerate(original["runs"])
                         if row["transfer"] == "kogut_susskind")]
    for position in positions:
        row = original["runs"][position]
        components = [
            ("constant_rate_fit", "rate_min"),
            ("constant_rate_fit", "rate_max"),
            ("constant_rate_fit", "relative_frobenius_residual"),
            ("fiber_dependent_rates", "rate_min"),
            ("fiber_dependent_rates", "rate_max"),
            ("fiber_dependent_rates", "spread_max_over_min"),
            ("fiber_dependent_rates", "offdiagonal_mass_single_flip"),
            ("dobrushin", "eta_star"),
            ("dobrushin", "unit_rate_floor_c_star_times_1_minus_eta"),
        ]
        if row["transfer"] == "wilson":
            components.append(("fiber_dependent_rates", "offdiagonal_mass_outside_single_flip"))
        else:
            components.append(("variable_rate_floor", "lower_bound_value"))
        for block, key in components:
            changed = copy.deepcopy(original)
            before = changed["runs"][position][block][key]
            changed["runs"][position][block][key] = before * 1.01 if before else 0.01
            _rehash(changed)
            with pytest.raises(verifier.ReceiptValidationError):
                _verify(specimen, changed)
        for index in range(row["n_links"]):
            changed = copy.deepcopy(original)
            changed["runs"][position]["constant_rate_fit"]["rates"][index] *= 1.01
            _rehash(changed)
            with pytest.raises(verifier.ReceiptValidationError):
                _verify(specimen, changed)
        key = "lambda_max" if row["transfer"] == "wilson" else "ground_energy"
        changed = copy.deepcopy(original)
        changed["runs"][position][key] *= 1.01
        _rehash(changed)
        with pytest.raises(verifier.ReceiptValidationError, match=key):
            _verify(specimen, changed)


def test_free_identity_is_checked_even_for_matching_forged_replays(specimen):
    changed = copy.deepcopy(specimen[0])
    changed["runs"][0]["spectral"]["gap_H"] *= 1.01
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="free identity"):
        _verify(specimen, changed, changed)


@pytest.mark.parametrize("field", ["relative_frobenius_residual", "eta_star", "outside"])
def test_roundoff_budget_is_limited_to_analytic_zero_diagnostics(specimen, field):
    reference = copy.deepcopy(specimen[0])
    changed = copy.deepcopy(reference)
    run = changed["runs"][0]
    if field == "relative_frobenius_residual":
        block, key = run["constant_rate_fit"], field
    elif field == "eta_star":
        block, key = run["dobrushin"], field
    else:
        block, key = run["fiber_dependent_rates"], "offdiagonal_mass_outside_single_flip"
    block[key] = 1e-14
    if field == "eta_star":
        block["unit_rate_floor_c_star_times_1_minus_eta"] = 1 - block[key]
    _rehash(changed)
    _verify(specimen, changed, reference)
    block[key] = 1e-3
    if field == "eta_star":
        block["unit_rate_floor_c_star_times_1_minus_eta"] = 1 - block[key]
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="roundoff budget"):
        _verify(specimen, changed, changed)


def test_interacting_influence_and_outside_mass_use_relative_comparison(specimen):
    for block, key in [
        ("dobrushin", "eta_star"),
        ("fiber_dependent_rates", "offdiagonal_mass_outside_single_flip"),
        ("constant_rate_fit", "relative_frobenius_residual"),
    ]:
        reference = copy.deepcopy(specimen[0])
        reference["runs"][1][block][key] = 1e-200
        if block == "dobrushin":
            reference["runs"][1][block]["dobrushin_condition_holds"] = True
            reference["runs"][1][block]["unit_rate_floor_c_star_times_1_minus_eta"] = 1.0
        _rehash(reference)
        changed = copy.deepcopy(reference)
        changed["runs"][1][block][key] = 0.0
        _rehash(changed)
        with pytest.raises(verifier.ReceiptValidationError, match=key):
            _verify(specimen, changed, reference)


def test_ks_floor_checks_minimum_not_only_boolean(specimen):
    changed = copy.deepcopy(specimen[0])
    run = next(row for row in changed["runs"] if row["transfer"] == "kogut_susskind")
    fiber = run["fiber_dependent_rates"]
    fiber["rate_min"] = 0.5 * run["variable_rate_floor"]["lower_bound_value"]
    fiber["spread_max_over_min"] = fiber["rate_max"] / fiber["rate_min"]
    assert run["variable_rate_floor"]["numerical_min_respects_bound"] is True
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="contradicts analytic bound"):
        _verify(specimen, changed, changed)


@pytest.mark.parametrize("quantity", ["residual", "influence", "outside"])
def test_bindings_enforce_known_free_zeros_without_replay(specimen, quantity):
    changed = copy.deepcopy(specimen[0])
    run = changed["runs"][0]
    if quantity == "residual":
        run["constant_rate_fit"]["relative_frobenius_residual"] = 0.5
    elif quantity == "influence":
        run["dobrushin"]["eta_star"] = 0.5
        run["dobrushin"]["unit_rate_floor_c_star_times_1_minus_eta"] = 0.5
    else:
        run["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = 100.0
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="free zero diagnostic"):
        verifier.verify_bindings(changed, specimen[1], expected_grid=specimen[2])


@pytest.mark.parametrize("nonzero", [5e-324, 1e-20, 1.0])
def test_bindings_require_structural_ks_outside_zero(specimen, nonzero):
    changed = copy.deepcopy(specimen[0])
    local = next(run for run in changed["runs"] if run["transfer"] == "kogut_susskind")
    local["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = nonzero
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="must be exactly zero"):
        verifier.verify_bindings(changed, specimen[1], expected_grid=specimen[2])


@pytest.mark.parametrize("mutation", [
    "empty_runs", "missing_run", "duplicate_run", "extra_run", "reordered_runs",
    "missing_field", "extra_field", "wrong_grid", "empty_grid", "duplicate_grid",
    "physical_claim", "universal_claim", "wrong_flag", "wrong_parameter",
    "bool_number", "number_bool", "wrong_dimension", "wrong_rate_length",
    "stale_source", "stale_runs",
])
def test_incomplete_and_malformed_receipts_fail_closed(specimen, mutation):
    changed = copy.deepcopy(specimen[0])
    if mutation == "empty_runs":
        changed["runs"] = []
    elif mutation == "missing_run":
        changed["runs"].pop()
    elif mutation == "duplicate_run":
        changed["runs"][2] = copy.deepcopy(changed["runs"][1])
    elif mutation == "extra_run":
        changed["runs"].append(copy.deepcopy(changed["runs"][-1]))
    elif mutation == "reordered_runs":
        changed["runs"][1:3] = changed["runs"][1:3][::-1]
    elif mutation == "missing_field":
        del changed["runs"][1]["spectral"]["gap_H"]
    elif mutation == "extra_field":
        changed["runs"][1]["spectral"]["certified"] = True
    elif mutation == "wrong_grid":
        changed["grid_scope"]["L"] = [3]
    elif mutation == "empty_grid":
        changed["grid_scope"]["L"] = []
    elif mutation == "duplicate_grid":
        changed["grid_scope"]["L"] = [2, 2]
    elif mutation == "physical_claim":
        changed["physical_clay_receipt"] = True
    elif mutation == "universal_claim":
        changed["grid_scope"]["universal_no_go"] = True
    elif mutation == "wrong_flag":
        changed["runs"][1]["doob_generator_offdiagonal_nonpositive"] = True
    elif mutation == "wrong_parameter":
        changed["runs"][1]["parameters"]["beta_t"] = 0.7
    elif mutation == "bool_number":
        changed["runs"][1]["spectral"]["gap_H"] = True
    elif mutation == "number_bool":
        changed["runs"][1]["doob_generator_rows_sum_zero"] = 1
    elif mutation == "wrong_dimension":
        changed["runs"][1]["n_links"] = 9
    elif mutation == "wrong_rate_length":
        changed["runs"][1]["constant_rate_fit"]["rates"].pop()
    elif mutation == "stale_source":
        changed["producer_sha256"] = "0" * 64
    elif mutation == "stale_runs":
        changed["runs"][1]["spectral"]["gap_H"] *= 1 + 1e-10
    if mutation != "stale_runs":
        _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError):
        _verify(specimen, changed)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), 10**1000])
def test_nonfinite_numbers_rejected_before_comparison(specimen, value):
    changed = copy.deepcopy(specimen[0])
    changed["runs"][1]["spectral"]["gap_H"] = value
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError, match="finite"):
        _verify(specimen, changed)


@pytest.mark.parametrize("payload", [
    '{"schema":"first","schema":"second"}',
    '{"runs":[{"gap":NaN}]}',
    '{"runs":[{"gap":Infinity}]}',
    '{"runs":[{"gap":-Infinity}]}',
    '{"runs":[{"gap":1e1000}]}',
    '[]',
])
def test_strict_json_parser(tmp_path, payload):
    path = tmp_path / "receipt.json"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(verifier.ReceiptValidationError):
        verifier.load_receipt(path)


def test_default_verifier_does_not_accept_a_self_declared_smaller_grid(specimen):
    receipt, producer, _ = specimen
    with pytest.raises(verifier.ReceiptValidationError, match="requested replay grid"):
        verifier.verify_receipt(receipt, receipt, producer)


def test_source_binding_is_to_raw_bytes(specimen):
    receipt, producer, _ = specimen
    producer.write_bytes(producer.read_bytes().replace(b"\n", b"\r\n"))
    with pytest.raises(verifier.ReceiptValidationError, match="producer_sha256"):
        _verify(specimen, receipt)


@pytest.mark.parametrize("token", ["1e-1000", "-1e-1000", "2e-324", "-2e-324"])
def test_strict_loader_refuses_nonzero_json_numbers_lost_to_zero(tmp_path, token):
    path = tmp_path / "underflow.json"
    path.write_text('{"original_number":' + token + '}', encoding="utf-8")
    with pytest.raises(verifier.ReceiptValidationError, match="underflow|rounds to zero"):
        verifier.load_receipt(path)


@pytest.mark.parametrize("token", ["0", "-0.0", "0e-1000", "5e-324", "-5e-324", "1e-310"])
def test_strict_loader_keeps_zeros_and_representable_subnormals(tmp_path, token):
    path = tmp_path / "representable.json"
    path.write_text('{"original_number":' + token + '}', encoding="utf-8")
    value = verifier.load_receipt(path)["original_number"]
    assert value == float(token)
    if token in ("5e-324", "-5e-324", "1e-310"):
        assert value != 0


def _raw_ks_mass_token(receipt, token):
    changed = copy.deepcopy(receipt)
    local = next(run for run in changed["runs"] if run["transfer"] == "kogut_susskind")
    # Bind the attacker's anticipated rounded payload.  The raw token must
    # nevertheless be checked before it can become an exact-zero claim.
    local["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = json.loads(token)
    _rehash(changed)
    local["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] = "RAW_NUMBER_TOKEN"
    return json.dumps(changed, sort_keys=True).replace('"RAW_NUMBER_TOKEN"', token)


@pytest.mark.parametrize("token,accepted", [("0", True), ("-0.0", True), ("0e-1000", True),
                                            ("1e-1000", False), ("-1e-1000", False),
                                            ("5e-324", False)])
def test_cli_checks_original_json_before_structural_zero_claim(tmp_path, token, accepted):
    original = HERE / "receipts/z2_finite_transfer_receipt.json"
    path = tmp_path / "candidate.json"
    path.write_text(_raw_ks_mass_token(verifier.load_receipt(original), token), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(HERE / "verify_z2_finite_transfer_receipt.py"), str(path), str(original)],
        capture_output=True, text=True, check=False,
    )
    if accepted:
        assert result.returncode == 0, result.stderr
        assert "portable replay verified" in result.stdout
    else:
        assert result.returncode != 0, result.stdout
        assert "ReceiptValidationError" in result.stderr


@pytest.mark.parametrize("mutation", ["probability_lower", "probability_upper", "fit_ceiling",
                                      "influence_ceiling", "heat_bath_ceiling", "free_mass",
                                      "fiber_mass_lower", "fiber_mass_upper"])
def test_bindings_reject_mathematically_impossible_scalars(specimen, mutation):
    changed = copy.deepcopy(specimen[0])
    run = changed["runs"][1]
    if mutation == "probability_lower":
        run["spectral"].update(pi_min=0.03, pi_max=0.5)
        assert 0.5 + 31 * 0.03 > 1
    elif mutation == "probability_upper":
        run["spectral"].update(pi_min=0.001, pi_max=1 / 32)
        assert 0.001 + 31 / 32 < 1
    elif mutation == "fit_ceiling":
        run["constant_rate_fit"]["relative_frobenius_residual"] = 2.0
    elif mutation == "influence_ceiling":
        run["dobrushin"].update(eta_star=100.0, dobrushin_condition_holds=False,
                                unit_rate_floor_c_star_times_1_minus_eta=0.0)
    elif mutation == "heat_bath_ceiling":
        run["spectral"]["gap_unit_rate_heat_bath"] = 100.0
    elif mutation == "free_mass":
        changed["runs"][0]["fiber_dependent_rates"]["offdiagonal_mass_single_flip"] *= 0.5
    else:
        fibers = run["fiber_dependent_rates"]
        pairs = run["n_orbits"] * run["n_links"] / 2
        fibers["offdiagonal_mass_single_flip"] = (
            0.5 * pairs * fibers["rate_min"] if mutation == "fiber_mass_lower"
            else 2 * pairs * fibers["rate_max"]
        )
    _rehash(changed)
    with pytest.raises(verifier.ReceiptValidationError):
        verifier.verify_bindings(changed, specimen[1], expected_grid=specimen[2])
    # This is a necessary scalar condition, independently of whether two
    # claimed interacting replays happen to agree with each other.
    with pytest.raises(verifier.ReceiptValidationError):
        _verify(specimen, changed, changed)


@pytest.mark.parametrize("boundary", ["probability_lower", "probability_upper", "fit_ceiling",
                                      "influence_ceiling", "heat_bath_ceiling"])
def test_necessary_scalar_bounds_allow_dimension_roundoff(specimen, boundary):
    # These exercise only the necessary scalar bounds. Matching inputs here
    # are not an independent interacting eigensolve or a production receipt.
    changed = copy.deepcopy(specimen[0])
    run = changed["runs"][1]
    allowance = 128 * sys.float_info.epsilon * run["n_orbits"]
    if boundary == "probability_lower":
        run["spectral"].update(pi_min=0.03, pi_max=1 - 31 * 0.03 + allowance / 4)
    elif boundary == "probability_upper":
        run["spectral"].update(pi_min=0.001, pi_max=(1 - 0.001 - allowance / 4) / 31)
    elif boundary == "fit_ceiling":
        run["constant_rate_fit"]["relative_frobenius_residual"] = 1 + allowance / 4
    elif boundary == "influence_ceiling":
        run["dobrushin"].update(eta_star=7 * (1 + allowance / 4),
                                dobrushin_condition_holds=False,
                                unit_rate_floor_c_star_times_1_minus_eta=0.0)
    else:
        run["spectral"]["gap_unit_rate_heat_bath"] = 8 * (1 + allowance / 4)
    _rehash(changed)
    _verify(specimen, changed, changed)

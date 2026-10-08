#!/usr/bin/env python3
"""Tests for the source-facing Thomson spectral transport gate."""

from __future__ import annotations

import pytest
from fractions import Fraction

from thomson_spectral_transport import (
    build_source_transport_interval_certificate,
    blocked_missing_source_transport,
    validate_source_transport_payload,
)


def _budget() -> dict:
    return {"bound_interval": {"lo": "0", "hi": "1e-30"}}


def _source_measure() -> dict:
    return {
        "artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
        "finite_volume_levels": [
            {
                "ensemble_id": "ens0",
                "channel": "U1_Q",
                "levels": [{"level_id": "rho0", "s": "1.0", "energy": "1.0", "weight": "1.0"}],
            }
        ],
        "ward_projected_residues": [
            {"level_id": "rho0", "residue": "0.25", "current_normalization": "Q=T3+Y"}
        ],
        "current_normalization": "Q=T3+Y",
        "rho_had_or_measure": {
            "representation": "primitive_spectral_measure",
            "support_variable": "s",
            "pushforward_rule": "certified finite-volume continuum pushforward",
            "positivity_status": "certified_positive",
        },
        "transport_moment_certificate": {
            "kernel": "mZ(P)^2/(3*pi*s*(s+mZ(P)^2))",
            "Delta_had_image": {"lo": "4.30", "hi": "4.40"},
            "quadrature_error_bound": "1e-40",
            "tail_bound": "1e-40",
        },
        "systematics": {
            "statistical_budget": _budget(),
            "continuum_budget": _budget(),
            "finite_volume_budget": _budget(),
            "chiral_budget": _budget(),
            "current_matching_budget": _budget(),
            "quadrature_budget": _budget(),
            "endpoint_remainder_budget": _budget(),
        },
        "guards": {
            "stable_channel_only": False,
            "surrogate_hadron_artifact": False,
            "compare_only_external_endpoint": False,
        },
    }


def _valid_payload() -> dict:
    return {
        "artifact": "oph_source_ward_projected_thomson_transport",
        "source_only": True,
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "scheme": {
            "same_subtraction_as_a0": True,
            "scheme_id": "d10_u1q_thomson_v1",
            "normalization_convention": "Q=T3+Y",
        },
        "interval_backend": {
            "theorem_grade": True,
            "directed_outward_rounding": True,
            "library": "arb",
            "certificate": "arb-directed-outward-proof",
        },
        "p_interval": {"lo": "1.6308", "hi": "1.6311"},
        "endpoint_map": {
            "alpha_inv_image": {"lo": "136", "hi": "138"},
            "derivative_abs_bound": "5",
            "transport_error_bound": "1e-30",
            "components": {
                "a0_image": {"lo": "128.3", "hi": "128.4"},
                "Delta_lep_image": {"lo": "4.30", "hi": "4.31"},
                "Delta_had_image": {"lo": "4.30", "hi": "4.40"},
                "Delta_EW_image": {"lo": "0", "hi": "0.01"},
                "R_Q_image": {"lo": "0.03", "hi": "0.05"},
            },
        },
        "source_measure": _source_measure(),
        "delta_EW": {"zero_theorem": "declared_scheme_zero"},
        "fixed_point_certificate": {
            "self_map_pass": True,
            "uniqueness_pass": True,
            "G_image": {"lo": "1.6308", "hi": "1.6311"},
            "contraction_kappa": "0.001",
        },
    }


def test_missing_source_transport_blocks_promotion() -> None:
    result = blocked_missing_source_transport()

    assert result["status"] == "blocked_source_spectral_measure_missing"
    assert result["promotion_allowed"] is False
    assert "required_field_missing:source_measure" in result["reasons"]


def test_source_transport_rejects_hidden_alpha_compare_keys() -> None:
    payload = _valid_payload()
    payload["compare_alpha_inv"] = "137.035999177"

    result = validate_source_transport_payload(payload)

    assert result.promotion_allowed is False
    assert any(reason.startswith("forbidden_external_or_compare_key") for reason in result.reasons)


def test_source_transport_requires_scheme_and_fixed_point_certificate() -> None:
    payload = _valid_payload()
    payload["scheme"] = {"same_subtraction_as_a0": False}
    payload["fixed_point_certificate"]["uniqueness_pass"] = False

    result = validate_source_transport_payload(payload)

    assert result.promotion_allowed is False
    assert "scheme_not_locked_to_a0" in result.reasons
    assert "fixed_point_uniqueness_missing" in result.reasons


def test_complete_source_transport_contract_does_not_replay_a_certificate() -> None:
    result = validate_source_transport_payload(_valid_payload())

    assert result.status == "source_transport_contract_satisfied_unverified"
    assert result.contract_satisfied is True
    assert result.promotion_allowed is False
    assert result.reasons == ("source_certificate_replay_not_implemented",)


def test_boolean_only_source_transport_no_longer_promotes() -> None:
    payload = {
        "artifact": "oph_source_ward_projected_thomson_transport",
        "source_only": True,
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "scheme": {"same_subtraction_as_a0": True},
        "rho_had": {
            "positivity_certificate": True,
            "threshold_support": True,
            "ope_tail_certificate": True,
            "quadrature_error_bound": "1e-40",
        },
        "delta_EW": {"zero_theorem": "declared_scheme_zero"},
        "fixed_point_certificate": {"self_map_pass": True, "uniqueness_pass": True},
    }

    result = validate_source_transport_payload(payload)

    assert result.promotion_allowed is False
    assert "required_field_missing:source_measure" in result.reasons
    assert "required_field_missing:interval_backend" in result.reasons


def test_interval_certificate_builder_uses_source_payload_not_compare_alpha() -> None:
    certificate = build_source_transport_interval_certificate(_valid_payload())

    assert certificate["promotion_allowed"] is False
    assert certificate["contract_satisfied"] is True
    assert certificate["certificate_replayed"] is False
    assert certificate["external_inputs_used"] is False
    assert certificate["certified_intervals"]["alpha_interval"]["lo"].startswith("0.0072")
    assert "CODATA" not in str(certificate)


def test_interval_certificate_rejects_non_self_map() -> None:
    payload = _valid_payload()
    payload["p_interval"] = {"lo": "1.0", "hi": "1.1"}

    certificate = build_source_transport_interval_certificate(payload)

    assert certificate["promotion_allowed"] is False
    assert "fixed_point_self_map_failed" in certificate["reasons"]


def test_disjoint_component_sum_cannot_supply_an_endpoint():
    payload = _valid_payload()
    payload["endpoint_map"]["components"]["a0_image"] = {"lo": "1000", "hi": "1001"}
    result = validate_source_transport_payload(payload)
    assert not result.promotion_allowed
    assert not result.contract_satisfied
    assert "endpoint_components_disjoint_from_total" in result.reasons


@pytest.mark.parametrize("key,value,reason", [
    ("finite_volume_levels", ["nonsense"], "invalid_source_measure_finite_volume_levels"),
    ("ward_projected_residues", ["nonsense"], "invalid_source_measure_residue"),
    ("ward_projected_residues", [{"residue": "-1"}], "negative_decimal:source_measure.residue"),
])
def test_malformed_or_negative_spectral_data_fail_contract(key, value, reason):
    payload = _valid_payload()
    payload["source_measure"][key] = value
    result = validate_source_transport_payload(payload)
    assert not result.promotion_allowed
    assert not result.contract_satisfied
    assert reason in result.reasons


@pytest.mark.parametrize("label", ["not_positive", "unproved", "not_certified", "positive_but_not_proved"])
def test_negated_positivity_labels_are_not_certificates(label):
    payload = _valid_payload()
    payload["source_measure"]["rho_had_or_measure"]["positivity_status"] = label
    result = validate_source_transport_payload(payload)
    assert not result.promotion_allowed
    assert not result.contract_satisfied
    assert "source_measure_positivity_not_certified" in result.reasons


def test_numeric_zero_derivative_is_a_valid_contract_bound():
    payload = _valid_payload()
    payload["endpoint_map"]["derivative_abs_bound"] = 0
    result = validate_source_transport_payload(payload)
    assert result.contract_satisfied
    assert not result.promotion_allowed


def test_user_metadata_cannot_enable_backend_replay():
    payload = _valid_payload()
    payload["interval_backend"].update({"certificate_replayed": True, "verified": True})
    result = validate_source_transport_payload(payload)
    assert result.contract_satisfied
    assert not result.promotion_allowed
    assert result.to_json()["certificate_replayed"] is False


def test_reciprocal_endpoints_enclose_original_exact_rationals():
    certificate = build_source_transport_interval_certificate(_valid_payload(), precision=20)
    interval = certificate["certified_intervals"]["alpha_interval"]
    assert Fraction(interval["lo"]) <= Fraction(1, 138)
    assert Fraction(interval["hi"]) >= Fraction(1, 136)


def test_hadronic_moment_must_overlap_its_endpoint_component():
    payload = _valid_payload()
    payload["source_measure"]["transport_moment_certificate"]["Delta_had_image"] = {"lo": "9", "hi": "10"}
    result = validate_source_transport_payload(payload)
    assert not result.contract_satisfied
    assert "hadronic_moment_disjoint_from_endpoint_component" in result.reasons


@pytest.mark.parametrize("moment", [None, "proof", []])
def test_malformed_moment_returns_a_blocked_contract(moment):
    payload = _valid_payload()
    payload["source_measure"]["transport_moment_certificate"] = moment
    result = validate_source_transport_payload(payload)
    assert not result.promotion_allowed
    assert not result.contract_satisfied
    assert "required_field_missing:source_measure.transport_moment_certificate" in result.reasons


def test_correlated_endpoint_may_be_narrower_than_the_component_interval_sum():
    payload = _valid_payload()
    payload["endpoint_map"]["alpha_inv_image"] = {"lo": "137", "hi": "137"}
    result = validate_source_transport_payload(payload)
    assert result.contract_satisfied
    assert not result.promotion_allowed

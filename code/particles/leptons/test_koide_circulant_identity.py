"""Tests for the exact C3-circulant Koide identity receipt."""

from __future__ import annotations

import json
from decimal import Decimal, localcontext

import derive_koide_circulant_identity as lane


def test_exact_balance_equivalence() -> None:
    balance = Decimal(1) / Decimal(2).sqrt()
    assert (
        abs(lane.koide_from_modulus_ratio(balance) - Decimal(2) / Decimal(3))
        < Decimal("1e-27")
    )
    assert (
        abs(
            lane.modulus_ratio_from_koide(Decimal(2) / Decimal(3))
            - balance
        )
        < Decimal("1e-27")
    )


def test_phase_independence_has_a_positive_chamber_boundary() -> None:
    balance = 2.0**-0.5
    for phase in (0.0, 0.1, 2.0 / 9.0):
        values = lane.roots(1.0, balance, phase)
        assert min(values) >= 0.0
        assert abs(lane.physical_koide(values) - 2.0 / 3.0) < 1.0e-14
    outside = lane.roots(1.0, balance, 0.4)
    assert min(outside) < 0.0
    assert abs(lane.physical_koide(outside) - 2.0 / 3.0) > 1.0e-3


def test_artifact_separates_identity_model_and_comparison() -> None:
    artifact = lane.build_artifact()
    assert artifact["checks_pass"] is True
    assert artifact["source_only_physical_prediction"] is False
    assert artifact["public_physical_promotion_allowed"] is False
    assert (
        artifact["identity"]["equivalence"]
        == "Q = 2/3 iff |b|/a = 1/sqrt(2)"
    )
    assert (
        artifact["compare_only_PDG_2026_central_coordinate"][
            "significance_claimed"
        ]
        is False
    )
    assert artifact["compare_only_PDG_2026_central_coordinate"][
        "masses_MeV"
    ] == lane.EXPECTED_PDG_2026_MASSES_MEV
    assert Decimal(
        artifact["compare_only_PDG_2026_central_coordinate"]["Q"]
    ) == lane.EXPECTED_PDG_2026_Q
    assert artifact["checks"]["PDG_2026_mass_coordinate_is_exact"] is True
    assert "historically target-informed" in artifact[
        "current_MCPR_coordinate"
    ]["provenance"]


def test_two_scale_free_ratios_use_displayed_uncertainties() -> None:
    comparison = lane.build_artifact()["compare_only_PDG_2026_ratio_envelope"]
    rows = comparison["independent_ratios"]
    assert set(rows) == {"m_mu_over_m_e", "m_tau_over_m_e"}
    mu = rows["m_mu_over_m_e"]
    tau = rows["m_tau_over_m_e"]
    assert Decimal("0.086") < Decimal(mu["relative_central_residual_ppm"]) < Decimal(
        "0.087"
    )
    assert Decimal("0.022") < Decimal(
        mu["conservative_component_interval_halfwidth_ppm"]
    ) < Decimal("0.023")
    assert mu["inside_componentwise_one_sigma_interval"] is False
    assert tau["inside_componentwise_one_sigma_interval"] is True
    assert (
        comparison[
            "common_scale_under_identity_mass_scheme_necessary_condition_passed"
        ]
        is False
    )
    assert "no covariance or joint confidence" in comparison["interval_basis"]

    mcpr = json.loads(lane.MCPR.read_text(encoding="utf-8"))[
        "conditional_prediction"
    ]
    masses = [Decimal(value) for value in mcpr["masses_over_E_star"]]
    with localcontext() as context:
        context.prec = 70
        assert abs(
            masses[1] / masses[0] - Decimal(mcpr["ratios"]["m_mu_over_m_e"])
        ) < Decimal("1e-54")
        assert abs(
            masses[2] / masses[0] - Decimal(mcpr["ratios"]["m_tau_over_m_e"])
        ) < Decimal("1e-54")


def test_ratio_interval_mutation_fails_closed_without_changing_candidate() -> None:
    mcpr = json.loads(lane.MCPR.read_text(encoding="utf-8"))
    reference = json.loads(lane.PDG_REFERENCES.read_text(encoding="utf-8"))
    original = lane.compare_independent_ratios(mcpr, reference)
    assert (
        original["independent_ratios"]["m_mu_over_m_e"][
            "inside_componentwise_one_sigma_interval"
        ]
        is False
    )

    # A wider reported muon uncertainty changes this componentwise comparison
    # while preserving the frozen candidate and the tau ratio.
    reference["entries"]["muon"]["value_text"] = "105.6583755+-0.0000100"
    widened = lane.compare_independent_ratios(mcpr, reference)
    assert (
        widened["independent_ratios"]["m_mu_over_m_e"][
            "inside_componentwise_one_sigma_interval"
        ]
        is True
    )
    assert widened["independent_ratios"]["m_tau_over_m_e"] == original[
        "independent_ratios"
    ]["m_tau_over_m_e"]


def test_committed_artifact_is_byte_exact() -> None:
    artifact = lane.build_artifact()
    assert lane.DEFAULT_OUT.read_bytes() == lane.canonical_bytes(artifact)
    stored = json.loads(lane.DEFAULT_OUT.read_text(encoding="utf-8"))
    assert stored["artifact"] == "oph_koide_circulant_identity"

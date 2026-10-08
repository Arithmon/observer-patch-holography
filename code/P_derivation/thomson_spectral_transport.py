#!/usr/bin/env python3
"""Source-facing Ward-projected Thomson endpoint validator.

This module is deliberately a gate, not a fitted endpoint solver.  It accepts a
declared spectral transport contract only if the payload carries the same
D10 family/scheme identifiers as the anchor, declares a Ward-projected
hadronic spectral measure, supplies an EW remainder statement, and keeps all
external fine-structure comparison keys out of the source path. Contract
consistency is not verification of its source data or interval proof. This
module has no backend-certificate replay, so it cannot authorize promotion.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from fractions import Fraction
from typing import Any

from interval_backend import Interval
from paper_math import decimal_pi


FORBIDDEN_SOURCE_KEYS = {
    "CODATA",
    "NIST",
    "codata",
    "nist",
    "compare_alpha",
    "compare_alpha_inv",
    "measured_alpha",
    "measured_alpha_inv",
    "P_C",
    "S_required",
    "c_Q",
    "residual_second_order_coefficient",
    "missing_source_transport_delta_alpha_inv",
    "required_transport_delta_alpha_inv",
}

SOURCE_SYSTEMATICS_BUDGETS = (
    "statistical_budget",
    "continuum_budget",
    "finite_volume_budget",
    "chiral_budget",
    "current_matching_budget",
    "quadrature_budget",
    "endpoint_remainder_budget",
)

ENDPOINT_COMPONENT_INTERVALS = (
    "a0_image",
    "Delta_lep_image",
    "Delta_had_image",
    "Delta_EW_image",
    "R_Q_image",
)

THEOREM_GRADE_REQUIRED_FIELDS = (
    "artifact",
    "source_only",
    "source_family_id",
    "current",
    "scheme.same_subtraction_as_a0",
    "scheme.scheme_id",
    "scheme.normalization_convention",
    "interval_backend.theorem_grade",
    "interval_backend.directed_outward_rounding",
    "interval_backend.library",
    "interval_backend.certificate",
    "p_interval",
    "endpoint_map.alpha_inv_image",
    "endpoint_map.derivative_abs_bound",
    "endpoint_map.components.a0_image",
    "endpoint_map.components.Delta_lep_image",
    "endpoint_map.components.Delta_had_image",
    "endpoint_map.components.Delta_EW_image",
    "endpoint_map.components.R_Q_image",
    "endpoint_map.transport_error_bound",
    "fixed_point_certificate.G_image",
    "fixed_point_certificate.contraction_kappa",
    "source_measure.artifact",
    "source_measure.finite_volume_levels",
    "source_measure.ward_projected_residues",
    "source_measure.current_normalization",
    "source_measure.rho_had_or_measure",
    "source_measure.transport_moment_certificate",
    "source_measure.systematics.*.bound_interval",
    "source_measure.guards",
    "delta_EW.zero_theorem_or_source_bound",
)


@dataclass(frozen=True)
class TransportValidation:
    status: str
    promotion_allowed: bool
    reasons: tuple[str, ...]
    contract_satisfied: bool = False
    interval_certificate: dict[str, Any] | None = None
    required_fields: tuple[str, ...] = THEOREM_GRADE_REQUIRED_FIELDS

    def to_json(self) -> dict[str, Any]:
        payload = {
            "status": self.status,
            "promotion_allowed": self.promotion_allowed,
            "reasons": list(self.reasons),
            "contract_satisfied": self.contract_satisfied,
            "certificate_replayed": False,
            "required_fields": list(self.required_fields),
        }
        if self.interval_certificate is not None:
            payload["interval_certificate"] = self.interval_certificate
        return payload


def _dedupe(items: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(items))


def _decimal_to_json(value: Decimal) -> str:
    return format(value, "f")


def _interval_to_json(value: Interval) -> dict[str, str]:
    return {"lo": _decimal_to_json(value.lo), "hi": _decimal_to_json(value.hi)}


def _parse_decimal(value: Any, field: str, reasons: list[str]) -> Decimal | None:
    if value is None:
        reasons.append(f"required_field_missing:{field}")
        return None
    try:
        parsed = Decimal(str(value))
    except (InvalidOperation, ValueError):
        reasons.append(f"invalid_decimal:{field}")
        return None
    if not parsed.is_finite():
        reasons.append(f"nonfinite_decimal:{field}")
        return None
    return parsed


def _parse_nonnegative_decimal(value: Any, field: str, reasons: list[str]) -> Decimal | None:
    parsed = _parse_decimal(value, field, reasons)
    if parsed is not None and parsed < 0:
        reasons.append(f"negative_decimal:{field}")
    return parsed


def _parse_interval(value: Any, field: str, reasons: list[str], *, positive: bool = False) -> Interval | None:
    if not isinstance(value, dict):
        reasons.append(f"required_interval_missing:{field}")
        return None
    if "lo" not in value or "hi" not in value:
        reasons.append(f"required_interval_endpoint_missing:{field}")
        return None
    lo = _parse_decimal(value["lo"], f"{field}.lo", reasons)
    hi = _parse_decimal(value["hi"], f"{field}.hi", reasons)
    if lo is None or hi is None:
        return None
    if lo > hi:
        reasons.append(f"invalid_interval_order:{field}")
        return None
    if positive and lo <= 0:
        reasons.append(f"nonpositive_interval:{field}")
        return None
    return Interval(lo, hi)


def _has_bound_interval(payload: Any) -> bool:
    return isinstance(payload, dict) and isinstance(payload.get("bound_interval"), dict)


def _walk_keys(payload: Any, prefix: str = "") -> set[str]:
    keys: set[str] = set()
    if isinstance(payload, dict):
        for key, value in payload.items():
            key_text = str(key)
            full = f"{prefix}.{key_text}" if prefix else key_text
            keys.add(key_text)
            keys.add(full)
            keys.update(_walk_keys(value, full))
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            keys.update(_walk_keys(value, f"{prefix}[{index}]"))
    return keys


def source_payload_forbidden_keys(payload: dict[str, Any]) -> set[str]:
    keys = _walk_keys(payload)
    return {key for key in keys if key in FORBIDDEN_SOURCE_KEYS}


def _validate_interval_backend(payload: dict[str, Any], reasons: list[str]) -> None:
    backend = payload.get("interval_backend")
    if not isinstance(backend, dict):
        reasons.append("required_field_missing:interval_backend")
        return

    if backend.get("theorem_grade") is not True:
        reasons.append("interval_backend_not_theorem_grade")
    if backend.get("directed_outward_rounding") is not True:
        reasons.append("interval_backend_not_directed_outward")

    library = str(backend.get("library", "")).strip()
    if not library:
        reasons.append("required_field_missing:interval_backend.library")
    elif library.lower() in {"decimal", "decimalintervalbackend", "decimal_interval_backend"}:
        reasons.append("interval_backend_library_not_certified")

    if not (backend.get("certificate") or backend.get("proof_artifact")):
        reasons.append("required_field_missing:interval_backend.certificate")


def _validate_source_measure(payload: dict[str, Any], reasons: list[str]) -> None:
    measure = payload.get("source_measure") or payload.get("spectral_measure")
    if not isinstance(measure, dict):
        reasons.append("required_field_missing:source_measure")
        return

    if measure.get("artifact") != "oph_qcd_ward_projected_hadronic_spectral_measure":
        reasons.append("source_measure_artifact_mismatch")

    if not isinstance(measure.get("finite_volume_levels"), list) or not measure["finite_volume_levels"]:
        reasons.append("required_field_missing:source_measure.finite_volume_levels")
    else:
        for ensemble in measure["finite_volume_levels"]:
            if not isinstance(ensemble, dict) or not isinstance(ensemble.get("levels"), list) or not ensemble["levels"]:
                reasons.append("invalid_source_measure_finite_volume_levels")
                continue
            for level in ensemble["levels"]:
                if not isinstance(level, dict):
                    reasons.append("invalid_source_measure_level")
                    continue
                support = _parse_decimal(level.get("s"), "source_measure.level.s", reasons)
                energy = _parse_decimal(level.get("energy"), "source_measure.level.energy", reasons)
                _parse_nonnegative_decimal(level.get("weight"), "source_measure.level.weight", reasons)
                if support is not None and energy is not None:
                    if support <= 0 or energy <= 0:
                        reasons.append("source_measure_level_support_not_positive")
                    elif Fraction(support) != Fraction(energy) ** 2:
                        reasons.append("source_measure_level_s_energy_inconsistent")

    residues = measure.get("ward_projected_residues")
    if not isinstance(residues, list) or not residues:
        reasons.append("required_field_missing:source_measure.ward_projected_residues")
    else:
        for residue in residues:
            if not isinstance(residue, dict):
                reasons.append("invalid_source_measure_residue")
                continue
            _parse_nonnegative_decimal(residue.get("residue"), "source_measure.residue", reasons)

    has_current_normalization = bool(measure.get("current_normalization"))
    if isinstance(residues, list) and residues:
        has_current_normalization = has_current_normalization or all(
            isinstance(item, dict) and bool(item.get("current_normalization")) for item in residues
        )
    if not has_current_normalization:
        reasons.append("required_field_missing:source_measure.current_normalization")

    rho = measure.get("rho_had_or_measure")
    if not isinstance(rho, dict):
        reasons.append("required_field_missing:source_measure.rho_had_or_measure")
    else:
        if rho.get("support_variable") != "s":
            reasons.append("source_measure_support_variable_mismatch")
        positivity = str(rho.get("positivity_status", "")).lower()
        if positivity not in {"certified_positive", "proved_positive"}:
            reasons.append("source_measure_positivity_not_certified")
        if not rho.get("pushforward_rule"):
            reasons.append("required_field_missing:source_measure.rho_had_or_measure.pushforward_rule")

    moment = measure.get("transport_moment_certificate")
    if not isinstance(moment, dict):
        reasons.append("required_field_missing:source_measure.transport_moment_certificate")
    else:
        if not moment.get("kernel"):
            reasons.append("required_field_missing:source_measure.transport_moment_certificate.kernel")
        _parse_interval(
            moment.get("Delta_had_image"),
            "source_measure.transport_moment_certificate.Delta_had_image",
            reasons,
        )
        _parse_nonnegative_decimal(
            moment.get("quadrature_error_bound"),
            "source_measure.transport_moment_certificate.quadrature_error_bound",
            reasons,
        )
        _parse_nonnegative_decimal(
            moment.get("tail_bound"),
            "source_measure.transport_moment_certificate.tail_bound",
            reasons,
        )

    systematics = measure.get("systematics")
    if not isinstance(systematics, dict):
        reasons.append("required_field_missing:source_measure.systematics")
    else:
        for budget_name in SOURCE_SYSTEMATICS_BUDGETS:
            budget = systematics.get(budget_name)
            if not _has_bound_interval(budget):
                reasons.append(f"required_field_missing:source_measure.systematics.{budget_name}.bound_interval")
            else:
                interval = _parse_interval(
                    budget["bound_interval"],
                    f"source_measure.systematics.{budget_name}.bound_interval",
                    reasons,
                )
                if interval is not None and interval.lo < 0:
                    reasons.append(f"negative_systematics_budget:{budget_name}")

    guards = measure.get("guards")
    if not isinstance(guards, dict):
        reasons.append("required_field_missing:source_measure.guards")
    else:
        if guards.get("stable_channel_only") is not False:
            reasons.append("source_measure_guard_stable_channel_only_not_false")
        if guards.get("surrogate_hadron_artifact") is not False:
            reasons.append("source_measure_guard_surrogate_not_false")
        if guards.get("compare_only_external_endpoint") is not False:
            reasons.append("source_measure_guard_compare_endpoint_not_false")


def _validate_delta_ew(payload: dict[str, Any], reasons: list[str]) -> None:
    delta_ew = payload.get("delta_EW")
    if not isinstance(delta_ew, dict):
        reasons.append("delta_EW_missing")
        return

    if delta_ew.get("zero_theorem"):
        return

    source_bound = delta_ew.get("source_bound")
    if source_bound is None:
        reasons.append("delta_EW_zero_or_bound_missing")
        return
    _parse_interval(source_bound, "delta_EW.source_bound", reasons)


def build_source_transport_interval_certificate(
    payload: dict[str, Any],
    *,
    precision: int = 96,
) -> dict[str, Any]:
    """Check a declared interval contract and compute Decimal diagnostics.

    Backend names, certificate strings and asserted derivative bounds are
    unverified inputs. The local Decimal calculation is not directed interval
    arithmetic. Even a consistent contract does not establish a theorem.
    """
    reasons: list[str] = []
    forbidden = sorted(source_payload_forbidden_keys(payload))
    if forbidden:
        reasons.append("forbidden_external_or_compare_key:" + ",".join(forbidden))

    _validate_interval_backend(payload, reasons)

    p_interval = _parse_interval(payload.get("p_interval"), "p_interval", reasons, positive=True)

    endpoint = payload.get("endpoint_map")
    if not isinstance(endpoint, dict):
        reasons.append("required_field_missing:endpoint_map")
        endpoint = {}

    alpha_inv_image = _parse_interval(
        endpoint.get("alpha_inv_image") or endpoint.get("A_T_image"),
        "endpoint_map.alpha_inv_image",
        reasons,
        positive=True,
    )
    derivative_bound = _parse_nonnegative_decimal(
        endpoint.get("derivative_abs_bound", endpoint.get("A_T_prime_abs_bound")),
        "endpoint_map.derivative_abs_bound",
        reasons,
    )
    _parse_nonnegative_decimal(endpoint.get("transport_error_bound"), "endpoint_map.transport_error_bound", reasons)

    components = endpoint.get("components")
    parsed_components = {}
    if not isinstance(components, dict):
        reasons.append("required_field_missing:endpoint_map.components")
    else:
        for field in ENDPOINT_COMPONENT_INTERVALS:
            parsed_components[field] = _parse_interval(components.get(field), f"endpoint_map.components.{field}", reasons)
        summands = [parsed_components.get(field) for field in ENDPOINT_COMPONENT_INTERVALS[:4]]
        if alpha_inv_image is not None and all(value is not None for value in summands):
            # Exact rational comparison avoids a rounded sum accepting disjoint
            # enclosures. Overlap is necessary, not sufficient: correlations
            # can make a valid total enclosure narrower than the interval sum.
            lower = sum(Fraction(value.lo) for value in summands)
            upper = sum(Fraction(value.hi) for value in summands)
            if lower > Fraction(alpha_inv_image.hi) or upper < Fraction(alpha_inv_image.lo):
                reasons.append("endpoint_components_disjoint_from_total")
        measure = payload.get("source_measure") or payload.get("spectral_measure") or {}
        moment = measure.get("transport_moment_certificate", {}) if isinstance(measure, dict) else {}
        if not isinstance(moment, dict):
            moment = {}
        hadronic = _parse_interval(moment.get("Delta_had_image"),
                                   "source_measure.transport_moment_certificate.Delta_had_image", reasons)
        hadronic_component = parsed_components.get("Delta_had_image")
        if any(value is not None and value.hi < 0 for value in (hadronic, hadronic_component)):
            reasons.append("positive_hadronic_measure_has_wholly_negative_moment")
        if hadronic is not None and hadronic_component is not None:
            if hadronic.lo > hadronic_component.hi or hadronic.hi < hadronic_component.lo:
                reasons.append("hadronic_moment_disjoint_from_endpoint_component")
        ew = payload.get("delta_EW")
        ew_component = parsed_components.get("Delta_EW_image")
        if isinstance(ew, dict):
            ew_bound = None
            if ew.get("source_bound") is not None:
                ew_bound = _parse_interval(ew["source_bound"], "delta_EW.source_bound", reasons)
            if ew.get("zero_theorem"):
                if any(value is not None and not value.contains(0) for value in (ew_bound, ew_component)):
                    reasons.append("delta_EW_zero_theorem_disjoint_from_bound")
            if ew_bound is not None and ew_component is not None:
                if ew_bound.lo > ew_component.hi or ew_bound.hi < ew_component.lo:
                    reasons.append("delta_EW_source_bound_disjoint_from_endpoint_component")

    fixed_point = payload.get("fixed_point_certificate")
    if not isinstance(fixed_point, dict):
        reasons.append("fixed_point_certificate_missing")
        fixed_point = {}

    claimed_g_image = _parse_interval(
        fixed_point.get("G_image"),
        "fixed_point_certificate.G_image",
        reasons,
        positive=True,
    )
    claimed_kappa = _parse_nonnegative_decimal(
        fixed_point.get("contraction_kappa"),
        "fixed_point_certificate.contraction_kappa",
        reasons,
    )

    computed_g_image: Interval | None = None
    computed_kappa: Decimal | None = None
    alpha_interval: Interval | None = None
    if p_interval is not None and alpha_inv_image is not None and derivative_bound is not None:
        with localcontext() as ctx:
            ctx.prec = precision
            pi = decimal_pi(precision + 8)
            sqrt_pi = +pi.sqrt()
            phi = +((Decimal(1) + Decimal(5).sqrt()) / Decimal(2))
            computed_g_image = Interval(
                +(phi + sqrt_pi / alpha_inv_image.hi),
                +(phi + sqrt_pi / alpha_inv_image.lo),
            )
            computed_kappa = +(sqrt_pi * derivative_bound / (alpha_inv_image.lo * alpha_inv_image.lo))
            # These reciprocal endpoints can be directed exactly with Decimal;
            # the separate pi/sqrt G calculation above remains only an audit.
            ctx.rounding = "ROUND_FLOOR"
            alpha_lo = Decimal(1) / alpha_inv_image.hi
            ctx.rounding = "ROUND_CEILING"
            alpha_hi = Decimal(1) / alpha_inv_image.lo
            alpha_interval = Interval(alpha_lo, alpha_hi)

    if computed_g_image is not None and claimed_g_image is not None:
        if not computed_g_image.subset_of(claimed_g_image):
            reasons.append("fixed_point_G_image_does_not_enclose_outer_map")
        if p_interval is not None and not claimed_g_image.subset_of(p_interval):
            reasons.append("fixed_point_self_map_failed")

    if computed_kappa is not None and claimed_kappa is not None:
        if claimed_kappa >= 1:
            reasons.append("fixed_point_uniqueness_failed")
        if computed_kappa > claimed_kappa:
            reasons.append("fixed_point_contraction_kappa_too_small")

    if fixed_point.get("self_map_pass") is not True:
        reasons.append("fixed_point_self_map_missing")
    if fixed_point.get("uniqueness_pass") is not True:
        reasons.append("fixed_point_uniqueness_missing")

    contract_satisfied = not reasons
    # No implementation here reads or replays the supplied backend certificate.
    # Do not turn a nonempty metadata label into a physical-alpha theorem.
    reasons.append("source_certificate_replay_not_implemented")
    return {
        "artifact": "oph_source_transport_interval_certificate",
        "status": (
            "source_interval_contract_satisfied_unverified"
            if contract_satisfied
            else "blocked_source_interval_certificate_failed"
        ),
        "promotion_allowed": False,
        "contract_satisfied": contract_satisfied,
        "certificate_replayed": False,
        "verification_scope": "declared interval consistency and nondirected Decimal diagnostics only",
        "certified_intervals_status": "source bounds unverified; only reciprocal alpha endpoints are directed locally",
        "external_inputs_used": False,
        "reasons": list(_dedupe(reasons)),
        "theorem": {
            "map": "G(P)=phi+sqrt(pi)/A_T(P)",
            "source_endpoint": "A_T(P)=a0(P)+Delta_lep_src(P)+Delta_had_src(P)+Delta_EW_src(P)",
            "self_map_condition": "G(I_P) subset I_P",
            "uniqueness_condition": "sqrt(pi) * sup|A_T'(P)| / inf(A_T(I_P))^2 < 1",
            "alpha_interval_consequence": "alpha_* in 1/A_T(I_P)",
        },
        "certified_intervals": {
            "I_P": _interval_to_json(p_interval) if p_interval is not None else None,
            "A_T_image": _interval_to_json(alpha_inv_image) if alpha_inv_image is not None else None,
            "claimed_G_image": _interval_to_json(claimed_g_image) if claimed_g_image is not None else None,
            "decimal_audit_G_image": _interval_to_json(computed_g_image) if computed_g_image is not None else None,
            "alpha_interval": _interval_to_json(alpha_interval) if alpha_interval is not None else None,
        },
        "bounds": {
            "A_T_prime_abs_bound": _decimal_to_json(derivative_bound) if derivative_bound is not None else None,
            "claimed_contraction_kappa": _decimal_to_json(claimed_kappa) if claimed_kappa is not None else None,
            "decimal_audit_contraction_kappa": (
                _decimal_to_json(computed_kappa) if computed_kappa is not None else None
            ),
        },
        "required_fields": list(THEOREM_GRADE_REQUIRED_FIELDS),
    }


def validate_source_transport_payload(payload: dict[str, Any]) -> TransportValidation:
    reasons: list[str] = []

    forbidden = sorted(source_payload_forbidden_keys(payload))
    if forbidden:
        reasons.append("forbidden_external_or_compare_key:" + ",".join(forbidden))

    if payload.get("artifact") != "oph_source_ward_projected_thomson_transport":
        reasons.append("artifact_mismatch")

    if payload.get("source_only") is not True:
        reasons.append("source_only_false_or_missing")

    if payload.get("source_family_id") != "d10_running_tree":
        reasons.append("source_family_id_mismatch")

    if payload.get("current") != "U1_Q":
        reasons.append("ward_current_mismatch")

    scheme = payload.get("scheme", {})
    if not isinstance(scheme, dict) or scheme.get("same_subtraction_as_a0") is not True:
        reasons.append("scheme_not_locked_to_a0")
    else:
        if not scheme.get("scheme_id"):
            reasons.append("required_field_missing:scheme.scheme_id")
        if not scheme.get("normalization_convention"):
            reasons.append("required_field_missing:scheme.normalization_convention")

    _validate_source_measure(payload, reasons)
    _validate_delta_ew(payload, reasons)
    interval_certificate = build_source_transport_interval_certificate(payload)
    reasons.extend(reason for reason in interval_certificate["reasons"]
                   if reason != "source_certificate_replay_not_implemented")
    reasons_tuple = _dedupe(reasons)

    if reasons_tuple:
        return TransportValidation(
            status="blocked_source_transport_contract_failed",
            promotion_allowed=False,
            reasons=reasons_tuple,
            interval_certificate=interval_certificate,
        )

    return TransportValidation(
        status="source_transport_contract_satisfied_unverified",
        promotion_allowed=False,
        contract_satisfied=True,
        reasons=("source_certificate_replay_not_implemented",),
        interval_certificate=interval_certificate,
    )


def blocked_missing_source_transport() -> dict[str, Any]:
    return TransportValidation(
        status="blocked_source_spectral_measure_missing",
        promotion_allowed=False,
        reasons=(
            "required_field_missing:source_measure",
            "required_field_missing:interval_backend",
            "required_field_missing:endpoint_map",
            "delta_EW_missing",
            "fixed_point_certificate_missing",
        ),
    ).to_json()

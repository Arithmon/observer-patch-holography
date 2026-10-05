from __future__ import annotations

from fractions import Fraction

import pytest

import carrier_scale_bound_diagnostic as cb


def test_receipt_builds_with_expected_status() -> None:
    receipt = cb.finalize(cb.build_receipt())
    assert receipt["status"] == (
        "EXPOSED_RETROSPECTIVE_CARRIER_SCALE_BOUND__DIAGNOSTIC_ONLY"
    )
    boundary = receipt["comparison_boundary"]
    assert boundary["public_measurement_read"] is True
    assert boundary["comparison_permitted"] is False
    assert boundary["scored"] is False
    assert boundary["fz11_untouched"] is True


def test_bound_value_matches_exact_identification() -> None:
    receipt = cb.finalize(cb.build_receipt())
    assert receipt["bound"]["carrier_scale_upper_bound_m"] == "8.825e-29"
    assert receipt["bound"]["planck_length_headroom"] == "5.460e+06"
    # independent recomputation with floats
    import math

    a_gev = math.sqrt(20) / 1.0e13
    a_m = a_gev * 1.97327e-16
    assert abs(a_m - 8.825e-29) / a_m < 1e-3
    assert abs(a_m / 1.616255e-35 - 5.460e6) / 5.460e6 < 1e-3


def test_unrestricted_radii_constrain_effective_scale_only() -> None:
    import sympy as sp

    # Independent raw icosahedron moments along z; no class-certificate import.
    phi = (1 + sp.sqrt(5)) / 2
    unit_z = [0] * 4 + [sign / sp.sqrt(1 + phi**2) for sign in (-1, 1)] * 2
    unit_z += [sign * phi / sp.sqrt(1 + phi**2) for sign in (-1, 1)] * 2
    radius = sp.Symbol("radius", positive=True)
    k2 = sum((radius * z) ** 2 / 2 for z in unit_z)
    k4 = -sum((radius * z) ** 4 / 24 for z in unit_z)
    assert sp.simplify(k4 / k2) == -radius**2 / 20

    # Continuum normalization leaves arbitrary r^2 in C4/a^2. Even pure
    # single-radius members can move the inferred a ceiling by six decades.
    for r in (sp.Rational(1, 10**6), sp.Integer(10**6)):
        coefficient = sp.simplify((k4 / k2).subs(radius, r))
        assert sp.sqrt((-sp.Rational(1, 20)) / coefficient) == 1 / r

    # The entire cosine symbol, not merely its leading coefficient, has the
    # exact rescaling degeneracy (a, r) -> (s a, r/s).
    a, k, scale = sp.symbols("a k scale", positive=True)
    symbol = sum(1 - sp.cos(a * k * radius * z) for z in unit_z) / (a**2 * k2)
    rescaled = symbol.subs({a: scale * a, radius: radius / scale}, simultaneous=True)
    assert sp.simplify(rescaled - symbol) == 0

    receipt = cb.finalize(cb.build_receipt())
    general = receipt["translation"]["general_positive_weight_class"]
    assert general["leading_coefficient"] == "C4 = -(a^2/20)(mu4/mu2)"
    assert general["constrained_length"] == "a_eff = a sqrt(mu4/mu2)"
    assert general["single_radius"] == "a_eff = a r"
    assert general["bound_on_a_without_radius_contract"] is False
    assert receipt["bound"]["carrier_scale_upper_bound_m"] == "8.825e-29"


def test_no_linear_term_statement_is_true_of_the_stencil() -> None:
    # the stencil expansion carries even powers of (a k . u) only
    import a5_multipole_fixed_point_certificate as base

    verts = base.cartesian_vertices()
    for k in (1, 3, 5, 7):
        assert base.p_is_zero(base.moment_sum(verts, k))


def test_committed_receipt_is_byte_exact() -> None:
    committed = cb.RECEIPT_PATH.read_bytes()
    rebuilt = cb.finalize(cb.build_receipt())
    assert committed == cb.base.canonical_json_bytes(rebuilt)


def test_parent_pin_drift_fails_closed(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    import json

    tampered = json.loads(cb.base.RECEIPT_PATH.read_text())
    tampered["kinetic_stencil_conditional"]["expansion"] = "k^2 - (a^2/10) k^4"
    path = tmp_path / "parent.json"
    path.write_text(json.dumps(tampered))
    monkeypatch.setattr(cb.base, "RECEIPT_PATH", path)
    with pytest.raises(cb.base.FingerprintError):
        cb.build_receipt()

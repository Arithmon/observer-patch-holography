#!/usr/bin/env python3
"""Diagnose the frozen conditional light-quark ratio without selecting a law.

This is a downstream sensitivity calculation.  It never changes the RSCC
producer or uses reference values to choose coefficients.  Each ablation
removes exactly one already declared quadratic term while retaining all other
terms; it is not a theory-error envelope or an alternative prediction.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / "code/particles/runs/flavor/conditional_quark_mass_replay.json"
REFERENCE = ROOT / "code/particles/data/pdg_2026_conditional_quark_comparison.json"
VERIFIER = Path(__file__).with_name("verify_conditional_quark_mass_replay.py")
OUTPUT = ROOT / "code/particles/runs/flavor/conditional_quark_ratio_diagnostic.json"
D = Decimal


def number(value: object) -> Decimal:
    result = D(str(value))
    if not result.is_finite():
        raise ValueError("finite decimal required")
    return result


def fmt(value: Decimal, places: int = 12) -> str:
    return format(value, f".{places}f")


def native_diagnostic(packet: dict) -> dict:
    """Reconstruct ratios from the frozen exponents, with no comparison data."""
    if (packet["schema"] != "oph.conditional_quark_mass_replay.v1"
            or packet["model"]["primary"] != "full_rscc"
            or packet["guards"]["source_selected_flavor_law"] is not False
            or packet["guards"]["formula_discovery_target_informed"] is not True):
        raise ValueError("frozen conditional interpretation changed")
    with localcontext() as ctx:
        ctx.prec = 70
        full = packet["full_rscc"]
        lower = packet["lower_order_control"]
        c = {name: number(value) for name, value in full["coefficients"].items()}
        linear = [number(value) for value in full["basis"]["linear"]]
        quadratic = [number(value) for value in full["basis"]["quadratic"]]
        masses = {name: number(full["native_mass_chart_GeV"][name]) for name in "uds"}
        if any(value <= 0 for value in masses.values()):
            raise ValueError("positive light-quark coordinates required")
        ray_up, ray_down = c["b_u_ray"] / c["a_u"], c["b_d_ray"] / c["a_d"]

        def ratios(coefficients: dict[str, Decimal]) -> tuple[Decimal, Decimal, Decimal]:
            bu = ray_up * coefficients["a_u"] + c["q_u"]
            bd = ray_down * coefficients["a_d"] + c["q_d"]
            eta = (coefficients["mean_u"] - coefficients["mean_d"]) / 2
            u_over_d = (2 * c["B"] * eta
                        + (coefficients["a_u"] - coefficients["a_d"]) * linear[0]
                        + (bu - bd) * quadratic[0]).exp()
            s_over_d = (coefficients["a_d"] * (linear[1] - linear[0])
                        + bd * (quadratic[1] - quadratic[0])).exp()
            return u_over_d, s_over_d, 2 * s_over_d / (1 + u_over_d)

        u_over_d, s_over_d, ratio = ratios(c)
        if any(abs(actual - expected) > D("1e-33") for actual, expected in (
            (u_over_d, masses["u"] / masses["d"]),
            (s_over_d, masses["s"] / masses["d"]),
            (ratio, number(full["light_ratios"]["s_over_mean_ud"])),
        )):
            raise ValueError("frozen exponents and masses disagree")
        w2 = c["w"] ** 2
        corrections = {
            "mean_u_w2_over_29": ("mean_u", w2 / 29),
            "mean_d_minus_w2_over_432": ("mean_d", -w2 / 432),
            "a_u_w2_over_22": ("a_u", w2 / 22),
            "a_d_w2_over_420": ("a_d", w2 / 420),
            "common_delta_g": ("delta_g", c["delta_g"]),
        }
        ablations = {}
        joint = c.copy()
        for label, (name, amount) in corrections.items():
            changed = c.copy()
            changed[name] -= amount
            joint[name] -= amount
            result = ratios(changed)[2]
            ablations[label] = {"ratio_without_term": fmt(result),
                                "change_from_full": fmt(result - ratio)}
        control = ratios(joint)[2]
        if abs(control - number(lower["light_ratios"]["s_over_mean_ud"])) > D("1e-33"):
            raise ValueError("joint ablation fails frozen lower-order control")
        if ratios({**c, "delta_g": D(0)})[2] != ratio:
            raise ValueError("common normalization failed to cancel")
        return {
            "full": {"u_over_d": fmt(u_over_d), "s_over_d": fmt(s_over_d),
                     "s_over_mean_ud": fmt(ratio),
                     "log_s_over_d_linear": fmt(c["a_d"] * (linear[1] - linear[0])),
                     "log_s_over_d_quadratic": fmt(c["b_d"] * (quadratic[1] - quadratic[0]))},
            "single_term_ablations": ablations,
            "joint_lower_order_control": {"s_over_mean_ud": fmt(control),
                                          "change_from_full": fmt(control - ratio)},
            "exact_cancellations": ["common_mass_normalization", "common_delta_g"],
        }


def reference_diagnostic(native: dict, reference: dict) -> dict:
    if (reference["schema"] != "oph.pdg_2026_conditional_quark_comparison.v1"
            or reference["source"]["url"]
            != "https://pdg.lbl.gov/2026/tables/rpp2026-sum-quarks.pdf"
            or [row["quark"] for row in reference["masses"][:3]] != list("uds")):
        raise ValueError("unexpected PDG reference")
    ratios = {row["quantity"]: row for row in reference["light_quantities"]}
    if set(ratios) != {"mu_over_md", "ms_over_mud", "mud_gev"}:
        raise ValueError("all light ratio references required")
    if (ratios["ms_over_mud"]["quoted_confidence"] != "90%"
            or ratios["ms_over_mud"]["definition"] != "m_ud=(m_u+m_d)/2"
            or any(row["scheme"] != "MSbar" or row["scale"] != "2 GeV"
                   for row in reference["masses"][:3])):
        raise ValueError("quoted ratio confidence, definition or light-quark chart changed")
    with localcontext() as ctx:
        ctx.prec = 70
        direct = ratios["ms_over_mud"]
        lower = number(direct["central"]) - number(direct["minus"])
        upper = number(direct["central"]) + number(direct["plus"])
        full = number(native["full"]["s_over_mean_ud"])
        control = number(native["joint_lower_order_control"]["s_over_mean_ud"])
        u, d, s = [number(row["central"]) for row in reference["masses"][:3]]
        listed_mud = number(ratios["mud_gev"]["central"])
        return {
            "pdg_direct_quoted_90_percent": {"lower": fmt(lower), "central": fmt(number(direct["central"])),
                                             "upper": fmt(upper)},
            "full_below_lower_by": fmt(lower - full),
            "lower_order_below_lower_by": fmt(lower - control),
            "minimum_fractional_increase_in_s_at_fixed_u_d": fmt(lower / full - 1),
            "minimum_fractional_decrease_in_mean_ud_at_fixed_s": fmt(1 - full / lower),
            "displayed_marginal_central_arithmetic_only": {
                "from_u_d_s": fmt(2 * s / (u + d)),
                "from_s_and_separately_listed_mean_ud": fmt(s / listed_mud)},
            "boundary": "PDG mass and ratio rows are separate marginals, not one joint fit; neither ablation nor the required shift supplies a physical matching law or a theory-error bound.",
        }


def build() -> dict:
    packet = json.loads(PACKET.read_text())
    reference = json.loads(REFERENCE.read_text())
    spec = importlib.util.spec_from_file_location("conditional_quark_ratio_source_verifier", VERIFIER)
    if spec is None or spec.loader is None:
        raise ValueError("conditional source verifier unavailable")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    verifier.verify_payload(packet)
    native = native_diagnostic(packet)
    return {
        "schema": "oph.conditional_quark_ratio_diagnostic.v1",
        "status": "retrospective_conditional_sensitivity_not_model_selection",
        "native": native,
        "comparison": reference_diagnostic(native, reference),
        "source_pins": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (PACKET, REFERENCE, VERIFIER)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if OUTPUT.read_text() != rendered:
            raise SystemExit("conditional ratio diagnostic drift")
        print("conditional ratio diagnostic verified")
    else:
        OUTPUT.write_text(rendered)
        print(OUTPUT)


if __name__ == "__main__":
    main()

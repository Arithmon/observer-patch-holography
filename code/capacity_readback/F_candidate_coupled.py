#!/usr/bin/env python3
"""Coupled candidate for the capacity readback map F (G2-GAP-1 coupling theorem, 2026-07-15).

Composition (derivation and step ledger: G2_GAP_1_COUPLING_THEOREM.md):

    seed:       pi, from the D6 radius identity N/pi = (r_CRC/ell_star)^2
                (R_N_global_repair_tick_certificate.json, proved_by_certificate)
    form:       port-load inversion F(N) = pi * exp(X_read(N))  [CP-2]
    re-emission: X_read(N) = (1 - lambda) * log(N/pi) + lambda * x_EW(P),
                lambda = 1/2 recorded free in (0,1)               [CP-3]
    coupling:   x_EW(P) = 6*pi / (P * alpha_U(P)), forced by the balance
                condition Gamma_EW = t_tr, equivalently
                Pi_EW = beta_EW * P                               [CP-1]

Fixed point (conditional on CP-1..CP-3): N_CRC = pi * exp(6*pi/(P*alpha_U(P))),
the electroweak bridge capacity of CL-3, independent of lambda.

Blindness status, recorded before any comparison: the evaluation cone of this
candidate contains the bridge exponent 6*pi/(P*alpha_U) by construction. That
expression is named in F_READBACK_SPEC.md Section 3 as barred from a blind
CL-7 candidate. This candidate is therefore conditional (theorem-coupled), it
cannot serve as a blind test of CL-3, and its A7 block below is an informative
comparison, never a landing verdict. Moves no ledger row.

Inputs: certified forward-closure P enclosure and certified alpha_U(P_fwd)
enclosure, both from
code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json.
No measured Lambda, no SL-4 estimate as input (it appears once, inside the
clearly marked informative-comparison block).
"""

from __future__ import annotations

import json
from pathlib import Path

from mpmath import iv, mp, mpf

from toy_readback import _endpoints, _interval_json

from F_candidate_capL import (
    P_HI,
    P_LO,
    contraction_certificate_centered,
    fixed_point_enclosure_centered,
)

ARTIFACT_NAME = "oph_capacity_readback_candidate_coupled_2026-07-15"
PRECISION = 60

# Certified alpha_U(P_fwd) enclosure, interval_diagnostics.alpha_u of the
# P interval contraction certificate (mpmath.iv, outward rounding, iv_dps 60).
ALPHA_U_LO = "0.041124247441816685140889933889659717292128290416193554171910567414"
ALPHA_U_HI = "0.041124247441816685140889933889659717292128290616193554171910678971"

ENCLOSURE_HALF_WIDTH = "1e-25"
LAMBDA = "0.5"


def build() -> dict:
    iv.dps = PRECISION
    mp.dps = PRECISION
    P = iv.mpf([P_LO, P_HI])
    alpha_u = iv.mpf([ALPHA_U_LO, ALPHA_U_HI])
    lam = iv.mpf(LAMBDA)
    one = iv.mpf(1)

    # D10-side load forced by CP-1: x_EW = 6*pi/(P*alpha_U).
    x_ew = 6 * iv.pi / (P * alpha_u)
    x_lo, x_hi = _endpoints(x_ew)

    # --- load-coordinate certification: y = log(N/pi), C(y) = (1-lam)*y + lam*x_EW
    def C(y):
        return (one - lam) * y + lam * x_ew

    def Cp(y):
        return (one - lam) + 0 * y

    interval_y = iv.mpf([mp.nstr(x_lo - 1, PRECISION), mp.nstr(x_hi + 1, PRECISION)])
    cert_y = contraction_certificate_centered(C, Cp, interval_y)
    enclosure_y = fixed_point_enclosure_centered(
        C, Cp, interval_y, ENCLOSURE_HALF_WIDTH, lipschitz_L=mpf(LAMBDA)
    )
    y_box = iv.mpf([enclosure_y["enclosure"]["lo"], enclosure_y["enclosure"]["hi"]])

    # exact-solution cross-check: the affine fixed point is x_EW itself
    exact_check = {
        "algebraic_fixed_point": "y = x_EW(P), independent of lambda",
        "x_ew_enclosure": _interval_json(x_ew),
        "x_ew_inside_certified_box": bool(
            mpf(enclosure_y["enclosure"]["lo"]) <= x_lo
            and x_hi <= mpf(enclosure_y["enclosure"]["hi"])
        ),
    }

    # --- N-coordinate certification around the fixed point
    def F(n):
        return iv.pi * iv.exp((one - lam) * iv.log(n / iv.pi) + lam * x_ew)

    def Fp(n):
        return (one - lam) * F(n) / n

    n_ew = iv.pi * iv.exp(x_ew)
    n_lo, n_hi = _endpoints(n_ew)
    interval_n = iv.mpf([mp.nstr(n_lo * mpf("0.9"), PRECISION), mp.nstr(n_hi * mpf("1.1"), PRECISION)])
    cert_n = contraction_certificate_centered(F, Fp, interval_n)

    # certified N enclosure from the load-coordinate box, outward exponentiation
    n_box = iv.pi * iv.exp(y_box)
    n_box_lo, n_box_hi = _endpoints(n_box)

    certificate_core = {
        "coupling": {
            "cp1_balance_condition": "Pi_EW(P,N) = beta_EW*P with Pi_EW = 24*pi/(alpha_U*log(N/pi)); equivalently Gamma_EW = t_tr; equivalently per-port load log(N/pi)/12 = pi/(2*P*alpha_U)",
            "forced_load": "x_EW(P) = 6*pi/(P*alpha_U(P))",
            "x_ew_enclosure": _interval_json(x_ew),
        },
        "load_coordinate": {
            "map": "C(y) = (1-lambda)*y + lambda*x_EW, y = log(N/pi)",
            "lambda": LAMBDA,
            "contraction_certificate": cert_y,
            "fixed_point": enclosure_y,
            "exact_solution_check": exact_check,
        },
        "capacity_coordinate": {
            "map": "F(N) = pi * exp((1-lambda)*log(N/pi) + lambda*x_EW)",
            "contraction_certificate": cert_n,
            "fixed_point_enclosure": {
                **_interval_json(n_box),
                "relative_width": mp.nstr((n_box_hi - n_box_lo) / n_box_lo, 8),
            },
        },
    }

    # spec property record P1-P5
    properties = {
        "P1_well_definedness": "the coupled map is total and C^1 on the certificate interval; well-definedness of the physical nf/Obs/Cap_read chain at finite cutoff is carried by CP-2/CP-3, unconstructed",
        "P2_monotone": "F' = (1-lambda)*F(N)/N > 0 on the interval; certified enclosure above",
        "P3_bounds": "bracketing pair supplied by the certified self-map interval; L = 1/2",
        "P4_count_density": "open obligation: the count representation of the coupled membership clause is unconstructed; recorded with CP-2/CP-3, no coherence claim",
        "P5_nontrivial": "F' = 1/2 on the load coordinate; the map is neither the identity nor constant for lambda in (0,1); the lambda = 1 constant reading is the excluded CAP-B branch",
    }

    # A6 blindness record, written before the comparison block below
    blindness = {
        "inputs": [
            "certified forward-closure P enclosure (source-side computed fixed point)",
            "certified alpha_U(P_fwd) enclosure (same interval certificate)",
            "declared structure integers 12/24/6/4 and pi",
            "lambda = 1/2 recorded free averaging weight",
        ],
        "reads_measured_lambda": False,
        "reads_sl4_estimate_as_input": False,
        "cone_contains_cl3_bridge_expression": True,
        "v08_verdict": "fails the spec Section 3 blindness bar as written: the bridge exponent 6*pi/(P*alpha_U) is in the cone by construction (it is the theorem content). This candidate is conditional on CP-1..CP-3 and cannot serve as a blind CL-7 landing test; if the premises discharge, CL-7 closes into CL-3 and the live test is the bridge value against the Lambda-located capacity.",
        "dependency_cone": [
            "mpmath",
            "toy_readback interval machinery",
            "F_candidate_capL centered certificate helpers",
            "P_derivation interval contraction certificate (P and alpha_U enclosures)",
        ],
    }

    # ------------------------------------------------------------------
    # Informative comparison (A7-shaped, never a landing verdict).
    # The SL-4 reference value appears here for the first and only time.
    # ------------------------------------------------------------------
    sl4_reference = mpf("3.31e122")
    n_mid = (n_box_lo + n_box_hi) / 2
    comparison = {
        "status": "informative_only_not_a_landing_verdict",
        "reason": "A6 fails by construction; the spec landing protocol does not apply to a theorem-coupled candidate",
        "sl4_lambda_located_capacity": "3.31e122",
        "fixed_point_over_sl4": mp.nstr(n_mid / sl4_reference, 10),
        "log10_offset": mp.nstr(mp.log10(n_mid / sl4_reference), 6),
        "relative_offset": mp.nstr(n_mid / sl4_reference - 1, 6),
        "note": "the offset restates the CL-3 residual (bridge value 6.6 percent above the Lambda-located capacity in N units); it is the single live test if CP-1..CP-3 discharge",
    }

    return {
        "artifact": ARTIFACT_NAME,
        "specification": "F_READBACK_SPEC.md",
        "derivation": "G2_GAP_1_COUPLING_THEOREM.md",
        "status": "conditional_on_CP-1_CP-2_CP-3",
        "family": "coupled port-inversion readback: F(N) = pi * exp((1-lambda)*log(N/pi) + lambda*6*pi/(P*alpha_U))",
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": PRECISION,
            "rounding": "mpmath_interval_outward",
            "promotion_backend_required": "arb_or_mpfi_directed_outward",
        },
        "inputs": {
            "P_certified_enclosure": {"lo": P_LO, "hi": P_HI},
            "alpha_U_certified_enclosure": {"lo": ALPHA_U_LO, "hi": ALPHA_U_HI},
            "source": "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json",
            "structure_integers": [12, 24, 6, 4],
        },
        "certificate": certificate_core,
        "spec_properties": properties,
        "blindness": blindness,
        "informative_comparison": comparison,
        "moves_cl7": False,
        "cl7_status": "open_reduced_to_CP-1_CP-2_CP-3",
    }


def main() -> int:
    certificate = build()
    out = Path(__file__).resolve().parent / "runtime" / "F_candidate_coupled_certificates.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "status": certificate["status"],
        "fixed_point_N": certificate["certificate"]["capacity_coordinate"]["fixed_point_enclosure"],
        "load_banach_pass": certificate["certificate"]["load_coordinate"]["contraction_certificate"]["banach_pass"],
        "capacity_banach_pass": certificate["certificate"]["capacity_coordinate"]["contraction_certificate"]["banach_pass"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

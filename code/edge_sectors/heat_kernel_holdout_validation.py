#!/usr/bin/env python3
"""Held-out diagnostics of the edge-sector heat-kernel law on finite gauge groups.

This script accompanies the "Numerical diagnostics of the heat-kernel law"
section of the synthesis paper and implements the finite comparison protocol
requested by paper-audit issue #540:

  * it separates symmetry-forced checks from overconstrained tests;
  * every substantive test fits the diffusion parameter t on a declared
    subset of spectral sectors and prints held-out residuals for the
    remaining sectors;
  * groups with degenerate nontrivial spectra (Z2, Z3) are labelled as
    implementation checks, because for them the multi-sector agreement is
    forced by charge conjugation plus eigenvalue degeneracy before any
    heat-kernel ansatz is imposed.

Models
------
Z_n  : 2x2 periodic lattice gauge theory (8 links) with Z_n link spaces and

           H = -K sum_p Re(B_p) - h sum_l Re(X_l) - Gamma sum_v Re(A_v),

       exactly the Hamiltonian displayed in the paper (K = 1, Gamma = 5).
       Its unique h>0 ground state is computed in the n^3-dimensional
       zero-divergence/zero-winding electric sector (see README.md).
       Region A consists of links whose tail has x = 0; the electric-center
       edge charge at a boundary vertex v is the restricted star
       Q_v = prod_{l in star(v) cap A} X_l^{+/-1}.

S3   : the exact single-plaquette reduction described in the paper.  With
       Gauss's law imposed at every vertex the physical space is the
       3-dimensional space of class functions of the plaquette holonomy,
       spanned by the normalized characters {chi_triv, chi_sign, chi_std}.
       The magnetic term is multiplication by Re chi_std(g); the electric
       term is diagonal in the character basis with the Cayley-graph
       Laplacian eigenvalues (transposition generating set)
       lambda_triv = 0, lambda_sign = 6, lambda_std = 3.

Extraction and held-out protocol
--------------------------------
The heat-kernel ansatz is p_R proportional to d_R exp(-t lambda_R).  For each
nontrivial sector, t_R = ln((p_0/d_0)/(p_R/d_R)) / lambda_R.  The fit sector
(declared below, always the lowest nonzero eigenvalue) determines t; every
other nontrivial sector's weight is then a parameter-free held-out
prediction, and the printed residual is

    residual(R) = ln(p_R_measured / p_R_predicted) / abs(ln(p_R_predicted/p_0))

(a relative log-scale error), together with the eigenvalue-ratio diagnostic
log(p_R/p_0 d_R) / log(p_fit/p_0 d_fit) versus lambda_R / lambda_fit.

The Hamiltonians and their electric normalization are fixed inputs. Changing
the electric term changes the ground state. Only a common rescaling of the
extraction eigenvalues can be absorbed into the fitted diffusion parameter.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import math
from numbers import Integral, Real

import mpmath
import numpy as np

if __package__ in (None, ""):
    from abelian_ground_state import zn_edge_distribution
    from nonabelian_ground_state import s3_diagnostics, s3_edge_distribution
else:
    from .abelian_ground_state import zn_edge_distribution
    from .nonabelian_ground_state import s3_diagnostics, s3_edge_distribution


# ----------------------------------------------------------------------
# Held-out fitting.
# ----------------------------------------------------------------------

def fit_t(p0, p_fit, d_fit, lam_fit):
    """Fit supplied weights, requiring resolved final binary64 conversion."""
    values = [_positive(value) for value in (p0, p_fit, d_fit, lam_fit)]
    base, weight, dimension, eigenvalue = values
    ratio = Fraction(base)*Fraction(dimension)/Fraction(weight)
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    if Fraction(1, 2) <= ratio <= 2:
        delta = ratio-1
        logarithm = ctx.log1p(ctx.mpf(delta.numerator)/delta.denominator)
    else:
        logarithm = ctx.log(ctx.mpf(ratio.numerator)/ratio.denominator)
    return _resolved_float(logarithm/ctx.mpf(eigenvalue), "fitted time", zero_allowed=ratio == 1)


def _resolved_float(value, name, *, zero_allowed=False):
    result = float(value)
    if (not math.isfinite(result) or (result == 0 and not zero_allowed)
            or (value and abs((value.context.mpf(result)-value)/value) > value.context.mpf("1e-12"))):
        raise ValueError(f"{name} exceeds the resolved reporting range")
    return result


def _binary_real(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError("finite real fit parameters required")
    try:
        converted = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError("fit parameters exceed binary64 range") from exc
    if not math.isfinite(converted):
        raise ValueError("finite real fit parameters required")
    if isinstance(value, Integral):
        original = Fraction(int(value))
    elif isinstance(value, Fraction):
        original = Fraction(int(value.numerator), int(value.denominator))
    else:
        numerator, denominator = value.as_integer_ratio()
        original = Fraction(int(numerator), int(denominator))
    if Fraction(converted) != original:
        raise ValueError("fit parameters must be exactly represented binary64 values")
    return converted


def _positive(value):
    converted = _binary_real(value)
    if converted <= 0:
        raise ValueError("positive real fit parameters required")
    return converted


def predict(p0, d, lam, t):
    """Predict without intermediate underflow or unresolved output rounding."""
    p0, d, lam = map(_positive, (p0, d, lam))
    t = _binary_real(t)
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    logarithm = ctx.log(ctx.mpf(p0))+ctx.log(ctx.mpf(d))-ctx.mpf(t)*ctx.mpf(lam)
    return _resolved_float(ctx.exp(logarithm), "positive prediction")


def report_zn(n, h_values):
    lam = [4.0 * math.sin(math.pi * q / n) ** 2 for q in range(n)]
    distinct = len({round(l, 12) for l in lam[1:]})
    if distinct < 2:
        kind = "SYMMETRY-FORCED IMPLEMENTATION CHECK"
    else:
        kind = "OVERCONSTRAINED HELD-OUT TEST"
    print(f"\n=== Z_{n} ({kind}) ===")
    print(f"eigenvalues: {[round(l, 6) for l in lam]}")
    if n == 3:
        print(
            "note: lambda_1 = lambda_2 and charge conjugation force "
            "p_1 = p_2 and t_{q=1} = t_{q=2} for ANY conjugation-invariant "
            "distribution; agreement below tests the implementation only."
        )
    header = "h      " + "  ".join(f"p_{q:<8d}" for q in range(n)) + "  t(fit q=1)  held-out residuals"
    print(header)
    for h in h_values:
        p = zn_edge_distribution(n, h)
        t = fit_t(p[0], p[1], 1.0, lam[1])
        if abs(t*lam[1]) <= 1e-8:
            raise ValueError("Z_n log-gap is unresolved for a normalized sector comparison")
        residuals = []
        for q in range(2, (n // 2) + 1):
            pred = predict(p[0], 1.0, lam[q], t)
            log_gap = fit_t(p[0], p[q], 1., 1.)
            res = (t*lam[q]-log_gap)/abs(t*lam[q])
            ratio = log_gap/(t*lam[1])
            residuals.append(
                f"q={q}: pred {pred:.3e} meas {p[q]:.3e} "
                f"res {res:+.2%} ratio {ratio:.4f} (target {lam[q]/lam[1]:.4f})"
            )
        row = f"{h:<6.2f} " + "  ".join(f"{p[q]:.4e}" for q in range(n)) + f"  {t:<10.4f}"
        print(row)
        for r in residuals:
            print(f"        {r}")


def report_s3(h_values):
    print("\n=== S_3 (OVERCONSTRAINED HELD-OUT TEST, nonabelian) ===")
    print("eigenvalues: triv 0, sign 6, std 3 (distinct nonzero pair)")
    print("fit sector: std; held-out sector: sign; target log-ratio = 2")
    print("h      p_triv      p_sign      p_std       t(std)   pred p_sign  res      log-ratio excess")
    for h in h_values:
        result = s3_diagnostics(h)
        p = result["probabilities"]
        print(
            f"{h:<6.2f} {p['triv']:.4e}  {p['sign']:.4e}  {p['std']:.4e}  "
            f"{result['fit_time']:<8.4f} {result['predicted_sign']:.3e}    "
            f"{result['normalized_log_residual']:+.3e}  "
            f"{result['log_ratio_excess']:+.6e}"
        )
        print(f"        log(measured/predicted)={result['log_discrepancy']:+.6e}; "
              f"positive diffusion fit: {result['diffusion_fit']}")



def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--groups",
        nargs="+",
        default=["Z2", "Z3", "Z5", "S3"],
        choices=["Z2", "Z3", "Z5", "S3"],
    )
    parser.add_argument(
        "--h", nargs="+", type=float, default=None, help="electric couplings to scan"
    )
    args = parser.parse_args()

    for group in args.groups:
        if group == "S3":
            report_s3(args.h or [0.5, 1.0, 2.0, 5.0, 12.0, 100.0])
        else:
            n = int(group[1:])
            default_h = [0.05, 0.1, 0.2, 0.5] if n == 5 else [0.2, 0.5, 1.0, 2.0]
            report_zn(n, args.h or default_h)

    print(
        "\nSummary: Z2/Z3 are symmetry-forced implementation checks. "
        "Z5 fits t at q=1 and compares q=2; its samples supply no limit theorem. "
        "For the specified S3 one-plaquette model, the analytic sign-sector "
        "measured/predicted ratio is strictly below one at every finite h>=0 "
        "and increases to one as h grows. The fitted diffusion time is positive "
        "only above (1+sqrt(3)-sqrt(2))/24. These are results for the supplied "
        "finite Hamiltonians, with no continuous-group or physical-model transfer."
    )


if __name__ == "__main__":
    main()

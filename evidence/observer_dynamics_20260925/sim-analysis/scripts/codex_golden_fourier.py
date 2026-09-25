"""Exact finite-product formula; bounded numerical readings of prescribed modes.

The source population is a Cartesian product of xi_b = frac(b*phi), b=0..q-1.
For integer cube Fourier modes n, wrapping xi changes no phase. Hence
F_q(n) = product_j [sum_b exp(2*pi*i*n_j*b*phi) / q]. The one-dimensional
sum is geometric. This is a finite identity, not a cosmological spectrum.

Generate::
  python3 sim-analysis/scripts/codex_golden_fourier.py \
    sim-analysis/data/codex_audit_20260925/golden_fourier.json

Use --verify to check the deterministic receipt and its producer hash.
"""

from __future__ import annotations

import argparse
import cmath
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
import hashlib
import itertools
import json
import math
from pathlib import Path


LEVELS = (55, 89, 144)
# Declared finite list, including equal-|n| direction pairs. No optimized bands.
MODES = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1),
         (3, 0, 0), (2, 2, 1), (5, 0, 0), (4, 3, 0),
         (7, 0, 0), (6, 3, 2))
DIRECTION_PAIRS = (((3, 0, 0), (2, 2, 1)),
                   ((5, 0, 0), (4, 3, 0)),
                   ((7, 0, 0), (6, 3, 2)))


def _phi():
    return (Decimal(1) + Decimal(5).sqrt()) / 2


def _sin_pi(value):
    """Reduce in Decimal before passing a small argument to libm."""
    integer = value.to_integral_value(rounding=ROUND_HALF_EVEN)
    small = float(value - integer)
    return (-1 if int(integer) % 2 else 1) * math.sin(math.pi * small)


def axis_amplitude(q, n):
    if q <= 0 or not isinstance(q, int) or not isinstance(n, int):
        raise ValueError("q must be a positive integer and n an integer")
    if n == 0:
        return complex(1)
    with localcontext() as ctx:
        ctx.prec = 80
        t = Decimal(n) * _phi()
        ratio = _sin_pi(q * t) / (q * _sin_pi(t))
        phase = float(((q - 1) * t / 2) % 1)
    return ratio * cmath.exp(2j * math.pi * phase)


def amplitude(q, mode):
    if len(mode) != 3:
        raise ValueError("a spatial mode has three coordinates")
    return math.prod(axis_amplitude(q, n) for n in mode)


def explicit_cloud_amplitude(q, mode):
    """Independent q^3 evaluation, deliberately bounded to tiny validation clouds."""
    if q < 1 or q > 13 or len(mode) != 3:
        raise ValueError("explicit validation is limited to 1 <= q <= 13")
    with localcontext() as ctx:
        ctx.prec = 80
        phi = _phi()
        positions = [float((b * phi) % 1) for b in range(q)]
    phases = [cmath.exp(2j * math.pi * math.fsum(n * x for n, x in zip(mode, point)))
              for point in itertools.product(positions, repeat=3)]
    return complex(math.fsum(z.real for z in phases), math.fsum(z.imag for z in phases)) / q**3


def fibonacci_index(q):
    a, b, index = 0, 1, 0
    while a < q:
        a, b, index = b, a + b, index + 1
    if a != q:
        raise ValueError("q is not Fibonacci")
    return index


def asymptotic_constant(mode):
    with localcontext() as ctx:
        ctx.prec = 80
        phi = _phi()
        return math.prod((math.pi * n) ** 2 / (5 * _sin_pi(n * phi) ** 2)
                         for n in mode if n)


def mode_row(q, mode):
    value = amplitude(q, mode)
    power = abs(value) ** 2
    h = sum(n != 0 for n in mode)
    constant = asymptotic_constant(mode)
    return {
        "mode": list(mode), "norm_squared": sum(n * n for n in mode),
        "nonzero_coordinates": h,
        "amplitude_real": value.real, "amplitude_imag": value.imag,
        "normalized_density_power": power,
        "structure_factor_N_times_power": q**3 * power,
        "fixed_count_uniform_expected_power": 1 / q**3 if h else 1,
        "ratio_to_fixed_count_uniform_expectation": q**3 * power if h else 1,
        "fixed_mode_fibonacci_asymptotic": {
            "predicted_power_exponent_in_q": -4 * h,
            "constant": constant,
            "observed_power_times_q_to_4h_over_constant": power * q**(4 * h) / constant,
            "fitted": False,
        },
    }


def make_receipt():
    levels = []
    for q in LEVELS:
        rows = [mode_row(q, mode) for mode in MODES]
        by_mode = {tuple(r["mode"]): r for r in rows}
        pairs = []
        for a, b in DIRECTION_PAIRS:
            assert sum(n * n for n in a) == sum(n * n for n in b)
            pairs.append({
                "mode_a": list(a), "mode_b": list(b),
                "shared_norm_squared": sum(n * n for n in a),
                "power_a_over_power_b": by_mode[a]["normalized_density_power"] / by_mode[b]["normalized_density_power"],
            })
        levels.append({"q": q, "population": q**3, "fibonacci_index": fibonacci_index(q),
                       "modes": rows, "equal_norm_direction_pairs": pairs})
    checks = []
    for q in (2, 3, 5, 8):
        for mode in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (2, -1, 0), (3, 2, 1)):
            explicit = explicit_cloud_amplitude(q, mode)
            factored = amplitude(q, mode)
            checks.append({"q": q, "mode": list(mode),
                           "complex_amplitude_absolute_error": abs(explicit - factored),
                           "power_absolute_error": abs(abs(explicit)**2 - abs(factored)**2)})
    return {
        "schema": "oph.analysis.golden-cartesian-fourier.v1",
        "definition": {
            "positions": "L*(frac(b1*phi), frac(b2*phi), frac(b3*phi)); each bj=0,...,q-1",
            "mode": "k=(2*pi/L)*n with n in Z^3",
            "amplitude": "F_q(n)=q^-3 sum_b exp(2*pi*i*n.dot(frac(b*phi)))",
            "power": "abs(F_q(n))^2; no radial average or bin fitting",
            "structure_factor": "S_q(n)=q^3*abs(F_q(n))^2; n=0 is conserved total, not a fluctuation",
            "uniform_control": "for N independent uniform points in the same cube, E|F(n)|^2=1/N at each nonzero integer mode",
        },
        "theorems": {
            "factorization": "F_q(n)=product_j a_q(n_j), a_q(0)=1, a_q(n)=exp(pi*i*n*phi*(q-1))*sin(pi*n*q*phi)/(q*sin(pi*n*phi)) for n!=0",
            "fibonacci_identity": "phi*F_m-F_(m+1)=(-1)^(m+1)*phi^-m",
            "fixed_mode_asymptotic": "for h nonzero coordinates and fixed n, q^(4h)*|F_q(n)|^2 -> product_(n_j!=0) (pi*n_j)^2/[5*sin^2(pi*n_j*phi)] along q=F_m",
            "status": "elementary algebraic derivations documented in the audit note; not Lean checked",
        },
        "levels": levels,
        "explicit_small_cloud_validation": checks,
        "maximum_complex_absolute_error": max(r["complex_amplitude_absolute_error"] for r in checks),
        "numerics": "Decimal80 evaluates phi and reduces phases; double-precision sin/exp evaluate the final bounded expressions. These are numerical readings, not certified interval enclosures.",
        "nonclaims": [
            "The fixed Cartesian source geometry is an input, not a structure generated by repair.",
            "Integer cube Fourier modes avoid the uniform-window mean, but do not remove all finite-window interpretation issues.",
            "Fixed-box, fixed-mode Fibonacci refinement is not an infinite-volume hyperuniformity limit.",
            "Fixed modes decay absolutely even when ratios between directions grow; a direction ratio alone does not establish finite continuum anisotropy.",
            "No primordial scalar tilt, CMB temperature spectrum, radiative transfer, or physical wavelength is fitted or inferred.",
            "The formula describes the stated mathematical population; raw geometry files were not downloaded or byte-verified for this calculation.",
        ],
        "provenance": {
            "producer_name": Path(__file__).name,
            "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "cloud_downloaded": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    receipt = make_receipt()
    if args.verify:
        if json.loads(args.output.read_text()) != receipt:
            raise SystemExit("Fourier receipt differs from current producer")
        print("Verified formula readings, explicit small-cloud checks and producer hash")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        print(args.output)


if __name__ == "__main__":
    main()

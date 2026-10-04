"""Independent replay of retained witnesses using a direct edge-orbit formula.

No producer import, empirical likelihood or promotion of a frozen protocol.
"""

import itertools
import json
from pathlib import Path
import re

import mpmath as mp

from .bounds import verify_budget

DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE]([+-]?[0-9]+))?")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def decimal(value):
    require(isinstance(value, str) and len(value) < 100, "invalid decimal field")
    match = DECIMAL.fullmatch(value)
    require(match is not None, "invalid decimal field")
    # Bound the representation before invoking arbitrary-precision arithmetic.
    # This comfortably includes the diagnostics from all allowed precisions.
    require(abs(int(match.group(1) or "0")) <= 1000, "decimal exponent outside report range")
    try:
        result = mp.mpf(value)
    except (ValueError, TypeError):
        raise ValueError("invalid decimal field") from None
    require(mp.isfinite(result), "nonfinite decimal")
    return result


def direct_orbit():
    phi = (1+mp.sqrt(5))/2
    edges = []
    for i in range(3):
        for sign in (-1, 1):
            edges.append(tuple(mp.mpf(sign if j == i else 0) for j in range(3)))
    for signs in itertools.product((-1, 1), repeat=3):
        v = (signs[0]/mp.mpf(2), signs[1]/(2*phi), signs[2]*phi/2)
        edges.extend((v, v[1:]+v[:1], v[2:]+v[:2]))
    return edges


def check(packet):
    verify_budget()
    require(type(packet) is dict and set(packet) == {
        "scope", "global_threshold_absolute_error_eV", "witnesses"}, "report schema")
    require(packet["scope"] == "fixed-hypothesis kinematic feasibility; no empirical verdict", "scope drift")
    require(packet["global_threshold_absolute_error_eV"] == "1e-13", "error claim drift")
    rows = packet["witnesses"]
    require(type(rows) is list and len(rows) == 30, "witness census")
    with mp.workdps(90):
        m = mp.mpf("510998.95069")
        a = mp.sqrt(mp.mpf("1.6309682094039593"))*mp.mpf("1.616255e-35")/mp.mpf("1.973269804e-7")
        d = a*a/20
        orbit = direct_orbit()
        require(mp.mpf("5e-58") < d < mp.mpf("6e-58"), "coefficient outside proof box")
        seen = set()
        for row in rows:
            require(type(row) is dict and set(row) == {
                "hard_momentum_eV", "variant", "direction", "leading_soft_eV",
                "collinear_soft_energy_eV", "outgoing_small_fraction", "energy_residual_eV"}, "row schema")
            k = decimal(row["hard_momentum_eV"])
            require(k in tuple(mp.mpf(s) for s in ("1e17", "3e17", "1e18", "1e19", "1e20")), "momentum census")
            variant, direction = row["variant"], row["direction"]
            require(variant in ("photon_only", "photon_electron_positron"), "variant census")
            require(type(direction) is list and all(type(v) is int for v in direction), "direction schema")
            require(tuple(direction) in ((1,0,0), (1,1,1), (1,2,3)), "direction census")
            key = (k, variant, tuple(direction))
            require(key not in seen, "duplicate witness")
            seen.add(key)
            norm = mp.sqrt(sum(v*v for v in direction))
            dots = [sum(v[i]*direction[i] for i in range(3))/norm for v in orbit]

            def energy(r, mass):
                if mass and variant == "photon_only":
                    return mp.sqrt(r*r+mass*mass)
                # Complex chord identity, independently assembled support.
                lam = mp.fsum(abs(mp.expm1(1j*a*r*z))**2 for z in dots)/(10*a*a)
                return mp.sqrt(mass*mass+lam)

            eps = decimal(row["collinear_soft_energy_eV"])
            share = decimal(row["outgoing_small_fraction"])
            reported_leading = decimal(row["leading_soft_eV"])
            require(0 < eps < 200 and 0 < share <= mp.mpf("0.5"), "witness domain")
            # At s<200, |Omega(s)-s|<d*200^3<5e-51 eV.
            p, q = share*(k-eps), (1-share)*(k-eps)
            residual = energy(p,m)+energy(q,m)-energy(k,0)-energy(eps,0)
            require(abs(residual) < mp.mpf("1e-29"), "energy conservation")
            # Internal solve diagnostic; the serialized witness is checked
            # independently above at its own (35-digit) rounding tolerance.
            require(abs(decimal(row["energy_residual_eV"])) < mp.mpf("1e-35"), "reported residual")
            slope_difference = mp.diff(lambda r: energy(r,m),p)-mp.diff(lambda r: energy(r,m),q)
            require(abs(slope_difference) < mp.mpf("1e-40"), "stationarity")
            curvature = mp.diff(lambda r: energy(r,m),p,2)+mp.diff(lambda r: energy(r,m),q,2)
            require(curvature > 0, "stationary point is not a minimum")
            if variant == "photon_only":
                expected = m*m/k+d*k**3/4
            else:
                transition = (16*m*m/(3*d))**mp.mpf("0.25")
                expected = m*m/k+3*d*k**3/16 if k <= transition else m*mp.sqrt(3*d)*k/2
            require(abs(reported_leading-expected) < abs(expected)*mp.mpf("1e-33"), "leading formula")
            require(abs(eps-expected) < mp.mpf("1e-13"), "analytic threshold enclosure")
    return True


def load(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    with Path(path).open("rb") as source:
        raw = source.read(50001)
    require(len(raw) <= 50000, "oversized report")
    return json.loads(raw, object_pairs_hook=unique)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", default=str(Path(__file__).with_name("report.json")))
    args = parser.parse_args()
    check(load(args.path))
    print("Exact rational global error budget and all 30 independent witnesses pass.")

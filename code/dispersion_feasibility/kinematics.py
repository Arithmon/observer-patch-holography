"""Exact cosine-shell witnesses; no flux, cross section or observational score.

The universal control means photons and e+/e- only. Its massive shell is an
explicit hypothesis, E_m(p)^2 = m^2 + Lambda(p), not an OPH derivation.
All constants below are nominal inputs, not zero-input predictions.
"""

from functools import lru_cache
from fractions import Fraction

import mpmath as mp

P = Fraction("1.6309682094039593")
PLANCK_M = Fraction("1.616255e-35")
HBARC_EV_M = Fraction("1.973269804e-7")
MASS_EV = Fraction("510998.95069")
A2 = P * (PLANCK_M / HBARC_EV_M) ** 2
D = A2 / 20
MOMENTA = ("1e17", "3e17", "1e18", "1e19", "1e20")
DIRECTIONS = ((1, 0, 0), (1, 1, 1), (1, 2, 3))
VARIANTS = ("photon_only", "photon_electron_positron")


def number(value):
    if isinstance(value, Fraction):
        return mp.mpf(value.numerator) / value.denominator
    return mp.mpf(value)


def validate(k, variant, precision):
    if type(precision) is not int or not 60 <= precision <= 120:
        raise ValueError("precision must be an integer from 60 to 120")
    if variant not in VARIANTS:
        raise ValueError("unknown lepton variant")
    if isinstance(k, bool):
        raise ValueError("invalid hard momentum")
    k = number(k)
    if not mp.isfinite(k) or not mp.mpf("1e17") <= k <= mp.mpf("1e20"):
        raise ValueError("hard momentum outside the declared audit range")
    return k


def leading(k, variant):
    """Leading minimum, with x the smaller outgoing momentum fraction."""
    m, d = number(MASS_EV), number(D)
    if variant == "photon_only":
        return m*m/k + d*k**3/4, mp.mpf("0.5")
    if variant != "photon_electron_positron":
        raise ValueError("unknown lepton variant")
    y = min(mp.mpf("0.25"), m / (mp.sqrt(3*d)*k*k))
    # Rationalized smaller root avoids cancellation as y tends to zero.
    x = 2*y / (1 + mp.sqrt(1-4*y))
    return (m*m/y + 3*d*k**4*y)/(4*k), x


@lru_cache(maxsize=4)
def _edge_vectors(precision):
    # Independent geometric construction from the twelve regular vertices.
    with mp.workdps(precision):
        phi = (1 + mp.sqrt(5))/2
        vertices = []
        for s in (-1, 1):
            for t in (-1, 1):
                vertices.extend(((0, s, t*phi), (s, t*phi, 0), (t*phi, 0, s)))
        edges = []
        for i, v in enumerate(vertices):
            for w in vertices[i+1:]:
                diff = tuple(w[j]-v[j] for j in range(3))
                if abs(sum(z*z for z in diff)-4) < mp.mpf(10)**(-precision+10):
                    edges.append(tuple(z/2 for z in diff))
        if len(edges) != 30:
            raise ValueError("invalid edge census")
        return tuple(edges)


def projections(direction, precision):
    if len(direction) != 3 or any(isinstance(x, bool) for x in direction):
        raise ValueError("direction must have three real coordinates")
    values = tuple(number(x) for x in direction)
    if not all(mp.isfinite(x) for x in values):
        raise ValueError("nonfinite direction")
    norm = mp.sqrt(sum(x*x for x in values))
    if not norm:
        raise ValueError("zero direction")
    n = tuple(x/norm for x in values)
    return tuple(sum(v[j]*n[j] for j in range(3)) for v in _edge_vectors(precision))


def shell(r, mass, a, dots):
    """Energy, energy minus r, and radial group velocity, without 1-cos loss."""
    lam = mp.fsum(2*mp.sin(a*r*z/2)**2 for z in dots)/(5*a*a)
    energy = mp.sqrt(mass*mass + lam)
    correction = (mass*mass + lam-r*r)/(energy+r) if energy+r else mp.mpf(0)
    slope = mp.fsum(z*mp.sin(a*r*z) for z in dots)/(10*a*energy) if energy else mp.mpf(1)
    return energy, correction, slope


def witness(k, variant, direction=(1, 0, 0), precision=70):
    """Solve exact collinear energy conservation and stationarity.

    This supplies an independent numerical witness, not a numerical proof of
    global minimization. The direction-uniform analytic bracket is separate.
    """
    with mp.workdps(precision if type(precision) is int and 60 <= precision <= 120 else 70):
        k = validate(k, variant, precision)
        a, m = mp.sqrt(number(A2)), number(MASS_EV)
        dots = projections(direction, precision)
        eps0, x0 = leading(k, variant)
        _, hard_correction, _ = shell(k, 0, a, dots)

        def lepton(r):
            if variant == "photon_only":
                e = mp.sqrt(r*r+m*m)
                return e, m*m/(e+r), r/e
            return shell(r, m, a, dots)

        def residual(s, p):
            q = k-s-p
            return -2*s + lepton(p)[1] + lepton(q)[1] - hard_correction - shell(s, 0, a, dots)[1]

        if x0 == mp.mpf("0.5"):
            s = mp.findroot(lambda s: residual(s, (k-s)/2), (eps0*mp.mpf("0.999"), eps0*mp.mpf("1.001")))
            p = (k-s)/2
        else:
            scale = m*m/(x0*k)**2 + number(D)*k*k

            def equations(u, v):
                s, p = eps0*u, x0*k*v
                q = k-s-p
                return residual(s, p)/eps0, (lepton(p)[2]-lepton(q)[2])/scale

            u, v = mp.findroot(equations, (1, 1), tol=mp.mpf(10)**(-precision+20))
            s, p = eps0*u, x0*k*v
        if not 0 < p <= (k-s)/2 or not 0 < s < 200:
            raise ArithmeticError("witness outside the declared domain")
        energy_residual = residual(s, p)
        if abs(energy_residual) > mp.mpf("1e-35"):
            raise ArithmeticError("unresolved energy equation")
        return {
            "hard_momentum_eV": str(k), "variant": variant,
            "direction": list(direction),
            "leading_soft_eV": mp.nstr(eps0, 35),
            "collinear_soft_energy_eV": mp.nstr(shell(s, 0, a, dots)[0], 35),
            "outgoing_small_fraction": mp.nstr(p/(k-s), 35),
            "energy_residual_eV": mp.nstr(energy_residual, 8),
        }


def report():
    return {
        "scope": "fixed-hypothesis kinematic feasibility; no empirical verdict",
        "global_threshold_absolute_error_eV": "1e-13",
        "witnesses": [witness(k, variant, n) for k in MOMENTA for variant in VARIANTS for n in DIRECTIONS],
    }


if __name__ == "__main__":
    import json
    print(json.dumps(report(), indent=2, sort_keys=True))

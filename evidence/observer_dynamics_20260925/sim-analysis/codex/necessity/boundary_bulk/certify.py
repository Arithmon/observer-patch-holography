"""Numerical certificates for exact boundary/bulk statements.

The all-k positivity and all-radii uniqueness claims have analytical proofs
in REPORT.md. Finite grids below are diagnostics, not proofs by sampling.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import eval_legendre, gammaln, spherical_jn

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def bessel_square_exact(ell: int, theta: float = 0) -> float:
    """Integral x^(-1-theta)*j_l(x)^2 dx over (0,infinity)."""
    if int(ell) != ell or ell < 0 or not (-2 < theta < 2*ell):
        raise ValueError("convergence requires integer ell>=0 and -2<theta<2*ell")
    return float(np.exp(0.5*math.log(math.pi)-math.log(4)
                        +gammaln(1+theta/2)-gammaln(1.5+theta/2)
                        +gammaln(ell-theta/2)-gammaln(ell+2+theta/2)))


def bessel_tail_envelope(ell: int, x: float) -> float:
    """B_l(x) with |j_l(y)|<=B_l(x)/y for all y>=x>0."""
    if int(ell) != ell or ell < 0 or x <= 0:
        raise ValueError("invalid envelope arguments")
    return sum(math.factorial(ell+m)/(math.factorial(m)*math.factorial(ell-m)*(2*x)**m)
               for m in range(ell+1))


def panel_nodes(a: float, b: float, *, order: int = 32, width: float = math.pi):
    count = int(math.ceil((b-a)/width))
    edges = np.linspace(a, b, count+1)
    z, w = leggauss(order)
    mid, half = (edges[:-1]+edges[1:])/2, np.diff(edges)/2
    return (mid[:, None]+half[:, None]*z).ravel(), (half[:, None]*w).ravel()


def bessel_check(ell: int, theta: float, *, xmax: float = 4096, order: int = 32):
    exact = bessel_square_exact(ell, theta)
    # Adaptive quadrature resolves any noninteger power near zero. The rest
    # uses fixed Gauss panels shorter than one oscillation period.
    first, first_err = quad(lambda x: x**(-1-theta)*spherical_jn(ell, x)**2,
                            0, math.pi, epsabs=2e-13, epsrel=2e-13)
    x, w = panel_nodes(math.pi, xmax, order=order)
    finite = first + float(np.dot(w, x**(-1-theta)*spherical_jn(ell, x)**2))
    tail_bound = bessel_tail_envelope(ell, xmax)**2*xmax**(-2-theta)/(2+theta)
    return {"ell": ell, "theta": theta, "exact": exact,
            "finite_quadrature": finite, "first_panel_error_estimate": first_err,
            "exact_minus_finite": exact-finite, "analytical_positive_tail_bound": tail_bound}


def bump(r):
    r = np.asarray(r, dtype=float)
    u = r-3
    return np.where((u >= 0) & (u <= 1), u*u*(1-u)**2, 0.0)


def bump_fourier(k):
    """3D radial transform H(k), computed without any bulk power prior."""
    k = np.asarray(k, dtype=float)
    if np.any(~np.isfinite(k)) or np.any(k < 0):
        raise ValueError("wave numbers must be nonnegative and finite")
    flat = k.ravel()
    result = np.empty_like(flat)
    small = flat < 4
    # Direct nonoscillatory quadrature at low k avoids endpoint cancellation.
    if np.any(small):
        z, w = leggauss(64)
        r, weights = 3.5+z/2, w/2
        result[small] = 4*np.pi*(np.sinc(flat[small, None]*r/np.pi) @ (weights*r*r*bump(r)))
    if np.any(~small):
        q = flat[~small]
        iq = 1j*q
        # Exact endpoint expansion of integral (r*h(r))*exp(ikr) dr.
        upper = 8/iq**3-54/iq**4+144/iq**5-120/iq**6
        lower = 6/iq**3+30/iq**4+24/iq**5-120/iq**6
        integral = np.exp(4j*q)*upper-np.exp(3j*q)*lower
        result[~small] = 4*np.pi/q*integral.imag
    return result.reshape(k.shape)


def bump_shell_cl(ell: int, radius: float):
    """Angular spectrum perturbation from the exact real-space chord kernel."""
    if ell < 0 or int(ell) != ell or radius <= 0:
        raise ValueError("invalid shell arguments")
    # Break where the support boundary reaches the shell, for accurate Gauss
    # quadrature even when only a small angular cap intersects the bump.
    points = sorted({-1.0, 1.0, *[1-r*r/(2*radius*radius) for r in (3, 4)
                                 if -1 < 1-r*r/(2*radius*radius) < 1]})
    total = 0.0
    for a, b in zip(points[:-1], points[1:]):
        value = quad(lambda mu: float(bump(radius*math.sqrt(2-2*mu)))*eval_legendre(ell, mu),
                     a, b, epsabs=1e-12, epsrel=1e-12)[0]
        total += value
    return 2*math.pi*total


def bump_shell_fourier_cl(ell: int, radius: float, *, kmax: float = 1024):
    k, w = panel_nodes(0, kmax, order=32, width=math.pi/4)
    estimate = 2/math.pi*float(np.dot(w, k*k*bump_fourier(k)*spherical_jn(ell, k*radius)**2))
    # Three integrations by parts bound |H(k)|<=304*pi/k^4:
    # |g''(3)|+|g''(4)|+integral |g'''| <= 6+8+30+12+20 = 76.
    tail = 608*bessel_tail_envelope(ell, kmax*radius)**2/(3*radius**2*kmax**3)
    return estimate, tail


def log_shell_kernel(mu, amplitude: float = 1):
    mu = np.asarray(mu, dtype=float)
    if np.any(mu < -1) or np.any(mu >= 1):
        raise ValueError("use -1<=mu<1; the coincident log kernel is singular")
    return -amplitude/2*(1+np.log((1-mu)/2))


def anisotropy_fixture(increment: float = 0.4):
    points = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], dtype=float)
    adjacency = (np.abs(points@points.T) < 0.5).astype(float)
    lap0 = np.diag(adjacency.sum(axis=1))-adjacency
    b = np.zeros(6); b[0], b[2] = 1, -1
    lap = lap0 + increment*np.outer(b,b)
    cov0, cov = np.linalg.pinv(lap0, hermitian=True), np.linalg.pinv(lap, hermitian=True)
    # A genuine octahedral rotation, cyclically permuting x,y,z.
    perm = np.array([2,3,4,5,0,1])
    rotation = np.eye(6)[:, perm]
    lhs = cov@rotation-rotation@cov
    rhs = -cov@(lap@rotation-rotation@lap)@cov
    dipoles = points/math.sqrt(2)
    return {"baseline_rotation_commutator_norm": float(np.linalg.norm(cov0@rotation-rotation@cov0)),
            "perturbed_precision_commutator_norm": float(np.linalg.norm(lap@rotation-rotation@lap)),
            "perturbed_covariance_commutator_norm": float(np.linalg.norm(lhs)),
            "inverse_commutator_identity_residual": float(np.linalg.norm(lhs-rhs)),
            "dipole_covariance_eigenvalues": np.linalg.eigvalsh(dipoles.T@cov@dipoles).tolist(),
            "covariance_annihilates_constant_residual": float(np.linalg.norm(cov@np.ones(6)))}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_receipt():
    spec = json.loads((HERE/"spec.json").read_text())
    checks = []
    for ell in spec["bessel_checks"]["ell"]:
        for theta in spec["bessel_checks"]["theta"]:
            row = bessel_check(ell, theta, xmax=spec["bessel_checks"]["x_max"], order=32)
            low_order = bessel_check(ell, theta, xmax=spec["bessel_checks"]["x_max"], order=16)
            row["quadrature_order_difference"] = abs(row["finite_quadrature"]-low_order["finite_quadrature"])
            checks.append(row)
    cf = spec["positive_counterexample"]
    grid = cf["diagnostic_k_grid"]
    k = np.geomspace(grid["min"], grid["max"], grid["points"])
    spectral_delta_fraction = k**3*bump_fourier(k)/(240*np.pi)
    shell_checks = []
    for radius in cf["shell_checks"]["radii"]:
        for ell in cf["shell_checks"]["ell"]:
            real = bump_shell_cl(ell, radius)
            fourier, tail = bump_shell_fourier_cl(ell, radius, kmax=cf["shell_checks"]["k_max"])
            shell_checks.append({"radius": radius, "ell": ell,
                                 "real_space_delta_cl_for_unit_bump": real,
                                 "finite_fourier_delta_cl_for_unit_bump": fourier,
                                 "absolute_difference": abs(real-fourier),
                                 "analytical_fourier_tail_bound": tail})
    legendre_checks = []
    for ell in (1,2,3,6,12):
        numerical = 2*np.pi*quad(lambda mu: float(log_shell_kernel(mu))*eval_legendre(ell, mu),
                                 -1, 1, epsabs=1e-11, epsrel=1e-11)[0]
        exact = 2*np.pi/(ell*(ell+1))
        legendre_checks.append({"ell": ell, "numerical_cl": numerical, "exact_cl": exact,
                                "absolute_error": abs(numerical-exact)})
    return {"schema": "oph.boundary-bulk-necessity.receipt.v1",
            "source_sha256": sha256(__file__), "spec_sha256": sha256(HERE/"spec.json"),
            "read_only_source_hashes": {p: sha256(ROOT/p) for p in spec["source_files_read_only"]},
            "bessel_identity_checks": checks, "log_kernel_legendre_checks": legendre_checks,
            "positive_counterexample": {"analytical_minimum_alternative_over_baseline": 0.5,
                                        "sampled_max_abs_spectral_perturbation_fraction": float(np.max(np.abs(spectral_delta_fraction))),
                                        "sampled_min_either_alternative_over_baseline": float(1-np.max(np.abs(spectral_delta_fraction))),
                                        "unit_bump_at_r_3p5": float(bump(3.5)),
                                        "covariance_perturbation_at_r_3p5_for_A1": math.pi/120*float(bump(3.5)),
                                        "shell_checks": shell_checks},
            "anisotropy_fixture": anisotropy_fixture(spec["anisotropy_fixture"]["edge_increment"]),
            "interpretation": "Analytical identities and counterexamples with numerical certificates. No cosmological fit, no native bulk-gluing certificate, no measured tilt."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    value = build_receipt()
    path = HERE/"receipt.json"
    if args.verify:
        if json.loads(path.read_text()) != value:
            raise SystemExit("receipt differs")
        print("Verified boundary/bulk receipt")
    else:
        path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
        print("Wrote boundary/bulk receipt")


if __name__ == "__main__":
    main()

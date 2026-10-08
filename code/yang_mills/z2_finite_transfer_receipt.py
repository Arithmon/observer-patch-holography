#!/usr/bin/env python3
"""Finite Z2 lattice-gauge diagnostic for the ground-state-transform receipt.

The Yang--Mills gap paper assumes a finite ground-state-transform and
cross-fiber receipt: a unitary ``U_r`` with ``U_r Omega_r = 1`` such that
``U_r H_r U_r^{-1} = sum_C D_C`` with each ``D_C = c_C (I - E_C)`` a heat-bath
collar generator whose rate ``c_C`` is independent of the repaired value.
This module evaluates that receipt numerically on two small gauge systems:
Z2 lattice gauge theory on an ``L x L`` periodic spatial torus in the
gauge-invariant (Gauss-law) sector.  The free control has an analytic exact
identity; interacting calculations use float64 arithmetic with explicit
resolution checks. These checks are numerical safeguards, not interval proofs.

Two transfer objects are tested.

* ``wilson``: the reflection-positive Wilson transfer matrix
  ``T = exp(beta_s P / 2) K exp(beta_s P / 2)`` with ``K(s, s') =
  prod_l exp(beta_t s_l s'_l)``, and ``H = -log(T / lambda_max)``.
* ``kogut_susskind``: the Hamiltonian ``H = -lam sum_l X_l - sum_p U_p``,
  shifted so that its ground energy is zero.

For each, ``Omega`` is the Perron ground state, ``pi = Omega^2`` the
stationary time-zero law, ``L = D_Omega^{-1} H D_Omega`` the ground-state
(Doob) transform, and ``E_l`` the ``pi``-preserving conditional expectation
on the single-link collar fiber ``{o, X_l o}``.

The script reports

* the relative Frobenius residual of the best constant-rate fit
  ``L ~ sum_l c_l (I - E_l)``;
* the exact fiber-dependent rates ``c_l(o)`` that reproduce ``L`` when ``L``
  is single-flip, and their spread ``max/min`` (the cross-fiber receipt
  requires spread 1);
* the Dobrushin influence ``eta_*`` of ``pi`` under single-link heat bath,
  the floor ``delta_* = c_*(1 - eta_*)``, and the exact spectral gaps of
  ``H`` and of the unit-rate heat-bath generator.

Everything here is a finite diagnostic on a toy gauge system.  It is not a
physical compact-simple-gauge receipt, and the output JSON says so.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg.lapack import dgejsv

SCHEMA = "oph.yang_mills.z2_finite_transfer_receipt.v2"
RESOLUTION_RTOL = 1e-7
EPS = np.finfo(float).eps


# ----------------------------------------------------------------------------
# Lattice geometry
# ----------------------------------------------------------------------------


def require_nondegenerate_periodic_lattice(L: int) -> None:
    """Reject the self-loop ``L = 1`` cellulation.

    The bit-mask implementation below stores an incidence set with bitwise OR.
    On the physical one-site periodic cellulation each link is incident to the
    same site twice and each plaquette traverses each link twice, so those
    incidences cancel over Z2.  The OR representation is faithful only for the
    non-self-loop tori used by the committed receipt.
    """
    if L < 2:
        raise ValueError(
            "the physical periodic Z2 receipt requires L >= 2; at L = 1 "
            "repeated self-loop incidences cancel and cannot be represented "
            "by these OR masks"
        )


def link_index(L: int, x: int, y: int, direction: int) -> int:
    """Index of the link leaving site ``(x, y)`` in ``direction`` (0 = x, 1 = y)."""
    return 2 * ((x % L) * L + (y % L)) + direction


def plaquette_masks(L: int) -> list[int]:
    require_nondegenerate_periodic_lattice(L)
    masks = []
    for x in range(L):
        for y in range(L):
            m = 0
            m |= 1 << link_index(L, x, y, 0)
            m |= 1 << link_index(L, x + 1, y, 1)
            m |= 1 << link_index(L, x, y + 1, 0)
            m |= 1 << link_index(L, x, y, 1)
            masks.append(m)
    return masks


def star_masks(L: int) -> list[int]:
    """Coboundary of each site: links with exactly that site as an endpoint."""
    require_nondegenerate_periodic_lattice(L)
    masks = []
    for x in range(L):
        for y in range(L):
            m = 0
            m |= 1 << link_index(L, x, y, 0)
            m |= 1 << link_index(L, x, y, 1)
            m |= 1 << link_index(L, x - 1, y, 0)
            m |= 1 << link_index(L, x, y - 1, 1)
            masks.append(m)
    return masks


def gauge_group_masks(L: int) -> np.ndarray:
    """All gauge transformations as XOR masks (one per subset of sites)."""
    stars = star_masks(L)
    n_sites = len(stars)
    masks = np.zeros(1 << n_sites, dtype=np.int64)
    for subset in range(1 << n_sites):
        m = 0
        for site in range(n_sites):
            if subset >> site & 1:
                m ^= stars[site]
        masks[subset] = m
    return np.unique(masks)


def popcount(values: np.ndarray) -> np.ndarray:
    return np.bitwise_count(values.astype(np.uint64)).astype(np.int64)


# ----------------------------------------------------------------------------
# Orbit space
# ----------------------------------------------------------------------------


class Z2GaugeOrbits:
    """Gauge orbits of Z2 link configurations on the L x L periodic torus."""

    def __init__(self, L: int) -> None:
        require_nondegenerate_periodic_lattice(L)
        self.L = L
        self.n_links = 2 * L * L
        self.n_configs = 1 << self.n_links
        self.plaquettes = plaquette_masks(L)
        self.gauge = gauge_group_masks(L)
        configs = np.arange(self.n_configs, dtype=np.int64)
        # canonical representative: minimum over the gauge orbit
        reps = configs.copy()
        for g in self.gauge:
            reps = np.minimum(reps, configs ^ int(g))
        self.rep_of_config = reps
        self.reps, inverse = np.unique(reps, return_inverse=True)
        self.orbit_of_config = inverse
        self.n_orbits = len(self.reps)
        self.orbit_size = np.bincount(inverse, minlength=self.n_orbits)
        if not np.all(self.orbit_size == self.orbit_size[0]):
            raise RuntimeError("gauge action is not free; weighted orbit basis needed")
        # members[o] = all configurations in orbit o
        self.members = (self.reps[:, None] ^ self.gauge[None, :]).astype(np.int64)
        # plaquette sum P(o) = sum_p U_p for the representative
        P = np.zeros(self.n_orbits, dtype=np.int64)
        for m in self.plaquettes:
            parity = popcount(self.reps & m) & 1
            P += 1 - 2 * parity
        self.plaquette_sum = P
        # single-link flip as an involution on orbits
        self.flip = np.zeros((self.n_links, self.n_orbits), dtype=np.int64)
        for l in range(self.n_links):
            self.flip[l] = self.orbit_of_config[self.reps ^ (1 << l)]
            if np.any(self.flip[l] == np.arange(self.n_orbits)):
                raise RuntimeError("single-link flip has a fixed orbit")

    # -- transfer objects ---------------------------------------------------

    def wilson_transfer(self, beta_s: float, beta_t: float) -> np.ndarray:
        """Symmetric gauge-invariant Wilson transfer matrix on orbit space."""
        n = self.n_orbits
        T = np.zeros((n, n))
        weight = np.exp(0.5 * beta_s * self.plaquette_sum)
        for o in range(n):
            d = popcount(self.members ^ self.reps[o])  # (n_orbits, |G|)
            kin = np.exp(beta_t * (self.n_links - 2 * d)).sum(axis=1)
            T[o] = weight[o] * kin * weight
        return T

    def kogut_susskind(self, lam: float) -> np.ndarray:
        n = self.n_orbits
        H = np.zeros((n, n))
        H[np.arange(n), np.arange(n)] = -self.plaquette_sum.astype(float)
        for l in range(self.n_links):
            H[np.arange(n), self.flip[l]] -= lam
        return H


# ----------------------------------------------------------------------------
# Receipt evaluation
# ----------------------------------------------------------------------------


def _symmetric_finite_matrix(matrix: np.ndarray) -> np.ndarray:
    if np.ma.isMaskedArray(matrix):
        raise ValueError("masked matrices do not specify a complete transfer operator")
    array = np.asarray(matrix)
    if (array.ndim != 2 or array.shape[0] != array.shape[1] or len(array) < 2
            or np.iscomplexobj(array) or array.dtype.kind not in "fiu"
            or not np.isfinite(array).all()):
        raise ValueError("expected a finite real square matrix of dimension >= 2")
    array = array.astype(float)
    scale = float(np.max(np.abs(array)))
    if scale == 0 or not np.isfinite(scale):
        raise ValueError("matrix has no finite nonzero scale")
    if np.max(np.abs(array / scale - array.T / scale)) > 8 * EPS:
        raise ValueError("matrix must be symmetric")
    return array


def _positive_perron(omega: np.ndarray) -> np.ndarray:
    if omega.sum() < 0:
        omega = -omega
    if not np.isfinite(omega).all() or np.any(omega <= 0):
        raise RuntimeError("Perron support is not resolved as strictly positive")
    return omega / np.linalg.norm(omega)


def symmetric_log_hamiltonian(T: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Logarithm of an already-rounded matrix, only on its resolved spectrum.

    Scaling prevents unit-dependent overflow. A positive computed eigenvalue
    alone does not resolve its logarithm; a resolved bottom eigenvalue alone
    does not resolve the Perron state. Source Wilson calculations below avoid
    rounding away kinetic eigenvalues before this operation in the first place.
    """
    T = _symmetric_finite_matrix(T)
    scale = float(np.max(np.abs(T)))
    normalized = T / scale
    w, V = np.linalg.eigh(normalized)
    error = 4 * len(T) * EPS * np.linalg.norm(normalized, ord=np.inf)
    if w[0] <= error or error / w[0] > RESOLUTION_RTOL:
        raise RuntimeError("transfer log spectrum is unresolved in float64 precision")
    separation = w[-1] - w[-2]
    omega = _positive_perron(V[:, -1])
    if separation <= 2 * error or 2 * error / separation > RESOLUTION_RTOL * omega.min():
        raise RuntimeError("Perron eigenspace is unresolved in float64 precision")
    H = (V * (np.log(w[-1]) - np.log(w))) @ V.T
    maximum = float(w[-1]) * scale
    if not math.isfinite(maximum) or maximum <= 0:
        raise ValueError("transfer normalization is outside float64 range")
    return H, omega, maximum


def _real_parameter(value: float, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a finite real parameter")
    original = value
    try:
        value = float(value)
    except OverflowError as error:
        raise ValueError(f"{name} must be a finite real parameter") from error
    if value == 0 and original != 0:
        raise RuntimeError(f"{name} underflows the supported float64 precision")
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    if isinstance(original, np.floating) and original != 0:
        # Promote the narrowed value back before comparing. A wider input can
        # round to a nonzero binary64 subnormal with large relative damage;
        # checking only whether it became zero misses that loss. Ordinary
        # rounding and exactly representable binary64 subnormals remain valid.
        restored = type(original)(value)
        if abs((restored - original) / original) > 4 * EPS:
            raise RuntimeError(f"{name} loses input precision when converted to float64")
    return value


def dual_coupling(beta_t: float) -> float:
    """log(coth(beta_t)), including the near-identity transfer regime."""
    beta_t = _real_parameter(beta_t, "beta_t", positive=True)
    rate = (-math.log(math.tanh(beta_t)) if beta_t < 1
            else 2 * math.atanh(math.exp(-2 * beta_t)))
    if rate < np.finfo(float).tiny:
        raise RuntimeError("dual coupling is unresolved in float64 range")
    return rate


def wilson_hamiltonian(orbits: Z2GaugeOrbits, beta_s: float,
                       beta_t: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute H from the source factor, without forming its ill-conditioned Gram matrix.

    For electric cycles z in the annihilator of the gauge group, the normalized
    character matrix F[o,z] = (-1)^(rep_o dot z)/sqrt(n) is orthogonal and
    K = (2 cosh(beta_t))^E F diag(tanh(beta_t)^|z|) F.T.
    Thus T is a scalar times A A.T, A = D_s F diag(tanh(beta_t)^(|z|/2)).
    LAPACK GEJSV with JOBA='E' retains relative singular-value accuracy under
    column scaling; its conditioning depends on D_s F, not on the tiny kinetic
    eigenvalues. See https://www.netlib.org/lapack/explore-html/d8/d78/
    group__gejsv_gaca7ba7f1e8002c7a1d5bffa4ccbb541f.html .
    """
    beta_s = _real_parameter(beta_s, "beta_s")
    beta_t = _real_parameter(beta_t, "beta_t", positive=True)
    rate = dual_coupling(beta_t)
    n = orbits.n_orbits
    log_scale = orbits.n_links * (beta_t + math.log1p(math.exp(-2 * beta_t)))
    if beta_s == 0:
        # Exact free law, not a threshold applied to nearby interacting inputs.
        identity = np.eye(n)
        H = (rate / 2) * sum(identity - identity[flip] for flip in orbits.flip)
        omega = np.full(n, 1 / math.sqrt(n))
    else:
        log_condition = 0.5 * abs(beta_s) * int(orbits.plaquette_sum.max() - orbits.plaquette_sum.min())
        if log_condition > math.log(RESOLUTION_RTOL / (8 * n * EPS)):
            raise RuntimeError("Wilson spatial factor condition is unresolved in float64")
        row_logs = 0.5 * beta_s * orbits.plaquette_sum
        shift = float(row_logs.max())
        # Gauge-invariant electric characters, enumerated before any floating
        # transfer arithmetic. No selection on diagnostic success occurs here.
        cycles = np.arange(orbits.n_configs, dtype=np.int64)
        for star in star_masks(orbits.L):
            cycles = cycles[(popcount(cycles & star) & 1) == 0]
        if len(cycles) != n:
            raise RuntimeError("electric character dimension disagrees with orbit space")
        F = (1 - 2 * (popcount(orbits.reps[:, None] & cycles[None, :]) & 1)) / math.sqrt(n)
        column_scale = np.exp(-0.5 * rate * popcount(cycles))
        if column_scale.min() < np.finfo(float).tiny ** 0.5:
            raise RuntimeError("Wilson kinetic factor is unresolved in float64 range")
        A = np.exp(row_logs - shift)[:, None] * F * column_scale[None, :]
        # SciPy enum mapping: E,U,V,R,N,N. Require all columns to survive the
        # restricted range; do not transpose the column-scaled factor or perturb it.
        s, U, V, work, rank, info = dgejsv(
            A, joba=1, jobu=0, jobv=0, jobr=1, jobt=0, jobp=0,
        )
        if (info != 0 or rank[0] != n or rank[1] != n or rank[2] != 0
                or not np.isfinite(s).all() or s[-1] <= 0
                or not np.isfinite(U).all() or not np.isfinite(V).all()
                or not np.isfinite(work[:3]).all() or np.any(work[:3] <= 0)):
            raise RuntimeError("Wilson factor singular spectrum is unresolved")
        singular_scale = work[0] / work[1]
        reconstructed = ((U * s) @ V.T) * singular_scale
        relative_columns = np.linalg.norm((reconstructed - A) / np.linalg.norm(A, axis=0), axis=0)
        if relative_columns.max() * math.exp(log_condition) > RESOLUTION_RTOL:
            raise RuntimeError("Wilson factor column reconstruction is unresolved")
        separation = (s[0] - s[1]) / s[0]
        if separation <= 8 * n * EPS / RESOLUTION_RTOL:
            raise RuntimeError("Wilson Perron separation is unresolved in float64")
        omega = _positive_perron(U[:, 0])
        energies = 2 * (np.log(s[0]) - np.log(s))
        H = (U * energies) @ U.T
        # A good normwise H is insufficient for the subsequent pointwise Doob
        # division. Given an entrywise uncertainty estimate eps*n*||H||_2,
        # its Frobenius amplification is ||1/omega||_2 because ||omega||_2=1.
        # This is a resolution policy, not a certified SVD forward-error bound.
        # Scale by min(omega) to evaluate the ratio without reciprocal overflow.
        smallest = float(omega.min())
        scaled_doob = H * omega[None, :] * (smallest / omega)[:, None]
        doob_uncertainty = (EPS * n * float(energies.max())
                            * np.linalg.norm(smallest / omega)
                            / np.linalg.norm(scaled_doob))
        if not math.isfinite(doob_uncertainty) or doob_uncertainty > RESOLUTION_RTOL:
            raise RuntimeError("Wilson Doob transform precision is unresolved on the Perron support")
        log_scale += 2 * (shift + math.log(s[0]) + math.log(singular_scale))
    if log_scale > math.log(np.finfo(float).max):
        raise ValueError("reported transfer normalization is outside float64 range")
    maximum = math.exp(log_scale)
    if not math.isfinite(maximum) or maximum <= 0:
        raise ValueError("reported transfer normalization is outside float64 range")
    return H, omega, maximum


def ground_state(H: np.ndarray) -> tuple[np.ndarray, float]:
    """Resolve the positive ground state before using its component ratios.

    Conservation alone does not distinguish accurate probabilities from an
    unresolved mixture of nearly degenerate ground-sector eigenvectors. This
    conservative float-resolution policy can refuse otherwise accurate inputs.
    """
    H = _symmetric_finite_matrix(H)
    scale = float(np.max(np.abs(H)))
    normalized = H / scale
    w, V = np.linalg.eigh(normalized)
    omega = _positive_perron(V[:, 0])
    error = 4 * len(H) * EPS * np.linalg.norm(normalized, ord=np.inf)
    separation = w[1] - w[0]
    if separation <= 2 * error or 2 * error / separation > RESOLUTION_RTOL * omega.min():
        raise RuntimeError("ground-state Perron support precision is unresolved")
    energy = float(w[0]) * scale
    if not math.isfinite(energy):
        raise ValueError("ground energy is outside float64 range")
    return omega, energy


def doob_transform(H: np.ndarray, omega: np.ndarray, e0: float) -> np.ndarray:
    result = (H - e0 * np.eye(len(omega))) * omega[None, :] / omega[:, None]
    if not np.isfinite(result).all():
        raise RuntimeError("Doob transform precision is unresolved")
    # H omega = e0 omega implies conservation exactly. Test every row relative
    # to its own magnitude, so a large unrelated row cannot hide a violation.
    # Scaling first also makes the test independent of generator units.
    row_scale = np.max(np.abs(result), axis=1, keepdims=True)
    scaled = np.divide(result, row_scale, out=np.zeros_like(result), where=row_scale != 0)
    if np.any(np.abs(scaled.sum(axis=1)) > RESOLUTION_RTOL * np.abs(scaled).sum(axis=1)):
        raise RuntimeError("Doob row conservation is unresolved at the numerical resolution")
    return result


def heat_bath_projectors(orbits: Z2GaugeOrbits, pi: np.ndarray) -> list[np.ndarray]:
    n = orbits.n_orbits
    projectors = []
    for l in range(orbits.n_links):
        partner = orbits.flip[l]
        denom = pi + pi[partner]
        E = np.zeros((n, n))
        idx = np.arange(n)
        E[idx, idx] = pi / denom
        E[idx, partner] = pi[partner] / denom
        projectors.append(E)
    return projectors


def constant_rate_fit(Lgen: np.ndarray, projectors: list[np.ndarray]) -> dict[str, Any]:
    n = Lgen.shape[0]
    basis = np.stack([(np.eye(n) - E).ravel() for E in projectors], axis=1)
    scale = float(np.max(np.abs(Lgen)))
    if not math.isfinite(scale) or scale == 0:
        raise ValueError("rate fit requires a finite nonzero generator")
    target = Lgen.ravel() / scale
    coeff, *_ = np.linalg.lstsq(basis, target, rcond=None)
    residual = target - basis @ coeff
    rel = float(np.linalg.norm(residual) / np.linalg.norm(target))
    return {
        "rates": [float(c * scale) for c in coeff],
        "rate_min": float(coeff.min() * scale),
        "rate_max": float(coeff.max() * scale),
        "relative_frobenius_residual": rel,
    }


def fiber_dependent_rates(
    Lgen: np.ndarray, orbits: Z2GaugeOrbits, pi: np.ndarray
) -> dict[str, Any]:
    """Exact collar rates ``c_l(o)`` from the off-diagonal of ``L``.

    For ``o' = X_l o`` the heat-bath form gives
    ``L(o, o') = -c_l(o) pi(o') / (pi(o) + pi(o'))``.  The returned spread is
    ``max c / min c`` over all links and orbits; the cross-fiber receipt needs
    spread exactly one.  Off-diagonal mass of ``L`` outside single-link flips is
    returned separately; it must vanish for the single-flip form to be exact.
    """
    n = orbits.n_orbits
    idx = np.arange(n)
    rates = np.zeros((orbits.n_links, n))
    single_flip_mask = np.zeros((n, n), dtype=bool)
    for l in range(orbits.n_links):
        partner = orbits.flip[l]
        off = Lgen[idx, partner]
        rates[l] = -off * (pi + pi[partner]) / pi[partner]
        single_flip_mask[idx, partner] = True
    off_diag = Lgen - np.diag(np.diag(Lgen))
    outside = off_diag[~single_flip_mask]
    return {
        "rate_min": float(rates.min()),
        "rate_max": float(rates.max()),
        "spread_max_over_min": float(rates.max() / rates.min()) if rates.min() > 0 else math.inf,
        "all_rates_positive": bool(np.all(rates > 0)),
        "offdiagonal_mass_outside_single_flip": float(np.abs(outside).sum()),
        "offdiagonal_mass_single_flip": float(np.abs(off_diag[single_flip_mask]).sum()),
    }


def dobrushin_influence(orbits: Z2GaugeOrbits, pi: np.ndarray) -> dict[str, Any]:
    """Total-variation influence matrix of the single-link conditional kernels."""
    n_links = orbits.n_links
    idx = np.arange(orbits.n_orbits)
    # conditional probability that link l is in its representative state at o:
    # kernel(o) = pi(o) / (pi(o) + pi(X_l o)), as a function of the other links.
    kernel = np.zeros((n_links, orbits.n_orbits))
    for l in range(n_links):
        kernel[l] = pi / (pi + pi[orbits.flip[l]])
    influence = np.zeros((n_links, n_links))
    for l in range(n_links):
        for u in range(n_links):
            if u == l:
                continue
            # changing link u: compare kernel at o and at X_u o.  The kernel is
            # stated for the representative state of link l, which X_u does not
            # change, so the TV distance is the absolute difference.
            diff = np.abs(kernel[l] - kernel[l][orbits.flip[u]])
            influence[l, u] = float(diff.max())
    row_sums = influence.sum(axis=1)
    return {
        "eta_star": float(row_sums.max()),
        "row_sums": [float(r) for r in row_sums],
    }


def spectral_gap(M: np.ndarray, pi: np.ndarray | None = None) -> float:
    """Smallest nonzero eigenvalue of a generator symmetrisable by ``pi``."""
    if pi is not None:
        s = np.sqrt(pi)
        M = M * s[:, None] / s[None, :]
    w = np.linalg.eigvalsh(0.5 * (M + M.T))
    w = np.sort(w)
    return float(w[1])


def evaluate(orbits: Z2GaugeOrbits, transfer: str, **params: float) -> dict[str, Any]:
    expected = {"wilson": {"beta_s", "beta_t"}, "kogut_susskind": {"lam"}}
    if transfer not in expected or set(params) != expected[transfer]:
        raise ValueError("transfer parameters must exactly match the selected operator")
    if transfer == "wilson":
        H, omega, lam_max = wilson_hamiltonian(orbits, params["beta_s"], params["beta_t"])
        e0 = 0.0
        extra = {"lambda_max": lam_max}
    elif transfer == "kogut_susskind":
        _real_parameter(params["lam"], "lam", positive=True)
        H = orbits.kogut_susskind(params["lam"])
        omega, e0 = ground_state(H)
        extra = {"ground_energy": e0}
    else:
        raise ValueError(transfer)
    pi = omega**2
    pi = pi / pi.sum()
    Lgen = doob_transform(H, omega, e0)
    projectors = heat_bath_projectors(orbits, pi)
    unit_heat_bath = sum(np.eye(orbits.n_orbits) - E for E in projectors)
    fit = constant_rate_fit(Lgen, projectors)
    fibre = fiber_dependent_rates(Lgen, orbits, pi)
    dob = dobrushin_influence(orbits, pi)
    eta = dob["eta_star"]
    # Subtracting almost equal conditionals can lose the influence even when
    # the matrix, populations and conservation are individually well resolved.
    # This dimension-scaled contrast floor is a numerical safeguard, not an
    # interval guarantee for every derived observable. Only the source-exact
    # free law justifies bypassing it with a known zero influence.
    source_free = transfer == "wilson" and params["beta_s"] == 0
    contrast_floor = 8 * EPS * orbits.n_orbits * orbits.n_links
    if not source_free and contrast_floor > RESOLUTION_RTOL * eta:
        raise RuntimeError("Dobrushin influence contrast is unresolved at the numerical precision")
    result: dict[str, Any] = {
        "transfer": transfer,
        "parameters": params,
        "n_orbits": int(orbits.n_orbits),
        "n_links": int(orbits.n_links),
        "doob_generator_rows_sum_zero": True,  # checked per row in doob_transform
        "doob_generator_offdiagonal_nonpositive": bool(
            np.all(Lgen - np.diag(np.diag(Lgen)) <= 1e-12)
        ),
        "constant_rate_fit": fit,
        "fiber_dependent_rates": fibre,
        "dobrushin": {
            "eta_star": eta,
            "dobrushin_condition_holds": bool(eta < 1),
            "unit_rate_floor_c_star_times_1_minus_eta": float(max(0.0, 1 - eta)),
        },
        "spectral": {
            "gap_H": spectral_gap(Lgen, pi),
            "gap_unit_rate_heat_bath": spectral_gap(unit_heat_bath, pi),
            "pi_min": float(pi.min()),
            "pi_max": float(pi.max()),
        },
    }
    if transfer == "kogut_susskind":
        analytic_floor = 2.0 * params["lam"]
        result["variable_rate_floor"] = {
            "identity": "c_l(o) = lambda * (r_l(o) + 1/r_l(o))",
            "analytic_lower_bound": "c_l(o) >= 2 * lambda by AM-GM",
            "lower_bound_value": analytic_floor,
            "numerical_min_respects_bound": bool(
                fibre["rate_min"] >= analytic_floor - 1e-10
            ),
            "scope": (
                "finite Kogut-Susskind Doob transform; quotient-space "
                "approximate tensorization and continuum transfer not proved"
            ),
        }
    result.update(extra)
    return result


def run(L_values: list[int], betas: list[float], lams: list[float]) -> dict[str, Any]:
    runs = []
    for L in L_values:
        orbits = Z2GaugeOrbits(L)
        runs.append({"L": L, "transfer": "wilson", "control": "beta_s_zero",
                     **evaluate(orbits, "wilson", beta_s=0.0, beta_t=0.5)})
        for beta in betas:
            runs.append({"L": L, **evaluate(orbits, "wilson", beta_s=beta, beta_t=beta)})
        for lam in lams:
            runs.append({"L": L, **evaluate(orbits, "kogut_susskind", lam=lam)})
    receipt = {
        "schema": SCHEMA,
        "scope": "finite_gauge_diagnostic",
        "physical_clay_receipt": False,
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "grid_scope": {
            "L": L_values,
            "wilson_diagonal_beta_s_eq_beta_t": betas,
            "kogut_susskind_lambda": lams,
            "universal_no_go": False,
        },
        "system": "Z2 lattice gauge theory, L x L periodic spatial torus, gauge-invariant sector",
        "receipt_under_test": (
            "finite ground-state-transform and cross-fiber receipt: "
            "U_r H_r U_r^{-1} = sum_C c_C (I - E_C) with c_C independent of the repaired value"
        ),
        "ground_state_transform": "Doob transform by the Perron vector, pi = Omega^2",
        "collars": "one collar per spatial link, fiber {o, X_l o}, pi-preserving heat bath",
        "runs": runs,
    }
    receipt["sha256_of_runs"] = hashlib.sha256(
        json.dumps(runs, sort_keys=True).encode("utf-8")
    ).hexdigest()
    return receipt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--L", type=int, nargs="+", default=[2, 3])
    ap.add_argument("--beta", type=float, nargs="+", default=[0.1, 0.3, 0.5, 0.7, 1.0])
    ap.add_argument("--lam", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    receipt = run(args.L, args.beta, args.lam)
    text = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)
    if args.output:
        args.output.write_bytes((text + "\n").encode("utf-8"))
    for r in receipt["runs"]:
        tag = r["transfer"] + (" control" if "control" in r else "")
        print(
            f"L={r['L']} {tag:24s} params={r['parameters']} "
            f"fit_resid={r['constant_rate_fit']['relative_frobenius_residual']:.3e} "
            f"spread={r['fiber_dependent_rates']['spread_max_over_min']:.4f} "
            f"outside_single_flip={r['fiber_dependent_rates']['offdiagonal_mass_outside_single_flip']:.3e} "
            f"eta*={r['dobrushin']['eta_star']:.4f} "
            f"gapH={r['spectral']['gap_H']:.4f} gapHB={r['spectral']['gap_unit_rate_heat_bath']:.4f}"
        )


if __name__ == "__main__":
    main()

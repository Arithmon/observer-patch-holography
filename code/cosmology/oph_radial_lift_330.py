"""Source-side radial-lift reference implementation for OPH issue #330.

The module implements the mathematical receipts proved in
``RADIAL_LIFT_THEOREMS_330.tex``:

* exact shell/window projection from a homogeneous isotropic 3-D curvature field;
* exact Mellin integral for the thin-shell power-law branch;
* source-derived amplitude conversion (never an angular-spectrum fit);
* a rigorous finite-window stability bound in a weighted Bessel Hilbert space;
* finite-basis SVD/null-space and prior-selected continuation diagnostics;
* the explicit dilation-intertwiner residual needed to derive the radial power law;
* fail-closed receipt construction and the TT/TE/EE promotion firewall.

This module does not run a Boltzmann solver and contains no observational target.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import solve_triangular
from scipy.special import digamma, gammaln, spherical_jn


class RadialLiftInputError(ValueError):
    """Raised when an input violates the theorem contract."""


@dataclass(frozen=True)
class PrimordialAmplitude:
    """Exact source-side thin-shell amplitude conversion."""

    A_q: float
    theta: float
    A_zeta: float
    Z_q: float
    R_star: float
    k_pivot: float
    kR_pivot: float
    conversion_factor: float


@dataclass(frozen=True)
class WindowBound:
    """Floating-point evaluation of the analytic finite-window bound."""

    ell: int
    theta: float
    I_ell_theta: float
    J_ell_theta: float
    eta: float
    shell_norm: float
    absolute_cl_bound: float
    relative_to_shell_cl_bound: float


@dataclass(frozen=True)
class NullSpaceReport:
    shape: tuple[int, int]
    rank: int
    nullity: int
    singular_values: list[float]
    effective_threshold: float
    condition_number_nonzero: float | None
    null_basis: list[list[float]]
    relative_cutoff: float
    rank_metric: str = "raw_operator_euclidean"


@dataclass(frozen=True)
class PriorContinuation:
    """Numerical minimum-Q continuation; rank refers to the whitened operator."""

    p: list[float]
    residual: list[float]
    residual_norm: float
    objective: float
    resolution: list[list[float]]
    null_projector: list[list[float]]
    effective_rank: int
    relative_cutoff: float
    rank_metric: str = "row_equilibrated_prior_whitened_operator"


@dataclass(frozen=True)
class DilationReceipt:
    """Finite safe-band test of the physical dilation-intertwiner law."""

    theta: float
    scale_ratios: list[float]
    max_absolute_log_residual: float
    rms_log_residual: float
    passed: bool
    tolerance: float
    evaluated_pairs: int


# ---------------------------------------------------------------------------
# Exact shell formulas
# ---------------------------------------------------------------------------


def mellin_spherical_bessel_square(ell: int, theta: float) -> float:
    r"""Return ``I_l(theta)=∫ dln(x) x^{-theta} j_l(x)^2`` exactly.

    The formula is

    .. math::

        I_\ell(\theta)=\frac{\sqrt\pi}{4}
        \frac{\Gamma(1+\theta/2)}{\Gamma(3/2+\theta/2)}
        \frac{\Gamma(\ell-\theta/2)}
             {\Gamma(\ell+2+\theta/2)}.

    Absolute convergence holds for ``-2 < theta < 2*ell``.  The OPH retained
    band starts at ``ell=2``, so a common sufficient interval is ``(-2, 4)``.
    """

    ell = _integer_at_least(ell, 0, "ell")
    theta = _finite(theta, "theta")
    if not (-2.0 < theta < 2.0 * ell):
        raise RadialLiftInputError(
            f"Mellin integral requires -2 < theta < 2*ell; got ell={ell}, theta={theta}"
        )
    log_value = (
        0.5 * math.log(math.pi)
        - math.log(4.0)
        + gammaln(1.0 + 0.5 * theta)
        - gammaln(1.5 + 0.5 * theta)
        + gammaln(ell - 0.5 * theta)
        - gammaln(ell + 2.0 + 0.5 * theta)
    )
    value = math.exp(log_value)
    if not (math.isfinite(value) and value > 0.0):
        raise RadialLiftInputError("Mellin integral evaluated to a nonpositive/nonfinite value")
    return value


def derivative_mellin_norm(ell: int, theta: float) -> float:
    r"""Return ``J_l(theta)=∫ dln(x) x^{2-theta} j_l'(x)^2`` exactly.

    Integration by parts with the spherical-Bessel equation gives

    .. math::

        J_\ell(\theta)=I_\ell(\theta-2)
        -\left[\ell(\ell+1)-\frac{\theta(\theta+1)}2\right]I_\ell(\theta).

    The derivative norm is finite for ``0 < theta < 2*ell``.
    """

    ell = _integer_at_least(ell, 1, "ell")
    theta = _finite(theta, "theta")
    if not (0.0 < theta < 2.0 * ell):
        raise RadialLiftInputError("derivative Mellin norm requires 0 < theta < 2*ell")
    # Gamma recurrence reduces the difference to one positive norm times a
    # rational factor.  Computing I(theta - 2) first rounds theta - 2 to -2
    # for small positive theta and needlessly loses the domain endpoint.
    value = mellin_spherical_bessel_square(ell, theta) * (
        ell * (ell + 1.0) / theta + (theta + 1.0) * (theta - 2.0) / 4.0
    )
    if not (math.isfinite(value) and value > 0.0):
        raise RadialLiftInputError("derived derivative norm is negative or nonfinite")
    return value


def screen_gamma_ratio_cl(ell: ArrayLike, A_q: float, theta: float) -> NDArray[np.float64] | float:
    """Return the exact gamma-ratio screen spectrum ``C_l^q``."""

    A_q = _positive(A_q, "A_q")
    theta = _finite(theta, "theta")
    ell_arr = np.asarray(ell, dtype=float)
    if np.any(~np.isfinite(ell_arr)) or np.any(ell_arr < 2.0):
        raise RadialLiftInputError("ell must be finite and at least 2")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("the common retained-band interval is -2 < theta < 4")
    log_values = (
        math.log(A_q)
        + gammaln(ell_arr - 0.5 * theta)
        - gammaln(ell_arr + 2.0 + 0.5 * theta)
    )
    values = np.exp(log_values)
    if np.isscalar(ell):
        return float(values)
    return np.asarray(values, dtype=float)


def primordial_amplitude_from_screen(
    A_q: float,
    theta: float,
    *,
    Z_q: float,
    R_star: float,
    k_pivot: float,
) -> PrimordialAmplitude:
    r"""Convert source amplitude ``A_q`` to ``A_zeta`` without fitting.

    For

    .. math::

        \Delta_\zeta^2(k)=A_\zeta(k/k_*)^{-\theta},\qquad
        q(\hat n)=Z_q\,\zeta(R_*\hat n),

    the exact result is

    .. math::

        A_\zeta=\frac{A_q}{\pi^{3/2} Z_q^2(k_*R_*)^\theta}
        \frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
    """

    A_q = _positive(A_q, "A_q")
    theta = _finite(theta, "theta")
    Z_q = _positive(Z_q, "Z_q")
    R_star = _positive(R_star, "R_star")
    k_pivot = _positive(k_pivot, "k_pivot")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("thin-shell retained-band formula requires -2 < theta < 4")
    kR = k_pivot * R_star
    log_factor = (
        -1.5 * math.log(math.pi)
        - 2.0 * math.log(Z_q)
        - theta * math.log(kR)
        + gammaln(1.5 + 0.5 * theta)
        - gammaln(1.0 + 0.5 * theta)
    )
    factor = math.exp(log_factor)
    return PrimordialAmplitude(
        A_q=A_q,
        theta=theta,
        A_zeta=A_q * factor,
        Z_q=Z_q,
        R_star=R_star,
        k_pivot=k_pivot,
        kR_pivot=kR,
        conversion_factor=factor,
    )


def screen_amplitude_from_primordial(
    A_zeta: float,
    theta: float,
    *,
    Z_q: float,
    R_star: float,
    k_pivot: float,
) -> float:
    """Inverse of :func:`primordial_amplitude_from_screen`."""

    A_zeta = _positive(A_zeta, "A_zeta")
    theta = _finite(theta, "theta")
    Z_q = _positive(Z_q, "Z_q")
    R_star = _positive(R_star, "R_star")
    k_pivot = _positive(k_pivot, "k_pivot")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("thin-shell retained-band formula requires -2 < theta < 4")
    log_factor = (
        1.5 * math.log(math.pi)
        + 2.0 * math.log(Z_q)
        + theta * math.log(k_pivot * R_star)
        + gammaln(1.0 + 0.5 * theta)
        - gammaln(1.5 + 0.5 * theta)
    )
    return A_zeta * math.exp(log_factor)


def thin_shell_powerlaw_cl(
    ell: ArrayLike,
    A_zeta: float,
    theta: float,
    *,
    Z_q: float,
    R_star: float,
    k_pivot: float,
) -> NDArray[np.float64] | float:
    """Exact thin-shell angular spectrum of the source power law."""

    A_q = screen_amplitude_from_primordial(
        A_zeta, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    return screen_gamma_ratio_cl(ell, A_q, theta)


def primordial_amplitude_log_sensitivity(theta: float, kR_pivot: float) -> float:
    r"""Return ``∂_theta ln(A_zeta/A_q)`` at fixed ``Z_q`` and ``k_*R_*``."""

    theta = _finite(theta, "theta")
    kR_pivot = _positive(kR_pivot, "kR_pivot")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("theta must lie in (-2, 4)")
    return float(
        -math.log(kR_pivot)
        + 0.5 * (digamma(1.5 + 0.5 * theta) - digamma(1.0 + 0.5 * theta))
    )


# ---------------------------------------------------------------------------
# Exact finite-window forward kernel and stability bound
# ---------------------------------------------------------------------------


def normalized_radial_window(radii: ArrayLike, weights: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate and normalize a positive discrete radial window."""

    r = np.asarray(radii, dtype=float)
    w = np.asarray(weights, dtype=float)
    if r.ndim != 1 or w.shape != r.shape or r.size == 0:
        raise RadialLiftInputError("radii and weights must be nonempty matching vectors")
    if np.any(~np.isfinite(r)) or np.any(r <= 0.0):
        raise RadialLiftInputError("all radii must be positive and finite")
    if np.any(~np.isfinite(w)) or np.any(w < 0.0):
        raise RadialLiftInputError("window weights must be finite and nonnegative")
    largest = float(np.max(w))
    if largest <= 0.0:
        raise RadialLiftInputError("window weights must have positive total")
    scaled = w / largest
    return r, scaled / float(np.sum(scaled))


def window_transfer(
    ell: int,
    k: ArrayLike,
    radii: ArrayLike,
    weights: ArrayLike,
) -> NDArray[np.float64]:
    r"""Return ``Psi_l(k)=sum_i w_i j_l(k r_i)`` for a discrete window."""

    ell = _integer_at_least(ell, 0, "ell")
    k_arr = np.asarray(k, dtype=float)
    if k_arr.ndim != 1 or k_arr.size == 0 or np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be a nonempty positive finite vector")
    r, w = normalized_radial_window(radii, weights)
    values = spherical_jn(ell, np.outer(k_arr, r)) @ w
    return np.asarray(values, dtype=float)


def window_powerlaw_cl_quadrature(
    ell: int,
    A_zeta: float,
    theta: float,
    *,
    Z_q: float,
    k_pivot: float,
    k: ArrayLike,
    dlnk_weights: ArrayLike,
    radii: ArrayLike,
    radial_weights: ArrayLike,
) -> float:
    """Finite quadrature of the exact radial-window forward projection."""

    ell = _integer_at_least(ell, 0, "ell")
    A_zeta = _positive(A_zeta, "A_zeta")
    theta = _finite(theta, "theta")
    Z_q = _positive(Z_q, "Z_q")
    k_pivot = _positive(k_pivot, "k_pivot")
    k_arr = np.asarray(k, dtype=float)
    wk = np.asarray(dlnk_weights, dtype=float)
    if k_arr.ndim != 1 or wk.shape != k_arr.shape or k_arr.size == 0:
        raise RadialLiftInputError("k and dlnk_weights must be matching vectors")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be positive and finite")
    if np.any(~np.isfinite(wk)) or np.any(wk <= 0.0):
        raise RadialLiftInputError("dlnk_weights must be positive and finite")
    psi = window_transfer(ell, k_arr, radii, radial_weights)
    delta = A_zeta * (k_arr / k_pivot) ** (-theta)
    value = 4.0 * math.pi * Z_q**2 * float(np.dot(wk, delta * psi**2))
    return _positive(value, "projected C_l")


def finite_window_stability_bound(
    ell: int,
    theta: float,
    *,
    A_zeta: float,
    Z_q: float,
    k_pivot: float,
    R_star: float,
    radii: ArrayLike,
    radial_weights: ArrayLike,
) -> WindowBound:
    r"""Evaluate the analytic finite-window deviation bound for a thin shell.

    In the Hilbert space with norm

    .. math::

        \|f\|_\theta^2=\int_0^\infty d\ln k\,k^{-\theta}|f(k)|^2,

    set ``f_r(k)=j_l(kr)`` and ``f_W=sum_i w_i f_{r_i}``.  The theorem proves

    .. math::

        \|f_W-f_{R_*}\|_\theta\le
        \eta=\frac{2\sqrt{J_\ell(\theta)}}{\theta}
        \sum_iw_i|r_i^{\theta/2}-R_*^{\theta/2}|,

    and then bounds the difference of squared norms.  This retains the Bessel
    ultraviolet decay and is integrable; it replaces the unsafe pointwise Taylor
    bound that can lose that decay after a supremum is taken. The returned
    binary64 values are numerical evaluations, not outward-rounded interval
    certificates; callers must account separately for numerical error.
    """

    ell = _integer_at_least(ell, 1, "ell")
    theta = _finite(theta, "theta")
    if not (0.0 < theta < 2.0 * ell):
        raise RadialLiftInputError("window stability theorem requires 0 < theta < 2*ell")
    A_zeta = _positive(A_zeta, "A_zeta")
    Z_q = _positive(Z_q, "Z_q")
    k_pivot = _positive(k_pivot, "k_pivot")
    R_star = _positive(R_star, "R_star")
    r, w = normalized_radial_window(radii, radial_weights)
    I = mellin_spherical_bessel_square(ell, theta)
    J = derivative_mellin_norm(ell, theta)
    a = 0.5 * theta
    # Subtracting nearly equal powers can produce an exactly zero claimed
    # bound for a nontrivial window.  Retain the radius displacement through
    # log1p/expm1, including the near-scale-invariant theta -> 0+ limit.
    with np.errstate(over="ignore", invalid="ignore"):
        displacement = (r - R_star) / R_star
    near = np.abs(displacement) < 0.5
    log_ratio = np.empty_like(r)
    log_ratio[near] = np.log1p(displacement[near])
    log_ratio[~near] = np.log(r[~near]) - math.log(R_star)
    reference_power = R_star**a
    power_difference = reference_power * np.abs(np.expm1(a * log_ratio))
    eta = (2.0 * math.sqrt(J) / theta) * float(np.dot(w, power_difference))
    shell_norm = reference_power * math.sqrt(I)
    prefactor = 4.0 * math.pi * Z_q**2 * A_zeta * k_pivot**theta
    abs_bound = prefactor * eta * (2.0 * shell_norm + eta)
    shell_cl = prefactor * shell_norm**2
    rel_bound = abs_bound / shell_cl
    if not all(math.isfinite(x) for x in (eta, shell_norm, abs_bound, rel_bound)):
        raise RadialLiftInputError("window bound exceeds the floating-point range")
    return WindowBound(
        ell=ell,
        theta=theta,
        I_ell_theta=I,
        J_ell_theta=J,
        eta=eta,
        shell_norm=shell_norm,
        absolute_cl_bound=abs_bound,
        relative_to_shell_cl_bound=rel_bound,
    )


# ---------------------------------------------------------------------------
# Finite radial operator, null-space and prior-selected continuation
# ---------------------------------------------------------------------------


def radial_projection_matrix(
    ell: ArrayLike,
    k: ArrayLike,
    dlnk_weights: ArrayLike,
    *,
    Z_q: float,
    radii: ArrayLike,
    radial_weights: ArrayLike,
) -> NDArray[np.float64]:
    r"""Build ``A_{ell,j}=4*pi*Z_q^2*w_j*|Psi_l(k_j)|^2``."""

    ell_arr = np.asarray(ell, dtype=int)
    k_arr = np.asarray(k, dtype=float)
    wk = np.asarray(dlnk_weights, dtype=float)
    Z_q = _positive(Z_q, "Z_q")
    if ell_arr.ndim != 1 or ell_arr.size == 0 or np.any(ell_arr < 0):
        raise RadialLiftInputError("ell must be a nonempty vector of nonnegative integers")
    if k_arr.ndim != 1 or wk.shape != k_arr.shape or k_arr.size == 0:
        raise RadialLiftInputError("k and dlnk_weights must be matching vectors")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be positive and finite")
    if np.any(~np.isfinite(wk)) or np.any(wk <= 0.0):
        raise RadialLiftInputError("dlnk_weights must be positive and finite")
    A = np.empty((ell_arr.size, k_arr.size), dtype=float)
    for i, lval in enumerate(ell_arr):
        psi = window_transfer(int(lval), k_arr, radii, radial_weights)
        A[i, :] = 4.0 * math.pi * Z_q**2 * wk * psi**2
    return A


_RADIAL_EPS = np.finfo(float).eps
_RADIAL_RESOLUTION_RTOL = 1e-7


def _radial_array(value: ArrayLike, name: str, ndim: int) -> NDArray[np.float64]:
    """Validate original scalars before a common dtype can discard information."""

    def reject_masks(part: Any) -> None:
        if np.ma.isMaskedArray(part) or part is np.ma.masked:
            raise RadialLiftInputError(f"{name} cannot contain masked data")
        if isinstance(part, (list, tuple)):
            for item in part:
                reject_masks(item)
        elif isinstance(part, np.ndarray) and part.dtype == object:
            for item in part.flat:
                reject_masks(item)

    reject_masks(value)
    try:
        original = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as error:
        raise RadialLiftInputError(f"{name} must be a rectangular real array") from error
    if original.ndim != ndim or not original.size:
        raise RadialLiftInputError(f"{name} must be a nonempty {ndim}-D array")
    result = np.empty(original.shape, dtype=float)
    for index, scalar in np.ndenumerate(original):
        if isinstance(scalar, (bool, np.bool_)) or not isinstance(scalar, (Real, Decimal)):
            raise RadialLiftInputError(f"{name} must contain real, non-Boolean numbers")
        try:
            exact = Fraction(int(scalar)) if isinstance(scalar, Integral) else Fraction(*scalar.as_integer_ratio())
            narrowed = float(scalar)
        except (ValueError, TypeError, OverflowError, AttributeError) as error:
            raise RadialLiftInputError(f"{name} must contain finite float64-representable numbers") from error
        if not math.isfinite(narrowed):
            raise RadialLiftInputError(f"{name} must contain finite numbers")
        converted = Fraction(narrowed)
        if exact.denominator == 1 and converted != exact:
            raise RadialLiftInputError(f"{name} loses integer precision in float64")
        if exact and abs(converted - exact) > abs(exact) * Fraction(4 * _RADIAL_EPS):
            raise RadialLiftInputError(f"{name} loses input precision in float64")
        result[index] = narrowed
    return result


def _radial_rtol(value: float, shape: tuple[int, int]) -> float:
    supplied = float(_radial_array([value], "rtol", 1)[0])
    if not 0 < supplied < 1:
        raise RadialLiftInputError("rtol must satisfy 0 < rtol < 1")
    # A numerical cutoff, not a certificate of algebraic rank or exact zeros.
    return max(supplied, max(shape) * _RADIAL_EPS)


def _radial_float(value: Fraction, name: str) -> float:
    """Round a reported scalar once, refusing overflow and unresolved underflow."""
    try:
        result = float(value)
    except OverflowError as error:
        raise RadialLiftInputError(f"{name} is outside float64 output range") from error
    if not math.isfinite(result) or (value and abs(Fraction(result) - value) > abs(value) * Fraction(4 * _RADIAL_EPS)):
        raise RadialLiftInputError(f"{name} is unresolved in float64 output precision")
    return result


def _radial_fractions(array: NDArray[np.float64]) -> list[list[Fraction]]:
    return [[Fraction(float(x)) for x in row] for row in array]


def _radial_dot(row: Sequence[Fraction], vector: Sequence[Fraction]) -> Fraction:
    return sum((a * b for a, b in zip(row, vector)), Fraction())


def _radial_norm_parts(vector: Sequence[float]) -> tuple[float, float]:
    scale = float(max(map(abs, vector), default=0.0))
    return (scale, math.hypot(*(float(x) / scale for x in vector))) if scale else (0.0, 0.0)


def _radial_norm(vector: Sequence[float], name: str) -> float:
    scale, unit_norm = _radial_norm_parts(vector)
    return _radial_float(Fraction(scale) * Fraction(unit_norm), name)


def _radial_svd(matrix: NDArray[np.float64], rtol: float, *, complete: bool = False):
    scale = float(np.max(np.abs(matrix)))
    normalized = matrix / scale if scale else matrix.copy()
    if scale:
        # A nonzero subnormal quotient can already have lost most of its
        # information. Rescaling the singular value later does not restore it.
        # Ordinary normalized values have standard division precision; check
        # the underflow boundary from the exact supplied binary64 entries.
        for index in zip(*np.nonzero((matrix != 0) & (np.abs(normalized) < np.finfo(float).tiny))):
            normalized[index] = _radial_float(
                Fraction(float(matrix[index])) / Fraction(scale), "operator normalization precision"
            )
    try:
        # A complete right basis only needs the full SVD for a wide matrix.
        u, singular, vh = np.linalg.svd(normalized, full_matrices=complete and matrix.shape[1] > matrix.shape[0])
    except np.linalg.LinAlgError as error:
        raise RadialLiftInputError("radial SVD did not converge") from error
    if not all(np.all(np.isfinite(x)) for x in (u, singular, vh)):
        raise RadialLiftInputError("radial SVD produced nonfinite output")
    threshold = rtol * float(singular[0])
    rank = int(np.count_nonzero(singular > threshold))
    return scale, u, singular, vh, threshold, rank


def radial_null_space_report(matrix: ArrayLike, *, rtol: float = 1e-12) -> NullSpaceReport:
    """Report the raw operator's numerical right kernel from one scaled SVD.

    Singular directions below the relative cutoff are unresolved, not proven
    algebraic zeros. The cutoff is at least ``max(shape)*eps``. At rank zero,
    the condition number on the retained subspace is undefined and is ``None``.
    """

    A = _radial_array(matrix, "matrix", 2)
    cutoff = _radial_rtol(rtol, A.shape)
    scale, _, singular, vh, threshold, rank = _radial_svd(A, cutoff, complete=True)
    return NullSpaceReport(
        shape=(int(A.shape[0]), int(A.shape[1])),
        rank=rank,
        nullity=int(A.shape[1] - rank),
        singular_values=[_radial_float(Fraction(float(x)) * Fraction(scale), "singular value") for x in singular],
        effective_threshold=_radial_float(Fraction(threshold) * Fraction(scale), "rank threshold"),
        condition_number_nonzero=float(singular[0] / singular[rank - 1]) if rank else None,
        null_basis=vh[rank:].tolist(),
        relative_cutoff=cutoff,
    )


def minimum_prior_continuation(
    matrix: ArrayLike,
    screen_cl: ArrayLike,
    *,
    prior_center: ArrayLike,
    prior_precision: ArrayLike,
    rtol: float = 1e-12,
) -> PriorContinuation:
    r"""Numerically minimize the prior under the supplied constraints ``A p=C``.

    It minimizes

    .. math::

        \frac12(p-p_0)^TQ(p-p_0)\quad\text{subject to}\quad Ap=C,

    and therefore returns

    .. math::

        p_*=p_0+Q^{-1}A^T(AQ^{-1}A^T)^+(C-Ap_0).

    The formula is evaluated using an equilibrated Cholesky factor of Q and
    one SVD of the row-equilibrated whitened operator, never an inverse of Q
    or a Gram matrix. Each nonzero row is divided by its maximum magnitude,
    together with the corresponding constraint, to remove equation units.
    The correction, resolution, null projector and effective_rank all use
    that SVD's retained subspace. This numerical rank can differ from the raw-A
    report when the prior changes conditioning. It is not an exact-rank claim.
    Incompatible constraints and unresolved returned quantities are refused.
    """

    A = _radial_array(matrix, "matrix", 2)
    C = _radial_array(screen_cl, "screen_cl", 1)
    p0 = _radial_array(prior_center, "prior_center", 1)
    Q = _radial_array(prior_precision, "prior_precision", 2)
    if C.shape != (A.shape[0],) or p0.shape != (A.shape[1],):
        raise RadialLiftInputError("matrix, screen_cl and prior_center dimensions do not match")
    if Q.shape != (A.shape[1], A.shape[1]):
        raise RadialLiftInputError("prior_precision must be a finite square matrix")
    if not np.array_equal(Q, Q.T):
        raise RadialLiftInputError("prior_precision must be symmetric")
    cutoff = _radial_rtol(rtol, A.shape)
    if np.any(np.diag(Q) <= 0):
        raise RadialLiftInputError("prior_precision must be positive definite")
    # Q = D L L.T D. Diagonal equilibration separates coordinate units from
    # the conditioning of its correlations, and avoids overflowing Q scales.
    diagonal = np.sqrt(np.diag(Q))
    q_exact = _radial_fractions(Q)
    d_exact = [Fraction(float(x)) for x in diagonal]
    correlation = np.array([
        [_radial_float(q_exact[i][j] / (d_exact[i] * d_exact[j]), "prior equilibration")
         for j in range(Q.shape[0])] for i in range(Q.shape[0])
    ])
    try:
        minimum_eigenvalue = float(np.linalg.eigvalsh(correlation)[0])
        # A conservative numerical policy on prior geometry, not an interval
        # eigenvalue certificate. A successful Cholesky alone cannot resolve
        # a tiny correlation eigenvalue, even if the constraints fit exactly.
        prior_norm = float(np.linalg.norm(correlation, ord=np.inf))
        prior_error_scale = Q.shape[0] * _RADIAL_EPS * prior_norm
        if minimum_eigenvalue <= 0 or prior_error_scale > _RADIAL_RESOLUTION_RTOL * minimum_eigenvalue:
            raise RadialLiftInputError("prior_precision correlation geometry is unresolved at float64 precision")
        prior_condition = prior_norm / minimum_eigenvalue
        prior_error_estimate = prior_error_scale / minimum_eigenvalue
        L = np.linalg.cholesky(correlation)
    except np.linalg.LinAlgError as error:
        raise RadialLiftInputError("prior_precision must be numerically positive definite") from error
    a_exact = _radial_fractions(A)
    c_exact = [Fraction(float(x)) for x in C]
    p0_exact = [Fraction(float(x)) for x in p0]
    rhs = [c - _radial_dot(row, p0_exact) for row, c in zip(a_exact, c_exact)]
    # Remove each equation's units before whitening, so independently tiny
    # and huge rows cannot erase one another in a global normalization.
    normalized = A.copy()
    for i, row in enumerate(a_exact):
        row_scale = max(map(abs, row), default=Fraction())
        if row_scale:
            normalized[i] = [_radial_float(x / row_scale, "operator row scaling") for x in row]
            rhs[i] /= row_scale
    equilibrated = np.array([
        [_radial_float(Fraction(float(x)) / d_exact[j], "whitened operator")
         for j, x in enumerate(row)] for row in normalized
    ])
    B = solve_triangular(L, equilibrated.T, lower=True, check_finite=False).T
    if not np.all(np.isfinite(B)):
        raise RadialLiftInputError("whitened operator is outside float64 range")
    row_scales = np.max(np.abs(B), axis=1)
    for i, row_scale in enumerate(row_scales):
        if row_scale:
            B[i] = [_radial_float(Fraction(float(x)) / Fraction(float(row_scale)), "equation equilibration")
                    for x in B[i]]
            rhs[i] /= Fraction(float(row_scale))
    b_scale, u, singular, vh, _, rank = _radial_svd(B, cutoff, complete=True)
    # Separate checks on the prior and inverse miss amplification of whitening
    # roundoff by the inverse. Treat the Cholesky factor as a nearby metric
    # (the prior term above), then propagate row/triangular-solve roundoff
    # through both kappa(L) and the retained kappa(B). Since H=L L.T,
    # sqrt(||H||inf/lambda_min(H)) estimates the first amplification.
    # This first-order conditioning policy is not a certified error bound.
    retained_condition = float(singular[0] / singular[rank - 1]) if rank else 0.0
    inverse_error_estimate = (
        max(A.shape) * _RADIAL_EPS * math.sqrt(prior_condition) * retained_condition
    )
    if prior_error_estimate + inverse_error_estimate > _RADIAL_RESOLUTION_RTOL:
        raise RadialLiftInputError("combined prior and radial inverse geometry is unresolved at float64 precision")
    V = vh[:rank].T
    # Back-transform only retained right singular vectors; no inverse of Q.
    W = solve_triangular(L.T, V, lower=False, check_finite=False)
    W = np.array([
        [_radial_float(Fraction(float(x)) / d_exact[i], "prior back-transform") for x in row]
        for i, row in enumerate(W)
    ]).reshape(A.shape[1], rank)
    rhs_scale = max(map(abs, rhs), default=Fraction())
    tolerance = 128 * max(A.shape) * _RADIAL_EPS
    if rhs_scale:
        rhs_unit = np.array([_radial_float(x / rhs_scale, "constraint scaling") for x in rhs])
        fitted = u[:, :rank] @ (u[:, :rank].T @ rhs_unit)
        if math.hypot(*(rhs_unit - fitted)) > tolerance * math.hypot(*rhs_unit):
            raise RadialLiftInputError("constraints are inconsistent with the retained numerical range")
        if not rank:
            raise RadialLiftInputError("nonzero constraints have no retained operator direction")
        coordinates = (u[:, :rank].T @ rhs_unit) / singular[:rank]
        correction_scale = rhs_scale / Fraction(b_scale)
        coord_exact = [Fraction(float(x)) for x in coordinates]
        correction = [_radial_dot(row, coord_exact) * correction_scale for row in _radial_fractions(W)]
    else:
        correction = [Fraction()] * A.shape[1]
    p = [_radial_float(center + delta, "continued spectrum") for center, delta in zip(p0_exact, correction)]
    p_exact = list(map(Fraction, p))
    residual_exact = [c - _radial_dot(row, p_exact) for row, c in zip(a_exact, c_exact)]
    # Check each supplied equation in its own units; no absolute unit floor.
    for row, target, error in zip(a_exact, c_exact, residual_exact):
        scale = abs(target) + sum((abs(a * x) for a, x in zip(row, p_exact)), Fraction())
        if abs(error) > Fraction(tolerance) * scale:
            raise RadialLiftInputError("returned continuation does not resolve every supplied constraint")
    residual = [_radial_float(x, "constraint residual") for x in residual_exact]
    displacement = [x - center for x, center in zip(p_exact, p0_exact)]
    # The final addition must retain the optimum's correction. A tiny target
    # cannot be declared fitted merely because the prior has large cancelling
    # components, nor may rounding the correction away produce objective zero.
    correction_energy = _radial_dot(correction, [_radial_dot(row, correction) for row in q_exact])
    rounding_error = [actual - intended for actual, intended in zip(displacement, correction)]
    rounding_energy = _radial_dot(rounding_error, [_radial_dot(row, rounding_error) for row in q_exact])
    if correction_energy < 0 or rounding_energy < 0 or rounding_energy > Fraction(_RADIAL_RESOLUTION_RTOL) ** 2 * correction_energy:
        raise RadialLiftInputError("prior correction is unresolved in the returned spectrum precision")
    objective_exact = _radial_dot(displacement, [_radial_dot(row, displacement) for row in q_exact]) / 2
    if objective_exact < 0:
        raise RadialLiftInputError("prior quadratic form is not positive at the returned displacement")
    I = np.eye(A.shape[1])
    if rank == A.shape[1]:
        resolution = I
        null_projector = np.zeros_like(I)
    elif not rank:
        resolution = np.zeros_like(I)
        null_projector = I
    else:
        def projector(basis):
            back = solve_triangular(L.T, basis, lower=False, check_finite=False)
            back = [[Fraction(float(x)) / d_exact[i] for x in row]
                    for i, row in enumerate(back)]
            dual = (basis.T @ L.T) * diagonal
            columns = _radial_fractions(dual.T)
            values = [[_radial_dot(row, column) for column in columns] for row in back]
            terms = [[sum((abs(a * b) for a, b in zip(row, column)), Fraction())
                      for column in columns] for row in back]
            return values, terms

        resolved, resolved_terms = projector(V)
        unresolved, unresolved_terms = projector(vh[rank:].T)
        resolution, null_projector = np.empty_like(I), np.empty_like(I)
        for i in range(A.shape[1]):
            for j in range(A.shape[1]):
                # Evaluate both complementary subspaces. For off-diagonals,
                # choose the sum with less cancellation; on the diagonal keep
                # the smaller component, avoiding subtraction from one there.
                choose_resolved = (abs(resolved[i][j]) <= abs(unresolved[i][j]) if i == j
                                   else resolved_terms[i][j] <= unresolved_terms[i][j])
                if choose_resolved:
                    r, n = resolved[i][j], Fraction(int(i == j)) - resolved[i][j]
                else:
                    n, r = unresolved[i][j], Fraction(int(i == j)) - unresolved[i][j]
                resolution[i, j] = _radial_float(r, "resolution projector")
                null_projector[i, j] = _radial_float(n, "null projector")
    if not np.all(np.isfinite(null_projector)):
        raise RadialLiftInputError("null projector is outside float64 output range")
    # Whitening/truncation must not turn an original, resolved constraint into
    # a claimed null direction. This is checked in every original row's units.
    null_columns = _radial_fractions(null_projector.T)
    for row in a_exact:
        row_scale = sum(map(abs, row), Fraction())
        for column in null_columns:
            column_scale = max(map(abs, column), default=Fraction())
            if abs(_radial_dot(row, column)) > Fraction(tolerance) * row_scale * column_scale:
                raise RadialLiftInputError("numerical null projector loses an original constraint")
    return PriorContinuation(
        p=p,
        residual=residual,
        residual_norm=_radial_norm(residual, "constraint residual norm"),
        objective=_radial_float(objective_exact, "prior objective"),
        resolution=resolution.tolist(),
        null_projector=null_projector.tolist(),
        effective_rank=rank,
        relative_cutoff=cutoff,
    )


def source_powerlaw(k: ArrayLike, A_zeta: float, theta: float, k_pivot: float) -> NDArray[np.float64]:
    """Return the one-dimensional source-derived radial family."""

    k_arr = np.asarray(k, dtype=float)
    if k_arr.ndim != 1 or k_arr.size == 0 or np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be a nonempty positive finite vector")
    return _positive(A_zeta, "A_zeta") * (k_arr / _positive(k_pivot, "k_pivot")) ** (-_finite(theta, "theta"))


def forward_residual(matrix: ArrayLike, spectrum: ArrayLike, screen_cl: ArrayLike) -> Mapping[str, Any]:
    """Return a non-fitting residual, retaining cancellation in supplied data.

    Products and sums use the exact values of the validated binary64 inputs
    before one output rounding. A zero target has relative residual zero only
    for a zero residual; otherwise its relative residual is undefined and this
    function refuses the report. No dimensional denominator floor is used.
    """

    A = _radial_array(matrix, "matrix", 2)
    p = _radial_array(spectrum, "spectrum", 1)
    C = _radial_array(screen_cl, "screen_cl", 1)
    if p.shape != (A.shape[1],) or C.shape != (A.shape[0],):
        raise RadialLiftInputError("matrix, spectrum and screen_cl dimensions do not match")
    p_exact = [Fraction(float(x)) for x in p]
    predicted_exact = [_radial_dot(row, p_exact) for row in _radial_fractions(A)]
    residual_exact = [Fraction(float(c)) - value for c, value in zip(C, predicted_exact)]
    predicted = [_radial_float(x, "forward prediction") for x in predicted_exact]
    residual = [_radial_float(x, "forward residual") for x in residual_exact]
    residual_scale, residual_unit = _radial_norm_parts(residual)
    target_scale, target_unit = _radial_norm_parts(C)
    if not target_scale:
        if residual_scale:
            raise RadialLiftInputError("nonzero residual has undefined relative error against a zero target")
        relative = 0.0
    else:
        relative = _radial_float(
            Fraction(residual_scale) * Fraction(residual_unit) /
            (Fraction(target_scale) * Fraction(target_unit)), "relative forward residual"
        )
    return {
        "predicted": predicted,
        "residual": residual,
        "absolute_l2_residual": _radial_norm(residual, "absolute forward residual"),
        "relative_l2_residual": relative,
    }


# ---------------------------------------------------------------------------
# Physical dilation-intertwiner receipt
# ---------------------------------------------------------------------------


def dilation_intertwiner_receipt(
    k: ArrayLike,
    delta_zeta_sq: ArrayLike,
    theta: float,
    *,
    scale_ratios: Sequence[float],
    tolerance: float,
) -> DilationReceipt:
    r"""Test ``Delta^2(b k)=b^{-theta} Delta^2(k)`` on a finite safe band.

    The theorem-level producer is the physical operator relation

    ``D_s^{-1} C_zeta D_s = exp(-theta*s) C_zeta``

    on the common ``dln(k)`` mode basis.  This finite function tests its diagonal
    consequence after interpolation; it does not substitute for the source DAG
    that proves the operator intertwiner.
    """

    k_arr = np.asarray(k, dtype=float)
    p = np.asarray(delta_zeta_sq, dtype=float)
    theta = _finite(theta, "theta")
    tolerance = _positive(tolerance, "tolerance")
    if k_arr.ndim != 1 or p.shape != k_arr.shape or k_arr.size < 3:
        raise RadialLiftInputError("k and delta_zeta_sq must be matching vectors with at least 3 points")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0) or np.any(np.diff(k_arr) <= 0.0):
        raise RadialLiftInputError("k must be strictly increasing, positive and finite")
    if np.any(~np.isfinite(p)) or np.any(p <= 0.0):
        raise RadialLiftInputError("delta_zeta_sq must be positive and finite")
    logk = np.log(k_arr)
    logp = np.log(p)
    residuals: list[float] = []
    ratios: list[float] = []
    for b_raw in scale_ratios:
        b = _positive(b_raw, "scale ratio")
        if math.isclose(b, 1.0):
            continue
        shifted = logk + math.log(b)
        mask = (shifted >= logk[0]) & (shifted <= logk[-1])
        if not np.any(mask):
            continue
        interpolated = np.interp(shifted[mask], logk, logp)
        res = interpolated - logp[mask] + theta * math.log(b)
        residuals.extend(float(x) for x in res)
        ratios.append(b)
    if not residuals:
        raise RadialLiftInputError("no scale-ratio pair lies inside the safe k band")
    arr = np.asarray(residuals, dtype=float)
    max_res = float(np.max(np.abs(arr)))
    rms = float(np.sqrt(np.mean(arr**2)))
    return DilationReceipt(
        theta=theta,
        scale_ratios=ratios,
        max_absolute_log_residual=max_res,
        rms_log_residual=rms,
        passed=bool(max_res <= tolerance),
        tolerance=tolerance,
        evaluated_pairs=int(arr.size),
    )


def approximate_dilation_shape_bound(
    k: ArrayLike,
    epsilon_log_slope: ArrayLike,
    *,
    k_pivot: float,
) -> NDArray[np.float64]:
    r"""Integrate a bound on ``|d ln Delta^2/d ln k + theta|``.

    Returns the cumulative upper bound

    .. math::

        \left|\ln\frac{\Delta^2(k)}{A_\zeta(k/k_*)^{-\theta}}\right|
        \le\left|\int_{\ln k_*}^{\ln k}\epsilon(e^t)dt\right|.
    """

    k_arr = np.asarray(k, dtype=float)
    eps = np.asarray(epsilon_log_slope, dtype=float)
    k_pivot = _positive(k_pivot, "k_pivot")
    if k_arr.ndim != 1 or eps.shape != k_arr.shape or k_arr.size < 2:
        raise RadialLiftInputError("k and epsilon_log_slope must be matching vectors")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0) or np.any(np.diff(k_arr) <= 0.0):
        raise RadialLiftInputError("k must be strictly increasing, positive and finite")
    if np.any(~np.isfinite(eps)) or np.any(eps < 0.0):
        raise RadialLiftInputError("epsilon_log_slope must be finite and nonnegative")
    logk = np.log(k_arr)
    # Include the pivot by interpolation in the cumulative trapezoid grid.
    if not (k_arr[0] <= k_pivot <= k_arr[-1]):
        raise RadialLiftInputError("k_pivot must lie inside the k grid")
    pivot_log = math.log(k_pivot)
    grid = np.unique(np.concatenate([logk, np.array([pivot_log])]))
    eps_grid = np.interp(grid, logk, eps)
    increments = 0.5 * (eps_grid[:-1] + eps_grid[1:]) * np.diff(grid)
    cumulative = np.concatenate([[0.0], np.cumsum(increments)])
    pivot_index = int(np.where(np.isclose(grid, pivot_log, rtol=0.0, atol=1e-14))[0][0])
    bound_grid = np.abs(cumulative - cumulative[pivot_index])
    return np.interp(logk, grid, bound_grid)


# ---------------------------------------------------------------------------
# Fail-closed receipts
# ---------------------------------------------------------------------------


RADIAL_RECEIPTS = {
    "SCR330_SOURCE_SHELL_EMBEDDING_RECEIPT",
    "SCR330_PHYSICAL_MODE_BASIS_RECEIPT",
    "SCR330_RADIAL_DILATION_INTERTWINER_RECEIPT",
    "SCR330_THIN_SHELL_MELLIN_LIFT_RECEIPT",
    "SCR330_FINITE_WINDOW_KERNEL_RECEIPT",
    "SCR330_RADIAL_NULL_REPORT",
    "SCR330_RADIAL_FORWARD_RESIDUAL_RECEIPT",
    "SCR330_RADIAL_TOMOGRAPHY_RECEIPT",
    "SCR330_RADIAL_PROMOTION_RECEIPT",
    "SCR330_TRANSFER_FIREWALL_RECEIPT",
}


def build_radial_receipt(
    *,
    receipt: str,
    passed: bool,
    claim_tier: str,
    source_dag: Mapping[str, Any],
    blockers: Iterable[str] = (),
    payload: Mapping[str, Any] | None = None,
    physical_tt_te_ee_claim: bool = False,
) -> dict[str, Any]:
    """Build a fail-closed radial-lift receipt."""

    if receipt not in RADIAL_RECEIPTS:
        raise RadialLiftInputError(f"unknown radial receipt: {receipt}")
    if claim_tier not in {"E0", "E1", "E2", "E3", "E4", "E5"}:
        raise RadialLiftInputError("claim_tier must be E0 through E5")
    blocker_set = set(str(x) for x in blockers)
    no_target_ancestor = not _dag_has_blacklisted_ancestor(source_dag)
    if passed and not no_target_ancestor:
        passed = False
        blocker_set.add("measurement_fit_or_likelihood_ancestor")
    if physical_tt_te_ee_claim and claim_tier != "E5":
        passed = False
        blocker_set.add("tt_te_ee_claim_before_E5")
    if receipt == "SCR330_RADIAL_PROMOTION_RECEIPT" and claim_tier != "E4":
        passed = False
        blocker_set.add("radial_primordial_promotion_requires_E4")
    canonical = json.dumps(source_dag, sort_keys=True, separators=(",", ":")).encode("utf-8")
    result: dict[str, Any] = {
        "schema_version": "scr330-radial-v2",
        "receipt": receipt,
        "passed": bool(passed),
        "claim_tier": claim_tier,
        "source_dag_hash": "sha256:" + hashlib.sha256(canonical).hexdigest(),
        "no_measurement_fit_likelihood_ancestor": bool(no_target_ancestor),
        "physical_tt_te_ee_claim": bool(physical_tt_te_ee_claim),
        "blockers": sorted(blocker_set),
    }
    if payload:
        result["payload"] = dict(payload)
    return result


def write_json(path: str | Path, value: Mapping[str, Any] | Sequence[Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def summary_payload(A_q: float, theta: float, *, Z_q: float, R_star: float, k_pivot: float) -> dict[str, Any]:
    """Return a compact source-side summary useful in simulator receipts."""

    amp = primordial_amplitude_from_screen(
        A_q, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    return {"amplitude": asdict(amp), "A_zeta_over_A_q": amp.conversion_factor}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _dag_has_blacklisted_ancestor(dag: Mapping[str, Any]) -> bool:
    nodes = dag.get("nodes")
    if not isinstance(nodes, list):
        return True
    for node in nodes:
        if not isinstance(node, Mapping):
            return True
        if any(
            bool(node.get(key))
            for key in (
                "measurement",
                "fit",
                "likelihood",
                "posterior",
                "planck_shape",
                "target_calibrated_proxy",
                "metadata_unknown",
            )
        ):
            return True
    return False


def _finite(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RadialLiftInputError(f"{name} must be finite")
    return result


def _positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise RadialLiftInputError(f"{name} must be positive")
    return result


def _integer_at_least(value: int, minimum: int, name: str) -> int:
    if isinstance(value, bool):
        raise RadialLiftInputError(f"{name} must be an integer")
    result = int(value)
    if result != value or result < minimum:
        raise RadialLiftInputError(f"{name} must be an integer at least {minimum}")
    return result


if __name__ == "__main__":
    theta = 1.630968209403959 / 48.0
    payload = summary_payload(1.0, theta, Z_q=1.0, R_star=1.0, k_pivot=1.0)
    print(json.dumps(payload, indent=2, sort_keys=True))

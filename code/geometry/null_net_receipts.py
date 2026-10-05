#!/usr/bin/env python3
"""Finite Gaussian null-net diagnostics, with explicit limits on their scope.

The supplied anti-periodic half-filled ring has a Fourier-projector covariance.
Arcs of at most half the ring have faithful reduced states. Site coverage and
finite relative commutant dimensions are algebraic facts, not scaling-limit
standardness.

One smeared bond expectation has an explicit limit and a uniform error bound.
This does not supply mixed-GNS convergence for all observables. Chiral packet
leakage is directional, but the FULL subspace's maximum leakage is equal at
opposite times: finite-dimensional unitarity forbids proper one-sided
subspace inclusion. The resummed commutator envelope is a scalar shape fit,
not Lie closure. All available even ranges are included by default, with
the old R=8 fit retained as a truncation control.

See NULL_NET_AUDIT.md for proofs, independent controls and the old failures.
Running this file refreshes runs/null_net_receipt_report.json (schema 2).
"""
from __future__ import annotations

from fractions import Fraction
import json
from numbers import Real
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import expm

if __package__:
    from .modular_clock_instrumentation import (
        _validate_ring, arc_entanglement_hamiltonian, ring_correlation_value,
    )
else:
    from modular_clock_instrumentation import (
        _validate_ring, arc_entanglement_hamiltonian, ring_correlation_value,
    )

HERE = Path(__file__).resolve().parent
REPORT_PATH = HERE / "runs" / "null_net_receipt_report.json"


def _rings(rings: tuple[int, ...]) -> tuple[int, ...]:
    rings = tuple(rings)
    if not rings:
        raise ValueError("rings must be nonempty")
    for n in rings:
        _validate_ring(n)
    if any(a >= b for a, b in zip(rings, rings[1:])):
        raise ValueError("rings must be strictly increasing")
    return rings


def _matrix(value: np.ndarray, *, hermitian: bool = False) -> np.ndarray:
    value = np.asarray(value)
    if (value.ndim != 2 or value.shape[0] != value.shape[1]
            or value.shape[0] == 0 or value.dtype.kind not in "fciu"
            or not np.all(np.isfinite(value))):
        raise ValueError("expected a finite numeric square matrix")
    try:
        with np.errstate(over="raise", invalid="raise"):
            value = np.asarray(value, dtype=complex if np.iscomplexobj(value) else float)
    except (FloatingPointError, OverflowError) as exc:
        raise ValueError("matrix is outside the binary64 range") from exc
    scale = float(np.max(np.abs(value)))
    if not np.isfinite(scale):
        raise ValueError("matrix magnitude is not representable")
    if hermitian and scale:
        normalized = value / scale
        if np.max(np.abs(normalized - normalized.conj().T)) > 64 * np.finfo(float).eps:
            raise ValueError("matrix must be Hermitian")
    return value


def embed(h_small: np.ndarray, n_total: int, offset: int = 0) -> np.ndarray:
    """Embed a finite square matrix without discarding complex entries."""
    h_small = _matrix(h_small)
    m = h_small.shape[0]
    if (type(n_total) is not int or type(offset) is not int
            or offset < 0 or n_total < offset + m):
        raise ValueError("embedding size and offset must contain the matrix")
    out = np.zeros((n_total, n_total), dtype=h_small.dtype)
    out[offset:offset + m, offset:offset + m] = h_small
    return out


def correlation_matrix(n_ring: int) -> np.ndarray:
    """Full anti-periodic covariance. Wrapping distance loses the seam sign."""
    _validate_ring(n_ring)
    with mp.workdps(120):
        vals = [float(ring_correlation_value(n_ring, d)) for d in range(n_ring)]
    return np.array([[vals[abs(i - j)] for j in range(n_ring)]
                     for i in range(n_ring)])


def nti_receipt(n_ring: int) -> dict:
    """Finite nested matrix-algebra commutant dimension; no limit receipt."""
    _validate_ring(n_ring)
    free_modes = n_ring // 2 - n_ring // 4
    return {"n_ring": n_ring, "relative_commutant_modes": free_modes,
            "relative_commutant_dim_log2": 2 * free_modes,
            "nontrivial": free_modes > 0}


def weak_additivity_receipt(n_ring: int) -> dict:
    """Site coverage on one finite ring, not cofinal weak additivity."""
    _validate_ring(n_ring)
    m = n_ring // 2
    covered, shifts = set(), []
    for s in range(0, n_ring, m // 2):
        covered |= {(s + j) % n_ring for j in range(m)}
        shifts.append(s)
        if len(covered) == n_ring:
            break
    return {"n_ring": n_ring, "translates_used": len(shifts),
            "covers_ring": len(covered) == n_ring}


def separating_modulus_receipt(n_ring: int) -> dict:
    """Numerical gap of a faithful finite arc; no uniform cyclicity bound."""
    _validate_ring(n_ring)
    m = n_ring // 2
    with mp.workdps(120):
        vals = [ring_correlation_value(n_ring, d) for d in range(m)]
        c_arc = mp.matrix([[vals[abs(i - j)] for j in range(m)] for i in range(m)])
        evals, _ = mp.eigsy(c_arc)
        gap = min(min(e, 1 - e) for e in evals)
        if not 0 < gap <= mp.mpf("0.5"):
            raise ValueError("occupation gap unresolved at the working precision")
        log_gap = float(mp.log10(gap))
    return {"n_ring": n_ring, "occupation_gap_log10": log_gap,
            "strictly_faithful": True, "working_precision_digits": 120}


def bond_limit_error_bound(n_ring: int) -> Fraction:
    """Exact conservative bound for the analytic expectation, not float error.

    |E_N - (2/pi) integral_0^1 f| <= (2/3)(3/N+a)/(1-a),
    a=(22/(7N))^2/6, for the fixed Gaussian bump below. See the proof.
    """
    _validate_ring(n_ring)
    a = Fraction(22, 7 * n_ring) ** 2 / 6
    return Fraction(2, 3) * (Fraction(3, n_ring) + a) / (1 - a)


def mixed_gns_cauchy_receipt(rings: tuple[int, ...]) -> dict:
    """Legacy entry point: only ONE fixed bond expectation is controlled.

    E_N = (1/N) sum_{j=0}^{N-2} f(j/N) <c_j^dag c_{j+1}+h.c.>.
    The seam bond is excluded, as in the original observable. The exact
    rational bound applies to the mathematical expectation; reported floats
    are numerical evaluations, not outward-rounded enclosures.
    """
    rings = _rings(rings)
    values = []
    bounds = [bond_limit_error_bound(n) for n in rings]
    with mp.workdps(120):
        sigma, center = mp.mpf(1) / 20, mp.mpf(1) / 4
        integral = sigma * mp.sqrt(mp.pi / 2) * (
            mp.erf((1 - center) / (mp.sqrt(2) * sigma))
            + mp.erf(center / (mp.sqrt(2) * sigma)))
        limit = float(2 * integral / mp.pi)
        for n in rings:
            average = mp.fsum(mp.exp(-((mp.mpf(j) / n - center) / sigma) ** 2 / 2)
                              for j in range(n - 1)) / n
            values.append(float(2 * ring_correlation_value(n, 1) * average))
    diffs = [abs(b - a) for a, b in zip(values, values[1:])]
    return {
        "scope": "one supplied smeared bond expectation; not mixed-GNS convergence",
        "rings": list(rings), "values": values, "limit_value": limit,
        "limit_error_bounds_exact": [str(b) for b in bounds],
        "pairwise_cauchy_bounds_exact": [str(a + b) for a, b in zip(bounds, bounds[1:])],
        "cauchy_differences": diffs,
        "cauchy_decreasing": len(diffs) > 1 and all(a > b for a, b in zip(diffs, diffs[1:])),
        "single_observable_limit_proved": True,
        "mixed_gns_certified": False,
    }


def modular_subspace_diagnostic(h_a: np.ndarray, h_b: np.ndarray,
                                t_mod: float, packet: np.ndarray) -> dict:
    """Compare packet leakage with the maximum over the entire B subspace.

    B is the first len(h_b) coordinates of A. The supplied nonzero B packet
    is normalized. Finite Hermitian matrices and finite positive model time
    are required. Unitary/SVD outputs are numerical diagnostics, never an
    exact inclusion certificate. Internal B evolution preserves B, so it
    cannot change the maximum leakage of exp(i h_A t) on B.
    """
    h_a, h_b = _matrix(h_a, hermitian=True), _matrix(h_b, hermitian=True)
    ma, mb = h_a.shape[0], h_b.shape[0]
    if not 0 < mb < ma:
        raise ValueError("B must be a nonempty proper subspace of A")
    if isinstance(t_mod, (bool, np.bool_)) or not isinstance(t_mod, Real):
        raise ValueError("modular time must be finite and positive")
    t_mod = float(t_mod)
    if not np.isfinite(t_mod) or t_mod <= 0:
        raise ValueError("modular time must be finite and positive")
    packet = np.asarray(packet)
    if (packet.shape != (mb,) or packet.dtype.kind not in "fciu"
            or not np.all(np.isfinite(packet))):
        raise ValueError("packet must be a finite vector in B")
    packet = np.asarray(packet, dtype=complex)
    size = float(np.max(np.abs(packet)))
    if not np.isfinite(size) or size == 0:
        raise ValueError("packet must have a representable nonzero norm")
    packet = packet / size
    packet /= np.linalg.norm(packet)
    h_b = embed(h_b, ma)
    rows = []
    for sign in (1, -1):
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                u = expm(1j * sign * t_mod * h_a) @ expm(-1j * sign * t_mod * h_b)
        except (FloatingPointError, OverflowError, ValueError) as exc:
            raise ValueError("modular flow is numerically unresolved") from exc
        if not np.all(np.isfinite(u)):
            raise ValueError("modular flow is numerically unresolved")
        defect = float(np.linalg.norm(u.conj().T @ u - np.eye(ma), 2))
        if defect > 1e-10:
            raise ValueError("modular flow failed the numerical unitarity check")
        block = u[mb:, :mb]
        rows.append((float(np.linalg.norm(block @ packet) ** 2),
                     float(np.linalg.norm(block, 2) ** 2), defect))
    lp, lm = rows[0][0], rows[1][0]
    lo, hi = sorted((lp, lm))
    # Zero leakage is not floored into a finite ratio. None means undefined.
    ratio = hi / lo if lo > 0 else None
    if ratio is not None and not np.isfinite(ratio):
        ratio = None
    return {
        "modular_time": t_mod, "leakage_plus": lp, "leakage_minus": lm,
        "lower_packet_leakage_sign": "+" if lp < lm else "-" if lm < lp else None,
        "asymmetry_ratio": ratio,
        "max_subspace_leakage_plus": rows[0][1],
        "max_subspace_leakage_minus": rows[1][1],
        "max_subspace_sign_difference": abs(rows[0][1] - rows[1][1]),
        "unitarity_defect_max": max(r[2] for r in rows),
        "proper_finite_subspace_inclusion_possible": False,
        "hsm_inclusion_certified": False,
    }


def hsm_compression_receipt(n_ring: int, t_mod: float = 0.12) -> dict:
    """Legacy name: chiral packet transport, NOT half-sided inclusion."""
    _validate_ring(n_ring)
    ma, mb = n_ring // 2, n_ring // 4
    h_a = arc_entanglement_hamiltonian(n_ring, ma)
    h_b = arc_entanglement_hamiltonian(n_ring, mb)
    xs = np.arange(mb)
    packet = (np.exp(-((xs - mb / 2) ** 2) / (2 * (mb / 6) ** 2))
              * np.exp(1j * np.pi * xs / 2))
    return {"n_ring": n_ring, **modular_subspace_diagnostic(h_a, h_b, t_mod, packet)}


def momentum_profile(comm: np.ndarray, rmax: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Even-range envelope and triangle bound for omitted ranges.

    The bound is for exact finite entries; floating evaluation is not an
    interval certificate. This map has a large kernel: it does not certify
    equality of operators. None sums every fitting even range.
    """
    comm = _matrix(comm)
    n = comm.shape[0]
    if np.iscomplexobj(comm) or n < 3:
        raise ValueError("commutator must be real antisymmetric, size at least three")
    size = float(np.max(np.abs(comm)))
    if size and np.max(np.abs(comm / size + comm.T / size)) > 64 * np.finfo(float).eps:
        raise ValueError("commutator must be real antisymmetric")
    if rmax is not None and (type(rmax) is not int or rmax < 2 or rmax % 2):
        raise ValueError("rmax cutoff must be None or a positive even integer")
    profile, tail = np.zeros(n - 2), np.zeros(n - 2)
    with np.errstate(over="raise", invalid="raise"):
        try:
            for j in range(n - 2):
                for r in range(2, n, 2):
                    i = j + 1 - r // 2
                    if i < 0 or i + r >= n:
                        continue
                    term = (r / 2) * comm[i, i + r]
                    if rmax is None or r <= rmax:
                        profile[j] += (-1) ** ((r - 2) // 2) * term
                    else:
                        tail[j] += abs(term)
        except FloatingPointError as exc:
            raise ValueError("momentum envelope is numerically unresolved") from exc
    return profile, tail


def _shape_fit(measured: np.ndarray, target: np.ndarray) -> dict:
    measured, target = np.asarray(measured), np.asarray(target)
    if (measured.ndim != 1 or measured.size < 2 or target.shape != measured.shape
            or measured.dtype.kind not in "fiu" or target.dtype.kind not in "fiu"
            or not np.all(np.isfinite(measured)) or not np.all(np.isfinite(target))):
        raise ValueError("shape fit requires at least two finite real paired samples")
    try:
        with np.errstate(over="raise", invalid="raise"):
            measured, target = measured.astype(float), target.astype(float)
    except (FloatingPointError, OverflowError) as exc:
        raise ValueError("shape samples are outside the binary64 range") from exc
    sm, st = float(np.max(np.abs(measured))), float(np.max(np.abs(target)))
    if sm == 0 or st == 0:
        raise ValueError("zero signal or target cannot support a shape fit")
    m, g = measured / sm, target / st
    alpha_scaled = float(np.dot(m, g) / np.dot(g, g))
    # Avoid an overflowing/underflowing scale ratio when the final alpha is
    # representable. Replay the RETURNED coefficient in normalized units.
    try:
        alpha = float(Fraction(sm) * Fraction(alpha_scaled) / Fraction(st))
        replay_weight = float(Fraction(alpha) * Fraction(st) / Fraction(sm))
    except OverflowError as exc:
        raise ValueError("shape normalization is numerically unresolved") from exc
    if alpha == 0 and alpha_scaled != 0:
        raise ValueError("shape normalization is numerically unresolved")
    residual = float(np.linalg.norm(m - replay_weight * g) / np.linalg.norm(m))
    if not np.isfinite(alpha) or not np.isfinite(residual):
        raise ValueError("shape fit is numerically unresolved")
    return {"relative_residual": residual, "normalization_alpha": alpha}


def lie_closure_receipt(n_ring: int, rmax: int | None = None) -> dict:
    """Legacy name: scalar momentum-envelope fit, NOT operator Lie closure."""
    _validate_ring(n_ring)
    if n_ring < 32 or n_ring % 8:
        raise ValueError("Lie-envelope rings must be divisible by eight and at least 32")
    if rmax is not None and (type(rmax) is not int or rmax < 2 or rmax % 2):
        raise ValueError("rmax cutoff must be None or a positive even integer")
    m, off_b = n_ring // 2, n_ring // 8
    h = arc_entanglement_hamiltonian(n_ring, m)
    h_a, h_b = embed(h, n_ring), embed(h, n_ring, off_b)
    comm = h_a @ h_b - h_b @ h_a
    measured, tail = momentum_profile(comm, rmax)
    historical, _ = momentum_profile(comm, 8)

    def beta_and_deriv(x, u, v):
        s = lambda y: np.sin(np.pi * y / n_ring)  # noqa: E731
        c = lambda y: np.cos(np.pi * y / n_ring)  # noqa: E731
        b = 2 * n_ring * s(x - u) * s(v - x) / s(v - u)
        db = 2 * np.pi * (c(x - u) * s(v - x) - s(x - u) * c(v - x)) / s(v - u)
        return b, db

    u_a, v_a, u_b, v_b = -0.5, m - 0.5, off_b - 0.5, off_b + m - 0.5
    xs = np.arange(n_ring - 2) + 1.0
    idx = np.where((xs > u_b) & (xs < v_a))[0][2:-2]
    ba, dba = beta_and_deriv(xs[idx], u_a, v_a)
    bb, dbb = beta_and_deriv(xs[idx], u_b, v_b)
    target = ba * dbb - dba * bb
    fitted = _shape_fit(measured[idx], target)
    control = _shape_fit(np.array([comm[j, j + 2] for j in idx]), target)
    old = _shape_fit(historical[idx], target)
    # i[real symmetric, real symmetric] is orthogonal to ALL real symmetric
    # Hermitian matrices, including the previously claimed generator span.
    operator_norm = float(np.linalg.norm(comm))
    if operator_norm == 0 or not np.isfinite(operator_norm):
        raise ValueError("commutator must be resolved and nonzero")
    return {
        "n_ring": n_ring, "scope": "normalized scalar envelope fit only",
        "rmax": rmax, "all_even_ranges": rmax is None,
        "fit_sample_count": len(idx), **fitted,
        "relative_residual_unresummed_control": control["relative_residual"],
        "relative_residual_r8_control": old["relative_residual"],
        "omitted_envelope_bound_max": float(np.max(tail[idx])),
        "commutator_frobenius_norm": operator_norm,
        "real_symmetric_span_relative_distance": float(np.linalg.norm((comm - comm.T) / 2) / operator_norm),
        "operator_lie_closure_certified": False,
    }


def instrument_null_net(rings: tuple[int, ...] = (16, 32, 64)) -> dict:
    rings = _rings(rings)
    if any(n % 8 for n in rings):
        raise ValueError("instrumented rings must be divisible by eight")
    nti = [nti_receipt(n) for n in rings]
    wa = [weak_additivity_receipt(n) for n in rings]
    sep = [separating_modulus_receipt(n) for n in rings]
    cauchy = mixed_gns_cauchy_receipt(rings)
    hsm = [hsm_compression_receipt(n) for n in rings]
    lie = [lie_closure_receipt(n) for n in rings if n >= 32]
    lie_res = [entry["relative_residual"] for entry in lie]
    signs = {x["lower_packet_leakage_sign"] for x in hsm}
    ratios = [x["asymmetry_ratio"] for x in hsm]
    verdicts = {
        "nti_all_stages": all(x["nontrivial"] for x in nti),
        "weak_additivity_all_stages": all(x["covers_ring"] for x in wa),
        "separating_all_stages": all(x["strictly_faithful"] for x in sep),
        "mixed_gns_cauchy_decreasing": cauchy["cauchy_decreasing"],
        "packet_asymmetry_min_ratio": min(ratios) if all(r is not None for r in ratios) else None,
        "packet_direction_consistent": len(signs) == 1 and None not in signs,
        "lie_closure_residuals": lie_res,
        "lie_closure_percent_level": bool(lie_res) and all(r < 0.02 for r in lie_res),
        "lie_closure_residuals_decreasing": len(lie_res) > 1 and all(a > b for a, b in zip(lie_res, lie_res[1:])),
        "lie_closure_rate_certified": False,
    }
    return {
        "artifact": "oph_null_net_receipt_instrumentation", "schema_version": 2,
        "object_id": "NullNetReceipts_Issue503", "issue": 503,
        "audit_issue": 1033,
        "scope": (
            "Supplied finite Gaussian family. Analytic finite faithfulness and "
            "one-observable convergence; numerical packet transport and scalar "
            "envelope fits. No mixed-GNS, HSM, operator Lie-closure, continuum "
            "standardness or physical identification is certified. Legacy verdict "
            "names containing 'lie_closure' refer only to envelope diagnostics."
        ),
        "rings": list(rings), "nti": nti, "weak_additivity": wa,
        "separating": sep, "mixed_gns_cauchy": cauchy,
        "hsm_compression": hsm, "lie_closure": lie, "verdicts": verdicts,
        "receipts_witnessed": {
            "nti": verdicts["nti_all_stages"],
            "weak_additivity": verdicts["weak_additivity_all_stages"],
            "separating_faithfulness": verdicts["separating_all_stages"],
            "single_bond_observable_limit": True,
            "mixed_gns_cauchy": False, "hsm_compression_one_particle": False,
            "modular_lie_closure_percent_level": False,
        },
        "receipts_pending": [
            "Mixed-GNS convergence for the required observable family",
            "Cyc limit clause and standard relative commutants on the common GNS space",
            "Half-sided modular inclusion and modular-intersection relations in the limit",
            "Operator Lie closure and any continuum convergence rate",
            "Source-causal physical placement, faithful embedding, continuum limit and UC/VR/scale identification",
        ],
    }


def main() -> None:
    report = instrument_null_net()
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report["verdicts"], indent=2, allow_nan=False))
    print(json.dumps(report["receipts_witnessed"], indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

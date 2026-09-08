#!/usr/bin/env python3
"""Finite free-fermion modular-profile diagnostics on a supplied collar family.

Half-filled anti-periodic hopping ground states (Fermi velocity one) supply
C_I and h = log((1-C_I)/C_I), evaluated with 120-digit eigendecomposition.
The real h matrix is then rounded to float64. On each half-ring arc, sum
EVERY available odd hopping range centered at each bond:

    beta(j) = -2 sum_r r (-1)^((r-1)/2) h[i,i+r], i=j+(1-r)/2.

This removes finite range truncation, not eigensolver or continuum error.
An optional odd cutoff R retains the historical diagnostic. The triangle
inequality gives |beta_full(j)-beta_R(j)| <=
2 sum_{r>R, fitting at j} r |h[i,i+r]|. This is an exact finite algebraic
bound for exact matrix entries; its reported float evaluation is numerical,
not an interval certificate. No fixed R has a proved vanishing tail here.

Compare with beta_CFT(x)=2N sin(pi(x-u)/N) sin(pi(v-x)/N)/sin(pi(v-u)/N),
using x=j+1/2 and cuts u=-1/2, v=m-1/2, and with the corresponding Mobius
cross-ratio. The long-range terms matter for the continuum reading; see
Eisler, Tonni and Peschel, https://arxiv.org/abs/1902.04474 and, for rings,
Eisler and Peschel, https://arxiv.org/abs/1805.00078. These references do not
certify this finite computation or a uniform error bound for its estimator.

The report preserves the R=9 control: its median profile residual improves
at N=16,32,64 but worsens at N=128. Finite profile agreement and decreasing
cross-ratio differences are diagnostics, not a KMS identity, an analytic
Bisognano-Wichmann identification, or a Cauchy-limit certificate. The
ground-state family, model time and velocity normalization are supplied.
Cap-interior, multiresolution and physical clock attachments remain separate.

Run:
    python3 code/geometry/modular_clock_instrumentation.py
writes code/geometry/runs/modular_clock_instrumentation_report.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
REPORT_PATH = HERE / "runs" / "modular_clock_instrumentation_report.json"

mp.mp.dps = 120
RMAX = 9  # historical fixed-depth control; not the default summation rule
DEFAULT_RINGS = (16, 32, 64, 128)


def _validate_ring(n_ring: int) -> None:
    # Multiples of four avoid a zero-energy occupation ambiguity under the
    # anti-periodic twist.
    if type(n_ring) is not int or n_ring < 8 or n_ring % 4:
        raise ValueError("n_ring must be an integer multiple of four, at least eight")


def _validate_profile_ring(n_ring: int) -> None:
    _validate_ring(n_ring)
    if n_ring % 8:
        raise ValueError("profile rings must be divisible by eight for exact quarter-arc readouts")


# ---------------------------------------------------------------------------
# extended-precision modular data
# ---------------------------------------------------------------------------

def ring_correlation_value(n_ring: int, d: int) -> mp.mpf:
    """C(d) = (1/N) sum_{filled k} cos(k d), k = (2t+1) pi / N (anti-periodic),
    filled = negative-energy modes of -cos k (half filling)."""
    _validate_ring(n_ring)
    if type(d) is not int:
        raise ValueError("d must be an integer")
    s = mp.mpf(0)
    for t in range(n_ring):
        k = (2 * t + 1) * mp.pi / n_ring
        if mp.cos(k) > 0:
            s += mp.cos(k * d)
    return s / n_ring


def arc_entanglement_hamiltonian(n_ring: int, m: int) -> np.ndarray:
    """h = log((1-C_I)/C_I) for the m-site arc, via mpmath eigendecomposition;
    returned as float64. The local context enforces the declared precision
    even if a caller has changed mpmath's process-wide default."""
    _validate_ring(n_ring)
    if type(m) is not int or not 0 < m <= n_ring // 2:
        raise ValueError("arc size must be an integer between one and half the ring")
    with mp.workdps(120):
        vals = {d: ring_correlation_value(n_ring, d) for d in range(m)}
        c_arc = mp.matrix(m, m)
        for i in range(m):
            for j in range(m):
                c_arc[i, j] = vals[abs(i - j)]
        evals, vecs = mp.eigsy(c_arc)
        diag = mp.matrix(m, m)
        for i in range(m):
            if not 0 < evals[i] < 1:
                raise ValueError("occupation spectrum unresolved at the working precision")
            diag[i, i] = mp.log((1 - evals[i]) / evals[i])
        h_arc = vecs * diag * vecs.T
        return np.array([[float(h_arc[i, j]) for j in range(m)] for i in range(m)])


def _profile_and_tail(h_arc: np.ndarray, rmax: int | None) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the finite sum and the triangle-inequality omitted-range bound."""
    h_arc = np.asarray(h_arc)
    if (h_arc.ndim != 2 or h_arc.shape[0] != h_arc.shape[1]
            or h_arc.shape[0] < 2 or not np.isrealobj(h_arc)
            or not np.all(np.isfinite(h_arc))):
        raise ValueError("h_arc must be a finite real square matrix of size at least two")
    if rmax is not None and (type(rmax) is not int or rmax < 1 or rmax % 2 != 1):
        raise ValueError("rmax must be None or a positive odd integer")
    m = h_arc.shape[0]
    beta = np.zeros(m - 1)
    tail = np.zeros(m - 1)
    for j in range(m - 1):
        total = 0.0
        for r in range(1, m, 2):
            i = j + (1 - r) // 2
            if i < 0 or i + r >= m:
                continue
            weight = r * (-1) ** (((r - 1) // 2) % 2)
            if rmax is None or r <= rmax:
                total += -2.0 * weight * h_arc[i, i + r]
            else:
                tail[j] += 2.0 * r * abs(h_arc[i, i + r])
        beta[j] = total
    if not np.all(np.isfinite(beta)) or not np.all(np.isfinite(tail)):
        raise ValueError("profile or omitted-range bound overflowed")
    return beta, tail


def resummed_profile(h_arc: np.ndarray, rmax: int | None = None) -> np.ndarray:
    """Sum all fitting odd ranges by default; an explicit cutoff is diagnostic."""
    return _profile_and_tail(h_arc, rmax)[0]


def truncation_tail_bound(h_arc: np.ndarray, rmax: int | None = None) -> np.ndarray:
    """Numerically evaluate the exact finite omitted-range triangle bound."""
    return _profile_and_tail(h_arc, rmax)[1]


def cft_profile(n_ring: int, x: np.ndarray, u: float, v: float) -> np.ndarray:
    s = lambda y: np.sin(np.pi * y / n_ring)  # noqa: E731
    return 2.0 * n_ring * s(x - u) * s(v - x) / s(v - u)


# ---------------------------------------------------------------------------
# receipts
# ---------------------------------------------------------------------------

def _stage_from_hamiltonian(n_ring: int, h_arc: np.ndarray, rmax: int | None) -> dict:
    m = n_ring // 2
    beta, tail = _profile_and_tail(h_arc, rmax)
    u, v = -0.5, m - 0.5
    xs = np.arange(m - 1) + 0.5
    beta_geo = cft_profile(n_ring, xs, u, v)
    return {"m": m, "beta": beta, "beta_geo": beta_geo, "u": u, "v": v,
            "truncation_tail_bound": tail}


def stage_data(n_ring: int, *, rmax: int | None = None) -> dict:
    _validate_profile_ring(n_ring)
    return _stage_from_hamiltonian(
        n_ring, arc_entanglement_hamiltonian(n_ring, n_ring // 2), rmax
    )


def kms_profile_receipt(data: dict, n_ring: int,
                        interior_fraction: float = 0.5) -> dict:
    beta, beta_geo = data["beta"], data["beta_geo"]
    if not 0 < interior_fraction <= 1 or not np.all(np.isfinite(beta)):
        raise ValueError("invalid interior window or nonfinite profile")
    m1 = len(beta)
    lo = int(m1 * (1 - interior_fraction) / 2)
    hi = m1 - lo
    rel = np.abs(beta[lo:hi] / beta_geo[lo:hi] - 1.0)
    rel_wrong = np.abs(beta[lo:hi] / (1.2 * beta_geo[lo:hi]) - 1.0)
    return {
        "n_ring": n_ring,
        "arc_sites": data["m"],
        "interior_bonds": int(hi - lo),
        "median_relative_residual_2pi": float(np.median(rel)),
        "max_relative_residual_2pi": float(np.max(rel)),
        "median_relative_residual_wrong_1p2": float(np.median(rel_wrong)),
        "separation_factor": float(np.median(rel_wrong) / np.median(rel)),
    }


def crossratio_receipt(data: dict, n_ring: int) -> dict:
    """Transport cross-ratio over the fixed angular quadruple (1/4, 3/4 of
    the arc, anchored at the cuts)."""
    m, beta = data["m"], data["beta"]
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0):
        raise ValueError("cross-ratio transport requires a strictly positive finite profile")
    # exact quarter/three-quarter fractions of the arc: a - u = m/4,
    # b - u = 3m/4, so the Moebius target is stage-independent
    a_idx, b_idx = m // 4 - 1, (3 * m) // 4 - 1
    t = float(np.sum(1.0 / beta[a_idx:b_idx]))
    cr_latt = float(np.exp(2.0 * np.pi * t))
    s = lambda y: np.sin(np.pi * y / n_ring)  # noqa: E731
    u, v = data["u"], data["v"]
    a, b = a_idx + 0.5, b_idx + 0.5
    cr_geo = float((s(b - u) * s(v - a)) / (s(v - b) * s(a - u)))
    return {
        "n_ring": n_ring,
        "quadruple_bonds": [a_idx, b_idx],
        "cr_lattice": cr_latt,
        "cr_moebius": cr_geo,
        "relative_error": abs(cr_latt / cr_geo - 1.0),
    }


def _summarize_stages(stages: dict[int, dict], rings: tuple[int, ...]) -> dict:
    kms = [kms_profile_receipt(stages[n], n) for n in rings]
    crs = [crossratio_receipt(stages[n], n) for n in rings]
    kms_res = [k["median_relative_residual_2pi"] for k in kms]
    cr_err = [c["relative_error"] for c in crs]
    cr_vals = [c["cr_lattice"] for c in crs]
    cauchy = [abs(cr_vals[i + 1] - cr_vals[i]) for i in range(len(cr_vals) - 1)]
    verdicts = {
        "kms_residual_decreasing": all(
            kms_res[i] > kms_res[i + 1] for i in range(len(kms_res) - 1)
        ),
        "kms_final_median_residual": kms_res[-1],
        "wrong_normalization_separated": all(
            k["separation_factor"] > 5.0 for k in kms
        ),
        "crossratio_error_decreasing": all(
            cr_err[i] > cr_err[i + 1] for i in range(len(cr_err) - 1)
        ),
        "crossratio_final_relative_error": cr_err[-1],
        "crossratio_cauchy_decreasing": all(
            cauchy[i] > cauchy[i + 1] for i in range(len(cauchy) - 1)
        ) if len(cauchy) > 1 else True,
    }
    diagnostics = {
        "profile_agreement_with_decreasing_finite_residual": bool(
            verdicts["kms_residual_decreasing"]
            and verdicts["wrong_normalization_separated"]
        ),
        "crossratio_agreement_with_decreasing_finite_differences": bool(
            verdicts["crossratio_error_decreasing"]
            and verdicts["crossratio_cauchy_decreasing"]
        ),
    }
    return {"kms_profile_receipt": kms, "crossratio_receipt": crs,
            "verdicts": verdicts, "finite_diagnostics": diagnostics}


def instrument_tower(rings: tuple[int, ...] = DEFAULT_RINGS) -> dict:
    if len(rings) < 3:
        raise ValueError("at least three increasing cutoffs are required for finite difference trends")
    for n in rings:
        _validate_profile_ring(n)
    if any(a >= b for a, b in zip(rings, rings[1:])):
        raise ValueError("ring cutoffs must be strictly increasing")
    stages, fixed_stages, tails = {}, {}, []
    for n in rings:
        h = arc_entanglement_hamiltonian(n, n // 2)
        stages[n] = _stage_from_hamiltonian(n, h, None)
        fixed_stages[n] = _stage_from_hamiltonian(n, h, RMAX)
        tail = fixed_stages[n]["truncation_tail_bound"]
        tails.append({"n_ring": n, "max_absolute_omitted_range_bound": float(np.max(tail)),
                      "max_relative_omitted_range_bound": float(np.max(tail / stages[n]["beta_geo"])),
                      "max_absolute_full_minus_fixed_profile": float(np.max(abs(stages[n]["beta"] - fixed_stages[n]["beta"])))})
    return {
        "artifact": "oph_modular_clock_instrumentation",
        "object_id": "ModularClockInstrumentation_Issue503",
        "issue": 503,
        "scope": (
            "boundary-collar instrumentation: declared Gaussian MaxEnt "
            "reference states (critical hopping, anti-periodic twist, half "
            "filling, v_F = 1) on the stage-r cap-boundary collar rings; "
            "modular data computed with 120-digit eigendecomposition then "
            "rounded to float64; every available odd hopping range is summed. "
            "Finite profile and transport diagnostics only: no verified KMS "
            "identity, analytic BW identification, continuum convergence, "
            "Cauchy-limit certificate or physical clock attachment."
        ),
        "rings": list(rings),
        "resummation_rule": "all_fitting_odd_ranges_at_each_bond",
        "resummation_rmax_by_ring": [n // 2 - 1 for n in rings],
        **_summarize_stages(stages, rings),
        "fixed_depth_control": {
            "resummation_rmax": RMAX,
            **_summarize_stages(fixed_stages, rings),
            "omitted_range_diagnostics": tails,
            "interpretation": "The fixed-depth finite trend may fail; its omitted tail has no proved vanishing refinement bound.",
        },
        "error_scope": {
            "finite_range_truncation_removed": True,
            "roundoff_interval_certified": False,
            "continuum_profile_error_certified": False,
            "tail_bound": "2 sum over omitted fitting odd r of r abs(h[i,i+r]); triangle inequality for exact entries, numerical evaluation here",
        },
        "receipts_witnessed": {
            "geometric_2pi_kms_boundary_collar": False,
            "modular_cross_ratio_boundary_collar": False,
        },
        "receipts_pending": [
            "uniform continuum profile and transport error bounds",
            "same-algebra KMS identity and BW identification",
            "cap-interior modular data (full 2D cap algebras)",
            "Cyc/NTI/weak-additivity/MI null-net receipts",
            "physical event/link and clock attachment, with the exact embedding "
            "or controlled continuum approximation required by the chosen route",
            "UC/VR/scale physical-identification receipts",
        ],
    }


def main() -> None:
    report = instrument_tower()
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_bytes((json.dumps(report, indent=2, allow_nan=False) + "\n").encode("utf-8"))
    print(json.dumps(report["verdicts"], indent=2))
    print(json.dumps(report["receipts_witnessed"], indent=2))


if __name__ == "__main__":
    main()

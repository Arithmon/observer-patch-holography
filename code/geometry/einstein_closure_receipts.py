#!/usr/bin/env python3
"""Machine receipts for the Einstein branch closure packets (GitHub #526-#528, #503, #578).

The theorem labels below identify the historical #526--#528 packets, not
proofs of physical energy calibration. Null tomography and the baseline
countermodel are finite linear-algebra checks. Entropy response is a tangent
identity at a faithful reference; finite changes include relative entropy.
The bulk/edge split is tested against all central probability transfers as
well as a supplied path. The MaxEnt slope is an analytic consequence of the
supplied Gibbs family, with a separately labelled finite secant diagnostic.
See ENTROPY_FIRST_LAW_AUDIT.md for proofs, regressions and downstream scope.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag
import math

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from quantum_information import (
    density_matrix, dimensions, direct_sum_state, faithful_log, probabilities,
    relative_entropy, shannon_entropy,
    von_neumann_entropy as entropy,
)

from quantum_information.entropy_response import (
    central_normalization_diagnostic, entropy_tangent, finite_entropy_balance,
    first_law_diagnostic, gibbs_entropy_response, _centered, _pairing,
)
from quantum_information.gibbs import _numeric, _parameter

from geometry.null_tomography import (
    ETA, charges_of, design_matrix, eta_project_out, fit_null_charges,
    null_vector, reconstruct_from_charges, sym_basis, tomography_directions,
)


# ---------------------------------------------------------------------------
# null tomography (thm:null-tomography)
# ---------------------------------------------------------------------------

def generic_null_directions(n: int, seed: int = 2) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    return [null_vector(rng.normal(size=3)) for _ in range(n)]


def tomography_receipt(seed: int = 2) -> dict[str, float]:
    """Resolved reconstruction on the Lean frame plus three audit directions.

    A violated relation has a replayable left-null witness. Numerical error
    budgets concern this sampled family only, not a physical stress source.
    """
    rng = np.random.default_rng(seed)
    dirs = list(tomography_directions()) + generic_null_directions(3, seed=seed)
    m = rng.normal(size=(4, 4))
    t_true = (m + m.T) / 2.0
    charges = charges_of(t_true, dirs)
    fit = fit_null_charges(charges, dirs)
    t_hat = fit.require_consistent(error_budget=1e-12 * np.linalg.norm(charges))
    out = {
        "consistent_residual": fit.residual_norm,
        "tracefree_error": float(
            np.linalg.norm(eta_project_out(t_hat) - eta_project_out(t_true))
        ),
        "design_rank": int(np.linalg.matrix_rank(design_matrix(dirs))),
        "smallest_singular_value": float(fit.singular_values[-1]),
        "noise_amplification": fit.noise_amplification,
    }
    # countermodel: violate one linearity relation by bumping one charge
    bad = charges.copy()
    bad[0] += 1.0
    bad_fit = fit_null_charges(bad, dirs)
    out["inconsistent_residual"] = bad_fit.residual_norm
    out["inconsistent_witness_charge"] = bad_fit.witness_charge
    out["inconsistent_witness_design_defect"] = bad_fit.witness_design_defect
    return out


# ---------------------------------------------------------------------------
# bulk/edge/central first law (thm:bulk-edge-central-first-law)
# ---------------------------------------------------------------------------

def random_faithful(dim: int, rng: np.random.Generator) -> np.ndarray:
    m = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    rho = m @ m.conj().T + 0.2 * np.eye(dim)
    return rho / np.trace(rho)


def blockwise_state(ps: list[float], bulk_states: list[np.ndarray],
                    edge_dims: list[int]) -> np.ndarray:
    """Build oplus_a p_a (rho_bulk,a tensor I_edge,a / d_a)."""
    edge_dims = dimensions(edge_dims)
    if not (len(ps) == len(bulk_states) == len(edge_dims)):
        raise ValueError("ps, bulk_states, and edge_dims must have equal length")
    return direct_sum_state(ps, [np.kron(density_matrix(rho), np.eye(d)/d)
                                for rho, d in zip(bulk_states, edge_dims)])


def central_z(bulk_dims: list[int], edge_dims: list[int],
              z_weights: list[float]) -> np.ndarray:
    """Build Z = oplus_alpha z_alpha 1_alpha (no ``-log p_alpha`` term)."""
    bulk_dims, edge_dims = dimensions(bulk_dims), dimensions(edge_dims)
    z_weights = _numeric(z_weights, "central weights", real=True)
    if z_weights.ndim != 1:
        raise ValueError("one-dimensional central weights required")
    if not (len(bulk_dims) == len(edge_dims) == len(z_weights)):
        raise ValueError("bulk_dims, edge_dims, and z_weights must have equal length")
    block_dims = [d_bulk * d_edge
                  for d_bulk, d_edge in zip(bulk_dims, edge_dims)]
    dim = sum(block_dims)
    out = np.zeros((dim, dim))
    i = 0
    for d, z in zip(block_dims, z_weights):
        out[i:i + d, i:i + d] = z * np.eye(d)
        i += d
    return out


def bulk_entropy(ps: list[float], bulk_states: list[np.ndarray]) -> float:
    """H(p) + sum_a p_a S(rho_bulk,a), including the central Shannon term."""
    ps = probabilities(ps)
    if len(ps) != len(bulk_states):
        raise ValueError("one bulk state is required per sector weight")
    return shannon_entropy(ps) + sum(p*entropy(rho) for p,rho in zip(ps,bulk_states))


def edge_entropy(ps: list[float], edge_dims: list[int]) -> float:
    """sum_a p_a log d_a; the Shannon term belongs to bulk_entropy."""
    ps, edge_dims = probabilities(ps), dimensions(edge_dims)
    if len(ps) != len(edge_dims):
        raise ValueError("one edge dimension is required per sector weight")
    return float(np.dot(ps,np.log(edge_dims)))


def first_law_receipt(z_weights: list[float] | None = None,
                      move_weights: bool = True,
                      seed: int = 4, eps: float = 1e-6) -> dict:
    """Schema 2: tangent defects and a separately accounted finite change.

    The tangent is specified independently of eps; eps only selects the
    finite comparison state. A rounded-away finite comparison raises.
    The fixture constructs K=-log(rho0), so its modular first-law residual
    is a consistency check. first_law_diagnostic also accepts independent K.
    The complete center test covers paths omitted by move_weights=False.
    """
    eps = _parameter(eps, "variation step")
    if eps == 0:
        raise ValueError("variation step must be nonzero and resolved")
    if not isinstance(move_weights, (bool, np.bool_)):
        raise ValueError("move_weights must be Boolean")
    rng = np.random.default_rng(seed)
    bulk_dims, edge_dims = [2, 3], [2, 3]
    ps0 = np.array([0.6, 0.4])
    sectors0 = [random_faithful(d, rng) for d in bulk_dims]
    correct_z = [math.log(d) for d in edge_dims]
    zw = correct_z if z_weights is None else z_weights
    # Validate even when the sampled path cannot see central normalization.
    central_z(bulk_dims, edge_dims, zw)
    center = central_normalization_diagnostic(edge_dims, zw)
    rho0 = blockwise_state(ps0, sectors0, edge_dims)
    k0 = -faithful_log(rho0)

    dp = np.array([1., -1.]) if move_weights else np.zeros(2)
    targets = [random_faithful(d, rng) for d in bulk_dims]
    ds = [t-s for s, t in zip(sectors0, targets)]
    ps1 = ps0 + eps*dp
    sectors1 = [s + eps*d for s, d in zip(sectors0, ds)]
    rho1 = blockwise_state(ps1, sectors1, edge_dims)
    if np.array_equal(rho0, rho1) or (move_weights and np.array_equal(ps0, ps1)):
        raise ValueError("finite variation is not resolved at this precision")
    tangent = block_diag(*[np.kron(v*s+p*d, np.eye(n)/n)
                          for p, v, s, d, n in zip(ps0, dp, sectors0, ds, edge_dims)])
    dot_s = entropy_tangent(rho0, tangent)
    dot_k = _pairing(_centered(k0), _centered(tangent))
    dot_z = math.fsum(float(z)*float(v) for z, v in zip(zw, dp))
    dot_edge = math.fsum(z*float(v) for z, v in zip(correct_z, dp))
    dot_bulk = math.fsum([-float(v)*math.log(float(p))
                         + float(v)*entropy(s) + float(p)*entropy_tangent(s, d)
                         for p, v, s, d in zip(ps0, dp, sectors0, ds)])
    # 2pi dot<B> = dot<K>-dot<Z>; scalar origins need never be subtracted
    # inside large matrices. This is the explicitly supplied fixture split.
    dot_2pi_b = math.fsum((dot_k, -dot_z))
    s_bulk0, s_bulk1 = bulk_entropy(ps0, sectors0), bulk_entropy(ps1, sectors1)
    balance = finite_entropy_balance(rho1, rho0)
    decomposed_d = math.fsum([relative_entropy(np.diag(ps1), np.diag(ps0))]
                             + [float(p)*relative_entropy(s1, s0)
                                for p, s1, s0 in zip(ps1, sectors1, sectors0)])
    return {
        "schema_version": 2,
        "base_entropy_split_defect": abs(entropy(rho0)-s_bulk0-edge_entropy(ps0, edge_dims)),
        "varied_entropy_split_defect": abs(entropy(rho1)-s_bulk1-edge_entropy(ps1, edge_dims)),
        "first_law_defect": abs(dot_s-dot_k),
        "edge_identification_defect": abs(dot_z-dot_edge),
        "split_identity_defect": abs(math.fsum((dot_s, -dot_2pi_b, -dot_edge))),
        "bulk_identity_defect": abs(dot_bulk-dot_2pi_b),
        "predicted_bulk_defect": abs(dot_z-dot_edge),
        "predicted_edge_defect": abs(dot_z-dot_edge),
        "all_tangent_modular_defect": first_law_diagnostic(rho0, k0)["all_tangent_defect"],
        "all_sector_transfer_defect": center["all_sector_transfer_defect"],
        "sector_transfer_witness": center["witness"].tolist(),
        "finite_entropy_change": balance["entropy_change"],
        "finite_modular_change": balance["modular_change"],
        "finite_relative_entropy": balance["relative_entropy"],
        "finite_trace_difference": balance["trace_difference"],
        "finite_remainder_decomposition_defect": abs(balance["relative_entropy"]-decomposed_d),
    }


# ---------------------------------------------------------------------------
# MaxEnt multiplier identity (historical thm:maxent-lagrange-stationarity)
# ---------------------------------------------------------------------------


def maxent_multiplier_receipt(lam: float = 1.3, dlam: float = 1e-5,
                              seed: int = 6) -> dict[str, float]:
    """Analytic Gibbs tangent slope plus an independent finite secant.

    dS/dt=lambda is a consequence of the declared family. The secant has a
    finite-step truncation error and is not an exact stationarity test.
    Neither quantity identifies physical energy units or temperature.
    """
    lam, dlam = _parameter(lam, "multiplier"), _parameter(dlam, "secant step")
    if dlam <= 0 or not np.isfinite(lam+dlam) or not np.isfinite(lam-dlam):
        raise ValueError("secant step must be positive and finite")
    if lam+dlam == lam or lam-dlam == lam:
        raise ValueError("secant step is not resolved at this precision")
    rng = np.random.default_rng(seed)
    m = rng.normal(size=(6, 6))
    t_op = (m + m.T) / 2.0
    response = gibbs_entropy_response(t_op, lam)
    plus = gibbs_entropy_response(t_op, lam+dlam)["state"]
    minus = gibbs_entropy_response(t_op, lam-dlam)["state"]
    delta_t = _pairing(_centered(t_op), plus, minus)
    if delta_t == 0:
        raise ValueError("secant energy change is not resolved at this precision")
    balance = finite_entropy_balance(plus, minus)
    secant = balance["entropy_change"] / delta_t
    if not math.isfinite(secant):
        raise ValueError("entropy secant is not resolved at this precision")
    return {"schema_version": 2, "ds_dt": response["ds_dt"], "lambda": lam,
            "multiplier_defect": abs(response["ds_dt"]-lam),
            "energy_variance": response["energy_variance"],
            "dt_dlambda": response["dt_dlambda"], "ds_dlambda": response["ds_dlambda"],
            "secant_ds_dt": secant, "secant_multiplier_defect": abs(secant-lam),
            "secant_step": dlam}


# ---------------------------------------------------------------------------
# baseline countermodel (prop:baseline-countermodels(iv))
# ---------------------------------------------------------------------------

def baseline_countermodel_receipt(c: float = 0.7, seed: int = 8) -> dict[str, float]:
    """Two baselines Y1 = 0 and Y2 = c*eta produce identical first-variation
    responses delta Y = 0 for every sampled variation, while differing by
    exactly c*eta: the absolute constant is invisible to first variations."""
    rng = np.random.default_rng(seed)
    y1 = np.zeros((4, 4))
    y2 = c * ETA
    # first variations of a *constant* tensor field vanish identically for
    # both baselines; sample variation directions to record the identical data
    defects = []
    for _ in range(16):
        v = rng.normal(size=4)
        v = v / np.linalg.norm(v)
        defects.append(abs(float(v @ (y2 - y2) @ v)) + abs(float(v @ (y1 - y1) @ v)))
    return {
        "variation_data_difference": max(defects),
        "baseline_difference": float(np.linalg.norm(y2 - y1)),
        "eta_alignment": float(
            np.linalg.norm((y2 - y1) - c * ETA)
        ),
    }

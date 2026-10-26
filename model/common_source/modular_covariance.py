"""Finite modular dilation bounds and the boundary of compressed readouts.

For Hermitian G and P and real kappa, cyclicity of the finite trace gives
  ||i[G,P]-kappa P||_HS^2 = ||[G,P]||_HS^2+kappa^2||P||_HS^2.
Thus a nonzero P cannot obey a nontrivial boost/translation Ward relation
on the full finite algebra. This is an algebraic identity, not a numerical
rank test. The native example uses G=log(sigma), so the parameter is the
existing dimensionless modular parameter, not native attempts or seconds.

Compression can evade this *full-algebra* bound through a boundary term:
  Q i[G,P] Q = i[QGQ, QPQ]
              + i(QG(1-Q)PQ-QP(1-Q)GQ).
The second term must be retained. A native two-state example matches one
compressed derivative through that term, while failing finite-parameter
dilation. It proves no continuum or physical null-momentum identification.
"""
from __future__ import annotations

import json
import math

import numpy as np

from ..native_modular_reference import local_preparation
from ..observer import native


def _hermitian(value, name):
    raw = np.asarray(value)
    if (np.ma.isMaskedArray(value) or raw.ndim != 2
            or raw.shape[0] != raw.shape[1] or not len(raw)
            or raw.dtype.kind not in "iufc" or not np.isfinite(raw).all()):
        raise ValueError(f"{name} must be a finite nonempty square matrix")
    matrix = raw.astype(complex)
    if not np.allclose(matrix, matrix.conj().T, rtol=0, atol=1e-12):
        raise ValueError(f"{name} must be Hermitian")
    return matrix


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite real scalar")
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _hs_squared(matrix):
    return float(np.vdot(matrix, matrix).real)


def ward_certificate(generator, translation, kappa):
    """Evaluate the exact trace identity on a declared finite matrix pair.

    Translation is a candidate Hermitian observable; positivity is not
    needed for this bound. The values returned here are numerical checks of
    the analytic identity in the module docstring, not interval certificates.
    """
    g = _hermitian(generator, "generator")
    p = _hermitian(translation, "translation")
    if g.shape != p.shape:
        raise ValueError("generator and translation must share one algebra")
    rate = _real(kappa, "dilation rate")
    commutator = g @ p-p @ g
    derivative = 1j*commutator
    defect = derivative-rate*p
    expected = _hs_squared(commutator)+rate**2*_hs_squared(p)
    denominator = abs(rate)*math.sqrt(_hs_squared(p))
    return {
        "dimension": len(g),
        "translation_hs_norm": math.sqrt(_hs_squared(p)),
        "commutator_hs_norm": math.sqrt(_hs_squared(commutator)),
        "defect_hs_norm": math.sqrt(_hs_squared(defect)),
        "squared_identity_residual": abs(_hs_squared(defect)-expected),
        "global_relative_defect_lower_bound": 1.0 if denominator else None,
        "global_relative_defect": math.sqrt(_hs_squared(defect))/denominator
            if denominator else None,
    }


def compressed_certificate(generator, translation, projection, kappa):
    """Check the projected Ward relation with its exact boundary contribution.

    A subspace invariant under G or P has zero boundary contribution and
    inherits the same finite-algebra obstruction. Low-energy matrix elements
    need an explicit domain and boundary budget; their labels alone do not
    provide an escape from it.
    """
    g = _hermitian(generator, "generator")
    p = _hermitian(translation, "translation")
    q = _hermitian(projection, "projection")
    if g.shape != p.shape or q.shape != g.shape:
        raise ValueError("generator, translation and projection must share one algebra")
    if not np.allclose(q @ q, q, rtol=0, atol=1e-12):
        raise ValueError("projection must be an orthogonal projector")
    rate = _real(kappa, "dilation rate")
    outside = np.eye(len(g))-q
    gq, pq = q @ g @ q, q @ p @ q
    internal = 1j*(gq @ pq-pq @ gq)
    boundary = 1j*(q @ g @ outside @ p @ q-q @ p @ outside @ g @ q)
    derivative = q @ (1j*(g @ p-p @ g)) @ q
    residual = derivative-rate*pq
    return {
        "projection_rank": int(round(np.trace(q).real)),
        "projected_translation_hs_norm": math.sqrt(_hs_squared(pq)),
        "internal_derivative_hs_norm": math.sqrt(_hs_squared(internal)),
        "boundary_term_hs_norm": math.sqrt(_hs_squared(boundary)),
        "compressed_defect_hs_norm": math.sqrt(_hs_squared(residual)),
        "decomposition_hs_error": math.sqrt(_hs_squared(derivative-internal-boundary)),
        "necessary_boundary_plus_defect_lower_bound": abs(rate)*math.sqrt(_hs_squared(pq)),
        "generator_leakage_hs_norm": math.sqrt(_hs_squared(outside @ g @ q)),
        "translation_leakage_hs_norm": math.sqrt(_hs_squared(outside @ p @ q)),
    }


def finite_dilation_certificate(generator, translation, parameter, kappa):
    """Full-algebra finite-parameter norm bound for unitary conjugation.

    Unitary conjugation preserves ||P||, hence
    ||U P U* - exp(kappa*s) P|| >= |exp(kappa*s)-1| ||P||.
    This also applies to bounded infinite-dimensional observables. Physical
    continuum translation generators can be unbounded, with their domains
    handled separately.
    """
    g = _hermitian(generator, "generator")
    p = _hermitian(translation, "translation")
    if g.shape != p.shape:
        raise ValueError("generator and translation must share one algebra")
    s, rate = _real(parameter, "modular parameter"), _real(kappa, "dilation rate")
    factor = math.exp(rate*s)
    eigenvalues, vectors = np.linalg.eigh(g)
    unitary = (vectors*np.exp(1j*s*eigenvalues)) @ vectors.conj().T
    p_norm = float(np.linalg.norm(p, 2))
    return {
        "modular_parameter": s,
        "dilation_factor": factor,
        "operator_defect": float(np.linalg.norm(unitary @ p @ unitary.conj().T-factor*p, 2)),
        "operator_defect_lower_bound": abs(factor-1)*p_norm,
    }


def native_certificate():
    """Attach the bounds to the exact same prepared 66-state reference.

    The positive two-state P and coherent probe Q are explicitly declared
    diagnostic observables, not source-selected null momentum or a vacuum.
    No success/failure preparation record is discarded into a claim about
    the continuing native history.
    """
    weights = local_preparation.reference()
    states = native.states(binary=True)
    low = next(i for i, state in enumerate(states)
               if sum(state[p] for p in local_preparation.CORE) == 0)
    high = next(i for i, state in enumerate(states)
                if sum(state[p] for p in local_preparation.CORE) == 2)
    ratio = weights[high]/weights[low]
    if ratio != 5:
        raise ArithmeticError("the source modular gap changed")
    g = np.diag(np.log([float(w) for w in weights]))
    gap = math.log(float(ratio))
    p = np.zeros((len(states), len(states)), complex)
    p[low, low] = p[high, high] = 1
    p[low, high], p[high, low] = .5j, -.5j
    probe = np.zeros(len(states))
    probe[low] = probe[high] = 1/math.sqrt(2)
    q = np.outer(probe, probe)
    rate = gap/2
    compressed = compressed_certificate(g, p, q, rate)
    if compressed["compressed_defect_hs_norm"] > 2e-14:
        raise ArithmeticError("the projected derivative witness failed")
    parameter = 1/gap
    unitary = np.diag(np.exp(1j*parameter*np.diag(g)))
    observed = float(np.vdot(probe, unitary @ p @ unitary.conj().T @ probe).real)
    target = math.exp(.5)
    if abs(observed-(1+.5*math.sin(1))) > 2e-14:
        raise ArithmeticError("native finite-parameter witness changed")
    return {
        "reference": "native local incident preparation; ideal faithful target",
        "dimension": len(states),
        "reference_weights_per_configuration": sorted(set(map(str, weights))),
        "parameter": "dimensionless modular s with U(s)=sigma^(i*s)",
        "gap": "log(5)",
        "declared_dilation_rate": "log(5)/2",
        "candidate_translation_positive_eigenvalues": [.5, 1.5],
        "global_ward": ward_certificate(g, p, rate),
        "compressed_ward": compressed,
        "compressed_finite_parameter": {
            "modular_parameter": parameter,
            "observed_expression": "1+sin(s*log(5))/2",
            "dilation_target_expression": "exp(s*log(5)/2)",
            "observed": observed,
            "dilation_target": target,
            "difference": observed-target,
        },
        "global_finite_parameter": finite_dilation_certificate(g, p, parameter, rate),
        "physical_boost_or_null_momentum_identified": False,
        "continuum_no_go_claimed": False,
    }


if __name__ == "__main__":
    print(json.dumps(native_certificate(), indent=2, sort_keys=True))

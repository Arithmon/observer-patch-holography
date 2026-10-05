"""Tangent entropy response and finite relative-entropy accounting, in nats.

These are diagnostics on supplied finite states and generators, not physical
energy calibration or exact-arithmetic certificates. Matrix logarithms and
spectra retain the support restrictions of the shared state library.
"""

from fractions import Fraction
import math

import numpy as np

from .algebras import resolved_state_spectrum
from .gibbs import _finite, _numeric, _observables, _parameter, _scaled, _thermal
from .states import (
    density_matrix, dimensions, faithful_density_matrix, faithful_log, relative_entropy,
)


def _rounded(value, name):
    """Round one exact binary64 accumulation, refusing lost nonzero results."""
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} exceeds finite numerical range") from exc
    if not math.isfinite(result) or (value and result == 0):
        raise ValueError(f"{name} is not resolved at this precision")
    return result


def _pairing(a, b, reference=None):
    """Re Tr(a* b), or Re Tr(a* (b-reference)), accumulated exactly.

    Only the supplied binary64 components are exact here, not their producer
    (e.g. a matrix logarithm). Subtract before rounding the final scalar.
    """
    result = Fraction(0)
    for x, y, z in zip(a.flat, b.flat,
                       np.zeros_like(b).flat if reference is None else reference.flat):
        for xr, yr, zr in ((x.real, y.real, z.real), (x.imag, y.imag, z.imag)):
            result += Fraction(float(xr)) * (Fraction(float(yr))-Fraction(float(zr)))
    return _rounded(result, "trace pairing")


def _centered(value):
    matrices, scales, _ = _observables([value])
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = matrices[0] * scales[0]
    _finite(result, "centered operator")
    if (np.any((matrices[0].real != 0) & (result.real == 0))
            or np.any((matrices[0].imag != 0) & (result.imag == 0))):
        raise ValueError("centered operator is not resolved at this precision")
    return result


def _norm(value):
    return _finite(math.hypot(*value.real.flat, *value.imag.flat), "operator norm")


def _tangent(value, size):
    raw = _numeric(value, "tangent")
    centered = _centered(raw)  # shared Hermiticity/shape validation
    if raw.shape != (size, size):
        raise ValueError("tangent and reference must use the same algebra")
    trace = sum((Fraction(float(x.real)) for x in np.diag(raw)), Fraction(0))
    if abs(trace) > Fraction(64*np.finfo(float).eps) * Fraction(_norm(raw)):
        raise ValueError("normalized-state tangent must have trace zero")
    # Project only the admitted trace roundoff. No finite state is normalized.
    return centered


def entropy_tangent(sigma, tangent):
    """dS_sigma[D] = -Tr(D log sigma) on normalized Hermitian tangents.

    Trace roundoff up to 64 eps ||D||_HS is explicitly projected out; larger
    violations raise, including small tangents with a large relative trace.
    """
    log_sigma = faithful_log(_numeric(sigma, "reference state"))
    direction = _tangent(tangent, len(log_sigma))
    return -_pairing(_centered(log_sigma), direction)


def first_law_diagnostic(sigma, generator):
    """Test an independently supplied K against *all* normalized tangents.

    R = traceless(K + log sigma); the worst response error over HS-unit
    tangents is ||R||_HS, attained by R/||R||_HS. A scalar K is invisible.
    Zero residual has no nonzero witness. Returned witness defects quantify
    floating-point trace/norm error; no threshold declares a theorem proved.
    """
    log_sigma = faithful_log(_numeric(sigma, "reference state"))
    k = _centered(generator)
    if k.shape != log_sigma.shape:
        raise ValueError("generator and reference must use the same algebra")
    with np.errstate(over="ignore", invalid="ignore"):
        residual = _centered(_finite(k + _centered(log_sigma), "response residual"))
    defect = _norm(residual)
    witness = _scaled(residual, defect) if defect else np.zeros_like(residual)
    return {
        "all_tangent_defect": defect,
        "witness": witness,
        "witness_response": _pairing(residual, witness),
        "witness_trace_defect": abs(math.fsum(float(x.real) for x in np.diag(witness))),
        "witness_norm_defect": abs(_norm(witness) - (1. if defect else 0.)),
    }


def finite_entropy_balance(rho, sigma):
    """Compute S(rho)-S(sigma) = < -log sigma, rho-sigma > - D(rho||sigma).

    This evaluates the identity, rather than independently verifying it.
    It avoids subtraction of two order-one entropies near equilibrium.
    sigma must be faithful; rho may have resolved zero eigenvalues. No
    indefinite source admitted by spectral-roundoff clipping is used.
    The trace difference is reported, not normalized away. D minus that
    difference is the nonnegative Bregman remainder for PSD inputs, even
    when accepted input traces differ within the state validator tolerance.
    """
    a = density_matrix(_numeric(rho, "source state"))
    resolved_state_spectrum(a)
    b = density_matrix(_numeric(sigma, "reference state"))
    log_b = faithful_log(b)
    if a.shape != b.shape:
        raise ValueError("source and reference must use the same algebra")
    linear = -_pairing(log_b, a, b)
    divergence = relative_entropy(a, b)
    trace_difference = _pairing(np.eye(len(a)), a, b)
    remainder = math.fsum((divergence, -trace_difference))
    if not np.array_equal(a, b) and remainder <= 0:
        raise ValueError("finite entropy remainder is not resolved at this precision")
    return {
        "entropy_change": math.fsum((linear, -divergence)),
        "modular_change": linear,
        "relative_entropy": divergence,
        "trace_difference": trace_difference,
        "bregman_remainder": remainder,
    }


def central_normalization_diagnostic(edge_dims, z_weights):
    """Complete center test for z_a = log(d_a) modulo a common constant.

    The spread of z-log(d) is the worst error for transferring one unit of
    probability between sectors (positive and negative masses each one).
    The returned max-to-min transfer witnesses it, even when the caller's
    sampled path holds sector weights fixed. A one-sector center is vacuous.
    """
    dims = dimensions(edge_dims)
    z = _numeric(z_weights, "central weights", real=True)
    if z.shape != (len(dims),):
        raise ValueError("one central weight is required per edge dimension")
    errors = [Fraction(float(v))-Fraction(math.log(d)) for v, d in zip(z, dims)]
    hi, lo = max(range(len(dims)), key=errors.__getitem__), min(range(len(dims)), key=errors.__getitem__)
    witness = np.zeros(len(dims))
    if errors[hi] != errors[lo]:
        witness[hi], witness[lo] = 1., -1.
    return {"all_sector_transfer_defect": _rounded(errors[hi]-errors[lo], "central spread"),
            "witness": witness}


def gibbs_entropy_response(observable, lam):
    """Analytic derivative of a supplied one-constraint Gibbs family.

    dt/dlambda=-Var(T), dS/dlambda=-lambda Var(T). The ratio is an analytic
    consequence of that family, not an independent numerical verification.
    A scalar observable has no energy coordinate and no slope dS/dt.
    Scalar energy origins are removed before the variance is evaluated.
    """
    lam = _parameter(lam, "multiplier")
    matrices, scales, _ = _observables([observable])
    c, scale = matrices[0], float(scales[0])
    # Diagonalize the observable, including at lambda=0. Probabilities must
    # not be recovered by rotating a nearly pure reconstructed matrix back:
    # that would replace small genuine masses by eigensolver roundoff.
    energies, vectors = np.linalg.eigh(c)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        coefficient = _finite(lam*scale, "scaled multiplier")
        ham = _finite(coefficient*energies, "thermal Hamiltonian")
    if lam != 0 and (coefficient == 0
            or np.any((energies != 0) & (ham == 0))):
        raise ValueError("thermal Hamiltonian is not resolved at this precision")
    diagonal, _, _, _ = _thermal(np.diag(ham))
    masses = [Fraction(float(p)) for p in np.diag(diagonal).real]
    rho = faithful_density_matrix((vectors*np.diag(diagonal).real) @ vectors.conj().T)
    # Pairwise variance avoids subtracting two large moments, and also the
    # tiny false variance caused by rounding a dominant eigenvalue's mean.
    variance = sum((masses[i]*masses[j]*(Fraction(float(energies[i]))
                                       - Fraction(float(energies[j])))**2
                    for i in range(len(energies)) for j in range(i)), Fraction(0))
    variance *= Fraction(scale)**2 / sum(masses)**2
    var = _rounded(variance, "energy variance")
    if var <= 0:
        raise ValueError("energy coordinate has no resolved nonzero variance")
    ds = _rounded(-Fraction(lam)*variance, "entropy derivative")
    return {"state": rho, "energy_variance": var, "dt_dlambda": -var,
            "ds_dlambda": ds, "ds_dt": _finite(ds/(-var), "entropy slope")}

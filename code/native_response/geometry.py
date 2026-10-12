"""Small-matrix N4 bare-load Dirichlet response at the tracial reference.

These floating-point functions evaluate the finite formulas. They are not
interval certificates. The all-angle bounds follow from the accompanying
analytic argument, not from sampling these functions.
"""
from __future__ import annotations

from numbers import Real

import numpy as np
from scipy.linalg import expm

from .source import checked_edge, edges, source_generators


def checked_theta(theta: Real) -> float:
    """Accept a finite real scalar, excluding booleans and masked values."""
    if (np.ma.isMaskedArray(theta) or isinstance(theta, (bool, np.bool_))
            or not isinstance(theta, Real)):
        raise ValueError("theta must be a finite real scalar")
    try:
        value = float(theta)
    except (OverflowError, ValueError) as exc:
        raise ValueError("theta must be a finite real scalar") from exc
    if not np.isfinite(value):
        raise ValueError("theta must be a finite real scalar")
    return value


def exchange_amplitude(edge) -> np.ndarray:
    """Return B=H-iG on the two endpoint private factors (dimension 36)."""
    a, b = checked_edge(edge)
    generators = source_generators()
    g = generators[b] - generators[a]
    h = generators[a] + generators[b] - generators.sum(axis=0) / 6
    identity = np.eye(6)
    local = h - 1j * g
    return np.kron(local, identity) + np.kron(identity, local)


def exchange_generator(edge) -> np.ndarray:
    """Return T=[[0,B],[-B dagger,0]] in load-major odd-block ordering."""
    amplitude = exchange_amplitude(edge)
    zero = np.zeros_like(amplitude)
    return np.block([[zero, amplitude], [-amplitude.conj().T, zero]])


def bare_contrast() -> np.ndarray:
    """Return Z=(N_a-N_b) tensor I_36 on |10>, |01>."""
    return np.kron(np.diag([1, -1]), np.eye(36))


def odd_unitary(edge, theta: Real) -> np.ndarray:
    """Evaluate the 72-dimensional odd-block unitary exp(theta T).

    The diagnostic numerical domain is |theta|<=1. The unitarity check is a
    consistency check, not a forward-error certificate for the exponential.
    The exact theorem separately covers every real theta.
    """
    angle = checked_theta(theta)
    if abs(angle) > 1:
        raise ValueError("numerical diagnostics require abs(theta) <= 1")
    generator = exchange_generator(edge)
    with np.errstate(over="ignore", invalid="ignore"):
        argument = angle * generator
    if not np.all(np.isfinite(argument)):
        raise ValueError("theta is outside the numerical exponential range")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        unitary = expm(argument)
    if (not np.all(np.isfinite(unitary)) or
            np.linalg.norm(unitary.conj().T @ unitary - np.eye(72), ord=2) > 1e-10):
        raise ValueError("theta is outside the resolved numerical unitary range")
    return unitary


def load_response(edge, theta: Real) -> dict:
    """Return the odd-block bare response shared by conjugated C and E.

    R=Tr_load(W dagger Z W)/2, output=Phi_e(Z)=W(I_2 tensor R)W dagger,
    alpha=1-Tr(R^2)/36. All traces defining alpha are on private dimension36.
    """
    unitary = odd_unitary(edge, theta)
    pullback = unitary.conj().T @ bare_contrast() @ unitary
    retained = (pullback[:36, :36] + pullback[36:, 36:]) / 2
    output = unitary @ np.kron(np.eye(2), retained) @ unitary.conj().T
    alpha = 1 - float(np.trace(retained @ retained).real / 36)
    return {"pullback": pullback, "retained_private": retained,
            "output": output, "alpha": alpha}


def edge_dirichlet(theta: Real) -> np.ndarray:
    """Compute K_pq=tau[N_p(I-Phi)N_q] by summing all thirty seams.

    tau is normalized trace on binary total-two loads and all twelve private
    factors. The factor1/396 includes the odd-shell probability10/33, the
    conversion from Z to endpoint occupation, and the uniform30-seam mix.
    """
    angle = checked_theta(theta)
    result = np.zeros((12, 12))
    for a, b in edges():
        weight = load_response((a, b), angle)["alpha"] / 396
        result[a, a] += weight
        result[b, b] += weight
        result[a, b] -= weight
        result[b, a] -= weight
    return result


def alpha_lower_bound(theta: Real) -> float:
    """Evaluate max(1/9,1-c theta^4), c=4(49+19sqrt(5))/3.

    This is a floating-point rendering of the proved bound, not a directed
    rounding enclosure. The branch avoids overflow for large finite theta.
    """
    angle = abs(checked_theta(theta))
    coefficient = 4 * (49 + 19 * np.sqrt(5)) / 3
    threshold = ((8 / 9) / coefficient) ** 0.25
    return 1 / 9 if angle >= threshold else 1 - coefficient * angle**4

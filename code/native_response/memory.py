"""Exact locality reduction of the full-network second-moment defect.

The formulas refer to normalized trace on the binary total-two shell and
twelve six-level private factors. Numerical evaluation only uses matrices
of dimension 72, 36, 6 and the thirty-seam Gram matrix. No private factor is
reset between the two events whose correlation is being computed.
"""
from __future__ import annotations

import numpy as np

from .geometry import bare_contrast, checked_theta, load_response
from .source import edges, laplacian


LOAD_COVARIANCE = 5 / 33


def _partial_endpoint(block: np.ndarray, endpoint: int) -> np.ndarray:
    """Normalized partial trace of a 36-dimensional private operator."""
    if np.shape(block) != (36, 36) or endpoint not in (0, 1):
        raise ValueError("expected a 36 by 36 block and endpoint 0 or 1")
    tensor = np.asarray(block).reshape(6, 6, 6, 6)
    return (np.einsum("abcb->ac", tensor) if endpoint == 0
            else np.einsum("abad->bd", tensor)) / 6


def incidence() -> np.ndarray:
    result = np.zeros((30, 12))
    for index, (a, b) in enumerate(edges()):
        result[index, a], result[index, b] = 1, -1
    return result


def closure_data(theta) -> dict:
    """Evaluate the full-network leakage Gram matrix at the supplied theta.

    R_e=Phi_e(Z_e)-beta Z_e, beta=1-alpha. R_e has zero private trace
    in each diagonal load block. Disjoint seams therefore have zero Gram
    entry. For seams sharing p, both must have odd load. The configuration
    with p occupied has nine completions in the other nine ports, while
    the complementary configuration has one. This gives weights 9/66 and
    1/66. Products of off-diagonal load blocks of distinct seams have zero
    load trace; same-seam coherences are included in the full local norm.

    With L=I-Phi and P the orthogonal projection onto centered bare means,
    leakage_matrix is the quadratic form of P L(I-P)L P. Its normalized
    version divides by the bare-mean covariance 5/33. These are also the
    excess second moments for Phi, since I preserves this subspace.
    """
    angle = checked_theta(theta)
    seams = edges()
    reduced, norms, alphas = [], [], []
    contrast = bare_contrast()
    for edge in seams:
        local = load_response(edge, angle)
        alpha = local["alpha"]
        residual = local["output"] - (1 - alpha) * contrast
        alphas.append(alpha)
        norms.append(float(np.vdot(residual, residual).real / 72))
        reduced.append([
            [_partial_endpoint(residual[:36, :36], endpoint),
             _partial_endpoint(residual[36:, 36:], endpoint)]
            for endpoint in (0, 1)
        ])
    gram = np.zeros((30, 30))
    for e, edge in enumerate(seams):
        gram[e, e] = (10 / 33) * norms[e]
        for f in range(e):
            shared = set(edge).intersection(seams[f])
            if not shared:
                continue
            port = next(iter(shared))
            a, b = edge.index(port), seams[f].index(port)
            left, right = reduced[e][a], reduced[f][b]
            value = (9 * np.trace(left[a] @ right[b])
                     + np.trace(left[1 - a] @ right[1 - b])).real / (66 * 6)
            gram[e, f] = gram[f, e] = value
    boundary = incidence()
    leakage = boundary.T @ gram @ boundary / 3600
    return {"theta": angle, "alpha": float(alphas[0]),
            "alphas": np.array(alphas), "edge_gram": gram,
            "leakage_matrix": leakage,
            "normalized_leakage": leakage / LOAD_COVARIANCE}


def band_diagnostics(theta) -> list[dict]:
    """Numerical per-band curvature and its invariance under clock rescaling."""
    data = closure_data(theta)
    values, vectors = np.linalg.eigh(laplacian())
    answer = []
    for value in (5 - np.sqrt(5), 6.0, 5 + np.sqrt(5)):
        basis = vectors[:, np.abs(values - value) < 1e-9]
        excess = np.linalg.eigvalsh(basis.T @ data["normalized_leakage"] @ basis)
        decay = data["alpha"] * value / 60
        answer.append({"laplacian_eigenvalue": float(value),
                       "multiplicity": basis.shape[1],
                       "initial_decay": float(decay),
                       "second_moment_excess": excess.tolist(),
                       "clock_invariant_excess": (excess / decay**2).tolist()})
    return answer

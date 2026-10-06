"""Finite Petz channels and state distances without implicit state repair.

Natural logs are used by the shared entropy API. Fidelity here is squared.
Support inversion has no positive-eigenvalue cutoff: unresolved dense support
raises. Numerical PSD roundoff follows states.py; it is not an exact rank proof.
"""

import math

import numpy as np

from .gibbs import _numeric
from .states import (
    ATOL, _spectrum, _unresolved_positive_spectrum, density_matrix, dimensions,
    partial_trace,
)


def validated_state(value):
    """Strict conversion followed by the shared finite-state validation."""
    return density_matrix(_numeric(value, "state"))


def matrix_sqrt_psd(rho):
    """Square root of a normalized state, retaining positive spectral mass."""
    _, values, vectors = _spectrum(_numeric(rho, "state"))
    return (vectors*np.sqrt(values)) @ vectors.conj().T


def _support_inverse(rho):
    a, values, vectors = _spectrum(_numeric(rho, "reference marginal"))
    kernel = vectors[:, values == 0]
    if (_unresolved_positive_spectrum(a, values)
            or (kernel.size and np.any(a @ kernel != 0))):
        raise ValueError("reference support is numerically unresolved")
    inverse = np.zeros_like(values)
    positive = values > 0
    inverse[positive] = 1/np.sqrt(values[positive])
    return (vectors*inverse) @ vectors.conj().T, kernel


def matrix_inv_sqrt_psd(rho, cutoff=0.):
    """Moore--Penrose inverse root; a positive support cutoff is forbidden."""
    cutoff = _numeric(cutoff, "support cutoff", real=True)
    if cutoff.ndim != 0 or cutoff != 0:
        raise ValueError("positive spectral mass cannot be removed by a cutoff")
    return _support_inverse(rho)[0]


def state_fidelity(rho, sigma):
    """Squared Uhlmann fidelity of supplied states, with no normalization."""
    a, b = validated_state(rho), validated_state(sigma)
    if a.shape != b.shape:
        raise ValueError("fidelity requires states on the same space")
    # Singular values avoid squaring the condition number in sqrt(a) b sqrt(a).
    singular = np.linalg.svd(matrix_sqrt_psd(a) @ matrix_sqrt_psd(b),
                             compute_uv=False)
    return math.fsum(singular)**2


def trace_distance(rho, sigma):
    """Half trace norm of the supplied normalized-state difference."""
    a, b = validated_state(rho), validated_state(sigma)
    if a.shape != b.shape:
        raise ValueError("trace distance requires states on the same space")
    return .5*math.fsum(abs(np.linalg.eigvalsh(a-b)))


def petz_recovery_kraus(reference_bc, dims):
    """Kraus operators for a fixed B -> BC Petz channel.

    On the support of tau_B this is the usual transpose channel of Tr_C.
    The explicit completion sends mass on ker(tau_B) to tau_BC. That choice
    makes the map trace preserving on all B inputs and does not affect an
    AB state with B marginal tau_B. No state-dependent renormalization occurs.
    """
    dims = dimensions(dims)
    if len(dims) != 2:
        raise ValueError("reference dimensions must name B and C")
    b, c = dims
    tau = validated_state(reference_bc)
    if tau.shape != (b*c, b*c):
        raise ValueError("reference and tensor dimensions do not match")
    marginal = partial_trace(tau, dims, [0])
    inverse, kernel = _support_inverse(marginal)
    root = matrix_sqrt_psd(tau)
    kraus = [root @ np.kron(np.eye(b), np.eye(c)[:, j:j+1]) @ inverse
             for j in range(c)]
    kraus.extend(np.outer(root[:, j], kernel[:, k].conj())
                 for k in range(kernel.shape[1]) for j in range(b*c))
    completeness = sum(k.conj().T @ k for k in kraus)
    if (not all(np.all(np.isfinite(k)) for k in kraus)
            or np.linalg.norm(completeness-np.eye(b), ord="fro") > ATOL*b):
        raise ValueError("Petz channel trace preservation is numerically unresolved")
    for k in kraus:
        k.setflags(write=False)
    return tuple(kraus)


def petz_recovery(rho_abc, dims=(2, 2, 2)):
    """Recover ABC from its AB marginal using its fixed BC reference."""
    dims = dimensions(dims)
    if len(dims) != 3:
        raise ValueError("recovery dimensions must name A, B and C")
    rho = validated_state(rho_abc)
    ab = partial_trace(rho, dims, [0, 1])
    bc = partial_trace(rho, dims, [1, 2])
    kraus = petz_recovery_kraus(bc, dims[1:])
    lifted = [np.kron(np.eye(dims[0]), k) for k in kraus]
    recovered = sum(k @ ab @ k.conj().T for k in lifted)
    # Validation cannot manufacture missing trace or change the normalization.
    return validated_state(recovered)

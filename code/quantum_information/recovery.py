"""Finite Petz channels and state distances without implicit state repair.

Natural logs are used by the shared entropy API. Fidelity here is squared.
Support is determined from the original Hermitian entries, without a positive
eigenvalue cutoff. Unresolved numerical functions raise explicitly.
"""

from fractions import Fraction

import numpy as np

from .gibbs import _numeric
from .states import (
    ATOL, _spectrum, _unresolved_positive_spectrum, density_matrix, dimensions,
)
from .state_factors import UnresolvedStatePrecision, source_state, fidelity_factors
from .state_functions import source_root, source_distance, source_reduction, entries_array, channel_entries
from .positive_rounding import _is_exact_psd


class UnresolvedPetzSupport(UnresolvedStatePrecision):
    """The supplied reference cannot be support-inverted at this precision."""


def validated_state(value):
    """Strict conversion followed by the shared finite-state validation."""
    return density_matrix(_numeric(value, "state"))


def matrix_sqrt_psd(rho):
    """Square root of a normalized state, retaining positive spectral mass."""
    return source_root(source_state(rho))


def _support_inverse(rho):
    return _source_inverse(source_state(rho))


def _source_inverse(source):
    a, values, vectors = _spectrum(source.raw)
    kernel = vectors[:, values == 0]
    if (np.count_nonzero(values > 0) != source.rank
            or _unresolved_positive_spectrum(a, values)
            or (kernel.size and np.any(a @ kernel != 0))):
        raise UnresolvedPetzSupport("reference support is numerically unresolved")
    return source_root(source, inverse=True), kernel


def matrix_inv_sqrt_psd(rho, cutoff=0.):
    """Moore--Penrose inverse root; a positive support cutoff is forbidden."""
    cutoff = _numeric(cutoff, "support cutoff", real=True)
    if cutoff.ndim != 0 or cutoff != 0:
        raise ValueError("positive spectral mass cannot be removed by a cutoff")
    return _support_inverse(rho)[0]


def state_fidelity(rho, sigma):
    """Squared Uhlmann fidelity of supplied states, with no normalization."""
    return fidelity_factors(source_state(rho), source_state(sigma))


def trace_distance(rho, sigma):
    """Half trace norm of the supplied normalized-state difference."""
    return source_distance(source_state(rho), source_state(sigma))


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
    return _source_kraus(source_state(reference_bc), dims)


def _source_kraus(tau, dims):
    b, c = dims
    if tau.raw.shape != (b*c, b*c):
        raise ValueError("reference and tensor dimensions do not match")
    marginal = source_reduction(tau, dims, [0])
    inverse, kernel = _source_inverse(marginal)
    root = source_root(tau)
    kraus = [root @ np.kron(np.eye(b), np.eye(c)[:, j:j+1]) @ inverse
             for j in range(c)]
    kraus.extend(np.outer(root[:, j], kernel[:, k].conj())
                 for k in range(kernel.shape[1]) for j in range(b*c))
    completeness = sum(k.conj().T @ k for k in kraus)
    if (not all(np.all(np.isfinite(k)) for k in kraus)
            or np.linalg.norm(completeness-np.eye(b), ord="fro") > ATOL*b):
        raise UnresolvedStatePrecision("Petz channel trace preservation is numerically unresolved")
    # Necessary accuracy check, not a proof of the channel's full action.
    # Correct rounding of tiny root components can still spoil the reference
    # when those binary64 components are multiplied to assemble a channel.
    reproduced = channel_entries(kraus, marginal.entries)
    tolerance = Fraction(ATOL) * b * c
    for i in range(b*c):
        for j in range(b*c):
            for actual, expected in zip(reproduced[i][j], tau.entries[i][j]):
                scale = abs(expected) if expected else tau.trace
                if abs(actual - expected) > tolerance * scale:
                    raise UnresolvedStatePrecision("Petz reference components are numerically unresolved")
    for k in kraus:
        k.setflags(write=False)
    return tuple(kraus)


def petz_recovery(rho_abc, dims=(2, 2, 2)):
    """Recover ABC from its AB marginal using its fixed BC reference."""
    dims = dimensions(dims)
    if len(dims) != 3:
        raise ValueError("recovery dimensions must name A, B and C")
    rho = source_state(rho_abc)
    ab = source_reduction(rho, dims, [0, 1])
    bc = source_reduction(rho, dims, [1, 2])
    kraus = _source_kraus(bc, dims[1:])
    # A PSD state with a rank-one A or C marginal factors across that cut.
    # Its Petz recovery is exactly itself (also for accepted trace roundoff).
    # Use this source certificate instead of manufacturing an indefinite
    # matrix by rounding a known fixed point. The channel-domain checks above
    # still apply; numerical near-purity cannot select this identity.
    if (source_reduction(rho, dims, [0]).rank == 1 or
            source_reduction(rho, dims, [2]).rank == 1):
        recovered = entries_array(rho.entries)
    else:
        lifted = [np.kron(np.eye(dims[0]), k) for k in kraus]
        recovered = entries_array(channel_entries(lifted, ab.entries))
    # A rounded matrix product can be indefinite even when its exact channel
    # expression is PSD. Do not feed that output to a different repaired state
    # distance or add a diagonal floor to manufacture a successful diagnostic.
    result = validated_state(recovered)
    if not _is_exact_psd(result):
        raise UnresolvedStatePrecision("recovered matrix positivity is numerically unresolved")
    return result

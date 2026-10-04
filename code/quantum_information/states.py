"""Validated finite-dimensional density matrices, reductions and entropy.

Natural logarithms and the ordinary matrix trace are used. ATOL bounds
floating-point Hermiticity, normalization and negative-eigenvalue roundoff;
it is not a cutoff on positive eigenvalues or a physical uncertainty.
No positive eigenvalue is discarded or replaced by a logarithm floor.
These numerical diagnostics do not certify exact rank from approximate data.
"""

from math import prod

import numpy as np

ATOL = 1e-12


def dimensions(dims):
    if not isinstance(dims, (list, tuple)) or not dims:
        raise ValueError("nonempty subsystem dimensions required")
    if any(isinstance(d, (bool, np.bool_)) or not isinstance(d, (int, np.integer))
           or d <= 0 for d in dims):
        raise ValueError("subsystem dimensions must be positive integers")
    return tuple(int(d) for d in dims)


def probabilities(values):
    raw = np.asarray(values)
    if raw.ndim != 1 or raw.size == 0 or raw.dtype.kind not in "iuf":
        raise ValueError("nonempty real probability vector required")
    p = raw.astype(float)
    if not np.all(np.isfinite(p)) or np.any(p < 0):
        raise ValueError("probabilities must be finite and nonnegative")
    if abs(float(p.sum()) - 1) > ATOL:
        raise ValueError("probabilities must sum to one; no implicit normalization")
    return p


def _spectrum(value):
    raw = np.asarray(value)
    if (raw.ndim != 2 or raw.shape[0] == 0 or raw.shape[0] != raw.shape[1]
            or raw.dtype.kind not in "iufc"):
        raise ValueError("nonempty square numeric density matrix required")
    a = raw.astype(complex)
    if not np.all(np.isfinite(a)):
        raise ValueError("density matrix must be finite")
    if np.linalg.norm(a - a.conj().T, ord="fro") > ATOL:
        raise ValueError("density matrix must be Hermitian")
    if abs(np.trace(a) - 1) > ATOL:
        raise ValueError("density matrix must have trace one")
    # Symmetrization is restricted to the declared roundoff tolerance.
    a = (a + a.conj().T) / 2
    eigenvalues, vectors = np.linalg.eigh(a)
    if eigenvalues[0] < -ATOL:
        raise ValueError("density matrix must be positive semidefinite")
    return a, np.maximum(eigenvalues, 0), vectors


def density_matrix(value):
    """Validate a normalized state; never normalize an invalid input."""
    return _spectrum(value)[0]


def _indices(indices, count):
    if not isinstance(indices, (list, tuple)):
        raise ValueError("subsystem indices must be a list or tuple")
    if any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer))
           or not 0 <= i < count for i in indices):
        raise ValueError("invalid subsystem index")
    if len(set(indices)) != len(indices):
        raise ValueError("duplicate subsystem index")
    return tuple(sorted(indices))


def partial_trace(rho, dims, keep):
    """Retain subsystems in their original order, including the empty set."""
    dims = dimensions(dims)
    keep = _indices(keep, len(dims))
    a = density_matrix(rho)
    if a.shape[0] != prod(dims):
        raise ValueError("subsystem dimensions do not match density matrix")
    tensor = a.reshape(dims + dims)
    current = len(dims)
    for i in reversed(range(len(dims))):
        if i not in keep:
            tensor = np.trace(tensor, axis1=i, axis2=i + current)
            current -= 1
    size = prod(dims[i] for i in keep)
    return tensor.reshape(size, size)


def von_neumann_entropy(rho):
    _, eigenvalues, _ = _spectrum(rho)
    positive = eigenvalues[eigenvalues > 0]
    return float(-np.sum(positive * np.log(positive)))


def shannon_entropy(values):
    p = probabilities(values)
    positive = p[p > 0]
    return float(-np.sum(positive * np.log(positive)))


def faithful_log(rho):
    """Logarithm on the full algebra; a singular state has no such log."""
    _, eigenvalues, vectors = _spectrum(rho)
    if eigenvalues[0] <= 0:
        raise ValueError("state is not faithful; full matrix logarithm undefined")
    return (vectors * np.log(eigenvalues)) @ vectors.conj().T


def relative_entropy(rho, sigma):
    """Umegaki D(rho || sigma), with +infinity on detected support escape.

Zero eigenvalues of sigma remain zero. A nonzero numerical component of
rho in its computed kernel gives infinity, even if small; there is no
support-leakage tolerance that can silently turn infinity into a finite
answer. Near singularity, roundoff can prevent reliable support decisions;
use exact or higher-precision support data for an exact-support theorem.
"""
    a, eig_a, _ = _spectrum(rho)
    b, eig_b, vec_b = _spectrum(sigma)
    if a.shape != b.shape:
        raise ValueError("relative entropy requires the same algebra")
    if np.array_equal(a, b):
        return 0.0
    kernel = vec_b[:, eig_b == 0]
    if kernel.size and np.any(a @ kernel != 0):
        return float("inf")
    positive_a = eig_a[eig_a > 0]
    positive_b = eig_b > 0
    v = vec_b[:, positive_b]
    diagonal = np.real(np.diag(v.conj().T @ a @ v))
    return float(np.sum(positive_a * np.log(positive_a))
                 - np.dot(diagonal, np.log(eig_b[positive_b])))


def _parts(dims, parts):
    dims = dimensions(dims)
    groups = [_indices(part, len(dims)) for part in parts]
    flat = [i for group in groups for i in group]
    if len(flat) != len(set(flat)):
        raise ValueError("information-theoretic subsystems must be disjoint")
    return dims, groups


def mutual_information(rho, dims, part_x, part_y):
    dims, (x, y) = _parts(dims, (part_x, part_y))
    entropy = lambda keep: von_neumann_entropy(partial_trace(rho, dims, keep))
    return entropy(x) + entropy(y) - entropy(x + y)


def conditional_mutual_information(rho, dims, part_a, part_b, part_c):
    dims, (a, b, c) = _parts(dims, (part_a, part_b, part_c))
    entropy = lambda keep: von_neumann_entropy(partial_trace(rho, dims, keep))
    return entropy(a + b) + entropy(b + c) - entropy(b) - entropy(a + b + c)


def direct_sum_state(weights, states):
    """Retain a classical sector label; this is not an unlabelled mixture."""
    p = probabilities(weights)
    if len(p) != len(states):
        raise ValueError("one normalized state is required per sector weight")
    blocks = [density_matrix(state) for state in states]
    out = np.zeros((sum(len(block) for block in blocks),) * 2, dtype=complex)
    offset = 0
    for weight, block in zip(p, blocks):
        end = offset + len(block)
        out[offset:end, offset:end] = weight * block
        offset = end
    return out


def one_sided_projection(operator, d_left, d_right):
    """HS projection onto M_left tensor 1 + 1 tensor M_right.

    This is a linear operator calculation, not a density-matrix operation.
    A multi-sector collar applies it separately inside each central block.
    """
    left, right = dimensions((d_left,d_right))
    a = np.asarray(operator)
    if (a.shape != (left*right,)*2 or a.dtype.kind not in "iufc"
            or not np.all(np.isfinite(a))):
        raise ValueError("finite operator with matching tensor dimensions required")
    a = a.astype(complex)
    tensor = a.reshape(left,right,left,right)
    l = np.trace(tensor,axis1=1,axis2=3)/right
    r = np.trace(tensor,axis1=0,axis2=2)/left
    return np.kron(l,np.eye(right))+np.kron(np.eye(left),r)-np.trace(a)/(left*right)*np.eye(left*right)

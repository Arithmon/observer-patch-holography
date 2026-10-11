"""The six-level response in SOURCE_RECORD_GLUING.tex, in its source chart.

The vertex order matches source_selection_model/verify_response.py. In this
chart (0, 2) is an edge and (0, 1) is not; integer labels from another carrier
chart must be transported before using these matrices.
"""
from __future__ import annotations

from itertools import combinations
from numbers import Integral

import numpy as np
from scipy.linalg import block_diag


def _frame(phi: float) -> np.ndarray:
    base = np.array([(0, s, t * phi) for s in (1, -1) for t in (1, -1)])
    return np.vstack((base, base[:, [1, 2, 0]], base[:, [2, 0, 1]]))


def vertices() -> np.ndarray:
    """Return the twelve unnormalized source calculation-frame vectors."""
    return _frame((1 + np.sqrt(5)) / 2)


def edges() -> tuple[tuple[int, int], ...]:
    """Infer nearest-neighbor seams from the same source calculation frame."""
    frame = vertices()
    phi = (1 + np.sqrt(5)) / 2
    return tuple((a, b) for a, b in combinations(range(12), 2)
                 if np.isclose(frame[a] @ frame[b], phi, rtol=0, atol=1e-12))


def checked_edge(edge) -> tuple[int, int]:
    """Validate a source edge, retaining the caller's endpoint orientation."""
    if np.ma.isMaskedArray(edge):
        raise ValueError("source edge cannot contain masked entries")
    try:
        endpoints = tuple(edge)
    except TypeError as exc:
        raise ValueError("source edge needs two integer ports") from exc
    if (len(endpoints) != 2 or any(isinstance(p, (bool, np.bool_)) or
                                 not isinstance(p, Integral) for p in endpoints)):
        raise ValueError("source edge needs two integer ports")
    a, b = map(int, endpoints)
    if tuple(sorted((a, b))) not in edges():
        raise ValueError("ports are not a seam in the source chart")
    return a, b


def _hat(vector: np.ndarray) -> np.ndarray:
    x, y, z = vector
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def source_generators() -> np.ndarray:
    """Reconstruct all twelve skew-Hermitian D_p from the source formula.

    Galois conjugation replaces sqrt(5) by -sqrt(5) in the real frame; it
    does not conjugate the complex unit in the upper quadrupole block.
    """
    frame = vertices()
    conjugate = _frame((1 - np.sqrt(5)) / 2)
    identity = np.eye(3)
    return np.array([
        block_diag(_hat(v / 2) + 1j * (np.outer(v, v) / 2
                                      - (v @ v / 6 - 1 / 12) * identity),
                   _hat(w / 2))
        for v, w in zip(frame, conjugate)
    ])


def laplacian() -> np.ndarray:
    """Return the integer carrier Laplacian in this same source chart."""
    result = np.zeros((12, 12), dtype=int)
    for a, b in edges():
        result[a, a] += 1
        result[b, b] += 1
        result[a, b] -= 1
        result[b, a] -= 1
    return result

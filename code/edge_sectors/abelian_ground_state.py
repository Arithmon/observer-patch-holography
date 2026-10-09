"""Resolve the existing 2x2 Z_n ground state in its electric symmetry sector.

For h > 0 the original link-basis Hamiltonian is irreducible with nonpositive
off-diagonal entries. Its unique positive ground state is fixed by the four
star permutations and the two noncontractible electric-cut permutations,
which commute with H. In the electric basis this means zero divergence and
zero winding. These fields are exactly plaquette boundaries over Z_n, also
when n is composite. Four plaquette potentials modulo their common constant
therefore give n**3 states instead of n**8.

Fix c00 = 0 and write (a,b,c) = (c01,c10,c11). In the original link order the
electric field is (-a,b,a,c-a,b-c,-b,c-b,a-c), modulo n. The restricted star
at (0,0) measures q = b-c. The magnetic term changes one potential by one;
changing c00 and restoring the chosen representative changes all three by
minus one. The Gauss term is the constant -20 in this sector.

The reduction is exact algebra. The eigensolve remains floating point: we
require a positive vector and check every component of its eigen-equation
against positive neighboring amplitudes before summing squared amplitudes.
These checks detect lost tails; they are not interval error bounds. A failed
check is an unresolved calculation, not a zero probability. At h = 0 the
original full model has degenerate topological sectors, so this API refuses
to choose one on the caller's behalf.
"""
from __future__ import annotations

import math
from fractions import Fraction
from numbers import Integral, Rational, Real

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import ArpackNoConvergence, eigsh, spsolve


_COMPONENT_RTOL = 1e-8


def _reduced_operator(n: int, h: float):
    """Return a scaled H + 20 + 8h, its positive hopping, and edge charges."""
    a, b, c = np.indices((n, n, n))
    scale = max(1.0, h)
    diagonal = 4.0 * (h / scale) * (
        np.sin(np.pi * a / n) ** 2 + np.sin(np.pi * b / n) ** 2
        + np.sin(np.pi * (b - c) / n) ** 2
        + np.sin(np.pi * (c - a) / n) ** 2
    )
    indices = np.arange(n**3).reshape((n, n, n))
    neighbors = []
    for sign in (-1, 1):
        neighbors.extend(np.roll(indices, sign, axis=axis).ravel() for axis in range(3))
        neighbors.append(np.roll(indices, (sign,) * 3, axis=(0, 1, 2)).ravel())
    # Duplicate translations at n=2 must add, since Re(B) then equals B.
    hopping = coo_matrix(
        (np.full(8 * n**3, 0.5 / scale),
         (np.tile(indices.ravel(), 8), np.concatenate(neighbors))),
        shape=(n**3, n**3),
    ).tocsr()
    diagonal = diagonal.ravel()
    return diags(diagonal) - hopping, diagonal, hopping, ((b - c) % n).ravel()


def _schur_ground(diagonal, hopping):
    """Positive tail solve when the nonzero-electric block is dominant.

    Write the shifted Hamiltonian as [[0,-b.T],[-b,A]]. For ground energy
    -s and origin amplitude 1, its other amplitudes are z=(A+sI)^-1 b and
    s=b.T z. Strict diagonal dominance makes A a positive M-matrix, so its
    inverse is nonnegative. The scalar right side decreases with s: the
    ground root lies between zero and b.T A^-1 b. Double the upper endpoint
    to separate its residual from zero, then scale the bracket to [0,1],
    retaining small roots without an absolute root-finder floor.
    """
    block = diags(diagonal[1:]) - hopping[1:, 1:]
    source = hopping[1:, 0].toarray().ravel()
    upper = 2.0 * float(source @ spsolve(block, source))
    if not math.isfinite(upper) or upper <= 0:
        raise ValueError("Z_n ground-state precision unresolved: Schur root underflow")

    def tails(ratio):
        return spsolve(block + diags(np.full(len(source), upper * ratio)), source)

    ratio = brentq(lambda x: x - float(source @ tails(x)) / upper,
                   0.0, 1.0, xtol=1e-300, rtol=8 * np.finfo(float).eps)
    return np.concatenate(([1.0], tails(ratio)))


def _sector_probabilities(ground, charges, n):
    """Sum exact squares of supplied binary64 amplitudes before conversion.

    A basis weight may underflow even when its sector's total is resolved.
    Fraction arithmetic preserves those weights through summation and
    normalization. The single final conversion must preserve each positive
    probability with at most 1e-12 relative error, including subnormals.
    This controls conversion only, not the preceding ground-state solve.
    """
    weights = [Fraction(0) for _ in range(n)]
    for charge, amplitude in zip(charges, ground, strict=True):
        exact_amplitude = Fraction(float(amplitude))
        weights[int(charge)] += exact_amplitude**2
    total = sum(weights)
    probabilities = []
    for weight in weights:
        exact_probability = weight / total
        probability = float(exact_probability)
        if probability <= 0:
            raise ValueError("Z_n ground-state precision unresolved: positive sector probability underflow")
        if abs(Fraction(probability) - exact_probability) > exact_probability / 10**12:
            raise ValueError("Z_n ground-state precision unresolved: sector probability conversion")
        probabilities.append(probability)
    return np.array(probabilities)


def zn_edge_distribution(n: int, h: float) -> np.ndarray:
    """Return all n electric-center probabilities, or refuse unresolved precision.

    n is an integer >= 2 and h is a finite, strictly positive real coupling
    exactly representable in binary64; conversion must not change its value.
    All probabilities of this finite connected model are strictly positive.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n < 2:
        raise ValueError("n must be an integer >= 2")
    if isinstance(h, (bool, np.bool_)) or not isinstance(h, Real):
        raise ValueError("h must be a finite, strictly positive real coupling")
    original_h = h
    try:
        h = float(original_h)
    except (OverflowError, ValueError) as exc:
        raise ValueError("h must be a finite, strictly positive real coupling") from exc
    if not math.isfinite(h):
        raise ValueError("h must be a finite, strictly positive real coupling")
    if isinstance(original_h, Rational):
        numerator, denominator = int(original_h.numerator), int(original_h.denominator)
    elif hasattr(original_h, "as_integer_ratio"):
        numerator, denominator = original_h.as_integer_ratio()
    else:
        raise ValueError("h must supply an exact value representable in binary64")
    float_numerator, float_denominator = h.as_integer_ratio()
    if numerator * float_denominator != float_numerator * denominator:
        raise ValueError("h must be exactly representable in binary64 without changing its value")
    if h <= 0:
        raise ValueError("h must be finite and strictly positive; h=0 has degenerate topological sectors")
    n = int(n)
    operator, diagonal, hopping, charges = _reduced_operator(n, h)
    if np.any(diagonal[1:] == 0):
        raise ValueError("Z_n ground-state precision unresolved: nonzero electric diagonal underflow")
    dim = n**3
    if np.min(diagonal[1:]) > np.max(np.asarray(hopping.sum(axis=1))):
        # Eigensolvers control absolute vector error, which can overwhelm a
        # representable rare sector. Here positive block solves resolve tails.
        ground = _schur_ground(diagonal, hopping)
    else:
        try:
            _, vectors = eigsh(operator, k=1, which="SA", v0=np.ones(dim),
                               maxiter=20000, tol=1e-13)
        except ArpackNoConvergence as exc:
            raise ValueError("Z_n ground-state precision unresolved: eigensolve did not converge") from exc
        ground = np.asarray(vectors[:, 0], dtype=float)
    if not np.isfinite(ground).all() or ground[0] == 0:
        raise ValueError("Z_n ground-state precision unresolved: missing positive ground vector")
    if ground[0] < 0:
        ground = -ground
    if np.any(ground <= 0):
        raise ValueError("Z_n ground-state precision unresolved: nonpositive amplitude")
    # The zero-electric row has zero diagonal. Recover -E there without
    # subtracting large electric and magnetic contributions to a Rayleigh sum.
    incoming = hopping @ ground
    minus_energy = incoming[0] / ground[0]
    left = (diagonal + minus_energy) * ground
    if (not math.isfinite(minus_energy) or minus_energy <= 0
            or not np.isfinite(left).all() or np.any(incoming <= 0)):
        raise ValueError("Z_n ground-state precision unresolved: vanished positive eigen-equation")
    relative_residual = np.abs(left - incoming) / np.maximum(left, incoming)
    if not np.isfinite(relative_residual).all() or np.max(relative_residual) > _COMPONENT_RTOL:
        raise ValueError("Z_n ground-state precision unresolved: componentwise eigen-equation residual")

    return _sector_probabilities(ground, charges, n)

"""Finite-precision Majorana readout without squaring the condition number.

Return ascending singular masses and a unitary U with U.T @ M @ U = diag(m)
to a normwise relative tolerance. This is a numerical residual check, not an
interval certificate, relative error bound for tiny masses, or choice of a
physical basis inside a degenerate eigenspace. A degenerate SVD block that is
not already diagonal under congruence requires a separate block solver and is
explicitly refused.
"""

import numpy as np

from quantum_information.gibbs import _numeric


def validate_takagi(matrix, masses, unitary):
    """Check a supplied readout in scaled units, including the real shortcut."""
    matrix = _numeric(matrix, "Majorana matrix")
    masses = _numeric(masses, "Takagi masses", real=True)
    unitary = _numeric(unitary, "Takagi columns")
    if (matrix.ndim != 2 or not matrix.shape[0] or matrix.shape[0] != matrix.shape[1]
            or unitary.shape != matrix.shape or masses.shape != (len(matrix),)):
        raise ValueError("Takagi readout shapes disagree")
    if np.any(masses < 0):
        raise ValueError("Takagi masses must be nonnegative")
    scale = max(float(np.max(np.abs(matrix.real))), float(np.max(np.abs(matrix.imag))))
    if scale == 0:
        scaled, relative_masses = matrix, masses
    else:
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            scaled = matrix.real / scale + 1j * (matrix.imag / scale)
            relative_masses = masses / scale
        if not np.all(np.isfinite(relative_masses)):
            raise ValueError("Takagi masses exceed their matrix scale")
        if (np.any((matrix.real != 0) & (scaled.real == 0))
                or np.any((matrix.imag != 0) & (scaled.imag == 0))):
            raise ValueError("Majorana normalization loses nonzero components")
    if np.max(np.abs(scaled - scaled.T)) > 1e-12:
        raise ValueError("Majorana matrix must be complex symmetric")
    tolerance = 1e-10 * float(np.max(np.abs(scaled)))
    if np.max(np.abs(unitary.T @ scaled @ unitary - np.diag(relative_masses))) > tolerance:
        raise ValueError("Takagi masses disagree with the positive congruence diagonal")
    if np.max(np.abs(unitary.conj().T @ unitary - np.eye(len(matrix)))) > 1e-12:
        raise ValueError("Takagi columns are not unitary")


def sorted_takagi(matrix):
    matrix = _numeric(matrix, "Majorana matrix")
    if matrix.ndim != 2 or not matrix.shape[0] or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Majorana matrix must be nonempty and square")
    scale = max(float(np.max(np.abs(matrix.real))), float(np.max(np.abs(matrix.imag))))
    if scale == 0:
        return np.zeros(len(matrix)), np.eye(len(matrix), dtype=complex)
    # Divide components separately: complex division can overflow its reciprocal
    # for small units. Refuse a scale change that discards supplied information.
    with np.errstate(under="ignore"):
        scaled = matrix.real / scale + 1j * (matrix.imag / scale)
    if (np.any((matrix.real != 0) & (scaled.real == 0))
            or np.any((matrix.imag != 0) & (scaled.imag == 0))):
        raise ValueError("Majorana normalization loses nonzero components")
    if np.max(np.abs(scaled - scaled.T)) > 1e-12:
        raise ValueError("Majorana matrix must be complex symmetric")

    # Diagonalizing M†M loses singular values below sqrt(epsilon) times the
    # largest mass and can underflow even when every mass is representable.
    # Right singular vectors directly supply the Takagi columns, up to phases
    # (and an additional congruence within unresolved degenerate blocks).
    _, singular, vh = np.linalg.svd(scaled)
    order = np.argsort(singular)
    singular = singular[order]
    unitary = vh.conj().T[:, order]
    congruence = unitary.T @ scaled @ unitary
    tolerance = 1e-10 * float(np.max(singular))
    offdiag = congruence - np.diag(np.diag(congruence))
    if np.max(np.abs(offdiag)) > tolerance:
        raise ValueError("Takagi eigenspaces require a degenerate-block congruence resolution")
    unitary = unitary @ np.diag(np.exp(-0.5j * np.angle(np.diag(congruence))))
    diagonalized = unitary.T @ scaled @ unitary
    if np.max(np.abs(diagonalized - np.diag(singular))) > tolerance:
        raise ValueError("Takagi masses disagree with the positive congruence diagonal")
    if np.max(np.abs(unitary.conj().T @ unitary - np.eye(len(matrix)))) > 1e-12:
        raise ValueError("Takagi columns are not unitary")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        masses = singular * scale
    if not np.all(np.isfinite(masses)) or np.any((singular != 0) & (masses == 0)):
        raise ValueError("Takagi masses exceed finite binary64 range")
    return masses, unitary

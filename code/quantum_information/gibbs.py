"""Finite Gibbs numerics shared by MaxEnt inference and direct-sum collars.

Scalar energy origins are removed before diagonalization. Relative sector
energies remain physical and are retained in the global partition function.
These are floating-point states with checked support, not exact certificates.
"""

from fractions import Fraction
import math

import numpy as np

from .states import faithful_density_matrix, finite_real_scalar


def _finite(value, name):
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} exceeds finite numerical range")
    return value


def _numeric(value, name, *, real=False):
    if np.ma.is_masked(value):
        raise ValueError(f"{name} must not contain masked or missing values")
    # Check mixed Python containers before NumPy coerces True to 1.
    if isinstance(value, (list, tuple)):
        for part in value:
            _numeric(part, name, real=real)
    raw = np.asarray(value)
    if raw.dtype.kind not in ("iuf" if real else "iufc"):
        raise ValueError(f"{name} must contain finite numeric values, not Booleans")
    _finite(raw, name)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = raw.astype(float if real else complex)
    _finite(result, name)
    if raw.dtype.kind in "iu":
        changed = any(int(x) != int(y) for x, y in zip(raw.flat, result.real.flat))
    else:
        changed = np.any(raw.real != result.real) or np.any(raw.imag != result.imag)
    if changed:
        raise ValueError(f"{name} loses information in float64 conversion")
    return result


def _max_component(a):
    return max(float(np.max(np.abs(a.real))), float(np.max(np.abs(a.imag))))


def _scaled(a, scale):
    # Real division avoids the complex reciprocal overflowing for tiny units.
    result = a.real/scale + 1j*(a.imag/scale)
    if (np.any((a.real != 0) & (result.real == 0))
            or np.any((a.imag != 0) & (result.imag == 0))):
        raise ValueError("observable normalization loses nonzero components")
    return _finite(result, "observable normalization")


def _observables(constraints):
    """Return S_a = (anchor_a + mean_a) I + scale_a C_a.

Subtract a diagonal anchor before scaling: a huge scalar must neither erase
off-diagonal entries nor dominate the Hermiticity check on the variation.
Retain the anchor and remaining mean separately until all scalar energies
are summed: other anchors can cancel and expose a much smaller mean.
"""
    if not isinstance(constraints, (list, tuple)) or not constraints:
        raise ValueError("nonempty constraint family required")
    centered, scales, offsets = [], [], []
    shape = None
    for value in constraints:
        a = _numeric(value, "constraints")
        if a.ndim != 2 or not a.shape[0] or a.shape[0] != a.shape[1]:
            raise ValueError("square constraint operators required")
        if shape is not None and a.shape != shape:
            raise ValueError("constraints must use the same algebra")
        shape = a.shape
        d = len(a)
        anchor = float(a[0, 0].real)
        b = a.copy()
        with np.errstate(over="ignore", invalid="ignore"):
            b[np.diag_indices(d)] -= anchor
        _finite(b, "observable differences")
        scale = _max_component(b)
        if scale == 0:
            centered.append(b)
            scales.append(1.)
            offsets.append((anchor, 0.))
            continue
        b = _scaled(b, scale)
        if np.linalg.norm(b-b.conj().T) > 1e-12:
            raise ValueError("constraints must be Hermitian")
        b = (b+b.conj().T)/2
        mean = float(np.trace(b).real/d)
        b[np.diag_indices(d)] -= mean
        centered.append(b)
        scales.append(scale)
        offsets.append((anchor, _finite(scale*mean, "observable offset")))
    return centered, np.array(scales), np.array(offsets)


def _thermal(ham):
    energies, vectors = np.linalg.eigh(ham)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        shifted = energies-energies[0]
        weights = np.exp(-shifted)
        probs = weights/weights.sum()
    if np.any(probs == 0) or not np.all(np.isfinite(probs)):
        raise ValueError("Gibbs spectrum underflow; faithful-state precision is insufficient")
    rho = faithful_density_matrix((vectors*probs) @ vectors.conj().T)
    return rho, math.log(weights.sum())-energies[0], probs, vectors


def _parameter(value, name):
    finite_real_scalar(value, name)
    return float(_numeric(value, name, real=True))


def gibbs_sectors(hamiltonians, central_energies, beta=1.):
    """Return (probability, conditional state) for each supplied sector.

    The full Hamiltonian is the direct sum of H_a + e_a I, with ordinary
    matrix trace and one common inverse temperature. Sector dimensions may
    differ. The scalar anchor, residual trace and central energy are kept
    separately; exact rational accumulation of those binary64 scalars removes
    a common origin before any rounding of relative log weights.

    Every sector and joint spectral direction must retain positive mass;
    unresolved dense-state support or underflow raises ValueError. beta=0
    returns the maximally mixed full state, so sector masses follow dimension.
    Negative beta is allowed for these finite spectra. No physical energy,
    temperature or center probabilities are inferred by this constructor.
    """
    if (not isinstance(hamiltonians, (list, tuple)) or not hamiltonians
            or not isinstance(central_energies, (list, tuple, np.ndarray))
            or (isinstance(central_energies, np.ndarray) and central_energies.ndim != 1)
            or len(hamiltonians) != len(central_energies)):
        raise ValueError("one central energy is required per nonempty sector")
    beta = _parameter(beta, "inverse temperature")
    states, spectra, log_weights = [], [], []
    for ham, energy in zip(hamiltonians, central_energies):
        energy = _parameter(energy, "central energy")
        centered, scales, offsets = _observables([ham])
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            coefficient = _finite(beta*scales[0], "scaled inverse temperature")
            variation = _finite(coefficient*centered[0], "thermal Hamiltonian")
        if beta != 0 and (coefficient == 0
                or np.any((centered[0].real != 0) & (variation.real == 0))
                or np.any((centered[0].imag != 0) & (variation.imag == 0))):
            raise ValueError("thermal Hamiltonian underflow; precision is insufficient")
        state, log_z, spectrum, _ = _thermal(variation)
        # Fractions act only on a few scalar offsets, not the eigensolver.
        # log_z already includes the sector's dimension and excitation gaps.
        scalar = sum((Fraction(float(x)) for x in offsets[0]), Fraction(energy))
        log_weights.append(Fraction(float(log_z))-Fraction(beta)*scalar)
        states.append(state)
        spectra.append(spectrum)

    origin = max(log_weights)
    try:
        differences = np.array([float(origin-value) for value in log_weights])
    except OverflowError as exc:
        raise ValueError("Gibbs sector underflow; relative energies exceed precision") from exc
    with np.errstate(under="ignore"):
        weights = np.exp(-differences)
        probabilities = weights/math.fsum(weights)
        if (np.any(probabilities == 0)
                or any(np.any(p*spectrum == 0)
                       for p, spectrum in zip(probabilities, spectra))):
            raise ValueError("normalized Gibbs sector underflow; precision is insufficient")
    return [(float(p), state) for p, state in zip(probabilities, states)]

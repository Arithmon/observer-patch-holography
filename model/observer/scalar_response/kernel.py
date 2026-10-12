"""Retained Weyl branches using all field modes and a small Gram matrix.

The field law and coherent preparation are the existing scalar adapter's
inputs. Free evolution is represented in the Heisenberg field directions;
the state keeps every branch amplitude and Weyl phase. Nothing is reset
between reads, and no oscillator occupation cutoff is introduced.

This is a bounded floating evaluator, not a roundoff certificate. The
separate continuation bounds concern the ideal mathematical field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import fsum
from numbers import Integral

import numpy as np

from .field import ScalarField, _array, _current, _outcome, _real


MAX_DIRECTIONS = 7
MAX_TERMS = 128
NORM_TOLERANCE = 1e-10


def _settings(times, couplings):
    if not isinstance(times, (list, tuple)) or not 1 <= len(times) <= MAX_DIRECTIONS:
        raise ValueError("one to seven ordered model times required")
    times = tuple(_real(x, "model time") for x in times)
    if any(t < 0 for t in times) or any(a > b for a, b in zip(times, times[1:])):
        raise ValueError("model times must be nonnegative and ordered")
    couplings = (1.,) * len(times) if couplings is None else couplings
    if not isinstance(couplings, (list, tuple)) or len(couplings) != len(times):
        raise ValueError("one real detector coupling per read required")
    return times, tuple(_real(x, "detector coupling") for x in couplings)


def _finite(value, name):
    if not np.all(np.isfinite(value)):
        raise ArithmeticError(f"nonfinite {name}")
    return value


@dataclass(frozen=True)
class WeylKernel:
    """Gaussian reference and inner products of the declared field directions.

    G_ab=<i v(t_a),i v(t_b)> and m_a=<X(t_a)> for unit preparation
    current. The signed preparation is current*m. A positive Hermitian G
    supplies the Gaussian CCR covariance; means are scalar translations.
    from_field binds these data to the actual existing scalar action.
    """

    covariance: np.ndarray
    unit_means: np.ndarray
    times: tuple
    current: int = 1
    mode_count: int = 0

    def __post_init__(self):
        times, _ = _settings(self.times, None)
        covariance = _array(self.covariance, "Weyl covariance", 2)
        means = _array(self.unit_means, "Weyl means", 1)
        n = len(times)
        if covariance.shape != (n, n) or means.shape != (n,) or np.any(means.imag != 0):
            raise ValueError("one real mean and square covariance per direction required")
        if not np.array_equal(covariance, covariance.conj().T):
            raise ValueError("Weyl covariance must be Hermitian")
        if np.any(covariance.diagonal().real < 0):
            raise ValueError("Weyl variances must be nonnegative")
        scale = max(float(np.max(np.abs(covariance.real), initial=0.)),
                    float(np.max(np.abs(covariance.imag), initial=0.)))
        if scale and np.linalg.eigvalsh(covariance / scale)[0] < -2e-12:
            raise ValueError("Weyl covariance must be positive semidefinite")
        if (isinstance(self.mode_count, (bool, np.bool_))
                or not isinstance(self.mode_count, Integral) or self.mode_count < 0):
            raise ValueError("nonnegative integer mode count required")
        covariance.setflags(write=False)
        real_means = means.real.copy()
        real_means.setflags(write=False)
        object.__setattr__(self, "covariance", covariance)
        object.__setattr__(self, "unit_means", real_means)
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "current", _current(self.current))
        object.__setattr__(self, "mode_count", int(self.mode_count))

    @classmethod
    def from_field(cls, field, times=(19/20, 1.), current=1, gamma=4., couplings=None):
        """Accumulate every mode in 2D slices; memory is O(q^2 + reads^2)."""
        if not isinstance(field, ScalarField):
            raise ValueError("the declared ScalarField is required")
        times, couplings = _settings(times, couplings)
        current, gamma = _current(current), _real(gamma, "gamma")
        n = len(times)
        g_sums = [[[] for _ in times] for _ in times]
        m_sums = [[] for _ in times]
        eigenvalues = field.laplacian_eigenvalues
        transverse_squared = np.outer(field.transverse_modes**2, field.transverse_modes**2)
        base = 1 + eigenvalues[:, None] + eigenvalues[None, :]
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            for lam, f, g in zip(eigenvalues, field.preparation_modes, field.detector_modes, strict=True):
                omega = np.sqrt(base + lam)
                variance = g*g*transverse_squared/(2*omega)
                for a, (ta, ca) in enumerate(zip(times, couplings)):
                    m_sums[a].append(float(np.sum(gamma*ca*f*g*transverse_squared*np.sin(omega*ta)/omega)))
                    for b in range(a, n):
                        phase = np.exp(1j*omega*(times[b]-ta))
                        g_sums[a][b].append(complex(np.sum(ca*couplings[b]*variance*phase)))
        covariance = np.zeros((n, n), complex)
        for a in range(n):
            for b in range(a, n):
                values = g_sums[a][b]
                value = complex(fsum(z.real for z in values), fsum(z.imag for z in values))
                covariance[a, b] = value
                covariance[b, a] = value.conjugate()
        means = np.array([fsum(values) for values in m_sums])
        _finite(covariance, "mode covariance")
        _finite(means, "mode mean")
        return cls(covariance, means, times, current, field.mode_count)

    @classmethod
    def from_modes(cls, omega, field_vector, coherent_initial, times, current=1, couplings=None):
        """Small explicit-mode control, independent of the tensor slice assembly."""
        times, couplings = _settings(times, couplings)
        omega = _array(omega, "frequencies", 1)
        beta = _array(field_vector, "field vector", 1)
        alpha = _array(coherent_initial, "coherent preparation", 1)
        if (not 1 <= len(omega) <= 32768 or beta.shape != omega.shape or alpha.shape != omega.shape
                or np.any(omega.imag != 0) or np.any(omega.real <= 0) or np.any(beta.imag != 0)):
            raise ValueError("bounded positive frequencies and matching real field vector required")
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            directions = np.array([1j*c*beta*np.exp(1j*omega.real*t)
                                   for t, c in zip(times, couplings)])
            covariance = directions.conj() @ directions.T
            means = 2*np.imag(directions @ alpha.conj())
        _finite(covariance, "mode covariance")
        _finite(means, "mode mean")
        # Fill the lower triangle from the same contractions, not separate
        # BLAS roundoff. Diagonal inner products are real by definition.
        for a in range(len(times)):
            covariance[a, a] = covariance[a, a].real
            for b in range(a+1, len(times)):
                covariance[b, a] = covariance[a, b].conjugate()
        return cls(covariance, means, times, current, len(omega))

    def for_current(self, current):
        return WeylKernel(self.covariance, self.unit_means, self.times, current, self.mode_count)

    def vacuum(self):
        """The prepared coherent reference, before any of the retained reads."""
        return WeylState(self, np.zeros((1, len(self.times)), dtype=int), np.array([1.]), -1)

    def statistics(self):
        return {"means": (self.current*self.unit_means).tolist(),
                "covariance_real": self.covariance.real.tolist(),
                "covariance_imag": self.covariance.imag.tolist(),
                "times": list(self.times), "mode_count": self.mode_count}

    def fingerprint(self):
        # A checkpoint belongs to this actual floating kernel. This is an
        # execution binding, not a portable numerical prediction digest.
        value = {**self.statistics(), "current": self.current,
                 "unit_means": self.unit_means.tolist()}
        return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()

    def two_read_law(self):
        """Analytic independent check; the experiment executes actual branches."""
        if len(self.times) != 2:
            raise ValueError("the analytic law requires exactly two reads")
        v1, v2 = self.covariance.diagonal().real
        cross = self.covariance[0, 1]
        m1, m2 = self.current*self.unit_means
        s1, s2 = np.exp(-v1/2)*np.sin(m1), np.exp(-v2/2)*np.sin(m2)
        attenuation = np.cos(cross.imag)
        parity = .5*(np.exp(-(v1+v2-2*cross.real)/2)*np.cos(m1-m2)
                     - np.exp(-(v1+v2+2*cross.real)/2)*np.cos(m1+m2))
        values = {(y, z): float((1+y*s1+z*attenuation*s2+y*z*parity)/4)
                  for y in (-1, 1) for z in (-1, 1)}
        if not all(np.isfinite(p) and 0 <= p <= 1 for p in values.values()):
            raise ArithmeticError("unresolved analytic joint probabilities")
        return values


@dataclass(frozen=True)
class WeylState:
    """Sum c_n D(i sum_a n_a v(t_a)/2)|alpha>, with all relative phases.

    Integer shifts encode half-field steps exactly. Coefficients and the
    small Gram matrix are floating diagnostics. A state admits at most
    seven ordered reads (128 terms); resource exhaustion is explicit.
    """

    kernel: WeylKernel
    shifts: np.ndarray
    coefficients: np.ndarray
    last_index: int = -1

    def __post_init__(self):
        if not isinstance(self.kernel, WeylKernel) or np.ma.isMaskedArray(self.shifts):
            raise ValueError("a Weyl kernel and unmasked integer shifts required")
        raw_shifts = np.asarray(self.shifts, dtype=object)
        if (raw_shifts.ndim != 2 or any(isinstance(x, (bool, np.bool_))
                or not isinstance(x, Integral) for x in raw_shifts.flat)):
            raise ValueError("half-shifts must be integers, not boolean or rounded symbols")
        shifts = np.asarray(self.shifts)
        coefficients = _array(self.coefficients, "Weyl coefficients", 1)
        if (shifts.ndim != 2 or shifts.dtype.kind not in "iu"
                or not 1 <= len(shifts) <= MAX_TERMS or shifts.shape[1] != len(self.kernel.times)
                or coefficients.shape != (len(shifts),) or np.any(shifts > 7) or np.any(shifts < -7)):
            raise ValueError("bounded integer half-shifts and matching coefficients required")
        if (isinstance(self.last_index, (bool, np.bool_)) or not isinstance(self.last_index, Integral)
                or not -1 <= self.last_index < len(self.kernel.times)):
            raise ValueError("invalid last-read index")
        if np.any(shifts[:, self.last_index+1:] != 0):
            raise ValueError("state contains shifts from future reads")
        shifts = shifts.astype(np.int64, copy=True)
        shifts.setflags(write=False)
        coefficients.setflags(write=False)
        object.__setattr__(self, "shifts", shifts)
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "last_index", int(self.last_index))

    def gram(self):
        half = self.shifts/2
        result = np.empty((len(half), len(half)), complex)
        means = self.kernel.current*self.kernel.unit_means
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            for a, left in enumerate(half):
                for b, right in enumerate(half):
                    difference = right-left
                    distance = float(difference @ self.kernel.covariance.real @ difference)
                    phase = float(left @ self.kernel.covariance.imag @ right)
                    mean_phase = float(difference @ means)
                    _finite((distance, phase, mean_phase), "Weyl Gram exponent")
                    exponent = -.5*distance + 1j*phase + 1j*mean_phase
                    result[a, b] = np.exp(exponent)
        return _finite(result, "Weyl Gram matrix")

    def norm_squared(self):
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            norm = np.vdot(self.coefficients, self.gram() @ self.coefficients)
        _finite(norm, "Weyl branch norm")
        if abs(norm.imag) > NORM_TOLERANCE or norm.real < 0:
            raise ArithmeticError("unresolved Weyl branch norm")
        return float(norm.real)

    def normalized(self):
        norm = self.norm_squared()
        if norm <= 1e-14:
            raise ArithmeticError("zero or unresolved Weyl branch cannot be normalized")
        state = WeylState(self.kernel, self.shifts, self.coefficients/np.sqrt(norm), self.last_index)
        if abs(state.norm_squared()-1) > NORM_TOLERANCE:
            raise ArithmeticError("unresolved normalized Weyl state")
        return state

    def instrument(self, index, outcome):
        """Apply the next retained binary instrument to this actual state."""
        outcome = _outcome(outcome)
        if (isinstance(index, (bool, np.bool_)) or not isinstance(index, Integral)
                or not self.last_index < index < len(self.kernel.times)):
            raise ValueError("read index must advance within the declared history")
        if 2*len(self.coefficients) > MAX_TERMS:
            raise ValueError("retained Weyl term budget exceeded")
        if abs(self.norm_squared()-1) > NORM_TOLERANCE:
            raise ValueError("normalized premeasurement Weyl state required")
        shifts, coefficients = [], []
        for sign, factor in ((-1, .5), (1, -.5j*outcome)):
            moved = self.shifts.copy()
            moved[:, index] += sign
            with np.errstate(over="ignore", invalid="ignore", under="ignore"):
                phases = np.exp(-.25j*sign*(self.shifts @ self.kernel.covariance.imag[index]))
            # h_new=sign/2 and h_old=n/2 give phase -i sign G_im*n/4.
            shifts.append(moved)
            coefficients.append(self.coefficients*factor*_finite(phases, "Weyl phase"))
        branch = WeylState(self.kernel, np.concatenate(shifts), np.concatenate(coefficients), int(index))
        probability = branch.norm_squared()
        if not 0 < probability <= 1:
            raise ArithmeticError("zero or unresolved instrument probability")
        return probability, branch.normalized()

    def checkpoint(self):
        return {"schema": "oph.scalar_weyl_state.v1", "kernel": self.kernel.fingerprint(),
                "last_index": self.last_index, "half_shifts": self.shifts.tolist(),
                "coefficients": [[float(z.real), float(z.imag)] for z in self.coefficients]}

    @classmethod
    def restore(cls, kernel, payload):
        keys = {"schema", "kernel", "last_index", "half_shifts", "coefficients"}
        if (not isinstance(kernel, WeylKernel) or not isinstance(payload, dict) or set(payload) != keys
                or payload["schema"] != "oph.scalar_weyl_state.v1"
                or payload["kernel"] != kernel.fingerprint()):
            raise ValueError("checkpoint does not match its actual field kernel")
        pairs = _array(payload["coefficients"], "checkpoint coefficients", 2)
        if pairs.shape[1:] != (2,) or np.any(pairs.imag != 0):
            raise ValueError("real/imaginary coefficient pairs required")
        state = cls(kernel, payload["half_shifts"], pairs[:, 0].real+1j*pairs[:, 1].real,
                    payload["last_index"])
        if abs(state.norm_squared()-1) > NORM_TOLERANCE:
            raise ValueError("checkpoint is not a normalized state")
        return state

"""Checked PDG mixing coordinates; historical schemas are adapters below.

Angles are in the first quadrant. A full Dirac-phase readout requires all
five PDG anchor entries to be nonzero; at a chart boundary it is undefined.
See README.md for the determinant identity and the numerical contract.
"""

from fractions import Fraction
import math

import numpy as np

from quantum_information.gibbs import _numeric


def _readout(matrix):
    u = _numeric(matrix, "mixing matrix")
    if u.shape != (3, 3):
        raise ValueError("mixing matrix must have shape (3, 3)")
    with np.errstate(over="ignore", invalid="ignore"):
        residual = np.linalg.norm(u.conj().T @ u - np.eye(3), ord="fro")
    if not np.isfinite(residual) or residual > 1e-12:
        raise ValueError("mixing matrix must be unitary to absolute tolerance 1e-12")
    anchors = (u[0, 0], u[0, 1], u[1, 2], u[2, 2], u[0, 2].conjugate())
    phase = 1 + 0j
    for z in (*anchors, np.linalg.det(u).conjugate()):
        scale = max(abs(z.real), abs(z.imag))
        if scale == 0:
            raise ValueError("Dirac phase is undefined at a mixing-chart boundary")
        # Normalize each factor separately: neither a tiny product nor a
        # complex reciprocal may underflow/overflow and erase the phase.
        z = complex(z.real / scale, z.imag / scale)
        phase *= z / abs(z)
    signed_delta = math.atan2(phase.imag, phase.real)
    angles = (math.atan2(abs(u[0, 1]), abs(u[0, 0])),
              math.atan2(abs(u[1, 2]), abs(u[2, 2])),
              math.atan2(abs(u[0, 2]), math.hypot(abs(u[0, 0]), abs(u[0, 1]))))
    # Exact quartet on the supplied binary64 entries, then one rounding.
    real, imag = Fraction(1), Fraction(0)
    for z in (u[0, 0], u[1, 1], u[0, 1].conjugate(), u[1, 0].conjugate()):
        a, b = Fraction(float(z.real)), Fraction(float(z.imag))
        real, imag = real*a - imag*b, real*b + imag*a
    jarlskog = float(imag)
    if imag and jarlskog == 0:
        raise ValueError("Jarlskog invariant underflows binary64")
    return (*angles, signed_delta), jarlskog


def mixing_parameters(matrix):
    """Nine canonical fields, with delta in [0, 2*pi); no input repair."""
    values, jarlskog = _readout(matrix)
    delta = values[3] % math.tau
    values = (*values[:3], 0.0 if delta == math.tau else delta)
    names = ("theta12", "theta23", "theta13", "delta")
    return {**{name + "_rad": value for name, value in zip(names, values)},
            **{name + "_deg": math.degrees(value) for name, value in zip(names, values)},
            "J": jarlskog}


def pmns_degrees(matrix):
    return {key: value for key, value in mixing_parameters(matrix).items() if not key.endswith("_rad")}


def ckm_parameters(matrix):
    result = mixing_parameters(matrix)
    return dict(zip(("theta_12", "theta_23", "theta_13", "delta_ckm", "jarlskog"),
                    (result[key] for key in ("theta12_rad", "theta23_rad", "theta13_rad", "delta_rad", "J"))))


def pmns_signed(matrix):
    values, jarlskog = _readout(matrix)
    return dict(zip(("theta_12", "theta_23", "theta_13", "delta_pmns", "jarlskog"), (*values, jarlskog)))

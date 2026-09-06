"""Read a model-internal clock from ordered charged configurations.

The consumer reads only q_exact from a separately authenticated software
observer episode. Neither velocity, action timestamp nor repair-cycle count
enters the clock integral. Geometry, action, energy and temporal gauge are
declared inputs. Polyline and floating-point errors are numerical diagnostics,
not a rigorous trajectory enclosure or a laboratory clock calibration.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE = HERE / "runtime/whitney_charged_instrument_receipt.json"
OUTPUT = HERE / "runtime/whitney_ephemeris_clock_receipt.json"
ENERGY = 493*np.sqrt(5)/120+1969/160
PINS = (
    "code/electromagnetism/whitney_ephemeris_clock.py",
    "code/electromagnetism/verify_whitney_ephemeris_clock.py",
    "code/electromagnetism/test_whitney_ephemeris_clock.py",
    "paper/tex_fragments/WHITNEY_EPHEMERIS_CLOCK.tex",
    "code/electromagnetism/verify_whitney_charged_dynamics.py",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
)


def action_at(q, tangent):
    """Full restricted temporal kinetic energy, including basis derivative."""
    alpha, cr, ci, br, bi = np.asarray(q, dtype=float)
    da, dcr, dci, dbr, dbi = np.asarray(tangent, dtype=float)
    c, b, dc, db = cr+1j*ci, br+1j*bi, dcr+1j*dci, dbr+1j*dbi
    phi = (1+np.sqrt(5))/2
    volume = 20*phi**2/3
    kappa = 3*volume/(2+3*phi)
    nodes, weights = leggauss(4)
    s = (nodes+1)/2
    weights = weights*1.5*volume*(1-s)**2
    rotation = np.exp(.25j*alpha)
    difference = rotation*c-b
    value = s*rotation*c+(1-s)*b
    derivative = s*rotation*dc+(1-s)*db+.25j*da*s*(1-s)*difference
    kinetic = .5*kappa*da**2+weights@abs(derivative)**2
    potential = kappa*abs(difference)**2+weights@(.5*abs(value)**2+.125*abs(value)**4)
    return float(kinetic), float(potential)


def durations(configurations, energy=ENERGY, order=24, warped=False):
    """Jacobi duration on each polygon segment; labels carry no duration."""
    q = np.asarray(configurations, dtype=float)
    if q.ndim != 2 or q.shape[1] != 5 or len(q) < 2 or not np.isfinite(q).all():
        raise ValueError("finite sequence of at least two five-coordinate configurations required")
    if not np.isfinite(energy):
        raise ValueError("finite supplied energy required")
    nodes, weights = leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    result = []
    for first, second in zip(q[:-1], q[1:], strict=True):
        delta = second-first
        if not np.any(delta):
            raise ValueError("nonregular constant segment")
        value = 0.0
        for u, weight in zip(nodes, weights, strict=True):
            # Smooth positive-speed regrading of the SAME polygon segment.
            position = (u+u*u)/2 if warped else u
            rate = .5+u if warped else 1.0
            kinetic, potential = action_at(first+position*delta, rate*delta)
            if energy <= potential or kinetic <= 0:
                raise ValueError("clock requires positive kinetic energy and an open Hill region")
            value += weight*np.sqrt(kinetic/(energy-potential))
        result.append(float(value))
    return result


def read_configurations(packet):
    """Deliberately no access to timestamps, velocities, or cycle counts."""
    return [[float(Fraction(value)) for value in frame["q_exact"]] for frame in packet["frames"]]


def build():
    import verify_whitney_charged_instrument as instrument
    source = instrument.load(SOURCE)
    instrument.verify(source)
    q = read_configurations(source)
    series = []
    for stride in (8, 4, 2, 1):
        increments = durations(q[::stride])
        series.append({"stride": stride, "configurations": len(q[::stride]),
                       "increments": increments, "duration": float(sum(increments))})
    fine = series[-1]["increments"]
    return {
        "schema": "oph.whitney_ephemeris_clock.v1",
        "scope": "ACTION_DERIVED_MODEL_CLOCK__NUMERICAL_POLYLINE_READOUT",
        "source": {"path": SOURCE.relative_to(ROOT).as_posix(),
                   "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                   "consumed_fields": ["frames[].q_exact"]},
        "source_pins": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
        "contract": {"energy_exact": "1969/160 + (493/120)*sqrt(5)",
                     "formula": "d_tau=sqrt(G(q)[dq,dq]/(2*(E-V(q))))",
                     "gauge": "temporal gauge of the decoded classical path",
                     "input_clock": False, "calibrated_physical_clock": False,
                     "quantum_clock_operator": False, "rigorous_numerical_enclosure": False,
                     "energy_source": "supplied initial condition and the same complete coupled action",
                     "curve": "piecewise linear interpolation of decoded configurations",
                     "quadrature": "24-point Gauss-Legendre per segment; 32-point independent-order controls"},
        "refinements": series,
        "controls": {"warped_duration": sum(durations(q, order=32, warped=True)),
                     "higher_order_duration": sum(durations(q, order=32)),
                     "wrong_energy_duration": sum(durations(q, energy=ENERGY+1)),
                     "solver_interval_for_comparison_only": 2.0},
        "cumulative_clock": [0.0]+np.cumsum(fine).tolist(),
    }


if __name__ == "__main__":
    packet = build()
    OUTPUT.write_bytes((json.dumps(packet, sort_keys=True, indent=2, allow_nan=False)+"\n").encode("utf-8"))
    print(OUTPUT)

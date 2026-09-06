"""Independent action-metric replay of a configuration-only model clock.

Reconstructs full local tetrahedral fields, not the clock producer's reduced
one-dimensional formulas. Provenance of the source is replayed separately.
Numerical tolerances do not certify floating-point or physical clock errors.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def peer(name):
    spec = importlib.util.spec_from_file_location("ephemeris_"+name, HERE/(name+".py"))
    if spec is None or spec.loader is None:
        raise ValueError("missing sibling verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fields = peer("verify_whitney_charged_dynamics")


def instrument_verifier():
    return peer("verify_whitney_charged_instrument")


OUTPUT = HERE / "runtime/whitney_ephemeris_clock_receipt.json"
SOURCE = HERE / "runtime/whitney_charged_instrument_receipt.json"
ENERGY = float(Fraction(1969, 160))+float(Fraction(493, 120))*np.sqrt(5)
PINS = {
    "code/electromagnetism/whitney_ephemeris_clock.py",
    "code/electromagnetism/verify_whitney_ephemeris_clock.py",
    "code/electromagnetism/test_whitney_ephemeris_clock.py",
    "paper/tex_fragments/WHITNEY_EPHEMERIS_CLOCK.tex",
    "code/electromagnetism/verify_whitney_charged_dynamics.py",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path=OUTPUT):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result
    def number(text):
        result = float(text)
        require(np.isfinite(result), "nonfinite JSON number")
        return result
    def constant(_):
        raise ValueError("nonfinite JSON constant")
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_float=number, parse_constant=constant)


def equal(value, expected, name):
    require(json.dumps(value, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), name)


def near(value, expected, name, tolerance=2e-10):
    def numeric(x):
        if isinstance(x, list):
            return all(numeric(v) for v in x)
        return type(x) in (int, float) and np.isfinite(x)
    require(numeric(value), "finite numeric "+name)
    a, b = np.asarray(value), np.asarray(expected)
    require(a.shape == b.shape and np.max(abs(a-b), initial=0) <= tolerance, name)


def action_at(q, tangent):
    mesh = fields.geometry(4)
    a, v = fields.expanded(q), fields.expanded(tangent)
    element = mesh["elements"][0]
    data = fields.element_fields(a, v, element)
    weight = 20*element["weights"]
    kinetic = .5*v[:42]@mesh["M"]@v[:42]+weight@abs(data["scalar_dot"])**2
    potential = .5*a[:42]@mesh["K"]@a[:42]
    potential += weight@(np.sum(abs(data["spatial"])**2, axis=1)
                         +.5*abs(data["scalar"])**2+.125*abs(data["scalar"])**4)
    return float(kinetic), float(potential)


def increments(configurations, energy=ENERGY, warped=False):
    result = []
    for left, right in zip(configurations[:-1], configurations[1:], strict=True):
        left, right = np.asarray(left), np.asarray(right)
        delta = right-left
        require(np.any(delta), "nonregular constant segment")
        def density(u):
            position = (u+u*u)/2 if warped else u
            speed = .5+u if warped else 1
            kinetic, potential = action_at(left+position*delta, speed*delta)
            require(energy > potential and kinetic > 0, "open Hill region and regularity")
            return np.sqrt(kinetic/(energy-potential))
        value, error = quad(density, 0, 1, epsabs=2e-13, epsrel=2e-13)
        require(error < 1e-10, "quadrature diagnostic")
        result.append(value)
    return result


def verify(packet):
    require(type(packet) is dict and set(packet) == {
        "schema", "scope", "source", "source_pins", "contract", "refinements",
        "controls", "cumulative_clock"}, "clock schema keys")
    equal(packet["schema"], "oph.whitney_ephemeris_clock.v1", "schema")
    equal(packet["scope"], "ACTION_DERIVED_MODEL_CLOCK__NUMERICAL_POLYLINE_READOUT", "scope")
    equal(packet["contract"], {
        "energy_exact": "1969/160 + (493/120)*sqrt(5)",
        "formula": "d_tau=sqrt(G(q)[dq,dq]/(2*(E-V(q))))",
        "gauge": "temporal gauge of the decoded classical path",
        "input_clock": False, "calibrated_physical_clock": False,
        "quantum_clock_operator": False, "rigorous_numerical_enclosure": False,
        "energy_source": "supplied initial condition and the same complete coupled action",
        "curve": "piecewise linear interpolation of decoded configurations",
        "quadrature": "24-point Gauss-Legendre per segment; 32-point independent-order controls"}, "clock contract")
    require(type(packet["source_pins"]) is dict and set(packet["source_pins"]) == PINS, "pin census")
    for path in PINS:
        equal(packet["source_pins"][path], hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), "pin "+path)
    equal(packet["source"], {"path": SOURCE.relative_to(ROOT).as_posix(),
          "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          "consumed_fields": ["frames[].q_exact"]}, "source attachment")
    instrument = instrument_verifier()
    source = instrument.load(SOURCE)
    instrument.verify(source)
    configurations = [[float(Fraction(x)) for x in frame["q_exact"]] for frame in source["frames"]]
    require(len(configurations) == 81, "configuration census")
    require(type(packet["refinements"]) is list and len(packet["refinements"]) == 4, "refinement census")
    durations = []
    for row, stride in zip(packet["refinements"], (8, 4, 2, 1), strict=True):
        require(type(row) is dict and set(row) == {"stride", "configurations", "increments", "duration"}, "refinement keys")
        equal(row["stride"], stride, "stride")
        equal(row["configurations"], 80//stride+1, "configuration count")
        values = increments(configurations[::stride])
        near(row["increments"], values, "independent segment durations")
        near(row["duration"], sum(values), "independent total duration")
        durations.append(sum(values))
    fine = values
    controls = packet["controls"]
    require(type(controls) is dict and set(controls) == {"warped_duration", "higher_order_duration",
            "wrong_energy_duration", "solver_interval_for_comparison_only"}, "control census")
    near(controls["warped_duration"], sum(increments(configurations, warped=True)), "reparameterization")
    near(controls["higher_order_duration"], durations[-1], "quadrature-order comparison")
    wrong = sum(increments(configurations, ENERGY+1))
    near(controls["wrong_energy_duration"], wrong, "wrong energy control")
    require(abs(wrong-durations[-1]) > .05, "nonvacuous energy dependence")
    equal(controls["solver_interval_for_comparison_only"], 2.0, "comparison-only interval")
    near(packet["cumulative_clock"], [0.0]+np.cumsum(fine).tolist(), "cumulative clock")
    errors = [abs(x-2) for x in durations]
    require(all(a > b for a, b in zip(errors[:-1], errors[1:])), "finite refinement control")
    return {"accepted": True, "duration": durations[-1], "comparison_error": errors[-1],
            "source_configurations": 81, "calibrated_physical_clock": False,
            "rigorous_numerical_enclosure": False, "quantum_clock_operator": False}


if __name__ == "__main__":
    print(json.dumps(verify(load()), sort_keys=True))

"""Conditional SI dictionary for a specified model clock and angular frequency.

Fixing c sets a length/time ratio. The seconds per model clock unit and the
identification E=hbar*omega are declared calibration inputs, not model results.
No repair count is treated as an angular frequency or a particle mass.
"""
from __future__ import annotations

import argparse
import json
import math

# Exact SI defining constants; floating evaluation below is approximate.
# https://www.bipm.org/en/measurement-units/si-defining-constants
C_M_PER_S = 299_792_458
H_J_S = 6.62607015e-34
HBAR_J_S = H_J_S / math.tau


def attachment(seconds_per_unit: float, signal_speed: float, *,
               rest_angular_frequency: float | None = None,
               kinetic_curvature: float | None = None) -> dict:
    """Convert a supplied dispersion in dimensionless model coordinates.

    omega is radians per model-clock unit; k is radians per model-length unit.
    curvature is d²omega/dk² at rest, independently extracted from dispersion.
    The dimensionless massive-shell check is omega(0)*curvature/signal_speed².
    """
    for name, value in (("seconds_per_unit", seconds_per_unit),
                        ("signal_speed", signal_speed)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if rest_angular_frequency is not None and (
            not math.isfinite(rest_angular_frequency) or rest_angular_frequency < 0):
        raise ValueError("rest angular frequency must be finite and nonnegative")
    if kinetic_curvature is not None and (
            not math.isfinite(kinetic_curvature) or kinetic_curvature < 0
            or rest_angular_frequency is None or rest_angular_frequency <= 0):
        raise ValueError("massive curvature check needs positive rest frequency and finite nonnegative curvature")
    energy_unit = HBAR_J_S / seconds_per_unit
    mass_unit = energy_unit / C_M_PER_S**2
    length_unit = C_M_PER_S * seconds_per_unit / signal_speed
    speed_squared = signal_speed * signal_speed
    if any(not math.isfinite(value) or value <= 0
           for value in (energy_unit, mass_unit, length_unit, speed_squared)):
        raise ValueError("unit scales exceed the finite floating arithmetic range")
    result = {
        "status": "CONDITIONAL_UNIT_ATTACHMENT",
        "constants_source": "https://www.bipm.org/en/measurement-units/si-defining-constants",
        "seconds_per_model_clock_unit_supplied": seconds_per_unit,
        "model_signal_speed_supplied": signal_speed,
        "metres_per_model_length_unit": length_unit,
        "joules_per_unit_angular_frequency": energy_unit,
        "kilograms_per_unit_rest_angular_frequency": mass_unit,
        "conditional_rest_energy_joules": None if rest_angular_frequency is None else energy_unit * rest_angular_frequency,
        "conditional_rest_mass_kg": None if rest_angular_frequency is None else mass_unit * rest_angular_frequency,
        "independent_mass_shell_ratio": None if kinetic_curvature is None else rest_angular_frequency * kinetic_curvature / speed_squared,
        "energy_frequency_identification_supplied": "E = hbar * omega; rest mass = E_rest / c^2",
        "physical_clock_calibration_derived": False,
        "candidate_particle_mass_established": False,
    }
    if any(isinstance(value, float) and not math.isfinite(value) for value in result.values()):
        raise ValueError("converted mode exceeds the finite floating arithmetic range")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds-per-unit", type=float, required=True)
    parser.add_argument("--signal-speed", type=float, required=True)
    parser.add_argument("--rest-angular-frequency", type=float)
    parser.add_argument("--kinetic-curvature", type=float)
    args = parser.parse_args()
    print(json.dumps(attachment(args.seconds_per_unit, args.signal_speed,
        rest_angular_frequency=args.rest_angular_frequency,
        kinetic_curvature=args.kinetic_curvature), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

"""Exact clock-scale ambiguity for the native transfer's active modular probe.

This is a local, target-free sanity check.  The source transfer determines
dimensionless modular gaps and the declared controller applies one modular
step per *incident event*.  Neither operation assigns seconds to that event.
The transfer gaps are not identified here with particle rest energies.
"""

from __future__ import annotations

from fractions import Fraction

from . import construction


def _positive_fraction(value: Fraction, name: str) -> Fraction:
    if type(value) is not Fraction or value <= 0:
        raise ValueError(f"{name} must be a positive exact Fraction")
    return value


def _scale_interval(interval: tuple[Fraction, Fraction], scale: Fraction):
    lo, hi = interval
    if not 0 < lo <= hi:
        raise ArithmeticError("a strictly positive certified modular gap is required")
    return lo / scale, hi / scale


def clock_scale_witness(
    seconds_per_event_a: Fraction = Fraction(1),
    seconds_per_event_b: Fraction = Fraction(2),
) -> dict[str, object]:
    """Exhibit two distinct physical-frequency scales for the same source law.

    For each positive event duration tau, the assigned angular frequency is
    omega_j=E_j/tau.  At the only native times t_n=n*tau, both assignments
    give the identical phase omega_j*t_n=n*E_j.  A conditional rest-energy
    interpretation would give hbar*omega_j and inherit the same ambiguity.
    This construction neither chooses tau nor claims such an interpretation.
    """

    tau_a = _positive_fraction(seconds_per_event_a, "seconds_per_event_a")
    tau_b = _positive_fraction(seconds_per_event_b, "seconds_per_event_b")
    if tau_a == tau_b:
        raise ValueError("the two event-duration assignments must differ")

    source = construction.spectrum_audit()
    zero, slow, fast = construction.energies()
    if zero != (Fraction(0), Fraction(0)):
        raise ArithmeticError("the source's stationary modular mode changed")
    if not (0 < slow[0] <= slow[1] < fast[0] <= fast[1]):
        raise ArithmeticError("the two source modular gaps lost ordering")

    ratio = (fast[0] / slow[1], fast[1] / slow[0])
    a = (_scale_interval(slow, tau_a), _scale_interval(fast, tau_a))
    b = (_scale_interval(slow, tau_b), _scale_interval(fast, tau_b))

    # A common multiplication changes both frequencies but cancels from the
    # ratio.  It cannot be selected from the event-indexed phase or current law.
    for original, rate_a, rate_b in zip((slow, fast), a, b, strict=True):
        for endpoint, omega_a, omega_b in zip(original, rate_a, rate_b, strict=True):
            if omega_a * tau_a != endpoint or omega_b * tau_b != endpoint:
                raise ArithmeticError("a clock embedding changed native phases")
            if omega_a * tau_a != omega_b * tau_b:
                raise ArithmeticError("clock embeddings no longer share event phases")
    return {
        "source_characteristic_polynomial": source["characteristic_polynomial"],
        "modular_gaps_per_declared_incident_event": (slow, fast),
        "dimensionless_gap_ratio_enclosure": ratio,
        "seconds_per_event": (tau_a, tau_b),
        "angular_frequency_per_second_enclosures": (a, b),
        "same_native_event_phases": True,
        "physical_energy_or_rest_mass_selected": False,
        "conditional_energy_formula": "E_j(SI)=hbar*gap_j/seconds_per_event",
        "source_clock_is_attempt_transfer_but_probe_ticking_is_incident_event": True,
    }

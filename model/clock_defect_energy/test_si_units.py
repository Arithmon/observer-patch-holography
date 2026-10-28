import math

import pytest

from .si_units import C_M_PER_S, attachment


def test_unit_attachment_does_not_invent_a_rest_frequency():
    result = attachment(1.0, 1.0)
    assert result["conditional_rest_mass_kg"] is None
    assert result["metres_per_model_length_unit"] == C_M_PER_S
    assert not result["candidate_particle_mass_established"]


def test_clock_unit_change_preserves_a_consistently_rescaled_physical_mode():
    # Same physical dispersion, with a clock unit ten times as long.
    first = attachment(0.01, 2.0, rest_angular_frequency=3.0, kinetic_curvature=4/3)
    second = attachment(0.1, 20.0, rest_angular_frequency=30.0, kinetic_curvature=40/3)
    for field in ("conditional_rest_mass_kg", "conditional_rest_energy_joules",
                  "metres_per_model_length_unit", "independent_mass_shell_ratio"):
        assert math.isclose(first[field], second[field], rel_tol=1e-14)
    assert math.isclose(first["independent_mass_shell_ratio"], 1)
    assert math.isclose(first["conditional_rest_energy_joules"],
                        first["conditional_rest_mass_kg"] * C_M_PER_S**2)


def test_unit_conversion_cannot_fix_a_failed_mass_shell_check():
    assert attachment(1, 2, rest_angular_frequency=3,
                      kinetic_curvature=1)["independent_mass_shell_ratio"] == 0.75


@pytest.mark.parametrize("seconds,speed", [(0,1),(-1,1),(1,0),(1,-1),(float('nan'),1),(1,float('inf'))])
def test_invalid_clock_or_signal_scales_are_rejected(seconds, speed):
    with pytest.raises(ValueError):
        attachment(seconds, speed)


def test_massive_curvature_needs_a_rest_gap():
    with pytest.raises(ValueError, match="positive rest frequency"):
        attachment(1, 1, kinetic_curvature=1)


@pytest.mark.parametrize("seconds,speed", [(1e308,1),(1,1e-320),(1,1e308)])
def test_unrepresentable_derived_scales_are_rejected(seconds, speed):
    with pytest.raises(ValueError, match="floating arithmetic range"):
        attachment(seconds, speed, rest_angular_frequency=1, kinetic_curvature=1)

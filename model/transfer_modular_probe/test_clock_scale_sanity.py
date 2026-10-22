"""The same native event history admits distinct SI clock embeddings."""

from fractions import Fraction as F
from math import log, sqrt

import pytest

from . import clock_scale_sanity, construction
from ..native_clock import clock


def test_source_gaps_and_native_record_law_do_not_select_seconds():
    witness = clock_scale_sanity.clock_scale_witness(F(1, 100), F(1, 50))
    gaps = witness["modular_gaps_per_declared_incident_event"]
    assert gaps == construction.energies()[1:]
    assert witness["source_characteristic_polynomial"] == "(x-1)(x^2-5x/3+311/450)"

    # The exact marked source law has no duration argument, including the
    # no-current branch.  Both embeddings thus predict this complete word law.
    reduced = construction.reduced()
    word = (1, 0, 1)
    native_probability = sum(clock.word_law(reduced, word))
    assert native_probability == F(53, 2700)
    assert witness["same_native_event_phases"] is True
    assert witness["physical_energy_or_rest_mass_selected"] is False

    rates_a, rates_b = witness["angular_frequency_per_second_enclosures"]
    for rate_a, rate_b in zip(rates_a, rates_b, strict=True):
        assert tuple(omega_a / omega_b for omega_a, omega_b in zip(rate_a, rate_b)) == (F(2), F(2))
    ratio_lo, ratio_hi = witness["dimensionless_gap_ratio_enclosure"]
    expected_ratio = -log((25 - sqrt(3)) / 30) / -log((25 + sqrt(3)) / 30)
    assert 1 < ratio_lo <= ratio_hi
    assert abs(float(ratio_lo) - expected_ratio) < 1e-14
    assert abs(float(ratio_hi) - expected_ratio) < 1e-14


@pytest.mark.parametrize("a,b", [(F(0), F(1)), (F(-1), F(1)), (F(1), F(1)), (1, F(2))])
def test_invalid_or_identical_clock_assignments_fail(a, b):
    with pytest.raises(ValueError):
        clock_scale_sanity.clock_scale_witness(a, b)

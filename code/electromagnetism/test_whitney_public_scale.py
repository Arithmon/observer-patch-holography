"""Original-input public Gaussian logs at binary64 range boundaries.

The logarithm can be representable even when an intermediate squared norm
is not. The oracle retains the original supplied floats before squaring or
dividing; it does not use a producer norm, metric, or density value.
"""
import warnings
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_interacting_quantum as quantum
import whitney_quantum_packet as packet


MINIMUM = float.fromhex("0x0.0000000000001p-1022")


def original_log(coordinates, sigma, curved):
    mp = mpmath.mp.clone()
    mp.dps = 100
    exponent = mp.fsum(mp.mpf(float(value))**2 for value in coordinates)
    exponent /= 4*mp.mpf(float(sigma))**2
    result = -14*(mp.log(2*mp.pi)+2*mp.log(mp.mpf(float(sigma))))-exponent
    if curved:
        # At zero charge the metric is independent of the matter coordinates.
        # This absolute determinant is independently derived from the golden
        # cone's exact Whitney/simplex masses in test_whitney_kinetic_precision.
        neutral_logdet = mp.mpf(
            "-59.981818227201022084204145712934299528443641467366830477835"
        )
        result -= neutral_logdet/4
    return result


def evaluate(coordinates, sigma, curved):
    if curved:
        return quantum.gaussian_state_log_amplitude(coordinates, sigma, charge=0)
    return quantum.gaussian_half_density_log(coordinates, sigma)


@pytest.mark.parametrize("curved", [False, True], ids=["half_density", "neutral_state"])
@pytest.mark.parametrize("values,sigma", [
    ([2.], 1.),
    ([2e154], 1.),
    ([1.5e154, 1.5e154], 1.),
    ([1e-169], MINIMUM),
    ([1e308], 1e154),
    ([1e308], 1e308),
    ([], MINIMUM),
], ids=["ordinary", "square_overflow", "sum_overflow", "subnormal_width",
        "large_balanced_units", "huge_width", "zero_at_small_width"])
def test_representable_original_gaussian_log_remains_evaluable(values, sigma, curved):
    coordinates = np.zeros(56)
    coordinates[30:30+len(values)] = values
    expected = float(original_log(coordinates, sigma, curved))
    assert np.isfinite(expected)
    # Capture the defective implementation's arithmetic warning so the
    # retained baseline fails on its wrong result, not a warning category.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = evaluate(coordinates, sigma, curved)
    assert np.isfinite(actual), "a representable original-input log became infinite"
    assert actual == pytest.approx(expected, rel=2e-14, abs=2e-12)
    assert not caught


@pytest.mark.parametrize("curved", [False, True], ids=["half_density", "neutral_state"])
@pytest.mark.parametrize("value,sigma", [(3e154, 1.), (1., MINIMUM)],
                         ids=["unreportable_log", "unreportable_ratio"])
def test_unreportable_gaussian_log_is_an_explicit_resolution_failure(value, sigma, curved):
    coordinates = np.zeros(56)
    coordinates[30] = value
    expected = original_log(coordinates, sigma, curved)
    assert abs(expected) > np.finfo(float).max
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="precision|range|represent"):
            evaluate(coordinates, sigma, curved)


@pytest.mark.parametrize("indices,values", [
    ([30, 31], [1.5e308, -1.5e308]),
    ([30, 43], [1.5e308, 1.5e308]),
], ids=["opposite_neutral_nodes", "finite_complex_components"])
def test_neutral_gaussian_does_not_evaluate_unused_charged_differences(indices, values):
    coordinates = np.zeros(56)
    coordinates[indices] = values
    sigma = 1.5e308
    expected = float(original_log(coordinates, sigma, True))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = quantum.gaussian_state_log_amplitude(coordinates, sigma, charge=0)
    assert actual == pytest.approx(expected, rel=0, abs=2e-10)
    assert not caught


@pytest.mark.parametrize("node,amplitude", [(1, 1.), (1, 1e308), (0, 4e307), (1, 1e-310)],
                         ids=["ordinary", "large_boundary", "large_center", "subnormal_tangent"])
def test_representable_packet_cotangent_does_not_overflow_its_resolution_check(node, amplitude):
    # At psi=0 the scalar kinetic block is exactly twice the nodal simplex
    # mass. Its center diagonal is V/5, and each boundary diagonal is V/20.
    coordinates, velocity = np.zeros(68), np.zeros(68)
    velocity[42+node] = amplitude
    mp = mpmath.mp.clone()
    mp.dps = 100
    volume = 10+10*mp.sqrt(5)/3
    expected = float(mp.mpf(amplitude)*volume/(5 if node == 0 else 20))
    assert np.isfinite(expected) and expected > 0
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = packet.phase_space(coordinates, velocity)
    assert np.isfinite(actual["momentum"]).all()
    assert actual["momentum"][30+node] == pytest.approx(expected, rel=2e-12, abs=0)
    assert actual["constant_moment_map"] == 0
    assert not caught, "a representable cotangent overflowed only the precision guard"


def test_unreportable_packet_cotangent_is_not_returned_as_infinity_and_nan():
    coordinates, velocity = np.zeros(68), np.zeros(68)
    velocity[42] = 6e307
    mp = mpmath.mp.clone()
    mp.dps = 100
    expected = mp.mpf(float(velocity[42]))*(10+10*mp.sqrt(5)/3)/5
    assert expected > np.finfo(float).max
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="precision|range|represent"):
            packet.phase_space(coordinates, velocity)


def test_packet_scale_uses_real_coordinates_without_overflowing_a_complex_modulus():
    coordinates, velocity = np.zeros(68), np.zeros(68)
    velocity[43] = velocity[56] = 1.6e308
    mp = mpmath.mp.clone()
    mp.dps = 100
    expected = float(mp.mpf(float(velocity[43]))*(10+10*mp.sqrt(5)/3)/20)
    assert np.isfinite(expected)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = packet.phase_space(coordinates, velocity)
    np.testing.assert_allclose(actual["momentum"][[31, 44]], expected, rtol=2e-12, atol=0)
    assert not caught

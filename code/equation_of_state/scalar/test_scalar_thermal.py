"""Independent oscillator partition controls and corrupt-receipt regressions."""
from copy import deepcopy
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from scalar_eos import thermal_point
import verify_scalar_eos as verifier


def partition_control(spectrum, scale, mass, temperature):
    """Differentiate log Z in dimensionless coordinates, from original inputs.

    No producer frequencies, occupations, stress components or support choices
    enter this control. A private context also checks global precision isolation.
    """
    ctx = mp.mp.clone()
    ctx.dps = 110
    scale, mass, temperature = map(ctx.mpf, (scale, mass, temperature))
    volume, beta = scale**3, 1 / temperature

    def log_partition(lam, b=beta, v=volume, m=mass):
        frequency = ctx.sqrt(m*m + lam / v**(ctx.mpf(2)/3))
        return -ctx.log1p(-ctx.exp(-b*frequency))

    modes = []
    for lam in map(ctx.mpf, spectrum):
        energy = -ctx.diff(lambda t: log_partition(lam, b=beta*ctx.exp(t)), 0)/beta
        pressure = temperature/volume * ctx.diff(
            lambda t: log_partition(lam, v=volume*ctx.exp(t)), 0)
        potential = -temperature/2 * ctx.diff(
            lambda t: log_partition(lam, m=mass*ctx.exp(t)), 0)
        frequency = ctx.sqrt(mass*mass + lam/scale**2)
        modes.append((frequency, energy/frequency, energy, energy/2,
                      3*pressure*volume/2, potential, pressure,
                      -temperature*log_partition(lam)))
    result = dict(zip(('frequencies', 'occupations', 'mode_thermal_energies',
                       'mode_K', 'mode_G', 'mode_U'),
                      (list(column) for column in zip(*modes))))
    energy = ctx.fsum(row[2] for row in modes)
    pressure = ctx.fsum(row[6] for row in modes)
    result.update(volume=volume, energy=energy, energy_density=energy/volume,
                  pressure=pressure, w=pressure*volume/energy,
                  helmholtz_free_energy=ctx.fsum(row[7] for row in modes))
    return result


def assert_partition_result(result, expected):
    for key, values in expected.items():
        for actual, reference in zip(np.atleast_1d(result[key]), np.atleast_1d(values)):
            if reference == 0:
                assert actual == 0, key
            else:
                assert np.isfinite(actual), key
                assert abs(mp.mpf(float(actual))/reference-1) < mp.mpf('2e-12'), key


@pytest.mark.parametrize('spectrum,scale,mass,temperature', [
    ([1., 4., 9.], 2., 3., 4.),
    ([1.], 1., 1e8, 1e8),  # K and U agree in binary64; pressure is positive.
    ([1.], 1., 0., .02),  # log(1-exp(-50)) must remain negative.
    ([1e200], 1e100, 0., 1.),  # no overflowing product before volume division
    ([1e-200], 1., 1e-100, 1e-100),  # E*lambda underflows before division
    ([0., 1., 4.], 1., 2., 3.),  # a massive zero-gradient mode is admissible
])
def test_thermal_components_follow_independent_partition(spectrum, scale, mass, temperature):
    expected = partition_control(spectrum, scale, mass, temperature)
    result = thermal_point(np.array(spectrum), scale, mass, temperature)
    assert_partition_result(result, expected)
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize('key', ['occupations', 'mode_thermal_energies', 'mode_K', 'mode_G'])
@pytest.mark.parametrize('replacement', [0., -1e-12])
def test_retained_thermal_verifier_must_not_erase_positive_modes(key, replacement):
    receipt = verifier.load(Path(__file__).with_name('scalar_eos_receipt.json'))
    changed = deepcopy(receipt)
    mode = min((value, j, i) for j, row in enumerate(changed['thermal_cases'])
               for i, value in enumerate(row[key]) if value > 0)
    _, case_index, mode_index = mode
    changed['thermal_cases'][case_index][key][mode_index] = replacement
    with mp.workdps(70), pytest.raises(ValueError, match=key):
        verifier.check_thermal(changed, verifier.tensor_spectrum())


@pytest.mark.parametrize('spectrum,scale,mass,temperature', [
    ([1.], 1., 0., -1.), ([1.], 1., 0., 0.),
    ([1.], 0., 0., 1.), ([1.], -1., 0., 1.),
    ([1.], 1., -1., 1.), ([-1.], 1., 2., 1.),
    ([0.], 1., 0., 1.), ([], 1., 0., 1.),
    ([[1.]], 1., 0., 1.), ([np.nan], 1., 0., 1.),
    ([1.], 1., 0., np.inf), ([1.], 1., True, 1.),
])
def test_invalid_thermal_inputs_are_refused(spectrum, scale, mass, temperature):
    with pytest.raises(ValueError):
        thermal_point(np.array(spectrum), scale, mass, temperature)

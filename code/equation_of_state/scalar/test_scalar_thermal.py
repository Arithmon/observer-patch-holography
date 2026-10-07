"""Independent oscillator partition controls and corrupt-receipt regressions."""
from copy import deepcopy
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from scalar_eos import thermal_point
import verify_scalar_eos as verifier


def partition_control(spectrum, scale, mass, temperature, *, digits=110):
    """Differentiate log Z in dimensionless coordinates, from original inputs.

    No producer frequencies, occupations, stress components or support choices
    enter this control. A private context also checks global precision isolation.
    """
    ctx = mp.mp.clone()
    ctx.dps = digits
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
    ctx = mp.mp.clone()
    ctx.dps = 110
    for key, values in expected.items():
        for actual, reference in zip(np.atleast_1d(result[key]), np.atleast_1d(values)):
            if reference == 0:
                assert actual == 0, key
            else:
                assert np.isfinite(actual), key
                assert abs(ctx.mpf(float(actual))/reference-1) < ctx.mpf('2e-12'), key


@pytest.mark.parametrize('spectrum,scale,mass,temperature', [
    ([1., 4., 9.], 2., 3., 4.),
    ([1.], 1., 1e8, 1e8),  # K and U agree in binary64; pressure is positive.
    ([1.], 1., 0., .02),  # log(1-exp(-50)) must remain negative.
    ([1e200], 1e100, 0., 1.),  # no overflowing product before volume division
    ([1e-240], 1., 1e-120, 1e-120),  # E*lambda underflows before division
    ([1e-300], 1e-100, 0., 1e-50),  # tiny volume with a finite density
    ([0., 1., 4.], 1., 2., 3.),  # a massive zero-gradient mode is admissible
])
def test_thermal_components_follow_independent_partition(spectrum, scale, mass, temperature):
    expected = partition_control(spectrum, scale, mass, temperature)
    result = thermal_point(np.array(spectrum), scale, mass, temperature)
    assert_partition_result(result, expected)
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize('key', ['occupations', 'mode_thermal_energies', 'mode_K', 'mode_G', 'mode_U'])
@pytest.mark.parametrize('replacement', [0., -1e-12, 'half'])
def test_retained_thermal_verifier_must_not_erase_positive_modes(key, replacement):
    receipt = verifier.load(Path(__file__).with_name('scalar_eos_receipt.json'))
    changed = deepcopy(receipt)
    mode = min((value, j, i) for j, row in enumerate(changed['thermal_cases'])
               for i, value in enumerate(row[key]) if value > 0)
    _, case_index, mode_index = mode
    changed['thermal_cases'][case_index][key][mode_index] = mode[0]/2 if replacement == 'half' else replacement
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


def test_massless_potential_is_exact_zero_in_retained_verification():
    receipt = verifier.load(Path(__file__).with_name('scalar_eos_receipt.json'))
    receipt['thermal_cases'][0]['mode_U'][0] = 1e-15
    with mp.workdps(70), pytest.raises(ValueError, match='mode_U exact zero'):
        verifier.check_thermal(receipt, verifier.tensor_spectrum())


def test_resolved_subnormal_components_remain_available():
    result = thermal_point([1.], 1., 0., 1/710)
    assert 0 < abs(result['helmholtz_free_energy']) < np.finfo(float).tiny
    assert_partition_result(result, partition_control([1.], 1., 0., 1/710))


@pytest.mark.parametrize('spectrum,scale,mass,temperature,field', [
    ([1.], 1., 0., 1/720, 'occupations|helmholtz_free_energy'),
    ([1.], 1., 0., 1e-300, 'occupations'),
    ([1.], 1e200, 1., 1., 'mode_G|volume'),
    ([1.], 1e-200, 0., 1e200, 'volume|energy_density|pressure'),
    ([0.], 1., 1e-308, 1e308, 'occupations'),
])
def test_unresolved_report_components_raise_instead_of_becoming_zero_or_infinite(
        spectrum, scale, mass, temperature, field):
    with pytest.raises(ValueError, match=field):
        thermal_point(spectrum, scale, mass, temperature)


def test_high_temperature_mode_with_subnormal_rest_mass():
    result = thermal_point([0.], 1., 1e-308, 1.)
    # x ~ 1e-308: Z ~ 1/x, E ~ T, F ~ T log(x), far below
    # any relevant correction to binary64 in these separately checked fields.
    assert result['frequencies'] == [1e-308]
    assert result['occupations'] == pytest.approx([1e308])
    assert result['energy'] == pytest.approx(1.)
    assert result['mode_K'] == pytest.approx([.5])
    assert result['mode_U'] == pytest.approx([.5])
    assert result['pressure'] == result['w'] == 0
    assert result['helmholtz_free_energy'] == pytest.approx(np.log(1e-308))


def test_resolved_pressure_far_below_the_rest_energy():
    # Resolve the volume dependence in the separate partition derivative,
    # rather than applying higher precision to the producer's rounded energy.
    expected = partition_control([1e-100], 1., 1e100, 1e100, digits=420)
    result = thermal_point([1e-100], 1., 1e100, 1e100)
    assert 0 < result['w'] < 1e-299
    assert_partition_result(result, expected)


@pytest.mark.parametrize('unit', [1e-120, 1., 1e120])
def test_energy_unit_change_preserves_each_thermal_component(unit):
    spectrum = [unit**2, 4*unit**2]
    result = thermal_point(spectrum, 2., 3*unit, 4*unit)
    expected = partition_control(spectrum, 2., 3*unit, 4*unit)
    assert_partition_result(result, expected)
    reference = thermal_point([1., 4.], 2., 3., 4.)
    for key in expected:
        factor = 1 if key in ('volume', 'occupations', 'w') else unit
        assert np.asarray(result[key])/factor == pytest.approx(reference[key], rel=2e-12, abs=0)


def test_thermal_evaluation_does_not_use_or_change_global_precision():
    reference = thermal_point([1.], 1., 0., .02)
    with mp.workdps(6):
        assert thermal_point([1.], 1., 0., .02) == reference
        assert mp.mp.dps == 6


@pytest.mark.parametrize('spectrum', [
    [1., True], np.array([True]), [1+0j], ['1'],
    np.ma.array([1., 4.], mask=[False, True]), [2**53+1],
])
def test_thermal_spectrum_does_not_coerce_invalid_or_lossy_data(spectrum):
    with pytest.raises(ValueError):
        thermal_point(spectrum, 1., 0., 1.)


def test_extended_precision_input_cannot_be_silently_rounded():
    if np.finfo(np.longdouble).nmant <= np.finfo(float).nmant:
        pytest.skip('platform has no precision beyond binary64')
    more_precise = np.longdouble(1) + np.finfo(np.longdouble).eps
    with pytest.raises(ValueError, match='loses information'):
        thermal_point([more_precise], 1., 0., 1.)


def test_thermal_mode_permutation_and_duplication():
    reference = thermal_point([1., 4., 9.], 1., 2., 3.)
    permuted = thermal_point([9., 1., 4.], 1., 2., 3.)
    doubled = thermal_point([1., 4., 9., 1., 4., 9.], 1., 2., 3.)
    for key in ('energy', 'energy_density', 'pressure', 'helmholtz_free_energy'):
        assert permuted[key] == reference[key]
        assert doubled[key] == 2*reference[key]
    assert doubled['w'] == reference['w']


def test_retained_spectrum_and_all_thermal_cases_agree_with_independent_verifier():
    receipt = verifier.load(Path(__file__).with_name('scalar_eos_receipt.json'))
    spectrum = receipt['spatial_operator']['generalized_eigenvalues']
    for case in receipt['thermal_cases']:
        new = thermal_point(spectrum, case['scale'], case['mass'], case['temperature'])
        assert set(new) == set(case)
        case.update(new)
    with mp.workdps(70):
        result = verifier.check_thermal(receipt, verifier.tensor_spectrum())
    assert result[1] < .06
    assert result[2] == pytest.approx(1/3)

"""Physical preparation, split isolation, alias and calibration checks."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest

from .data import load_observations
from . import executable_calibration as c
from ...observer.scalar_response.experiment import compare_receipt


def synthetic_axes(alpha, thermal, extent=4):
    # Independent formula expressed as amplitude times trigonometric phase.
    values = {}
    for i, j in c.AXIS_KEYS:
        u, v = i*extent/10, j*extent/10
        phase = 2*alpha*v
        values[i, j] = math.exp(-(thermal+.5)*(u*u+v*v))*complex(
            math.cos(phase), math.sin(phase))
    return values


@pytest.mark.parametrize('alpha,thermal', [(-1.1, .04), (.8, .2), (-2.4, .7)])
def test_calibration_recovers_independent_physical_preparations(alpha, thermal):
    fit = c.calibrate_state(synthetic_axes(alpha, thermal))
    assert fit['parameters']['alpha_real'] == pytest.approx(alpha, abs=2e-9)
    assert fit['parameters']['thermal_occupation'] == pytest.approx(thermal, abs=2e-9)
    assert fit['calibration_spin_mean_rmse'] < 1e-10
    assert fit['scope']['normalizable_positive_state']
    assert not fit['scope']['parameter_confidence_region_certified']


def test_only_axes_are_accepted_and_validation_mutation_cannot_change_fit():
    observations = load_observations()
    axis = observations['calibration']
    original = deepcopy(axis)
    fit = c.calibrate_state(axis)
    observations['validation'].update({key: 1-1j for key in observations['validation']})
    assert c.calibrate_state(axis) == fit
    assert axis == original
    with pytest.raises(ValueError, match='exactly the 21'):
        c.calibrate_state({**axis, **observations['validation']})
    assert fit['scope']['withheld_data_previously_examined_in_earlier_analysis']
    assert not fit['scope']['prospective_prediction']


def test_two_coordinate_conventions_are_separate_physical_calibrations():
    axis = load_observations()['calibration']
    primary = c.calibrate_state(axis, extent=4)
    alternate = c.calibrate_state(axis, extent=5)
    assert primary['parameters']['alpha_real'] == pytest.approx(.96341309, abs=2e-7)
    assert primary['parameters']['thermal_occupation'] == pytest.approx(.12777123, abs=2e-7)
    assert alternate['parameters']['alpha_real'] == pytest.approx(.76974579, abs=2e-7)
    assert alternate['parameters']['thermal_occupation'] == 0
    assert alternate['thermal_boundary_active']
    assert primary['calibration_input_sha256'] == alternate['calibration_input_sha256']
    assert not primary['scope']['absolute_coordinate_extent_identified']
    assert not alternate['scope']['absolute_coordinate_extent_identified']


def test_noisy_origin_does_not_renormalize_or_supply_visibility():
    axis = synthetic_axes(.7, .1)
    axis[0, 0] = .2+.7j
    fit = c.calibrate_state(axis)
    assert fit['parameters']['alpha_real'] == pytest.approx(.7, abs=2e-9)
    assert fit['parameters']['thermal_occupation'] == pytest.approx(.1, abs=2e-9)
    origin = next(row for row in fit['calibration_rows'] if row['index'] == [0, 0])
    assert origin['observed_real'] == .2
    assert origin['observed_imag'] == .7
    assert origin['fitted_real'] == 1
    assert origin['fitted_imag'] == 0
    assert fit['calibration_spin_mean_rmse'] == pytest.approx(math.sqrt((.8**2+.7**2)/42))


def test_phase_alias_and_unknown_direction_are_reported():
    extent = 4
    period = math.pi*10/extent
    axis = synthetic_axes(.8+period, .1)
    fit = c.calibrate_state(axis)
    assert fit['parameters']['alpha_real'] == pytest.approx(.8, abs=2e-9)
    assert fit['alpha_alias_period'] == pytest.approx(period)
    assert not fit['scope']['physical_coordinate_direction_identified']
    assert not fit['optimizer']['global_optimality_certified']


@pytest.mark.parametrize('extent', [0, -1, True, float('inf'), float('nan'), '4'])
def test_invalid_coordinate_extent_rejected(extent):
    with pytest.raises(ValueError, match='extent'):
        c.calibrate_state(synthetic_axes(.7, .1), extent=extent)


@pytest.mark.parametrize('value', [1.01, -1.01j, float('nan'), complex(0, float('inf')), True, '1'])
def test_invalid_spin_means_rejected(value):
    axis = synthetic_axes(.7, .1)
    axis[0, 0] = value
    with pytest.raises(ValueError, match='Spin means|spin component'):
        c.calibrate_state(axis)


def test_independent_components_may_exceed_unit_disk():
    axis = synthetic_axes(.7, .1)
    axis[0, 0] = .9+.9j
    fit = c.calibrate_state(axis)
    assert fit['parameters']['alpha_real'] == pytest.approx(.7, abs=2e-9)


@pytest.mark.parametrize('key', [(True, 0), (1., 0), (0, 0, 0), '0,0'])
def test_invalid_axis_keys_rejected(key):
    axis = synthetic_axes(.7, .1)
    axis.pop((0, 0))
    axis.pop((1, 0))
    axis[key] = 1
    with pytest.raises(ValueError, match='coordinates'):
        c.calibrate_state(axis)


def test_analytic_jacobian_fit_does_not_require_finite_difference_noise():
    fit = c.calibrate_state(synthetic_axes(-1.9, .09))
    singular = fit['optimizer']['jacobian_singular_values']
    assert len(singular) == 2 and min(singular) > 0
    assert fit['optimizer']['gradient_optimality_diagnostic'] < 1e-10
    n = fit['parameters']['thermal_occupation']
    ratio = n/(1+n)
    probabilities = (1-ratio)*ratio**np.arange(100)
    assert np.all(probabilities >= 0)
    assert np.sum(probabilities) == pytest.approx(1)


@pytest.mark.parametrize('extent', [4, 5])
def test_published_calibration_has_json_native_types_and_strict_receipt_roundtrip(extent):
    calibration = c.calibrate_state(load_observations()['calibration'], extent=extent)

    def audit(value):
        if type(value) is dict:
            assert all(type(key) is str for key in value)
            for item in value.values():
                audit(item)
        elif type(value) is list:
            for item in value:
                audit(item)
        else:
            assert type(value) in (str, float, int, bool, type(None))

    audit(calibration)
    decoded = json.loads(json.dumps(calibration, allow_nan=False))
    compare_receipt(decoded, calibration)
    assert decoded == calibration

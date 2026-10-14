"""Physical preparation, split isolation, alias and calibration checks."""
from copy import deepcopy
import json
import math
from types import SimpleNamespace

import mpmath as mp
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


def _refine(axis, seeds, extent=4):
    keys, values = c._axes(axis)
    coordinates = np.asarray(keys, dtype=float)*extent/10
    return c._refine_parameters(values, coordinates, seeds, math.pi*5/extent)


def test_platform_solver_seeds_and_start_diagnostics_do_not_define_preparation(monkeypatch):
    """Actual Linux/Windows basin seeds must declare the same binary source input."""
    axis = load_observations()['calibration']
    old_seed = [.9634130901023388, .12777123190717807]
    new_seed = [.9634130900677007, .12777123192307055]

    def calibrate(seed, successful_start):
        calls = iter(range(c.STARTS))
        predicted = synthetic_axes(*seed)
        cost = .5*sum(abs(predicted[key]-axis[key])**2 for key in axis)

        def solver(*args, **kwargs):
            selected = next(calls)
            return SimpleNamespace(
                x=np.asarray(seed), success=successful_start is None or selected == successful_start,
                cost=cost, optimality=1e-9 if successful_start is None else 2e-9)

        monkeypatch.setattr(c, 'least_squares', solver)
        return c.calibrate_state(axis)

    first, second = calibrate(old_seed, None), calibrate(new_seed, 17)
    assert first == second
    assert first['parameters'] == {
        'alpha_real': float.fromhex('0x1.ed447b0eb4eb9p-1'),
        'thermal_occupation': float.fromhex('0x1.05acec6756e4bp-3')}
    assert not first['optimizer']['solver_start_and_iteration_diagnostics_retained']
    assert not first['optimizer']['global_optimality_certified']


@pytest.mark.parametrize('seed,success,cost', [
    ([.8, .2], False, 0.),
    ([float('nan'), .2], True, 0.),
    ([.8, float('inf')], True, 0.),
    ([.8, .2], True, float('nan')),
    ([.8, .2], True, float('inf')),
    ([.8, .2], True, -1.),
    ([.8], True, 0.),
    ([.8, -.1], True, 0.),
    ([4., .2], True, 0.),
    ([True, False], True, 0.),
])
def test_failed_or_malformed_solver_results_cannot_supply_preparation(monkeypatch, seed, success, cost):
    monkeypatch.setattr(c, 'least_squares', lambda *args, **kwargs:
                        SimpleNamespace(x=np.asarray(seed), success=success, cost=cost))
    with pytest.raises(ArithmeticError, match='No axis-only'):
        c.calibrate_state(synthetic_axes(.8, .2))


@pytest.mark.parametrize('sign', [-1, 1])
@pytest.mark.parametrize('distance', [0., 1e-8])
@pytest.mark.parametrize('thermal', [0., .1])
def test_valid_alpha_boundary_and_neighbor_remain_supported(sign, distance, thermal):
    bound = math.pi*5/4
    alpha = sign*(bound-distance)
    fit = c.calibrate_state(synthetic_axes(alpha, thermal))
    actual = fit['parameters']
    assert -bound <= actual['alpha_real'] <= bound
    # The two interval endpoints are aliases; retain the original means as
    # the oracle rather than presupposing which endpoint wins its tiny tie.
    displacement = actual['alpha_real']-alpha
    assert min(abs(displacement), abs(displacement-2*bound),
               abs(displacement+2*bound)) < 1e-12
    assert actual['thermal_occupation'] == pytest.approx(thermal, abs=2e-15)
    assert fit['calibration_spin_mean_rmse'] < 1e-14
    assert fit['optimizer']['gradient_optimality_diagnostic'] < 1e-10


def test_thermal_boundary_sign_and_tiny_neighbor_are_checked_from_original_axes():
    axis = load_observations()['calibration']
    parameters, diagnostics = _refine(axis, [(.7697457903270853, 0.)], extent=5)
    assert parameters[1] == 0
    assert diagnostics['thermal_boundary_gradient'] > .6
    for thermal in (0., 5e-13):
        parameters, diagnostics = _refine(synthetic_axes(.8, thermal), [(.8, 0.)])
        assert parameters[1] == 0
        assert diagnostics['gradient_optimality_diagnostic'] < c.BOUNDARY_KKT_ALLOWANCE
    # A spurious optimizer boundary success must not erase an interior
    # occupation when the original-axis derivative resolves its direction.
    for thermal in (1e-8, .2):
        with pytest.raises(ArithmeticError, match='stationary refinement'):
            _refine(synthetic_axes(.8, thermal), [(.8, 0.)])
    parameters, _ = _refine(synthetic_axes(.8, 2e-12), [(.8, .01)])
    assert parameters[1] == pytest.approx(2e-12, abs=1e-16)


def test_original_input_objective_and_derivative_match_independent_mp_oracle():
    axis = load_observations()['calibration']
    parameters, diagnostics = _refine(axis, [(.9634130900677007, .12777123192307055)])
    ctx = mp.mp.clone()
    ctx.dps = 70

    def exact(value):
        numerator, denominator = float(value).as_integer_ratio()
        return ctx.mpf(numerator)/denominator

    def objective(alpha, thermal):
        total = ctx.mpf(0)
        for (i, j), observed in axis.items():
            u, v = exact(i*4/10), exact(j*4/10)
            predicted = ctx.exp(-(thermal+ctx.mpf('.5'))*(u*u+v*v)+2j*alpha*v)
            total += abs(predicted-ctx.mpc(exact(observed.real), exact(observed.imag)))**2
        return total

    alpha, thermal = map(exact, parameters)
    loss = objective(alpha, thermal)
    derivatives = [ctx.diff(lambda a: objective(a, thermal), alpha),
                   ctx.diff(lambda n: objective(alpha, n), thermal)]
    assert diagnostics['sum_of_squared_real_residuals'] == float(loss)
    assert diagnostics['gradient_optimality_diagnostic'] == pytest.approx(
        float(max(map(abs, derivatives))), rel=1e-12, abs=1e-28)


def test_refinement_isolates_mp_context_and_seed_order():
    axis = load_observations()['calibration']
    seeds = [(.9634130901023388, .12777123190717807),
             (.9634130900677007, .12777123192307055), (.95, .13)]
    expected, diagnostics = _refine(axis, seeds)
    with mp.workdps(7):
        before = (mp.mp.prec, mp.mp.dps, mp.mp.trap_complex)
        actual, changed_order = _refine(axis, list(reversed(seeds)))
        assert (mp.mp.prec, mp.mp.dps, mp.mp.trap_complex) == before
    assert np.array_equal(actual, expected)
    assert changed_order == diagnostics


def test_distinct_equal_cost_basins_use_canonical_parameter_tie_break():
    # Erasing the imaginary means makes the original-axis objective exactly
    # even in alpha. Its two nonzero stationary basins are equally admissible;
    # start order must not select a different declared source preparation.
    axis = {key: complex(value.real) for key, value in synthetic_axes(1.8, .1).items()}
    negative, negative_diagnostics = _refine(axis, [(-1.8, .1)])
    positive, positive_diagnostics = _refine(axis, [(1.8, .1)])
    assert negative[0] < 0 < positive[0]
    assert np.array_equal(negative, positive*np.array([-1., 1.]))
    assert negative_diagnostics == positive_diagnostics
    for seeds in ([(-1.8, .1), (1.8, .1)], [(1.8, .1), (-1.8, .1)]):
        selected, diagnostics = _refine(axis, seeds)
        assert np.array_equal(selected, negative)
        assert diagnostics == negative_diagnostics


def test_interior_refinement_success_requires_original_axis_stationarity(monkeypatch):
    context = mp.mp.clone()
    monkeypatch.setattr(c.mp.mp, 'clone', lambda: context)
    monkeypatch.setattr(context, 'findroot', lambda *args, **kwargs:
                        (context.mpf('.8'), context.mpf('.2')))
    with pytest.raises(ArithmeticError, match='stationary refinement'):
        _refine(synthetic_axes(.8, .1), [(.8, .2)])

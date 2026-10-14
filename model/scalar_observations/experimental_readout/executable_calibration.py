"""Axis-only calibration of one normalizable displaced thermal state.

The analytic characteristic function is used to calibrate two preparation
parameters, never to produce the executable detector's off-axis predictions.
The API deliberately accepts only the 21 calibration records, not a data grid
or a validation dictionary. This is a retrospective follow-up on already
examined public data; the split alone does not make it a prospective test.
"""
from collections.abc import Mapping
import hashlib
import json
import math
import numbers

import numpy as np
from scipy.optimize import least_squares


AXIS_KEYS = frozenset({(i, 0) for i in range(11)} | {(0, j) for j in range(11)})
STARTS = 33


def _axes(axis_means):
    if not isinstance(axis_means, Mapping):
        raise ValueError('Calibration requires the 21 axis means as a mapping')
    for key in axis_means:
        if (not isinstance(key, tuple) or len(key) != 2
                or any(type(index) is not int for index in key)):
            raise ValueError('Calibration coordinates must be pairs of integer indices')
    if set(axis_means) != AXIS_KEYS:
        raise ValueError('Calibration requires exactly the 21 nonnegative axis points')
    keys = sorted(AXIS_KEYS)
    values = []
    for key in keys:
        value = axis_means[key]
        if isinstance(value, bool) or not isinstance(value, numbers.Complex):
            raise ValueError('Spin means must be finite numeric complex values')
        value = complex(value)
        if (not math.isfinite(value.real) or not math.isfinite(value.imag)
                or abs(value.real) > 1 or abs(value.imag) > 1):
            raise ValueError('Each observed spin component must lie in [-1,1]')
        values.append(value)
    return keys, np.asarray(values, dtype=complex)


def _axis_model(parameters, coordinates):
    """Calibrating curve; coordinates are checked at the public API boundary."""
    alpha, occupation = parameters
    real, imaginary = coordinates.T
    return np.exp(-(occupation+.5)*(real*real+imaginary*imaginary)
                  + 2j*alpha*imaginary)


def calibrate_state(axis_means, *, extent=4.0):
    """Fit real displacement and nonnegative thermal occupation using axes only.

    Every complex mean contributes two equally weighted real residuals. There
    is no fitted gain, origin renormalization, offset or independent axis
    width. The physical state is D(alpha) rho_thermal(n) D(alpha)^dagger.

    Uniform imaginary-coordinate spacing aliases alpha with period pi/step;
    the numerical fit is restricted to the principal interval and is not an
    identification of the laboratory displacement or coordinate sign.
    """
    keys, values = _axes(axis_means)
    if (isinstance(extent, bool) or not isinstance(extent, numbers.Real)
            or not math.isfinite(extent) or extent <= 0):
        raise ValueError('Coordinate extent must be a positive finite real number')
    extent = float(extent)
    coordinates = np.asarray(keys, dtype=float)*extent/10
    radial = np.sum(coordinates**2, axis=1)
    imaginary = coordinates[:, 1]

    def residual(parameters):
        difference = _axis_model(parameters, coordinates)-values
        return np.concatenate((difference.real, difference.imag))

    def jacobian(parameters):
        mean = _axis_model(parameters, coordinates)
        derivatives = np.column_stack((2j*imaginary*mean, -radial*mean))
        return np.vstack((derivatives.real, derivatives.imag))

    # One full phase-alias interval. Starts are interior points, fixed without
    # looking at validation values; this is numerical multistart, not a
    # certificate of a global nonlinear optimum.
    alpha_bound = math.pi*5/extent
    fits = [least_squares(residual, [initial, .04], jac=jacobian,
                          bounds=([-alpha_bound, 0], [alpha_bound, np.inf]),
                          xtol=1e-13, ftol=1e-13, gtol=1e-13,
                          max_nfev=1000)
            for initial in np.linspace(-alpha_bound*.999, alpha_bound*.999, STARTS)]
    valid = [index for index, fit in enumerate(fits)
             if fit.success and np.isfinite(fit.x).all() and np.isfinite(fit.cost)]
    if not valid:
        raise ArithmeticError('No axis-only displaced-thermal fit converged')
    selected = min(valid, key=lambda index: fits[index].cost)
    fit = fits[selected]
    parameters = fit.x.copy()
    # A boundary optimum is exactly the vacuum thermal component. This tiny
    # rounding is disclosed and its axis residuals are recomputed below.
    if parameters[1] < 1e-12:
        parameters[1] = 0.
    prediction = _axis_model(parameters, coordinates)
    errors = residual(parameters)
    singular = np.linalg.svd(jacobian(parameters), compute_uv=False)
    calibration_rows = [
        {'index': list(key), 'real_coordinate': float(coordinate[0]),
         'imaginary_coordinate': float(coordinate[1]),
         'observed_real': float(observed.real), 'observed_imag': float(observed.imag),
         'fitted_real': float(fitted.real), 'fitted_imag': float(fitted.imag)}
        for key, coordinate, observed, fitted in zip(keys, coordinates, values, prediction)]
    bound_input = [{'index': list(key), 'real': float(value.real), 'imag': float(value.imag)}
                   for key, value in zip(keys, values)]
    return {
        'schema': 'oph.scalar_observations.axis_state_calibration.v1',
        'state_family': 'real_displaced_thermal_single_oscillator',
        'characteristic_convention': 'exp(-(n+1/2)*abs(beta)^2 + 2i*alpha*imag(beta))',
        'parameters': {'alpha_real': float(parameters[0]),
                       'thermal_occupation': float(parameters[1])},
        'coordinate_extent': extent,
        'alpha_alias_period': 2*alpha_bound,
        'alpha_fitted_principal_interval': [-alpha_bound, alpha_bound],
        'parameter_count': 2,
        'calibration_complex_points': len(keys),
        'calibration_real_statistics': 2*len(keys),
        'calibration_input_sha256': hashlib.sha256(
            json.dumps(bound_input, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        'calibration_rows': calibration_rows,
        'calibration_spin_mean_rmse': float(np.sqrt(np.mean(errors**2))),
        'calibration_max_component_residual': float(np.max(np.abs(errors))),
        'thermal_boundary_active': bool(parameters[1] == 0),
        'boundary_rounding_threshold': 1e-12,
        'optimizer': {
            'method': 'bounded least_squares with analytic Jacobian and fixed phase multistarts',
            'number_of_starts': STARTS,
            'successful_starts': len(valid),
            'selected_start': selected,
            'sum_of_squared_real_residuals': float(errors@errors),
            'gradient_optimality_diagnostic': float(fit.optimality),
            'jacobian_singular_values': [float(value) for value in singular],
            'global_optimality_certified': False},
        'scope': {
            'calibration_uses_withheld_values': False,
            'withheld_data_previously_examined_in_earlier_analysis': True,
            'prospective_prediction': False,
            'normalizable_positive_state': True,
            'arbitrary_axis_lookup_used_as_state': False,
            'origin_renormalization_or_fitted_visibility': False,
            'physical_coordinate_direction_identified': False,
            'absolute_coordinate_extent_identified': False,
            'thermal_parameter_is_independently_measured_temperature': False,
            'parameter_confidence_region_certified': False,
            'apparatus_error_budget_authenticated': False,
            'calibration_residual_is_sampling_confidence_interval': False},
    }

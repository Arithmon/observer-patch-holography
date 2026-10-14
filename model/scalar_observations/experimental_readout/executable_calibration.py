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

import mpmath as mp
import numpy as np
from scipy.optimize import least_squares


AXIS_KEYS = frozenset({(i, 0) for i in range(11)} | {(0, j) for j in range(11)})
STARTS = 33
REFINEMENT_DIGITS = 80
BOUNDARY_ROUNDING = 1e-12
BOUNDARY_KKT_ALLOWANCE = 1e-10
BASIN_OBJECTIVE_ALLOWANCE = 1e-12


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


def _refine_parameters(values, coordinates, seeds, alpha_bound):
    """Refine numerical basin seeds from original axis values, then select a point.

    This is numerical stationary-point refinement, not a global-optimum or
    parameter-confidence certificate. Only final binary64 parameters enter
    the source program and its separately exact occupation certificate.
    """
    ctx = mp.mp.clone()
    ctx.dps = REFINEMENT_DIGITS

    def exact(value):
        numerator, denominator = float(value).as_integer_ratio()
        return ctx.mpf(numerator)/denominator

    points = [(exact(u), exact(v), exact(z.real), exact(z.imag))
              for (u, v), z in zip(coordinates, values, strict=True)]

    def gradient(alpha, thermal):
        first = second = ctx.mpf(0)
        for u, v, real, imaginary in points:
            radius = u*u+v*v
            amplitude = ctx.exp(-(thermal+ctx.mpf('.5'))*radius)
            sine, cosine = ctx.sin(2*alpha*v), ctx.cos(2*alpha*v)
            first += 4*v*amplitude*(real*sine-imaginary*cosine)
            second += 2*radius*amplitude*(real*cosine+imaginary*sine-amplitude)
        return first, second

    def objective(alpha, thermal):
        result = ctx.mpf(0)
        for u, v, real, imaginary in points:
            amplitude = ctx.exp(-(thermal+ctx.mpf('.5'))*(u*u+v*v))
            result += (amplitude*ctx.cos(2*alpha*v)-real)**2
            result += (amplitude*ctx.sin(2*alpha*v)-imaginary)**2
        return result

    def boundary(alpha):
        step = min(ctx.mpf('.01'), exact(alpha_bound)/100)
        return ctx.findroot(lambda a: gradient(a, ctx.mpf(0))[0],
                            (alpha-step, alpha+step), tol=ctx.mpf('1e-60'), maxsteps=80)

    bound = exact(alpha_bound)

    candidates = {}
    for seed in seeds:
        alpha, thermal = map(exact, seed)
        try:
            if thermal < exact(BOUNDARY_ROUNDING):
                alpha, thermal = boundary(alpha), ctx.mpf(0)
            else:
                alpha, thermal = ctx.findroot(gradient, (alpha, thermal),
                                             tol=ctx.mpf('1e-60'), maxsteps=80)
            # An unconstrained stationary point can lie just outside a valid
            # box optimum. Refine the free coordinate on that active face;
            # both gradient signs are checked below, including its corners.
            alpha_face = abs(alpha) > bound
            thermal_face = thermal < exact(BOUNDARY_ROUNDING)
            if alpha_face:
                alpha = bound if alpha > 0 else -bound
                if thermal_face:
                    thermal = ctx.mpf(0)
                else:
                    thermal = ctx.findroot(lambda n: gradient(alpha, n)[1], thermal,
                                           tol=ctx.mpf('1e-60'), maxsteps=80)
                    if thermal < exact(BOUNDARY_ROUNDING):
                        thermal = ctx.mpf(0)
            elif thermal_face:
                alpha, thermal = boundary(alpha), ctx.mpf(0)
                if abs(alpha) > bound:
                    alpha = bound if alpha > 0 else -bound
            if (not ctx.isfinite(alpha) or not ctx.isfinite(thermal)
                    or abs(alpha) > bound or thermal < 0):
                continue
            da, dn = gradient(alpha, thermal)
            if alpha == bound:
                alpha_gradient = max(ctx.mpf(0), da)
            elif alpha == -bound:
                alpha_gradient = max(ctx.mpf(0), -da)
            else:
                alpha_gradient = abs(da)
            alpha_allowance = exact(BOUNDARY_KKT_ALLOWANCE) if abs(alpha) == bound else ctx.mpf('1e-40')
            if alpha_gradient > alpha_allowance:
                continue
            if thermal == 0:
                # A materially negative derivative points into the feasible
                # thermal domain, so that boundary is not an admissible fit.
                if dn < -exact(BOUNDARY_KKT_ALLOWANCE):
                    continue
            elif abs(dn) > ctx.mpf('1e-40'):
                continue
            point = (float(alpha), float(thermal))
            if not all(math.isfinite(value) for value in point):
                continue
            # Reevaluate at the actual declared binary point. The same point
            # has the same objective regardless of its optimizer start.
            a, n = map(exact, point)
            loss = objective(a, n)
            if loss > objective(*map(exact, seed))+ctx.mpf('1e-20'):
                continue
            da, dn = gradient(a, n)
            projected_alpha = max(ctx.mpf(0), da) if a == bound else (
                max(ctx.mpf(0), -da) if a == -bound else abs(da))
            projected = max(projected_alpha, max(ctx.mpf(0), -dn) if n == 0 else abs(dn))
            candidates[point] = (float(loss), float(projected), float(dn), float(da))
        except (ValueError, ZeroDivisionError, OverflowError):
            continue
    if not candidates:
        raise ArithmeticError('No admissible original-axis stationary refinement converged')
    selected = min(candidates, key=lambda point: (candidates[point][0], point))
    loss, projected, thermal_gradient, alpha_gradient = candidates[selected]
    return np.array(selected), {
        'sum_of_squared_real_residuals': loss,
        'gradient_optimality_diagnostic': projected,
        'thermal_boundary_gradient': thermal_gradient if selected[1] == 0 else None,
        'alpha_boundary_gradient': alpha_gradient if abs(selected[0]) == alpha_bound else None,
    }


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
    seeds = []
    for fit in fits:
        value = np.asarray(fit.x)
        if (not fit.success or value.shape != (2,) or value.dtype.kind not in 'iuf'
                or not np.isfinite(value).all() or not math.isfinite(fit.cost)
                or fit.cost < 0 or abs(value[0]) > alpha_bound or value[1] < 0):
            continue
        seeds.append(value)
    if not seeds:
        raise ArithmeticError('No axis-only displaced-thermal fit converged')
    scores = [float(residual(seed)@residual(seed)) for seed in seeds]
    best = min(scores)
    seeds = [seed for seed, score in zip(seeds, scores, strict=True)
             if score <= best+BASIN_OBJECTIVE_ALLOWANCE*max(1., best)]
    parameters, refinement = _refine_parameters(values, coordinates, seeds, alpha_bound)
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
        'boundary_rounding_threshold': BOUNDARY_ROUNDING,
        'optimizer': {
            'method': 'bounded least_squares basin seeds followed by original-input high-precision stationary refinement',
            'number_of_starts': STARTS,
            'refinement_decimal_digits': REFINEMENT_DIGITS,
            'selection_rule': 'least recomputed objective at distinct refined binary points; canonical parameter tuple breaks ties',
            'basin_objective_allowance': BASIN_OBJECTIVE_ALLOWANCE,
            'solver_start_and_iteration_diagnostics_retained': False,
            'boundary_kkt_gradient_allowance': BOUNDARY_KKT_ALLOWANCE,
            **refinement,
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

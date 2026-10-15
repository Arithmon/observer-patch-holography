"""Preparation fit accepting only the explicitly labelled zero-duration scan."""
from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.optimize import least_squares

from ..scalar_observations.experimental_readout.executable_calibration import (
    _refine_parameters,
)


def calibrate(probe_coordinates, spin_real, spin_imag):
    """Fit D(alpha) rho_thermal(nbar) D(alpha)^dagger on the t=0 P scan.

    Coordinates are literal Im(beta), not the displacement of the readout
    pulse (beta/2). There is no fitted visibility, offset or coordinate scale.
    This accepts one scan, not the time series or a validation dictionary.
    """
    arrays = [np.asarray(value) for value in (probe_coordinates, spin_real, spin_imag)]
    if any(a.ndim != 1 or a.dtype.kind not in 'iuf' or not np.isfinite(a).all()
           for a in arrays):
        raise ValueError('three finite real vectors are required')
    beta, real, imag = [np.asarray(a, dtype=float) for a in arrays]
    if (len(beta) < 3 or real.shape != beta.shape or imag.shape != beta.shape
            or not np.all(np.diff(beta) > 0) or beta[0] < 0
            or np.max(np.abs(real)) > 1 or np.max(np.abs(imag)) > 1):
        raise ValueError('ordered nonnegative coordinates and binary spin means required')
    observations = real + 1j*imag
    coordinates = np.column_stack((np.zeros(len(beta)), beta))

    def curve(parameters):
        alpha, thermal = parameters
        return np.exp(-(thermal+.5)*beta**2 + 2j*alpha*beta)

    def residual(parameters):
        difference = curve(parameters)-observations
        return np.r_[difference.real, difference.imag]

    def jacobian(parameters):
        values = curve(parameters)
        derivatives = np.column_stack((2j*beta*values, -beta**2*values))
        return np.vstack((derivatives.real, derivatives.imag))

    fits = [least_squares(residual, (start, .04), jac=jacobian,
                          bounds=([-2., 0.], [2., np.inf]),
                          ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=1000)
            for start in np.linspace(-1.998, 1.998, 33)]
    seeds = [fit.x for fit in fits if fit.success and np.isfinite(fit.x).all()]
    if not seeds:
        raise ArithmeticError('no preparation fit converged')
    objectives = [float(residual(seed)@residual(seed)) for seed in seeds]
    best = min(objectives)
    seeds = [seed for seed, loss in zip(seeds, objectives, strict=True)
             if loss <= best+1e-12*max(1., best)]
    parameters, refinement = _refine_parameters(observations, coordinates, seeds, 2.)
    fitted = curve(parameters)
    rows = [dict(beta_imag=float(b), observed_real=float(r), observed_imag=float(i),
                 fitted_real=float(z.real), fitted_imag=float(z.imag))
            for b, r, i, z in zip(beta, real, imag, fitted, strict=True)]
    input_rows = [[float(b), float(r), float(i)]
                  for b, r, i in zip(beta, real, imag, strict=True)]
    return {
        'parameters': {'alpha_real': float(parameters[0]),
                       'thermal_occupation': float(parameters[1])},
        'calibration_complex_points': len(beta), 'parameter_count': 2,
        'calibration_spin_mean_rmse': float(np.sqrt(np.mean(residual(parameters)**2))),
        'calibration_input_sha256': hashlib.sha256(json.dumps(
            input_rows, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
        'rows': rows, 'stationary_refinement': refinement,
        'method': '33 bounded least-squares starts; 80-digit original-input stationary refinement',
        'alpha_domain': [-2., 2.], 'thermal_domain': 'nonnegative',
        'zero_duration_only': True, 'visibility_offset_scale_fitted': False,
        'parameter_uncertainty_propagated': False,
        'global_optimality_certified': False,
    }

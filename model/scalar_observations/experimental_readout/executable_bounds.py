"""Rational occupation certificates for the executable oscillator readout.

The declared preparation caps a geometric thermal distribution at n=d-1,
then executes a finite-code displacement. The finite smooth pointer dilation
is compared with the full oscillator operation, not with compressed effects.
No fitted-parameter uncertainty, floating implementation error or apparatus
error is included in this mathematical occupation certificate.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache

from ...observer.source_scalar.bounds import root_upper
from ...observer.source_scalar.cutoff import exponential_upper


def _rational(value, name, *, nonnegative=False):
    if type(value) is not F or (nonnegative and value < 0):
        raise ValueError(f"{name} must be an exact Fraction" +
                         (" greater than or equal to zero" if nonnegative else ""))
    return value


def _dimension(value):
    if type(value) is not int or value < 2:
        raise ValueError("occupation dimension must be an integer of at least two")
    return value


def _digits(value):
    if type(value) is not int or value < 10:
        raise ValueError("directed precision must be an integer of at least ten digits")
    return value


@lru_cache(maxsize=4096)
def thermal_moment(weight, thermal, displacement_squared, *, digits=40):
    """Upper bound Tr[rho_displaced * weight**N], with directed exp.

    The exact expression is exp((q-1)*|z|**2/h)/h, h=1-nbar*(q-1).
    It follows by integrating the positive coherent-state representation of
    a thermal state; nbar=0 follows by the coherent-state Poisson series.
    """
    weight = _rational(weight, "exponential weight")
    thermal = _rational(thermal, "thermal occupation", nonnegative=True)
    displacement_squared = _rational(displacement_squared, "squared displacement",
                                    nonnegative=True)
    digits = _digits(digits)
    denominator = 1-thermal*(weight-1)
    if weight <= 1 or denominator <= 0:
        raise ValueError("weight must obey 1 < q and nbar*(q-1) < 1")
    exponent = (weight-1)*displacement_squared/denominator
    return exponential_upper(exponent, digits=digits)/denominator


def _weight_candidates(dimension, thermal, radius):
    # The derivative of log(M(q)/q**d) is strictly increasing. Its sign
    # selects two neighboring weights on a fixed rational mesh. Selection
    # improves tightness only: every individual candidate is a valid bound.
    ceiling = -(-thermal.numerator//thermal.denominator)
    denominator = 32*(1+ceiling)
    maximum = 512
    if thermal:
        exclusive = F(denominator, 1)/thermal
        maximum = min(maximum, -(-exclusive.numerator//exclusive.denominator)-1)

    def derivative(numerator):
        q = 1+F(numerator, denominator)
        h = 1-thermal*(q-1)
        return radius/(h*h)+thermal/h-F(dimension, 1)/q

    lower, upper = 1, maximum
    while lower < upper:
        middle = (lower+upper)//2
        if derivative(middle) < 0:
            lower = middle+1
        else:
            upper = middle
    indices = sorted({max(1, lower-1), lower})
    return tuple(1+F(index, denominator) for index in indices)


@lru_cache(maxsize=4096)
def thermal_tail(dimension, thermal, displacement_squared, *, digits=40):
    """Certified Chernoff bound for n>=d and separately for n=d.

    The same number bounds a supremum over all displacements of squared
    radius no greater than the supplied rational radius. It need not be a
    tight estimate of the actual photon-number distribution.
    """
    dimension = _dimension(dimension)
    thermal = _rational(thermal, "thermal occupation", nonnegative=True)
    displacement_squared = _rational(displacement_squared, "squared displacement",
                                    nonnegative=True)
    digits = _digits(digits)
    candidates = []
    for weight in _weight_candidates(dimension, thermal, displacement_squared):
        moment = thermal_moment(weight, thermal, displacement_squared, digits=digits)
        candidates.append((moment/weight**dimension, weight, moment))
    bound, weight, moment = min(candidates)
    # The undisplaced distribution is exactly geometric, including vacuum.
    # Keeping this additional bound avoids an artificial numerical floor.
    if displacement_squared == 0:
        bound = min(bound, (thermal/(1+thermal))**dimension)
    return {
        "exponential_weight": weight,
        "squared_displacement_upper": displacement_squared,
        "number_moment_upper": moment,
        "occupation_tail_upper": min(F(1), bound),
        "unclipped_chernoff_tail_upper": moment/weight**dimension,
        "undisplaced_geometric_bound_used": displacement_squared == 0,
    }


def cutoff_certificate(dimension, *, thermal, alpha_real, beta_real, beta_imag,
                       alpha_imag=F(0), digits=40):
    """Bound the actual capped-thermal, finite-displacement pointer readout.

    alpha and beta denote the displacement in D(alpha) and the detector's
    characteristic point D(beta). The latter has coherent pointer arms
    D(+beta/2) and D(-beta/2). Number rotations and pointer-only rotations
    preserve the projection and need no further occupation budget.

    All inputs are exact rational declared parameters. A caller may supply
    Fraction.from_float(value) to certify the exact binary value used in its
    supplied preparation; this does not certify the floating source pulses.
    """
    dimension = _dimension(dimension)
    digits = _digits(digits)
    thermal = _rational(thermal, "thermal occupation", nonnegative=True)
    alpha_real = _rational(alpha_real, "real preparation displacement")
    alpha_imag = _rational(alpha_imag, "imaginary preparation displacement")
    beta_real = _rational(beta_real, "real detector displacement")
    beta_imag = _rational(beta_imag, "imaginary detector displacement")
    preparation_radius = alpha_real**2+alpha_imag**2
    beta_squared = beta_real**2+beta_imag**2
    cross_term = abs(alpha_real*beta_real+alpha_imag*beta_imag)
    pointer_radius = preparation_radius+beta_squared/4+cross_term
    preparation = thermal_tail(dimension, thermal, preparation_radius, digits=digits)
    pointer = thermal_tail(dimension, thermal, pointer_radius, digits=digits)
    scale = 10**digits
    preparation_residual = root_upper(
        preparation_radius*dimension*preparation["occupation_tail_upper"], scale=scale)
    pointer_residual = root_upper(
        beta_squared*dimension*pointer["occupation_tail_upper"]/4, scale=scale)
    projected_error = preparation_residual+pointer_residual
    common_low_vector_error = root_upper(
        projected_error**2+pointer["occupation_tail_upper"], scale=scale)
    ratio = thermal/(1+thermal)
    capped_mass = ratio**dimension
    # The common initial n<d mixture has squared norm 1-r**d. Keeping
    # this factor sharpens the standard pure-state trace-distance bound.
    common_low_trace_error = root_upper(
        (1-capped_mass)*common_low_vector_error**2, scale=scale)
    trace_error = min(F(1), capped_mass+common_low_trace_error)
    # At beta=0 both components are pointer-only fixed identities, whatever
    # field preparation was supplied. Their probability error is exactly 0.
    probability_error = F(0) if beta_squared == 0 else trace_error
    return {
        "schema": "oph.experimental_readout.executable_cutoff.v1",
        "occupation_dimension": dimension,
        "thermal_occupation": thermal,
        "preparation_displacement_real": alpha_real,
        "preparation_displacement_imag": alpha_imag,
        "detector_displacement_real": beta_real,
        "detector_displacement_imag": beta_imag,
        "rounding_decimal_digits": digits,
        "thermal_geometric_ratio": ratio,
        "capped_thermal_probability_mass": capped_mass,
        "preparation_tail_certificate": preparation,
        "pointer_tail_certificate": pointer,
        "preparation_projected_residual_upper": preparation_residual,
        "pointer_projected_residual_upper": pointer_residual,
        "projected_common_low_path_error_upper": projected_error,
        "common_low_vector_error_upper": common_low_vector_error,
        "common_low_trace_error_upper": common_low_trace_error,
        "complete_instrument_trace_distance_upper": trace_error,
        "probability_error_upper": probability_error,
        "spin_mean_error_upper": 2*probability_error,
        "scope": {
            "capped_geometric_thermal_preparation": True,
            "finite_displacement_preparation_included": True,
            "finite_smooth_pointer_dilation_included": True,
            "both_pointer_outcomes_and_terminal_field_included": True,
            "unrestricted_displaced_thermal_reference": True,
            "preparation_postselection": False,
            "field_rotation_commutes_with_occupation_projection": True,
            "rational_directed_bound_arithmetic": True,
            "calibrated_parameter_uncertainty_included": False,
            "floating_thermal_weight_rounding_included": False,
            "numerical_source_execution_error_included": False,
            "physical_apparatus_error_included": False,
        },
    }

"""Infer local scalar coefficients from whole-field product-curve records.

The preparation and detector family is supplied; no isolated generator
service enters calibration.  Every onsite and bond term participates in
the same fixed product schedule at each duration.  The old inverse's linear
algebra is reused, with its derivative bound justified here by the sum of
centered factor-generator norms rather than one isolated Hamiltonian norm.
"""
from __future__ import annotations

from numbers import Integral

from ..interacting_reconstruction import recovery as inverse


def coefficient_from_records(probabilities, first_time, second_time, *,
                             product_generator_norm_sum_cap, probability_error=0.):
    """Estimate one derivative of a fixed real-generator product curve.

    Let U(t)=prod_l exp(-i t alpha_l H_l), with a fixed list of real symmetric
    H_l and real coefficients alpha_l.  Real product-number preparation and
    imaginary-coherence readout imply s(-t)=-s(t).  Its derivative at zero
    is the relevant matrix element of sum_l alpha_l H_l, after spectator
    contraction.  If sum_l |alpha_l| ||H_l-c_l I|| <= S, multinomial
    differentiation gives |s^(5)(t)| <= (2S)^5/2.  Thus the same Richardson
    radius as the earlier inverse applies with B replaced by S.

    The result identifies the derivative of the actual product curve; it
    does not approximate that curve by exp(-it sum_l alpha_l H_l).  There
    is no additional calibration splitting error.  Floating/source errors
    remain excluded unless covered by an explicit probability_error margin.
    """
    result = inverse.coefficient_from_records(probabilities, first_time, second_time,
        norm_bound=product_generator_norm_sum_cap, probability_error=probability_error)
    result['product_generator_norm_sum_cap'] = result.pop('generator_norm_bound')
    return result


def recover(records_first, records_second, clock_first, clock_second, *,
            omega, product_generator_norm_sum_cap, probability_error=0., edge_threshold=.1,
            require_positive_interaction=True):
    """Recover coefficients and conditional cube geometry without isolation.

    Input shapes match the earlier inverse.  Onsite probes now start from
    the full product blank, or from one addressed |1> with other sites blank.
    The local imaginary-coherence instruments marginalize every spectator.
    Pair probes start from the full blank and read joint |00>/<11| coherence,
    again marginalizing the rest.  Zero spectator position means eliminate
    external bonds from the first derivative only; higher-order paths remain
    in the records and are covered by the whole-schedule remainder bound.

    The supplied polynomial family, basis/preparations, common fixed product
    schedule, calibrated clock and zero-or-separated-positive topology family
    are still premises.  No general Hamiltonian tomography or continuum
    recovery follows from these two-time derivative estimates.
    """
    if not isinstance(records_first, dict) or not isinstance(records_second, dict):
        raise ValueError('two whole-field calibration record dictionaries are required')
    first_steps = records_first.get('fixed_symmetric_steps')
    second_steps = records_second.get('fixed_symmetric_steps')
    if (isinstance(first_steps, bool) or not isinstance(first_steps, Integral)
            or isinstance(second_steps, bool) or not isinstance(second_steps, Integral)
            or first_steps < 1 or first_steps != second_steps):
        raise ValueError('both records must declare the same positive fixed symmetric step count')
    if any(record.get('complete_field_schedule_in_every_trial') is not True
           or record.get('independently_isolated_action_services') is not False
           for record in (records_first, records_second)):
        raise ValueError('both records must declare whole-field calibration without isolated services')
    try:
        dimensions = tuple(records_first.get('occupation_dimensions', ()))
        second_dimensions = tuple(records_second.get('occupation_dimensions', ()))
    except TypeError as exc:
        raise ValueError('the same occupation dimensions are required at both durations') from exc
    if (len(dimensions) != 8 or dimensions != second_dimensions
            or any(isinstance(d, bool) or not isinstance(d, Integral) or d < 4 for d in dimensions)):
        raise ValueError('the same eight occupation dimensions, each at least four, are required at both durations')
    result = inverse.recover(records_first, records_second, clock_first, clock_second,
        omega=omega, norm_bound_onsite=product_generator_norm_sum_cap,
        norm_bound_bond=product_generator_norm_sum_cap, probability_error=probability_error,
        edge_threshold=edge_threshold, require_positive_interaction=require_positive_interaction)
    for row in result['onsite_coefficient_details']:
        for detail in row:
            detail['product_generator_norm_sum_cap'] = detail.pop('generator_norm_bound')
    for detail in result['bond_coefficient_details']:
        detail['product_generator_norm_sum_cap'] = detail.pop('generator_norm_bound')
    del result['isolated_generator_family_is_supplied']
    result.update({
        'whole_field_schedule_is_supplied': True,
        'fixed_symmetric_steps': int(first_steps),
        'calibration_occupation_dimensions': tuple(int(d) for d in dimensions),
        'single_generator_isolation_used': False,
        'spectator_outcomes_marginalized': True,
        'product_generator_norm_sum_cap': float(product_generator_norm_sum_cap),
        'calibration_curve': 'fixed real-generator product scaled by calibrated model duration',
        'calibration_bounds_compare_product_derivative_at_zero': True,
        'separate_calibration_splitting_error_required': False,
    })
    return result

"""Calibration bounds for a fixed product family of real source generators.

The sum of centered factor norms controls derivatives of the executed product
as its supplied model-coordinate duration varies.  A bound only on the sum
Hamiltonian is insufficient.  No semigroup or isolated-generator evolution
is assumed, and no separate splitting error is needed to identify the product
family's derivative at zero.  Clock, floating and physical service errors
remain separate from these ideal finite-dimensional certificates.
"""
from __future__ import annotations

from ..interacting_reconstruction.bounds import (
    _rational, _upper, richardson_radius,
)
from ..relativistic_scalar.bounds import generator_sum_bound


def product_coefficient_radius(first_time, second_time,
                               product_generator_norm_sum_cap=96,
                               probability_error=0):
    """Bound the derivative reconstructed from two local coherence laws.

    Let U(t)=product_l exp(-i t alpha_l H_l), with a factor sequence and real
    symmetric H_l fixed independently of t, and a real initial state.  For
    an imaginary coherence observable O of norm 1/2, s(t)=<U* O U> is odd.
    If B >= sum_l |alpha_l| ||H_l-c_l I||, the product rule gives
    ||U^(r)(t)|| <= B^r after removing a scalar phase.  Consequently
    |s^(5)(t)| <= (2B)^5/2 and the ordinary generalized Richardson remainder
    applies.  The first derivative depends on sum_l alpha_l H_l, although
    finite-duration U(t) need not equal that Hamiltonian's exponential.

    The displayed interval controls this derivative, not finite-time
    continuous evolution.  The same factor sequence and scaled coefficients
    must be used at both durations; changing the splitting step count is
    outside this argument.  A supplied per-law error delta propagates as in
    the isolated inverse, and does not automatically include clock error.
    """
    cap = _rational(product_generator_norm_sum_cap,
                    "centered product-generator norm-sum cap")
    result = richardson_radius(first_time, second_time, cap, probability_error)
    return {**result,
            "schedule_centered_norm_sum_bound": _upper(cap),
            "fixed_scaled_factor_sequence_required": True,
            "continuous_semigroup_assumed": False,
            "independent_splitting_error_added": False,
            "clock_floating_and_physical_service_error_included": False}


def product_generator_cap_certificate(dimensions=(4,)*8, coupling=1,
                                      declared_cap=96):
    """Later backend audit for the supplied canonical cube source schedule.

    Uses the existing exact default-action generator bound: centered onsite
    polynomials plus all twelve bond norm bounds.  In a full symmetric step,
    the two onsite half steps together have weight one and every bond has
    weight one.  Repeating n equally scaled steps leaves the norm sum
    unchanged.  This result is for the default m=c=a=1, Omega=2 action with
    coupling zero or one; it must not be fed back as fitted inverse data.
    """
    cap = _rational(declared_cap, "declared product norm-sum cap")
    actual_upper = generator_sum_bound(dimensions, coupling)
    return {
        "declared_schedule_centered_norm_sum_bound": _upper(cap),
        "exact_schedule_centered_norm_sum_upper": str(actual_upper),
        "schedule_centered_norm_sum_upper": _upper(actual_upper),
        "declared_cap_covers_source_schedule": actual_upper <= cap,
        "backend_audit_only": True,
        "supplied_action": "canonical cube; m=c=a=1, Omega=2; lambda zero or one",
        "global_hamiltonian_norm_alone_used": False,
        "floating_or_physical_service_error_included": False,
    }

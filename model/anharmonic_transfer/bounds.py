"""Ideal finite-occupation channel and selected-readout error certificates.

The certificate compares a symmetric dephasing/unitary split to the unsplit
finite Lindblad equation, with the same initial state and apparatus. It does
not certify occupation removal, inferred parameters or floating arithmetic.
All algebraic constants use exact rationals and an outward square root.
"""
from __future__ import annotations

from fractions import Fraction as F
import math
from numbers import Integral, Real

from ..relativistic_scalar.exact import sqrt_upper


def _integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f'{name} must be an integer at least {minimum}')
    return int(value)


def _rational(value, name, *, positive=False):
    if (isinstance(value, bool) or not isinstance(value, (Real, F))
            or not math.isfinite(value) or value < 0 or (positive and value == 0)):
        raise ValueError(f'{name} must be finite and nonnegative, positive where required')
    return F(value)


def _upper(value):
    value = F(value)
    result = float(value)
    if not math.isfinite(result):
        raise ValueError('bound is too large to represent')
    return math.nextafter(result, math.inf) if F(result) < value else result


def _field(name, value):
    return {name: _upper(value), f'exact_{name}': str(value)}


def strang_bound(dimension, steps, subdivisions, *, delta=2*math.pi*500,
                 omega=math.pi/150e-6, pulse_duration=50e-6,
                 dephasing_rate=18.):
    """Bound the complete drive schedule's half diamond distance.

    A=-i ad_H, D=-gamma ad_N^2/2, H=delta N-sign Omega X(a+a†)/2.
    One subpulse is exp(h D/2) exp(h A) exp(h D/2), h=tau/subdivisions.
    There are 4*steps*subdivisions such subpulses. Arbitrary exact CPTP
    interleavings, including spin rotations and commuting spin/free-number
    pulses, do not amplify this bound. Additional split noncommuting pulses
    need their own budgets. The finite [a,a†]=I-d|d-1><d-1| term is retained.

    A binary endpoint probability differs by at most the returned radius
    before selection. Use conditional_readout_bound after terminal selection.
    Float inputs denote their exact stored binary values, not parameter
    confidence intervals. Source-word arithmetic is outside this certificate.
    """
    d = _integer(dimension, 'occupation dimension', 2)
    k = _integer(steps, 'step count', 0)
    subdivisions = _integer(subdivisions, 'subdivisions', 1)
    delta = _rational(delta, 'harmonic frequency')
    omega = _rational(omega, 'drive frequency')
    tau = _rational(pulse_duration, 'pulse duration', positive=True)
    gamma = _rational(dephasing_rate, 'dephasing rate')
    n = d-1
    j = omega*sqrt_upper(F(n))
    # ||ad_L|| <= 2 delta*j + Omega²*d/2, L=-i[H,-i[H,N]].
    # Its last term is the actual finite-occupation boundary contribution.
    aa = gamma*((2*delta*j+omega*omega*d/2)*n+4*j*j)
    dd = 2*gamma*gamma*j*n*n
    h = tau/subdivisions
    count = 4*k*subdivisions
    per_pulse = h**3*(dd/24+aa/12)/2
    radius = min(F(1), count*per_pulse)
    return {
        'dimension': d, 'steps': k, 'subdivisions_per_pulse': subdivisions,
        'drive_subpulses': count,
        **_field('harmonic_frequency_radians_per_second', delta),
        **_field('drive_frequency_radians_per_second', omega),
        **_field('dephasing_rate_per_second', gamma),
        **_field('pulse_duration_seconds', tau),
        **_field('subpulse_duration_seconds', h),
        **_field('total_driven_duration_seconds', 4*k*tau),
        **_field('hamiltonian_double_commutator_diamond_upper', aa),
        **_field('dephasing_double_commutator_diamond_upper', dd),
        **_field('one_subpulse_half_diamond_upper', per_pulse),
        **_field('half_diamond_error_upper', radius),
        'finite_CCR_boundary_retained': True,
        'norm_is_uniform_over_inputs_and_reference_systems': True,
        'comparison': 'ideal symmetric split versus exact finite Lindblad drive schedule',
        'time_parameter': 'supplied simulated Schrodinger-frame duration in seconds',
        'scope': {
            'finite_occupation_dimension': d,
            'identical_initial_state_and_instruments': True,
            'quadrature_error_included': False,
            'occupation_removal_error_included': False,
            'floating_arithmetic_error_included': False,
            'physical_parameter_or_apparatus_error_included': False,
        },
    }


def channel_history_bound(applications):
    """Telescope (invocation_count, half_diamond_radius) channel contracts.

    Input radii may be exact Fractions reconstructed from a channel's exact
    certificate string. Arbitrary intervening CPTP instruments and retained
    records are permitted. Conditioning is performed only after this bound.
    """
    total = F(0)
    terms = []
    for count, radius in applications:
        count = _integer(count, 'channel invocation count', 0)
        radius = _rational(radius, 'channel half diamond radius')
        if radius > 1:
            raise ValueError('channel half diamond radius must not exceed one')
        total += count*radius
        terms.append({'invocations': count, **_field('per_invocation_radius', radius)})
    return {'terms': terms, **_field('half_diamond_error_upper', min(F(1), total)),
            'all_instrument_branches_retained': True,
            'conditioning_included': False}


def conditional_readout_bound(trace_distance_upper, selected_probability):
    """Bound every selected detector law using a retained selection mass.

    For two complete normalized states/laws at trace distance <=epsilon,
    conditional distance is <=epsilon/q, where q is either law's positive
    selection mass. Any positive lower bound for q may be supplied instead.
    Keeping the rejected branch gives this sharper bound without a factor
    two. The reference selection mass is guaranteed positive only if q>epsilon.
    The q<=epsilon case therefore reports a vacuous radius and no guarantee
    that a reference conditional law exists.

    A numerical acceptance value is interpreted as exact here; a caller with
    a certified arithmetic interval should supply its positive lower endpoint.
    """
    epsilon = _rational(trace_distance_upper, 'trace distance upper bound')
    q = _rational(selected_probability, 'selected probability', positive=True)
    if epsilon > 1 or q > 1:
        raise ValueError('trace distance and selected probability must not exceed one')
    probability = min(F(1), epsilon/q)
    return {
        **_field('selected_probability_or_lower_bound', q),
        **_field('unconditional_trace_distance_upper', epsilon),
        **_field('reference_selected_probability_lower', max(F(0), q-epsilon)),
        **_field('conditional_detector_tv_upper', probability),
        **_field('conditional_binary_probability_error_upper', probability),
        **_field('conditional_spin_mean_error_upper', 2*probability),
        'reference_selection_guaranteed_positive': q > epsilon,
        'both_selection_branches_retained': True,
        'acceptance_arithmetic_error_included': False,
    }


def execution_error_budget(dimension, steps, subdivisions, selected_probability,
                           *, environment_applications=(), **pulse_parameters):
    """Combine splitting and implemented-environment ideal channel errors."""
    splitting = strang_bound(dimension, steps, subdivisions, **pulse_parameters)
    environment = channel_history_bound(environment_applications)
    total = min(F(1), F(splitting['exact_half_diamond_error_upper'])
                +F(environment['exact_half_diamond_error_upper']))
    return {'splitting': splitting, 'environment': environment,
            **_field('complete_trace_distance_upper', total),
            'selected_readout': conditional_readout_bound(total, selected_probability),
            'occupation_removal_error_included': False,
            'floating_arithmetic_error_included': False,
            'physical_parameter_or_apparatus_error_included': False}

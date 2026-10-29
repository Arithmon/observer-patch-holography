"""Exact ideal-compilation certificates; no floating or hardware error claim."""
from fractions import Fraction as F
from math import factorial


def rational(x, name, *, positive=False):
    if type(x) not in (int, F) or x < 0 or (positive and not x):
        raise ValueError('nonnegative exact rational '+name+' required')
    return F(x)


def norm_sum(speed, interaction):
    return 6*rational(speed, 'speed')+rational(interaction, 'interaction')


def history(segments, speed=F(1, 8), interaction=F(1)):
    b = norm_sum(speed, interaction)
    total = F()
    for time, steps in segments:
        time = rational(time, 'model time')
        if type(steps) is not int or steps < 1:
            raise ValueError('positive integer steps required')
        total += b**3*time**3/(3*steps**2)
    return {'split_generator_norm_sum_upper': str(b),
            'complete_history_half_trace_upper': str(min(F(1), total)),
            'bound_includes_quantum_output_and_untouched_reference': True,
            'branch_count_multiplier': 1, 'uniform_in_box_volume': True,
            'floating_and_physical_errors_included': False}


def boundary(radius, total_time, speed=F(1, 8), interaction=F(1)):
    if type(radius) is not int or radius < 1:
        raise ValueError('positive radius required')
    t = rational(total_time, 'model time')
    b = norm_sum(speed, interaction)
    error = 2*(b*t)**radius/factorial(radius)
    return {'first_potentially_different_Taylor_order': radius,
            'complete_history_half_trace_upper': str(min(F(1), error)),
            'comparison': 'same fixed-spacing N=2 field on larger cubes or infinite Z^3',
            'apparatus_support': 'origin and its +x neighbor; no remote intervention',
            'spatial_continuum_claimed': False}


def interaction_witness(speed=F(1, 8), interaction=F(1), steps=16):
    c, u = rational(speed, 'speed'), rational(interaction, 'interaction', positive=True)
    b = norm_sum(c, u)
    time = F(1, 32)/b
    lower = 11*u*time/30
    split = F(history(((time, steps),), c, u)['complete_history_half_trace_upper'])
    free_split = F(history(((time, steps),), c, 0)['complete_history_half_trace_upper'])
    return {'model_time': str(time), 'first_derivative': str(u/2),
            'continuous_interaction_minus_free_probability_lower': str(lower),
            'ideal_executed_interaction_minus_free_probability_lower': str(lower-split-free_split),
            'all_box_radii_and_infinite_lattice': True,
            'observable': 'occupied down mode in (origin+i neighbor)/sqrt(2)',
            'preparation': 'origin up occupied; down=(origin+neighbor)/sqrt(2)',
            'proof': 'Duhamel commutator telescoping; 2Bt=1/16 and exp(1/16)-1<=1/15',
            'floating_or_hardware_error_certified': False}


def low_momentum(halfwidth, duration, speed=F(1, 8)):
    k = rational(halfwidth, 'coordinate momentum cutoff')
    t, c = rational(duration, 'model time'), rational(speed, 'speed')
    # Triangle inequality applied to three Pauli terms: no hidden numerical root.
    error = c*k**3/2
    return {'component_momentum_cutoff': str(k), 'one_body_operator_difference_upper': str(error),
            'free_two_particle_unitary_difference_upper': str(2*t*error),
            'interaction_projection_covered': False,
            'reason': 'local contact scattering can leave the selected momentum patch',
            'no_taste_discarded_in_execution': True}

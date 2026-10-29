"""Use the proved spectral enclosure to reduce a declared splitting budget.

This changes only the number of source evolution steps. The full N=2 action,
supplied global reference, native count decoder and complete detector remain.
"""
from fractions import Fraction as F
from math import isqrt

from ..band_clock import clock
from ..interacting_field import source
from ..local_repair_family import inference
from ..pair_cyclic_control.spectral_bound import canonical_clock_bound


def error_budget(parameter=F(-1), steps=21):
    # Reuse the original domain and statistical/splitting certificate.
    original = clock.error_budget(parameter, steps)
    improved = canonical_clock_bound(F(parameter))
    result = dict(original)
    result['centered_full_H_operator_norm_upper'] = str(improved['centered_norm_upper'])
    result['rate_error_halfdiamond_upper_on_confidence_event'] = str(
        improved['rate_error_halfdiamond_upper'])
    total = (F(original['ideal_split_halfdiamond_upper'])
             + improved['rate_error_halfdiamond_upper']
             + F(original['statistical_failure_upper']))
    result['complete_average_halfdiamond_upper'] = str(min(F(1), total))
    result['spectral_certificate'] = 'model/pair_cyclic_control/spectral_bound.py'
    return result


def minimum_steps(target=None, parameter=F(-1)):
    """Smallest positive integer satisfying this exact bound, not a global optimum.

    The old 32-step uniform ceiling is the default. No empirical output or
    target particle mass enters the choice. Other splitting methods, tighter
    estimates and shorter actual durations can give different optima.
    """
    clock.error_budget(parameter, 1)
    if target is None:
        target = F(clock.error_budget()['complete_average_halfdiamond_upper'])
    if type(target) not in (int, F) or not 0 < target < 1:
        raise ValueError('exact rational target strictly between zero and one required')
    target = F(target)
    fixed = (canonical_clock_bound(F(parameter))['rate_error_halfdiamond_upper']
             + inference.fixed_budget()['statistical_failure_upper'])
    coefficient = F(7, 4)**3*(F(5, 4)*abs(F(parameter)))**3/3
    allowance = target-fixed
    if allowance < 0 or (allowance == 0 and coefficient):
        raise ValueError('target is below the retained rate/statistical floor')
    if coefficient == 0:
        return 1
    ratio = coefficient/allowance
    steps = max(1, isqrt(ratio.numerator//ratio.denominator))
    if steps*steps*ratio.denominator < ratio.numerator:
        steps += 1
    return steps


def execute(counts, state, parameter=F(-1), space=None):
    """Execute the reduced schedule with every original stop/read output."""
    steps = minimum_steps(parameter=parameter)
    return clock.execute(counts, state, parameter, steps=steps, space=space)


def resources(space, duration, steps):
    """Full evolution and both complete occupation-pointer services."""
    result = source.resources(space, duration, steps, segments=1)
    _, pointer = source.occupation()
    for key in source.COUNT_KEYS:
        result[key] += 2*pointer['resources'][key]
    # The old floating residual diagnostic is evolution-only at this point;
    # expose counted resources only, never label it a hardware certificate.
    return {key: result[key] for key in (*source.COUNT_KEYS, 'routing_FSWAP_calls',
        'completed_hopping_blocks', 'completed_contact_blocks')} | {
        'fresh_pointer_blanks': 2, 'committed_pointer_bits': 2,
        'field_mode_register_size': space.modes,
        'scope': 'evolution plus two pointers; preparation and count acquisition separate',
        'physical_duration_or_energy_claimed': False}

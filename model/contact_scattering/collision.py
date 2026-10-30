"""Complete source-operated finite-time interaction echo and occupation read.

This is a controlled local interaction experiment, not an asymptotic S matrix.
Initial separated fermions, blank pointer and backward free evolution are
supplied services. Every first-failure export and both read outcomes remain.
"""
from __future__ import annotations

from fractions import Fraction as F
import math
import numpy as np
from scipy.linalg import expm

from ..field_thermal.thermal import execute_with_exports, failure_probability
from ..full_fock_field.field import FieldGeometry, FullFock, source_word
from ..full_fock_field.local_simulation import split_bound
from ..observer.source_scalar import compiler


def geometry(axis=0, interaction=1.):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('spatial axis 0,1,2 required')
    sides = [1, 1, 1]; sides[axis] = 2
    return FieldGeometry(tuple(sides), speed=.125, interaction=interaction)


def preparation():
    """Two separated fermions entangled with an untouched reference qubit."""
    value = np.zeros((16, 2), complex)
    value[9, 0] = 1/math.sqrt(2)       # modes0,3
    value[5, 1] = 1j/math.sqrt(2)      # modes0,2
    return value


def exact_echo(axis=0, duration=4., interaction=1., initial=None):
    g = geometry(axis, interaction); free = geometry(axis, 0.)
    space = FullFock(4)
    initial = preparation() if initial is None else space.state(initial)
    if type(duration) not in (int, float) or not math.isfinite(duration):
        raise ValueError('finite model duration required')
    return (expm(1j*duration*space.hamiltonian(free).toarray())@
            expm(-1j*duration*space.hamiltonian(g).toarray())@initial)


def execute(axis=0, *, duration=F(4), steps=32, interaction=1., initial=None):
    if type(duration) not in (int, F) or type(steps) is not int or not 1 <= steps <= 64:
        raise ValueError('exact rational duration and split count1..64 required')
    if abs(duration) > 8:
        raise ValueError('bounded source echo requires |model duration|<=8')
    g = geometry(axis, interaction); free = geometry(axis, 0.)
    space = FullFock(4)
    state = preparation() if initial is None else space.state(initial)
    if state.ndim == 1:
        state = state[:, None]
    if state.shape[1] > 2 or abs(float(np.vdot(state, state).real)-1) > 1e-12:
        raise ValueError('normalized input with at most two reference columns required')
    word = (*source_word(g, float(duration), steps),
            *source_word(free, -float(duration), steps),
            *compiler.cnot(0, 4))
    initial_pointer = np.zeros((32, state.shape[1]), complex)
    initial_pointer[:16] = state
    success, failures = execute_with_exports(word, 5, initial_pointer)
    target = exact_echo(axis, float(duration), interaction, state)
    expected = (target*((space.addresses%2 == 0)[:, None]),
                target*((space.addresses%2 == 1)[:, None]))
    outputs = (success[:16], success[16:])
    failure_mass = failure_probability(failures)
    vector_error = math.sqrt(sum(float(np.linalg.norm(a-b)**2)
                                 for a, b in zip(outputs, expected, strict=True))+failure_mass)
    ideal_bound = (split_bound(g.sides, F(duration), steps, interaction=F(str(g.interaction)), refined=True)+
                   split_bound(free.sides, F(duration), steps, interaction=F(0), refined=True))
    probability = [float(np.vdot(x, x).real) for x in outputs]
    ref = success.conj().T@success
    for event in failures:
        for branch in event:
            ref += branch.conj().T@branch
    gram_error = float(np.linalg.norm(ref-state.conj().T@state))
    if (not np.isfinite([failure_mass, vector_error, gram_error, *probability]).all() or
            gram_error > 2e-8 or vector_error > float(ideal_bound)+2e-8):
        raise ArithmeticError('complete source interaction echo failed admission')
    return {'outputs':outputs, 'first_failure_events':failures,
            'report':{
                'axis':axis, 'model_duration':str(duration), 'steps':steps,
                'interaction':g.interaction, 'source_resources':compiler.resources(word),
                'initialized_field_bits':4, 'initialized_pointer_bits':1,
                'supplied_untouched_reference_dimension':state.shape[1],
                'pointer_probabilities':probability,
                'failure_probability':failure_mass, 'reference_gram_residual':gram_error,
                'complete_vector_error':vector_error,
                'ideal_word_half_diamond_bound':str(ideal_bound),
                'floating_admission_tolerance':2e-8,
                'exact_empty_probability':float(np.linalg.norm(expected[0])**2),
                'initial_empty_probability':float(np.linalg.norm(state[space.addresses%2 == 0])**2),
                'all_outputs_retained':True,
                'initialization_and_physical_clock_supplied':True,
                'asymptotic_scattering_or_physical_cross_section_claimed':False}}


def collision_certificate():
    """Exact rational enclosure for the canonical z-bond experiment.

    Its active two-level block is [[0,ic],[-ic,U]]. For c=1/8,U=1,t=4,
    p_empty=(1-Re alpha)/4 with the expression proved in README.md.
    Taylor remainders and the square-root bracket are rational, not fitted
    to the numerical source output.
    """
    scale=10**12
    lower=F(math.isqrt(17*scale*scale),scale)
    upper=lower+F(1,scale)
    if not lower*lower <= 17 < upper*upper or lower < 4:
        raise ArithmeticError('square-root enclosure failed')
    def trig(x, odd):
        return sum(((-1)**((n-odd)//2)*x**n/math.factorial(n)
                    for n in range(odd,40,2)),F())
    arguments=(F(1,2),F(2),lower/2)
    remainder=max(abs(x)**40/math.factorial(40) for x in arguments)
    if remainder >= F(1,100):
        raise ArithmeticError('Taylor domain failed')
    ch,ct,cq=(trig(x,0) for x in arguments)
    sh,st,sq=(trig(x,1) for x in arguments)
    probability=(1-ct*(ch*cq+sh*sq/lower)-st*4*ch*sq/lower)/4
    error=100*remainder+(upper-lower)
    raw_lower,raw_upper=probability-error,probability+error
    display_scale=10**6
    lo=F((raw_lower*display_scale).__floor__(),display_scale)
    hi=F((raw_upper*display_scale).__ceil__(),display_scale)
    g=geometry(2)
    split=split_bound(g.sides,F(4),32,refined=True)
    if not F(46,1000)<lo<hi<F(48,1000) or lo-split <= F(36,1000):
        raise ArithmeticError('interaction signal enclosure failed')
    return {'axis':2,'model_duration':'4','speed':'1/8','interaction':'1',
            'steps':32,'exact_empty_probability_interval':[str(lo),str(hi)],
            'ideal_word_half_diamond_bound':str(split),
            'ideal_free_empty_probability':'0',
            'ideal_word_interaction_signal_lower':str(lo-split),
            'preparation_supplied':True,'floating_and_hardware_errors_included':False}

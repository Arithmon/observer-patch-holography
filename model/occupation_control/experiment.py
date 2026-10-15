"""End-to-end validation of the supplied continuous apparatus law.

This evaluator generates polynomial paths and checks their full defects; it
does not execute native observer words or read empirical response values.
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np

from . import intervals as iv
from .solver import _validate_state_values, evolve

PARAMETERS = {
    'alpha_real': -1.0081292326480422,
    'thermal': .13314125172459518,
    'delta': 2*math.pi*500,
    'omega': math.pi/150e-6,
    'pulse_duration': 50e-6,
    'theta': .8,
    'gamma': 18.,
    'omega_carrier': math.pi/35e-6,
    'probe_h': .4,
}
CASES = {'a': (0., 'left'), 'b': (math.pi/20, 'left'),
         'c': (math.pi/20, 'right'), 'd': (math.pi/20, 'center'),
         'e': (math.pi/10, 'left')}


def thermal_state(dimension, thermal):
    """Capped geometric mixture, with exact thermal-tail and rounding bounds."""
    if (type(dimension) is not int or dimension < 2 or isinstance(thermal, bool)
            or not math.isfinite(thermal) or thermal < 0):
        raise ValueError('valid dimension and nonnegative finite thermal input required')
    thermal = Fraction(thermal)
    ratio = thermal/(1+thermal)
    weights = [(1-ratio)*ratio**n for n in range(dimension-1)]+[ratio**(dimension-1)]
    rounded = list(map(float, weights))
    rho = np.zeros((2*dimension, 2*dimension), complex)
    rho[np.arange(dimension), np.arange(dimension)] = rounded
    rounding = sum(abs(Fraction(x)-y) for x, y in zip(rounded, weights, strict=True))/2
    return rho, {'thermal_tail': iv.rational(ratio**dimension)[1],
                 'thermal_rounding': iv.rational(rounding)[1]}


def _add(budgets):
    return {key: iv.sum_up([budget.get(key, 0.) for budget in budgets])
            for key in sorted({key for budget in budgets for key in budget})}


def _total(budget):
    return iv.sum_up(list(budget.values()))


def detector(rho, h, *, cells=8, order=12):
    """Pointer initially |+>, opposite half displacements, then exact -Y read.

    The full two-spin-branch binary detector law is retained. There is no
    postselection in propagation; selection is only the final ratio of records.
    """
    if isinstance(h, bool) or not math.isfinite(h) or h <= 0:
        raise ValueError('positive finite literal detector coordinate required')
    _validate_state_values(rho)
    joint = np.block([[rho*.5, rho*.5], [rho*.5, rho*.5]])
    target_half = iv.complex_scale(iv.complex_exact(rho), iv.exact(.5))
    block_error = iv.complex_sub(iv.complex_exact(rho*.5), target_half)
    embedding_error = float(iv.up(2*iv.trace_upper(*block_error)))
    # D(ih)=exp(ihQ). Pointer off-diagonal trace is chi(ih)/2.
    joint, budget = evolve(joint, cells=cells, order=order, sectors=4,
                           force=-h/2, pointer_sign=True)
    budget['pointer_embedding_rounding'] = embedding_error
    d = len(rho)//2
    masses, signals = [], []
    for spin in (0, 1):
        first = range(spin*d, (spin+1)*d)
        second = range((2+spin)*d, (3+spin)*d)
        masses.append(sum((Fraction(joint[i, i].real)+Fraction(joint[j, j].real)
                           for i, j in zip(first, second, strict=True)), Fraction(0)))
        signals.append(2*sum((Fraction(joint[i, j].imag)
                             for i, j in zip(first, second, strict=True)), Fraction(0)))
    exact_reads = [[(q+s)/2, (q-s)/2] for q, s in zip(masses, signals, strict=True)]
    if masses[0] <= 0 or any(value < 0 for row in exact_reads for value in row):
        raise ArithmeticError('approximate detector does not give nonnegative records')
    reads = [list(map(float, row)) for row in exact_reads]
    mean = signals[0]/masses[0]
    mean_float = float(mean)
    estimator_float = mean_float/(math.sqrt(2)*h)
    denominator = iv.mul(iv.square_root(iv.exact(2.)), iv.exact(h))
    estimator_enclosure = sorted([mean/Fraction(float(endpoint)) for endpoint in denominator])
    estimator_rounding = max(abs(Fraction(estimator_float)-endpoint) for endpoint in estimator_enclosure)
    conditional = [value/masses[0] for value in exact_reads[0]]
    conditional_rounding = max(abs(Fraction(float(value))-value) for value in conditional)
    return {'spin_detector_joint_probabilities': reads,
            'acceptance': float(masses[0]),
            'acceptance_exact_binary_sum': str(masses[0]),
            'selected_detector_probabilities': list(map(float, conditional)),
            'selected_mean': mean_float,
            'selected_estimator': estimator_float,
            'selected_probability_evaluation_error_upper': iv.rational(conditional_rounding)[1],
            'selected_mean_evaluation_error_upper': iv.rational(abs(Fraction(mean_float)-mean))[1],
            'estimator_evaluation_error_upper': iv.rational(estimator_rounding)[1],
            'mass_drift': float(abs(sum(masses)-1))}, budget


def run(*, dimension=32, case='a', steps=50, pulse_cells=5, order=12,
        record_steps=None):
    if (case not in CASES or type(steps) is not int or steps < 0
            or type(pulse_cells) is not int or pulse_cells < 1):
        raise ValueError('known case and nonnegative integer steps required')
    if record_steps is None:
        record_steps = [0, steps] if steps else [0]
    if (not isinstance(record_steps, (list, tuple)) or not record_steps
            or list(record_steps) != sorted(set(record_steps)) or
            any(type(s) is not int or not 0 <= s <= steps for s in record_steps)):
        raise ValueError('sorted unique record steps within horizon required')
    p = PARAMETERS
    phi, side = CASES[case]
    alpha = {'left': p['alpha_real'], 'right': -p['alpha_real'], 'center': 0.}[side]
    state, initial = thermal_state(dimension, p['thermal'])
    state, preparation = evolve(state, cells=8, order=order, force=alpha, force_imaginary=True)
    prior = [initial, preparation]
    evolution = []
    records = []
    # Preserve the source's binary gap-duration value as the declared premise.
    gap = abs(phi)/p['omega_carrier']
    for step in range(steps+1):
        if step in record_steps:
            read, read_budget = detector(state, p['probe_h'], order=order)
            budget = _add(prior+evolution+[read_budget])
            epsilon = _total(budget)
            finite_error = _total({key: value for key, value in budget.items()
                                   if key not in ('thermal_tail', 'occupation_defect')})
            cutoff_only = iv.sum_up([epsilon, finite_error])
            q = Fraction(read['acceptance_exact_binary_sum'])
            # The approximate path is Hermitian but not guaranteed positive or
            # trace one. The general conditional estimate uses 2 epsilon/q.
            conditional = estimator_radius = mean_radius = None
            if q > 2*Fraction(epsilon):
                conditional = min(1., iv.rational(2*Fraction(epsilon)/q+
                                  Fraction(read['selected_probability_evaluation_error_upper']))[1])
                mean_radius = iv.rational(4*Fraction(epsilon)/q+
                                         Fraction(read['selected_mean_evaluation_error_upper']))[1]
                denominator_lower = iv.mul(iv.square_root(iv.exact(2.)), iv.exact(p['probe_h']))[0]
                estimator_radius = iv.rational(4*Fraction(epsilon)/q/Fraction(float(denominator_lower))+
                                              Fraction(read['estimator_evaluation_error_upper']))[1]
            records.append({'steps': step, 'published_coordinate_ms': step/5,
                            'drive_elapsed_seconds': step*4*(p['pulse_duration']+gap),
                            **read, 'half_trace_error_upper': epsilon,
                            'exact_finite_path_half_trace_error_upper': finite_error,
                            'unrestricted_vs_exact_finite_half_trace_error_upper': cutoff_only,
                            'selected_probability_error_upper': conditional,
                            'selected_mean_error_upper': mean_radius,
                            'estimator_error_upper': estimator_radius,
                            'acceptance_positive_certified': q > 2*Fraction(epsilon),
                            'budget': budget,
                            'preparation_budget': _add(prior),
                            'evolution_budget': _add(evolution),
                            'detector_budget': read_budget})
        if step == steps:
            break
        for sign, multiple, spin_angles in ((-1, 1, (-phi, phi)),
                                            (1, 2, (phi, -phi)),
                                            (-1, 1, None)):
            state, budget = evolve(state, cells=pulse_cells*multiple, order=order,
                                   duration=Fraction(p['pulse_duration'])*multiple,
                                   delta=p['delta'], force=-sign*p['omega']/2,
                                   force_spin_flip=True, gamma=p['gamma'])
            evolution.append(budget)
            if spin_angles is not None:
                for index, angle in enumerate(spin_angles):
                    if angle:
                        # A normalized control parameter integrates the finite
                        # carrier interval; oscillator/noise rates include gap.
                        state, budget = evolve(state, cells=2, order=order,
                                               delta=Fraction(p['delta'])*Fraction(gap),
                                               gamma=Fraction(p['gamma'])*Fraction(gap),
                                               spin_x=angle/2)
                        evolution.append(budget)
                    if index == 0:
                        state, budget = evolve(state, cells=2, order=order, spin_z=p['theta']/2)
                        evolution.append(budget)
    return {'schema': 'oph.occupation_control.v1', 'dimension': dimension,
            'case': case, 'steps': steps, 'pulse_cells': pulse_cells, 'order': order,
            'parameters': {**p, 'phi': phi, 'alpha_real': alpha, 'gap_duration': gap},
            'records': records,
            'certificate': {'state_error_includes_occupation': True,
                            'state_error_includes_time_discretization': True,
                            'state_error_includes_evaluator_roundoff': True,
                            'preparation_and_detector_displacements_included': True,
                            'native_source_words_executed': False,
                            'native_compiler_error_included': False,
                            'apparatus_or_calibration_uncertainty_included': False,
                            'empirical_outcomes_loaded': False,
                            'arithmetic_assumptions': 'IEEE-754 binary64 elementary operations and sqrt, round-to-nearest, gradual underflow',
                            'target': 'declared exact binary parameters in the conservative unrestricted spin-oscillator Lindblad apparatus',
                            'method': 'piecewise polynomial full-defect integral, outward intervals, CPTP trace-norm contraction'}}

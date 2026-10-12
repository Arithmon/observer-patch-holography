"""A readable finite modular oscillator on the actual interacting history.

The unconditional terminal d=2 scalar state is rho=(m I-delta X)/2,
where m is its retained total weight (one up to numerical execution error).
The two eigenvalues p_minus,p_plus are obtained from that state, never fitted
to a desired clock rate. Its modular Hamiltonian is alpha I+beta X with
2 beta=log(p_minus/p_plus).

An emitted Z(+pi/4) pulse perturbs the existing branch ensemble. The emitted
X(+beta*s) pulse implements the Schroedinger dual of sigma_s(A)=rho^(is) A
rho^(-is). A source CNOT pointer then reads the scalar number bit. Its full
law is P(1;s)=(m+delta sin(2 beta*s))/2. Each setting extends a fresh copy of
the same complete history; no outcome, hidden working state, or scalar
density is substituted as a fresh preparation. All old records survive.

This is a probability-based clock on a declared monotone phase window, not
a single-shot time reading. Its period aliases, the unperturbed reference
is stationary, and the source pulse schedule and analogue angles are supplied.
Neither model-time nor physical-time calibration follows from this operation.
"""
from __future__ import annotations

import math
from numbers import Real

import numpy as np

from ..observer.source_scalar import compiler as c
from ..observer.source_scalar.experiment import encoded
from . import experiment as e
from . import operational as o


_X = np.array([[0., 1.], [1., 0.]])
_NUMERICAL_TOLERANCE = 2e-10


def _finite(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f'{name} must be a finite real number')
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number')
    return value


def parameters(result):
    """Extract this history's clock; reject a trivial or incompatible state."""
    if type(result['dimension']) is not int or result['dimension'] != 2:
        raise ValueError('this emitted clock uses the actual d=2 scalar code')
    rho, mass = o.scalar_state(result)
    o.density_spectrum(rho)
    delta = -float(np.trace(rho@_X).real)
    discrepancy = float(np.linalg.norm(rho-(mass*np.eye(2)-delta*_X)/2, 2))
    if discrepancy > _NUMERICAL_TOLERANCE:
        raise ValueError('terminal state is not the declared coordinate-sign readout state')
    p_minus, p_plus = (mass+delta)/2, (mass-delta)/2
    if min(p_minus, p_plus) <= o.RANK_TOLERANCE:
        raise ValueError('this full-algebra modular clock requires a faithful scalar state')
    if abs(delta) <= _NUMERICAL_TOLERANCE:
        raise ValueError('the actual scalar state has no nontrivial modular clock contrast')
    signed_gap = math.log(p_minus/p_plus)
    frequency, contrast = abs(signed_gap), abs(delta)
    return {'mass': mass, 'p_minus': p_minus, 'p_plus': p_plus,
            'signed_contrast': delta, 'contrast': contrast,
            'signed_modular_gap': signed_gap, 'modular_frequency': frequency,
            'state_form_operator_error': discrepancy,
            'monotone_parameter_window': [0., math.pi/(2*frequency)],
            'bounded_slope_parameter_window': [0., math.pi/(3*frequency)],
            'minimum_slope_on_bounded_window': contrast*frequency/4,
            'inverse_lipschitz_on_bounded_window': 4/(contrast*frequency),
            'period_in_modular_parameter': 2*math.pi/frequency}


def inferred_parameter(probability, clock):
    """Invert a population frequency on the explicitly chosen first window.

    This does not turn one observed bit into an exact time. Finite sampling,
    state estimation and hardware errors need their own probability bounds.
    """
    probability = _finite(probability, 'outcome probability')
    value = (2*probability-clock['mass'])/clock['contrast']
    if value < -1e-10 or value > 1+1e-10:
        raise ValueError('frequency lies outside the chosen monotone clock window')
    # Only endpoint roundoff is clipped, never a state or outcome weight.
    return math.asin(min(1., max(0., value)))/clock['modular_frequency']


def sample_setting(result, parameter, prepared=True):
    """Execute one complete setting on all four retained quantum branches."""
    parameter = _finite(parameter, 'modular parameter')
    if type(prepared) is not bool:
        raise ValueError('the preparation control must be Boolean')
    clock = parameters(result)
    vectors = o.terminal_vectors(result)
    word = [c.pulse('Z', math.pi/4 if prepared else 0., [2]),
            c.pulse('X', clock['signed_modular_gap']*parameter/2, [2])]
    unitary, guard = c.stabilize_executed_unitary(c.execute(word, 4))
    pointer, pointer_meta = e.pointer(2, 2)
    kraus = tuple(branch@unitary for branch in pointer)
    completeness_error = float(np.linalg.norm(
        sum(branch.conj().T@branch for branch in kraus)-np.eye(16), 2))
    rows, marginal = [], [0., 0.]
    old_marginal_error = 0.
    for record, vector in sorted(vectors.items()):
        subtotal = 0.
        for bit, branch in enumerate(kraus):
            terminal = branch@vector
            probability = float(np.vdot(terminal, terminal).real)
            rows.append({'old_record': list(record), 'clock_bit': bit,
                         'probability': probability,
                         'unnormalized_terminal_state': encoded(terminal)})
            marginal[bit] += probability
            subtotal += probability
        old_marginal_error = max(old_marginal_error,
                                float(abs(subtotal-np.vdot(vector, vector).real)))
    predicted = (clock['mass']+clock['signed_contrast']*math.sin(
        clock['signed_modular_gap']*parameter))/2 if prepared else clock['mass']/2
    probability_error = abs(marginal[1]-predicted)
    if max(probability_error, old_marginal_error, completeness_error) > _NUMERICAL_TOLERANCE:
        raise AssertionError('source clock word or complete retained instrument disagrees')
    source_counts = c.resources(word)
    appended = {key: source_counts[key]+pointer_meta['resources'][key]
                for key in e.COUNT_KEYS}
    inherited = result['resources_per_realized_history']
    total = {key: inherited[key]+appended[key] for key in e.COUNT_KEYS}
    return {'modular_parameter': parameter,
            'phase': clock['modular_frequency']*parameter,
            'probe_prepared': prepared,
            'clock_bit_marginal': marginal,
            'analytic_number_one_weight': predicted,
            'analytic_probability_error': probability_error,
            'old_record_marginal_error': old_marginal_error,
            'full_input_kraus_completeness_error': completeness_error,
            'branches': rows,
            'execution': {'probe_and_modular_word': word, 'unitary_guard': guard,
                          'clock_pointer': pointer_meta},
            'appended_source_resources_per_realized_history': {
                **appended, 'system_qubits': 4, 'pointer_qubits_per_read': 1,
                'fresh_pointer_blanks': 1, 'new_committed_classical_bits': 1},
            'total_source_resources_per_realized_history': {
                **total, 'system_qubits': 4, 'pointer_qubits_per_read': 1,
                'fresh_pointer_blanks': inherited['fresh_pointer_blanks']+1,
                'committed_classical_bits': inherited['committed_classical_bits']+1},
            'all_outcomes_retained': True}


def receipt(result):
    """A same-history finite clock instrument and its negative controls."""
    clock = parameters(result)
    phases = (0., math.pi/12, math.pi/6, math.pi/4, math.pi/3)
    settings = [sample_setting(result, phase/clock['modular_frequency']) for phase in phases]
    inversion_errors = [abs(inferred_parameter(row['clock_bit_marginal'][1], clock)
                            -row['modular_parameter']) for row in settings]
    stationary = [sample_setting(result, row['modular_parameter'], prepared=False)
                  for row in settings]
    selected = settings[2]
    wrapped = sample_setting(result, selected['modular_parameter']
                             +clock['period_in_modular_parameter'])
    mirrored = sample_setting(result, (math.pi-selected['phase'])
                              /clock['modular_frequency'])
    stationary_error = max(abs(row['clock_bit_marginal'][1]-clock['mass']/2)
                           for row in stationary)
    wrapped_error = abs(wrapped['clock_bit_marginal'][1]-selected['clock_bit_marginal'][1])
    mirrored_error = abs(mirrored['clock_bit_marginal'][1]-selected['clock_bit_marginal'][1])
    probabilities = [row['clock_bit_marginal'][1] for row in settings]
    checks = {'same_history_source_clock_law': all(
                  row['analytic_probability_error'] < _NUMERICAL_TOLERANCE for row in settings),
              'strictly_increasing_on_declared_window': all(
                  x < y for x, y in zip(probabilities, probabilities[1:])),
              'frequency_inverts_on_chosen_window': max(inversion_errors) < 1e-8,
              'unperturbed_reference_does_not_tick': stationary_error < _NUMERICAL_TOLERANCE,
              'period_wrap_is_indistinguishable': wrapped_error < _NUMERICAL_TOLERANCE,
              'outside_window_mirror_is_indistinguishable': mirrored_error < _NUMERICAL_TOLERANCE,
              'every_setting_preserves_old_records_and_all_outcomes': all(
                  len(row['branches']) == 8 and row['old_record_marginal_error'] < _NUMERICAL_TOLERANCE
                  for row in settings+stationary+[wrapped, mirrored])}
    if not all(checks.values()):
        raise AssertionError(checks)
    return {'schema': 'oph.common_source.modular_clock.v1',
            'dimension': 2, 'steps_per_segment': result['steps_per_segment'],
            'native_current': result['native_current'],
            'feedback_enabled': result['feedback_enabled'],
            'state_access': 'unconditional terminal site-0 scalar algebra; old records retained by the experiment',
            'clock_parameters': clock,
            'convention': 'sigma_s(A)=rho^(is) A rho^(-is); probe unitary exp(+is K), K=-log(rho), with scalar global phase omitted',
            'probe': 'emitted Z(+pi/4) on each actual retained working branch, followed by emitted modular X rotation and number-bit pointer',
            'clock_law': 'P(clock_bit=1;s)=(mass+delta*sin(s*log(p_minus/p_plus)))/2',
            'settings': settings, 'unprepared_stationary_control': stationary,
            'wrapped_parameter_control': wrapped,
            'outside_window_mirror_control': mirrored,
            'maximum_parameter_inversion_error': max(inversion_errors),
            'checks': checks, 'all_pass': all(checks.values()),
            'scope': {
                'one_binary_outcome_is_an_exact_time_reading': False,
                'same_system_sequential_use_of_different_settings': False,
                'fresh_complete_history_for_each_setting': True,
                'clock_rate_and_contrast_fitted_to_target': False,
                'source_pulse_schedule_and_analogue_angles_supplied': True,
                'physical_time_calibration_established': False,
                'original_hamiltonian_model_time_identification_established': False,
                'noise_or_reference_state_smoothing_added': False,
                'numerical_interval_certificate': False,
                'state_estimation_sampling_service_and_spatial_errors_certified': False,
                'meaning': 'a readable finite modular-phase instrument; no autonomous tick counter, universal clock, BW cap identification or continuum physics claim'}}

"""Execute the exact native completion attachment in the common CQ runtime.

The experiment reads a terminal binary pointer after every hidden native
repair record has been discarded. Both pointer outcomes and their complete
quantum/reference outputs remain. This is a native-channel semantic execution,
not emitted source pulses, an apparatus identification or an empirical fit.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np

from . import algebra, completion, runtime


TOLERANCE = 2e-12
EVENT_COUNTS = (1, 2, 4, 8)


def _bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or value < 0 or value > TOLERANCE:
        raise ArithmeticError(f'{name} failed the complete native-channel comparison')
    return value


def _preparation():
    a, u, _ = completion.witness_indices()
    matrix = np.zeros((completion.DIMENSION, completion.DIMENSION), complex)
    matrix[np.ix_((a, u), (a, u))] = .5
    return matrix


def _read():
    plus = _preparation()
    return runtime.Instrument(
        kraus=((plus,), (np.eye(completion.DIMENSION)-plus,)), visible=True,
        binding={'instrument': 'complete_Luders_binary_read',
                 'plus_effect': 'projector onto (|a>+|u>)/sqrt(2)',
                 'minus_effect': 'full identity minus plus projector',
                 'preparation_and_read_interface_supplied': True})


def _program(events, mode):
    if type(events) is not int or events not in EVENT_COUNTS:
        raise ValueError('one of the four declared native event counts is required')
    if mode not in ('iterated', 'compressed', 'completion'):
        raise ValueError('a declared completion comparison mode is required')
    maps = (completion.compressed_kraus(events) if mode == 'compressed' else
            completion.binary_kraus('E' if mode == 'completion' else 'C'))
    hidden = runtime.Instrument(kraus=(maps,), binding={
        'kind': 'native_binary_channel_semantics', 'mode': mode,
        'fixed_seam': [0, 1], 'all_native_records_discarded': True,
        'emitted_native_pulse_words': False})
    return runtime.Program(
        f'fixed_seam_{mode}_{events}', completion.DIMENSION, _preparation(),
        (runtime.Event('hidden_fixed_seam_repair', sequence=(hidden,),
                       repeat=events if mode == 'iterated' else 1),
         runtime.Event('complete_binary_read', instruments=(_read(),))),
        binding={'reference': 'I66/66 for the entropy-completion theorem',
                 'preparation': 'specified coherent pure witness, not the tracial reference',
                 'represented_C_events': events if mode != 'completion' else None,
                 'comparison_to_E': mode == 'completion',
                 'parameter': 'fixed-seam event count, not physical time',
                 'native_pulse_execution': False})


def _analytic_terminal(events):
    """Independent dyadic coefficients of the prepared C**n density."""
    a, u, v = completion.witness_indices()
    rho = np.zeros((completion.DIMENSION, completion.DIMENSION), complex)
    rho[a, a], rho[u, u], rho[v, v] = .5, .25, .25
    rho[a, u] = rho[u, a] = float(F(1, 2**(events+1)))
    return rho


def _terminal_error(history, rho):
    projections = (_preparation(), np.eye(completion.DIMENSION)-_preparation())
    expected = tuple(p@rho@p for p in projections)
    if history.records != ((0,), (1,)):
        raise ArithmeticError('the complete binary record dictionary changed')
    return _bounded(sum(float(np.linalg.norm(branch.density-target, 'fro'))
                        for branch, target in zip(history.branches, expected, strict=True)),
                    'independent terminal density')


def _public_plus(history):
    # Exercise the same public observer algebra as the other joined programs.
    value = algebra.public_expectation(algebra.observer_view(history), {(0,): 1., (1,): 0.})
    if not (math.isfinite(value.real) and math.isfinite(value.imag)) or abs(value.imag) > TOLERANCE:
        raise ArithmeticError('public binary read is not a finite real probability')
    return float(value.real)


def _reference_comparison():
    a, u, _ = completion.witness_indices()
    r = 2
    p, q = r*a, r*u+1
    initial = np.zeros((r*completion.DIMENSION,)*2, complex)
    initial[p, p] = initial[q, q] = .5
    initial[p, q], initial[q, p] = -.5j, .5j
    base = _program(4, 'iterated')
    program = runtime.Program(
        'reference_fixed_seam_checkpoint', completion.DIMENSION, base.initial_density,
        (runtime.Event('first_two_hidden_repairs', sequence=base.events[0].sequence, repeat=2),
         runtime.Event('later_two_hidden_repairs', sequence=base.events[0].sequence, repeat=2),
         base.events[1]), binding=base.binding)
    full = runtime.execute(program, initial=initial)
    compressed = runtime.execute(_program(4, 'compressed'), initial=initial)
    prefix = runtime.execute(program, initial=initial, stop_after=1)
    serialized = runtime.checkpoint_to_json(prefix)
    restored = runtime.checkpoint_from_json(program, serialized)
    continued = runtime.execute(program, checkpoint=restored)
    compression_error = _bounded(runtime.complete_distance(full, compressed), 'reference compression')
    continuation_error = _bounded(runtime.complete_distance(full, continued), 'reference checkpoint')
    # The untouched reference's marginal is computed from complete quantum
    # output, not reconstructed from its two visible pointer probabilities.
    marginal = sum(np.trace(branch.density.reshape(66, r, 66, r), axis1=0, axis2=2)
                   for branch in full.branches)
    reference_error = _bounded(np.linalg.norm(marginal-np.eye(r)/2, 'fro'),
                               'untouched reference marginal')
    if full.reference_dimension != r or any(branch.density.shape != (132, 132) for branch in full.branches):
        raise ArithmeticError('the complete reference carrier was reduced')
    return {'native_events': 4, 'reference_dimension': r,
            'preparation': '(|a,0>+i|u,1>)/sqrt(2)', 'visible_records': 2,
            'checkpoint_after_native_events': 2, 'native_events_after_checkpoint': 2,
            'quantum_output_dimension_per_record': 132,
            'compressed_complete_cq_distance_diagnostic': compression_error,
            'checkpoint_complete_cq_distance_diagnostic': continuation_error,
            'reference_marginal_residual_diagnostic': reference_error,
            'checkpoint_includes_full_quantum_state': True,
            'checkpoint_json_independently_replayed': True}


def run():
    """Verify exact mathematics and independently execute every retained law."""
    exact = completion.audit()
    baseline = runtime.execute(_program(1, 'completion'))
    baseline_plus = _public_plus(baseline)
    _bounded(abs(baseline_plus-3/8), 'completion probability')
    rows = []
    for events in EVENT_COUNTS:
        program = _program(events, 'iterated')
        direct = runtime.execute(program)
        compressed = runtime.execute(_program(events, 'compressed'))
        expected_distance = completion.half_diamond_distance(events)
        expected_plus = F(3, 8)+expected_distance
        plus = _public_plus(direct)
        compressed_plus = _public_plus(compressed)
        distance = runtime.complete_distance(direct, baseline)
        diagnostics = {
            'iterated_plus_probability': plus,
            'compressed_plus_probability': compressed_plus,
            'completion_plus_probability': baseline_plus,
            'probability_residual': _bounded(max(abs(plus-float(expected_plus)),
                                               abs(compressed_plus-float(expected_plus))),
                                            'analytic probability'),
            'complete_cq_distance_to_E': distance,
            'distance_residual': _bounded(abs(distance-float(expected_distance)),
                                         'saturated complete-output distance'),
            'compressed_complete_cq_distance': _bounded(runtime.complete_distance(direct, compressed),
                                                        'complete compression'),
            'independent_terminal_density_residual': _terminal_error(direct, _analytic_terminal(events)),
            'complete_probability_mass_residual': _bounded(max(abs(direct.probability_mass-1),
                                                               abs(compressed.probability_mass-1)),
                                                            'complete probability mass'),
        }
        if any(not math.isfinite(value) for value in diagnostics.values()):
            raise ArithmeticError('nonfinite completion execution diagnostic')
        rows.append({'native_events': events, 'complete_records': [[0], [1]],
                     'ideal_C_power_probabilities': [str(expected_plus), str(1-expected_plus)],
                     'ideal_E_probabilities': ['3/8', '5/8'],
                     'exact_total_variation_and_half_diamond': str(expected_distance),
                     'iterated_hidden_channel_evaluations': events,
                     'compressed_hidden_channel_evaluations': 1,
                     'diagnostics': diagnostics})
    return {'schema': 'oph.primitive_observer.executed_completion.v1', 'all_pass': True,
            'exact_attachment': exact, 'runs': rows,
            'entangled_reference_and_checkpoint': _reference_comparison(),
            'scope': {'executor': 'shared finite classical/quantum runtime',
                      'native_rule_semantics_executed': True,
                      'native_pulse_words_executed': False,
                      'all_intermediate_native_records_forgotten': True,
                      'both_terminal_read_outcomes_retained': True,
                      'physical_event_duration_reduced': False,
                      'same_completion_does_not_select_finite_step_law': True,
                      'preparation_read_and_binary_domain_supplied': True,
                      'empirical_data_received_or_fitted': False,
                      'roundoff_or_physical_service_certificate': False}}

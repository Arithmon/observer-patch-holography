"""Independent endpoint replay of the joint quantum clock/response experiment.

This module does not use the producer's operator blocks, resolvent, transformed
cone or Kraus matrices. It builds C directly from native endpoint and current
equality. Fractions preserve the complete unnormalized source/reference state.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache

from ..observer import native, quantum_mixing

STATES = native.states(binary=True)
INDEX = {state: i for i, state in enumerate(STATES)}
REGIONS = ((0, 1), (0, 5))
UNION = (0, 1, 5)


def _clean(values):
    return {key: value for key, value in values.items() if value}


def _add(target, source):
    for key, value in source.items():
        target[key] = target.get(key, F(0)) + value


def _mark(seam, current):
    return (tuple(int(any(port in region for port in seam)) for region in REGIONS),
            tuple(current * (int(port == seam[1]) - int(port == seam[0]))
                  for port in UNION))


@lru_cache(maxsize=None)
def endpoint_table(seam):
    """Two independent native endpoint calls per basis state, no matrix import."""
    return tuple(tuple((INDEX[output], record[1]) for output, record in
                       (native.native(state, seam, coin) for state in STATES))
                 for coin in (0, 1))


@lru_cache(maxsize=None)
def mark_menu(seam=None):
    edges = (seam,) if seam is not None else quantum_mixing.seams()
    marks = set()
    for edge in edges:
        incident = any(port in UNION for port in edge)
        for table in endpoint_table(edge):
            for _, current in table:
                marks.add(_mark(edge, current) if incident else None)
    return tuple(sorted(marks, key=repr))


def native_step(rho, seam=None):
    """One C attempt, retaining the declared incident mark and forgetting its coin.

    A C Kraus operator groups only equal native currents at a fixed seam and
    coin. Different current records cannot interfere. Reference indices r,s
    are copied unchanged; off-diagonal reference blocks are never traced away.
    """
    edges = (seam,) if seam is not None else quantum_mixing.seams()
    denominator = 2 * len(edges)
    out = {mark: defaultdict(F) for mark in mark_menu(seam)}
    for edge in edges:
        incident = any(port in UNION for port in edge)
        for table in endpoint_table(edge):
            for (i, j, r, s), value in rho.items():
                left, jl = table[i]
                right, jr = table[j]
                if jl == jr:
                    mark = _mark(edge, jl) if incident else None
                    out[mark][left, right, r, s] += value / denominator
    return {mark: _clean(matrix) for mark, matrix in out.items()}


def interval(rho, cap=2):
    if type(cap) is not int or cap < 1:
        raise ValueError('positive integer attempt cap required')
    out = {mark: {} for mark in mark_menu() if mark is not None}
    waiting = dict(rho)
    for _ in range(cap):
        step = native_step(waiting)
        for mark in out:
            _add(out[mark], step[mark])
        waiting = step[None]
    return {mark: _clean(matrix) for mark, matrix in out.items()}, waiting


def preparation(reference=False, dephased=False):
    if type(reference) is not bool or type(dephased) is not bool:
        raise ValueError('Boolean preparation switches required')
    def address(ports):
        return INDEX[tuple(int(port in ports) for port in range(12))]
    left, right = address((0, 1)), address((0, 5))
    vectors = ((left, 0), (right, int(reference)))
    return {(i, j, r, s): F(1, 2)
            for i, r in vectors for j, s in vectors
            if not dephased or (i, r) == (j, s)}


@lru_cache(maxsize=None)
def apparatus_columns(gamma):
    """R_(1,5) after (-1)^(gamma*n0*n1), with exact rational R."""
    if type(gamma) is not int or gamma not in (0, 1):
        raise ValueError('gamma must be zero or one')
    columns = []
    for i, state in enumerate(STATES):
        phase = -1 if gamma and state[0] and state[1] else 1
        pair = state[1], state[5]
        if pair[0] == pair[1]:
            columns.append(((i, F(phase)),))
        else:
            swapped = list(state)
            swapped[1], swapped[5] = swapped[5], swapped[1]
            other = INDEX[tuple(swapped)]
            columns.append(((i, phase * F(3, 5)),
                            (other, phase * F(4 if pair == (1, 0) else -4, 5))))
    return tuple(columns)


def apparatus(rho, gamma):
    out = defaultdict(F)
    columns = apparatus_columns(gamma)
    for (i, j, r, s), value in rho.items():
        for a, u in columns[i]:
            for b, v in columns[j]:
                out[a, b, r, s] += u * v * value
    return _clean(out)


def finish(outputs, word, status, rho):
    """Both occupation-projector outcomes, including zero-density outcomes."""
    for bit in (0, 1):
        outputs[word, status, bit] = {
            key: value for key, value in rho.items()
            if STATES[key[0]][1] == bit and STATES[key[1]][1] == bit}


@lru_cache(maxsize=None)
def _execute(gamma, cap, reference, dephased):
    first, timeout = interval(preparation(reference, dephased), cap)
    output = {}
    finish(output, (), 'timeout_first', timeout)
    for mark, rho in first.items():
        word = (('scheduled', mark),)
        feedback_seam = (0, 1) if mark[1][0] > 0 else (0, 5)
        for feedback_mark, after_feedback in native_step(rho, feedback_seam).items():
            prefix = word + (('feedback', feedback_mark),)
            second, timeout = interval(apparatus(after_feedback, gamma), cap)
            finish(output, prefix, 'timeout_second', timeout)
            for second_mark, final in second.items():
                finish(output, prefix + (('scheduled', second_mark),), 'complete', final)
    return output


def execute(gamma=0, cap=2, reference=False, dephased=False):
    """Independent direct protocol, with copy-safe returned dictionaries."""
    if type(gamma) is not int or gamma not in (0, 1):
        raise ValueError('gamma must be zero or one')
    if type(cap) is not int or not 1 <= cap <= 3:
        raise ValueError('attempt cap from one through three required')
    if type(reference) is not bool or type(dephased) is not bool:
        raise ValueError('Boolean preparation switches required')
    return {key: dict(rho) for key, rho in _execute(gamma, cap, reference, dephased).items()}


def trace(rho):
    return sum((value for (i, j, r, s), value in rho.items()
                if i == j and r == s), F(0))


def reference_trace(rho):
    out = defaultdict(F)
    for (i, j, r, s), value in rho.items():
        if r == s:
            out[i, j, 0, 0] += value
    return _clean(out)


def clock_joint(outputs):
    """Complete (status, tickA, tickB, final occupation) distribution."""
    out = defaultdict(F)
    for (word, status, bit), rho in outputs.items():
        a = sum(mark[0][0] for _, mark in word)
        b = sum(mark[0][1] for _, mark in word)
        out[status, a, b, bit] += trace(rho)
    return dict(sorted(out.items()))


def matched_marginal_control(joint):
    """Break clock/field dependence within status, preserving both marginals.

    Timeouts stay explicit. Conditional independence is imposed only within
    each terminal status, so even all status-conditioned marginal fits agree.
    """
    mass, clocks, bits = defaultdict(F), defaultdict(F), defaultdict(F)
    for (status, a, b, bit), value in joint.items():
        mass[status] += value
        clocks[status, a, b] += value
        bits[status, bit] += value
    return {(status, a, b, bit): clocks[status, a, b] * bits[status, bit] / mass[status]
            if mass[status] else F(0)
            for status, a, b, bit in joint}


def total_variation(first, second):
    return sum((abs(first.get(key, F(0)) - second.get(key, F(0)))
                for key in first.keys() | second.keys()), F(0)) / 2


def clock_detector_covariance(joint, status='complete', ticks=(1, 3)):
    """Conditional indicator covariance from one complete recorded joint law."""
    mass = event = detector = both = F(0)
    for (label, a, b, bit), value in joint.items():
        if label == status:
            mass += value
            if (a, b) == ticks:
                event += value
            if bit == 1:
                detector += value
            if (a, b) == ticks and bit == 1:
                both += value
    if not mass:
        raise ValueError('the conditioning status has zero probability')
    return (both*mass - event*detector)/mass**2


def retained_entanglement():
    """Independent remote-wait replay of a source/reference entangled timeout."""
    def address(ports):
        return INDEX[tuple(int(port in ports) for port in range(12))]
    a, b = address((0, 2)), address((0, 3))
    initial = {(a, a, 0, 0): F(1, 2), (b, b, 1, 1): F(1, 2),
               (a, b, 0, 1): F(1, 2), (b, a, 1, 0): F(1, 2)}
    _, timeout = interval(initial, 2)
    state = {key: value for key, value in timeout.items()
             if STATES[key[0]][1] == 0 == STATES[key[1]][1]}
    witness = (state.get((a, a, 1, 1), F(0)) + state.get((b, b, 0, 0), F(0))
               - state.get((a, b, 0, 1), F(0)) - state.get((b, a, 1, 0), F(0)))/2
    return state, witness


def _serial_joint(joint):
    return [{'status': status, 'ticks': [a, b], 'bit': bit, 'probability': str(value)}
            for (status, a, b, bit), value in sorted(joint.items())]


def audit():
    """Exact response and matched-marginal obstruction; no fitted time scale."""
    # Import only at this comparison boundary. The independent execution above
    # uses no producer matrices, resolvent, Kraus groups or instrument compiler.
    from . import construction as producer

    outputs = [execute(gamma) for gamma in (0, 1)]
    joint = [clock_joint(value) for value in outputs]
    if any(sum(value.values(), F(0)) != 1 for value in joint):
        raise ArithmeticError('complete output law must have mass one')
    controls = [matched_marginal_control(value) for value in joint]
    clock_marginals = []
    bit_marginals = []
    for value in joint:
        clocks, bits = defaultdict(F), defaultdict(F)
        for (status, a, b, bit), probability in value.items():
            clocks[status, a, b] += probability
            bits[status, bit] += probability
        clock_marginals.append(dict(clocks))
        bit_marginals.append(dict(bits))
    if clock_marginals[0] != clock_marginals[1]:
        raise ArithmeticError('the declared coherent setting unexpectedly changed the clock marginal')
    incoherent = [execute(gamma, dephased=True) for gamma in (0, 1)]
    if incoherent[0] != incoherent[1]:
        raise ArithmeticError('a diagonal preparation cannot see the supplied diagonal phase')
    reference = [execute(gamma, reference=True) for gamma in (0, 1)]
    for gamma in (0, 1):
        if {key: reference_trace(rho) for key, rho in reference[gamma].items()} != incoherent[gamma]:
            raise ArithmeticError('discarding the untouched reference must give the diagonal mixture')
    for attached, histories in ((False, outputs), (True, reference)):
        for gamma, expected in enumerate(histories):
            for effective in (False, True):
                if producer.execute(gamma, reference=attached, effective=effective) != expected:
                    raise ArithmeticError('independent endpoint replay differs from complete output')
    nonzero_reference_response = sum(
        abs(reference[0][key].get(entry, F(0)) - reference[1][key].get(entry, F(0)))
        for key in reference[0]
        for entry in reference[0][key].keys() | reference[1][key].keys())
    entangled_timeout, expectation = retained_entanglement()
    expected_timeout = producer.execute(initial=producer.entangled_preparation())[
        (), 'timeout_first', 0]
    if entangled_timeout != expected_timeout or expectation >= 0:
        raise ArithmeticError('independent retained-timeout entanglement witness failed')
    return {
        'classification': 'exact native quantum clock-response correspondence control',
        'parameter': 'recorded incident-attempt counts; no seconds calibration',
        'complete_output_words': len(outputs[0]),
        'independent_native_endpoint_oracle': {
            'complete_source_or_effective_comparisons': 8,
            'all_unnormalized_source_reference_entries_equal': True,
            'zero_probability_words_and_both_detector_outcomes_retained': True},
        'joint_clock_detector_laws': [_serial_joint(value) for value in joint],
        'same_clock_marginal_both_settings': True,
        'full_joint_response_tv': str(total_variation(joint[0], joint[1])),
        'matched_status_clock_and_detector_marginals_control_tv': [
            str(total_variation(value, control)) for value, control in zip(joint, controls)],
        'complete_status_clock_1_3_and_detector_one_covariance': [
            str(clock_detector_covariance(value)) for value in joint],
        'clock_relabeling_scope': 'invariant under bijective relabeling of recorded clock pairs',
        'invariance_under_clock_coarsening_claimed': False,
        'dephased_input_zero_phase_response': True,
        'reference_trace_equals_dephased_process': True,
        'reference_full_matrix_entry_l1_response_not_trace_norm': str(nonzero_reference_response),
        'independent_timeout_reference_entanglement': {
            'branch_mass': str(trace(entangled_timeout)),
            'unnormalized_partial_transpose_expectation': str(expectation),
            'conditional_partial_transpose_expectation': str(expectation/trace(entangled_timeout))},
        'physical_clock_or_field_identification_established': False,
    }

"""Independent endpoint, reference, prefix and phase tests of the joined law.

The oracle below enumerates loads and performs floor/ceiling averaging. It
does not use the joined channel's target maps or the native update function.
These are floating falsification checks of the exact proof in QUANTUM_NATIVE,
not outward numerical or physical apparatus certificates.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
import math

import numpy as np
import pytest

from ..observer.source_discrimination import circuit
from . import composition, protocol
from . import quantum_native as joint


def _states():
    configurations = []
    for first, second in combinations_with_replacement(range(12), 2):
        value = [0]*12
        value[first] += 1
        value[second] += 1
        configurations.append(tuple(value))
    return tuple(sorted(configurations))


STATES = _states()
CODES = tuple(circuit.encode(state) for state in STATES)


def _event(state, seam, coin, law):
    left, right = seam
    low, odd = divmod(state[left]+state[right], 2)
    output = list(state)
    output[left], output[right] = low+coin*odd, low+(1-coin)*odd
    current = state[left]-output[left]
    current_code = {-1: 2, 0: 0, 1: 1}[current]
    raw = current_code+(4*odd if law == 'E' else 0)
    return tuple(output), current, raw


@lru_cache(maxsize=None)
def _native_map(law, first_coin, first_raw, selected, last_coin, last_raw):
    result = np.zeros((512, 78), complex)
    if selected is None:
        return result
    seam = (0, 5) if selected == 0 else (0, 1)
    for column, state in enumerate(STATES):
        middle, _, raw1 = _event(state, (0, 1), first_coin, law)
        final, _, raw2 = _event(middle, seam, last_coin, law)
        if (raw1, raw2) == (first_raw, last_raw):
            result[circuit.encode(final), column] = .5
    return result


def _oracle_map(law, record):
    c1, raw1, bits, c2, raw2 = record
    return _native_map(law, c1, raw1, bits[2] if bits else None, c2, raw2)


def _current(law, coin, raw):
    # The E high bit records endpoint-total parity. The low two bits have
    # the same current meaning as C, independently of the contrast decoder.
    return {0: 0, 1: 1, 2: -1}[raw % 4]


@pytest.fixture(scope='module')
def histories():
    return {law: joint.execute(law, 2, 16) for law in ('C', 'E')}


@pytest.fixture(scope='module')
def work_targets():
    return {current: protocol.execute(protocol.continuation_program(
                2, 16, native_current=current), backend='product')
            for current in (-1, 0, 1)}


@pytest.mark.parametrize('law', ['C', 'E'])
def test_every_joint_branch_has_the_independent_native_and_work_operator(
        histories, work_targets, law):
    history = histories[law]
    assert STATES == joint.STATES
    assert len(STATES) == 78 and len(set(CODES)) == 78
    seen_native = set()
    all_bits = set()
    for branch in history.branches:
        c1, raw1, bits, _, _ = branch.record
        # Native matrices are shared across quantum prefixes with the same
        # second smooth bit. Check all columns of each executed matrix once.
        native_key = (id(branch.native_map), c1, raw1, bits[2] if bits else None,
                      *branch.record[3:])
        if native_key not in seen_native:
            np.testing.assert_allclose(branch.native_map,
                                       _oracle_map(law, branch.record),
                                       atol=3e-11, rtol=0)
            seen_native.add(native_key)
        if not bits:
            assert branch.abort
            np.testing.assert_array_equal(branch.work_map, np.eye(16))
            continue
        all_bits.add(bits)
        expected = next(b for b in work_targets[_current(law, c1, raw1)].branches
                        if b.record == bits)
        np.testing.assert_allclose(branch.work_map, expected.kraus,
                                   atol=3e-11, rtol=0)
    assert len(all_bits) == 16


def test_joint_entanglement_reference_and_first_prefix_survive_complete_continuation(
        histories, work_targets):
    """One input entangles all 78 native basis states, work, and a reference."""
    rng = np.random.default_rng(10990073)
    initial = rng.normal(size=(78, 16, 2))+1j*rng.normal(size=(78, 16, 2))
    initial /= np.linalg.norm(initial)
    expected_reference = np.einsum('ijr,ijs->rs', initial.conj(), initial)
    actual_reference = np.zeros((2, 2), complex)
    prefix_actual, prefix_expected = {}, {}
    error_squared, masses = [], []
    maps = {current: {b.record: b.kraus for b in h.branches}
            for current, h in work_targets.items()}
    for record, actual in joint.output_on_joint_input(histories['C'], initial):
        c1, raw1, bits, _, _ = record
        prefix_actual[c1, raw1] = prefix_actual.get((c1, raw1), 0.)+float(np.vdot(actual, actual).real)
        actual_reference += np.einsum('ijr,ijs->rs', actual.conj(), actual)
        if bits:
            native = _oracle_map('C', record)
            work = maps[_current('C', c1, raw1)][bits]
            expected = np.einsum('ai,bj,ijr->abr', native, work, initial,
                                 optimize=True)
        else:
            expected = np.zeros_like(actual)
        error_squared.append(float(np.linalg.norm(actual-expected)**2))
        masses.append(float(np.vdot(actual, actual).real))
    # Before any work read, the complete first-native record effect is the
    # diagonal current-group projector times its coin probability 1/2.
    for coin in (0, 1):
        for column, state in enumerate(STATES):
            _, _, raw = _event(state, (0, 1), coin, 'C')
            prefix_expected[coin, raw] = prefix_expected.get((coin, raw), 0.)+float(
                np.vdot(initial[column], initial[column]).real)/2
    for key, probability in prefix_actual.items():
        assert abs(probability-prefix_expected.get(key, 0.)) < 3e-11
    assert abs(math.fsum(masses)-1) < 3e-11
    assert math.sqrt(math.fsum(error_squared)) < 3e-10
    np.testing.assert_allclose(actual_reference, expected_reference, atol=3e-11, rtol=0)


def test_uniform_binary_diagonal_native_input_recovers_the_retained_classical_benchmark(histories):
    binary = [index for index, state in enumerate(STATES) if max(state) == 1]
    assert len(binary) == 66
    loads = np.zeros(512)
    for state, code in zip(STATES, CODES, strict=True):
        loads[code] = state[0]
    contributions = {key: [] for key in (
        'probability_mass', 'mean_final_port0_load',
        'mean_final_port0_load_disconnected_actuator',
        'last_smooth_plus_probability', 'clock_bit_one_probability')}
    for branch in histories['C'].branches:
        c1, _, bits, c2, _ = branch.record
        columns = abs(branch.native_map[:, binary])**2
        work_weight = float(np.vdot(branch.work_state, branch.work_state).real)/66
        probability = work_weight*float(columns.sum())
        contributions['probability_mass'].append(probability)
        contributions['mean_final_port0_load'].append(
            work_weight*float(np.sum(loads[:, None]*columns)))
        if not bits:
            continue  # Its mass was retained; the halted path has no later bit.
        contributions['last_smooth_plus_probability'].append(probability*(bits[2] == 0))
        contributions['clock_bit_one_probability'].append(probability*bits[3])
        # The old benchmark also publishes a disconnected-actuator control.
        # Its observable depends on the original diagonal input and the two
        # coins, so evaluate that counterfactual without normalizing columns.
        disconnected = []
        for column in binary:
            middle, _, _ = _event(STATES[column], (0, 1), c1, 'C')
            final, _, _ = _event(middle, (0, 1), c2, 'C')
            disconnected.append(final[0])
        contributions['mean_final_port0_load_disconnected_actuator'].append(
            work_weight*float(np.dot(np.sum(columns, axis=0), disconnected)))
    programs = {current: protocol.continuation_program(2, 16, native_current=current)
                for current in (-1, 0, 1)}
    old = composition.native_experiment(programs)['summary']
    for key, values in contributions.items():
        assert abs(math.fsum(values)-old[key]) < 3e-11, key


def _phase_input(theta):
    state = np.zeros(78, complex)
    state[CODES.index(0)] = 1/math.sqrt(2)
    state[CODES.index(1)] = np.exp(1j*theta)/math.sqrt(2)
    return state


def _plus_probability(history, theta):
    # Direct rank-one effect on the native output, without importing the
    # implementation's interference formula or fine-readout calculation.
    initial = _phase_input(theta)
    probabilities = []
    for branch in history.branches:
        native = branch.native_map@initial
        if not branch.abort:
            probabilities.append(abs((native[0]+native[1])/math.sqrt(2))**2
                                 *float(np.vdot(branch.work_state, branch.work_state).real))
    return math.fsum(probabilities)


def test_surviving_native_phase_obeys_the_second_smooth_feedback_formula(
        histories, work_targets):
    q = math.fsum(b.probability for b in work_targets[0].branches if b.record[2] == 0)
    values = {law: [_plus_probability(history, phase)
                    for phase in (0, math.pi, math.pi/2)]
              for law, history in histories.items()}
    assert abs(values['C'][0]-values['C'][1]-(1+q)/4) < 3e-11
    assert abs(values['E'][0]-values['E'][1]) < 3e-11
    for law in ('C', 'E'):
        assert abs(values[law][2]-(values[law][0]+values[law][1])/2) < 3e-11
    assert abs(values['C'][2]-values['E'][2]) < 3e-11
    assert abs(values['C'][0]-values['E'][0]-(1+q)/8) < 3e-11
    # The actual emitted Y rotation plus full 512-code read realizes the
    # same single effect while retaining its strictly finer outcome menu.
    readout = joint.readout(histories['C'], _phase_input(0))
    assert abs(readout['plus_code_zero_probability']-values['C'][0]) < 3e-11
    assert abs(readout['probability_mass']-1) < 3e-11
    valid_count = sum(not b.abort for b in histories['C'].branches)
    abort_count = sum(b.abort for b in histories['C'].branches)
    assert readout['complete_outcome_count'] == 512*valid_count+abort_count

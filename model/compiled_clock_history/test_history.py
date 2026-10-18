"""Independent source-chart, complete-output and continuation controls."""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
from sympy.polys.domains import QQ_I

from ..observer import native, quantum_mixing
from ..quantum_observer_correspondence import clock_response as independent
from . import history as h, source


UNION = (0, 1, 5)
OUTSIDE = (2, 3, 4, 6, 7, 8, 9, 10, 11)
STATES = native.states(binary=True)


def chart(state):
    bits = tuple(state[p] for p in OUTSIDE)
    shell = sorted(row for row in product((0, 1), repeat=9) if sum(row) == sum(bits))
    return sum(state[p] << bit for bit, p in enumerate(UNION))+8*shell.index(bits)


CODES = tuple(map(chart, STATES))


def scalar(value):
    return QQ_I(F(float(value.real)), F(float(value.imag)))


def lift(outputs):
    return {word: {(CODES[i], CODES[j], r, s): QQ_I(value)
                   for (i, j, r, s), value in state.items()}
            for word, state in outputs.items()}


def entry_error(first, second):
    total = F(0)
    for record in first.keys() | second.keys():
        left, right = first.get(record, {}), second.get(record, {})
        for key in left.keys() | right.keys():
            difference = left.get(key, QQ_I.zero)-right.get(key, QQ_I.zero)
            for component in (difference.x, difference.y):
                total += F(abs(int(component.numerator)), int(component.denominator))
    return total/2


def test_chart_is_fixed_on_all_native_states_and_all_padding_is_refused():
    assert source.codes() == CODES and len(set(CODES)) == 66
    for state, code in zip(STATES, CODES, strict=True):
        assert source.encode(state) == code
        assert source.decode(code) == state
    for code in set(range(512))-set(CODES):
        with pytest.raises(ValueError):
            source.decode(code)
    # A local ancilla at temporary bit3 must not replace a spectator-shell bit.
    assert source.resources('readout')['blank_pointer_qubits'] == 1
    assert source.resources('preparation', reference=True)['blank_reference_qubits'] == 1


@pytest.mark.parametrize('reference', (False, True))
def test_preparation_uses_original_binary_components_and_explicit_reference_columns(reference):
    vector = source.prepared(reference)
    actual = h.prepare(reference)
    expected = {}
    for i, x in enumerate(vector):
        for j, y in enumerate(vector):
            a, b = scalar(x), scalar(y)
            value = a*QQ_I(b.x, -b.y)
            if value:
                expected[i % 8, j % 8, i//8, j//8] = value
    assert actual == expected
    assert abs(h.trace(actual)-1) < F(1, 10**12)
    # Structural equality above checks every returned component even on a
    # platform where one compiler cancellation happens to be exactly zero.


def test_every_native_partial_permutation_uses_original_endpoint_and_current():
    expected = {}
    for seam in quantum_mixing.seams():
        groups = {}
        for coin in (0, 1):
            for state in STATES:
                after, record = native.native(state, seam, coin)
                assert record[0] == coin
                groups.setdefault(record, {})[chart(state)] = chart(after)
        expected[seam] = tuple((record[1], tuple(sorted(mapping.items())))
                               for record, mapping in sorted(groups.items()))
    assert dict(h._native_maps()) == expected


def test_native_extension_keeps_complement_failure_and_destroys_only_cross_code_terms():
    good, bad = CODES[0], 7  # 7 has three union occupations and is outside the code.
    state = {(i, j, 0, 0): QQ_I(F(1, 2)) for i, j in product((good, bad), repeat=2)}
    completed, waiting, aborts = h.interval(state, 2)
    assert aborts[0] == {(bad, bad, 0, 0): QQ_I(F(1, 2))}
    assert aborts[1] == {}
    assert sum(map(h.trace, (*completed.values(), waiting, *aborts))) == 1
    for output in (*completed.values(), waiting):
        assert all(i in CODES and j in CODES for i, j, _, _ in output)
    assert h.trace(waiting) == F(9, 50)


def test_local_map_lifts_spectator_identity_without_rounding_or_normalization():
    matrix = np.eye(8, dtype=complex)
    matrix[0, 0] = 1.+2.**-20
    matrix[1, 0] = 1j/8
    # Spectator coherences and a nontrivial reference matrix unit both survive.
    state = {(8, 16, 1, 0): QQ_I(1, 2)}
    expected = {}
    for a in (0, 1):
        for b in (0, 1):
            x, y = scalar(matrix[a, 0]), scalar(matrix[b, 0])
            expected[8+a, 16+b, 1, 0] = x*QQ_I(1, 2)*QQ_I(y.x, -y.y)
    assert h.local(matrix, state) == expected
    positive = h.local(matrix, {(8, 8, 0, 0): QQ_I.one})
    assert h.trace(positive) == F(1+2.**-20)**2+F(1, 64) > 1


@pytest.fixture(scope='module')
def executed():
    return {(gamma, reference): h.execute(gamma=gamma, cap=1, reference=reference)
            for gamma, reference in product((0, 1), (False, True))}


@pytest.mark.parametrize('gamma,reference', product((0, 1), (False, True)))
def test_all_raw_outputs_agree_with_independent_native_history(executed, gamma, reference):
    actual = executed[gamma, reference]
    target = lift(independent.execute(gamma=gamma, cap=1, reference=reference))
    normal = {record for record in actual if not record[1].startswith('code_abort')}
    assert normal == set(target) and len(normal) == 1442
    assert len(actual) == 1564
    error = entry_error(actual, target)
    assert error == h.distance_upper(actual, target) < F(1, 10**11)
    assert abs(sum(map(h.trace, actual.values()))-1) < F(1, 10**12)
    assert sum(h.trace(state) for record, state in actual.items()
               if record[1].startswith('code_abort')) > 0
    first = sum(h.trace(state) for record, state in actual.items() if record[1] == 'timeout_first')
    second = sum(h.trace(state) for record, state in actual.items() if record[1] == 'timeout_second')
    assert abs(first-F(3, 5)) < F(1, 10**12)
    assert abs(second-F(6, 25)) < F(1, 10**12)


def test_interaction_is_seen_in_joint_outputs_even_when_clock_marginals_agree(executed):
    def joint(output):
        result = defaultdict(F)
        for (word, status, bit), state in output.items():
            ticks = tuple(sum(mark[0][i] for _, mark in word) for i in range(2))
            result[status, ticks, bit] += h.trace(state)
        return result
    a, b = (joint(executed[gamma, False]) for gamma in (0, 1))
    tv = sum(abs(a.get(key, 0)-b.get(key, 0)) for key in a.keys() | b.keys())/2
    assert tv > F(1, 100)
    for status, ticks, _ in a:
        assert abs(sum(a.get((status, ticks, bit), 0)-b.get((status, ticks, bit), 0)
                       for bit in (0, 1))) < F(1, 10**12)
    assert entry_error(executed[0, True], executed[1, True]) > F(1, 100)


def test_zero_erasing_padding_does_not_equal_the_executed_complete_instrument(executed):
    actual = executed[1, False]
    erased = {word: state for word, state in actual.items() if not word[1].startswith('code_abort')}
    assert h.distance_upper(actual, erased) > 0
    missing_timeout = {word: state for word, state in actual.items() if word[1] == 'complete'}
    assert h.distance_upper(actual, missing_timeout) > F(2, 5)


def test_code_failure_still_executes_both_terminal_pointer_outcomes(monkeypatch):
    ket = np.zeros((512, 1), complex)
    ket[7, 0] = 1
    monkeypatch.setattr(source, 'preparation', lambda reference=False: ket)
    actual = h.execute(cap=1)
    expected = []
    for matrix in source.local_readout():
        column = tuple(map(scalar, matrix[:, 7]))
        expected.append({(i, j, 0, 0): x*QQ_I(y.x, -y.y)
                         for i, x in enumerate(column) for j, y in enumerate(column) if x and y})
    for bit in (0, 1):
        assert actual[(), 'code_abort_first_0', bit] == expected[bit]
    assert all(not state for record, state in actual.items() if record[1] != 'code_abort_first_0')
    assert abs(sum(map(h.trace, actual.values()))-1) < F(1, 10**12)


@pytest.mark.parametrize('component', ('x', 'y'))
@pytest.mark.parametrize('numerator,denominator', ((1.5, 1), (1, 2.5), (True, 1),
    (1, True), (-1, -1), (2, 2), (1, 0)))
def test_mutated_gaussian_rational_slots_cannot_create_false_zero_bounds(component, numerator, denominator):
    value = QQ_I(F(123, 127), F(17, 19))  # Fresh components; never mutate domain constants.
    part = getattr(value, component)
    part.numerator, part.denominator = numerator, denominator
    state = {(0, 0, 0, 0): value}
    for operation in (lambda: h.conjugate(value), lambda: h.trace(state),
                      lambda: h.distance_upper({'record': state}, {})):
        with pytest.raises(ValueError):
            operation()


def test_malformed_exact_scalar_is_rejected_by_complete_output_acceptance():
    from . import receipt
    actual = {record: {} for record in h.record_menu(1)}
    word = ((), 'timeout_first', 0)
    damaged = QQ_I(F(123, 127))
    damaged.x.numerator, damaged.x.denominator = 1.5, 1
    actual[word] = {(0, 0, 0, 0): damaged}
    target = {word: {(0, 0, 0, 0): QQ_I.one}}
    with pytest.raises(ValueError):
        receipt.audit_output(actual, target, cap=1)


def test_exact_entrywise_bound_keeps_honest_small_real_and_imaginary_parts():
    state = {(0, 0, 0, 0): QQ_I(F(-1, 2), F(1, 3))}
    assert h.distance_upper({'record': state}, {}) == F(5, 12)


def test_serialized_continuation_reuses_the_exact_retained_prefix():
    stage = h.first_stage(gamma=1, cap=1, reference=True)
    payload = h.checkpoint_json(stage)
    resumed = h.execute(gamma=1, cap=1, reference=True, checkpoint=payload)
    whole = h.execute(gamma=1, cap=1, reference=True)
    assert resumed == whole and len(whole) == 1564


def test_highest_declared_cap_composes_executor_and_independent_target():
    actual = h.execute(cap=h.MAX_CAP)
    target = h.ideal_history(cap=h.MAX_CAP)
    independently_executed = lift(independent.execute(gamma=1, cap=h.MAX_CAP))
    assert target == independently_executed
    assert h.distance_upper(actual, target) < F(1, 10**11)
    timeout = sum(h.trace(state) for (_, status, _), state in target.items()
                  if status.startswith('timeout'))
    assert timeout == 1-(1-F(3, 5)**h.MAX_CAP)**2


@pytest.mark.parametrize('mutation', ('state', 'record', 'resource', 'binding', 'reference', 'cap'))
def test_checkpoint_corruption_is_rejected(mutation):
    stage = h.first_stage(gamma=1, cap=1, reference=True)
    value = json.loads(h.checkpoint_json(stage))
    if mutation == 'state':
        value['waiting'][0][1] = str(F(value['waiting'][0][1])+F(1, 7))
    elif mutation == 'record':
        value['completed'][0][0][0][0] = True
    elif mutation == 'resource':
        name = next(iter(value['preparation_resources']))
        value['preparation_resources'][name] += 1
    elif mutation == 'binding':
        name = next(iter(value['source_binding']))
        value['source_binding'][name] = '0'*64
    elif mutation == 'reference':
        value['parameters'][2] = 1
    else:
        value['scheduled_attempt_cap'] = True
    payload = json.dumps(value, sort_keys=True, separators=(',', ':'))
    with pytest.raises(ValueError):
        h.execute(gamma=1, cap=1, reference=True, checkpoint=payload)


@pytest.mark.parametrize('kwargs', ({'gamma': True}, {'gamma': 2}, {'cap': True},
    {'cap': 0}, {'cap': 4}, {'cap': 5}, {'reference': 1}, {'checkpoint': {}}))
def test_original_protocol_types_are_checked(kwargs):
    with pytest.raises(ValueError):
        h.execute(**kwargs)


@pytest.mark.parametrize('bad', (np.full((8, 8), np.nan), np.eye(8, dtype=bool),
    np.eye(8, dtype=np.int64), np.ma.array(np.eye(8), mask=False)))
def test_local_array_boundary_rejects_nonbinary_or_masked_data(bad):
    with pytest.raises(ValueError):
        h.local(bad, {(0, 0, 0, 0): QQ_I.one})


@pytest.mark.skipif(np.dtype(np.longdouble).itemsize <= 8, reason='no extended binary dtype on this platform')
def test_local_array_never_rounds_extended_precision_before_validation():
    with pytest.raises(ValueError):
        h.local(np.eye(8, dtype=np.longdouble), {(0, 0, 0, 0): QQ_I.one})

"""Native clock source execution, exact stopped laws and adverse controls."""
from fractions import Fraction as F
from itertools import product

import mpmath
import numpy as np
import pytest

from ..native_relational_readout import construction as relational
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime
from . import clock


def _sum(*states):
    result = {}
    for state in states:
        for key, value in state.items():
            result[key] = result.get(key, F(0))+value
    return {key: value for key, value in result.items() if value}


@pytest.mark.parametrize('law,count,activity', tuple(product(('C', 'E'), (1, 4, 16, 64, 128), (F(1), F(1, 2)))))
def test_complete_exact_signed_distribution_and_resource_moments(law, count, activity):
    result = clock.statistics(law, count, activity=activity)
    assert tuple(result['probabilities']) == clock.records(count)
    assert len(result['probabilities']) == 2*count+1
    assert sum(result['probabilities'].values()) == 1
    assert 0 < result['timeout_probability'] <= result['uniform_tail_bound'] < 1
    assert result['expected_native_attempts'] == activity*result['expected_controller_slots']
    assert result['expected_native_attempts'] < F(33, 5)
    assert result['expected_controller_slots'] <= count
    for record in result['probabilities']:
        cost = dict(clock.outcome_resources(count, record, activity=activity))
        assert cost['controller_slots'] == record[0]
        if activity == 1:
            assert cost['native_attempts'] == cost['native_current_readouts'] == record[0]
        else:
            assert 'native_attempts' not in cost
            assert cost['native_attempts_upper'] == record[0]


@pytest.mark.parametrize('law,sign', tuple(product(('C', 'E'), (-1, 1))))
def test_signed_outputs_recover_existing_exact_absolute_current_instrument(law, sign):
    count = 4
    marked = relational.source(law)
    actual = clock.exact_outputs(law, count, sign=sign)
    for n in range(1, count+1):
        assert _sum(actual[n, -1], actual[n, 1]) == relational.exact_stopped_state(marked, sign=sign, first_at=n)
    assert actual[count, 0] == relational.exact_tail(marked, count, sign)
    eventual = clock.exact_eventual_outputs(law, sign=sign)
    assert _sum(*eventual.values()) == relational.exact_stopped_state(marked, sign=sign)
    # Exact first-sign probabilities are asymmetric for the declared source
    # seam orientation; replacing current by an invented fair sign fails.
    assert clock.statistics(law, 1, sign=sign)['probabilities'] == {
        (1, -1): F(1, 20), (1, 1): F(31, 300), (1, 0): F(127, 150)}


@pytest.mark.parametrize('law', ('C', 'E'))
def test_actual_numerical_outputs_agree_with_independent_latched_runtime(law):
    count = 3
    marked = relational.source(law)
    # This earlier implementation groups absolute currents and executes every
    # latch word through the shared runtime, independent of our stopped loop.
    program = relational.program(marked, count)
    history = runtime.execute(program, stop_after=count)
    actual = clock.execute_native(law, count)
    for branch in history.branches:
        flags = branch.record
        if not all(a <= b for a, b in zip(flags, flags[1:])):
            np.testing.assert_array_equal(branch.density, np.zeros((66, 66)))
            continue
        if 1 in flags:
            n = flags.index(1)+1
            expected = actual[n, -1]+actual[n, 1]
        else:
            expected = actual[count, 0]
        np.testing.assert_allclose(branch.density, expected, atol=3e-15, rtol=0)
    bound = clock.compare_outputs(law, count, actual)
    assert 0 < bound['distance_upper'] < F(1, 10**12)
    assert 1 < bound['reported_nuclear_norm_upper'] < 1+F(1, 10**11)
    assert all(not value.flags.writeable for value in actual.values())


def test_numerical_execution_does_not_substitute_exact_target(monkeypatch):
    clock.source('C')
    def forbidden(*args, **kwargs):
        raise AssertionError('exact target was substituted for numerical native execution')
    monkeypatch.setattr(clock, 'exact_outputs', forbidden)
    monkeypatch.setattr(relational, 'exact_preparation', forbidden)
    actual = clock.execute_native('C', 2)
    assert tuple(actual) == clock.records(2)
    assert abs(sum(np.trace(value) for value in actual.values())-1) < 3e-15


@pytest.mark.parametrize('law', ('C', 'E'))
def test_native_numeric_idling_matches_exact_targets(law):
    arrays = clock.execute_native(law, 4, activity=F(1, 2))
    answer = clock.compare_outputs(law, 4, arrays, activity=F(1, 2))
    assert 0 < answer['distance_upper'] < F(1, 10**12)
    assert clock.compare_outputs(law, 4, arrays)['distance_upper'] > F(1, 100)


def test_C_E_timing_agrees_while_retained_native_coherence_differs():
    c, e = [clock.exact_outputs(law, 4) for law in ('C', 'E')]
    assert clock.statistics('C', 4)['probabilities'] == clock.statistics('E', 4)['probabilities']
    assert c[1, -1] == e[1, -1] and c[1, 1] == e[1, 1]
    assert c[2, -1] != e[2, -1] or c[2, 1] != e[2, 1]
    assert c[4, 0] != e[4, 0]


@pytest.mark.parametrize('law', ('C', 'E'))
def test_cap_extension_preserves_earlier_stops_and_keeps_full_tail(law):
    small, large = clock.exact_outputs(law, 4), clock.exact_outputs(law, 16)
    for record, state in small.items():
        if record[1]:
            assert large[record] == state
    unresolved_mass = sum((clock._trace(state) for (n, j), state in large.items()
                           if n > 4 or j == 0), F(0))
    assert unresolved_mass == clock._trace(small[4, 0])
    # Extending a cap advances the old timeout state; it does not preserve
    # that state's complete quantum output or replace it by a reset input.
    assert large[16, 0] != small[4, 0]


def test_numeric_prefixes_reuse_stops_and_retain_actual_earlier_timeout(monkeypatch):
    original = runtime._apply
    calls = []
    def counted(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(runtime, '_apply', counted)
    prefixes = clock.execute_prefixes('C', 4, activity=F(1, 2))
    assert len(calls) == 4
    for cap, arrays in prefixes.items():
        assert tuple(arrays) == clock.records(cap)
        assert arrays[1, -1] is prefixes[1][1, -1]
        assert clock.compare_outputs('C', cap, arrays, activity=F(1, 2))['distance_upper'] < F(1, 10**12)
    assert not np.array_equal(prefixes[1][1, 0], prefixes[4][4, 0])
    separately = clock.execute_native('C', 2, activity=F(1, 2))
    for record, array in separately.items():
        np.testing.assert_array_equal(array, prefixes[2][record])


@pytest.mark.parametrize('law', ('C', 'E'))
def test_signed_instrument_retains_arbitrary_reference_and_uniform_tail_bound(law):
    marked, original = clock.source(law), relational.source(law)
    entangled = quantum.reference_preparation()
    for n in range(1, 4):
        signed = runtime._apply(marked.instrument, entangled, 2)
        absolute = runtime._apply(original.instrument, entangled, 2)
        np.testing.assert_allclose(signed[0]+signed[2], absolute[1], atol=2e-16, rtol=0)
        np.testing.assert_allclose(signed[1], absolute[0], atol=2e-16, rtol=0)
        entangled = signed[1]
        assert np.trace(entangled).real <= float(F(13, 15)**n)+2e-15
    assert np.linalg.norm(entangled.reshape(66, 2, 66, 2)[:, 0, :, 1]) > 0


@pytest.mark.parametrize('law', ('C', 'E'))
def test_native_eventual_idle_invariance_does_not_imply_work_idle_invariance(law):
    assert clock.exact_eventual_outputs(law) == clock.exact_eventual_outputs(law, activity=F(1, 2))
    # A work qubit toggled once per controller slot is a simple exact work
    # observable. Native-only idling changes its complete stopped output.
    values = []
    for activity in (F(1), F(1, 2)):
        probabilities = clock.statistics(law, 64, activity=activity)['probabilities']
        values.append(sum(((-1)**n*probability for (n, _), probability in probabilities.items()), F(0)))
    assert abs(values[0]-values[1]) > F(1, 50)
    # Joint idling instead cancels on every native/work matrix unit. For
    # this toggle the no-event eigenvalue on work coherence is -w.
    q = F(1, 2)
    for row in clock.source(law).no_current_counts:
        for count in row:
            w = F(count, 60)
            assert -q/(60*(1-((1-q)-q*w))) == -F(1, 60)/(1+w)


def test_directed_comparison_preserves_nonhermitian_defect_and_large_norm_gain():
    arrays = clock.execute_native('C', 1)
    bad = dict(arrays)
    changed = arrays[1, -1].copy()
    changed[0, 1] += 3j
    bad[1, -1] = changed
    answer = clock.compare_outputs('C', 1, bad)
    row = answer['per_record'][1, -1]
    target = relational.dense(clock.exact_outputs('C', 1)[1, -1])
    measured = np.linalg.svd(changed-target, compute_uv=False).sum()/2
    assert row['distance_upper'] > F(3, 2)
    assert measured <= float(row['distance_upper'])
    assert np.linalg.svd(changed, compute_uv=False).sum() <= float(row['reported_nuclear_norm_upper'])
    # Independently evaluate the exact Frobenius comparison at high precision.
    ctx = mpmath.mp.clone(); ctx.dps = 90
    target = clock.exact_outputs('C', 1)[1, -1]
    squared = ctx.mpf(0)
    for i in range(66):
        for j in range(66):
            value = target.get((i, j), F(0))
            real = F(float(changed[i, j].real))-value
            imag = F(float(changed[i, j].imag))
            squared += (ctx.mpf(real.numerator)/real.denominator)**2+(ctx.mpf(imag.numerator)/imag.denominator)**2
    exact_squared = row['frobenius_error_squared']
    assert abs(squared-ctx.mpf(exact_squared.numerator)/exact_squared.denominator) < ctx.mpf('1e-85')
    upper = row['distance_upper']
    assert ctx.sqrt(66*squared)/2 <= ctx.mpf(upper.numerator)/upper.denominator


def test_sparse_comparison_keeps_smallest_binary64_defect_outside_target_support():
    arrays = clock.execute_native('C', 1)
    original = clock.compare_outputs('C', 1, arrays)['per_record'][1, -1]['frobenius_error_squared']
    target = clock.exact_outputs('C', 1)[1, -1]
    index = next((i, j) for i in range(66) for j in range(66)
                 if (i, j) not in target and arrays[1, -1][i, j] == 0)
    changed = arrays[1, -1].copy()
    tiny = np.nextafter(0., 1.)
    changed[index] = tiny
    arrays[1, -1] = changed
    actual = clock.compare_outputs('C', 1, arrays)['per_record'][1, -1]['frobenius_error_squared']
    assert actual-original == F(float(tiny))**2 > 0


@pytest.mark.parametrize('count', (True, False, 0, -1, 129, 1., F(1), np.int64(1)))
def test_invalid_caps_are_rejected(count):
    for function in (clock.exact_outputs, clock.execute_native, clock.execute_prefixes, clock.statistics):
        with pytest.raises(ValueError):
            function('C', count)
    with pytest.raises(ValueError):
        clock.records(count)


@pytest.mark.parametrize('activity', (True, 1, 1., 0., F(0), F(-1), F(2), np.float64(.5)))
def test_nonfraction_or_invalid_activity_is_rejected(activity):
    with pytest.raises(ValueError):
        clock.exact_outputs('C', 1, activity=activity)
    with pytest.raises(ValueError):
        clock.execute_native('C', 1, activity=activity)
    with pytest.raises(ValueError):
        clock.exact_eventual_outputs('C', activity=activity)


@pytest.mark.parametrize('warm', (False, True))
def test_original_labels_and_sign_are_validated_before_cached_source(warm):
    class Alias(str):
        pass
    clock._source.cache_clear()
    if warm:
        clock.source('C')
    with pytest.raises(ValueError):
        clock.source(Alias('C'))
    for sign in (True, False, 1., np.int64(1), 0, 2):
        with pytest.raises(ValueError):
            clock.exact_outputs('C', 1, sign=sign)


def test_bad_complete_record_menus_and_original_numeric_arrays_are_rejected():
    arrays = clock.execute_native('C', 1)
    for records in ({(1, -1): arrays[1, -1]}, {**arrays, (2, 0): arrays[1, 0]},
                    {(True, j): value for (_, j), value in arrays.items()}):
        with pytest.raises(ValueError):
            clock.compare_outputs('C', 1, records)
    for value in ([[]], np.zeros((2, 2)), np.zeros((66, 66), bool),
                  np.zeros((66, 66), object), np.ma.array(np.zeros((66, 66))),
                  np.full((66, 66), np.inf), np.full((66, 66), np.nan)):
        with pytest.raises(ValueError):
            clock.compare_outputs('C', 1, {**arrays, (1, -1): value})
    for record in ((True, 1), (1, 1.), (1, 0), (5, 1), (0, -1), [1, -1]):
        with pytest.raises(ValueError):
            clock.outcome_resources(4, record)


def test_fresh_exact_outputs_cannot_poison_later_targets():
    states = clock.exact_outputs('C', 1)
    states[1, -1].clear()
    assert clock.exact_outputs('C', 1)[1, -1]

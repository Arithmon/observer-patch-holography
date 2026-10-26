"""Independent sharp-capacity witness and complete source execution checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import math

import numpy as np
import pytest

from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..repair_scheduler import coins
from . import native_routing as routing, record_code as code


def _oracle(law, seam, coin, state):
    a, b = seam
    low, odd = divmod(state[a]+state[b], 2)
    final = list(state)
    final[a], final[b] = low+coin*odd, low+(1-coin)*odd
    raw = {0: 0, 1: 1, -1: 2}[state[a]-final[a]]
    return circuit.encode(tuple(final)), raw+(4*odd if law == 'E' else 0)


def _integer_permutation(edges):
    images = list(range(512))
    for target, pattern in edges:
        mask = 511 ^ (1 << target)
        images = [v ^ (1 << target) if v & mask == pattern else v for v in images]
    return tuple(images)


def test_six_matchings_partition_nonedge_graph_and_attain_sharp_arbitrary_code_bound():
    certificate = code.exact_certificate()
    assert certificate['maximum_exact_record_correctable_code_dimension'] == 6
    assert certificate['upper_bound_native_carrier_dimension'] == 78
    assert certificate['full_carrier_binary_support_derivation']['binary_support_derived_not_assumed']
    assert certificate['arbitrary_complex_codes_and_references_in_upper_bound']
    pairs = tuple(pair for row in code.MATCHINGS for pair in row)
    assert len(set(pairs)) == 36
    assert set(pairs).isdisjoint(code.SEAMS)
    assert set(pairs) | set(code.SEAMS) == set(combinations(range(12), 2))
    assert Counter(p for pair in pairs for p in pair) == Counter({p: 6 for p in range(12)})
    for row in code.MATCHINGS:
        assert Counter(p for pair in row for p in pair) == Counter(range(12))


@pytest.mark.parametrize('law', ('C', 'E'))
def test_full78_singleton_effects_follow_independent_floor_arithmetic_on_every_raw_record(law):
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    assert len(states) == 78
    doubles = {port: index[tuple(2*int(i == port) for i in range(12))] for port in range(12)}
    witnessed = Counter()
    ranks = Counter()
    for seam in code.SEAMS:
        a, b = seam
        expected = {(coin, raw): [] for coin in (0, 1)
                    for raw in range(1 << circuit.RECORD_BITS[law])}
        # Independent arithmetic: floor quotients, no implementation's
        # divmod, native endpoint call, support certificate or branch images.
        for coin in (0, 1):
            for i, state in enumerate(states):
                total = state[a]+state[b]
                left, right = (total+coin)//2, (total+1-coin)//2
                final = list(state)
                final[a], final[b] = left, right
                current = state[a]-left
                raw = 0 if current == 0 else 1 if current == 1 else 2
                if law == 'E':
                    raw += 4*(total-2*(total//2))
                expected[coin, raw].append((i, index[tuple(final)]))
        actual = {key: [] for key in expected}
        for coin, label, mapping in coins.coin_maps(law, seam, binary=False):
            raw = {value: raw for raw, value in circuit.record_values(law, coin).items()}[label]
            actual[coin, raw] = list(mapping)
            ranks[len(mapping)] += 1
        assert actual == expected
        singleton_records = ((0, 2, b), (1, 1, a)) if law == 'C' else tuple(
            (coin, raw, port) for coin in (0, 1) for raw, port in ((1, a), (2, b)))
        for coin, raw, port in singleton_records:
            assert [i for i, _ in actual[coin, raw]] == [doubles[port]]
            witnessed[port] += 1
        if law == 'C':
            # Dropping or reversing the recorded coin destroys the claimed
            # rank-one effect: the opposite coin has eleven inputs.
            assert len(actual[0, 1]) == len(actual[1, 2]) == 11
        else:
            assert len(actual[0, 0]) == len(actual[1, 0]) == 56
    expected_count = 5 if law == 'C' else 10
    assert witnessed == Counter({port: expected_count for port in range(12)})
    assert ranks == (Counter({1: 60, 11: 60, 66: 60}) if law == 'C'
                     else Counter({1: 120, 10: 120, 56: 60}))


def test_double_occupation_is_detected_in_a_two_dimensional_code_not_assumed_away():
    states = native.states()
    double = states.index((2,)+(0,)*11)
    other = next(i for i, state in enumerate(states) if state[0] == 0 and max(state) == 1)
    columns = np.zeros((78, 2), dtype=int)
    columns[double, 0], columns[other, 1] = 1, 1
    # On edge(0,1), C coin1/current+1 (raw1) isolates exactly double0.
    mapping = next(mapping for coin, label, mapping in coins.coin_maps('C', (0, 1), binary=False)
                   if coin == 1 and label == 1)
    effect = np.zeros(78, dtype=int)
    effect[[i for i, _ in mapping]] = 1
    gram = columns.T@(effect[:, None]*columns)
    np.testing.assert_array_equal(gram, np.diag([1, 0]))
    assert not np.array_equal(2*gram, np.trace(gram)*np.eye(2, dtype=int))


def test_full78_certificate_rejects_silently_omitted_double_occupation_effect(monkeypatch):
    original = coins.coin_maps

    def missing_singleton(law, seam, binary=True):
        rows = original(law, seam, binary=binary)
        if not binary and law == 'C' and seam == (0, 1):
            return tuple((coin, label, () if coin == 1 and label == 1 else mapping)
                         for coin, label, mapping in rows)
        return rows

    monkeypatch.setattr(coins, 'coin_maps', missing_singleton)
    with pytest.raises(ArithmeticError, match='full78 source effects'):
        code.exact_certificate()


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_raw_branches_have_exact_scalar_gram_with_independent_endpoint_oracle(law):
    for seam in code.SEAMS:
        total = F(0)
        for coin in (0, 1):
            for raw in range(1 << circuit.RECORD_BITS[law]):
                images = tuple(tuple(output for state in row
                                     for output, actual in [_oracle(law, seam, coin, state)]
                                     if actual == raw) for row in code.SUPPORTS)
                assert code.branch_images(law, seam, coin, raw) == images
                sizes = tuple(map(len, images))
                assert len(set(sizes)) == 1
                assert len(set(x for row in images for x in row)) == sum(sizes)
                probability = F(sizes[0], 12)
                assert code.branch_weight(law, seam, coin, raw) == probability
                total += probability
        assert total == 1


@pytest.mark.parametrize('law', ('C', 'E'))
def test_exact_source_effects_supply_upper_bound_operator_identities(law):
    states = native.states(binary=True)
    occupations = np.asarray(states, dtype=int).T
    np.testing.assert_array_equal(occupations.sum(axis=0), np.full(66, 2))
    nonedge = np.asarray([tuple(i for i, n in enumerate(x) if n) not in code.SEAMS for x in states])
    assert np.count_nonzero(nonedge) == 36
    assert set(occupations[:, nonedge].sum(axis=1)) == {6}
    for seam in code.SEAMS:
        a, b = seam
        effects = []
        for _, _, mapping in coins.coin_maps(law, seam):
            effect = np.zeros(66, dtype=int)
            effect[[i for i, _ in mapping]] = 1
            effects.append(effect)
        ab = occupations[a]*(1-occupations[b])
        ba = occupations[b]*(1-occupations[a])
        assert sum(ab) == sum(ba) == 10
        assert any(np.array_equal(ab, d) for d in effects)
        assert any(np.array_equal(ba, d) for d in effects)
        np.testing.assert_array_equal(ab-ba, occupations[a]-occupations[b])
        assert np.count_nonzero(occupations[a]-ab) == 1


@pytest.mark.parametrize('law', ('C', 'E'))
def test_declared_same66_recovery_is_complete_unitary_for_every_native_branch(law):
    states = native.states(binary=True)
    source_rows = [circuit.encode(state) for state in states]
    v = code.ideal_code_matrix()[source_rows]
    source_index = {value: i for i, value in enumerate(source_rows)}
    for seam in code.SEAMS:
        for coin in (0, 1):
            for raw in circuit.record_values(law, coin):
                recovery = code.declared_recovery_matrix66(law, seam, coin, raw)
                np.testing.assert_allclose(recovery.conj().T@recovery, np.eye(66), atol=3e-15, rtol=0)
                branch = np.zeros((66, 6), complex)
                for logical, row in enumerate(code.branch_images(law, seam, coin, raw)):
                    branch[[source_index[value] for value in row], logical] = 1/math.sqrt(12)
                expected = math.sqrt(float(code.branch_weight(law, seam, coin, raw)))*v
                np.testing.assert_allclose(recovery@branch, expected, atol=5e-16, rtol=0)


def test_every_source_route_and_branch_permutation_has_its_exact_complete_gray_word():
    # Source routing has no three-seam restriction; independently execute
    # its emitted integer gate schedule on ALL512 basis addresses.
    for seam in code.SEAMS:
        route = code._route_permutation(seam)
        assert _integer_permutation(routing.permutation_edges(route)) == route
        inverse = _integer_permutation(routing.permutation_edges(route, inverse=True))
        assert all(inverse[route[i]] == i for i in range(512))
        source_codes = set(circuit.SOURCE_CODES)
        assert all(route[i] == i for i in range(512) if i not in source_codes)
        for state in native.states():
            routed = circuit.decode(route[circuit.encode(state)])
            assert routed[0] == state[seam[0]] and routed[1] == state[seam[1]]
    for law in ('C', 'E'):
        for seam in code.SEAMS:
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    if code.branch_weight(law, seam, coin, raw):
                        perm = code._branch_permutation(law, seam, coin, raw)
                        assert _integer_permutation(routing.permutation_edges(perm)) == perm


def test_actual_all_seam_schedules_decode_independent_native_images_to_declared_logical_slots():
    report = code._compiled_layout_certificate()
    assert report['inverse_recovery_schedules_checked'] == 300
    assert report['independent_native_image_logical_slots_checked'] == 4320
    assert report['native_routes_checked_on_all78_states_and_padding'] == 30
    assert report['all_emitted_integer_schedules_match_independent_native_logical_layout']


def test_unexecuted_seam_wrong_but_unitary_recovery_fails_every_producer_gate(monkeypatch):
    original = code._branch_permutation

    def swapped_logical_columns(law, seam, coin, raw):
        images = list(original(law, seam, coin, raw))
        if (law, seam, coin, raw) == ('C', (0, 5), 0, 1):
            images[0], images[1] = images[1], images[0]
        return tuple(images)

    monkeypatch.setattr(code, '_branch_permutation', swapped_logical_columns)
    wrong = swapped_logical_columns('C', (0, 5), 0, 1)
    # This old consistency check passes: the wrong schedule is a valid
    # unitary and implements its own supplied target permutation exactly.
    assert _integer_permutation(routing.permutation_edges(wrong)) == wrong
    monkeypatch.setattr(code, 'code_matrix', lambda: pytest.fail('invalid schedule reached source execution'))
    # Simulate an already cached numerical budget: fresh semantic schedule
    # checking must still precede returning it to a public caller.
    monkeypatch.setattr(code, '_source_composition_report', lambda: {'cached': True})
    for producer in (code.exact_certificate, code.execution_report, code.source_composition_report):
        with pytest.raises(ArithmeticError, match='compiled recovery schedule'):
            producer()


@pytest.mark.parametrize('part', ('encoder', 'native_route'))
def test_wrong_but_unitary_compiled_encoder_or_native_route_is_rejected(part, monkeypatch):
    if part == 'encoder':
        original = code._encoder_permutation

        def wrong():
            images = list(original())
            images[0], images[1] = images[1], images[0]
            return tuple(images)

        monkeypatch.setattr(code, '_encoder_permutation', wrong)
    else:
        original = code._route_permutation

        def wrong(seam):
            images = list(original(seam))
            if seam == (0, 5):
                # Change the selected endpoint law, rather than merely
                # relabeling spectator coordinates undone by the inverse.
                first = circuit.encode((2,)+(0,)*11)
                second = circuit.encode((0,)*5+(2,)+(0,)*6)
                images[first], images[second] = images[second], images[first]
            return tuple(images)

        monkeypatch.setattr(code, '_route_permutation', wrong)
    with pytest.raises(ArithmeticError, match='compiled'):
        code.exact_certificate()


def test_source_encoder_and_decoder_preserve_complex_reference_columns_and_raw_rows():
    actual, expected = code.code_matrix(), code.ideal_code_matrix()
    np.testing.assert_allclose(actual, expected, atol=1e-9, rtol=0)
    rng = np.random.default_rng(61266)
    logical = rng.normal(size=(6, 3))+1j*rng.normal(size=(6, 3))
    encoded = code.encoding(logical)
    np.testing.assert_allclose(encoded, expected@logical, atol=3e-9, rtol=0)
    decoded = code.decoding(encoded)
    target = np.zeros((512, 3), complex)
    target[:6] = logical
    np.testing.assert_allclose(decoded, target, atol=6e-9, rtol=0)
    assert decoded.shape == (512, 3)
    for value in (actual, expected, encoded, decoded):
        with pytest.raises(ValueError):
            value.setflags(write=True)


@pytest.mark.parametrize('law,seam', (('C', (0, 1)), ('E', (2, 10))))
def test_source_all_raw_corrected_branch_columns_and_arbitrary_references(law, seam):
    branches = code.repair_branches(law, seam, code.code_matrix())
    assert len(branches) == 2*(1 << circuit.RECORD_BITS[law])
    rng = np.random.default_rng(60666)
    reference = rng.normal(size=(6, 4))+1j*rng.normal(size=(6, 4))
    gram = np.zeros((6, 6), complex)
    for branch in branches:
        probability = branch['code_branch_weight']
        expected = math.sqrt(float(probability))*code.ideal_code_matrix()
        assert branch['amplitude'].shape == (512, 6)
        np.testing.assert_allclose(branch['amplitude'], expected, atol=2e-9, rtol=0)
        np.testing.assert_allclose(branch['amplitude']@reference, expected@reference, atol=1e-8, rtol=0)
        gram += branch['amplitude'].conj().T@branch['amplitude']
        assert branch['recovery_applied'] == (not branch['abort'])
    np.testing.assert_allclose(gram, np.eye(6), atol=1e-8, rtol=0)


def test_source_recovery_is_unitary_on_unknown_offcode_padding_and_keeps_zero_records():
    rng = np.random.default_rng(512612)
    state = rng.normal(size=(512, 2))+1j*rng.normal(size=(512, 2))
    state /= np.linalg.norm(state)
    recovered = code.apply_recovery('E', (2, 10), 0, 5, state)
    assert recovered.shape == state.shape
    np.testing.assert_allclose(recovered.conj().T@recovered, state.conj().T@state,
                               atol=2e-9, rtol=0)
    # Valid doubled-occupation record has zero code weight; its complete
    # off-code action is identity, not zero or a code projection.
    assert code.branch_weight('E', (2, 10), 0, 1) == 0
    np.testing.assert_array_equal(code.apply_recovery('E', (2, 10), 0, 1, state), state)
    with pytest.raises(ValueError, match='halt'):
        code.apply_recovery('E', (2, 10), 0, 3, state)


def test_complete_unrecovered_raw_instrument_retains_unnormalized_padding_references():
    state = np.zeros((512, 3), complex)
    state[511, 0], state[400, 1], state[0, 2] = 2+.3j, -.7j, .6
    branches = code.repair_branches('E', (0, 1), state, recover=False)
    assert len(branches) == 16
    assert all(b['amplitude'].shape == state.shape and not b['recovery_applied'] for b in branches)
    gram = sum((b['amplitude'].conj().T@b['amplitude'] for b in branches),
               np.zeros((3, 3), complex))
    np.testing.assert_allclose(gram, state.conj().T@state, atol=1e-10, rtol=0)
    assert np.trace(gram).real > 4


def test_channel_labels_are_operationally_distinct_and_costs_match_existing_source_counter():
    assert code.branch_weight('C', (0, 1), 0, 0) == F(5, 12)
    assert code.branch_weight('E', (0, 1), 0, 0) == F(1, 3)
    assert code.branch_images('C', (0, 1), 0, 0) != code.branch_images('E', (0, 1), 0, 0)
    for permutation in (code._encoder_permutation(), code._route_permutation((2, 10))):
        assert code._permutation_resources(permutation) == routing.permutation_resources(permutation)


def test_uniform_support_source_words_and_resource_counts_cover_all_recovery_sizes():
    for count in (1, 4, 5, 6):
        word = compiler.compile_real_state([1/math.sqrt(count)]*count+[0.]*(8-count))
        actual = compiler.execute(word, 3)
        inverse = compiler.execute(compiler.inverse(word), 3)
        target = np.asarray([1/math.sqrt(count)]*count+[0.]*(8-count))
        np.testing.assert_allclose(actual[:, 0], target, atol=1e-12, rtol=0)
        np.testing.assert_allclose(inverse@actual, np.eye(8), atol=2e-12, rtol=0)
    encoding = code.encoding_resources()
    assert encoding['blank_code_qubits'] == 6
    assert encoding['proper_code_service_events'] == 2*encoding['native_pulses']
    for law in ('C', 'E'):
        for coin in (0, 1):
            for raw in range(1 << circuit.RECORD_BITS[law]):
                row = code.recovery_resources(law, (0, 1), coin, raw)
                assert row['blank_code_qubits'] == row['new_classical_records'] == 0
                assert row['archive_record_rewrites'] == 0
                assert row['proper_code_service_events'] == 2*row['native_pulses']
                assert bool(row['native_pulses']) == bool(code.branch_weight(law, (0, 1), coin, raw))


@pytest.mark.parametrize('bad', [True, np.ones(6, dtype=bool), [1., False, 0., 0., 0., 0.],
                                np.zeros(5), np.empty((6, 0)), [math.nan]*6,
                                np.ma.array(np.zeros(6), mask=False)])
def test_invalid_logical_inputs_fail_before_source_execution(bad, monkeypatch):
    monkeypatch.setattr(code, 'code_matrix', lambda: pytest.fail('invalid input executed source'))
    with pytest.raises(ValueError):
        code.encoding(bad)


@pytest.mark.parametrize('law,seam,coin,raw', [
    ('Q', (0, 1), 0, 0), ('C', (1, 0), 0, 0), ('C', [0, 1], 0, 0),
    ('C', (False, 1), 0, 0), ('C', (0, 2), 0, 0), ('C', (0, 1), True, 0),
    ('E', (0, 1), 0, False), ('C', (0, 1), 0, 4)])
def test_malformed_channel_labels_are_not_coerced(law, seam, coin, raw):
    for operation in (code.branch_images, code.branch_weight, code.recovery_resources):
        with pytest.raises(ValueError):
            operation(law, seam, coin, raw)


@pytest.mark.parametrize('bad', [np.zeros(511), np.empty((512, 0)),
                                np.ones(512, dtype=bool), [math.nan]*512])
def test_invalid_source_inputs_fail_before_word_execution(bad, monkeypatch):
    monkeypatch.setattr(code, '_apply_permutation', lambda *a, **k: pytest.fail('invalid source execution'))
    with pytest.raises(ValueError):
        code.apply_recovery('C', (0, 1), 0, 0, bad)
    with pytest.raises(ValueError):
        code.decoding(bad)


def test_every_source_schedule_is_covered_by_the_uniform_composition_budget():
    report = code.source_composition_report()
    assert report['numerical_source_composition_bound_not_interval_or_hardware_certificate']
    assert report['unrestricted_reference_extension_has_same_operator_bound']
    assert report['raw_and_stabilized_macro_conventions_both_covered']
    for direction, inverse in (('forward', False), ('inverse', True)):
        _, guard = routing._executed_macro(inverse)
        change = report[f'{direction}_full512_macro_polar_change']
        residual = report[f'{direction}_full512_stabilized_macro_operator_error']
        assert change == guard['projection_change']
        assert report[f'{direction}_full512_macro_operator_error'] == residual+change
        assert guard['dimension'] == 512 and change > 0
    assert 0 < report['encoder_operator_error_bound'] < 1e-7
    assert 0 < report['decoder_operator_error_bound'] < 1e-7
    assert set(report['mixing_angle_synthesis_diagnostics']) == {'1', '4', '5', '6'}
    for row in report['mixing_angle_synthesis_diagnostics'].values():
        assert row['literal_syndrome_wire_support_validated']
        assert row['literal_inverse_word_is_exact_reversed_negative_angle_word']
        assert row['requested_inverse_adjoint_error'] < 1e-10
        assert row['requested_inverse_uniform_to_zero_error'] < 1e-10
    for law in ('C', 'E'):
        row = report['laws'][law]
        assert row['all_thirty_seam_schedules_and_raw_recoveries_counted']
        assert 0 < row['maximum_corrected_attempt_operator_error'] < 1e-6
        counts = row['maximum_corrected_attempt_resources']
        assert counts['proper_code_service_events'] == 2*counts['native_pulses']
    report['laws']['C']['maximum_corrected_attempt_operator_error'] = 100
    assert code.source_composition_report()['laws']['C']['maximum_corrected_attempt_operator_error'] < 1e-6


@pytest.mark.parametrize('mutation', ('logical_wires', 'inverse_only', 'nominal_angle'))
def test_literal_mixing_words_reject_wrong_unitaries_before_cached_composition(mutation, monkeypatch):
    original = code._mixing_word

    def wrong(count, inverse=False):
        word = original(count, inverse)
        if count != 4:
            return word
        if mutation == 'logical_wires':
            return compiler.remap(word, (3, 4, 5, 0, 1, 2, 6, 7, 8))
        if mutation == 'inverse_only':
            return original(count, False)
        forward = original(count, False)
        forward.append(compiler.pulse('Y', .25, [3]))
        return compiler.inverse(forward) if inverse else forward

    # Independent full-register action exposes the unitarity-preserving
    # wire error on the actual uniform four-term syndrome input.
    if mutation == 'logical_wires':
        state = np.zeros(512, complex)
        state[[0, 8, 16, 24]] = .5
        actual = compiler.execute(wrong(4, True), 9, state)
        assert np.linalg.norm(actual-np.eye(512)[:, 0]) > 1.
    monkeypatch.setattr(code, '_mixing_word', wrong)
    monkeypatch.setattr(code, '_source_composition_report',
                        lambda *args: pytest.fail('invalid literal word reached cached composition'))
    with pytest.raises(ArithmeticError, match='mixing word|angle synthesis'):
        code.source_composition_report()

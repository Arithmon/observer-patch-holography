"""Independent complete-output and hostile checks for transported clocks."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import gzip
import json
import math
import sys

import numpy as np
import pytest

from . import history_clock as clock


STATES = clock.native.states(binary=True)
INDEX = {state: i for i, state in enumerate(STATES)}


def endpoint_maps(law, edge):
    """Direct endpoint balancing, independent of every production Kraus map."""
    groups = {}
    for coin in (0, 1):
        for i, state in enumerate(STATES):
            a, b = edge
            low, parity = divmod(state[a]+state[b], 2)
            out = list(state)
            out[a], out[b] = low+coin*parity, low+(1-coin)*parity
            current = state[a]-out[a]
            label = current if law == 'C' else state[a]-state[b]
            groups.setdefault((coin, label, current), []).append((i, INDEX[tuple(out)]))
    return tuple((key, tuple(pairs)) for key, pairs in sorted(groups.items()))


def independent_reference():
    weights = (F(2, 225), F(1, 36), F(2, 45))
    return np.array([float(weights[state[0]+state[1]]) for state in STATES])


def apply_partial(density, pairs, ref):
    out = np.zeros_like(density)
    src = [i*ref+r for i, _ in pairs for r in range(ref)]
    dst = [j*ref+r for _, j in pairs for r in range(ref)]
    out[np.ix_(dst, dst)] = density[np.ix_(src, src)]/2
    return out


def pair_rotation(pair, sign):
    out = np.eye(66, dtype=complex)
    out[np.ix_(pair, pair)] = np.array([[1, sign], [-sign, 1]])/math.sqrt(2)
    return out


def oracle(density, law, placement):
    """Independent ideal full-output instrument in native66 coordinates.

    Native maps and composed supports are obtained only by integer endpoint
    balancing. No production policy, support transport, ratios, clock, analyzer
    or native-step routine is called.
    """
    ref = len(density)//66
    pair = tuple(INDEX[tuple(int(p in occupied) for p in range(12))]
                 for occupied in ((0, 1), (0, 2)))
    sigma = independent_reference()
    r = np.kron(pair_rotation(pair, 1), np.eye(ref))
    work = r@density@r.conj().T
    if placement == 'before':
        phases = np.repeat(np.exp(-4j*np.log(60*sigma)), ref)
        work = phases[:, None]*work*phases.conj()[None, :]
    result = {}
    for (coin, label, current), first in endpoint_maps(law, (0, 5)):
        h1 = (0, 5, coin, label)
        edge = (0, 7) if current > 0 else (1, 8)
        middle = apply_partial(work, first, ref)
        for (coin2, label2, _), second in endpoint_maps(law, edge):
            h = (h1, (*edge, coin2, label2))
            final = apply_partial(middle, second, ref)
            follow = dict(second)
            composed = {i: follow[j] for i, j in first if j in follow}
            if placement == 'transported' and composed:
                ratios = np.ones(66)
                for i, j in composed.items():
                    ratios[j] = 60*sigma[i]
                phases = np.repeat(np.exp(-4j*np.log(ratios)), ref)
                final = phases[:, None]*final*phases.conj()[None, :]
            if all(i in composed for i in pair):
                after_pair = tuple(composed[i] for i in pair)
                analyzer = np.kron(pair_rotation(after_pair, -1), np.eye(ref))
                final = analyzer@final@analyzer.conj().T
            for (coin3, label3, current3), read in endpoint_maps(law, (1, 5)):
                result['read', h, coin3, label3, current3] = apply_partial(final, read, ref)
    return result


def assert_complete_outputs(actual, expected):
    reads = {key: value for key, value in actual.items() if key[0] == 'read'}
    assert set(reads) == set(expected)
    for key in reads:
        np.testing.assert_allclose(reads[key], expected[key], rtol=0, atol=3e-16)
    for key, value in actual.items():
        if key[0] == 'abort':
            np.testing.assert_array_equal(value, np.zeros_like(value))


@pytest.fixture(scope='module')
def ideal_outputs():
    density = np.diag(independent_reference()).astype(complex)
    return {law: {mode: clock.execute(density, law, mode=mode, raw=False)
                  for mode in ('before', 'transported', 'fixed', 'zero')}
            for law in ('C', 'E')}


@pytest.mark.parametrize('law', ['C', 'E'])
def test_independent_complete_mixed_history_and_distinguishing_controls(law, ideal_outputs):
    density = np.diag(independent_reference()).astype(complex)
    rows = ideal_outputs[law]
    for mode in ('before', 'transported'):
        assert_complete_outputs(rows[mode], oracle(density, law, mode))
    assert set(rows['before']) == set(rows['transported'])
    for key in rows['before']:
        np.testing.assert_allclose(rows['before'][key], rows['transported'][key], rtol=0, atol=4e-16)
    assert sum(np.trace(block).real for block in rows['transported'].values()) == pytest.approx(1, abs=2e-15)
    probabilities = {mode: {k: np.trace(v).real for k, v in blocks.items()} for mode, blocks in rows.items()}
    for wrong in ('fixed', 'zero'):
        distance = sum(abs(probabilities[wrong][k]-probabilities['transported'][k])
                       for k in probabilities[wrong])/2
        assert distance > 1e-5


@pytest.mark.parametrize('law', ['C', 'E'])
def test_independent_complex_entangled_history_keeps_both_reference_columns(law):
    rng = np.random.default_rng(1151)
    vector = rng.normal(size=132)+1j*rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    before = clock.execute(density, law, mode='before', raw=False)
    after = clock.execute(density, law, raw=False)
    assert_complete_outputs(after, oracle(density, law, 'transported'))
    for key in before:
        np.testing.assert_allclose(before[key], after[key], rtol=0, atol=4e-16)
    marginal = sum((block.reshape(66, 2, 66, 2).trace(axis1=0, axis2=2)
                    if key[0] == 'read' else block for key, block in after.items()), np.zeros((2, 2), complex))
    np.testing.assert_allclose(marginal, density.reshape(66, 2, 66, 2).trace(axis1=0, axis2=2),
                               rtol=0, atol=2e-15)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_faithful_extension_has_four_core_locations_and_preserves_every_corner(law):
    """Derive moved cores directly from endpoint currents, not port_permutation."""
    cores, expected_by_history = set(), {}
    for (coin, label, current), first in endpoint_maps(law, (0, 5)):
        core = {0, 1}
        if current:
            core = {5 if p == 0 else 0 if p == 5 else p for p in core}
        edge = (0, 7) if current > 0 else (1, 8)
        for (coin2, label2, current2), second in endpoint_maps(law, edge):
            following = dict(second)
            composed = {i: following[j] for i, j in first if j in following}
            if not composed:
                continue
            moved = ({edge[1] if p == edge[0] else edge[0] if p == edge[1] else p
                      for p in core} if current2 else core)
            cores.add(tuple(sorted(moved)))
            weights = (F(2, 225), F(1, 36), F(2, 45))
            expected = tuple(weights[sum(state[p] for p in moved)] for state in STATES)
            expected_by_history[((0, 5, coin, label), (*edge, coin2, label2))] = (expected, composed)
    assert cores == {(0, 1), (0, 8), (1, 5), (5, 8)}
    codes = clock.full.support()[0]
    words = {clock.ratios()}
    assert {b.history for b in clock.branches(law)} == set(expected_by_history)
    for branch in clock.branches(law):
        expected, composed = expected_by_history[branch.history]
        actual = clock.ratios(branch, law=law)
        words.add(actual)
        assert tuple(actual[c]/60 for c in codes) == expected
        assert all(actual[c] == 1 for c in set(range(512))-set(codes))
        for source, target in composed.items():
            assert expected[target] == (F(2, 225), F(1, 36), F(2, 45))[
                STATES[source][0]+STATES[source][1]]
    assert len(words) == 4


@pytest.mark.parametrize('law', ['C', 'E'])
def test_native_reference_is_unprobed_posterior_while_probe_retains_finite_coherence(law):
    sigma = np.diag(independent_reference()).astype(complex)
    pair = clock.probe_pair()
    rotation = pair_rotation(pair, 1)
    probed = rotation@sigma@rotation.conj().T
    admitted = {b.history: b for b in clock.branches(law)}
    coherent_branches = 0
    for (coin, label, current), first in endpoint_maps(law, (0, 5)):
        edge = (0, 7) if current > 0 else (1, 8)
        for (coin2, label2, _), second in endpoint_maps(law, edge):
            h = ((0, 5, coin, label), (*edge, coin2, label2))
            if h not in admitted:
                continue
            branch = admitted[h]
            unprobed = apply_partial(apply_partial(sigma, first, 1), second, 1)
            prepared = apply_partial(apply_partial(probed, first, 1), second, 1)
            reference = clock.transport.branch_reference(law, branch, clock.transport.reference())
            target = np.zeros((66, 66), complex)
            for i, weight in zip(reference.support, reference.unnormalized, strict=True):
                target[i, i] = float(weight)
            np.testing.assert_allclose(unprobed, target, rtol=0, atol=2e-18)
            if all(i in branch.inputs for i in pair):
                output = tuple(branch.outputs[branch.inputs.index(i)] for i in pair)
                assert unprobed[output] == 0
                assert abs(prepared[output]) == pytest.approx(float(F(1, 480)), rel=0, abs=2e-18)
                coherent_branches += 1
    assert coherent_branches > 0


def raw_key_oracle(law, word):
    def decode(edge, slot):
        coin, raw = slot
        labels = {}
        for (c, label, current), _ in endpoint_maps(law, edge):
            if c == coin:
                parity = abs(label) % 2 if law == 'E' else 0
                code = {0: 0, 1: 1, -1: 2}[current]+4*parity
                labels[code] = (label, current)
        return labels.get(raw)
    first = decode((0, 5), word[0])
    if first is None:
        return None
    edge = (0, 7) if first[1] > 0 else (1, 8)
    second, read = decode(edge, word[1]), decode((1, 5), word[2])
    if second is None or read is None:
        return None
    h = ((0, 5, word[0][0], first[0]), (*edge, word[1][0], second[0]))
    return 'read', h, word[2][0], read[0], read[1]


@pytest.mark.parametrize('law,raw_size,semantic_count', [('C', 4, 64), ('E', 8, 216)])
def test_full_raw_alphabet_lazy_zero_extension_and_missing_attained_blocks(law, raw_size, semantic_count, ideal_outputs):
    outputs = ideal_outputs[law]['transported']
    slots = tuple(product(range(2), range(raw_size)))
    counts, mass = 0, 0.
    for word in product(slots, repeat=3):
        expected = raw_key_oracle(law, word)
        assert clock.raw_read_key(law, word) == expected
        block = clock.raw_output(outputs, law, word)
        assert block.shape == (66, 66)
        if expected is None:
            np.testing.assert_array_equal(block, np.zeros((66, 66)))
        else:
            counts += 1
            np.testing.assert_array_equal(block, outputs[expected])
        mass += float(np.trace(block).real)
    assert counts == semantic_count
    assert mass == pytest.approx(1, abs=2e-15)
    certificate = clock.raw_alphabet_certificate(law, outputs)
    assert certificate['three_attempt_raw_words'] == (2*raw_size)**3
    assert certificate['attained_semantic_record_slots'] == semantic_count
    key = next(key for key in outputs if key[0] == 'read')
    word = next(word for word in product(slots, repeat=3) if raw_key_oracle(law, word) == key)
    missing = dict(outputs); del missing[key]
    with pytest.raises((KeyError, ValueError, ArithmeticError)):
        clock.raw_output(missing, law, word)
    with pytest.raises(ArithmeticError):
        clock.raw_alphabet_certificate(law, missing)
    detached = clock.raw_output(outputs, law, word)
    detached[0, 0] += 1
    assert detached[0, 0] != outputs[key][0, 0]


@pytest.mark.parametrize('damage', ['mapping', 'current', 'duplicate', 'coin_bool', 'map_float'])
def test_native_certificate_rejects_correct_matrices_with_wrong_metadata(monkeypatch, damage):
    maps = clock.interface.fixed_maps('C', clock.FIRST_EDGE)
    first = maps[0]
    if damage == 'mapping':
        pairs = list(first.mapping)
        pairs[0] = pairs[0][0], pairs[1][1]
        first = replace(first, mapping=tuple(pairs))
    elif damage == 'current':
        first = replace(first, current=1)
    elif damage == 'coin_bool':
        first = replace(first, coin=False)
    elif damage == 'map_float':
        first = replace(first, mapping=tuple((float(i), j) for i, j in first.mapping))
    damaged = (*maps, first) if damage == 'duplicate' else (first, *maps[1:])
    monkeypatch.setattr(clock.interface, 'fixed_maps', lambda law, edge: damaged)
    with pytest.raises((ArithmeticError, ValueError)):
        clock.certify_native('C', clock.FIRST_EDGE)


def test_native_certificate_rejects_nonfinite_matrix_and_charges_wrong_unitary_support(monkeypatch):
    maps = clock.interface.fixed_maps('E', clock.FIRST_EDGE)
    matrix = maps[0].matrix.copy(); matrix[0, 0] = float('nan')
    monkeypatch.setattr(clock.interface, 'fixed_maps', lambda law, edge: (replace(maps[0], matrix=matrix), *maps[1:]))
    with pytest.raises(ValueError):
        clock.certify_native('E', clock.FIRST_EDGE)


def test_full512_rotation_certificate_checks_wrong_sign_and_unused_padding():
    pair = clock.probe_pair()
    ideal = clock.rotation(pair, 1, False)
    assert clock.certify_rotation(pair, 1, ideal) < clock.ROTATION_GATE
    wrong = clock.rotation(pair, -1, False)
    assert clock.certify_rotation(pair, 1, wrong) > clock.ROTATION_GATE
    padding = next(i for i in range(512) if i not in clock.full.support()[0])
    wrong = ideal.copy(); wrong[padding, padding] = -1
    np.testing.assert_allclose(wrong.conj().T@wrong, np.eye(512), rtol=0, atol=3e-16)
    assert clock.certify_rotation(pair, 1, wrong) > 1
    wrong = ideal.copy(); wrong[padding, 0] = float('nan')
    with pytest.raises(ValueError):
        clock.certify_rotation(pair, 1, wrong)


@pytest.mark.parametrize('kwargs', [dict(law=True), dict(mode='wait'), dict(raw=1), dict(s=True),
    dict(s=4.), dict(s=F(17)), dict(raw=True, s=F(3))])
def test_execute_rejects_wrong_original_domains(kwargs):
    density = np.eye(66)/66
    args = {'law': 'C', 'raw': False, **kwargs}
    with pytest.raises(ValueError):
        clock.execute(density, **args)


@pytest.mark.parametrize('density', [np.eye(65)/65, np.eye(198)/198,
    np.eye(66)*float('nan'), -np.eye(66)/66, np.eye(66), np.eye(66).astype(int)])
def test_execute_rejects_wrong_or_nonfinite_density(density):
    with pytest.raises(ValueError):
        clock.execute(density, 'C', raw=False)


@pytest.mark.parametrize('word', [(), ((False, 0),)*3, ((0, True),)*3, ((0, -1),)*3,
    ((0, 4),)*3, [(0, 0)]*3])
def test_raw_word_validation(word):
    with pytest.raises(ValueError):
        clock.raw_read_key('C', word)


@pytest.mark.parametrize('numerator,denominator', [(True, 2), (1, True), (2, 4), (1, -2), (1, 0)])
def test_parameter_rejects_malformed_fraction_slots(numerator, denominator):
    value = object.__new__(F)
    value._numerator, value._denominator = numerator, denominator
    with pytest.raises(ValueError):
        clock.clock(value, clock.ratios(), False)


def test_clock_and_rotation_domain_and_immutable_cached_source_outputs():
    pair = clock.probe_pair()
    for bad_pair, bad_sign in ((pair, True), ((True, 2), 1), ((0, 0), 1), ((0, 66), 1)):
        with pytest.raises(ValueError):
            clock.rotation(bad_pair, bad_sign, False)
    with pytest.raises(ValueError):
        clock.clock(F(4), list(clock.ratios()), False)
    with pytest.raises(ValueError):
        clock.clock(F(4), (F(0),)*512, False)
    values, parameter = clock.ratios(), F(4)
    actual = clock.clock(parameter, values, True)
    saved = actual.copy()
    with pytest.raises(ValueError):
        actual.flags.writeable = True
    parameter._numerator = 3
    values[clock.full.support()[0][0]]._numerator += 1
    np.testing.assert_array_equal(clock.clock(F(4), clock.ratios(), True), saved)
    rotation = clock.rotation(pair, 1, True)
    with pytest.raises(ValueError):
        rotation.flags.writeable = True


def receipt_structure():
    """Small custody fixture; the packet fixture separately executes the source."""
    return {
        'schema': clock.SCHEMA, 'parameter': '4',
        'source_sha256': {'model/native_modular_reference/history_clock.py': 'a'*64},
        'source_budget': {'complete_history_halfdiamond_upper': '1/1000000'},
        'complete_preparation_history_halftrace_upper': '1/1000000',
        'histories': [{'law': 'C', 'complete_prepared_mass_diagnostic': 1.,
            'records': [{'record': ['preparation', 'herald_failure'], 'probability': 1.}],
            'wrong_fixed_reference_record_variation_diagnostic': 0.001,
            'nonzero_clock_record_variation_diagnostic': 0.002,
            **{key: {'lower': '0', 'upper': clock.NUMERICAL_RESIDUAL_CEILING}
               for key in clock.RESIDUAL_INTERVAL_FIELDS}}]}


def test_cli_gzip_is_deterministic_roundtrips_and_rejects_corruption(monkeypatch, tmp_path, capsys):
    representative = receipt_structure()
    monkeypatch.setattr(clock, 'run', lambda: deepcopy(representative))
    first, second = tmp_path/'first.json.gz', tmp_path/'different-name.json.gz'
    for target in (first, second):
        monkeypatch.setattr(sys, 'argv', ['history_clock', '--write', str(target)])
        clock.main()
    encoded = first.read_bytes()
    assert encoded == second.read_bytes()
    assert encoded[:2] == b'\x1f\x8b'
    assert encoded[3] & 8 == 0  # No embedded output filename.
    assert encoded[4:8] == bytes(4)  # Fixed mtime.
    assert gzip.decompress(encoded).decode('utf-8') == clock.full.corner.canonical(representative)
    monkeypatch.setattr(sys, 'argv', ['history_clock', '--check', str(first)])
    clock.main()
    assert 'verified' in capsys.readouterr().out
    first.write_bytes(b'not a gzip archive')
    with pytest.raises(OSError):
        clock.main()
    first.write_bytes(gzip.compress(b'{"schema":"first","schema":"duplicate"}', mtime=0))
    with pytest.raises(ValueError, match='duplicate'):
        clock.main()


@pytest.mark.parametrize('key', ['source', 'preparation'])
@pytest.mark.parametrize('bound', ['0', '-1', '1/10000000000', '1/0', 'not-a-bound', 1e-6, True])
def test_identical_receipts_reject_budgets_incompatible_with_retained_intervals(key, bound):
    bad = receipt_structure()
    if key == 'source':
        bad['source_budget']['complete_history_halfdiamond_upper'] = bound
    else:
        bad['complete_preparation_history_halftrace_upper'] = bound
    with pytest.raises(ValueError, match='scientific bound'):
        clock.compare_receipt(bad, deepcopy(bad))
    valid = receipt_structure()
    valid['source_budget']['complete_history_halfdiamond_upper'] = clock.NUMERICAL_RESIDUAL_CEILING
    valid['complete_preparation_history_halftrace_upper'] = clock.NUMERICAL_RESIDUAL_CEILING
    clock.compare_receipt(valid, deepcopy(valid))


@pytest.fixture(scope='module')
def packet():
    # One complete actual source replay; no monkeypatched shortcuts or reused
    # stored producer report supplies the reference answer.
    return json.loads(clock.full.corner.canonical(clock.run()))


def test_complete_source_producer_and_receipt_roundtrip(packet, tmp_path):
    target = tmp_path/'history.json'
    target.write_text(json.dumps(packet, allow_nan=False), encoding='utf-8')
    clock.compare_receipt(json.loads(target.read_text(encoding='utf-8')), packet)
    assert len(packet['histories']) == 2
    for row in packet['histories']:
        assert row['wrong_fixed_reference_record_variation_diagnostic'] > 1e-5
        assert row['nonzero_clock_record_variation_diagnostic'] > 1e-5
        assert abs(sum(record['probability'] for record in row['records'])-1) < 2e-10
        for field in clock.RESIDUAL_INTERVAL_FIELDS:
            assert row[field] == {'lower': '0', 'upper': '1/5000000000'}
    budget = packet['source_budget']
    assert budget['distinct_clock_words_verified'] == 4
    assert budget['reference_core_location_bits_upper'] == 8
    for row in packet['histories']:
        assert row['reference_extension_cores'] == [[0, 1], [0, 8], [1, 5], [5, 8]]
    assert budget['maximum_coherent_load_export_service_calls'] == 2*budget['maximum_source_pulses_per_successful_history']
    assert budget['native_repair_attempts_after_preparation'] == 2
    assert budget['native_detector_attempts'] == 1
    assert budget['native_export_checks'] == 3
    assert packet['scope']['floating_per_pulse_apparatus_error_bound_claimed'] is False
    assert packet['scope']['faithful_extension_is_a_port_permutation_of_native_reference'] is True
    assert packet['scope']['extension_is_executed_routing_or_native_graph_symmetry'] is False


@pytest.mark.parametrize('damage', ['hash', 'budget', 'missing_record', 'zero_signal', 'nan', 'negative_probability'])
def test_receipt_rejects_forged_bounds_custody_and_complete_outputs(packet, damage):
    bad = deepcopy(packet)
    if damage == 'hash':
        bad['source_sha256'][next(iter(bad['source_sha256']))] = '0'*64
    elif damage == 'budget':
        bad['source_budget']['complete_history_halfdiamond_upper'] = '0'
    elif damage == 'missing_record':
        bad['histories'][0]['records'].pop()
    elif damage == 'zero_signal':
        bad['histories'][0]['nonzero_clock_record_variation_diagnostic'] = 0.
    elif damage == 'nan':
        bad['histories'][0]['nonzero_clock_record_variation_diagnostic'] = float('nan')
    else:
        bad['histories'][0]['records'][0]['probability'] = -1e-100
    with pytest.raises(ValueError):
        clock.compare_receipt(bad, packet)


@pytest.mark.parametrize('field', clock.RESIDUAL_INTERVAL_FIELDS)
@pytest.mark.parametrize('replacement', [0., {'lower': '0', 'upper': '0'},
    {'lower': '0', 'upper': '1/10000000000'}, {'lower': '0', 'upper': '1/1000000000'},
    {'lower': '0', 'upper': 2e-10}])
def test_residual_intervals_cannot_claim_false_zero_or_changed_ceiling(packet, field, replacement):
    bad = deepcopy(packet)
    bad['histories'][0][field] = deepcopy(replacement)
    with pytest.raises(ValueError):
        clock.compare_receipt(bad, packet)
    # Even coordinated equal reports cannot change the canonical admission.
    with pytest.raises(ValueError):
        clock.compare_receipt(bad, deepcopy(bad))


@pytest.mark.parametrize('residual_index', range(5))
@pytest.mark.parametrize('value', [math.nextafter(2e-10, math.inf), -1e-30, float('nan')])
def test_producer_checks_each_original_residual_before_interval_projection(monkeypatch, residual_index, value):
    # Isolate the admission path: source execution itself has the separate
    # real packet replay above. No fabricated fixture is retained as evidence.
    monkeypatch.setattr(clock, 'source_budget', lambda: {
        'complete_history_halfdiamond_upper': '1/1000000', 'channel_norm_upper': '1'})
    monkeypatch.setattr(clock, '_joined_budget', lambda budget: ({}, F(1, 10**6)))
    monkeypatch.setattr(clock, 'execute', lambda density, *a, **k: {('terminal',): density})
    monkeypatch.setattr(clock, 'raw_alphabet_certificate', lambda *a: {})
    monkeypatch.setattr(clock, '_variation', lambda *a: 0.001)
    monkeypatch.setattr(clock.preparation, 'population_maps', lambda **k: {'success': np.eye(66)})
    monkeypatch.setattr(clock, 'prepared_target', lambda *a: {})
    values = [0.]*5
    values[residual_index] = value
    observations = iter(values)
    monkeypatch.setattr(clock.previous, 'distance', lambda *a: next(observations))
    with pytest.raises(ArithmeticError, match='fresh clock residual'):
        clock.run()


def test_residual_projection_checks_the_separate_scientific_bound_and_detaches():
    with pytest.raises(ArithmeticError, match='scientific bound'):
        clock._residual_interval(1e-12, F(1, 10**13))
    retained = clock._residual_interval(1e-12, F(1, 10**8))
    retained['upper'] = '0'
    assert clock._residual_interval(0., F(1, 10**8)) == {
        'lower': '0', 'upper': '1/5000000000'}


def test_probability_replay_limits_total_variation_not_only_each_slot(packet):
    bad = deepcopy(packet)
    rows = bad['histories'][0]['records']
    positive = [i for i, row in enumerate(rows) if row['probability'] > 1e-8][:8]
    assert len(positive) == 8
    for k, i in enumerate(positive):
        rows[i]['probability'] += 1e-10 if k % 2 else -1e-10
    # The old componentwise comparator accepts this mass-preserving mutation.
    clock.full.corner.compare_receipt(bad, packet)
    expected = packet['histories'][0]['records']
    assert math.fsum(abs(a['probability']-b['probability'])
                     for a, b in zip(rows, expected, strict=True))/2 > 3.9e-10
    with pytest.raises(ValueError, match='total-variation'):
        clock.compare_receipt(bad, packet)
    # A small portable probability difference remains admitted.
    near = deepcopy(packet)
    near_rows = near['histories'][0]['records']
    near_rows[positive[0]]['probability'] += 1e-11
    near_rows[positive[1]]['probability'] -= 1e-11
    clock.compare_receipt(near, packet)

"""Hostile source, quantum-reference and complete-history preparation audits."""
from collections import defaultdict
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, repeat
import json

import numpy as np
import pytest

from . import preparation as p


def test_native_target_selector_has_exact_uniform_subnormalized_columns():
    words = tuple(product((0, 1), repeat=7))
    addresses = [p.coins.decode_coin_word(word, 66) for word in words]
    assert addresses.count(None) == 62
    assert all(addresses.count(t) == 1 for t in range(66))
    selector = p.coins.selector_certificate(66, 32)
    assert F(selector['timeout_probability']) == F(31, 64)**32
    assert 66*F(selector['accepted_probability_per_address'])+F(selector['timeout_probability']) == 1


def test_all_5148_reset_columns_are_exact_complete_and_label_correlated():
    table = p.preparation_table()
    assert len(table.inputs) == 78 and len(table.targets) == 66
    per_target = (1-table.selector_timeout)/66
    assert all(len(row) == 66 and all(type(v) is F and 0 <= v <= per_target for v in row)
               for row in table.probabilities)
    bound = F(p.preparation_certificate()['complete_halfdiamond_upper'])
    assert 0 < max(table.failure_by_input) <= bound < F(1, 10**7)
    assert p.preparation_certificate()['ideal_initial_label_probabilities'] == {
        '00': '15/22', '01': '5/33', '10': '5/33', '11': '1/66'}
    # Exact uniformity after success conditioning is false, even though the
    # unconditioned labelled replacement is accurate and reference-complete.
    assert any(len(set(row)) > 1 for row in table.probabilities)
    with pytest.raises(TypeError):
        table.probabilities[0][0] = F(0)
    with pytest.raises(AttributeError):
        table.source_digest = 'fake'


def test_arbitrary_entangled_reference_complete_deficit_not_conditioned_success():
    table = p.preparation_table(1, 2, 2)
    failure = table.failure_by_input
    pair = min(range(78), key=failure.__getitem__), max(range(78), key=failure.__getitem__)
    vector = np.zeros(156, complex)
    vector[2*pair[0]], vector[2*pair[1]+1] = .6, .8j
    density = np.outer(vector, vector.conj())
    output = p.apply_preparation(density, 1, 2, 2)
    ideal = output['reference']/66
    success_deficit = np.array([ideal-block for block in output['success']])
    assert min(np.linalg.eigvalsh(block).min() for block in success_deficit) >= -1e-15
    expected_failure = .36*float(failure[pair[0]])+.64*float(failure[pair[1]])
    halftrace = (sum(np.linalg.norm(v, 'nuc') for v in success_deficit)
                 +np.linalg.norm(output['target_timeout'], 'nuc')
                 +np.linalg.norm(output['reset_timeout'], 'nuc'))/2
    assert abs(halftrace-expected_failure) < 1e-13
    np.testing.assert_allclose(output['success'].sum(axis=0)+output['target_timeout']+output['reset_timeout'],
                               output['reference'], rtol=0, atol=1e-14)
    # Renormalized success changes the reference. The theorem does not silently
    # replace it by the original marginal or retain the private resolved x.
    conditioned = output['success'].sum(axis=0)
    conditioned /= np.trace(conditioned)
    assert np.linalg.norm(conditioned-output['reference']) > 1e-5


@pytest.mark.parametrize('law', ('C', 'E'))
def test_source_coarse_coin_identity_with_correlated_reuse_and_adaptive_work(law):
    assert p.coins.effect_certificate()['coarse_coin_effect'] == 'I/2'
    result = p.correlated_coin_audit(law)
    assert result['complete_coin_words'] == 8
    assert result['array_reduced_effect_audit_pass']
    assert F(result['array_roundoff_allowance']) == F(64, 2**52)
    # Checking a current subbranch instead of the coarse coin fails: its
    # effect does not contain all source columns and it can steer R.
    maps = p.coins.coin_maps(law, (0, 1), False)
    assert any(0 < len(mapping) < 78 for _, _, mapping in maps)


@lru_cache(maxsize=2)
def _independent_stopped_columns(cap, rounds):
    states = tuple(p.native.states(binary=True))
    seams = tuple(p.reset.quantum.seams())
    alpha, tau = 1-F(1, 16**rounds), F(1, 16**rounds)
    result = []
    for original in states:
        waiting, stopped, rejected = {original: F(1)}, defaultdict(F), defaultdict(F)
        for _ in range(cap):
            following = defaultdict(F)
            for state, mass in waiting.items():
                rejected[state] += mass*tau
                for seam, coin in product(seams, (0, 1)):
                    final, _ = p.native.native(state, seam, coin)
                    target = stopped if 0 in seam or 1 in seam else following
                    target[final] += mass*alpha/60
            waiting = following
        result.append((stopped, rejected, waiting))
    return result


@pytest.mark.parametrize('cap', (1, 2))
def test_stopped_complete_kernels_match_independent_native_path_enumeration(cap):
    kernel = p.stopped_kernel(cap, 1)
    independent = _independent_stopped_columns(cap, 1)
    for column, branches in enumerate(independent):
        for matrix, actual in zip((kernel.success, kernel.selector_timeout, kernel.incident_timeout), branches, strict=True):
            assert all(F(matrix[row][column], kernel.denominator) == actual.get(state, 0)
                       for row, state in enumerate(p.core.source().states))
    assert all(sum(matrix[row][column] for matrix in
                   (kernel.success, kernel.selector_timeout, kernel.incident_timeout) for row in range(66))
               == kernel.denominator for column in range(66))
    if cap == 2:
        # A B ordering would alter terminal populations for selected inputs.
        a, b = p.core.source().remote_counts, p.core.source().incident_counts
        assert p._left(a, b) != p._left(b, a)


def _independent_short_fields(initial_index, feedback):
    table = p.preparation_table(1, 2, 2)
    initial = table.probabilities[initial_index]
    stopped = _independent_stopped_columns(2, 1)
    indices = {state: index for index, state in enumerate(table.targets)}
    branches = {(tuple(state[:2]), ()): defaultdict(F) for state in table.targets}
    for state, mass in zip(table.targets, initial, strict=True):
        branches[tuple(state[:2]), ()][state] += mass
    terminal = defaultdict(lambda: defaultdict(F))
    for cursor, port in enumerate((0, 1)):
        following = defaultdict(lambda: defaultdict(F))
        for (label, word), population in branches.items():
            for state, weight in population.items():
                stopped_branch = stopped[indices[state]]
                for flag, output in zip(('selector_timeout', 'incident_timeout'), stopped_branch[1:], strict=True):
                    for final, mass in output.items():
                        terminal[flag, label, word][final] += weight*mass
                for final, mass in stopped_branch[0].items():
                    bit = final[port]
                    if feedback and cursor == 0 and bit:
                        for coin in (0, 1):
                            changed, _ = p.native.native(final, (0, 1), coin)
                            following[label, (*word, bit)][changed] += weight*mass/2
                    else:
                        following[label, (*word, bit)][final] += weight*mass
        branches = following
    for (label, word), values in branches.items():
        terminal['complete', label, word] = values
    return terminal


@pytest.mark.parametrize('feedback', (False, True))
def test_short_complete_field_histories_match_independent_endpoint_oracle(feedback):
    law = p.finite_field_law(7, feedback=feedback, seed_scans=1, retries=2,
                             target_rounds=2, selector_rounds=1, cap=2, read_ports=(0, 1))
    independent = _independent_short_fields(7, feedback)
    assert sum(map(sum, law.counts)) == law.denominator
    assert law.records == p.record_alphabet((0, 1))
    for key, counts in zip(law.records, law.counts, strict=True):
        if key[0] in ('target_timeout', 'reset_timeout'):
            assert len(counts) == 1  # explicitly discarded work, not a fabricated state
        else:
            assert len(counts) == 66
            assert all(F(value, law.denominator) == independent[key].get(state, 0)
                       for value, state in zip(counts, p.core.source().states, strict=True))
    assert any(v == 0 for row in law.counts for v in row)  # retain impossible terminal states


def _prefix_to_field(seed=20261025):
    initial, auxiliary = p.native.states()[0], p.native.states()[3]
    history = p.execute_history(initial, auxiliary, p._bits(seed))
    assert history.flag == 'complete'
    first = next(i for i, row in enumerate(history.native_records) if row[1] == 'schedule')
    return initial, auxiliary, history, tuple(row[3] for row in history.native_records[:first])


def test_microexecution_reuses_auxiliary_keeps_every_timeout_and_true_terminal_state():
    initial, auxiliary, history, prefix = _prefix_to_field()
    assert p.replay_history(history, initial, auxiliary)
    target_timeout = p.execute_history(initial, auxiliary, repeat(1), target_rounds=2)
    assert target_timeout.flag == 'target_timeout' and target_timeout.native_attempts == 14
    assert target_timeout.terminal_work == initial and target_timeout.initial_label is None
    selector_timeout = p.execute_history(initial, auxiliary, (*prefix, *([1]*5)), selector_rounds=1)
    assert selector_timeout.flag == 'selector_timeout' and selector_timeout.word == ()
    remote = p.core.source().remote_edges[0]
    remote_bits = tuple((remote >> shift) & 1 for shift in reversed(range(5)))
    waiting_timeout = p.execute_history(initial, auxiliary, (*prefix, *remote_bits, 0), cap=1)
    assert waiting_timeout.flag == 'incident_timeout' and waiting_timeout.word == ()
    assert waiting_timeout.native_attempts == len(prefix)+6
    # No fresh tracial auxiliary is inserted at the boundary or each schedule.
    actual_aux = auxiliary
    for carrier, _, seam, coin, _ in history.native_records:
        if carrier == 'auxiliary':
            actual_aux, _ = p.native.native(actual_aux, seam, coin)
    assert history.terminal_auxiliary == actual_aux
    for law in ('C', 'E'):
        failed = p.execute_history(initial, auxiliary, repeat(0), law=law, seed_scans=1, retries=1)
        assert failed.flag == 'reset_timeout'
        assert failed.reset_reason in ('seed_timeout', 'flip_timeout', 'probe_timeout', 'route_timeout')
        assert failed.terminal_work in p.native.states(binary=True)


def test_zero_current_incident_attempts_count_and_feedback_consumes_a_real_tick():
    initial, auxiliary, _, prefix = _prefix_to_field()
    # Select seam0=(0,1), retain all four incident attempts including zeros.
    history = p.execute_history(initial, auxiliary, (*prefix, *([0]*24)))
    assert history.flag == 'complete' and history.read_clock == (1, 2, 3, 4)
    scheduled = [row for row in history.native_records if row[1] == 'scheduled']
    assert len(scheduled) == 4 and any(row[4] == 0 for row in scheduled)
    for seed in range(30):
        changed = p.execute_history(initial, auxiliary, p._bits(seed), feedback=True)
        if changed.flag == 'complete' and changed.word[0]:
            break
    assert changed.word[0] == 1
    assert changed.read_clock == (1, 3, 4, 5)
    assert sum(row[1] == 'feedback' for row in changed.native_records) == 1


def test_microscopic_fine_records_and_private_terminal_states_cannot_be_forged():
    initial, auxiliary, history, _ = _prefix_to_field()
    bad_records = list(history.native_records)
    carrier, phase, seam, coin, mark = bad_records[0]
    bad_records[0] = carrier, phase, seam, coin, mark+1
    mutations = [replace(history, native_records=tuple(bad_records)),
                 replace(history, read_clock=(0, 1, 2, 3)),
                 replace(history, initial_label=(1, 1)),
                 replace(history, terminal_auxiliary=p.native.states()[10]),
                 replace(history, native_records=history.native_records[:-1])]
    for bad in mutations:
        with pytest.raises(ValueError):
            p.replay_history(bad, initial, auxiliary)


@pytest.fixture(scope='module')
def receipt():
    return p.run()


def test_default_complete_source_field_connection_is_nonvacuous(receipt):
    bound = F(receipt['resources']['complete_field_record_TV_upper'])
    assert bound < F(9, 10**9)
    assert receipt['resources']['all_native_attempts_max'] == 13399
    assert F(receipt['feedback_separation']['uniform_all_input_TV_separation_lower']) > F(1, 20)
    assert receipt['feedback_separation']['existing_ideal_TV'] == '5/99'
    for row in receipt['finite_field_laws']:
        assert row['records_including_impossible_zero_branches'] == 186
        assert row['mass_exactly_one']
        assert F(row['complete_success_probability_lower']) > 1-F(1, 10**7)
        assert F(row['complete_record_TV_to_existing_ideal_upper']) <= bound
    assert len(receipt['microscopic_native_executions']) == 4
    assert all(row['full78_inputs_executed'] == 78 for row in receipt['microscopic_native_executions'])
    kernel = p.stopped_kernel()
    z = p._completion_probability(p.CAP, p.SELECTOR_ROUNDS)
    assert all(F(sum(row[column] for row in kernel.success), kernel.denominator) == z
               for column in range(66))
    exact = max(1-(1-f)*z**4 for f in p.preparation_table().failure_by_input)
    assert F(receipt['resources']['exact_native_continuation_halfdiamond_enclosure'][0]) <= exact <= bound
    assert bound <= F(receipt['resources']['complete_field_record_TV_union_upper'])


@pytest.mark.parametrize('text', ('{"schema": "first", "schema": "second"}', '{"all_pass": NaN}'))
def test_cli_receipt_loader_rejects_duplicate_keys_and_nonfinite_json(receipt, monkeypatch, tmp_path, text):
    path = tmp_path/'hostile.json'
    path.write_text(text, encoding='utf-8')
    monkeypatch.setattr(p, 'run', lambda: deepcopy(receipt))
    monkeypatch.setattr('sys.argv', ['preparation', '--check', str(path)])
    with pytest.raises((ValueError, SystemExit)):
        p.main()


@pytest.mark.parametrize('mutate', (
    lambda r: r.update(schema='wrong'),
    lambda r: r.update(all_pass=1),
    lambda r: r['parameters'].update(retries=True),
    lambda r: r['resources'].update(complete_field_record_TV_upper='0'),
    lambda r: r['preparation'].update(input_columns=66),
    lambda r: r['preparation'].update(source_success_table_sha256='0'*64),
    lambda r: r['scope'].update(physical_time_or_full_federation_claimed=True),
    lambda r: r.update(microscopic_native_executions=[]),
    lambda r: r.update(finite_field_laws=[]),
    lambda r: r['finite_field_laws'][0]['successful_record_probability_enclosures'].pop(),
    lambda r: r['finite_field_laws'][0].update(complete_success_probability_lower='1'),
    lambda r: r['source_sha256'].clear(),
))
def test_receipt_comparator_rejects_self_forgery(receipt, monkeypatch, mutate):
    # Canonical source replay was executed once by the fixture. A supplied
    # alleged actual result must never authorize its own matching forgery.
    monkeypatch.setattr(p, 'run', lambda: deepcopy(receipt))
    bad = deepcopy(receipt)
    mutate(bad)
    with pytest.raises(ValueError):
        p.compare_receipt(bad, bad)


def test_receipt_exact_types_nan_and_aliases_fail_closed(receipt, monkeypatch):
    monkeypatch.setattr(p, 'run', lambda: deepcopy(receipt))
    assert p.compare_receipt(deepcopy(receipt), deepcopy(receipt))
    bad = deepcopy(receipt)
    bad['reference_execution']['complete_failure_probability'] = float('nan')
    with pytest.raises(ValueError):
        p.compare_receipt(bad)
    bad = deepcopy(receipt)
    bad['parameters']['seed_scans'] = 4.0
    with pytest.raises(ValueError):
        p.compare_receipt(bad)


@pytest.mark.parametrize('kwargs', (
    {'initial_index': True}, {'initial_index': 78}, {'feedback': 1}, {'cap': 0},
    {'cap': 65}, {'selector_rounds': True}, {'selector_rounds': 9},
    {'read_ports': [0, 1]}, {'read_ports': (True,)}, {'read_ports': (2,)},
    {'read_ports': ()}, {'target_rounds': 33}, {'retries': 0}, {'seed_scans': 1.0},
))
def test_exact_law_contract_rejects_ambiguous_or_unsupported_inputs(kwargs):
    with pytest.raises(ValueError):
        p.finite_field_law(**kwargs)


@pytest.mark.parametrize('density', (np.eye(66), np.full((78, 78), np.nan),
                                   np.eye(78, dtype=int), [[1]], np.eye(312)))
def test_quantum_reference_execution_refuses_invalid_arrays(density):
    with pytest.raises(ValueError):
        p.apply_preparation(density)


def test_source_fault_does_not_pass_exact_kraus_audit(monkeypatch):
    original = p.native.native
    def bad(state, seam, coin):
        final, (c, current) = original(state, seam, coin)
        return final, (c, -current)
    monkeypatch.setattr(p.native, 'native', bad)
    with pytest.raises(ArithmeticError):
        p.reset.source_audit()

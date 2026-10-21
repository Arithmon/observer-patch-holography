"""Complete-record scheduler histories and strict evidence boundaries."""
import copy
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from model.closed_prediction_protocol import construction as old
from model.native_interface import construction as interface
from model.observer import native
from model.observer.scalar_response.experiment import compare_receipt
from . import construction as build


# Independent spectral and Kraus implementations may sum in different orders.
# This remains four orders below the complete-mass execution guard.
ORACLE_ATOL = 1e-14


def _reference_input():
    left, right = interface.phase_pairs()[0]
    vector = np.zeros(132, complex)
    vector[2*left], vector[2*right+1] = .6, .8j
    return np.outer(vector, vector.conj())


def _full_input():
    states = native.states()
    collision = next(i for i, x in enumerate(states) if max(x) == 2)
    binary = build._binary_indices()[0]
    vector = np.zeros(78, complex)
    vector[collision], vector[binary] = .6, .8j
    return np.outer(vector, vector.conj())


def _oracle(density, law, attempts, rounds):
    """Existing unbounded program and independent sparse prefix propagation."""
    q = F(16)**(-rounds)
    program = interface.program(law, attempts)
    unbounded = interface.oracle(program, law, attempts, initial=density)
    result = {(0, 0, bit): unbounded[0, bit] for bit in (0, 1)}
    ref = len(density)//66
    herald, _ = interface.herald_and_pointer(law)
    mixer = np.kron(interface.mixer(), np.eye(ref))
    for h in (1, 2):
        selected = sum((np.kron(k, np.eye(ref))@density@np.kron(k.conj().T, np.eye(ref))
                        for k in herald.kraus[h]), np.zeros_like(density))
        prefix = mixer@selected@mixer.conj().T
        for stage in range(1, attempts+1):
            # A timeout retains the prefix before that work attempt.
            result[h, stage, 0] = float(q*(1-q)**(stage-1))*old.prepare(prefix, law, stage-1)
        for bit in (0, 1):
            result[h, attempts+1, bit] = float((1-q)**attempts)*unbounded[h, bit]
    return result


def test_default_complete_alphabet_retains_all_twelve_records():
    expected = {(0, 0, 0), (0, 0, 1)}
    expected |= set(product((1, 2), (1, 2, 3), (0,)))
    expected |= set(product((1, 2), (4,), (0, 1)))
    assert len(build.public_records()) == 12
    assert set(build.public_records()) == expected


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('attempts, rounds', [(1, 1), (3, 1), (3, 7)])
def test_complete_outputs_match_independent_unbounded_and_timeout_prefix_laws(law, attempts, rounds):
    density = _reference_input()
    actual = build.trial_outputs(density, law, attempts=attempts, rounds=rounds)
    expected = _oracle(density, law, attempts, rounds)
    assert actual.keys() == expected.keys() == set(build.public_records(attempts))
    for record in actual:
        np.testing.assert_allclose(actual[record], expected[record], rtol=0, atol=ORACLE_ATOL)
    np.testing.assert_array_equal(actual[0, 0, 1], np.zeros_like(density))
    reduced = sum(np.trace(value.reshape(66, 2, 66, 2), axis1=0, axis2=2)
                  for value in actual.values())
    np.testing.assert_allclose(reduced,
        np.trace(density.reshape(66, 2, 66, 2), axis1=0, axis2=2), rtol=0, atol=1e-15)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_timeout_output_is_before_unperformed_work_and_cannot_be_reset(law):
    density = _reference_input()
    expected = _oracle(density, law, 3, 1)
    actual = build.trial_outputs(density, law, attempts=3, rounds=1)
    distances = []
    for h in (1, 2):
        timeout = expected[h, 1, 0]
        incorrectly_advanced = old.prepare(timeout, law, 1)
        distances.append(np.linalg.norm(timeout-incorrectly_advanced, 'nuc')/2)
        np.testing.assert_allclose(actual[h, 1, 0], timeout, rtol=0, atol=ORACLE_ATOL)
    assert max(distances) > 1e-5


@pytest.mark.parametrize('law', ['C', 'E'])
def test_uniform_step_preserves_entangled_reference_and_matches_existing_full_channel(law):
    density = _reference_input()
    actual = build.uniform_native_step(density, law)
    np.testing.assert_allclose(actual, old.prepare(density, law, 1), rtol=0, atol=ORACLE_ATOL)
    np.testing.assert_allclose(np.trace(actual.reshape(66, 2, 66, 2), axis1=0, axis2=2),
        np.trace(density.reshape(66, 2, 66, 2), axis1=0, axis2=2), rtol=0, atol=ORACLE_ATOL)


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('scale', [0., .125])
def test_zero_and_subnormalized_branches_are_never_normalized(law, scale):
    full = _full_input()
    prepared = build.prepare_full(full, law, 1)
    small = build.prepare_full(scale*full, law, 1)
    np.testing.assert_allclose(small, scale*prepared, rtol=0, atol=1e-16)
    first = build.trial_outputs(prepared, law, rounds=1)
    scaled = build.trial_outputs(small, law, rounds=1)
    for record in first:
        np.testing.assert_allclose(scaled[record], scale*first[record], rtol=0, atol=1e-16)
    assert sum(np.trace(value).real for value in scaled.values()) == pytest.approx(scale, abs=1e-15)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_two_trial_execution_retains_all_144_records_and_actual_subnormalized_inputs(law, monkeypatch):
    # One sweep deliberately does not establish the production accuracy.
    # Relax that assertion only to exercise complete composition cheaply.
    monkeypatch.setattr(build, 'PREPARATION_ALLOWANCE', F(1))
    original = build.prepare_full
    inputs, outputs = [], []

    def traced(density, candidate, count):
        inputs.append(density.copy())
        result = original(density, candidate, count)
        outputs.append(result.copy())
        return result

    monkeypatch.setattr(build, 'prepare_full', traced)
    report = build.execute_history(law, 1)
    records = build.public_records()
    assert len(report['complete_two_trial_law']) == 144
    assert {tuple(row['record']) for row in report['complete_two_trial_law']} == {
        (*first, *second) for first, second in product(records, repeat=2)}
    assert len(inputs) == 13
    np.testing.assert_array_equal(inputs[0], _full_input())
    first = _oracle(outputs[0], law, 3, 7)
    expected = {}
    for i, record in enumerate(records, start=1):
        np.testing.assert_allclose(inputs[i], build._embed_binary(first[record]), rtol=0, atol=ORACLE_ATOL)
        assert np.trace(outputs[i]).real == pytest.approx(np.trace(inputs[i]).real, abs=1e-15)
        for second, value in _oracle(outputs[i], law, 3, 7).items():
            expected[*record, *second] = float(np.trace(value).real)
    for row in report['complete_two_trial_law']:
        assert row['probability'] == pytest.approx(expected[tuple(row['record'])], abs=1e-15)
    assert sum(row['probability'] for row in report['complete_two_trial_law']) == pytest.approx(1., abs=2e-14)
    assert report['work_replaced_between_trials'] is False


def _counts(record=(0, 0, 0)):
    result = {r: 0 for r in build.public_records()}
    result[record] = build.TRIALS
    return result


@pytest.mark.parametrize('record, expected', [
    ((0, 0, 0), 'C'), ((0, 0, 1), 'outside_both_declared_models'),
    ((1, 4, 0), 'C'), ((1, 4, 1), 'E'),
    ((1, 1, 0), 'scheduler_timeout_inconclusive'),
    ((2, 3, 0), 'scheduler_timeout_inconclusive'),
])
def test_classifier_keeps_failed_heralds_impossible_records_and_timeouts(record, expected):
    assert build.classify(_counts(record)) == expected


def test_one_timeout_cannot_be_dropped_or_replaced_by_accepted_trial_total():
    counts = _counts((1, 4, 0))
    counts[1, 4, 0] -= 1
    counts[2, 2, 0] = 1
    assert build.classify(counts) == 'scheduler_timeout_inconclusive'
    counts[2, 2, 0] = 0
    with pytest.raises(ValueError):
        build.classify(counts)
    counts[1, 4, 0] += 1
    counts.pop((2, 2, 0))
    with pytest.raises(ValueError):
        build.classify(counts)


@pytest.mark.parametrize('key', [(False, 0, 0), (0., 0, 0), (np.int64(0), 0, 0),
                                (F(0), 0, 0), (0, 0), (3, 0, 0)])
def test_classifier_rejects_malformed_or_numeric_alias_record_keys(key):
    counts = _counts()
    value = counts.pop((0, 0, 0))
    counts[key] = value
    with pytest.raises(ValueError):
        build.classify(counts)


@pytest.mark.parametrize('value', [True, 0., np.int64(0), F(0), -1, None, '0'])
def test_classifier_rejects_nonoriginal_or_negative_counts(value):
    counts = _counts()
    counts[0, 0, 1] = value
    with pytest.raises(ValueError):
        build.classify(counts)


@pytest.mark.parametrize('attempts', [True, False, 1., F(1), np.int64(1), 0, 9, None])
def test_trial_rejects_invalid_attempt_cap(attempts):
    with pytest.raises(ValueError):
        build.trial_outputs(np.eye(66)/66, 'C', attempts=attempts)


@pytest.mark.parametrize('rounds', [True, False, 1., F(1), np.int64(1), 0, 33, None])
def test_trial_rejects_invalid_selector_cap(rounds):
    with pytest.raises(ValueError):
        build.trial_outputs(np.eye(66)/66, 'C', rounds=rounds)


@pytest.mark.parametrize('method', ['trial_outputs', 'uniform_native_step', 'prepare_full'])
@pytest.mark.parametrize('law', [True, 'C0', 'all', np.str_('C'), None])
def test_execution_rejects_invalid_candidate_before_any_shortcut(method, law):
    dimension = 78 if method == 'prepare_full' else 66
    arguments = (np.zeros((dimension, dimension), complex), law)
    with pytest.raises(ValueError):
        getattr(build, method)(*arguments, *(() if method != 'prepare_full' else (1,)))


@pytest.mark.parametrize('method', ['trial_outputs', 'uniform_native_step', 'prepare_full'])
@pytest.mark.parametrize('kind', ['list', 'integer', 'boolean', 'masked', 'nan', 'infinite',
                                 'vector', 'nonsquare', 'wrong_size', 'four_references'])
def test_execution_rejects_malformed_density_without_implicit_casting(method, kind):
    d = 78 if method == 'prepare_full' else 66
    density = np.eye(d)/d
    density = {'list': density.tolist(), 'integer': np.eye(d, dtype=int),
               'boolean': np.eye(d, dtype=bool), 'masked': np.ma.array(density, mask=False),
               'nan': np.full((d, d), np.nan), 'infinite': np.full((d, d), np.inf),
               'vector': np.zeros(d), 'nonsquare': np.zeros((d, d+1)),
               'wrong_size': np.zeros((d+1, d+1)),
               'four_references': np.eye(4*d)/(4*d)}[kind]
    with pytest.raises(ValueError):
        getattr(build, method)(density, 'C', *(() if method != 'prepare_full' else (1,)))


@pytest.mark.parametrize('count', [True, False, 1., F(1), np.int64(1), -1, 65, None])
def test_preparation_rejects_invalid_sweep_count(count):
    with pytest.raises(ValueError):
        build.prepare_full(np.eye(78)/78, 'C', count)


@pytest.mark.parametrize('method', ['trial_outputs', 'prepare_full'])
@pytest.mark.parametrize('kind', ['nonhermitian', 'negative', 'supernormalized'])
def test_density_entrypoints_reject_nonphysical_inputs(method, kind):
    d = 78 if method == 'prepare_full' else 66
    density = np.eye(d)/d
    if kind == 'nonhermitian':
        density[0, 1] = .1
    elif kind == 'negative':
        density[0, 0], density[1, 1] = -.1, .1+2/d
    else:
        density *= 2
    with pytest.raises(ValueError):
        getattr(build, method)(density, 'C', *(() if method != 'prepare_full' else (1,)))


def test_finite_pointer_mass_loss_is_rejected(monkeypatch):
    original = build.runtime._apply

    def apply(instrument, density, reference):
        outputs = original(instrument, density, reference)
        if len(outputs) == 2:
            outputs[0][0, 0] += .01
        return outputs

    monkeypatch.setattr(build.runtime, '_apply', apply)
    with pytest.raises(ArithmeticError, match='probability mass'):
        build.trial_outputs(np.eye(66)/66, 'C', attempts=1, rounds=1)


@pytest.mark.parametrize('corruption', [np.nan, np.inf, -np.inf])
def test_nonfinite_terminal_pointer_cannot_evade_complete_mass_check(monkeypatch, corruption):
    original = build.runtime._apply
    corrupted = []

    def apply(instrument, density, reference):
        outputs = original(instrument, density, reference)
        if len(outputs) == 2:  # Only the final binary pointer, after every work step.
            outputs[0][0, 0] = corruption
            corrupted.append(True)
        return outputs

    monkeypatch.setattr(build.runtime, '_apply', apply)
    with pytest.raises(ArithmeticError, match='nonfinite'):
        build.trial_outputs(np.eye(66)/66, 'C', attempts=1, rounds=1)
    assert corrupted


@pytest.mark.parametrize('side', [0, 1])
@pytest.mark.parametrize('corruption', [np.nan, np.inf, -np.inf])
def test_distance_rejects_nonfinite_complete_output_before_norm(side, corruption, monkeypatch):
    values = [np.eye(2)/2, np.eye(2)/2]
    values[side][0, 0] = corruption

    def unreachable(*args, **kwargs):
        raise AssertionError('nonfinite inputs must fail before the norm backend')

    monkeypatch.setattr(np.linalg, 'norm', unreachable)
    with pytest.raises(ArithmeticError, match='nonfinite'):
        build._distance(*values)


@pytest.mark.parametrize('corruption', [np.nan, np.inf])
def test_distance_rejects_nonfinite_norm_result(monkeypatch, corruption):
    monkeypatch.setattr(np.linalg, 'norm', lambda *args, **kwargs: corruption)
    with pytest.raises(ArithmeticError, match='nonfinite'):
        build._distance(np.eye(2)/2, np.eye(2)/2)


@pytest.fixture
def small_receipt():
    outputs = build.trial_outputs(np.eye(66)/66, 'C', rounds=1)
    return {'schema': 'oph.repair_scheduler.test.v1',
            'source_sha256': {'source.py': 'a'*64},
            'selector': build.coins.selector_certificate(rounds=1),
            'records': [{'record': list(key), 'probability': float(np.trace(value).real)}
                        for key, value in outputs.items()],
            'all_outputs_retained': True}


@pytest.mark.parametrize('mutation', [
    lambda r: r['records'].pop(),
    lambda r: r['records'][0]['record'].__setitem__(0, True),
    lambda r: r['records'][0].update(probability=float('nan')),
    lambda r: r['records'][0].update(probability=.5),
    lambda r: r['selector'].update(timeout_probability='0'),
    lambda r: r['selector'].update(maximum_auxiliary_native_attempts=True),
    lambda r: r['source_sha256'].clear(),
    lambda r: r.update(all_outputs_retained=False),
    lambda r: r.update(unproved_claim=True),
])
def test_complete_receipt_mutations_fail_shared_comparator(small_receipt, mutation):
    corrupted = copy.deepcopy(small_receipt)
    mutation(corrupted)
    with pytest.raises(ValueError):
        compare_receipt(corrupted, small_receipt)


def test_cli_checks_complete_receipt_and_rejects_duplicate_json(small_receipt, monkeypatch, tmp_path):
    from . import __main__ as command
    monkeypatch.setattr(command, 'run', lambda: copy.deepcopy(small_receipt))
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(small_receipt))
    monkeypatch.setattr(sys, 'argv', ['repair_scheduler', '--check', str(path)])
    command.main()
    path.write_text('{"schema":"first","schema":"second"}')
    with pytest.raises(ValueError):
        command.main()


@pytest.fixture
def retained_receipt():
    # Exercise the real compact schema without repeating its expensive producer.
    return json.loads(Path(__file__).with_name('results.json').read_text(encoding='utf-8'))


def _corrupt_receipt_domain(receipt, mutation):
    history = receipt['continuing_histories'][0]
    joint = history['joint_probabilities']
    auxiliary = receipt['auxiliary_executions'][0]
    if mutation == 'negative_impossible':
        joint[1][1] = -1e-10
    elif mutation == 'negative_possible':
        joint[2][2] = -1e-12
    elif mutation == 'positive_impossible':
        joint[1][1] = 1e-12
    elif mutation == 'negative_error':
        receipt['reference_histories'][0]['untouched_reference_marginal_error_diagnostic'] = -1e-12
    elif mutation == 'joint_mass':
        for i, row in enumerate(joint):
            for j in range(len(row)):
                if i != 1 and j != 1:
                    row[j] += 1e-10
    elif mutation == 'auxiliary_mass':
        for row in auxiliary['probabilities_by_round_then_address']:
            for i in range(len(row)):
                row[i] += 1e-10
    else:
        raise AssertionError(mutation)


@pytest.mark.parametrize('mutation', ['negative_impossible', 'negative_possible',
                                     'positive_impossible', 'negative_error',
                                     'joint_mass', 'auxiliary_mass'])
def test_cli_rejects_invalid_probability_laws_even_within_replay_tolerance(
        retained_receipt, mutation, monkeypatch, tmp_path):
    from . import __main__ as command
    corrupted = copy.deepcopy(retained_receipt)
    _corrupt_receipt_domain(corrupted, mutation)
    # This is precisely the gap that domain validation closes locally.
    compare_receipt(corrupted, retained_receipt)
    monkeypatch.setattr(command, 'run', lambda: copy.deepcopy(retained_receipt))
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(corrupted), encoding='utf-8')
    monkeypatch.setattr(sys, 'argv', ['repair_scheduler', '--check', str(path)])
    with pytest.raises(ValueError):
        command.main()


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf'), -1e-300, 1.0000000001, True])
def test_probability_domain_rejects_nonfinite_out_of_range_and_alias_values(retained_receipt, value):
    from . import __main__ as command
    retained_receipt['continuing_histories'][0]['joint_probabilities'][0][0] = value
    with pytest.raises(ValueError):
        command._validate_receipt_domains(retained_receipt)


@pytest.mark.parametrize('mutation', [
    lambda r: r.update(continuing_histories=None),
    lambda r: r['continuing_histories'][0].pop('joint_probabilities'),
    lambda r: r['continuing_histories'][0].update(joint_probabilities=[]),
    lambda r: r['continuing_histories'][0]['joint_probabilities'].pop(),
    lambda r: r['continuing_histories'][0]['joint_probabilities'][0].pop(),
    lambda r: r['continuing_histories'][0].pop('public_record_order'),
    lambda r: r['auxiliary_executions'][0].pop('probabilities_by_round_then_address'),
    lambda r: r['auxiliary_executions'][0].pop('timeout_probability_diagnostic'),
])
def test_domain_validation_rejects_malformed_grids_with_value_error(retained_receipt, mutation):
    from . import __main__ as command
    mutation(retained_receipt)
    with pytest.raises(ValueError):
        command._validate_receipt_domains(retained_receipt)


def test_legitimate_floating_rounding_remains_valid(retained_receipt, monkeypatch, tmp_path):
    from . import __main__ as command
    rounded = copy.deepcopy(retained_receipt)
    joint = rounded['continuing_histories'][0]['joint_probabilities']
    joint[0][0] = float(np.nextafter(joint[0][0], np.inf))
    joint[0][-1] = float(np.nextafter(joint[0][-1], 0.))
    auxiliary = rounded['auxiliary_executions'][0]['probabilities_by_round_then_address']
    auxiliary[0][0] = float(np.nextafter(auxiliary[0][0], np.inf))
    rounded['reference_histories'][0]['untouched_reference_marginal_error_diagnostic'] = 1e-15
    compare_receipt(rounded, retained_receipt)
    command._validate_receipt_domains(rounded)
    monkeypatch.setattr(command, 'run', lambda: copy.deepcopy(retained_receipt))
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(rounded), encoding='utf-8')
    monkeypatch.setattr(sys, 'argv', ['repair_scheduler', '--check', str(path)])
    command.main()


def test_write_cannot_emit_fresh_invalid_probability_diagnostics(retained_receipt, monkeypatch, tmp_path):
    from . import __main__ as command
    _corrupt_receipt_domain(retained_receipt, 'negative_impossible')
    monkeypatch.setattr(command, 'run', lambda: retained_receipt)
    path = tmp_path/'receipt.json'
    monkeypatch.setattr(sys, 'argv', ['repair_scheduler', '--write', str(path)])
    with pytest.raises(ValueError):
        command.main()
    assert not path.exists()

"""Success-conditioned accuracy, phase freedom and fail-closed census clocks."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import census_clock as c

WEIGHTS = (F(2, 225), F(1, 36), F(2, 45))


@pytest.fixture(scope='module')
def source_output():
    return c.run()


@pytest.fixture(scope='module')
def generated_receipt(source_output):
    return c.strict_json.loads(c.full.preparation.canonical(source_output))


@pytest.fixture
def receipt(generated_receipt, monkeypatch):
    # Generate once from the real producer, then isolate receipt mutations
    # from its expensive source compilation. The CLI always replays itself.
    monkeypatch.setattr(c, 'run', lambda: deepcopy(generated_receipt))
    return generated_receipt


def _counts():
    # An exact integer, fixed-budget record fixture near the native law.
    values = {key: 0 for key in c.census.LABELS}
    for q, (multiplicity, weight) in enumerate(zip((45, 20, 1), WEIGHTS)):
        values[f'success_q{q}'] = int(c.TRAINING_TRIALS*F(5, 33)*multiplicity*weight)
    values['herald_failure'] = c.TRAINING_TRIALS-sum(values.values())
    return values


def test_oscillation_uses_distinct_coordinates_and_covers_every_corner():
    box = ((WEIGHTS[0]*F(999, 1000), WEIGHTS[0]*F(1001, 1000)),
           (WEIGHTS[1]*F(999, 1000), WEIGHTS[1]*F(1001, 1000)),
           (WEIGHTS[2]*F(9, 10), WEIGHTS[2]*F(11, 10)))
    for s in (F(-2), F(), F(1), F(3, 2)):
        bound = c.inference_bound(WEIGHTS, box, s)
        for endpoints in product((0, 1), repeat=3):
            truth = tuple(box[q][endpoints[q]] for q in range(3))
            delta = float(s)*np.log(np.array(list(map(float, WEIGHTS)))/np.array(list(map(float, truth))))
            delta -= (delta.max()+delta.min())/2
            assert np.max(np.abs(np.exp(1j*delta)-1)) <= float(bound)+1e-15
    assert c.inference_bound(WEIGHTS, box) < F(3, 5)*c.learned.inference_bound(WEIGHTS, box)
    assert c.inference_bound(WEIGHTS, tuple((w, w) for w in WEIGHTS)) == 0
    scale = F(7, 13)
    assert c.inference_bound(WEIGHTS, box) == c.inference_bound(
        tuple(scale*w for w in WEIGHTS), tuple(tuple(scale*x for x in interval) for interval in box))


@pytest.mark.parametrize('weights,box', [
    ((F(), F(1), F(1)), ((F(1), F(1)),)*3),
    ((True, F(1), F(1)), ((F(1), F(1)),)*3),
    (WEIGHTS, ((F(), F(1)),)*3),
    (WEIGHTS, ((F(1), F(2)),)*3),
    (WEIGHTS, ((0., 1.),)*3),
    (WEIGHTS, [(F(1), F(2))]*3),
])
def test_invalid_weight_or_box_does_not_create_a_clock(weights, box):
    with pytest.raises(ValueError):
        c.inference_bound(weights, box)


@pytest.mark.parametrize('slot,value', [('_denominator', -1), ('_denominator', 0),
    ('_denominator', True), ('_denominator', 1.0), ('_numerator', True), ('_numerator', F(2))])
@pytest.mark.parametrize('position', ['weight', 'lower', 'upper', 'parameter'])
def test_malformed_fraction_slots_cannot_certify_a_clock(slot, value, position):
    weights = (F(2), F(1), F(1))
    intervals = ((F(2), F(2)), (F(1), F(1)), (F(1), F(1)))
    bad = F(2)
    setattr(bad, slot, value)
    parameter = F(1)
    if position == 'weight':
        weights = (bad, *weights[1:])
    elif position == 'parameter':
        parameter = bad
    else:
        intervals = (((bad, F(2)) if position == 'lower' else (F(2), bad)), *intervals[1:])
    with pytest.raises(ValueError):
        c.inference_bound(weights, intervals, parameter)


def test_exact_positive_fraction_inputs_are_canonicalized_before_comparison():
    # A valid positive value with unreduced slots has the same source value;
    # copying its integer pair canonicalizes it before equality/arithmetic.
    unreduced = F(2)
    unreduced._denominator = 2
    weights = (unreduced, F(1), F(1))
    boxes = tuple((value, value) for value in weights)
    assert c.inference_bound(weights, boxes, F(-1)) == 0
    assert c.inference_bound(weights, boxes, F(0)) == 0


def test_successful_admission_has_no_target_lookup_or_herald_dilution(monkeypatch):
    c.prior._record_contract.cache_clear()
    def forbidden(*args, **kwargs):
        raise AssertionError('target weights must not enter clock admission')
    monkeypatch.setattr(c.local, 'reference', forbidden)
    monkeypatch.setattr(c.local, 'budget', forbidden)
    decision = c.admission(_counts())
    assert decision['status'] == 'accept'
    budget = decision['confidence_event_bounds']
    assert F(1, 20) < budget['successful_two_read_channel_halfdiamond_upper'] < F(1, 10)
    assert budget['successful_reference_history_halftrace_upper'] <= c.SUCCESSFUL_ERROR_CAP
    assert budget['complete_flagged_history_halftrace_upper'] < budget['successful_reference_history_halftrace_upper']/5
    assert decision['complete_joint_protocol_halftrace_upper'] == (
        decision['fixed_complete_confidence_event_error_cap']
        +decision['inference']['confidence_failure_upper']*(c.learned.source_bounds()[1]+1)/2)
    assert decision['fixed_complete_confidence_event_error_cap'] < F(1, 66)
    assert decision['complete_joint_protocol_halftrace_upper'] < F(22, 1000)


def test_cap_is_fixed_even_when_observed_bound_changes_or_budget_is_wrong():
    baseline = c.admission(_counts())
    counts = _counts()
    counts['herald_failure'] += counts['success_q2']-100
    counts['success_q2'] = 100
    wide = c.admission(counts)
    assert wide['inference']['status'] == 'faithful'
    assert wide['status'] == 'stop'
    assert wide['reason'] == 'successful_history_confidence_box_too_wide'
    assert wide['complete_joint_protocol_halftrace_upper'] == baseline['complete_joint_protocol_halftrace_upper']
    assert c.admission({key: 2*value for key, value in _counts().items()})['reason'] == 'fixed_trial_budget_required'


@pytest.mark.parametrize('failure', ['unresolved', 'wide', 'compilation'])
def test_rejected_records_preserve_the_terminal_training_carrier(monkeypatch, failure):
    counts = _counts()
    if failure == 'unresolved':
        counts = {key: c.TRAINING_TRIALS if key == 'herald_failure' else 0 for key in counts}
    if failure == 'wide':
        counts['herald_failure'] += counts['success_q2']-100
        counts['success_q2'] = 100
    def certificate(*args, **kwargs):
        if failure != 'compilation':
            raise AssertionError('rejected training must not compile')
        raise ArithmeticError('source array missed its fixed gate')
    def prepared(*args, **kwargs):
        raise AssertionError('rejected training must not start future preparation')
    monkeypatch.setattr(c.learned, 'certificate', certificate)
    monkeypatch.setattr(c.learned, 'prepared', prepared)
    decision, outputs = c.prepared_from_records(counts, 'C', initial=37)
    assert decision['status'] == 'stop'
    assert len(outputs) == 1
    np.testing.assert_array_equal(next(iter(outputs.values())), np.diag([int(i == 37) for i in range(66)]))


def test_successful_channel_bound_on_entangled_input_retains_feedback_and_aborts():
    inferred = c.census.infer_counts(_counts())
    weights, intervals = inferred['central_subnormalized_weights'], inferred['subnormalized_weight_intervals']
    rng = np.random.default_rng(4128)
    vector = rng.normal(size=132)+1j*rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    raw = c.learned.execute(density, weights, 'C')
    ideal = c.learned.execute(density, WEIGHTS, 'C', raw=False)
    assert c.history.distance(raw, ideal) <= float(c.bounds(weights, intervals)['successful_two_read_channel_halfdiamond_upper'])
    assert all(a.shape == ((2, 2) if key[0] == 'abort' else (132, 132)) for key, a in raw.items())
    assert abs(sum(c.history.probabilities(raw).values())-1) < 1e-9
    assert c.history.distance(ideal, c.learned.execute(density, WEIGHTS, 'C', raw=False, feedback=False)) > 1e-4


@pytest.mark.parametrize('mutation', ['boolean', 'count', 'cap', 'source', 'mass', 'negative', 'success_error'])
def test_retained_receipt_rejects_proof_and_complete_history_mutations(mutation, receipt):
    actual = receipt
    changed = deepcopy(actual)
    if mutation == 'boolean':
        changed['training_counts']['warmup_timeout'] = False
    elif mutation == 'count':
        changed['training_counts']['success_q2'] += 1
    elif mutation == 'cap':
        changed['acceptance_protocol']['fixed_successful_error_cap'] = '1/9'
    elif mutation == 'source':
        key = next(iter(changed['source_sha256']))
        changed['source_sha256'][key] = '0'*64
    elif mutation == 'mass':
        changed['histories'][0]['record_probabilities'] = {
            key: value+1.9e-10 for key, value in changed['histories'][0]['record_probabilities'].items()}
    elif mutation == 'negative':
        key = next(iter(changed['histories'][0]['record_probabilities']))
        changed['histories'][0]['record_probabilities'][key] = -1e-12
    else:
        changed['histories'][0]['successful_output_halftrace_diagnostic'] = .2
    with pytest.raises(ValueError):
        c.compare_receipt(changed, actual)


def test_retained_receipt_passes_its_independent_consistency_checks(receipt):
    c.compare_receipt(receipt, receipt)


@pytest.mark.parametrize('mutation', ['bound', 'count', 'confidence', 'resource', 'conditioning', 'source_bound', 'source_path'])
def test_verifier_reconstructs_proof_fields_even_against_the_same_forged_receipt(mutation, receipt):
    changed = deepcopy(receipt)
    if mutation == 'bound':
        changed['acceptance_protocol']['fixed_complete_confidence_event_error_cap'] = '0'
    elif mutation == 'count':
        changed['training_counts']['success_q2'] += 1
        changed['training_counts']['herald_failure'] -= 1
    elif mutation == 'confidence':
        changed['statistical_failure_upper'] = '0'
    elif mutation == 'resource':
        changed['training_native_attempts_upper'] -= 1
    elif mutation == 'source_bound':
        changed['source_certificate']['source_history_halfdiamond_upper'] = '0'
    elif mutation == 'source_path':
        changed['source_sha256'].pop(next(iter(changed['source_sha256'])))
    else:
        changed['acceptance_protocol']['bound_conditioned_on_acceptance'] = True
    with pytest.raises(ValueError):
        c.compare_receipt(changed, changed)


@pytest.mark.parametrize('mutation', ['empty', 'duplicate_law', 'abort', 'extra_branch',
                                    'parameter', 'schema', 'scope', 'initial', 'stream_hash'])
def test_self_comparison_requires_complete_histories_and_declared_scope(mutation, receipt):
    changed = deepcopy(receipt)
    if mutation == 'empty':
        changed['histories'] = []
    elif mutation == 'duplicate_law':
        changed['histories'][1]['law'] = 'C'
    elif mutation == 'abort':
        changed['histories'][0]['record_probabilities'].pop("('abort', 0)")
    elif mutation == 'extra_branch':
        changed['histories'][0]['record_probabilities']['invented'] = 0.
    elif mutation == 'parameter':
        changed['parameter'] = '2'
    elif mutation == 'schema':
        changed['schema'] = 'oph.native_modular_reference.census_clock.v2'
    elif mutation == 'scope':
        changed['scope']['native_event_count_identified_with_modular_parameter_or_physical_seconds'] = True
    elif mutation == 'initial':
        changed['continuation_initial_native_state'] = True
    else:
        changed['training_record_and_state_sha256'] = 'g'*64
    with pytest.raises(ValueError):
        c.compare_receipt(changed, changed)


@pytest.mark.parametrize('mutation', ['stream_hash', 'normalized_records'])
def test_matching_fabricated_history_requires_source_replay(receipt, mutation):
    changed = deepcopy(receipt)
    if mutation == 'stream_hash':
        changed['training_record_and_state_sha256'] = '0'*64
    else:
        probabilities = changed['histories'][0]['record_probabilities']
        # Preserve success/preparation masses, the complete alphabet, and
        # positivity while fabricating which successful records occurred.
        keys = [key for key, mass in probabilities.items() if key.startswith("('read',") and mass > 1e-5]
        first, second = keys[:2]
        moved = probabilities[first]/2
        probabilities[first] -= moved
        probabilities[second] += moved
    c._compare_receipt(changed, changed)
    with pytest.raises(ValueError):
        c.compare_receipt(changed, changed)
    with pytest.raises(ValueError):
        c.compare_receipt(receipt, changed)
    c.compare_receipt(receipt)


def test_public_verifier_owns_replay_for_each_comparison(generated_receipt, monkeypatch):
    calls = []
    def replay():
        calls.append('source replay')
        return deepcopy(generated_receipt)
    monkeypatch.setattr(c, 'run', replay)
    c.compare_receipt(generated_receipt)
    c.compare_receipt(generated_receipt, deepcopy(generated_receipt))
    assert calls == ['source replay', 'source replay']


def test_public_bounds_do_not_alias_fixed_caps_or_cached_source_certificate(source_output, generated_receipt):
    decision = c.admission(_counts())
    source = source_output['source_certificate']
    weights = source_output['inferred_subnormalized_weights']
    expected_source = c.full.preparation.canonical(c.learned.certificate(weights))
    expected_bounds = c.full.preparation.canonical(decision)
    exposed = [decision['fixed_successful_error_cap'],
               source_output['acceptance_protocol']['fixed_successful_error_cap'],
               *(value for value in source.values() if type(value) is F)]
    original = [(value, value.numerator, value.denominator) for value in exposed]
    try:
        for value in exposed:
            value._numerator = 0
        assert c.SUCCESSFUL_ERROR_CAP == F(1, 10)
        assert c.full.preparation.canonical(c.admission(_counts())) == expected_bounds
        assert c.full.preparation.canonical(c.learned.certificate(weights)) == expected_source
    finally:
        for value, numerator, denominator in original:
            value._numerator, value._denominator = numerator, denominator

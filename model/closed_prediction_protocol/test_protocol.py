"""Independent source, dependence, checkpoint and locked-design controls."""
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest

from ..native_interface import construction as interface
from ..native_quantum_history import construction as quantum
from ..observer_prediction_audit import prospective
from ..primitive_observer import runtime
from . import construction as c


@pytest.fixture(scope='module')
def result():
    return c.run()


def test_exact_preparation_bound_and_minimal_certified_attempt_count():
    for error in (F(1, 1000), F(1, 1570800), F(1, 10**9)):
        count = c.preparation_count(error)
        assert F(33)*F(229, 240)**count <= error
        assert F(33)*F(229, 240)**(count-1) > error
    assert c.preparation_count(F(1, 1570800)) == 379
    for error in (1e-5, F(0), F(1), F(-1), True):
        with pytest.raises(ValueError):
            c.preparation_count(error)


def test_locked_design_includes_every_failure_and_accounts_for_all_native_costs():
    result = c.design()
    selected = result['selected']
    assert selected['unread_attempts'] == 3
    assert selected['fixed_trials'] == 3927
    assert selected['maximum_native_attempts'] == 1507968
    assert selected['maximum_native_attempts'] == min(row['maximum_native_attempts'] for row in result['locked_design_search'])
    for row in result['locked_design_search']:
        n = row['fixed_trials']
        certificate = row['likelihood_certificate']
        assert n*F(certificate['log_inverse_affinity_lower']) >= F(certificate['log_inverse_error_upper'])
        assert n*F(row['preparation_trace_distance_upper']) == F(1, 400)
        assert row['maximum_native_attempts'] == n*(row['discarded_preparation_attempts_per_trial']+row['unread_attempts']+2)
    budget = result['error_budget']
    assert F(budget['exact_source_model_decision_upper']) == F(3, 400)
    assert F(budget['additional_conditional_trial_TV_allowance'])*3927 == F(1, 400)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_sparse_execution_matches_actual_native_kraus_on_complex_reference(law):
    # A dense entangled input exercises all nonzero matrix entries, rather
    # than only the source's selected phase-pair subspace. This sample is
    # not an independent all-matrix-unit identity proof.
    rng = np.random.default_rng(217)
    vector = rng.normal(size=132)+1j*rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    instrument = quantum.native_instrument(law)
    expected = density.copy()
    for _ in range(3):
        expected, = runtime._apply(instrument, expected, 2)
    actual = c.prepare(density, law, 3)
    assert np.linalg.norm(actual-expected, 'nuc') < 2e-12
    assert abs(np.trace(actual)-1) < 2e-12
    assert np.linalg.eigvalsh(actual)[0] > -2e-12
    assert np.array_equal(c.prepare(density, law, 0), density)


def test_preparation_retains_unnormalized_branches_and_zero_records():
    rho = np.eye(66, dtype=complex)/66
    assert np.linalg.norm(c.prepare(rho/7, 'C', 8)-c.prepare(rho, 'C', 8)/7) < 1e-14
    assert not np.any(c.prepare(np.zeros((66, 66), complex), 'E', 379))
    for bad in (True, -1, 10001, .5):
        with pytest.raises(ValueError):
            c.prepare(rho, 'C', bad)
    for bad in (np.zeros((3, 3)), np.full((66, 66), np.nan)):
        with pytest.raises(ValueError):
            c.prepare(bad, 'C', 1)


def test_checkpoint_replays_complete_density_and_rejects_mutation():
    rho = c.prepare(interface.reference_input(), 'C', 3)
    payload = c._serialize(rho)
    assert np.array_equal(c._restore(payload, rho), rho)
    mutated = json.loads(payload)
    mutated['real'][0][0] += .01
    with pytest.raises(ValueError):
        c._restore(json.dumps(mutated), rho)
    mutated = json.loads(payload)
    mutated['shape'] = [66, 66]
    with pytest.raises(ValueError):
        c._restore(json.dumps(mutated), rho)


def test_without_preparation_a_previous_record_changes_the_next_trial():
    first = c.trial_outputs(np.eye(66)/66, 'C', 3)
    probabilities = prospective.exact_law('C', 3)
    # A nonzero failed herald already changes the next preparation yield.
    failed = first[0, 0]
    following = c.trial_outputs(failed, 'C', 3)
    mass = float(np.trace(failed).real)
    conditional = [float(np.trace(following[key]).real)/mass for key in c.RECORDS]
    assert sum(abs(p-float(q)) for p, q in zip(conditional, probabilities))/2 > .01


def test_conditional_perturbations_need_no_independence_for_transcript_bound():
    p, q = prospective.exact_law('C', 3), prospective.exact_law('E', 3)
    epsilon = F(1, 1000)
    for base, other, candidate in ((p, q, 'C'), (q, p, 'E')):
        # History-dependent errors are explicitly correlated: the second
        # perturbation depends on the complete first-trial outcome.
        first = list(base); first[0] -= epsilon; first[2] += epsilon
        joint = {}
        for a, b in product(range(6), repeat=2):
            conditional = list(base)
            shift = epsilon if a % 2 else -epsilon
            conditional[0] += shift; conditional[2] -= shift
            joint[a, b] = first[a]*conditional[b]
        assert sum(joint.values()) == 1
        tv = sum(abs(joint[a, b]-base[a]*base[b]) for a, b in joint)/2
        assert tv <= 2*epsilon
        error = F(0)
        for word, probability in joint.items():
            counts = tuple(word.count(index) for index in range(6))
            if prospective.classify(counts, p, q, 2) != candidate:
                error += probability
        affinity_upper = F(prospective.trial_design(base, other)['affinity_upper'])
        assert error <= affinity_upper**2+2*epsilon


def test_fixed_future_classifier_rejects_incomplete_or_impossible_records():
    n = c.design()['selected']['fixed_trials']
    for law in ('C', 'E'):
        probabilities = prospective.exact_law(law, 3)
        counts = [int(n*p) for p in probabilities]
        counts[0] += n-sum(counts)
        assert c.classify(tuple(counts)) == law
    assert c.classify((n-1, 1, 0, 0, 0, 0)) == 'outside_both_declared_models'
    with pytest.raises(ValueError):
        c.classify((n-1, 0, 0, 0, 0, 0))


def test_actual_two_trial_and_reference_execution_obey_the_uniform_bound(result):
    epsilon = F(result['design']['selected']['preparation_trace_distance_upper'])
    for row in result['sequential_execution']:
        assert len(row['complete_two_trial_law']) == 36
        assert row['initial_preparation_trace_distance_diagnostic'] <= float(epsilon)+2e-10
        assert row['two_trial_TV_from_ideal_product_diagnostic'] <= 2*float(epsilon)+2e-10
        assert row['maximum_second_conditional_law_TV_diagnostic'] <= float(epsilon)+2e-10
        assert row['normalization_residual_diagnostic'] < 2e-10
        assert not row['initial_state_replaced_between_trials']
    for row in result['entangled_reference_execution']:
        assert row['prepared_distance_from_product_diagnostic'] <= float(epsilon)+2e-10
        assert row['untouched_reference_marginal_distance_diagnostic'] < 2e-10
    assert result['resource_budget']['maximum_source_pulses'] == 267036
    assert result['fresh_coherent_source']['resources_per_mixer']['native_pulses'] == 34


def test_receipt_keeps_exact_error_budget_and_source_comparison_scope(result):
    assert json.loads(json.dumps(result)) == result
    altered = json.loads(json.dumps(result))
    altered['design']['error_budget']['exact_source_model_decision_upper'] = '0'
    with pytest.raises((ValueError, AssertionError)):
        c.compare(altered, result)
    assert not result['fresh_coherent_source']['normalized_physical_instrument_error_certified']
    assert not result['theorem_attachment']['iid_assumed_for_executed_process']

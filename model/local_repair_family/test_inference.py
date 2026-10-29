"""Independent probability, shared-reference and hostile decoder controls."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations
import json
import math

import numpy as np
import pytest

from . import inference as inference


def counts(weights=(25, 25)+(8,)*10):
    n = 1 << 20
    # Only this test/acquisition-side oracle knows the true port weights.
    return {(a, b): (n, (n*weights[a]*2//(weights[a]+weights[b])+1)//2)
            for a, b in inference.SEAMS}


def test_fixed_budget_uses_independent_high_precision_exponential_and_phase_constant():
    result = inference.fixed_budget()
    assert result['native_attempts'] == 31_457_280
    assert result['clock_halfdiamond_upper_on_confidence_event'] == F(3, 32)
    with localcontext() as context:
        context.prec = 90
        exponent = Decimal(2)*Decimal(1 << 20)*(Decimal(3)/1000)**2
        independent = 60*(-exponent).exp()
        actual = result['statistical_failure_upper']
        exact_decimal = Decimal(actual.numerator)/Decimal(actual.denominator)
        assert independent <= exact_decimal < Decimal(4)/10**7
        assert exact_decimal-independent < Decimal(1)/10**12
    assert result['complete_inference_halfdiamond_upper'] < F(15, 100)
    assert result['bound_conditioned_on_admission'] is False


def test_certificate_exact_menu_and_no_hidden_weight_input():
    report = inference.exact_certificate()
    assert json.loads(json.dumps(report)) == report
    # Enumerate by reduced rational denominator rather than ordered (a,b).
    menu = {F(a, denominator) for denominator in range(2, 65)
            for a in range(1, denominator)
            if a <= 32 and denominator-a <= 32
            and F(1, 5) <= F(a, denominator) <= F(4, 5)}
    assert len(menu) == 489
    assert min(abs(a-b) for a in menu for b in menu if a != b) == F(1, 3843)
    assert report['empirical_counts_treated_as_exact_probabilities'] is False
    with pytest.raises(TypeError):
        inference.infer_counts(counts(), true_weights=(25, 25)+(8,)*10)


@pytest.mark.parametrize('weights', [(25, 25)+(8,)*10, (1,)*12,
                                    (1, 3)*6, tuple(range(8, 20))])
def test_same_true_reference_in_every_interval_and_exact_clock_bound(weights):
    result = inference.infer_counts(counts(weights))
    assert result['status'] == 'accept'
    assert len(result['held_out_checks']) == 19
    assert all(row['compatible'] for row in result['held_out_checks'])
    estimated = result['relative_weights']
    for actual, interval in zip(weights, result['weight_intervals'], strict=True):
        truth = F(actual, weights[0])
        assert interval[0] <= truth <= interval[1]
    # Actual H66 eigenvalues, with ONE common phase removed. This is not a
    # collection of independently chosen two-port reference witnesses.
    errors = [math.log(float(estimated[a]*estimated[b]/F(weights[a]*weights[b], weights[0]**2)))
              for a, b in combinations(range(12), 2)]
    radius = (max(errors)-min(errors))/2
    assert radius <= float(result['budget']['clock_halfdiamond_upper_on_confidence_event'])
    centered = np.array(errors)-(max(errors)+min(errors))/2
    # Tensor with a complex reference and test a coherent superposition of
    # extremal energy differences, independent of the producer's bound code.
    i, j = int(np.argmin(centered)), int(np.argmax(centered))
    vector = np.zeros((66, 2), complex)
    vector[i, 0] = 1/math.sqrt(2); vector[j, 1] = 1j/math.sqrt(2)
    rotated = np.exp(1j*centered)[:, None]*vector
    distance = math.sqrt(max(0., 1-abs(np.vdot(vector, rotated))**2))
    assert distance <= float(result['budget']['clock_halfdiamond_upper_on_confidence_event'])
    witness = result['shared_feasibility_witness']
    for row in result['observations']:
        a, b = row['seam']; lo, hi = row['probability_interval']
        assert lo <= witness[a]/(witness[a]+witness[b]) <= hi


def test_all_coin_columns_are_state_independent_including_collisions_and_reference():
    # Independent floor/ceiling and record partition; no inference or local
    # channel helper supplies this completeness calculation.
    from ..observer import native
    states = native.states()
    for weights in ((25, 25)+(8,)*10, (1, 3)*6):
        for a, b in inference.SEAMS:
            for coin in (0, 1):
                probability = F(weights[a if coin else b], weights[a]+weights[b])
                for law in ('C', 'E'):
                    groups = {}
                    for i, state in enumerate(states):
                        total = state[a]+state[b]
                        output = list(state)
                        output[a] = (total+coin)//2
                        output[b] = (total+1-coin)//2
                        value = state[a]-output[a] if law == 'C' else state[a]-state[b]
                        groups.setdefault(value, []).append((i, tuple(output)))
                    effect = [F()]*78
                    for pairs in groups.values():
                        assert len({output for _, output in pairs}) == len(pairs)
                        for i, _ in pairs:
                            effect[i] += probability
                    assert effect == [probability]*78


def test_held_out_edges_validate_but_do_not_fit_the_estimate_or_error_budget():
    original = counts()
    base = inference.infer_counts(original)
    edge = base['design']['held_out_seams'][0]
    changed = dict(original)
    n, k = changed[edge]; changed[edge] = (n, k+100)
    replay = inference.infer_counts(changed)
    assert replay['status'] == 'accept'
    assert replay['relative_weights'] == base['relative_weights']
    assert replay['weight_intervals'] == base['weight_intervals']
    assert replay['budget'] == base['budget']
    changed[edge] = (n, n//5)
    stopped = inference.infer_counts(changed)
    assert stopped['status'] == 'stop'
    assert stopped['stop_reason'] == 'held_out_cycle_incompatible'
    assert stopped['relative_weights'] is None
    assert stopped['budget'] == base['budget']


@pytest.mark.parametrize('pattern', range(6))
def test_joint_confidence_event_cannot_be_rejected_even_at_adversarial_interval_edges(pattern):
    weights = ((1, 4)*6 if pattern < 3 else (25, 25)+(8,)*10)
    data = counts(weights)
    displacement = (3*(1 << 20))//1000-1
    for number, (edge, (n, k)) in enumerate(data.items()):
        sign = 1 if ((number+pattern) % (pattern+2)) else -1
        data[edge] = (n, k+sign*displacement)
        a, b = edge
        assert abs(F(data[edge][1], n)-F(weights[a], weights[a]+weights[b])) < F(3, 1000)
    result = inference.infer_counts(data)
    assert result['status'] == 'accept'
    estimates = result['relative_weights']
    with localcontext() as context:
        context.prec = 70
        errors = []
        for a, b in combinations(range(12), 2):
            ratio = estimates[a]*estimates[b]/F(weights[a]*weights[b], weights[0]**2)
            errors.append((Decimal(ratio.numerator)/Decimal(ratio.denominator)).ln())
        assert (max(errors)-min(errors))/2 < Decimal(3)/32


def test_separate_edge_intersections_do_not_replace_one_common_feasibility_witness():
    # On a triangle every interval contains the positive ratio 11/10. Each
    # coordinate is locally valid, but their oriented product cannot be one.
    intervals = {(0, 1): (F(11, 10), F(6, 5)),
                 (1, 2): (F(11, 10), F(6, 5)),
                 (0, 2): (F(4, 5), F(10, 11))}
    assert inference._feasible(intervals) is None
    valid = {(0, 1): (F(1), F(6, 5)), (1, 2): (F(1), F(6, 5)),
             (0, 2): (F(1), F(36, 25))}
    witness = inference._feasible(valid)
    assert witness is not None
    assert all(lo <= witness[a]/witness[b] <= hi for (a, b), (lo, hi) in valid.items())


def test_domain_incompatibility_is_a_retained_stop_not_a_missing_or_zero_weight_clock():
    data = counts()
    data[inference.SEAMS[0]] = (1 << 20, 0)
    result = inference.infer_counts(data)
    assert result['status'] == 'stop'
    assert result['stop_reason'] == 'probability_domain_incompatible'
    assert result['relative_weights'] is None
    assert result['budget']['bound_conditioned_on_admission'] is False
    assert len(result['observations']) == 30


@pytest.mark.parametrize('case', ['missing', 'extra', 'reversed', 'bool_key',
                                 'list_value', 'bool_count', 'float_count',
                                 'negative_count', 'over_count', 'wrong_budget', 'subclass'])
def test_original_complete_count_domain_rejected_cold_and_after_valid_execution(case):
    class EvilDict(dict):
        pass
    data = counts(); edge = inference.SEAMS[0]; n, k = data[edge]
    if case == 'missing': del data[edge]
    elif case == 'extra': data[(0, 2)] = (n, k)
    elif case == 'reversed': data[edge[::-1]] = data.pop(edge)
    elif case == 'bool_key': data[(False, 1)] = data.pop(edge)
    elif case == 'list_value': data[edge] = [n, k]
    elif case == 'bool_count': data[edge] = (n, True)
    elif case == 'float_count': data[edge] = (n, float(k))
    elif case == 'negative_count': data[edge] = (n, -1)
    elif case == 'over_count': data[edge] = (n, n+1)
    elif case == 'wrong_budget': data[edge] = (n+1, k)
    elif case == 'subclass': data = EvilDict(data)
    for _ in range(2):
        with pytest.raises(ValueError): inference.infer_counts(data)
        assert inference.infer_counts(counts())['status'] == 'accept'


@pytest.mark.parametrize('bad', [True, 1., float('nan'), '1', F(1001, 1000), -2])
def test_parameter_domain_cold_and_warm(bad):
    for _ in range(2):
        with pytest.raises(ValueError): inference.fixed_budget(bad)
        assert inference.fixed_budget(1)['parameter'] == 1


def test_malformed_fraction_and_public_result_mutations_cannot_poison_later_results():
    malformed = F(1)
    object.__setattr__(malformed, '_denominator', -1)
    with pytest.raises(ValueError): inference.fixed_budget(malformed)
    result = inference.infer_counts(counts())
    baseline = inference.infer_counts(counts())
    object.__setattr__(result['relative_weights'][1], '_numerator', 999)
    object.__setattr__(result['weight_intervals'][2][0], '_numerator', 999)
    object.__setattr__(result['budget']['probability_radius'], '_numerator', 999)
    result['status'] = 'forged'
    assert inference.infer_counts(counts()) == baseline


def test_poisoned_inherited_logarithm_cache_cannot_change_any_inference_output():
    from ..native_modular_reference import common_caps
    baseline = inference.infer_counts(counts())
    # This ratio was evaluated by the removed optional log metadata path.
    ratio = baseline['relative_weights'][1]/baseline['weight_intervals'][1][1]
    cached = common_caps.logarithm_interval(ratio)
    originals = {id(x): (x, x.numerator, x.denominator) for x in cached}
    try:
        for x, n, d in originals.values():
            object.__setattr__(x, '_numerator', n+100*d)
        assert common_caps.logarithm_interval(ratio)[0] > 90
        assert inference.infer_counts(counts()) == baseline
        assert 'relative_log_error_intervals' not in baseline
    finally:
        for x, n, d in originals.values():
            object.__setattr__(x, '_numerator', n)
            object.__setattr__(x, '_denominator', d)
    assert common_caps.logarithm_interval(ratio)[0] < 0


def test_exact_identification_is_unique_up_to_primitive_scale_and_checks_all_edges():
    weights = (25, 25)+(8,)*10
    probabilities = {(a, b): F(weights[a], weights[a]+weights[b]) for a, b in inference.SEAMS}
    assert inference.exact_identification(probabilities) == weights
    scaled = {(a, b): F(2*weights[a], 2*weights[a]+2*weights[b]) for a, b in inference.SEAMS}
    assert inference.exact_identification(scaled) == weights
    edge = next(e for e in inference.SEAMS if e not in inference.TREE)
    probabilities[edge] += F(1, 10000)
    with pytest.raises(ValueError): inference.exact_identification(probabilities)
    with pytest.raises(ValueError):
        inference.exact_identification({e: F(k, n) for e, (n, k) in counts().items()})


def test_native_graph_relabeling_changes_only_the_recovered_port_labels():
    # An independently found exact graph automorphism; the decoder continues
    # to use its fixed public tree, not a truth-selected acquisition design.
    # Reflection through antipodal pairing is an actual graph automorphism.
    adjacency = {i: {b if a == i else a for a, b in inference.SEAMS if i in (a, b)}
                 for i in range(12)}
    permutation = []
    for i in range(12):
        distance = {i: 0}; todo = [i]
        for a in todo:
            for b in adjacency[a]:
                if b not in distance: distance[b] = distance[a]+1; todo.append(b)
        permutation.append(next(j for j in range(12) if distance[j] == 3))
    assert {tuple(sorted((permutation[a], permutation[b]))) for a, b in inference.SEAMS} == set(inference.SEAMS)
    weights = (1, 3)*6
    moved = tuple(weights[permutation[i]] for i in range(12))
    actual = inference.infer_counts(counts(moved))
    assert actual['relative_weights'] == tuple(F(w, moved[0]) for w in moved)


def test_zero_and_negative_parameters_share_confidence_failure_without_fake_conditioning():
    zero = inference.fixed_budget(0)
    negative = inference.fixed_budget(F(-1, 2))
    assert zero['clock_halfdiamond_upper_on_confidence_event'] == 0
    assert negative['clock_halfdiamond_upper_on_confidence_event'] == F(3, 64)
    assert zero['complete_inference_halfdiamond_upper'] == zero['statistical_failure_upper']

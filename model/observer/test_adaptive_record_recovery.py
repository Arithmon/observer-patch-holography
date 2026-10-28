"""Independent endpoint and literal-history controls for adaptive tape decay."""
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import adaptive_record_recovery as adaptive
from . import native, record_recovery


def endpoint(state, edge, coin, law):
    """Independent integer endpoint balancing and complete raw encoding."""
    a, b = edge
    low, parity = divmod(state[a]+state[b], 2)
    out = list(state)
    out[a], out[b] = low+coin*parity, low+(1-coin)*parity
    current = state[a]-out[a]
    raw = {0: 0, 1: 1, -1: 2}[current]
    return tuple(out), raw+(4*parity if law == 'E' else 0)


def literal_pair_histories(controller, left, right, attempts):
    """No backward recurrence: retain each surviving complete record word."""
    histories = [((), controller.initial_memory, left, right, F(1))]
    for _ in range(attempts):
        following = []
        for history, memory, x, y, weight in histories:
            for action, (law, edge) in enumerate(adaptive.actions()):
                probability = controller.probabilities[memory][action]
                if not probability:
                    continue
                for coin in (0, 1):
                    out_x, raw_x = endpoint(x, edge, coin, law)
                    out_y, raw_y = endpoint(y, edge, coin, law)
                    if raw_x == raw_y:
                        next_memory = controller.updates[memory][action][coin][raw_x]
                        following.append((history+((action, coin, raw_x),), next_memory,
                                          out_x, out_y, weight*probability/2))
        histories = following
    return histories


def one_law_controller(law, eta, preferred=0):
    probabilities = tuple(eta/30+(1-eta if i % 30 == preferred else 0)
                          if name == law else F(0)
                          for i, (name, _) in enumerate(adaptive.actions()))
    updates = tuple(tuple(tuple(0 for _ in range(4 if name == 'C' else 8))
                          for _ in (0, 1)) for name, _ in adaptive.actions())
    return adaptive.Controller(eta, (probabilities,), (updates,))


def test_all_raw_maps_match_independent_endpoint_oracle_and_every_zero_outcome():
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    rows = adaptive.branch_maps()
    assert len(rows) == 720
    assert sum(len(row.inputs) for row in rows) == 9360
    assert any(not row.inputs for row in rows)
    for row in rows:
        law, edge = adaptive.actions()[row.action]
        pairs = [(i, index[out]) for i, state in enumerate(states)
                 for out, raw in (endpoint(state, edge, row.coin, law),) if raw == row.raw]
        assert tuple(zip(row.inputs, row.outputs, strict=True)) == tuple(pairs)
        assert len(row.outputs) == len(set(row.outputs))


def test_all_pair_seam_agreement_counts_and_exact_maxima():
    states = native.states()
    for law in ('C', 'E'):
        totals = np.zeros((78, 78), dtype=int)
        for edge in record_recovery.seams():
            actual = np.zeros((78, 78), dtype=int)
            for coin in (0, 1):
                labels = np.array([endpoint(state, edge, coin, law)[1] for state in states])
                actual += labels[:, None] == labels[None, :]
            a, b = edge
            contrasts = np.array([state[a]-state[b] for state in states])
            delta = contrasts[:, None]-contrasts[None, :]
            expected = 2-np.minimum(abs(delta), 2) if law == 'C' else 2*(delta == 0)
            np.testing.assert_array_equal(actual, expected)
            totals += actual
        np.testing.assert_array_equal(np.diag(totals), np.full(78, 60))
        assert max(totals[i, j] for i in range(78) for j in range(78) if i != j) == (
            50 if law == 'C' else 42)
    certificate = adaptive.exact_certificate()
    assert certificate['law_pair_seam_checks'] == 365040
    assert certificate['uniform_offdiagonal_agreement_maxima'] == {'C': '5/6', 'E': '7/10'}
    certificate['uniform_offdiagonal_agreement_maxima']['C'] = '1'
    assert adaptive.exact_certificate()['uniform_offdiagonal_agreement_maxima']['C'] == '5/6'


@pytest.mark.parametrize('binary', [False, True])
def test_switching_controller_dp_matches_complete_literal_histories(binary):
    controller = adaptive.example_controller()
    states = native.states(binary=binary)
    kernel = adaptive.finite_controller_kernel(controller, 2, binary=binary)
    assert all(kernel[i][i] == 1 for i in range(len(states)))
    assert kernel == tuple(zip(*kernel, strict=True))
    # Different collision, binary and neighboring inputs; both memory states
    # and both record laws occur. The oracle keeps every surviving full tape.
    pairs = ((0, 1), (0, len(states)-1), (len(states)-2, len(states)-1))
    law_words = set()
    for i, j in pairs:
        histories = literal_pair_histories(controller, states[i], states[j], 2)
        assert kernel[i][j] == sum((row[-1] for row in histories), F(0))
        law_words.update(tuple(adaptive.actions()[h[0]][0] for h in row[0]) for row in histories)
    assert ('C', 'C') in law_words and ('C', 'E') in law_words
    assert any(0 < value < 1 for i, row in enumerate(kernel) for j, value in enumerate(row) if i != j)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_uniform_controller_recovers_existing_iid_all_pair_kernel(law):
    controller = one_law_controller(law, F(1))
    actual = adaptive.finite_controller_kernel(controller, 2)
    expected = record_recovery.iid_rates(law)
    assert actual == tuple(tuple(rate**2 for rate in row) for row in expected)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_persistent_exploration_rate_is_attained_and_remains_nonzero(law):
    controller = one_law_controller(law, F(1, 2))
    states = native.states()
    kernel = adaptive.finite_controller_kernel(controller, 2)
    upper = adaptive.bounds(F(1, 2), 2, law_scope=law)['pair_agreement_upper']
    assert max(value for i, row in enumerate(kernel) for j, value in enumerate(row) if i != j) == upper
    # A common positive-weight tape has rank at least two: finite-time EB is
    # not a consequence of the quantitative bound, even with E records.
    i, j = next((i, j) for i in range(78) for j in range(i) if kernel[i][j] == upper)
    histories = literal_pair_histories(controller, states[i], states[j], 2)
    assert histories and all(x != y for _, _, x, y, _ in histories)
    assert all(weight > 0 for *_, weight in histories)


def test_complete_reference_outputs_have_the_retained_coherence_not_a_postselected_one():
    controller = adaptive.example_controller()
    states = native.states()
    i, j = 0, 77
    histories = literal_pair_histories(controller, states[i], states[j], 2)
    kernel = adaptive.finite_controller_kernel(controller, 2)
    # Bell work/reference input on the two chosen carrier states. Each common
    # raw flag has an unnormalized off-diagonal 2x2 block of weight/2.
    # Distinct flags are summed in trace norm, never coherently or normalized.
    distance = 0.
    for *_, weight in histories:
        block = np.array([[0, 1j*float(weight)/2], [-1j*float(weight)/2, 0]])
        distance += np.abs(np.linalg.eigvalsh(block)).sum()/2
    assert distance == pytest.approx(float(kernel[i][j]/2), rel=2e-13, abs=0)
    assert 0 < kernel[i][j] < 1
    # Initial populations and their untouched reference labels retain total
    # probability one, including histories that have destroyed the coherence.
    for k in (i, j):
        complete = literal_pair_histories(controller, states[k], states[k], 2)
        assert sum((row[-1] for row in complete), F(0)) == 1


def test_exact_lifetime_thresholds_and_extremely_small_surviving_coherence():
    for law, count in (('C', 58), ('E', 31)):
        bound = adaptive.bounds(F(1, 2), count, law_scope=law)
        previous = adaptive.bounds(F(1, 2), count-1, law_scope=law)
        assert bound['complete_reference_trace_distance_upper'] <= F(1, 4)
        assert previous['complete_reference_trace_distance_upper'] > F(1, 4)
        assert bound['encoded_final_recovery_trace_error_lower'] >= F(1, 4)
    tiny = adaptive.bounds(F(1, 10**300), 2)
    assert 0 < tiny['one_step_agreement_upper'] < 1
    assert tiny['pair_agreement_upper'] < tiny['one_step_agreement_upper']
    late = adaptive.bounds(F(1, 2), 10000)
    assert 0 < late['pair_agreement_upper'] < F(1, 10**300)
    assert late['full_input_optimal_entanglement_fidelity_upper'] > F(1, 78)
    assert late['encoded_final_recovery_trace_error_lower'] < F(1, 2)


def test_zero_event_identity_and_report_is_exact_json_safe():
    kernel = adaptive.finite_controller_kernel(adaptive.example_controller(), 0)
    assert all(value == 1 for row in kernel for value in row)
    packet = adaptive.report()
    json.dumps(packet, allow_nan=False)
    assert packet['executed_exact_controller']['attempts'] == 3
    assert not packet['scope']['finite_horizon_entanglement_breaking_claimed']
    assert not packet['scope']['whole_program_or_refinement_no_go']
    assert set(adaptive.source_bindings()) >= {
        'model/observer/adaptive_record_recovery.py',
        'model/observer/test_adaptive_record_recovery.py',
        'model/observer/native.py', 'model/observer/record_recovery.py'}


@pytest.mark.parametrize('eta', [True, 1, .5, float('nan'), F(0), F(-1, 2), F(3, 2)])
def test_rejects_noncanonical_or_nonexploring_original_parameters(eta):
    with pytest.raises(ValueError):
        adaptive.bounds(eta, 2)
    with pytest.raises(ValueError):
        adaptive.example_controller(eta)


@pytest.mark.parametrize('numerator,denominator', [(True, 2), (1, True), (2, 4), (1, -2), (1, 0)])
def test_rejects_forged_fraction_internals(numerator, denominator):
    value = object.__new__(F)
    value._numerator, value._denominator = numerator, denominator
    with pytest.raises(ValueError):
        adaptive.bounds(value, 2)


@pytest.mark.parametrize('kwargs', [dict(attempts=True), dict(attempts=-1), dict(attempts=2.),
    dict(binary=1), dict(law_scope=True), dict(law_scope='coherent'),
    dict(logical_dimension=True), dict(logical_dimension=1), dict(logical_dimension=79)])
def test_rejects_wrong_bound_domains(kwargs):
    arguments = {'attempts': 2, **kwargs}
    with pytest.raises(ValueError):
        adaptive.bounds(F(1, 2), **arguments)


def test_controller_rejects_hidden_missing_or_unfair_menus_and_invalid_memory():
    controller = adaptive.example_controller()
    bad = [replace(controller, initial_memory=True), replace(controller, initial_memory=2),
           replace(controller, eta=F(0)), replace(controller, probabilities=list(controller.probabilities))]
    # Static-seam control fails exploration, even though it is a valid native
    # experiment with surviving spectator coherences outside this theorem.
    row = (F(1),)+(F(0),)*59
    bad.append(replace(controller, probabilities=(row, row)))
    wrong = list(controller.updates[0])
    wrong[0] = (wrong[0][0][:-1], wrong[0][1])
    bad.append(replace(controller, updates=(tuple(wrong), controller.updates[1])))
    row = list(controller.probabilities[0]); row[0] = float(row[0])
    bad.append(replace(controller, probabilities=(tuple(row), controller.probabilities[1])))
    for damaged in bad:
        with pytest.raises(ValueError):
            adaptive.finite_controller_kernel(damaged, 1)
    with pytest.raises(ValueError):
        adaptive.finite_controller_kernel(lambda history: (0, 1), 1)


def test_excluded_static_seam_really_preserves_a_quantum_spectator_pair():
    states = native.states()
    left = tuple(int(i in (2, 3)) for i in range(12))
    right = tuple(int(i in (2, 4)) for i in range(12))
    assert left in states and right in states
    for law in ('C', 'E'):
        for coin in (0, 1):
            assert endpoint(left, (0, 1), coin, law) == (left, 0)
            assert endpoint(right, (0, 1), coin, law) == (right, 0)
    # An existing exact one-event code is another positive control. Its
    # repeated intermediate recovery is expressly outside the passive grammar.
    for law in ('C', 'E'):
        assert record_recovery.one_event_code(law)['exact_correctability_after_one_event']


def test_returned_branch_records_cannot_poison_owned_native_maps_or_certificate():
    certificate = adaptive.exact_certificate()
    controller = adaptive.example_controller()
    before = adaptive.finite_controller_kernel(controller, 2)
    rows = adaptive.branch_maps()
    row = next(row for row in rows if len(row.outputs) >= 2)
    original = row.outputs
    row.__dict__['outputs'] = (original[1], original[0], *original[2:])
    # Exercise the second bypass of frozen dataclass assignment too.
    object.__setattr__(rows[-1], 'inputs', (0,))
    assert adaptive.branch_maps()[0].outputs == original
    assert adaptive.finite_controller_kernel(controller, 2) == before
    assert adaptive.exact_certificate() == certificate


def test_returned_rationals_and_controller_tables_have_detached_custody():
    controller = adaptive.example_controller()
    controller.eta._denominator = 3
    assert adaptive.example_controller().eta == F(1, 2)
    supplied = F(1, 2)
    result = adaptive.bounds(supplied, 2)
    result['eta']._denominator = 3
    assert supplied == F(1, 2)
    original = adaptive.example_controller()
    checked = adaptive.validate_controller(original)
    before = original.probabilities[0][0]
    expected = F(before.numerator, before.denominator)
    checked.probabilities[0][0]._numerator += 1
    assert original.probabilities[0][0] == expected
    # Shared construction zeros are detached individually too.
    checked = adaptive.validate_controller(original)
    checked.probabilities[0][30]._numerator = 1
    assert checked.probabilities[0][31] == 0
    assert original.probabilities[0][30] == 0

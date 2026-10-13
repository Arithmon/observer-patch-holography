"""Independent endpoint histories, all-input kernels and hostile records."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

import pytest

from . import record_recovery as recovery


def floor_step(law, state, edge, coin):
    a, b = edge
    total = state[a] + state[b]
    output = list(state)
    output[a], output[b] = (total // 2, (total+1) // 2)[::1 if coin == 0 else -1]
    value = state[a] - output[a] if law == "C" else state[a] - state[b]
    return tuple(output), (a, b, coin, value)


def enumerate_histories(law, states, word=None):
    grouped = defaultdict(dict)
    length = 2 if word is None else len(word)
    for i, initial in enumerate(states):
        for coins in product((0, 1), repeat=length):
            state, history = initial, []
            for tick, coin in enumerate(coins):
                if word is not None:
                    edge = word[tick]
                elif tick == 0:
                    edge = (0, 1)
                else:
                    edge = (0, 5) if initial[0] - first[0] > 0 else (0, 1)
                state, record = floor_step(law, state, edge, coin)
                if tick == 0:
                    first = state
                history.append(record)
            grouped[tuple(history)][i] = state
    return grouped


@pytest.mark.parametrize("law", ("C", "E"))
@pytest.mark.parametrize("binary", (False, True))
def test_complete_observer_and_inverse_on_every_input(law, binary):
    states = recovery.carrier(binary)
    expected = enumerate_histories(law, states)
    actual = recovery.observer_instrument(law, binary=binary)
    assert {b.history: dict(zip(b.inputs, (states[i] for i in b.outputs), strict=True))
            for b in actual} == expected
    kernel = recovery.observer_kernel(law, binary=binary)
    for i, initial in enumerate(states):
        for history, images in expected.items():
            if i in images:
                assert recovery.restore_load(law, images[i], history) == initial
        for j in range(len(states)):
            assert kernel[i][j] == F(sum(i in images and j in images
                                        for images in expected.values()), 4)


@pytest.mark.parametrize("law,expected", (("C", "1301/2028"), ("E", "3329/6084")))
def test_existing_observer_optimal_recovery(law, expected):
    report = recovery.observer_report(law)
    assert report["optimal_maximally_mixed_entanglement_fidelity"] == expected
    assert report["full_input_correctable_block_sizes"] == [56, 10, 9, 1, 1, 1]
    assert report["full_input_correctable_algebra_dimension"] == 3320


@pytest.mark.parametrize("law", ("C", "E"))
@pytest.mark.parametrize("word", ((), ((0, 1),), ((0, 1), (0, 5)), ((1, 0), (0, 5), (0, 1))))
def test_dynamic_pair_kernel_matches_complete_coin_enumeration(law, word):
    states = recovery.carrier()
    groups = enumerate_histories(law, states, word)
    actual = recovery.fixed_word_kernel(law, word)
    expected = [[F(0)] * len(states) for _ in states]
    for images in groups.values():
        for i in images:
            for j in images:
                expected[i][j] += F(1, 2**len(word))
    assert actual == tuple(map(tuple, expected))


@pytest.mark.parametrize("law", ("C", "E"))
def test_iid_rates_and_two_event_recursion_are_native(law):
    states = recovery.carrier()
    rates = recovery.iid_rates(law)
    one = [[F(0)] * len(states) for _ in states]
    two = [[F(0)] * len(states) for _ in states]
    index = {state: i for i, state in enumerate(states)}
    for edge in recovery.seams():
        for coin in (0, 1):
            images = [floor_step(law, state, edge, coin) for state in states]
            for i, (x, rx) in enumerate(images):
                for j, (y, ry) in enumerate(images):
                    if rx == ry:
                        one[i][j] += F(1, 60)
                        two[i][j] += rates[index[x]][index[y]] / 60
    assert tuple(map(tuple, one)) == rates
    assert two == [[q*q for q in row] for row in rates]


@pytest.mark.parametrize("law,rate", (("C", "5/6"), ("E", "7/10")))
@pytest.mark.parametrize("binary", (False, True))
def test_iid_fidelity_includes_every_matrix_unit(law, rate, binary):
    report = recovery.iid_report(law, binary=binary)
    assert sum(report["rate_multiplicities"].values()) == report["dimension"]**2
    assert report["largest_offdiagonal_rate"] == rate
    assert report["optimal_entanglement_fidelity"]["0"] == "1"
    assert F(report["optimal_entanglement_fidelity"]["16"]) > F(1, report["dimension"])


@pytest.mark.parametrize("law", ("C", "E"))
def test_logical_qubit_has_equal_norm_orthogonal_native_branch_images(law):
    report = recovery.one_event_code(law)
    counts = report["branch_probability_counts"]
    assert sum(F(p)*n for p, n in counts.items()) == 1
    assert report["logical_dimension"] == 2
    assert report["compiled_logical_decoder"] is False


def test_forged_contrast_does_not_invert_even_when_current_is_zero():
    state = tuple(int(v in (2, 3)) for v in range(12))
    with pytest.raises(ValueError, match="incompatible"):
        recovery.restore_load("E", state, [(0, 1, 0, -1)])


@pytest.mark.parametrize("bad", (None, True, 1, "c", "pass", [], {}))
def test_invalid_candidate_rejected(bad):
    with pytest.raises(ValueError):
        recovery.observer_instrument(bad)


@pytest.mark.parametrize("bad", (True, -1, 0.5, "2", None))
def test_invalid_count_rejected(bad):
    with pytest.raises(ValueError):
        recovery.attempts(bad)


@pytest.mark.parametrize("record", (None, (), (0, 1, True, 0), (0, 1, 2, 0),
                                   (0, 1, 0, 2), (0, 2, 0, 0)))
def test_invalid_current_record_rejected(record):
    with pytest.raises(ValueError):
        recovery.record_current("C", record)


def test_noninjective_native_mutation_is_rejected(monkeypatch):
    original = recovery.native.native

    def corrupt(state, edge, coin):
        output, (c, _) = original(state, edge, coin)
        return output, (c, 0)

    monkeypatch.setattr(recovery.native, "native", corrupt)
    with pytest.raises((ValueError, ArithmeticError)):
        recovery.observer_instrument("C")


@pytest.mark.parametrize("law", ("C", "E"))
@pytest.mark.parametrize("corruption", ("coin", "boolean_coin", "current"))
def test_native_emitted_record_is_checked_before_relabelling(monkeypatch, law, corruption):
    original = recovery.native.native

    def corrupt(state, edge, coin):
        output, (c, current) = original(state, edge, coin)
        if corruption == "coin":
            c = 1-c
        elif corruption == "boolean_coin":
            c = bool(c)
        else:
            current += 1
        return output, (c, current)

    monkeypatch.setattr(recovery.native, "native", corrupt)
    with pytest.raises(ArithmeticError, match="emitted record changed"):
        recovery.observer_instrument(law)

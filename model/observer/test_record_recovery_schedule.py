"""Independent endpoint, component, stopping and bounded-policy controls."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product

import pytest

from . import record_recovery as recovery
from . import record_recovery_schedule as schedule


def _step(law, state, edge, coin):
    # Independent floor/ceiling primitive; no production step or decoder.
    a, b = edge
    out = list(state)
    total = state[a]+state[b]
    out[a] = total//2 if coin == 0 else (total+1)//2
    out[b] = total-out[a]
    value = state[a]-out[a] if law == "C" else state[a]-state[b]
    return tuple(out), (a, b, coin, value)


def test_component_criterion_against_every_subset_of_both_forests():
    for binary in (False, True):
        edges = schedule.forest(binary=binary)
        states = recovery.carrier(binary)
        assert len(edges) == (7 if binary else 8)
        for included in product((False, True), repeat=len(edges)):
            word = tuple(edge for edge, keep in zip(edges, included, strict=True) if keep)
            # Original contrast signatures determine E history fibers for
            # every coin word; do not call the production constraints here.
            signatures = {tuple(state[a]-state[b] for a, b in word) for state in states}
            assert schedule.separates(word, binary=binary) == (len(signatures) == len(states))
        assert schedule.separates(edges, binary=binary)
    assert not schedule.separates(schedule.forest(binary=True), binary=False)


def test_fixed_word_minima_follow_from_all_component_size_partitions():
    def partitions(total, minimum=1):
        if total == 0:
            yield ()
        for first in range(minimum, total+1):
            for tail in partitions(total-first, first):
                yield (first,)+tail
    for sizes in partitions(12):
        if sum(size <= 2 for size in sizes) <= 1:
            assert len(sizes) <= 4  # At least 12-4=8 visited edges.
        if sizes.count(1) <= 1 and sizes.count(2) <= 1:
            assert len(sizes) <= 5  # At least 12-5=7 visited edges.


@pytest.mark.parametrize("binary", (False, True))
def test_every_input_decodes_under_fixed_E_and_stopped_C_without_source_read(binary):
    edges = schedule.forest(binary=binary)
    count = len(edges)
    for law in ("C", "E"):
        for pattern in range(3):
            coins = (tuple((i+pattern) % 2 for i in range(count)) if law == "E"
                     else tuple(coin for i in range(count)
                                for coin in (((i+pattern) % 2,)*(1+(i+pattern) % 4)
                                             + (1-(i+pattern) % 2,))))
            for initial in recovery.carrier(binary):
                result = schedule.execute(law, initial, coins, binary=binary)
                assert result["complete"] is True
                assert result["unused_coins"] == ()
                assert result["coin_prefix_probability"] == F(1, 2**len(coins))
                state = initial
                for record, coin in zip(result["history"], coins, strict=True):
                    state, expected = _step(law, state, record[:2], coin)
                    assert record == expected
                assert result["terminal_load"] == state
                # The reconstruction call receives no terminal or initial load.
                assert schedule.reconstruct(law, result["history"], binary=binary) == initial


def test_initial_constraints_also_cover_arbitrary_prefixes_and_reversed_seams():
    words = ((0, 1), (1, 0), (0, 5), (5, 0), (6, 8))
    for law in ("C", "E"):
        histories = defaultdict(set)
        for initial in recovery.carrier():
            for coins in product((0, 1), repeat=len(words)):
                state, history = initial, []
                for edge, coin in zip(words, coins, strict=True):
                    state, record = _step(law, state, edge, coin)
                    history.append(record)
                histories[tuple(history)].add(initial)
        for history, expected in histories.items():
            assert set(schedule.decode(law, history)) == expected


def test_opposite_coin_identity_on_every_native_seam_input_and_first_coin():
    for state in recovery.carrier():
        for a, b in recovery.seams():
            for coin in (0, 1):
                after, first = _step("C", state, (a, b), coin)
                _, same = _step("C", after, (a, b), coin)
                _, opposite = _step("C", after, (a, b), 1-coin)
                assert same[-1] == 0
                assert 2*first[-1]+opposite[-1] == state[a]-state[b]


@pytest.mark.parametrize("binary", (False, True))
def test_stopping_distribution_by_independent_coin_state_recursion(binary):
    count = len(schedule.forest(binary=binary))
    active = {(0, None): F(1)}
    stopped = F(0)
    assert schedule.stopping_mass(0, binary=binary) == 0
    assert schedule.stopping_tail(0, binary=binary) == 1
    for attempts in range(1, 6*count+1):
        following, newly_stopped = defaultdict(F), F(0)
        for (edge, first), weight in active.items():
            for coin in (0, 1):
                if first is None:
                    following[edge, coin] += weight/2
                elif first == coin:
                    following[edge, first] += weight/2
                elif edge+1 == count:
                    newly_stopped += weight/2
                else:
                    following[edge+1, None] += weight/2
        active = following
        stopped += newly_stopped
        assert newly_stopped == schedule.stopping_mass(attempts, binary=binary)
        assert sum(active.values()) == schedule.stopping_tail(attempts, binary=binary)
        assert stopped+sum(active.values()) == 1


def test_unfinished_runs_and_unused_coins_are_explicit():
    initial = recovery.carrier()[0]
    incomplete = schedule.execute("C", initial, (0,)*40)
    assert incomplete["complete"] is False
    assert incomplete["attempts"] == 40
    assert incomplete["coin_prefix_probability"] == F(1, 2**40)
    assert initial in incomplete["compatible_initial_loads"]
    assert schedule.next_seam(incomplete["history"]) == schedule.forest()[0]
    finished = schedule.execute("C", initial, (0, 1)*8+(1, 0, 0))
    assert finished["complete"] is True
    assert finished["attempts"] == 16
    assert finished["unused_coins"] == (1, 0, 0)


def test_every_nonadjacent_difference_has_a_matching_C_coin_on_every_seam():
    edges = recovery.seams()
    for u, v in combinations(range(12), 2):
        if (u, v) in edges:
            continue
        for spectator in range(12):
            left, right = [0]*12, [0]*12
            left[u] += 1
            left[spectator] += 1
            right[v] += 1
            right[spectator] += 1
            difference = tuple(a-b for a, b in zip(left, right, strict=True))
            for edge in edges:
                matches = 0
                for coin in (0, 1):
                    out_left, rec_left = _step("C", left, edge, coin)
                    out_right, rec_right = _step("C", right, edge, coin)
                    if rec_left == rec_right:
                        matches += 1
                        assert tuple(a-b for a, b in zip(out_left, out_right, strict=True)) == difference
                assert matches >= 1


def test_bounded_record_only_adaptive_policies_retain_an_entangled_branch():
    policies = (
        lambda history: (0, 1),
        lambda history: schedule.forest()[len(history) % 8],
        lambda history: (0, 5) if history and history[-1][-1] > 0 else (0, 1),
        lambda history: None if len(history) == 5 else (6, 8),
    )
    for policy in policies:
        witness = schedule.common_c_history(policy, 40)
        left, right = witness["initial_pair"]
        for a, b, coin, value in witness["history"]:
            left, rec_left = _step("C", left, (a, b), coin)
            right, rec_right = _step("C", right, (a, b), coin)
            assert rec_left == rec_right == (a, b, coin, value)
        assert left != right
        assert witness["terminal_pair"] == (left, right)
        assert witness["coin_prefix_probability"] >= F(1, 2**40)


@pytest.mark.parametrize("law,history", (
    ("E", ((0, 1, 0, 2), (0, 1, 0, 2))),
    ("C", ((0, 1, 0, 1), (0, 1, 0, 1))),
    ("C", ((0, 1, True, 0),)),
    ("E", ((0, 2, 0, 0),)),
    ("E", ((0, 1, 0, 3),)),
    ("C", ((0, 1, 0),)),
))
def test_false_or_malformed_histories_fail_closed(law, history):
    with pytest.raises(ValueError):
        schedule.decode(law, history)


@pytest.mark.parametrize("value", (True, -1, 1.0, None, "3"))
def test_stopping_attempts_are_strict_integers(value):
    for function in (schedule.stopping_mass, schedule.stopping_tail):
        with pytest.raises(ValueError):
            function(value)


def test_scope_validation_and_ambiguous_input():
    with pytest.raises(ValueError):
        schedule.forest(binary=1)
    with pytest.raises(ValueError):
        schedule.execute("C", recovery.carrier()[0], (0, True))
    with pytest.raises(ValueError):
        schedule.reconstruct("E", ())
    with pytest.raises(ValueError):
        schedule.next_seam(((0, 5, 0, 0),))
    assert schedule.audit()["carriers"]["full"]["C_stopping_mean_attempts"] == 24
    assert schedule.audit()["carriers"]["binary"]["C_stopping_mean_attempts"] == 21

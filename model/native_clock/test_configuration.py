"""Exact source enumeration, impossible records and complete-instrument controls."""
from collections import defaultdict, deque
from dataclasses import replace
from fractions import Fraction as F
import inspect
import json

import numpy as np
import pytest

from . import configuration as c


def _tape(initial, *, full=False, mode='measure', target=None, fail=None,
          seed_scans=1, retries=2):
    """Source-side test driver; its hidden state is never sent to the policy."""
    controller = c.policy(seed_scans, retries, full=full, mode=mode, target=target)
    request, state, coins, used_failure = next(controller), initial, [], False
    while True:
        currents = []
        force_zero = request.phase == fail and not used_failure
        if force_zero:
            used_failure = True
        for _ in range(request.limit):
            options = [c.native.native(state, request.seam, coin) for coin in (0, 1)]
            moving = [coin for coin, (_, record) in enumerate(options) if record[1]]
            zero = [coin for coin, (_, record) in enumerate(options) if not record[1]]
            coin = zero[0] if force_zero else moving[0] if moving else 0
            state, (_, current) = options[coin]
            coins.append(coin)
            currents.append(current)
            if current:
                break
        try:
            request = controller.send(tuple(currents))
        except StopIteration as done:
            return tuple(coins), done.value, state


def _block(state, seam, retries):
    """Independent complete source enumeration of a capped mixed/equal probe."""
    active, terminal = {state: F(1)}, defaultdict(F)
    for _ in range(retries):
        following = defaultdict(F)
        for old, weight in active.items():
            for coin in (0, 1):
                out, (_, current) = c.native.native(old, seam, coin)
                if current:
                    terminal[out, current] += weight/2
                else:
                    following[out] += weight/2
        active = following
    for out, weight in active.items():
        terminal[out, 0] += weight
    assert sum(terminal.values()) == 1
    return tuple((weight, out, current) for (out, current), weight in terminal.items())


def _opposite_path(known, start):
    paths = deque([(start,)])
    visited = {start}
    while paths:
        path = paths.popleft()
        if known[path[-1]] != known[start]:
            return path
        for edge in c.quantum.seams():
            if path[-1] in edge:
                other = next(x for x in edge if x != path[-1])
                if other in known and other not in visited:
                    visited.add(other)
                    paths.append(path+(other,))
    raise AssertionError('no opposite bit')


def _oracle(initial, retries, scans, target):
    """All first-probe, flip, second-probe and terminal outcomes, independently.

    It reads simulator loads to enumerate the complete ensemble, as an oracle
    may; production policy does not. Nothing uses the production success table,
    compressed exponent paths, or its known-bit update implementation.
    """
    success, failed, seeds = F(), F(), []
    waiting = F(1)
    for _ in range(scans):
        for seam in c.quantum.seams():
            stationary = F()
            for weight, out, current in _block(initial, seam, 1):
                if current:
                    seeds.append((waiting*weight, out, {v: out[v] for v in seam}))
                else:
                    stationary += weight
            waiting *= stationary
    failed += waiting

    def grow(weight, state, known):
        nonlocal success, failed
        if len(known) == 12:
            for seam in c.route(state, target):
                alternatives = _block(state, seam, retries)
                moving = [(p, out) for p, out, j in alternatives if j]
                failed += weight*sum(p for p, _, j in alternatives if not j)
                assert len(moving) == 1
                p, state = moving[0]
                weight *= p
            success += weight
            return
        seam = next(edge for edge in c.quantum.seams()
                    if sum(v in known for v in edge) == 1)
        vertex = next(v for v in seam if v in known)
        for probability, out, current in _block(state, seam, retries):
            mass = weight*probability
            if current:
                grow(mass, out, {v: out[v] for v in set(known).union(seam)})
                continue
            path = _opposite_path(known, vertex)
            for left, right in reversed(tuple(zip(path[:-1], path[1:], strict=True))):
                alternatives = _block(out, tuple(sorted((left, right))), retries)
                failed += mass*sum(p for p, _, j in alternatives if not j)
                moving = [(p, after) for p, after, j in alternatives if j]
                assert len(moving) == 1
                p, out = moving[0]
                mass *= p
            for p, after, j in _block(out, seam, retries):
                if j:
                    grow(mass*p, after, {v: after[v] for v in set(known).union(seam)})
                else:
                    failed += mass*p
    for entry in seeds:
        grow(*entry)
    assert success+failed == 1
    return success


def test_source_kraus_gram_current_inverse_and_mutation_isolation():
    first = c.source_audit()
    assert first['C']['fine_branches'] == 180
    assert first['E']['fine_branches'] == 300
    first['C']['exact_gram'] = 'corrupted'
    assert c.source_audit()['C']['exact_gram'] == 'I'
    for mapping in c._current_maps().values():
        assert len(mapping) == len(set(mapping.values()))


def test_all_configuration_routes_are_source_swaps_and_diameter_four():
    assert c.routing_diameter() == 4
    for target in c.native.states(binary=True):
        for initial in c.native.states(binary=True):
            state = initial
            for seam in c.route(initial, target):
                options = [out for coin in (0, 1)
                           for out, (_, j) in (c.native.native(state, seam, coin),) if j]
                assert len(options) == 1
                state = options[0]
            assert state == target


@pytest.mark.parametrize('retries,scans', [(1, 1), (3, 2)])
def test_exact_success_table_matches_independent_complete_source_oracle(retries, scans):
    target = c.native.states(binary=True)[0]
    actual = c.success_probabilities(scans, retries, mode='reset', target=target)
    expected = tuple(_oracle(x, retries, scans, target) for x in c.native.states(binary=True))
    assert actual == expected
    assert max(1-p for p in actual) <= c.budget(scans, retries, mode='reset', target=target)['failure_upper']


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('mode', ['measure', 'restore', 'reset'])
def test_every_binary_input_has_actual_rank_one_success_words(law, mode):
    target = c.native.states(binary=True)[0] if mode == 'reset' else None
    for initial in c.native.states(binary=True):
        coins, expected, state = _tape(initial, mode=mode, target=target)
        actual = c.execute(initial, coins, 1, 2, mode=mode, target=target, law=law)
        assert actual.outcome == expected
        assert actual.outcome.success and actual.outcome.initial == initial
        assert actual.output == state == actual.outcome.final
        assert c.word_support(actual) == {initial: state}
        if mode != 'measure':
            assert state == (initial if mode == 'restore' else target)
        assert len(coins) <= c.budget(1, 2, mode=mode, target=target)['native_attempts_max']


@pytest.mark.parametrize('law', ['C', 'E'])
def test_collision_readout_retains_covering_sweep_and_resets_exactly(law):
    target = c.native.states(binary=True)[0]
    for initial in (x for x in c.native.states() if max(x) == 2):
        coins, expected, _ = _tape(initial, full=True, mode='reset', target=target)
        actual = c.execute(initial, coins, 1, 2, full=True, mode='reset', target=target, law=law)
        assert actual.outcome == expected
        assert expected.success and expected.initial == initial
        assert actual.output == target
        assert c.word_support(actual, full=True) == {initial: target}
        assert any(j for _, j in expected.currents[:30])
        # Dropping the covering currents would reconstruct the absorbed state,
        # not the original collision. This is the forbidden reset/read shortcut.
        without_cover = list(target)
        for (a, b), j in expected.currents[30:]:
            without_cover[a] += j
            without_cover[b] -= j
        assert tuple(without_cover) != initial


def test_full78_success_law_uses_the_complete_absorption_distribution():
    target = c.native.states(binary=True)[0]
    binary = c.success_probabilities(2, 3, mode='reset', target=target)
    full = c.success_probabilities(2, 3, full=True, mode='reset', target=target)
    assert len(full) == 78
    for row, actual in zip(c.absorption_counts(), full, strict=True):
        assert sum(row) == 2**30
        assert actual == sum(F(n, 2**30)*p for n, p in zip(row, binary, strict=True))
    for API in (c.success_probabilities, c.budget, c.certificate):
        with pytest.raises(ValueError, match='collisions'):
            API(full=True, mode='restore')


def test_all_zero_seed_timeout_does_not_guess_or_discard_failed_work():
    initial = c.native.states(binary=True)[0]
    # A source-valid no-movement coin exists at every binary seam.
    coins = tuple(next(coin for coin in (0, 1)
                       if c.native.native(initial, seam, coin)[1][1] == 0)
                  for seam in c.quantum.seams())
    actual = c.execute(initial, coins, 1, 1)
    assert actual.output == initial
    assert not actual.outcome.success
    assert actual.outcome.reason == 'seed_timeout'
    assert actual.outcome.initial is actual.outcome.final is None
    assert actual.outcome.attempts == 30
    assert len(actual.records) == 30


def test_failed_mixed_probe_cannot_be_reinterpreted_as_an_equal_bit():
    # Find a path whose first probe is mixed and force its complete zero reply.
    for initial in c.native.states(binary=True):
        controller, state = c.policy(1, 2), initial
        request = next(controller)
        while request.phase == 'seed':
            options = [c.native.native(state, request.seam, coin) for coin in (0, 1)]
            state, (_, j) = next((x for x in options if x[1][1]), options[0])
            request = controller.send((j,))
        if state[request.seam[0]] != state[request.seam[1]]:
            coins, outcome, _ = _tape(initial, fail='probe')
            actual = c.execute(initial, coins, 1, 2)
            assert actual.outcome == outcome
            assert not outcome.success and outcome.reason == 'probe_timeout'
            assert outcome.initial is None
            return
    raise AssertionError('no mixed-probe control found')


def test_controller_inputs_are_only_current_words_and_impossible_records_reject():
    assert 'initial' not in inspect.signature(c.policy).parameters
    for bad in ((True,), (0.0,), (2,), (), (0, 0)):
        controller = c.policy(1, 1)
        next(controller)
        with pytest.raises(ValueError):
            controller.send(bad)
    # First movement on (0,1) reveals two known endpoints. Repeated source
    # currents cannot send one beyond binary support, regardless of hidden load.
    tracker = c._Tracker(False)
    request = c.Request(c.quantum.seams()[0], 1, 'test')
    tracker.accept(request, (1,))
    with pytest.raises(ValueError, match='impossible'):
        tracker.accept(request, (1,))


@pytest.mark.parametrize('bad', [True, False, 0, -1, 1.0, F(1), 129])
def test_parameters_are_strict_before_cache_lookup(bad):
    c.success_probabilities(1, 1)
    with pytest.raises(ValueError):
        c.success_probabilities(bad, 1)
    with pytest.raises(ValueError):
        c.success_probabilities(1, bad)


def test_original_flags_targets_coins_and_forged_execution_laws_reject():
    initial = c.native.states(binary=True)[0]
    for full in (0, 1, 'false'):
        with pytest.raises(ValueError):
            c.success_probabilities(full=full)
    for target in (list(initial), tuple(bool(x) for x in initial), None, (0,)*12):
        with pytest.raises(ValueError):
            c.success_probabilities(mode='reset', target=target)
    for coins in ((), (True,), (0.0,), (2,)):
        with pytest.raises(ValueError):
            c.execute(initial, coins)
    coins, _, _ = _tape(initial)
    execution = c.execute(initial, coins, 1, 2)
    for law in ('bogus', True, 0, None):
        with pytest.raises(ValueError):
            c.word_support(replace(execution, law=law))
    for record in (((0, 1), True, 0), ((0, 1), 0, True), ((0.0, 1), 0, 0)):
        with pytest.raises(ValueError):
            c.word_support(replace(execution, records=(record,)))


def test_reference_safe_reset_requires_discarding_private_input_label():
    target = c.native.states(binary=True)[0]
    values = c.success_probabilities(1, 2, mode='reset', target=target)
    # A maximally entangled input over two different configuration labels.
    weights = np.array([float(values[0]), float(values[1])])/2
    success_reference = np.diag(weights)
    assert np.trace(success_reference) < 1
    missing = np.eye(2)/2-success_reference
    failure_trace = float(np.trace(missing))
    # Complete flag comparison: half of the missing ideal success trace plus
    # half of the actual failure trace, irrespective of failure-work coherence.
    complete_half_trace = (float(np.linalg.norm(missing, ord='nuc'))+failure_trace)/2
    assert complete_half_trace <= float(1-min(values))+1e-15
    assert values[0] != values[1]  # Conditioning actually changes this reference.
    assert not np.allclose(success_reference/np.trace(success_reference), np.eye(2)/2)
    # Even the ideal reset has correlated private label/reference. Discarding
    # the label yields the original marginal, while retaining it cannot factor.
    labelled = np.diag([0.5, 0, 0, 0.5])
    assert np.linalg.norm(labelled-np.eye(4)/4, ord='nuc')/2 == 0.5


def test_default_certificates_are_exact_and_receipt_is_current():
    payload = c.run()
    assert payload['routing_diameter'] == 4
    assert payload['binary_reset']['native_attempts_max'] == 2648
    assert payload['full78_reset']['native_attempts_max'] == 2678
    for key in ('binary_measurement', 'binary_restoration', 'binary_reset', 'full78_reset'):
        row = payload[key]
        assert 0 < F(row['exact_max_failure']) <= F(row['failure_upper']) < F(1, 10**7)
    assert json.loads(c.OUTPUT.read_text(encoding='utf-8')) == payload

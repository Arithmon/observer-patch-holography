"""Independent native-coin oracle for complete directed field histories."""
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import warnings

import numpy as np
import pytest

from ..observer import native
from . import core, federation as fed


def _source_image(occupied, seam, coin):
    """Call the original endpoint rule, independently of sector assembly."""
    a, b = seam
    output, (_, current) = native.native((int(a in occupied), int(b in occupied)),
                                         (0, 1), coin)
    target = set(occupied)-{a, b}
    target.update(port for port, value in zip(seam, output, strict=True) if value)
    return tuple(sorted(target)), current


def _source_step(population, seams, denominator):
    output = defaultdict(F)
    for occupied, mass in population.items():
        for seam in seams:
            for coin in (0, 1):
                target, _ = _source_image(occupied, seam, coin)
                output[target] += mass/denominator
    return dict(output)


def _finite_oracle(data, initial, *, cap, ticks, feedback):
    """Enumerate source alternatives; no sparse matrix or grid is consulted."""
    label = tuple(int(port in initial) for port in data.read_ports)
    live = {(label, ()): {initial: F(1)}}
    terminal = {record: {} for record in fed.records(ticks)}
    for tick in range(ticks):
        following = {}
        for (label, word), state in live.items():
            waiting, stopped = state, defaultdict(F)
            for _ in range(cap):
                for target, mass in _source_step(waiting, data.incident, data.denominator).items():
                    stopped[target] += mass
                waiting = _source_step(waiting, data.remote, data.denominator)
            terminal[label, word, 'timeout'] = waiting
            for bit in (0, 1):
                selected = {occupied: mass for occupied, mass in stopped.items()
                            if int(data.read_ports[tick % 2] in occupied) == bit}
                if feedback and tick == 0 and bit:
                    selected = _source_step(selected, (data.read_ports,), 2)
                following[label, (*word, bit)] = selected
        live = following
    terminal.update({(label, word, 'complete'): state for (label, word), state in live.items()})
    return terminal


@pytest.fixture(scope='module')
def data():
    return fed.source()


def test_integer_sector_operators_match_original_coins_and_exclusion(data):
    remote, incident, repair = fed.operators()
    index = {state: i for i, state in enumerate(data.states)}
    outside = tuple(port for port in range(data.ports) if port not in data.region)
    examples = (data.read_ports,
                tuple(sorted((data.read_ports[0], outside[0]))),
                tuple(sorted(outside[:2])),
                data.remote[0])
    for initial in examples:
        for matrix, seams in ((remote, data.remote), (incident, data.incident),
                              (repair, (data.read_ports,))):
            expected = Counter(_source_image(initial, seam, coin)[0]
                               for seam in seams for coin in (0, 1))
            column = matrix[:, index[initial]].toarray().ravel()
            observed = {data.states[i]: int(column[i]) for i in np.flatnonzero(column)}
            assert observed == dict(expected)
            assert sum(observed.values()) == 2*len(seams)
    # Occupied endpoints cannot acquire a third particle or coalesce into one.
    assert len(data.states) == 28680
    assert all(len(pair) == len(set(pair)) == 2 for pair in data.states)


def test_downward_grid_action_matches_python_integer_coin_sum(data):
    scale = 1 << 51
    index = {state: i for i, state in enumerate(data.states)}
    seeds = (data.read_ports, data.remote[0])
    weights = (scale//3, scale-scale//3)
    masses = np.zeros((len(data.states), 1), dtype=np.int64)
    for seed, weight in zip(seeds, weights, strict=True):
        masses[index[seed], 0] = weight
    for matrix, seams, denominator in zip(fed.operators(),
            (data.remote, data.incident, (data.read_ports,)),
            (data.denominator, data.denominator, 2), strict=True):
        expected = defaultdict(int)
        for seed, weight in zip(seeds, weights, strict=True):
            for seam in seams:
                for coin in (0, 1):
                    expected[_source_image(seed, seam, coin)[0]] += weight
        actual = fed._floor_apply(matrix, masses, denominator).ravel()
        assert set(map(int, np.flatnonzero(actual))) <= {index[pair] for pair in expected}
        for pair, numerator in expected.items():
            assert int(actual[index[pair]]) == numerator//denominator


@pytest.mark.parametrize('feedback', (False, True))
def test_every_terminal_lower_bound_and_total_deficit_against_native_oracle(data, feedback):
    outside = next(port for port in range(data.ports) if port not in data.region)
    initial = tuple(sorted((data.read_ports[0], outside)))
    result = fed.execute(cap=1, ticks=2, grid_bits=40, feedback=feedback, initial=initial)
    oracle = _finite_oracle(data, initial, cap=1, ticks=2, feedback=feedback)
    assert tuple(oracle) == result.records
    assert sum((sum(state.values(), F()) for state in oracle.values()), F()) == 1
    index = {state: i for i, state in enumerate(data.states)}
    scale = 1 << result.grid_bits
    for column, record in enumerate(result.records):
        expected = oracle[record]
        positive = set(map(int, np.flatnonzero(result.counts[:, column])))
        assert positive <= {index[pair] for pair in expected}
        for pair, exact in expected.items():
            assert F(int(result.counts[index[pair], column]), scale) <= exact
    proof = fed.budget(result)
    stored = F(sum(map(int, result.counts.ravel())), scale)
    assert proof['complete_l1_error'] == 1-stored
    assert proof['finite_half_nuclear_error'] == (1-stored)/2
    timeout = sum((sum(population.values(), F()) for record, population in oracle.items()
                   if record[2] == 'timeout'), F())
    assert timeout == proof['exact_timeout_probability']
    # Full remote evolution is visible in the first timeout, even when the
    # tracked endpoint values are unchanged. An identity wait would fail.
    label = tuple(int(port in initial) for port in data.read_ports)
    assert len(oracle[label, (), 'timeout']) > 1


def test_complete_word_menu_and_initial_disclosure_are_not_postselected(data):
    result = fed.execute(cap=1, ticks=1, grid_bits=40)
    assert len(result.records) == 12
    assert {label for label, _, _ in result.records} == set(product((0, 1), repeat=2))
    scale = 1 << result.grid_bits
    defect = fed.budget(result)['complete_l1_error']
    exact_counts = Counter(tuple(int(p in state) for p in data.read_ports) for state in data.states)
    for label, cardinality in exact_counts.items():
        lower = F(sum(int(result.counts[:, c].sum()) for c, record in enumerate(result.records)
                      if record[0] == label), scale)
        exact = F(cardinality, len(data.states))
        assert lower <= exact <= lower+defect
    with pytest.raises(ValueError):
        fed.budget(replace(result, records=result.records[:-1]))
    mutated = result.counts.copy()
    mutated[0, 0] += 1
    with pytest.raises(ValueError):
        fed.budget(replace(result, counts=mutated))


@pytest.mark.parametrize('reader', (fed.budget, fed.field_statistics, fed.resources, fed.summary))
def test_resealed_same_mass_output_requires_actual_source_replay(reader):
    result = fed.execute(cap=1, ticks=1, grid_bits=40)
    changed = result.counts.copy()
    source = tuple(np.argwhere(changed > 0)[0])
    destination = ((source[0]+1) % changed.shape[0], source[1])
    changed[source] -= 1
    changed[destination] += 1
    fabricated = replace(result, counts=changed)
    fabricated = replace(fabricated, seal=fed._seal(fabricated))
    # The output is nonnegative, dimensionally correct, equally normalized,
    # and has a valid public integrity digest. Only source replay rejects it.
    assert int(changed.sum()) == int(result.counts.sum())
    fed._compute.cache_clear()
    with pytest.raises(ValueError, match='source replay'):
        reader(fabricated)


def test_resealed_parameters_cannot_claim_an_unexecuted_policy():
    result = fed.execute(cap=1, ticks=1, grid_bits=40)
    fabricated = replace(result, feedback=True)
    fabricated = replace(fabricated, seal=fed._seal(fabricated))
    with pytest.raises(ValueError, match='source replay'):
        fed.resources(fabricated)


def test_returned_arrays_cannot_poison_cached_execution_metadata():
    result = fed.execute(cap=1, ticks=1, grid_bits=40)
    shape, seal = result.counts.shape, result.seal
    with pytest.raises(ValueError):
        result.counts.setflags(write=True)
    # NumPy 2.5 deprecates this still-available hostile metadata mutation.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Setting the shape on a NumPy array',
                                category=DeprecationWarning)
        result.counts.shape = (result.counts.size,)
    replay = fed.execute(cap=1, ticks=1, grid_bits=40)
    assert replay.counts.shape == shape
    assert replay.seal == seal
    assert fed.validate_execution(replay).counts is not replay.counts
    with pytest.raises(ValueError):
        fed.budget(result)


def test_record_counter_uses_both_labels_and_actual_feedback_current(data):
    outside = next(port for port in range(data.ports) if port not in data.region)
    initial = tuple(sorted((data.read_ports[0], outside)))
    initial_label = tuple(int(port in initial) for port in data.read_ports)
    touching = next(seam for seam in data.incident
                    if data.read_ports[0] in seam and seam != data.read_ports)
    tape = (data.remote[0], data.read_ports, touching, data.read_ports, data.remote[-1])
    for coins in product((0, 1), repeat=len(tape)):
        state, records = initial, []
        for event, (seam, coin) in enumerate(zip(tape, coins, strict=True), 1):
            state, current = _source_image(state, seam, coin)
            if seam in data.incident:
                records.append((event, *seam, coin, current))
            assert fed.readback(initial_label, tuple(records)) == tuple(
                int(port in state) for port in data.read_ports)
    assert fed.read_clock_values((1, 0, 1, 0), feedback=True) == (1, 3, 4, 5)
    assert core.read_clock((1, 0, 1, 0), feedback=True) == (1, 3, 4, 5)


@pytest.mark.parametrize('override', ({'cap': True}, {'ticks': True}, {'grid_bits': True},
                                    {'feedback': 1}, {'grid_bits': 52}, {'initial': (True, 82)},
                                    {'initial': (82, 1)}))
def test_controller_rejects_aliases_and_unbounded_integer_arithmetic(override):
    arguments = {'cap': 1, 'ticks': 1, 'grid_bits': 40}
    arguments.update(override)
    with pytest.raises(ValueError):
        fed.execute(**arguments)

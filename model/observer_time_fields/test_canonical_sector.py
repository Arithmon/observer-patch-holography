"""Full native configuration oracles for conserved-sector field readouts."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
import json
import math
from pathlib import Path

import pytest

from . import canonical_sector as c
from . import refinement as ref
from .test_refinement import _coin_step, _interval, _records

PORTS = (0, 1, 2, 3)
EDGES = ((0, 1), (1, 2), (2, 3))
INCIDENT = ((0, 1), (1, 2))


def _oracle(population, occupied, *, ticks, cap, feedback, incident=INCIDENT):
    configurations = tuple(tuple(int(p in selected) for p in range(population))
                           for selected in combinations(range(population), occupied))
    live = {(label, None): {state: F(1, len(configurations)) for state in configurations if state[:2] == label}
            for label in product((0, 1), repeat=2)}
    terminal = {record: {} for record in _records(ticks)}
    for completed in range(ticks):
        following = {}
        for (label, first), distribution in live.items():
            arrived, waiting = _interval(distribution, EDGES, incident, cap)
            terminal[label, first, None, completed, 'timeout'] = waiting
            selected = {first: arrived} if completed else {
                bit: {state: weight for state, weight in arrived.items() if state[0] == bit}
                for bit in (0, 1)}
            for bit, distribution in selected.items():
                if completed == 0 and feedback and bit:
                    distribution = _coin_step(distribution, ((0, 1),), 2)
                if completed+1 == ticks:
                    for last in (0, 1):
                        terminal[label, bit, last, ticks, 'complete'] = {
                            state: weight for state, weight in distribution.items() if state[1] == last}
                else:
                    following[label, bit] = distribution
        live = following
    return {record: sum(distribution.values(), F(0)) for record, distribution in terminal.items()}


@pytest.mark.parametrize('population,occupied', ((4, 0), (4, 1), (4, 2), (4, 4), (6, 3)))
@pytest.mark.parametrize('feedback', (False, True))
@pytest.mark.parametrize('cap,bits', ((1, 8), (3, 40)))
def test_all_canonical_success_and_timeout_intervals_match_full_native_oracle(population, occupied, feedback, cap, bits):
    oracle = _oracle(population, occupied, ticks=3, cap=cap, feedback=feedback)
    actual = c.law(PORTS, EDGES, INCIDENT, population=population, occupied=occupied,
                   ticks=3, cap=cap, bits=bits, feedback=feedback)
    assert set(actual.intervals()) == set(oracle)
    assert sum(oracle.values()) == 1
    for record, (low, high) in actual.intervals().items():
        assert 0 <= low <= oracle[record] <= high <= 1
    assert sum(oracle[record]-low for record, low in zip(actual.records, actual.lower_mass, strict=True)) == 1-sum(actual.lower_mass)
    assert actual.timeout == 1-(1-F(1, 3)**cap)**3


@pytest.mark.parametrize('occupied', range(5))
def test_uniform_sector_is_native_stationary_and_connected(occupied):
    states = tuple(tuple(int(p in selected) for p in PORTS) for selected in combinations(PORTS, occupied))
    uniform = {state: F(1, len(states)) for state in states}
    assert _coin_step(uniform, EDGES, 6) == uniform
    reached = {states[0]}
    for _ in states:
        reached |= {target for state in tuple(reached) for target in _coin_step({state: F(1)}, EDGES, 6)}
    assert reached == set(states)
    assert all(_coin_step({state: F(1)}, EDGES, 6).get(state, F(0)) >= F(1, 2) for state in states)


def test_disconnected_graph_does_not_determine_uniform_global_sector():
    seams = ((0, 1), (2, 3))
    concentrated = {(1, 1, 0, 0): F(1)}
    assert _coin_step(concentrated, seams, 4) == concentrated
    assert concentrated != {state: F(1, 6) for state in product((0, 1), repeat=4) if sum(state) == 2}


def test_same_density_does_not_identify_product_or_wrong_population_preparations():
    canonical = c.law(PORTS, EDGES, INCIDENT, population=6, occupied=3, ticks=2, cap=3, bits=44)
    local_sector = c.law(PORTS, EDGES, INCIDENT, population=4, occupied=2, ticks=2, cap=3, bits=44)
    independent = ref.law(PORTS, EDGES, INCIDENT, rho=F(1, 2), ticks=2, cap=3, bits=44)
    def marginal(value):
        return sum(weight for record, weight in zip(value.records, value.lower_mass, strict=True)
                   if record[0] == (1, 1))
    assert abs(marginal(canonical)-F(1, 5)) < F(1, 10**9)
    assert abs(marginal(local_sector)-F(1, 6)) < F(1, 10**9)
    assert abs(marginal(independent)-F(1, 4)) < F(1, 10**9)


def test_exact_integer_draw_sampler_is_uniform_without_replacement():
    class Draws:
        def __init__(self, values):
            self.values = iter(values)
        def randrange(self, bound):
            value = next(self.values)
            assert 0 <= value < bound
            return value
    counts = Counter(c.draw_sector(Draws(values), 4, 2)
                     for values in product(range(4), range(3), range(2), range(1)))
    assert len(counts) == math.comb(4, 2)
    assert set(counts.values()) == {4}
    assert all(sum(state) == 2 for state in counts)


@pytest.mark.parametrize('feedback', (False, True))
def test_normalized_two_screen_distance_keeps_finite_number_corrections(feedback):
    laws = [c.law(PORTS, EDGES, EDGES, population=population, occupied=occupied,
                  ticks=3, cap=1, bits=44, feedback=feedback) for population, occupied in ((4, 2), (6, 3))]
    exact = [_oracle(population, occupied, ticks=3, cap=1, feedback=feedback, incident=EDGES)
             for population, occupied in ((4, 2), (6, 3))]
    distance = sum(abs(exact[0][key]-exact[1][key]) for key in exact[0])/2
    low, high = c.screen_distance(*laws, F(0))
    assert 0 < low <= distance <= high
    assert high-low < F(1, 10**8)
    expanded = c.screen_distance(*laws, F(1, 1000))
    assert expanded[0] == max(F(0), low-F(1, 500))
    assert expanded[1] == min(F(1), high+F(1, 500))


def test_law_and_export_reject_forged_density_count_records_and_arithmetic():
    good = c.law(PORTS, EDGES, INCIDENT, population=6, occupied=3, ticks=2, cap=2, bits=20)
    first = good.records[0]
    for changed in (good._replace(occupied=True), good._replace(population=4),
                    good._replace(lower_mass=(good.lower_mass[0]+F(1, 2**40), *good.lower_mass[1:])),
                    good._replace(records=(((False, first[0][1]), *first[1:]), *good.records[1:]))):
        with pytest.raises(ValueError):
            c.export(changed, F(0))


@pytest.mark.parametrize('changes', ({'population': True}, {'occupied': 2.0}, {'occupied': -1},
    {'population': 3}, {'occupied': 5}, {'feedback': 1}, {'ticks': True}))
def test_invalid_canonical_inputs_rejected_before_shared_coefficient_cache(changes):
    params = dict(population=4, occupied=2, ticks=2, cap=2, bits=16)
    c.law(PORTS, EDGES, INCIDENT, **params)
    with pytest.raises(ValueError):
        c.law(PORTS, EDGES, INCIDENT, **dict(params, **changes))


def test_full_source_fixed_sector_continuation_uses_saved_observer_memory():
    full = c.execute(1, feedback=True)
    first_slot = next(record[1] for record in full.source.records[0] if record[0] == 'scheduled')
    prefix = c.execute(1, feedback=True, stop_after_slots=first_slot)
    restored = c.checkpoint_from_json(c.checkpoint_to_json(prefix))
    continued = c.execute(1, feedback=True, checkpoint=restored)
    assert continued == full
    assert tuple(map(sum, full.source.states)) == (5120, 81920)
    assert full.source.scheduled_hits == 4
    for side in (0, 1):
        suffix = full.source.records[side][len(restored.source.records[side]):]
        assert ref._counter_view(restored.counters[side], suffix, restored.first_fields[side]) == (
            full.counters[side], full.first_fields[side])
    with pytest.raises(ValueError):
        c.execute(1, feedback=True, counts=(5121, 81920), checkpoint=restored)


def test_serialized_canonical_checkpoint_rejects_changed_source_and_memory():
    prefix = c.execute(1, feedback=True, stop_after_slots=100)
    payload = json.loads(c.checkpoint_to_json(prefix))
    for section, key, value in (('source', 'initial_totals', [5121, 81920]),
                                ('source', 'source_digest', '0'*64),
                                (None, 'initial_bits', [[False, 0], [0, 0]]),
                                (None, 'counters', [[1, 1], [1, 1]])):
        changed = json.loads(json.dumps(payload))
        (changed if section is None else changed[section])[key] = value
        with pytest.raises(ValueError):
            c.checkpoint_from_json(json.dumps(changed))
    with pytest.raises(ValueError):
        c.checkpoint_from_json('{"schema":"duplicate",'+c.checkpoint_to_json(prefix)[1:])


def test_canonical_receipt_replays_and_retains_separate_screen_restrictions():
    from tools import strict_json
    from ..observer.scalar_response.experiment import compare_receipt
    result = c.run()
    compare_receipt(strict_json.load(Path(__file__).with_name('canonical_sector_results.json')), result)
    assert not result['scope']['independent_Bernoulli_product_preparation_required']
    rows = result['screens']
    assert rows[0]['initial_distinct_port_covariance'] != rows[1]['initial_distinct_port_covariance']
    for row in rows:
        assert all(len(protocol['records']) == 44 for protocol in row['field_protocols'])
        assert F(row['field_protocols'][0]['full_screen_eventual_connected_four'][0]) > 0
        assert F(row['field_protocols'][1]['full_screen_eventual_connected_four'][1]) < 0
        assert F(row['feedback_increase_of_final_field'][0]) > F(8, 100)

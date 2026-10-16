"""Independent scheduler, source-population and hostile receipt controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import math
import subprocess
import sys

import numpy as np
import pytest

from tools import strict_json
from ..observer import native
from ..observer.scalar_response.experiment import compare_receipt
from . import incident_refinement as c


def _scheduler_reward(word, clocks, total, ticks):
    # Bellman recursion for the number of matching subsequences. This does
    # not use a binomial coefficient or assume independent waiting times
    # for occurrences of the word. Skip and match are distinct choices.
    values = {}
    for remaining in range(ticks+1):
        values[remaining, len(word)] = F(1)
    for suffix in range(len(word)):
        values[0, suffix] = F(0)
    for remaining in range(1, ticks+1):
        for suffix in reversed(range(len(word))):
            reward = F(0)
            for letter in range(total):
                next_ticks = remaining-int(letter in clocks)
                if letter in clocks:
                    reward += values[next_ticks, suffix]/total
                if letter == word[suffix]:
                    reward += values[next_ticks, suffix+1]/total
            values[remaining, suffix] = reward/F(len(clocks), total)
    return values[ticks, 0]


def test_exact_stopped_word_formula_against_independent_scheduler_recursion():
    for total, clocks in ((3, {0}), (3, {0, 1}), (5, {0, 1})):
        for length in range(1, 5):
            for word in product(range(3), repeat=length):
                for ticks in (1, 2, 3):
                    expected = _scheduler_reward(word, clocks, total, ticks)
                    actual = c.word_expectation(ticks, len(clocks), length,
                                                sum(e in clocks for e in word[:-1]))
                    assert actual == expected
    # Two clock letters cannot both occur before the first clock tick;
    # ignoring the strict gap would give an incorrect positive expectation.
    assert c.word_expectation(1, 2, 2, 1) == 0
    assert c.word_expectation(1, 2, 2, 0) == F(1, 4)


def test_first_exit_walks_independently_enumerated_not_all_walks_or_dead_boundary():
    ports = (0, 1, 2)
    edges = ((0, 1), (1, 2))
    counts = c.exit_counts(ports, edges, (0,), (0, 0, 2), 8)
    adjacency = {0: (1,), 1: (0, 2), 2: (1, 3, 4)}
    paths = [(0,)]
    oracle = []
    for _ in range(8):
        following, escapes = [], 0
        for path in paths:
            for q in adjacency[path[-1]]:
                if q > 2:
                    escapes += 1
                else:
                    following.append((*path, q))
        paths = following
        oracle.append(escapes)
    assert counts == tuple(oracle) == (0, 0, 2, 0, 4, 0, 8, 0)


@pytest.fixture(scope='module')
def source_data():
    return c.geometry()


def test_source_support_and_longer_positive_series_fit_inside_certified_tail(source_data):
    ports, edges, incident, observed, exits, mapping = source_data
    assert len(observed) == 15 and len(incident) == 33 and sum(exits) == 45
    assert len(mapping) == 222
    counts = c.exit_counts(ports, edges, observed, exits, 64)
    assert counts[:10] == (0,)*8+(384, 4992)
    for ticks in (1, 2, 4, 8, 16):
        bound, remainder = c.exit_bound(counts[:48], ticks, 33)
        more, _ = c.exit_bound(counts, ticks, 33)
        partial = sum((n*c.word_expectation(ticks, 33, j)
                       for j, n in enumerate(counts, 1)), F(0))
        assert 0 < partial <= more <= bound < F(1, 10000)
        assert bound-more <= remainder
    # In particular the eight-tick *pair* comparison is below 4.3e-7.
    assert 2*c.exit_bound(counts[:48], 8, 33)[0] < F(43, 100000000)


def _native_population_return(edges, incident, cap, ticks):
    # Original endpoint/placement rule on ALL 16 binary configurations,
    # rather than the producer's one-particle (I+swap)/2 density reduction.
    rho = F(1, 3)
    states = {bits: rho**sum(bits)*(1-rho)**(4-sum(bits))*bits[0]
              for bits in product((0, 1), repeat=4)}
    for _ in range(ticks):
        waiting, completed = states, {}
        for _ in range(cap):
            next_wait = {}
            for bits, mass in waiting.items():
                for edge in edges:
                    for coin in (0, 1):
                        out, _ = native.native(bits, edge, coin)
                        target = completed if edge in incident else next_wait
                        target[out] = target.get(out, F(0))+mass/F(2*len(edges))
            waiting = next_wait
        states = completed
    return sum((mass*bits[0] for bits, mass in states.items()), F(0))


def test_density_covariance_against_all_configuration_native_endpoint_oracle():
    ports, edges, incident = (0, 1, 2, 3), ((0, 1), (1, 2), (2, 3)), ((0, 1),)
    result = c.density_response(ports, edges, incident, ticks=2, cap=4, bits=24)
    rho = F(1, 3)
    for row in result:
        ticks = row['scheduled_ticks']
        lo, hi = map(F, row['return_probability_interval'])
        timeout = F(row['unresolved_wait_allowance'])
        joint = _native_population_return(edges, incident, 4, ticks)
        # The independent finite-cap column has mass 1-timeout. Remove
        # exactly that cap probability before comparing the density column.
        finite_return = (joint-rho*rho*(1-timeout))/(rho*(1-rho))
        arithmetic = F(row['rounding_allowance'])
        assert lo <= finite_return <= lo+arithmetic
        assert hi-lo >= timeout
        assert arithmetic > 0


def test_directed_floor_each_component_against_python_integer_oracle():
    ports = (0, 1, 2, 3)
    edges = ((0, 1), (1, 2))
    raw = (2**52-19, 7, 8, 4)
    vector = np.array(raw, dtype=np.int64)
    expected = []
    for p in ports:
        value = 0
        for a, b in edges:
            swapped = raw[b] if p == a else raw[a] if p == b else raw[p]
            value += raw[p]+swapped
        expected.append(value//14)
    assert list(c._density_step(vector, edges, ports, 7)) == expected
    assert sum(expected) <= sum(raw)


def test_first_native_clock_field_response_has_closed_form(source_data):
    ports, edges, incident, _, _, _ = source_data
    result = c.density_response(ports, edges, incident, 1)[0]
    lo, hi = map(F, result['return_probability_interval'])
    # Remote seams never touch root zero before the first incident. Six of
    # 33 incident seams touch it and the source coin moves with probability1/2.
    exact = 1-F(6, 2*33)
    assert lo <= exact <= hi
    assert hi-lo < F(1, 10**11)


@pytest.mark.parametrize('bad', (True, -1, 0, 1.0, '1', None))
def test_public_clock_counts_refuse_coerced_or_empty_inputs(bad):
    with pytest.raises(ValueError):
        c.word_expectation(bad, 33, 9)
    with pytest.raises(ValueError):
        c.exit_bound((1,), bad, 33)
    with pytest.raises(ValueError):
        c.density_response((0, 1), ((0, 1),), ((0, 1),), bad)


def test_density_graph_domain_and_tail_refusals():
    for edges, incident in ((((0, 1), (0, 1)), ((0, 1),)),
                            (((0, 1),), ((0, True),)),
                            (((0, 1),), ((0, 2),))):
        with pytest.raises(ValueError):
            c.density_response((0, 1), edges, incident, 1)
    with pytest.raises(ValueError):
        c.exit_bound((1,), 8, 1)
    with pytest.raises(ValueError):
        c.exit_bound((True,), 1, 33)


def test_external_cutoff_bound_encloses_exact_small_binomial_tails():
    for slots, ticks, incident, total in product((1, 4, 12, 50), (1, 2, 4), (1, 2), (2, 5)):
        p = F(incident, total)
        exact = sum((math.comb(slots, j)*p**j*(1-p)**(slots-j)
                     for j in range(min(ticks, slots+1))), F(0))
        assert exact <= c.slot_timeout_bound(slots, ticks, incident, total) <= 1
    assert c.slot_timeout_bound(500000, 4, 33, 684879) < F(1, 10000000)


def test_retained_rational_receipt_and_material_forgery_rejected():
    actual = c.run()
    retained = strict_json.load(Path(__file__).with_name('incident_results.json'))
    compare_receipt(retained, actual)
    def exact_types(value):
        # The shared float-diagnostic tolerance must not certify altered
        # reported errors above this packet's fixed 1e-12 allowance.
        assert type(value) in (dict, list, str, int, bool)
        if type(value) is dict:
            for item in value.values():
                exact_types(item)
        elif type(value) is list:
            for item in value:
                exact_types(item)
    exact_types(actual)
    for mutate in (
        lambda r: r['geometry'].__setitem__('record_and_field_support_ports', 12),
        lambda r: r['complete_stopped_quantum_history_comparison'][2].__setitem__(
            'two_screen_complete_output_allowance', '0'),
        lambda r: r['finite_density_fields']['responses'][3].__setitem__(
            'either_full_screen_connected_field_interval', ['0', '0']),
        lambda r: r['scope'].__setitem__('field_covariance_execution_is_diagonal_sector', 1),
        lambda r: r['complete_quantum_clock_controls'][0].__setitem__(
            'fixed_probe_half_nuclear_error_upper', '1/10000000000'),
    ):
        wrong = deepcopy(retained)
        mutate(wrong)
        with pytest.raises(ValueError):
            compare_receipt(wrong, actual)


def test_real_cli_rejects_modified_clock_certificate(tmp_path):
    wrong = strict_json.load(Path(__file__).with_name('incident_results.json'))
    wrong['geometry']['incident_seams'] = 30
    path = tmp_path/'wrong.json'
    path.write_text(c.canonical(wrong), encoding='utf-8')
    process = subprocess.run([sys.executable, '-m', 'model.native_action.incident_refinement',
                              '--check', str(path)], text=True, capture_output=True)
    assert process.returncode != 0
    assert 'Receipt mismatch' in process.stderr and 'incident_seams' in process.stderr

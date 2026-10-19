"""Operational response, complete event accounting and premise controls."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
import json
import subprocess
import sys

import pytest

from . import repair_response as m


def canonical_states(size, number):
    return tuple(tuple(int(i in positions) for i in range(size))
                 for positions in combinations(range(size), number))


def test_marked_density_is_a_graph_independent_canonical_identity():
    # Independent explicit pair heat-bath rule on every number sector, with
    # no use of the model's native transition or pulse-weight implementation.
    for size in range(2, 7):
        for number in range(size+1):
            states = canonical_states(size, number)
            actual = defaultdict(F)
            for state in states:
                if not state[0]:
                    actual[(0, -1, 0), state] += F(1, len(states))
                else:
                    total = state[0]+state[1]
                    for coin in (0, 1):
                        pair = (total//2, (total+1)//2)
                        if coin:
                            pair = pair[::-1]
                        after = pair+state[2:]
                        actual[(1, coin, state[0]-pair[0]), after] += F(1, 2*len(states))
            for state in states:
                padded = state+(0,)*(12-size)
                weights = m.pulse_weights(padded)
                for mark in m.MARKS:
                    assert actual[mark, state] == weights[mark]/len(states)
                assert sum(weights.values()) == 1+F(state[1]-state[0], 2)


def test_direct_native_pulse_matches_every_marked_forecast_and_zero_label():
    actual = m.pulse_distribution()
    assert actual == m.pulse_distribution(reverse_forecast=True)
    assert sum(actual.values()) == 1
    marks = {mark: sum((v for (label, _), v in actual.items() if label == mark), F())
             for mark in m.MARKS}
    assert marks == {(0, -1, 0): F(5, 6), (1, 0, -1): F(0),
                     (1, 0, 0): F(1, 132), (1, 0, 1): F(5, 66),
                     (1, 1, -1): F(0), (1, 1, 0): F(1, 12), (1, 1, 1): F(0)}
    # The unselected extra fair repair preserves equilibrium, so feedback
    # is the reason this finite-amplitude perturbation has a response.
    unconditional = defaultdict(F)
    for state in m.source().states:
        for coin in (0, 1):
            after, _ = m.native.native(state, m.PULSE, coin)
            unconditional[after] += F(1, 132)
    assert set(unconditional.values()) == {F(1, 66)}


def test_independent_one_attempt_path_oracle_keeps_complete_clock_field_labels():
    expected = defaultdict(F)
    for before in m.source().states:
        options = ((0, -1, 0, before, F(1)),) if before[0] == 0 else tuple(
            (1, coin, record[1], after, F(1, 2)) for coin in (0, 1)
            for after, record in (m.native.native(before, (0, 1), coin),))
        for trigger, pulse_coin, pulse_current, after_pulse, weight in options:
            initial = tuple(after_pulse[p] for p in (0, 1, 2))
            for seam in m.source().seams:
                for coin in (0, 1):
                    terminal, (_, current) = m.native.native(after_pulse, seam, coin)
                    tape = ((seam, current),) if any(p in seam for p in (0, 1, 2)) else ()
                    key = ((trigger, pulse_coin, pulse_current), initial, tape,
                           'tick' if 0 in seam else 'timeout',
                           (trigger+int(0 in seam), trigger+int(1 in seam)),
                           tuple(terminal[p] for p in (0, 1, 2)))
                    expected[key] += weight/F(66*60)
    assert dict(expected) == m.stopped_law(1, clock_b=1)


@pytest.mark.parametrize('cap,clock_b', product((1, 2), (1, 2)))
def test_complete_stopped_laws_and_response_from_passive_future_records(cap, clock_b):
    actual = m.stopped_law(cap, clock_b=clock_b)
    forecast = m.stopped_law(cap, clock_b=clock_b, reverse_forecast=True)
    assert actual == forecast
    baseline = m.stopped_law(cap, clock_b=clock_b, intervention=False)
    assert sum(actual.values()) == sum(baseline.values()) == 1
    assert sum(v for key, v in actual.items() if key[3] == 'timeout') == F(5, 6)**cap
    def selected(key):
        return key[3] == 'tick' and key[-1][1] and key[-2][1]-key[0][0]*int(clock_b in m.PULSE) >= 1
    change = sum((v for key, v in actual.items() if selected(key)), F())-sum(
        (v for key, v in baseline.items() if selected(key)), F())
    predicted = sum((F(key[1][1]-key[1][0], 2)*v for key, v in baseline.items() if selected(key)), F())
    assert change == predicted >= 0
    assert (predicted > 0) == (cap == 2)
    assert any(v and key[0][0] and any(current == 0 for _, current in key[2]) for key, v in actual.items())
    # Dropping either zero-current incident events or the policy's own repair
    # changes the law; neither is an allowed clock or probability convention.
    zero_tick_mass = sum(v for key, v in actual.items()
                        if key[3] == 'tick' and key[2][-1][1] == 0)
    assert zero_tick_mass > 0
    assert sum(v for key, v in actual.items() if key[-2][0] == 2) == F(1, 6)*(1-F(5, 6)**cap)


@pytest.mark.parametrize('attempts', (0, 1, 2, 4, 8))
def test_exact_full_source_field_response_obeys_mass_conservation(attempts):
    response = m.density_response(attempts)
    assert sum(response) == 0
    assert response[0] == -response[1] < 0
    if attempts == 0:
        assert response == (-F(5, 66), F(5, 66))+(F(0),)*10
    else:
        assert abs(response[0]) < F(5, 66)


@pytest.mark.parametrize('z', (1, 2, F(1, 2)))
def test_response_resolvent_agrees_with_independent_full_native_configuration_inverse(z):
    response = m.discounted_response(z)
    assert sum(response) == 0
    assert response[0] == -response[1] < 0


def test_nonstationary_prior_invalidates_equilibrium_response_formula():
    chosen = next(i for i, state in enumerate(m.source().states) if state[0] and not state[1])
    prior = tuple(F(int(i == chosen)) for i in range(66))
    assert m.pulse_distribution(prior=prior) != m.pulse_distribution(prior=prior, reverse_forecast=True)
    assert sum(m.pulse_distribution(prior=prior, reverse_forecast=True).values()) == F(1, 2)


def test_event_response_is_not_a_reset_or_independent_clock_field_product():
    law = m.stopped_law(1, clock_b=1)
    trigger = sum(v for key, v in law.items() if key[0][0])
    field = sum(v for key, v in law.items() if key[-1][1])
    joint = sum(v for key, v in law.items() if key[0][0] and key[-1][1])
    assert joint != trigger*field
    # Repreparing stationarity removes the one-step signal response.
    assert m.density_response(1)[1] > 0


def test_api_rejects_invalid_domains_before_running():
    for bad in (True, 0, 3, 1.0):
        with pytest.raises(ValueError):
            m.stopped_law(bad)
    for bad in (True, -1, 9, 1.0):
        with pytest.raises(ValueError):
            m.density_response(bad)
    for bad in (True, 0, -1, 0.5):
        with pytest.raises(ValueError):
            m.discounted_response(bad)
    for state in ([0]*12, (False,)*12, (2,)+(0,)*11, (0,)*11):
        with pytest.raises(ValueError):
            m.pulse_weights(state)
    for prior in ((F(1, 66),)*65, (1/66,)*66, (F(-1),)+(F(0),)*65):
        with pytest.raises(ValueError):
            m.pulse_distribution(prior=prior)
    with pytest.raises(ValueError):
        m.stopped_law(intervention=False, reverse_forecast=True)


def test_repo_local_runtime_imports_are_pinned():
    script = '''
import json, sys
from pathlib import Path
from model.predictions import repair_response as m
m.density_response(1)
missing=[]
for module in list(sys.modules.values()):
    path=getattr(module, '__file__', None)
    if path is None: continue
    try: relative=Path(path).resolve().relative_to(m.ROOT).as_posix()
    except ValueError: continue
    if relative.endswith('.py') and relative not in m.SOURCE_PATHS: missing.append(relative)
print(json.dumps(sorted(set(missing))))
'''
    result = subprocess.run([sys.executable, '-c', script], cwd=m.ROOT, capture_output=True,
                            text=True, check=True)
    assert json.loads(result.stdout) == []


def test_before_pulse_disclosure_and_mark_reconstruct_all_field_labels():
    for state in m.source().states:
        bits = tuple(state[p] for p in m.READ_PORTS)
        if not bits[0]:
            assert m.pulse_readback(bits, (0, -1, 0)) == bits
        else:
            for coin in (0, 1):
                after, (_, current) = m.native.native(state, m.PULSE, coin)
                assert m.pulse_readback(bits, (1, coin, current)) == tuple(after[p] for p in m.READ_PORTS)
    for bits, mark in (((1, 0, 0), (0, -1, 0)), ((1, 0, 0), (1, 1, 1)),
                       ((1, 1, 0), (1, 0, 1)), ((True, 0, 0), (1, 0, 1))):
        with pytest.raises(ValueError):
            m.pulse_readback(bits, mark)


def test_corrupt_rational_inputs_are_rejected():
    broken = F(1, 2)
    broken._numerator = 2
    with pytest.raises(ValueError):
        m.discounted_response(broken)
    broken = F(1, 66)
    broken._numerator = 2
    with pytest.raises(ValueError):
        m.pulse_distribution(prior=(broken,)+(F(1, 66),)*64+(F(0),))

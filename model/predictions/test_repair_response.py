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


PROFILE = (F(1), F(-1))+(F(0),)*10


def test_control_changes_execution_probability_and_preserves_native_endpoint():
    matrices = m.controlled_matrices(PROFILE, F(1, 8))
    reference = m.controlled_matrices(PROFILE)
    import sympy as sp
    native = sp.Matrix(m.source().counts)/60
    assert reference['embedded'] == native
    assert reference['slot'] == (sp.eye(66)+native)/2
    # Direct seam and coin enumeration is an independent stationary-flux
    # check. The execution subkernel is not itself a Markov kernel.
    stationary_flux = defaultdict(F)
    for state in m.source().states:
        for a, b in m.source().seams:
            p = (1+F(1, 8)*(PROFILE[b]-PROFILE[a])*(state[a]-state[b]))/2
            assert F(1, 4) <= p <= F(3, 4)
            for coin in (0, 1):
                pair = ((state[a]+state[b])//2, (state[a]+state[b]+1)//2)
                if coin:
                    pair = pair[::-1]
                after = list(state); after[a], after[b] = pair
                stationary_flux[tuple(after)] += p/F(66*60)
    assert set(stationary_flux.values()) == {F(1, 132)}
    assert min(matrices['activity']) < F(1, 2) < max(matrices['activity'])
    assert reference['slot_derivative_at_zero'] != reference['embedded_derivative_at_zero']
    shifted = tuple(h+F(13, 7) for h in PROFILE)
    assert m.controlled_matrices(shifted, F(1, 8)) == matrices


@pytest.mark.parametrize('count', (0, 1, 2, 4))
def test_controlled_response_uses_actual_executed_clock(count):
    executed = m.controlled_response(PROFILE, count)
    slots = m.controlled_response(PROFILE, 2*count, clock='slot')
    assert sum(executed['terminal_probability_derivative']) == 0
    assert sum(executed['occupation_mean_derivative']) == 0
    if count:
        assert executed['occupation_mean_derivative'][0] > 0
        assert executed != slots
    else:
        assert set(executed['occupation_mean_derivative']) == {F(0)}
    assert m.controlled_response(tuple(-h for h in PROFILE), count) == {
        key: tuple(-v for v in values) for key, values in executed.items()}


def test_complete_controlled_paths_have_exact_likelihood_and_current_readback():
    initial = next(s for s in m.source().states if s[0] and not s[1])
    baseline = m.controlled_history(PROFILE, F(0), initial, ticks=1, slot_cap=2)
    epsilon = F(1, 8)
    changed = m.controlled_history(PROFILE, epsilon, initial, ticks=1, slot_cap=2)
    assert baseline.keys() == changed.keys()
    for key, mass in baseline.items():
        tape, status, count, terminal = key
        current_state = list(initial); likelihood = F(1)
        for (a, b), active, coin, current in tape:
            score = (PROFILE[b]-PROFILE[a])*(current_state[a]-current_state[b])
            likelihood *= 1+epsilon*score*(2*active-1)
            if not active:
                assert coin == -1 and current == 0
            else:
                assert coin in (0, 1)
            current_state[a] -= current; current_state[b] += current
        assert tuple(current_state) == terminal
        assert sum(row[1] for row in tape) == count
        assert (status == 'complete') == (count == 1)
        assert changed[key] == mass*likelihood
    assert sum(changed.values()) == 1
    assert any(key[1] == 'timeout' and mass for key, mass in changed.items())
    assert any(any(row[1] and not row[3] for row in key[0]) for key in changed)


@pytest.mark.parametrize('ticks,cap', ((1, 1), (1, 4), (2, 8)))
def test_stopped_likelihood_keeps_correlated_timeout_and_uniform_error(ticks, cap):
    control = m.likelihood_control(PROFILE, F(1, 8), ticks=ticks, slot_cap=cap)
    assert control['likelihood_mean'] == 1
    assert control['baseline_score_mean'] == 0
    assert control['linear_likelihood_remainder_second_moment'] <= control[
        'linear_likelihood_remainder_polynomial_upper'] <= control[
        'linear_likelihood_remainder_rational_upper']
    assert control['field_nonlinear_error']**2 <= control[
        'baseline_terminal_field_variance']*control['linear_likelihood_remainder_second_moment']
    assert control['biased_timeout']**2 <= control[
        'likelihood_second_moment']*control['baseline_timeout']
    if cap > 1:
        # A borrowed unbiased watchdog probability would miss this correlation.
        assert control['biased_timeout'] != control['baseline_timeout']
        assert control['field_nonlinear_error'] != 0
    assert control['capped_field_derivative'] != control['unbounded_clock_field_derivative']
    assert (control['capped_field_derivative']-control[
        'unbounded_clock_field_derivative'])**2 <= control['capped_derivative_error_squared_upper']


def test_native_coin_service_is_charged_and_does_not_tick_designated_dynamical_clock():
    state = m.source().states[0]
    for coins in ((), (0,), (1,), (0, 1, 1, 0)):
        result = m.native_service_coin_word(state, coins)
        assert result['terminal_load'] == state
        assert result['native_service_events'] == len(coins)
        assert result['dynamical_clock_ticks'] == 0
        assert result['branch_probability'] == F(1, 2**len(coins))
        assert len(result['records']) == len(coins)
    for bits in (1, 3, 11):
        result = m.dyadic_activation(F(2, 7), bits)
        assert 0 <= F(2, 7)-result['probability'] < result['absolute_error_upper']
        assert result['native_service_events_if_equal_seam_realized'] == bits
        assert result['uniform_seam_selection_service_is_additional']


@pytest.mark.parametrize('level', (0, 1, 2, 8, 100))
def test_uniform_finite_amplitude_continuum_signal_and_resource_budget(level):
    row = m.continuum_response_budget(level)
    assert row['executed_dynamical_repairs_at_tau_one'] == 1260*16**level
    assert row['selection_slot_cap'] == 4*row['executed_dynamical_repairs_at_tau_one']
    assert row['score_variance_upper'] < 36
    assert row['liminf_mean_response_lower'] == F(244359, 1024750000000000) > 0
    assert not row['complete_nonlinear_continuum_law_identified']
    assert not row['clock_is_a_primitive_port_incident_clock']
    bits = row['dyadic_activation_bits_per_slot']
    denominator = 16*row['ports']*row['selection_slot_cap']**2
    assert 2**(bits-1) < denominator <= 2**bits
    assert row['all_bounds_include_terminal_fields_on_timeout']


def test_controlled_inputs_reject_changed_domains_and_unsafe_activation():
    state = m.source().states[0]
    for profile in (list(PROFILE), PROFILE[:-1], (True,)+PROFILE[1:]):
        with pytest.raises(ValueError):
            m.activation(state, m.PULSE, profile, F(1, 8))
    for bad_state in ([*state], (False,)+state[1:], (0,)*12):
        with pytest.raises(ValueError):
            m.activation(bad_state, m.PULSE, PROFILE, F(1, 8))
    with pytest.raises(ValueError):
        m.activation(state, m.PULSE, PROFILE, F(1))
    with pytest.raises(ValueError):
        m.controlled_response(PROFILE, 1, clock='port-incident')
    for bad in (True, -1, 101, 1.0):
        with pytest.raises(ValueError):
            m.continuum_response_budget(bad)
    for bad in (F(0), F(-1), F(1)):
        with pytest.raises(ValueError):
            m.continuum_response_budget(0, bad)

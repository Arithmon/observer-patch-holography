"""One-history confidence, complete unknown outcomes, and source controls."""
from fractions import Fraction as F
import json

import pytest
import sympy as sp

from . import clock_bridge
from . import trajectory_readout as t


def test_poisson_certificate_lifts_to_independent_full_native_stopped_chain():
    native, blocks, _, _ = t.bits.source()
    full = clock_bridge.marked_source(native, 0)
    q = full['Q0']+full['Q1']
    class_of = {i: b for b, block in enumerate(blocks) for i in block}
    f = [state[0] for state in native.states]
    for name, c in t.source()[3].items():
        g = [c['poisson_solution'][class_of[i]] for i in range(66)]
        h = lambda i, j: (f[i] if name == 'occupation' else f[i]*f[j]
                         if name == 'adjacent_product' else abs(f[j]-f[i]))
        for i in range(66):
            assert sum(q[i, j]*(h(i, j)-c['stationary_mean']+g[j]-g[i])
                       for j in range(66)) == 0
            row = [h(i, j)+g[j] for j in range(66) if q[i, j]]
            assert max(row)-min(row) <= c['increment_range']
    assert t.source()[3]['occupation']['boundary_span'] == 3
    assert t.source()[3]['occupation']['increment_range'] == 3


def test_stationary_joint_law_is_not_independent_bit_sampling():
    c = t.source()[3]
    assert c['occupation']['stationary_mean'] == F(1, 6)
    assert c['adjacent_product']['stationary_mean'] == F(1, 11)
    assert c['movement']['stationary_mean'] == F(5, 33)
    assert F(1, 11) != F(1, 6)**2
    # Continuity fixes this relation, not independence between the moments.
    assert 2*(F(1, 6)-F(1, 11)) == F(5, 33)


def test_martingale_decomposition_on_an_actual_native_hidden_path():
    import random
    native, blocks, _, _ = t.bits.source()
    class_of = {native.states[i]: b for b, block in enumerate(blocks) for i in block}
    state = native.states[0]
    hidden = [class_of[state]]
    currents = []
    rng = random.Random(1249)
    while len(currents) < 80:
        seam = rng.choice(native.seams)
        state, (_, current) = t.bits.transfer.native.native(state, seam, rng.randrange(2))
        if 0 in seam:
            currents.append(-current if seam[0] == 0 else current)
            hidden.append(class_of[state])
    _, _, rewards, certs = t.source()
    result = t.estimate(currents)
    assert result['resolved']
    for name, c in certs.items():
        h, g, mu = rewards[name], c['poisson_solution'], c['stationary_mean']
        total = sum(h[a, b]-mu for a, b in zip(hidden, hidden[1:]))
        martingale = sum(h[a, b]-mu+g[b]-g[a] for a, b in zip(hidden, hidden[1:]))
        assert total == g[hidden[0]]-g[hidden[-1]]+martingale
        assert F(result['observables'][name]['observed_mean']) == (total+80*mu)/80


def test_rational_radius_has_the_declared_exponential_budget():
    for n in (1, 2, 13, 32768, 10**6):
        for s in (1, 12, 30):
            for name, c in t.source()[3].items():
                error = t.radius(n, s, observable=name)
                ticks = (n*error-c['boundary_span'])/c['increment_range']
                assert ticks.denominator == 1
                assert 2*ticks*ticks >= n*(s+1)
                assert ticks == 0 or 2*(ticks-1)**2 < n*(s+1)
    assert t.radius(10**6, 12) < F(1, 100)


def test_unknown_outcomes_are_kept_and_fixed_seam_control_is_detectable():
    result = t.estimate([0]*32768)
    assert not result['resolved']
    for name in ('occupation', 'adjacent_product'):
        assert result['observables'][name]['observed_mean'] is None
        assert result['observables'][name]['confidence_interval'] == ['0', '1']
    assert not result['observables']['movement']['prediction_in_interval']
    assert result['simultaneous_failure_upper'] == '3/4096'


def test_no_disclosure_is_needed_and_malformed_tapes_fail():
    result = t.estimate([0, -1, 0, 1, 0])
    assert result['observables']['occupation']['observed_mean'] == '3/5'
    assert result['observables']['adjacent_product']['observed_mean'] == '2/5'
    assert result['observables']['movement']['observed_mean'] == '2/5'
    for tape in ([], [True], [0.0], [1, 1], [-1, -1], '010'):
        with pytest.raises(ValueError):
            t.estimate(tape)
    for n, s in ((0, 12), (True, 12), (10, 0), (10, 2.0)):
        with pytest.raises(ValueError):
            t.radius(n, s)


def test_retained_nonstationary_execution_and_receipt():
    result = t.run()
    assert result == json.loads(t.OUTPUT.read_text(encoding='utf-8'))
    assert result['baseline']['all_predictions_retained']
    assert not result['fixed_empty_seam_control']['all_predictions_retained']
    assert result['ports_checked'] == 12

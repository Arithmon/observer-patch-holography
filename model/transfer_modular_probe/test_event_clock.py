"""Native event triggers, complete continuing outputs and finite clock reading."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import math
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

from . import construction as c
from . import event_clock as e


@pytest.mark.parametrize('cap', (1, 2))
def test_every_native_stopped_terminal_map_matches_independent_endpoint_paths(cap):
    _, maps, timeout, denominator = e.stopped_tables(cap)
    direct, direct_timeout = e.endpoint_oracle(cap)
    np.testing.assert_array_equal(timeout, direct_timeout)
    for mark in maps:
        np.testing.assert_array_equal(maps[mark], direct[mark])
    assert np.all(sum(maps.values()).sum(axis=1) == denominator-50**cap)
    assert np.all(timeout.sum(axis=1) == 50**cap)
    # Dropping timeout and renormalizing would change the declared experiment.
    assert F(int(sum(maps.values())[0].sum()), denominator) < 1


def test_every_incident_record_including_zero_current_triggers_same_derived_word():
    _, _, marks = e.native_tables()
    assert len(marks) == 15
    for mark in marks:
        assert e.triggered_word(mark) == c.diagonal_word(1)
    with pytest.raises(ValueError):
        e.triggered_word(((2, 3), 0))  # Not in the clock-carrier region.
    with pytest.raises(ValueError):
        e.triggered_word(('probe_source_pulse', 0))
    with pytest.raises(ValueError):
        e.triggered_word(((0, 1), False))


@pytest.mark.parametrize('mode', (1, 2))
@pytest.mark.parametrize('current', (-1, 0, 1))
def test_probe_continuation_matches_validated_shared_source_program(mode, current):
    program = c.program(mode, (1, 1))
    if current > 0:
        intervention = c.runtime.Event('native_current_triggered_balanced',
                                      instruments=(e._gates(mode)[0],))
        program = c.runtime.Program('current_controlled_modular_clock', 4,
            program.initial_density, (*program.events[:5], intervention, *program.events[5:]))
    history = c.factor_execute(program)
    _, first, second = e.probe_outputs(mode, current=current)
    for branch in history.branches:
        vector = second[branch.record]
        np.testing.assert_allclose(branch.density, np.outer(vector, vector.conj()), atol=2e-15, rtol=0)
    independent = e.density_oracle(mode, current=current)
    for mark, vector in second.items():
        np.testing.assert_allclose(np.outer(vector, vector.conj()), independent[2][mark], atol=5e-14, rtol=0)
    assert len(first) == 4 and len(second) == 16


def test_complete_native_probe_output_retains_padding_timeouts_and_terminal_fields():
    actual = e.complete_output(2)
    oracle = e.complete_output(2, ideal=True, oracle=True)
    assert len(actual) == 3661
    assert e._joint_diagnostic(actual, oracle) < 1e-13
    assert e.density_oracle_diagnostic(actual, 2) < 1e-13
    assert abs(sum(map(e._probability, actual))-1) < 1e-13
    assert any(e._probability(b) > 0 for b in actual if len(b[0]) == 4 and b[0][-1] == 3)
    assert actual[0][0] == ('timeout',)
    assert F(int(actual[0][1].sum()), actual[0][2]) == F(25, 36)
    # Local native records remain temporally correlated; resetting the carrier
    # after each tick would replace this conditional terminal vector.
    _, maps, _, denominator = e.stopped_tables(2)
    mark = next(k for k in maps if k[1] == 1)
    first = np.ones(66, dtype=np.int64)@maps[mark]
    assert len(set(map(int, first))) > 1


def test_reference_entanglement_survives_unread_source_tick_and_error_covers_all_columns():
    gate = e._gates(2)[2].kraus[0][0]
    vector = np.zeros(8, complex); vector[0] = vector[5] = 1/math.sqrt(2)
    out = np.kron(gate, np.eye(2))@vector
    reduced = out.reshape(4, 2).conj().T@out.reshape(4, 2)
    np.testing.assert_allclose(reduced, np.eye(2)/2, atol=5e-15, rtol=0)
    assert c.certify(gate, 'phase', s=1) <= c.ERROR_GATE
    assert c._sequence_allowance(8) < e.COMPLETE_ALLOWANCE
    damaged = gate.copy(); damaged[2, 1] += 1e-4
    assert c.certify(damaged, 'phase', s=1) > c.ERROR_GATE


def test_no_reset_and_native_current_intervention_have_distinct_complete_laws():
    ordinary = e.complete_output(2, feedback=False, ideal=True)
    reset = e.complete_output(2, feedback=False, ideal=True, reset=True)
    feedback = e.complete_output(2, feedback=True, ideal=True)
    select = lambda b: len(b[0]) == 4 and b[0][-1] == 0
    p = lambda rows: sum(e._probability(b) for b in rows if select(b))
    assert p(reset)-p(ordinary) > .001
    assert abs(p(feedback)-p(ordinary)) > .0001
    # This remains an unconditional complete-event contrast: timeouts stay in
    # the output and are not divided out to inflate the result.
    assert sum(e._probability(b) for b in ordinary if len(b[0]) == 4) == pytest.approx(F(121, 1296))


def test_passive_waiting_and_movement_only_triggering_miss_zero_current_modular_advance():
    _, first, _ = e.probe_outputs(2, current=0, feedback=False, ideal=True)
    probability = float(np.vdot(first[2], first[2]).real)
    assert probability > .01
    # With no U(1), the balanced preparation/analysis gives outcome0 exactly.
    b = np.eye(4); b[0, 0] = b[0, 2] = b[2, 0] = math.sqrt(.5); b[2, 2] = -math.sqrt(.5)
    passive = b@b@np.array([1., 0, 0, 0])
    assert abs(passive[2]) < 1e-15


def test_conjugation_proof_has_nonunit_norm_and_distinct_positive_roots():
    x, y = (25+sp.sqrt(3))/30, (25-sp.sqrt(3))/30
    assert 0 < y < x < 1
    assert sp.simplify(x*y) == sp.Rational(311, 450)
    assert sp.simplify(x.subs(sp.sqrt(3), -sp.sqrt(3))-y) == 0
    proof = e.multiplicative_independence()
    assert proof['integer_relation_lambda_plus_power_m_times_lambda_minus_power_n_equals_one'] == 'm=n=0'
    # A rational-frequency control violates the claimed joint injectivity.
    assert abs(np.exp(1j*2*np.pi)-1) < 1e-14
    assert abs(np.exp(1j*2*2*np.pi)-1) < 1e-14


@pytest.mark.parametrize('ticks', (0, 1, 64, 272, 4096))
def test_outward_phase_power_intervals_against_independent_high_precision(ticks):
    with localcontext() as context:
        context.prec = 85
        for mode, root in enumerate(((25+Decimal(3).sqrt())/30, (25-Decimal(3).sqrt())/30), 1):
            # Decimal Taylor after a high-precision 2pi reduction is only an
            # independent diagnostic; the certificate uses dyadic interval powering.
            angle = -ticks*root.ln()
            pi = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628')
            angle %= 2*pi
            cos = Decimal(1)+sum((-1)**k*angle**(2*k)/Decimal(math.factorial(2*k)) for k in range(1, 100))
            sin = sum((-1)**k*angle**(2*k+1)/Decimal(math.factorial(2*k+1)) for k in range(100))
            for actual, (lo, hi) in zip((cos, sin), e.phase_power(ticks, mode), strict=True):
                assert Decimal(lo.numerator)/Decimal(lo.denominator) <= actual
                assert actual <= Decimal(hi.numerator)/Decimal(hi.denominator)


def test_joint_phase_window_is_finite_and_has_a_certified_later_near_recurrence():
    audit = e.resolution_audit()
    assert audit['tick_window'] == [0, 64]
    assert audit['certified_approximate_recurrence_ticks'] == 272
    assert F(audit['approximate_recurrence_max_chord_upper']) < F(1, 8)
    assert not audit['unbounded_time_perfect_resolution_claimed']
    assert not e.multiplicative_independence()['single_measurement_determines_elapsed_parameter']


@pytest.mark.parametrize('bad', (True, -1, 0, 3, 2., '2'))
def test_invalid_attempt_budget_is_not_silently_a_clock_phase(bad):
    with pytest.raises(ValueError):
        e.stopped_tables(bad)


def test_invalid_phase_power_window_rejected():
    for bad in (True, -1, 2**20+1, 1.):
        with pytest.raises(ValueError):
            e.phase_power(bad, 1)


def test_reset_control_retains_reset_state_even_when_following_interval_times_out():
    initial, _, _ = e.probe_outputs(2, ideal=True)
    controlled = e.complete_output(2, feedback=False, ideal=True, reset=True)
    for record, _, _, vector in controlled:
        if len(record) == 3:
            weight = float(np.vdot(vector, vector).real)
            np.testing.assert_allclose(np.outer(vector, vector.conj()),
                weight*np.outer(initial, initial.conj()), atol=1e-15, rtol=0)


def test_runtime_source_closure_is_retained_in_a_clean_process():
    script = '''
import sys
from pathlib import Path
from model.transfer_modular_probe import event_clock as e
e.run()
missing = []
for module in list(sys.modules.values()):
    path = getattr(module, '__file__', None)
    if path is None:
        continue
    try:
        relative = Path(path).resolve().relative_to(e.ROOT).as_posix()
    except ValueError:
        continue
    if relative.endswith('.py') and relative not in e.SOURCE_PATHS:
        missing.append(relative)
assert not missing, sorted(set(missing))
'''
    subprocess.run([sys.executable, '-c', script], cwd=e.ROOT, check=True, timeout=120)

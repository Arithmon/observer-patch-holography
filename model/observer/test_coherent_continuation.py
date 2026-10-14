"""Independent finite checks and hostile controls for the N1 continuation audit."""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
import json
import subprocess
import sys

import pytest

from . import coherent_continuation as c


@pytest.fixture(autouse=True)
def fresh_caches():
    for function in (c.branches, c.reversal, c._kernel_counts):
        function.cache_clear()
    yield
    for function in (c.branches, c.reversal, c._kernel_counts):
        function.cache_clear()


def endpoint(state, seam, coin):
    """Independent floor/ceiling rule; never calls native.native."""
    a, b = seam
    total = state[a] + state[b]
    low, high = total // 2, (total + 1) // 2
    out = list(state)
    out[a], out[b] = (high, low) if coin else (low, high)
    return tuple(out), state[a] - out[a]


def basis(index):
    return tuple(int(i == index) for i in range(len(c.STATES)))


def inner(left, right):
    return sum((a*b for a, b in zip(left, right, strict=True)), F(0)) / c.EVENTS


def test_every_event_reversal_and_record_dictionary_against_endpoint_oracle():
    expected_kernel = [[0] * len(c.STATES) for _ in c.STATES]
    for arc, (output, edge, coin, current) in enumerate(c.branches()):
        x = c.STATES[arc // c.EVENTS]
        y, j = endpoint(x, c.quantum_mixing.seams()[edge], coin)
        assert (c.STATES[output], current) == (y, j)
        reverse_coin = coin if current == 0 else 1-coin
        restored, reverse_j = endpoint(y, c.quantum_mixing.seams()[edge], reverse_coin)
        assert (restored, reverse_j) == (x, -current)
        reversed_arc = output * c.EVENTS + edge * 2 + reverse_coin
        assert c.reversal()[arc] == reversed_arc
        assert c.decode_arc(reversed_arc) == (output, edge, coin, current)
        assert c.encode_record(c.decode_arc(arc)) == arc
        expected_kernel[arc // c.EVENTS][output] += 1
    assert c.kernel_counts() == tuple(map(tuple, expected_kernel))
    assert sum(i == j for i, j in enumerate(c.reversal())) == 3360
    assert len(set(c.branches())) == c.DIMENSION == 3960


def test_actual_coherent_native_event_and_discriminant_for_all_columns():
    # Disjoint native supports prove orthogonality for every pair of columns.
    # A nontrivial rational superposition verifies amplitudes, including signs.
    f = tuple(F((i % 7) - 3, 11) for i in range(len(c.STATES)))
    actual = [F(0)] * c.DIMENSION
    for arc, record in enumerate(c.branches()):
        actual[c.encode_record(record)] = f[arc // c.EVENTS]
    assert tuple(actual) == c.advance(c.embed(f), -1) == c.advance(c.embed(f), 1)
    assert c.compress(actual) == c.kernel_apply(f)
    assert c.norm_squared(actual) == sum(x*x for x in f)
    for x in range(len(c.STATES)):
        # The compression of each column is just its independently counted
        # classical branch multiplicities, not a replacement amplitude rule.
        direct = [0] * len(c.STATES)
        for seam in c.quantum_mixing.seams():
            for coin in (0, 1):
                y, _ = endpoint(c.STATES[x], seam, coin)
                direct[c.INDEX[y]] += 1
        assert c.kernel_apply(basis(x)) == tuple(F(n, 60) for n in direct)


def test_all_load_columns_and_full_amplitude_recurrence():
    for port in range(12):
        f = tuple(F(6*x[port]-1) for x in c.STATES)
        closed = []
        for x in c.STATES:
            value = F(6*x[port]-1)
            for a, b in c.quantum_mixing.seams():
                if port == a:
                    value -= F(6*(x[a]-x[b]), 60)
                elif port == b:
                    value -= F(6*(x[b]-x[a]), 60)
            closed.append(value)
        assert c.kernel_apply(f) == tuple(closed)
    f = tuple(F((i % 5)-2) for i in range(len(c.STATES)))
    current = c.embed(f)
    previous = None
    expected = f
    for step in range(7):
        assert c.compress(current) == expected
        following = c.kernel_apply(expected)
        if previous is not None:
            following = tuple(2*a-b for a,b in zip(following, previous, strict=True))
        previous, expected = expected, following
        current = c.advance(current)


def test_same_first_native_operation_different_later_complete_read():
    initial = c.embed(basis(0))
    assert c.advance(initial, -1) == c.advance(initial, 1)
    assert c.advance(c.advance(initial, 1), 1) == initial
    returns = c.compress(c.advance(c.advance(initial)))
    assert returns[0] == F(38, 75)
    # Return probability is the square; Ramsey probabilities are affine in
    # the interference amplitude. Neither is the classical K^2 return.
    classical_return = c.kernel_apply(c.kernel_apply(basis(0)))[0]
    assert classical_return == F(113, 150)
    assert returns[0]**2 == F(1444, 5625) != classical_return != returns[0]
    assert tuple(c.norm_squared(v) for v in c.ramsey_outputs(basis(0), 2, -1)) == (F(113,150), F(37,150))
    assert tuple(c.norm_squared(v) for v in c.ramsey_outputs(basis(0), 2, 1)) == (1, 0)
    for phase in (-1, 1):
        assert tuple(c.norm_squared(v) for v in c.ramsey_outputs(basis(0), 1, phase)) == (F(14,15), F(1,15))


def test_complex_input_and_spectator_reference_gram_preserved():
    # Two possibly nonorthogonal reference columns, each with exact real and
    # imaginary parts; their complete Gram matrix determines the reduced
    # spectator state. No complex floating arithmetic is used.
    f0 = tuple(F(3,5)*x for x in basis(0)), tuple(F(4,5)*x for x in basis(65))
    f1 = tuple(F(4,5)*x for x in basis(0)), tuple(F(3,5)*x for x in basis(65))
    columns = (f0, f1)
    def complex_inner(a, b):
        ar, ai = a; br, bi = b
        return inner(ar, br)+inner(ai, bi), inner(ar, bi)-inner(ai, br)
    prepared = tuple(tuple(c.embed(v) for v in col) for col in columns)
    for phase in (-1, 1):
        outputs = tuple(tuple(c.ramsey_outputs(v, 3, phase) for v in col) for col in columns)
        for i in range(2):
            for j in range(2):
                terms = [complex_inner((outputs[i][0][r], outputs[i][1][r]),
                                       (outputs[j][0][r], outputs[j][1][r])) for r in (0,1)]
                assert tuple(sum(term[k] for term in terms) for k in (0,1)) == complex_inner(prepared[i], prepared[j])


def test_same_committed_instrument_different_discriminant():
    counts, twisted = c.kernel_counts(), c.kernel_counts(True)
    signs_by_record = {}
    signs_by_contrast_record = {}
    for arc, (output, edge, coin, current) in enumerate(c.branches()):
        sign = -1 if c.quantum_mixing.seams()[edge] == (0,1) and current == 1 else 1
        signs_by_record.setdefault((edge,coin,current), set()).add(sign)
        a, b = c.quantum_mixing.seams()[edge]
        x = c.STATES[arc // c.EVENTS]
        signs_by_contrast_record.setdefault((edge,coin,x[a]-x[b]), set()).add(sign)
    assert all(len(signs) == 1 for signs in signs_by_record.values())
    assert all(len(signs) == 1 for signs in signs_by_contrast_record.values())
    # Constant branch signs cancel from every L_r rho L_r^*, even for a
    # reference-entangled input; the coherent arc overlaps retain their sign.
    changed = [(i,j) for i in range(66) for j in range(i+1,66) if counts[i][j] != twisted[i][j]]
    assert len(changed) == 10
    assert all(twisted[i][j] == -counts[i][j] == -1 for i,j in changed)
    # Independent sparse triangle enumeration for trace D^3 - trace K^3.
    difference = sum(twisted[i][j]*twisted[j][k]*twisted[k][i] - counts[i][j]*counts[j][k]*counts[k][i]
                     for i in range(66) for j in range(66) if counts[i][j]
                     for k in range(66) if counts[j][k] and counts[k][i])
    assert F(difference, 60**3) == -F(1,900)


def test_full_78_carrier_has_one_way_transition():
    x = (2,) + (0,)*11
    y, current = endpoint(x, (0,1), 0)
    assert current == 1 and y == (1,1)+(0,)*10
    assert all(endpoint(y,seam,coin)[0] != x for seam in c.quantum_mixing.seams() for coin in (0,1))


@pytest.mark.parametrize("bad", [(), (0,)*65, (False,)+(0,)*65, (0.5,)+(0,)*65, (1j,)+(0,)*65])
def test_inexact_or_wrong_dimension_input_rejected(bad):
    with pytest.raises(ValueError): c.embed(bad)


@pytest.mark.parametrize("phase", [0, 2, True, F(-1), -1.0, 1j])
def test_bad_phase_rejected(phase):
    with pytest.raises(ValueError): c.advance(c.embed(basis(0)), phase)
    with pytest.raises(ValueError): c.ramsey_outputs(basis(0), 0, phase)


@pytest.mark.parametrize("steps", [-1, True, F(1), 1.0])
def test_bad_steps_rejected(steps):
    with pytest.raises(ValueError): c.ramsey_outputs(basis(0), steps)


@pytest.mark.parametrize("record", [(0,0,0,2), (0,0,True,0), (66,0,0,0), (0,30,0,0), (0,0,0), [0,0,0,0]])
def test_bad_record_rejected(record):
    with pytest.raises(ValueError): c.encode_record(record)


def test_record_sign_mutant_rejected(monkeypatch):
    original = c.native.native
    def wrong(state, seam, coin):
        output, (coin, current) = original(state,seam,coin)
        return output, (coin,-current)
    monkeypatch.setattr(c.native, "native", wrong)
    with pytest.raises(ArithmeticError, match="current"): c.branches()


def test_noninvolutive_shift_mutant_rejected(monkeypatch):
    original = c.encode_record
    monkeypatch.setattr(c, "encode_record", lambda record: (original(record)+1) % c.DIMENSION)
    with pytest.raises(ArithmeticError, match="involution"): c.reversal()


def test_missing_reflection_mutant_rejected(monkeypatch):
    original = c.advance
    monkeypatch.setattr(c, "advance", lambda values, phase=-1: original(values, 1))
    with pytest.raises(ArithmeticError, match="recurrence"): c.audit()


def test_receipt_strict_types_scope_and_custody(tmp_path):
    result = c.run()
    assert result["audit"]["all_pass"]
    c.compare_receipt(deepcopy(result), result)
    for mutate in (
        lambda forged: forged["audit"].__setitem__("all_pass", 1),
        lambda forged: forged["scope"].__setitem__("source_pulses_emitted", True),
        lambda forged: forged["source_sha256"].__setitem__(c.SOURCE_PATHS[0], "0"*64),
    ):
        forged = deepcopy(result); mutate(forged)
        with pytest.raises(ValueError, match="mismatch"): c.compare_receipt(forged, result)
    for content in ('{"x":1,"x":2}', '{"x":NaN}'):
        path = tmp_path/"invalid.json"; path.write_text(content)
        with pytest.raises(ValueError): c.read_receipt(path)


def test_cli_exclusive_write_and_optimized_hostile_check(tmp_path):
    path = tmp_path/"receipt.json"
    command = [sys.executable, "-m", "model.observer.coherent_continuation"]
    first = subprocess.run(command + ["--write", str(path)], cwd=c.ROOT, capture_output=True, text=True)
    assert first.returncode == 0, first.stderr
    repeat = subprocess.run(command + ["--write", str(path)], cwd=c.ROOT, capture_output=True, text=True)
    assert repeat.returncode != 0
    forged = json.loads(path.read_text()); forged["audit"]["all_pass"] = 1
    path.write_text(c.canonical(forged))
    checked = subprocess.run([sys.executable, "-O", *command[1:], "--check", str(path)],
                             cwd=c.ROOT, capture_output=True, text=True)
    assert checked.returncode != 0 and "receipt mismatch" in checked.stderr


@pytest.mark.parametrize("bad", [0, 1, F(1), None, "true"])
def test_kernel_switch_rejects_inexact_types_even_after_cache_warmup(bad):
    c.kernel_counts(True)
    c.kernel_counts(False)
    with pytest.raises(ValueError):
        c.kernel_counts(bad)


def test_structural_selection_audit_keeps_the_two_control_scopes_distinct():
    result = c.selection_audit()
    endpoints = result["phase_endpoint_criteria"]
    assert endpoints["coins"] == [{"phase": 1, "diagonal": "1", "off_diagonal": "0"},
                                  {"phase": -1, "diagonal": "-29/30", "off_diagonal": "1/30"}]
    assert endpoints["uniform_coherent_state_fixed"]
    assert result["proper_covariance"]["arc_checks_per_reversal"] == 237600
    assert not result["conditional_coin_classification"]["full_event_symmetry_is_a_source_premise"]
    signed = result["moving_record_sign_control"]
    assert signed["each_committed_C_and_E_CP_branch_unchanged"]
    assert not signed["same_first_coherent_event"]
    assert not signed["uniform_coherent_state_fixed"]
    assert signed["complete_uniform_return_probabilities"] == ["529/1089", "560/1089"]
    assert signed["complete_same_Ramsey_probabilities"] == ["28/33", "5/33"]
    assert signed["trace_D_cubed_minus_trace_K_cubed"] == "-1/90"
    assert not result["A1_A2_A3_source_admission_proved"]


def test_full_coin_classification_and_common_reversal_on_complex_columns():
    # Simultaneous permutation covariance makes every diagonal entry equal
    # and every off-diagonal entry equal. The two orthogonal solutions fix
    # the uniform vector and are checked against the actual runtime coin.
    reverse = lambda v: tuple(v[i] for i in c.reversal())
    for phase in (-1, 1):
        diagonal = F(1+59*phase,60); off = F(1-phase,60)
        assert (diagonal-off)**2 == 1 and diagonal+59*off == 1
        for event in (0, 37, 59):
            v = tuple(F(i == event) for i in range(c.DIMENSION))
            actual = reverse(c.advance(v,phase))
            expected = tuple(diagonal if i == event else off if i < 60 else F(0)
                             for i in range(c.DIMENSION))
            assert actual == expected
        # For complex v=u+i w, Theta=S conjugation acts on both parts with
        # opposite signs; the same S works for both continuation rules.
        for v in (tuple(F(i%3-1,7) for i in range(c.DIMENSION)),
                  tuple(F(i%5-2,11) for i in range(c.DIMENSION))):
            assert c.advance(reverse(c.advance(v,phase)),phase) == reverse(v)
            before = c.compress(v)
            after = reverse(c.advance(v,phase))
            assert c.compress(after) == before
            assert sum((x-before[i//60])**2 for i,x in enumerate(v)) == sum((x-before[i//60])**2 for i,x in enumerate(after))


def test_all_moving_arc_signs_by_independent_endpoint_oracle():
    signs = []
    signed = [[0]*66 for _ in range(66)]
    branches_by_label = ({}, {})
    for i, state in enumerate(c.STATES):
        for edge, seam in enumerate(c.quantum_mixing.seams()):
            for coin in (0,1):
                output, current = endpoint(state,seam,coin)
                sign = 1 if output == state else -1
                signs.append(sign)
                signed[i][c.INDEX[output]] += sign
                for groups, label in zip(branches_by_label,
                                         (current,state[seam[0]]-state[seam[1]]),strict=True):
                    groups.setdefault((edge,coin,label), []).append((i,sign))
    assert signs.count(-1) == 600
    for groups in branches_by_label:
        # Every matrix unit surviving each CP branch gets the same phase
        # on its ket and bra; reference matrix units factor identically.
        assert all(a*b == 1 for columns in groups.values() for _,a in columns for _,b in columns)
    counts = c.kernel_counts()
    assert all(signed[i][j] == (counts[i][j] if i == j else -counts[i][j])
               for i in range(66) for j in range(66))
    return_amplitude = F(sum(signs),c.DIMENSION)
    assert return_amplitude == F(23,33)
    assert return_amplitude**2 == F(529,1089)
    assert c.norm_squared(tuple(F(1+s,2) for s in signs))/66 == F(28,33)
    assert c.norm_squared(tuple(F(1-s,2) for s in signs))/66 == F(5,33)


def test_actual_record_adaptive_two_event_histories_hide_scalar_phases():
    # A supplied current-controlled second seam is enough to check that
    # commitment plus adaptive feedback does not recover lost branch phases.
    # Equality of the CP map for each step gives the all-history extension.
    for candidate in ("C","E"):
        phases = {}
        for initial in c.STATES:
            for first_coin in (0,1):
                first, j1 = endpoint(initial,(0,1),first_coin)
                first_label = j1 if candidate == "C" else initial[0]-initial[1]
                seam = (0,5) if j1 > 0 else (0,1)
                for second_coin in (0,1):
                    _, j2 = endpoint(first,seam,second_coin)
                    second_label = j2 if candidate == "C" else first[seam[0]]-first[seam[1]]
                    label = (first_coin,first_label,seam,second_coin,second_label)
                    sign = (-1 if j1 else 1)*(-1 if j2 else 1)
                    phases.setdefault(label,set()).add(sign)
        assert all(len(signs) == 1 for signs in phases.values())


def test_duplicate_proper_rotations_rejected(monkeypatch):
    from ..common_source.slow_geometry import source
    monkeypatch.setattr(source(), "rotations", lambda: (tuple(range(12)),)*60)
    with pytest.raises(ArithmeticError, match="distinct proper rotations"):
        c.selection_audit()


def test_selection_audit_rejects_missing_reflection_on_nonuniform_inputs(monkeypatch):
    original = c.advance
    monkeypatch.setattr(c, "advance", lambda values, phase=-1: original(values,1))
    with pytest.raises(ArithmeticError, match="selection criteria"):
        c.selection_audit()

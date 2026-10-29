"""Independent endpoint, absorption, quantum-reference and hostile checks."""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import experiment, transport as t


def _endpoint(pair, seam, coin, law):
    """Separate occupation arithmetic and raw-address encoding."""
    a,b = seam
    x = [int(i in pair) for i in range(12)]
    low, odd = divmod(x[a]+x[b], 2)
    left = low+coin*odd
    current = x[a]-left
    raw = {0:0, 1:1, -1:2}[current]+(4*odd if law == "E" else 0)
    x[a],x[b] = left,low+(1-coin)*odd
    return tuple(x), raw, current


def test_all_native_flags_and_logical_independence_from_separate_endpoint_arithmetic():
    for center in t.LIVE:
        for law in ("C", "E"):
            for seam in t.SEAMS:
                for coin in (0,1):
                    seen = [Counter(),Counter()]
                    for logical,column in enumerate(t.supports(center)):
                        for pair in column:
                            output,raw,current = _endpoint(pair,seam,coin,law)
                            seen[logical][raw] += 1
                            images,target,abort = t.branch(center,law,seam,coin,raw)
                            assert t.circuit.encode(output) in images[logical]
                            assert not abort
                            moved = tuple(sorted(seam[1] if k == seam[0] else seam[0] if k == seam[1] else k for k in center)) if current else center
                            assert target == moved
                    assert seen[0] == seen[1]
                    for raw in range(1 << t.circuit.RECORD_BITS[law]):
                        images,_,_ = t.branch(center,law,seam,coin,raw)
                        assert len(images[0]) == seen[0][raw]


def _geometric_kernel(weights):
    """Independent core jump rates: every core port has occupation one-half.

    All boundary neighbors have zero code occupation. Thus current through
    a core-to-exterior seam has rate theta_b/[60(theta_a+theta_b)].
    """
    result = []
    for center in t.LIVE:
        row = Counter()
        for a in center:
            for edge in t.SEAMS:
                if a not in edge:
                    continue
                b = edge[0] if edge[1] == a else edge[1]
                if b in center:
                    continue
                target = tuple(sorted((set(center)-{a})|{b}))
                row[target] += F(weights[b],60*(weights[a]+weights[b]))
        row[center] = 1-sum(row.values())
        result.append(tuple(row[c] for c in t.CORES))
    return tuple(result)


@pytest.mark.parametrize("p", (t.parameters(), t.local.Parameters()))
def test_complete_killed_kernel_matches_geometric_jump_rates(p):
    expected = _geometric_kernel(p.port_weights)
    assert t.kernel(p,"C") == expected == t.kernel(p,"E")
    assert all(sum(row) == 1 for row in expected)
    # Renormalizing surviving rows would erase a positive stopping event.
    i = t.LIVE.index((0,2))
    assert sum(expected[i][t.CORES.index(c)] for c in t.LIVE) < 1


def test_exact_first_exit_equations_and_uniform_class_lumping():
    result = t.stopping_certificate()
    direct = _geometric_kernel(t.parameters().port_weights)
    times = list(map(F,result["expected_stopping_attempts"]))
    for i in range(60):
        assert times[i] == 1+sum(direct[i][t.CORES.index(c)]*times[j] for j,c in enumerate(t.LIVE))
    assert F(result["minimum_two_attempt_stop_probability"]) == F(49,245025)
    assert F(result["initial_movement_probability"]) == F(16,495)
    uniform = _geometric_kernel((1,)*12)
    for i,center in enumerate(t.LIVE):
        classes = Counter()
        for target,probability in zip(t.CORES,uniform[i],strict=True):
            classes[len(t.screening.safe_exterior(target))] += probability
        assert classes == ({4:F(29,30),2:F(1,30),0:F()} if center in t.SEAMS
                           else {4:F(1,30),2:F(19,20),0:F(1,60)})
        expected_time = F(150 if center in t.SEAMS else 120)
        assert expected_time == 1+sum(uniform[i][t.CORES.index(c)]*(150 if c in t.SEAMS else 120) for c in t.LIVE)


def test_holds_slow_exit_and_do_not_change_escape_class():
    p = t.parameters()
    slow = t.replace(p,repair_probability=F(1,2))
    base,changed = t.kernel(p),t.kernel(slow)
    for i,center in enumerate(t.LIVE):
        assert changed[i] == tuple(v/2+F(target==center,2) for target,v in zip(t.CORES,base[i],strict=True))


def test_same_reference_exterior_witness_and_modular_leakage():
    reference = t.local.reference(t.parameters())
    states = t.local.basis()
    exterior = t.exterior_certificate()["rows"]
    for row in exterior:
        center = tuple(row["core"])
        # Independent full exterior marginal of pi and its q1 condition.
        out = [i for i in range(12) if i not in center]
        sigma = Counter()
        conditional = Counter()
        mass = F()
        for x,w in zip(states,reference,strict=True):
            marginal = tuple(x[i] for i in out)
            sigma[marginal] += w
            if sum(marginal) == 1:
                conditional[marginal] += w
                mass += w
        distance = sum(abs(v-conditional[key]/mass) for key,v in sigma.items())/2
        assert distance == F(row["sharp_q1_full_exterior_distance_minimum"]) == 1-mass
    assert exterior[0]["sharp_q1_full_exterior_distance_minimum"] == "701/1501"
    for row in t.clock_compatibility()["rows"]:
        center = tuple(row["core"])
        v = t.ideal_code(center)
        evolved = v.copy()
        for x in states:
            theta = math.prod(w**n for w,n in zip(t.parameters().port_weights,x,strict=True))
            evolved[t.circuit.encode(x)] *= np.exp(-.37j*math.log(theta))
        loss = np.eye(2)-evolved.conj().T@v@v.conj().T@evolved
        expected = float(F(row["modular_code_leakage_sine_squared_coefficient"]))*math.sin(.37*math.log(25/8)/2)**2
        assert np.linalg.norm(loss-expected*np.eye(2)) < 1e-14
    assert t.clock_compatibility()["invariant_core_count"] == 19


@pytest.fixture(scope="module")
def diagnostics():
    return t.source_replay_diagnostics()


def test_source_compiled_motion_retains_complete_complex_reference_and_stop(diagnostics):
    assert all(diagnostics[k] < 1e-8 for k in diagnostics if k != "changing_core_C_E_path")
    path = diagnostics["changing_core_C_E_path"]
    assert [r["status"] for r in path] == ["continue","opposite_core_stop"]
    assert path[0]["core_before"] != path[0]["core_after"] == path[1]["core_before"]
    assert path[0]["core_after"] != list(t.INITIAL)
    assert path[1]["cumulative_flag_weight_given_selected_seams"] == "16/1089"
    assert [r["retained_raw_slots"] for r in path] == [8,16]
    experiment.admit(diagnostics)


def test_stop_branch_is_identity_service_not_hidden_reset_or_renormalization():
    center,law,seam,coin,raw = t.witness_path()[1]
    state = np.zeros((512,2),complex)
    state[7,0],state[23,1] = .3j,.4
    _,target,_ = t.branch(center,law,seam,coin,raw)
    assert target in t.STOP
    assert np.array_equal(t.recover(center,law,seam,coin,raw,state),state)
    assert all(v == 0 for v in t.recovery_resources(center,law,seam,coin,raw).values())


def test_layout_corruption_is_detected(monkeypatch):
    original = t.routing.permutation_edges
    def damaged(images,**kwargs):
        edges = original(images,**kwargs)
        return edges[1:] if edges else edges
    monkeypatch.setattr(t.routing,"permutation_edges",damaged)
    with pytest.raises(ArithmeticError,match="Gray schedule"):
        t.compiled_layout_certificate()


def test_complete_source_work_includes_stopping_attempt_and_encoding():
    result = t.resource_certificate()
    costs = list(map(F,result["equal_C_E_expected_source_pulses_until_and_including_stop"]))
    independent = _geometric_kernel(t.parameters().port_weights)
    for i,row in enumerate(result["rows"]):
        mean = sum(F(row["laws"][law]["uniform_seam_expected_source_pulses_per_attempt"]) for law in ("C","E"))/2
        assert costs[i] == mean+sum(independent[i][t.CORES.index(c)]*costs[j] for j,c in enumerate(t.LIVE))
    assert F(result["initial_expected_source_pulses_including_encoder_and_stopping_attempt"]) == result["initial_encoder"]["native_pulses"]+costs[t.LIVE.index(t.INITIAL)]
    assert result["maximum_source_pulses_per_attempt_all_live_cores_and_laws"] == 864887


def test_source_recovery_wrong_target_cannot_pass_replay(monkeypatch):
    original = t.recover
    def reset(center,law,seam,coin,raw,value):
        result = original(center,law,seam,coin,raw,value)
        images,target,_ = t.branch(center,law,seam,coin,raw)
        if not images[0] or target not in t.LIVE or target == t.INITIAL:
            return result
        # A genuine fixed-core-reset countercontrol after the correct map.
        decoded = t.record_code._apply_permutation(t.encoder_permutation(target),result,True)
        decoded = t.compiler.execute(t.fixed._mixing_word(len(t.supports(target)[0]),True),9,decoded)
        return t._encode_register(t.INITIAL,decoded)
    monkeypatch.setattr(t,"recover",reset)
    # All encoders still pass; this fails the complete branch comparison.
    with pytest.raises(ArithmeticError,match="fixed numerical"):
        experiment.admit(t.source_replay_diagnostics())


@pytest.mark.parametrize("bad", (float("nan"),float("inf"),-.1,1e-7,True))
def test_source_numeric_admission_is_fail_closed(diagnostics,bad):
    changed = deepcopy(diagnostics)
    changed["complete_path_max_branch_operator_error"] = bad
    with pytest.raises(ArithmeticError):
        experiment.admit(changed)


def test_receipt_rejects_boolean_integer_swap_and_fabricated_stop(diagnostics):
    admitted = experiment.admit(diagnostics)
    altered = deepcopy(admitted)
    altered["changing_core_C_E_path"][0]["coin"] = False
    with pytest.raises(ValueError):
        experiment.compare(admitted,altered)
    altered = deepcopy(diagnostics)
    altered["changing_core_C_E_path"][1]["status"] = "continue"
    with pytest.raises(ValueError):
        experiment.admit(altered)


@pytest.mark.parametrize("center", ((1,2),(1,0),(False,1),(0,0)))
def test_invalid_or_opposite_core_rejected(center):
    with pytest.raises(ValueError):
        t.supports(center)

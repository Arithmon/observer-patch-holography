"""Independent endpoint, exterior, exchange-fiber and source-interface checks."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import combinations, product
import json

import numpy as np
import pytest

from . import sector as s
from ..observer import native, record_recovery as records
from ..observer.source_discrimination import circuit
from ..native_modular_reference import history_transport
from ..common_source import encoded_memory as memory


def endpoint(x, edge, coin, law):
    """Independent integer endpoint oracle, including the E raw contrast."""
    a, b = edge
    y = list(x)
    y[a], y[b] = (x[a]+x[b]+coin)//2, (x[a]+x[b]+1-coin)//2
    value = x[a]-y[a] if law == "C" else x[a]-x[b]
    return tuple(y), (*edge, coin, value)


@pytest.mark.parametrize("law", ("C", "E"))
def test_record_lift_preserves_q_on_every_core_and_actual_endpoint(law):
    for edge in records.seams():
        for coin in (0, 1):
            for x in s.BASIS:
                y, rec = endpoint(x, edge, coin, law)
                assert records.recorded_step(law, x, edge, coin) == (y, rec)
                for c in s.CORES:
                    assert s.charge(x, c) == s.charge(y, s.moved_core(c, law, rec))


def test_sharp_sector_bounds_are_attained_on_actual_exterior_configurations():
    result = s.sector_certificate()
    assert [row["full_exterior_distance_minimum"] for row in result["sectors"]] == ["3/5", "4/9", "43/45"]
    assert [row["source_core_dimension"] for row in result["sectors"]] == [2970, 1320, 66]
    # Calculate reduced exterior probabilities, independently of total-number
    # grouping. Conditioned references attain the SAME full-exterior bound.
    for q in range(3):
        sigma, conditioned = {}, {}
        p = sum(w for x, w in zip(s.BASIS, s.reference()) if x[0]+x[1] == q)
        for x, w in zip(s.BASIS, s.reference()):
            outside = x[2:]
            sigma[outside] = sigma.get(outside, F()) + w
            conditioned[outside] = conditioned.get(outside, F()) + (w/p if x[0]+x[1] == q else 0)
        assert sum(abs(sigma[k]-conditioned[k]) for k in sigma)/2 == 1-p
    # The entire mixture is the reference: the sharp-sector assumption matters.
    assert sum(row["fixed_core_dimension"]*s.WEIGHTS[row["q"]] for row in result["sectors"]) == 1


def test_quantum_sector_superpositions_keep_the_number_witness_but_can_have_extra_distance():
    indices = [i for i, x in enumerate(s.BASIS) if x[0]+x[1] == 1]
    psi = np.zeros(66, complex)
    psi[indices[0]], psi[indices[-1]] = 1/np.sqrt(2), 1j/np.sqrt(2)
    rho = np.outer(psi, psi.conj())
    projection = np.diag([int(x[0]+x[1] == 1) for x in s.BASIS])
    assert np.trace(projection@rho).real == pytest.approx(1)
    sigma = np.diag(list(map(float, s.reference())))
    assert np.abs(np.linalg.eigvalsh(rho-sigma)).sum()/2 >= float(F(4, 9))


@pytest.mark.parametrize("law", ("C", "E"))
def test_complete_triangle_histories_and_independent_ordered_lifts_have_singleton_odd_fibers(law):
    # All native triangular faces, both orientations, and all 66 starting
    # states. This is a finite hostile witness suite; the general proof is
    # the common-record spectator invariant, not an exhaustive horizon claim.
    edges = set(records.seams())
    triangles = [t for t in combinations(range(12), 3)
                 if all(tuple(sorted(e)) in edges for e in combinations(t, 2))]
    assert len(triangles) == 20
    count = 0
    for a, b, c in triangles:
        for word in (((a, c), (a, b), (b, c)), ((b, c), (a, b), (a, c))):
            groups = {}
            for i, x in enumerate(s.BASIS):
                for coins in product((0, 1), repeat=3):
                    y, labels, tape = x, list(s.PAIRS[i]), []
                    for edge, coin in zip(word, coins):
                        z, record = endpoint(y, edge, coin, law)
                        changed = [p for p in range(12) if y[p] != z[p]]
                        if changed:
                            old = next(p for p in changed if y[p] == 1)
                            new = next(p for p in changed if z[p] == 1)
                            labels[labels.index(old)] = new
                        tape.append(record)
                        y = z
                    key = tuple(tape)
                    groups.setdefault(key, []).append((i, y == x and tuple(labels) == s.PAIRS[i][::-1]))
            branches = history_transport.word_instrument(law, word)
            assert set(groups) == {b.history for b in branches}
            for b in branches:
                expected = tuple(i for i, odd in groups[b.history] if odd)
                assert s.odd_closed_inputs(law, b) == expected
                if expected:
                    assert len(b.inputs) == 1
                    count += 1
    assert count == 40


def test_open_movement_can_have_large_fiber_and_erased_tape_is_not_the_same_channel():
    b = next(b for b in history_transport.word_instrument("C", ((0, 5),)) if b.history[0][-1] == 1)
    assert len(b.inputs) == 10 and s.odd_closed_inputs("C", b) == ()
    # Hiding the spectator-compatible inputs is rejected by exact native replay.
    forged = replace(b, inputs=b.inputs[:1], outputs=b.outputs[:1])
    with pytest.raises(ValueError, match="every actual compatible"):
        s.odd_closed_inputs("C", forged)


def test_rank_one_odd_flag_destroys_reference_entanglement_and_its_scalar_sign_is_unobservable():
    b = next(b for b in history_transport.word_instrument("C", s.EXCHANGE_WORD) if s.odd_closed_inputs("C", b))
    k = np.zeros((66, 66))
    k[b.outputs[0], b.inputs[0]] = np.sqrt(float(b.weight))
    other = (b.inputs[0]+1) % 66
    psi = np.zeros((66, 2), complex)
    psi[b.inputs[0], 0], psi[other, 1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    out = k@psi
    assert np.linalg.matrix_rank(out, tol=1e-12) == 1
    np.testing.assert_array_equal(np.outer(out.ravel(), out.ravel().conj()),
                                  np.outer((-out).ravel(), (-out).ravel().conj()))
    compression = psi.conj().T@k.conj().T@k@psi
    assert np.linalg.det(compression) == 0 and np.trace(compression).real > 0


def test_occupation_rank_gives_sharp_capacity_and_actual_all_seam_code_attains_it():
    result = s.code_capacity_certificate()
    assert result["optimal_exact_all_seam_correctable_dimension"] == 2
    assert set(result["exterior_occupation_ranks"].values()) == {2}
    # Occupation differences are the directed moving effects, on each seam.
    for a, b in records.seams():
        for x in s.BASIS:
            assert x[a]*(1-x[b])-x[b]*(1-x[a]) == x[a]-x[b]
    selected = [x for x in s.BASIS if x[0]+x[1] == 1]
    assert all(sum(x[2:]) == 1 for x in selected)


def test_source_executed_recovery_and_core_reset_control():
    result = s.execution_certificate()
    assert result["laws"]["C"]["after_recovery_without_core_reset_q1_probability"] == "7/8"
    assert result["laws"]["E"]["after_recovery_and_explicit_core_reset_q1_probability"] == "1"
    assert result["services"]["extra_core_reset_archive_bits_per_entry"] == 15
    assert result["services"]["fresh_quantum_bits_per_recovery"] == 0


def test_complete_terminal_read_retains_reference_and_padding_and_does_not_preserve_quantum_memory():
    a = np.zeros((512, 3), complex)
    a[0, 0], a[511, 1], a[400, 2] = .3j, 2., .7+.2j
    rows = s.terminal_rows(a)
    assert len(rows) == 512 and [i for i, _ in rows] == list(range(512))
    np.testing.assert_allclose(sum((r.conj().T@r for _, r in rows), np.zeros((3, 3))), a.conj().T@a)
    assert s.charge_read_masses(a)["outside_binary_shell"] == 4
    # A computational boundary read destroys the code's cross-row phases.
    plus = memory.ideal_code_matrix()@np.array([1., 1.])/np.sqrt(2)
    minus = memory.ideal_code_matrix()@np.array([1., -1.])/np.sqrt(2)
    np.testing.assert_array_equal(np.abs(plus)**2, np.abs(minus)**2)


@pytest.mark.parametrize("bad", ((False, 1), (1, 0), (0, 0), (0, 12), [0, 1], (0., 1)))
def test_core_validation_does_not_coerce_invalid_metadata(bad):
    with pytest.raises(ValueError):
        s.reference(bad)


@pytest.mark.parametrize("bad", (np.zeros(511), np.ones(512, bool), [float("nan")]*512))
def test_terminal_read_rejects_malformed_source_before_decoding(bad, monkeypatch):
    monkeypatch.setattr(memory, "decoding", lambda _: pytest.fail("invalid input reached decoder"))
    with pytest.raises(ValueError):
        s.terminal_rows(bad, setting="logical")


def test_type_sensitive_duplicate_safe_receipt_check(tmp_path):
    path = tmp_path/"receipt.json"
    for text in ('{"x":false}', '{"x":0,"x":0}'):
        path.write_text(text)
        with pytest.raises(ValueError):
            s.check_receipt(path, {"x": 0})
    path.write_text(json.dumps({"x": 0}))
    s.check_receipt(path, {"x": 0})


@pytest.mark.parametrize("law", ("C", "E"))
def test_same_code_activity_uses_all_native_attempts_and_is_reference_dependent(law):
    row = s.code_capacity_certificate()["complete_C_E_checks"][law]["uniform_seam_activity"]
    activity = []
    for x in s.BASIS:
        activity.append(F(sum(endpoint(x, edge, coin, law)[0] != x
                              for edge in records.seams() for coin in (0, 1)), 60))
    assert sum(activity)/66 == F(5, 33)
    code_support = {(0, 2), (1, 3), (0, 3), (1, 2)}
    assert {activity[i] for i, pair in enumerate(s.PAIRS) if pair in code_support} == {F(1, 6)}
    assert row["excess_over_uniform_source_per_attempt"] == "1/66"
    assert F(row["ideal_fixed_core_reference_background"]) == sum(a*w for a, w in zip(activity, s.reference()))
    assert F(row["ideal_fixed_core_reference_background"]) != F(5, 33)


def test_campaign_rejects_gate_promotion_and_missing_source_custody(tmp_path, monkeypatch):
    original = json.loads((s.HERE/"campaign.json").read_text())
    monkeypatch.setattr(s, "HERE", tmp_path)
    for key, value in (("integrated_gate_accepted", True), ("physical_comparison_count_added", False),
                       ("source_sha256", {})):
        changed = dict(original)
        changed[key] = value
        (tmp_path/"campaign.json").write_text(json.dumps(changed))
        with pytest.raises(ValueError):
            s.campaign_binding()

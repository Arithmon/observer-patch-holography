"""Independent Fock, full-record, arbitrary-reference and readout controls."""
from collections import Counter
from dataclasses import replace
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import statistics as s
from ..observer import native


def fock_hop(pair, destination, origin):
    """Independent bit-string annihilation then creation in the full Fock basis."""
    mask = sum(1 << p for p in pair)
    if not mask & (1 << origin):
        return None, 0
    phase = (-1)**((mask & ((1 << origin)-1)).bit_count())
    mask ^= 1 << origin
    if mask & (1 << destination):
        return None, 0
    phase *= (-1)**((mask & ((1 << destination)-1)).bit_count())
    mask |= 1 << destination
    return tuple(p for p in range(12) if mask & (1 << p)), phase


def apply(row, state):
    out = np.zeros_like(state)
    for i, j, phase in zip(row.inputs, row.outputs, row.signs, strict=True):
        out[j] = phase*state[i]
    return out


def half_distance(left, right):
    """Pure source/reference output distance without a large eigensolve."""
    left, right = left.reshape(-1, 1), right.reshape(-1, 1)
    # Nonzero eigenvalues of LL* - RR* are those of the compressed form.
    factors = np.column_stack((left, right))
    q, r = np.linalg.qr(factors)
    signs = np.diag([1]*left.shape[1]+[-1]*right.shape[1])
    return float(np.abs(np.linalg.eigvalsh(r@signs@r.conj().T)).sum()/2)


def test_all_bilinears_equal_independent_full_fock_restriction():
    for a, b in product(range(12), repeat=2):
        for pair, (j, phase) in zip(s.PAIRS, s.car_bilinear(a, b), strict=True):
            target, expected = fock_hop(pair, a, b)
            assert phase == expected
            assert (s.PAIRS[j] if phase else None) == target
    report = s.car_certificate()
    assert report["bilinear_commutator_columns"] == 12**4*66
    assert report["moving_CAR_columns"] == 1200
    assert report["standalone_creation_annihilation_CAR_on_H66"] is False


@pytest.mark.parametrize("law", ("C", "E"))
def test_complete_slots_endpoint_oracle_signs_and_effects(law):
    rows = s.branches(law)
    assert len(rows) == (240 if law == "C" else 480)
    assert len({(r.seam, r.coin, r.raw) for r in rows}) == len(rows)
    for row in rows:
        a, b = row.seam
        oracle = []
        for i, pair in enumerate(s.PAIRS):
            x = tuple(int(p in pair) for p in range(12))
            y = list(x)
            total = x[a]+x[b]
            y[a], y[b] = (total+row.coin)//2, (total+1-row.coin)//2
            current = x[a]-y[a]
            raw = {0: 0, 1: 1, -1: 2}[current]+(4*(total % 2) if law == "E" else 0)
            if raw != row.raw:
                continue
            target = tuple(p for p, n in enumerate(y) if n)
            if current:
                car_target, phase = fock_hop(pair, b if current == 1 else a, a if current == 1 else b)
                assert target == car_target
            else:
                phase = 1
            oracle.append((i, s.INDEX[target], phase))
        assert tuple(zip(row.inputs, row.outputs, row.signs)) == tuple(oracle)
        assert len(set(row.outputs)) == len(row.outputs)
        positive = s.branch(law, row.seam, row.coin, row.raw, 1)
        assert positive.inputs == row.inputs and positive.outputs == row.outputs
        assert set(positive.signs) <= {1}
    result = s.complete_certificate()[law]
    assert result["input_columns"] == 3960
    assert result["mixed_sign_slots"] == 42
    assert result["characters_checked"] == (1, -1)
    assert result["counts_are_per_character"] is True
    assert result["mixed_sign_count_character"] == -1
    assert result["positive_lift_matches_complete_native_instrument"] is True


def test_all_two_attempt_classical_feedback_histories_with_entangled_reference():
    rng = np.random.default_rng(1281)
    psi = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    psi /= np.linalg.norm(psi)
    total, different = 0.0, 0
    for coin, raw in product((0, 1), range(4)):
        first = ("C", (0, 5), coin, raw)
        # Classical feedback chooses the next actual seam from committed data.
        second_seam = (0, 11) if raw == 1 else (0, 1)
        for second_coin, second_raw in product((0, 1), range(8)):
            history = (first, ("E", second_seam, second_coin, second_raw))
            plus, minus = s.history_map(history, 1), s.history_map(history, -1)
            assert plus.inputs == minus.inputs and plus.outputs == minus.outputs
            assert plus.squared_amplitude == minus.squared_amplitude == F(1, 4)
            left, right = apply(plus, psi)/2, apply(minus, psi)/2
            # Every retained history and terminal population, including its
            # untouched reference state, agrees separately. No normalization.
            for i in range(66):
                assert np.allclose(np.outer(left[i], left[i].conj()),
                                   np.outer(right[i], right[i].conj()), atol=1e-16)
            assert np.allclose(left.conj().T@left, right.conj().T@right, atol=1e-16)
            total += float(np.linalg.norm(left)**2)
            different += half_distance(left, right) > 1e-8
    assert total == pytest.approx(1, abs=2e-15)
    assert different > 0
    identity = s.history_map(())
    assert identity.inputs == identity.outputs == tuple(range(66))
    assert identity.squared_amplitude == 1 and set(identity.signs) == {1}


@pytest.mark.parametrize("law", ("C", "E"))
def test_sharp_channel_distance_all_seams_and_reference_upper_bound(law):
    rng = np.random.default_rng(508)
    entangled = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    entangled /= np.linalg.norm(entangled)
    report = s.channel_distance_certificate()
    for entry in report["rows"]:
        seam = tuple(entry["seam"])
        expected = F(entry["laws"][law]["half_diamond_distance"])
        states = [entangled]
        if entry["attaining_input_pairs"]:
            state = np.zeros((66, 1), complex)
            for pair in entry["attaining_input_pairs"]:
                state[s.INDEX[tuple(pair)], 0] = 1/np.sqrt(2)
            states.append(state)
        for n, state in enumerate(states):
            distance = 0.0
            for coin, raw in product((0, 1), range(4 if law == "C" else 8)):
                plus = s.branch(law, seam, coin, raw, 1)
                minus = s.branch(law, seam, coin, raw, -1)
                distance += half_distance(apply(plus, state), apply(minus, state))/2
            assert distance <= float(expected)+2e-14
            if n == 1 or expected == 0:
                assert distance == pytest.approx(float(expected), abs=2e-14)
    assert Counter(row["laws"][law]["half_diamond_distance"] for row in report["rows"]) == {"1/2": 21, "0": 9}


@pytest.mark.parametrize("law", ("C", "E"))
def test_complete_witness_measurement_and_scalar_phase_dephasing_controls(law):
    indices = tuple(s.INDEX[p] for p in s.WITNESS_INPUT_PAIRS)
    psi = np.zeros((66, 1), complex)
    psi[list(indices), 0] = 1/np.sqrt(2)
    totals = {1: np.zeros(2), -1: np.zeros(2)}
    for coin, raw in product((0, 1), range(4 if law == "C" else 8)):
        row = s.branch(law, s.WITNESS_SEAM, coin, raw)
        positions = {i: j for i, j in zip(row.inputs, row.outputs, strict=True)}
        for character in (1, -1):
            actual = s.branch(law, s.WITNESS_SEAM, coin, raw, character)
            state = apply(actual, psi)/np.sqrt(2)
            components = [state[positions[i], 0] if i in positions else 0 for i in indices]
            read = np.array([components[0]+components[1], components[0]-components[1]])/np.sqrt(2)
            totals[character] += np.abs(read)**2
            # One scalar sign on the whole committed branch is invisible.
            assert np.allclose(np.outer(read, read.conj()), np.outer(-read, (-read).conj()))
            # Delete input coherence: both branch phase conventions give equal
            # +/- counts; total branch mass is still retained.
            for i in indices:
                singleton = np.zeros((66, 1), complex)
                singleton[i, 0] = 1/np.sqrt(2)
                out = apply(actual, singleton)/np.sqrt(2)
                entries = [out[positions[k], 0] if k in positions else 0 for k in indices]
                assert abs(entries[0]+entries[1])**2 == abs(entries[0]-entries[1])**2
    assert totals[1] == pytest.approx([1, 0], abs=1e-15)
    assert totals[-1] == pytest.approx([0.5, 0.5], abs=1e-15)
    report = s.witness()
    assert len([r for r in report["all_raw_outputs"] if r["law"] == law]) == (8 if law == "C" else 16)
    assert report["postselection_used"] is False


def test_chart_covariance_and_scalar_record_gauge_obstruction():
    result = s.covariance_certificate()
    assert result["proper_rotations"] == 60
    assert result["flat_local_cells"] == 495
    assert result["ordinary_group_product_columns"] == 60*60*66
    assert result["rotation_lift_is_projective_spin_half"] is False
    report = s.gauge_obstruction_certificate()
    one, exchange = report["witnesses"]
    assert [one["sign"], exchange["sign"]] == [1, -1]
    for law in ("C", "E"):
        assert Counter(one["directed_record_labels"][law]) == Counter(exchange["directed_record_labels"][law])
    # Independent arbitrary complex chart phases transport both state and read.
    rng = np.random.default_rng(610)
    chart = np.exp(1j*rng.normal(size=66))
    for row in s.branches("C"):
        values = [chart[j]*phase*chart[i].conjugate()
                  for i, j, phase in zip(row.inputs, row.outputs, row.signs, strict=True)]
        assert np.allclose(np.abs(values), 1)
    for loop in (one, exchange):
        phase = 1
        for left, right in zip(loop["pairs"], loop["pairs"][1:]):
            i, j = s.INDEX[tuple(left)], s.INDEX[tuple(right)]
            phase *= chart[j]*s.sign(tuple(left), tuple(right))*chart[i].conjugate()
        assert phase == pytest.approx(loop["sign"], abs=1e-14)


@pytest.mark.parametrize("bad", [False, True, 0, 1.0, -1.0, "-1", None])
def test_public_character_validation_is_strict_before_any_computation(bad):
    with pytest.raises(ValueError):
        s.branch("C", (0, 5), 0, 1, bad)
    with pytest.raises(ValueError):
        s.history_map((), bad)


@pytest.mark.parametrize("args", [
    ("C", (False, 5), 0, 1), ("C", (0.0, 5), 0, 1),
    ("C", (0, 5), False, 1), ("C", (0, 5), 0, True),
    ("C", (0, 5), 0, 1.0), ("C", [0, 5], 0, 1),
    ("C", (0, 2), 0, 1), ("E", (0, 5), 0, 8),
])
def test_public_complete_label_validation(args):
    s.branch("C", (0, 5), 0, 1)
    with pytest.raises(ValueError):
        s.branch(*args)


def test_census_rejects_changed_native_support_and_changed_signed_phase(monkeypatch):
    original = s.sign
    monkeypatch.setattr(s, "sign", lambda source, target, character=-1: 1)
    with pytest.raises(ArithmeticError, match="phase extension"):
        s.complete_certificate()
    monkeypatch.setattr(s, "sign", original)
    monkeypatch.setattr(native, "native", lambda state, seam, coin: (state, (coin, 0)))
    with pytest.raises(ArithmeticError, match="native branch support"):
        s.complete_certificate()


@pytest.mark.parametrize("damage", ("drop_zero", "duplicate_zero"))
def test_zero_and_abort_slots_cannot_be_dropped_or_duplicated(monkeypatch, damage):
    original = s.branches
    def changed(law, character=-1):
        rows = list(original(law, character))
        zero = next(i for i, row in enumerate(rows) if row.abort)
        if damage == "drop_zero":
            rows.pop(zero)
        else:
            another = next(i for i, row in enumerate(rows) if i != zero and row.abort)
            rows[another] = rows[zero]
        return tuple(rows)
    monkeypatch.setattr(s, "branches", changed)
    with pytest.raises(ArithmeticError, match="complete unique raw branch grid"):
        s.complete_certificate()


def test_full_certificate_rejects_positive_character_silently_using_negative_lift(monkeypatch):
    original = s.sign
    monkeypatch.setattr(s, "sign", lambda a, b, character=-1: original(a, b, -1))
    # All default negative-character checks still pass this mutation.  The
    # executable receipt must separately bind the claimed native + family.
    assert set(s.branch("C", (0, 5), 0, 1, 1).signs) == {-1, 1}
    with pytest.raises(ArithmeticError, match="positive lift differs"):
        s.certificate()


@pytest.mark.parametrize("damage", ("drop_zero", "wrong_support", "wrong_abort"))
def test_positive_lift_has_independent_complete_grid_support_and_metadata_checks(monkeypatch, damage):
    original = s.branches
    def changed(law, character=-1):
        rows = list(original(law, character))
        if character == 1:
            if damage == "drop_zero":
                rows.pop(next(i for i, row in enumerate(rows) if row.abort))
            elif damage == "wrong_support":
                index = next(i for i, row in enumerate(rows) if row.inputs)
                row = rows[index]
                rows[index] = replace(row, inputs=row.inputs[1:], outputs=row.outputs[1:], signs=row.signs[1:])
            else:
                index = next(i for i, row in enumerate(rows) if row.abort)
                rows[index] = replace(rows[index], abort=False)
        return tuple(rows)
    monkeypatch.setattr(s, "branches", changed)
    with pytest.raises(ArithmeticError, match="positive lift"):
        s.certificate()

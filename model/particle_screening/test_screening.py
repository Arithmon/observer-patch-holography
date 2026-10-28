"""Independent cap enumeration, source endpoint, and hostile code checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json

import numpy as np
import pytest

from . import screening as s
from ..exchange_sector import sector
from ..observer import record_recovery as records
from ..observer.source_discrimination import circuit


def _outside_number(pair, cap):
    return len(set(pair) & set(cap))


def test_all_cap_profiles_by_independent_basis_enumeration():
    """Enumerate actual cap outcomes, without invoking the hypergeom formula."""
    rows = s.cap_certificate()["rows"]
    sigma = sector.reference()
    for row in rows:
        k = row["cap_size"]
        for cap in combinations(range(2, 12), k):
            ref = Counter()
            for pair, weight in zip(sector.PAIRS, sigma, strict=True):
                ref[_outside_number(pair, cap)] += weight
            assert [str(ref[r]) for r in range(min(k, 2)+1)] == row["reference_exterior_number_distribution"]
            for q in range(3):
                pairs = [pair for pair in sector.PAIRS if len(set(pair) & set(s.CORE)) == q]
                conditional = Counter()
                for pair in pairs:
                    conditional[_outside_number(pair, cap)] += F(1, len(pairs))
                distance = sum(abs(conditional[r]-ref[r]) for r in range(min(k, 2)+1))/2
                assert str(distance) == row["sharp_q_minimax_trace_distances"][q]


def test_cap_minimax_is_not_a_bound_on_each_specified_cap():
    # A sharp q=1 basis state can be much closer to the reference on one
    # selected one-port cap than on another.
    pair = (0, 2)
    ref = F(61, 450)
    first = abs(F(int(2 in pair))-ref)
    second = abs(F(int(3 in pair))-ref)
    assert first == F(389, 450)
    assert second == ref
    assert s.cap_certificate()["rows"][0]["sharp_q_minimax_trace_distances"][1] == "8/225"


def test_coherent_sharp_state_can_have_more_cap_distance_than_number_witness():
    # Two source configurations with the same traced-out configuration retain
    # a genuinely quantum off-diagonal term on this two-port cap.
    cap = (2, 3)
    psi = {(0, 2): 1/np.sqrt(2), (0, 3): 1j/np.sqrt(2)}
    rho = np.zeros((4, 4), complex)
    sigma = np.zeros((4, 4), complex)

    def cap_address(pair):
        return sum(int(port in pair) << bit for bit, port in enumerate(cap))

    for left, amplitude in psi.items():
        for right, other in psi.items():
            if set(left)-set(cap) == set(right)-set(cap):
                rho[cap_address(left), cap_address(right)] += amplitude*other.conjugate()
    for pair, weight in zip(sector.PAIRS, sector.reference(), strict=True):
        sigma[cap_address(pair), cap_address(pair)] += float(weight)
    assert rho[1, 2] == pytest.approx(-.5j)
    trace_distance = np.abs(np.linalg.eigvalsh(rho-sigma)).sum()/2
    assert trace_distance > float(F(14, 225))


@pytest.mark.parametrize("center,expected", [((0, 1), 4), ((0, 2), 2), ((0, 3), 0)])
def test_safe_exterior_core_classes_and_one_port_occupation(center, expected):
    safe = s.safe_exterior(center)
    assert len(safe) == expected
    columns = s.sharp_code(center)
    if not safe:
        assert columns is None
        return
    assert columns is not None
    for column in columns:
        assert len(column) == expected
        assert all(len(set(pair) & set(center)) == 1 for pair in column)
        assert all(pair not in s.SEAMS for pair in column)
        frequencies = Counter(port for pair in column for port in pair if port not in center)
        assert set(frequencies) == set(safe)
        assert set(frequencies.values()) == {1}


def test_crossed_four_port_code_corrects_independent_endpoint_records_and_reference():
    columns = s.sharp_code(s.CORE)
    assert columns is not None
    assert columns == (
        ((0, 2), (1, 3), (0, 4), (1, 6)),
        ((0, 3), (1, 2), (0, 6), (1, 4)),
    )
    for law in ("C", "E"):
        for seam in records.seams():
            total = F()
            for coin in (0, 1):
                raw_of = {value: raw for raw, value in circuit.record_values(law, coin).items()}
                groups = [Counter(), Counter()]
                outputs = [set(), set()]
                for logical, column in enumerate(columns):
                    for pair in column:
                        x = tuple(int(port in pair) for port in range(12))
                        y, record = records.recorded_step(law, x, seam, coin)
                        raw = raw_of[record[-1]]
                        groups[logical][raw] += 1
                        outputs[logical].add((raw, y))
                assert groups[0] == groups[1]
                assert outputs[0].isdisjoint(outputs[1])
                total += sum(F(count, 8) for count in groups[0].values())
            assert total == 1
    # For any logical superposition V* n_j V is scalar, so the local
    # one-port state is input independent, including an untouched reference.
    v = np.zeros((66, 2), complex)
    for logical, column in enumerate(columns):
        for pair in column:
            v[sector.PAIRS.index(pair), logical] = .5
    np.testing.assert_array_equal(v.conj().T@v, np.eye(2))
    ref = F(61, 450)
    deviations = []
    for port in range(2, 12):
        effect = np.diag([int(port in pair) for pair in sector.PAIRS])
        compression = v.conj().T@effect@v
        expected = F(1, 4) if port in s.safe_exterior(s.CORE) else F()
        np.testing.assert_array_equal(compression, float(expected)*np.eye(2))
        deviations.append(abs(expected-ref))
    assert max(deviations) == F(61, 450)


def test_naive_ten_port_screening_breaks_a_real_native_flag():
    # Uniform exterior occupation would attain the unrestricted one-port
    # minimax, but core-exterior edge records expose this attempted code.
    outside = tuple(range(2, 12))
    groups = tuple(zip(outside[::2], outside[1::2], strict=True))
    columns = []
    for logical in (0, 1):
        column = []
        for u, v in groups:
            column.extend(((0, u if logical == 0 else v),
                           (1, v if logical == 0 else u)))
        columns.append(tuple(column))
    unequal = []
    for law in ("C", "E"):
        for seam in records.seams():
            for coin in (0, 1):
                groups_by_raw = [Counter(), Counter()]
                raw_of = {value: raw for raw, value in circuit.record_values(law, coin).items()}
                for logical, column in enumerate(columns):
                    for pair in column:
                        x = tuple(int(port in pair) for port in range(12))
                        _, record = records.recorded_step(law, x, seam, coin)
                        groups_by_raw[logical][raw_of[record[-1]]] += 1
                if groups_by_raw[0] != groups_by_raw[1]:
                    unequal.append((law, seam, coin))
    assert ("C", (0, 5), 0) in unequal
    assert ("E", (0, 5), 0) in unequal


@pytest.mark.parametrize("law", ("C", "E"))
def test_q2_edge_first_exit_matches_three_step_independent_dp(law):
    initial = s.CORE
    masses = {initial: F(1)}
    for t in range(1, 4):
        next_masses = Counter()
        for pair, mass in masses.items():
            x = tuple(int(port in pair) for port in range(12))
            for seam in records.seams():
                for coin in (0, 1):
                    y, _ = records.recorded_step(law, x, seam, coin)
                    out = tuple(port for port, n in enumerate(y) if n)
                    if out in s.SEAMS:
                        next_masses[out] += mass/F(60)
        masses = dict(next_masses)
        assert sum(masses.values()) == F(14, 15)**t
    assert s.transport_certificate()["uncontrolled_q2_edge_pair"][law]["edge_pair_mean_first_exit_attempts"] == "15"


def test_passive_error_horizon_is_strictly_after_previous_attempt_count():
    result = s.passive_memory_certificate()["rows"]
    for law, threshold in (("C", 23), ("E", 12)):
        assert result[law]["first_attempt_count_with_strictly_positive_final_decoder_error_lower"] == threshold
        before = s.adaptive.bounds(F(1), threshold-1, law_scope=law, binary=True,
                                   logical_dimension=2)
        after = s.adaptive.bounds(F(1), threshold, law_scope=law, binary=True,
                                  logical_dimension=2)
        assert before["encoded_final_recovery_trace_error_lower"] == 0
        assert after["encoded_final_recovery_trace_error_lower"] > 0


def test_receipt_rejects_type_swap_and_duplicate_keys(tmp_path):
    path = tmp_path/"receipt.json"
    for malformed in ('{"x":false}', '{"x":0,"x":0}'):
        path.write_text(malformed)
        with pytest.raises(ValueError):
            s.check_receipt(path, {"x": 0})
    path.write_text(json.dumps({"x": 0}))
    s.check_receipt(path, {"x": 0})


def test_campaign_pins_predecessors_and_science_hash(tmp_path, monkeypatch):
    original = json.loads((s.HERE/"campaign.json").read_text())
    monkeypatch.setattr(s, "HERE", tmp_path)
    for key, value in (("integrated_gate_accepted", True),
                       ("physical_comparison_count_added", False),
                       ("science_result_sha256", "0"*64)):
        candidate = dict(original)
        candidate[key] = value
        (tmp_path/"campaign.json").write_text(json.dumps(candidate))
        with pytest.raises(ValueError):
            s.campaign_binding(s.science_result())

"""Independent finite-state checks for the native defect decision experiment.

The oracle balances the two endpoint loads directly. It does not call the
producer's transition enumerator or the native update to obtain its counts.
This tests classical configuration invariants on the declared finite domain,
not all quantum, record-encoded or enlarged-source sectors.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from model.observer import native, native_loop_readout
from model.common_source import native_defect_decision as decision, record_code


@pytest.mark.parametrize("option", ["--verify", "--check"])
def test_receipt_cli_accepts_exact_evidence_and_rejects_physical_promotion(tmp_path, option):
    receipt = tmp_path/"receipt.json"
    packet = decision.report()
    receipt.write_text(json.dumps(packet))
    command = [sys.executable, "-m", "model.common_source.native_defect_decision", option, str(receipt)]
    root = Path(decision.__file__).resolve().parents[2]
    accepted = subprocess.run(command, cwd=root, capture_output=True, text=True)
    assert accepted.returncode == 0, accepted.stderr
    packet["scope"]["physical_energy_or_mass_derived"] = True
    damaged = json.dumps(packet)
    receipt.write_text(damaged)
    rejected = subprocess.run(command, cwd=root, capture_output=True, text=True)
    assert rejected.returncode != 0
    assert "native defect decision or source custody changed" in rejected.stderr
    assert receipt.read_text() == damaged


def _balance(loads, seam, coin):
    a, b = seam
    low, parity = divmod(loads[a] + loads[b], 2)
    result = list(loads)
    result[a] = low + coin * parity
    result[b] = low + (1 - coin) * parity
    return tuple(result)


def _counts(states, seams):
    indices = {state: i for i, state in enumerate(states)}
    rows = [[0] * len(states) for _ in states]
    for i, loads in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                rows[i][indices[_balance(loads, seam, coin)]] += 1
    return tuple(map(tuple, rows))


def _reachable(counts, start):
    reached = {start}
    pending = [start]
    for i in pending:
        for j, count in enumerate(counts[i]):
            if count and j not in reached:
                reached.add(j)
                pending.append(j)
    return reached


def test_all78_native_branch_counts_match_independent_endpoint_balancing():
    states = native.states()
    seams, _ = native_loop_readout.carrier()
    actual = decision.transitions()
    assert actual == _counts(states, seams)
    assert len(actual) == 78
    assert all(sum(row) == 60 for row in actual)
    assert sum(map(sum, actual)) == 4680


def test_binary_transfer_matches_existing_complete_native_enumerator():
    seams, _ = native_loop_readout.carrier()
    counts = decision.transitions(binary=True)
    existing, _ = native_loop_readout.binary_transfer(seams)
    assert counts == tuple(map(tuple, existing))
    assert counts == _counts(native.states(binary=True), seams)
    assert all(counts[i][j] == counts[j][i] for i in range(66) for j in range(66))
    assert all(counts[i][i] > 0 for i in range(66))


def test_complete_ce_record_maps_retain_every_zero_and_nonzero_outcome():
    states = native.states()
    indices = {state: i for i, state in enumerate(states)}
    seams, _ = native_loop_readout.carrier()
    rows = decision.complete_record_maps()
    actual = {(row["law"], tuple(row["seam"]), row["coin"], row["raw"]):
              tuple(map(tuple, row["pairs"])) for row in rows}
    assert len(rows) == len(actual) == 720
    for law, raw_count in (("C", 4), ("E", 8)):
        for seam in seams:
            for coin in (0, 1):
                expected = {raw: [] for raw in range(raw_count)}
                for i, state in enumerate(states):
                    after = _balance(state, seam, coin)
                    current = state[seam[0]]-after[seam[0]]
                    raw = {0: 0, 1: 1, -1: 2}[current]
                    if law == "E":
                        raw += 4*((state[seam[0]]+state[seam[1]]) % 2)
                    expected[raw].append((i, indices[after]))
                for raw, pairs in expected.items():
                    assert actual[law, seam, coin, raw] == tuple(pairs)
                    # Each recorded branch is a partial injection. With
                    # amplitude 1/sqrt(2), the two complete coin menus give
                    # the identity Kraus sum on the full 78-state domain.
                    assert len({target for _, target in pairs}) == len(pairs)
    assert sum(map(len, actual.values())) == 9360
    assert any(not pairs for pairs in actual.values())


def test_exact_reachability_has_only_one_closed_class_at_fixed_load_two():
    states = native.states()
    counts = decision.transitions()
    binary = {i for i, state in enumerate(states) if max(state) == 1}
    collisions = set(range(78)) - binary
    for i in binary:
        assert _reachable(counts, i) == binary
    for i in collisions:
        assert _reachable(counts, i) == binary | {i}
        assert counts[i][i] == 50
        assert sum(counts[i][j] for j in binary) == 10
    components = {frozenset(c) for c in decision.strongly_connected_components(counts)}
    assert components == {frozenset(binary)} | {frozenset({i}) for i in collisions}
    # A harmonic configuration readout is constant on the recurrent binary
    # class; f_i=(5/6)f_i+(1/6)c then fixes every collision value to the same c.
    assert all(F(counts[i][i], 60) != 1 for i in collisions)


def test_scc_algorithm_retains_directed_transient_components():
    counts = ((1, 1, 0, 0, 0), (1, 1, 1, 0, 0),
              (0, 0, 1, 0, 0), (0, 0, 1, 0, 1), (0, 0, 0, 1, 0))
    assert decision.strongly_connected_components(counts) == ((0, 1), (2,), (3, 4))


def test_collision_lifetime_and_potential_drop_use_the_same_transition():
    packet = decision.report()
    candidate = packet["collision_candidate"]
    survival = F(candidate["one_attempt_survival"])
    assert survival == F(5, 6)
    assert F(candidate["mean_attempts_to_split"]) == 1 / (1 - survival) == 6
    assert {int(k): F(v) for k, v in candidate["survival_by_attempt"].items()} == {
        n: survival ** n for n in range(13)}
    assert F(candidate["excess_quadratic_potential_before"]) == 1
    assert F(candidate["after_split"]) == 0
    seams, _ = native_loop_readout.carrier()
    for state in native.states():
        if max(state) != 2:
            continue
        before = F(sum(n*n for n in state), 2) - 1
        for seam in seams:
            for coin in (0, 1):
                after = _balance(state, seam, coin)
                drop = before - (F(sum(n*n for n in after), 2) - 1)
                assert drop == (1 if max(after) == 1 else 0)


def test_stationary_record_activity_is_not_quadratic_maintenance_cost():
    seams, _ = native_loop_readout.carrier()
    states = native.states(binary=True)
    activity = 0
    edge_changes = set()
    edge_cost = lambda loads: sum((loads[a]-loads[b])**2 for a, b in seams)
    for state in states:
        for seam in seams:
            for coin in (0, 1):
                after = _balance(state, seam, coin)
                activity += abs(state[seam[0]]-after[seam[0]])
                assert sum(n*n for n in state) == sum(n*n for n in after)
                edge_changes.add(edge_cost(after)-edge_cost(state))
    assert F(activity, 60*66) == F(5, 33)
    assert edge_changes == {-2, 0, 2}
    packet = decision.report()
    assert F(packet["transport_cost"]["stationary_nonzero_current_probability"]) == F(5, 33)
    assert list(map(F, packet["transport_cost"]["binary_quadratic_drop_values"])) == [F(0)]
    assert set(packet["transport_cost"]["binary_edge_cost_changes"]) == edge_changes
    assert packet["scope"]["classical_configuration_invariants_exhausted"]
    assert not packet["scope"]["all_quantum_or_record_encoded_sectors_excluded"]
    assert not packet["scope"]["physical_energy_or_mass_derived"]


def test_restricting_schedule_changes_sectors_and_is_not_the_original_protocol():
    seams, _ = native_loop_readout.carrier()
    retained = tuple(edge for edge in seams if 0 not in edge)
    states = native.states(binary=True)
    counts = decision.transitions(binary=True, seams=retained)
    assert counts == _counts(states, retained)
    occupied = {i for i, state in enumerate(states) if state[0] == 1}
    vacant = set(range(66))-occupied
    assert len(occupied) == 11 and len(vacant) == 55
    for part in (occupied, vacant):
        for i in part:
            assert _reachable(counts, i) == part
    assert sorted(map(len, decision.strongly_connected_components(counts))) == [11, 55]
    assert decision.report()["schedule_control"]["changed_source_schedule"]


def test_unrecovered_matching_code_escape_by_full_carrier_kraus_compression():
    states = native.states()
    indices = {state: i for i, state in enumerate(states)}
    # Integer columns C represent the declared code isometry C/sqrt(6).
    # The mixed code state is Pcode/6=C C^T/36 on the full 78-state carrier.
    c = np.zeros((78, 6), dtype=np.int64)
    for logical, support in enumerate(record_code.SUPPORTS):
        for state in support:
            c[indices[state], logical] = 1
    np.testing.assert_array_equal(c.T @ c, 6*np.eye(6, dtype=np.int64))
    retained_numerator = {"C": 0, "E": 0}
    for row in decision.complete_record_maps():
        m = np.zeros((78, 78), dtype=np.int64)
        for source, target in row["pairs"]:
            m[target, source] = 1
        compression = c.T @ m @ c
        retained_numerator[row["law"]] += int(np.sum(compression*compression))
    # K=M/sqrt(60) includes both the fair coin and uniform seam choice.
    # Tr(Pcode K (Pcode/6) K^T)=||C^T M C||_F^2/(60*6^3).
    escape = {law: 1-F(numerator, 60*6**3)
              for law, numerator in retained_numerator.items()}
    assert escape == {"C": F(31, 108), "E": F(55, 108)}
    packet = decision.corrected_memory_control()
    recorded = packet["unrecovered_uniform_C_E_code_escape_probability_from_Pcode_over_6"]
    assert {law: F(value) for law, value in recorded.items()} == escape
    assert packet["conditional_persistence_requires_existing_record_selected_recovery"]
    assert not packet["autonomous_defect_or_mass_selected"]


def test_quantum_collision_effect_and_entangled_reference_survival():
    states = native.states()
    collision = np.array([max(state) == 2 for state in states])
    projector = np.diag(collision.astype(np.int64))
    c = np.flatnonzero(collision)
    b = np.flatnonzero(~collision)
    # Integer complex amplitudes include coherence between collision sites,
    # collision/binary coherence, and entanglement with an untouched qubit.
    psi = np.zeros((78, 2), dtype=complex)
    psi[c[0], 0], psi[c[1], 0], psi[c[2], 1] = 1, 1j, 2
    psi[b[0], 0], psi[b[1], 1], psi[b[2], 1] = 1, 1j, 1
    reference = psi.conj().T @ psi
    norm = int(np.sum(abs(psi)**2))
    collision_norm = int(np.sum(abs(psi[collision])**2))
    assert norm == 9 and collision_norm == 6
    assert np.linalg.det(reference).real > 0
    effects = {law: np.zeros((78, 78), dtype=np.int64) for law in ("C", "E")}
    reference_sums = {law: np.zeros((2, 2), dtype=complex) for law in ("C", "E")}
    collision_sums = {"C": 0, "E": 0}
    branch_counts = {"C": 0, "E": 0}
    for row in decision.complete_record_maps():
        law = row["law"]
        m = np.zeros((78, 78), dtype=np.int64)
        for source, target in row["pairs"]:
            m[target, source] = 1
        selected = m[collision]
        # This is P_collision M (I-P_collision)=0, checked per raw outcome.
        assert not np.any(selected[:, ~collision])
        # Full 78-dimensional effect, including every off-diagonal entry.
        effects[law] += selected.T @ selected
        continued = m @ psi
        reference_sums[law] += continued.conj().T @ continued
        collision_sums[law] += int(np.sum(abs(continued[collision])**2))
        branch_counts[law] += 1
    witness = decision.report()["quantum_collision_survival"]["C_E"]
    for law in ("C", "E"):
        np.testing.assert_array_equal(effects[law], 50*projector)
        np.testing.assert_array_equal(reference_sums[law], 60*reference)
        assert F(collision_sums[law], 60*norm) == F(5, 6)*F(collision_norm, norm) == F(5, 9)
        nonzero = np.argwhere(effects[law])
        entries = [[int(i), int(j), int(effects[law][i, j])] for i, j in nonzero]
        assert witness[law]["integer_effect_entries"] == entries
        assert witness[law]["effect_denominator"] == 60
        assert witness[law]["nonzero_off_diagonal_entries"] == 0
        assert witness[law]["binary_to_collision_branch_addresses"] == 0
    assert branch_counts == {"C": 240, "E": 480}

"""Independent record projectors, complete source branches and hostile controls."""
from fractions import Fraction as F
from itertools import combinations_with_replacement

import numpy as np
import pytest

from . import record_recovery_source as recovery
from .source_discrimination import circuit
from .source_scalar import compiler as source


def _states():
    result = []
    for left, right in combinations_with_replacement(range(12), 2):
        state = [0] * 12
        state[left] += 1
        state[right] += 1
        result.append(tuple(state))
    return tuple(sorted(result))


def _label(candidate, state, coin):
    if candidate == "E":
        return state[0] - state[1]
    return state[0] - (sum(state[:2]) + coin) // 2


@pytest.mark.parametrize("candidate", ["C", "E"])
def test_all_input_record_projectors_have_exact_inverse_correspondence(candidate):
    report = recovery.exact_certificate(candidate)
    for coin in (0, 1):
        expected_survivors = sum(_label(candidate, x, coin) == _label(candidate, y, coin)
                                 for x in _states() for y in _states())
        assert report["coins"][coin] == {
            "coin": coin, "native_inputs": 78, "matrix_units_checked": 6084,
            "surviving_matrix_units": expected_survivors,
            "branch_amplitude_squared": "1/2"}
    assert report["archive_unchanged"]
    assert report["scratch_returns_blank"]
    assert report["arbitrary_reference_extension"]
    assert not report["cross_record_coherence_restored"]


@pytest.mark.parametrize("candidate, inverse_pulses, bits", [("C", 5370, 2), ("E", 5373, 3)])
def test_existing_inverse_words_copy_services_and_complete_raw_branches(candidate, inverse_pulses, bits):
    report = recovery.compiled_certificate(candidate)
    assert report["scratch_qubits"] == bits
    assert report["retained_classical_bits"] == bits + 1
    assert report["copy_cnot_words"] == bits
    assert report["blank_scratch_initializations"] == bits
    assert report["archive_record_rewrites"] == 0
    assert not report["certified_roundoff_or_physical_error_bound"]
    for coin in report["coins"]:
        assert coin["inverse_native_pulses"] == inverse_pulses
        assert coin["copy_native_pulses"] == 3 * bits
        assert coin["all_raw_branches_checked"] == 78 * 2**bits
        assert coin["retained_padding_abort_branches_checked"] == 78 * (1 if candidate == "C" else 3)
        assert coin["recovery_word_resources"] == {
            "native_pulses": inverse_pulses + 3 * bits,
            "proper_code_load_events": inverse_pulses + 3 * bits,
            "proper_code_export_events": inverse_pulses + 3 * bits}


@pytest.mark.parametrize("candidate", ["C", "E"])
def test_complete_executed_reference_branches_restore_only_history_projected_input(candidate):
    rng = np.random.default_rng(1096)
    states = _states()
    amplitudes = rng.normal(size=(78, 2)) + 1j * rng.normal(size=(78, 2))
    amplitudes /= np.linalg.norm(amplitudes)
    bits = circuit.RECORD_BITS[candidate]
    records = 2**bits
    total_mass = 0.
    for coin in (0, 1):
        forward = source.execute(circuit.repair_word(candidate, coin), 3 + bits)
        for raw in range(records):
            admitted = raw in circuit.record_values(candidate, coin)
            # Axes are archive, scratch, spectator code, endpoint code,
            # reference.  All branches and their unnormalized amplitudes
            # remain separate; no final marginal can hide a wrong block.
            actual = np.zeros((records, records, 64, 8, 2), complex)
            expected = np.zeros_like(actual)
            aborted = np.zeros((64, 8, 2), complex)
            for i, state in enumerate(states):
                code = circuit.encode(state)
                endpoint = forward[raw * 8:(raw + 1) * 8, code % 8]
                if not admitted:
                    # Different input endpoints with the same spectator
                    # contribute coherently to this single abort branch.
                    aborted[code // 8] += endpoint[:, None] * amplitudes[i] / np.sqrt(2)
                    continue
                restored = recovery.recover_local(candidate, coin, raw, endpoint)
                actual[:, :, code // 8] += (restored.reshape(records, records, 8, 1)
                                             * amplitudes[i] / np.sqrt(2))
                if circuit.record_values(candidate, coin)[raw] == _label(candidate, state, coin):
                    expected[raw, 0, code // 8, code % 8] += amplitudes[i] / np.sqrt(2)
            if admitted:
                assert np.linalg.norm(actual - expected) < 2e-11
                total_mass += np.linalg.norm(actual)**2
            else:
                # Keep numerical padding as an abort rather than normalizing
                # it away or quietly assigning it a valid decoded record.
                abort_mass = float(np.linalg.norm(aborted)**2)
                assert abort_mass < 1e-20
                total_mass += abort_mass
    assert abs(total_mass - 1) < 2e-11


def test_orthogonal_phase_pair_stays_identical_after_commitment_and_recovery():
    left = (1, 0, 1) + (0,) * 9
    right = (0, 1, 1) + (0,) * 9
    for candidate in ("C", "E"):
        for coin in (0, 1):
            assert _label(candidate, left, coin) != _label(candidate, right, coin)
        # Independent recorded projector action erases the only entries
        # in which the plus and minus density matrices differ.
        for sign in (-1, 1):
            restored = {}
            for coin in (0, 1):
                for i, x in enumerate((left, right)):
                    for j, y in enumerate((left, right)):
                        if _label(candidate, x, coin) == _label(candidate, y, coin):
                            restored[coin, _label(candidate, x, coin), i, j] = F(sign**(i + j), 4)
            if sign == -1:
                minus = restored
            else:
                assert restored == minus


def test_current_branch_keeps_a_phase_that_extra_parity_commitment_erases():
    x = (0, 0, 1, 1) + (0,) * 8
    z = (0, 1, 1) + (0,) * 9
    assert _label("C", x, 0) == _label("C", z, 0)
    assert _label("E", x, 0) != _label("E", z, 0)
    assert recovery.exact_certificate("C")["coins"][0]["surviving_matrix_units"] > (
        recovery.exact_certificate("E")["coins"][0]["surviving_matrix_units"])


@pytest.mark.parametrize("candidate", [None, True, 1, [], "c", "unknown"])
def test_invalid_candidates_rejected_before_cache_lookup(candidate):
    for operation in (lambda: recovery.inverse_word(candidate, 0),
                      lambda: recovery.copy_word(candidate),
                      lambda: recovery.exact_certificate(candidate),
                      lambda: recovery.compiled_certificate(candidate),
                      lambda: recovery.recover_local(candidate, 0, 0, np.zeros(8))):
        with pytest.raises(ValueError):
            operation()


@pytest.mark.parametrize("coin", [True, False, -1, 2, 0., None, "0"])
def test_invalid_coins_rejected_before_cache_lookup(coin):
    with pytest.raises(ValueError):
        recovery.inverse_word("C", coin)
    with pytest.raises(ValueError):
        recovery.recover_local("C", coin, 0, np.zeros(8))


@pytest.mark.parametrize("record", [True, False, -1, 3, 4, None, "0", 0.])
def test_invalid_current_records_are_not_recovered_as_valid_branches(record):
    with pytest.raises(ValueError):
        recovery.recover_local("C", 0, record, np.zeros(8))


@pytest.mark.parametrize("state", [[], np.zeros(7), np.zeros((8, 1)), [float("nan")] * 8,
                                    [float("inf")] * 8, ["0"] * 8,
                                    np.ma.array(np.zeros(8), mask=[True] + [False] * 7)])
def test_invalid_amplitudes_rejected(state):
    with pytest.raises(ValueError):
        recovery.recover_local("C", 0, 0, state)


@pytest.mark.parametrize("part", ["forward", "inverse", "copy"])
def test_relative_phase_corruption_is_rejected_even_if_probabilities_match(monkeypatch, part):
    original = recovery._executed
    position = {"forward": 0, "inverse": 1, "copy": 2}[part]
    def corrupted(candidate, coin):
        matrices = [matrix.copy() for matrix in original(candidate, coin)]
        matrices[position][:, 0] *= -1
        return tuple(matrices)
    monkeypatch.setattr(recovery, "_executed", corrupted)
    with pytest.raises(ArithmeticError, match="phase-exact"):
        recovery.compiled_certificate("C")


def test_wrong_inverse_word_is_rejected_and_no_stale_cache_masks_it(monkeypatch):
    recovery._executed.cache_clear()
    monkeypatch.setattr(recovery, "inverse_word", lambda candidate, coin: [])
    try:
        with pytest.raises(ArithmeticError, match="phase-exact"):
            recovery.compiled_certificate("C")
    finally:
        recovery._executed.cache_clear()


def test_public_state_output_cannot_mutate_cached_source_matrices():
    initial = np.eye(8)[0]
    first = recovery.recover_local("C", 0, 0, initial)
    first[:] = 99
    second = recovery.recover_local("C", 0, 0, initial)
    assert np.linalg.norm(second) == pytest.approx(1, abs=1e-11)
    for matrix in recovery._executed("C", 0):
        assert not matrix.flags.writeable


def test_audit_scope_excludes_coherent_tape_erasure_and_unexecuted_seams():
    result = recovery.audit()
    assert result["executed_seam"] == [0, 1]
    assert "abstract port relabelling only" in result["other_seams"]
    assert not result["new_protected_quantum_memory"]
    assert not result["new_hardware_execution"]

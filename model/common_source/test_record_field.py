"""Independent all-sector transport, finite records and state-update checks."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import record_field as field
from ..observer.source_scalar import compiler


def test_source_walk_matches_signed_exterior_lift_on_every_fock_sector():
    actual = field.source_walk()
    expected = field.second_quantization(field.permutation())
    assert np.max(abs(actual-expected)) < 2e-11
    assert np.max(abs(actual.conj().T@actual-np.eye(64))) < 4e-11
    # The source program transports every charged annihilator with its signs,
    # rather than merely reproducing the one-particle occupation probabilities.
    for source, destination in enumerate(field.permutation()):
        assert np.max(abs(actual@field.annihilator(source).T@actual.conj().T
                             - field.annihilator(destination).T)) < 4e-11


def test_ring_edge_swap_retains_intervening_fermionic_parity():
    permutation = list(range(6))
    permutation[0], permutation[5] = permutation[5], permutation[0]
    expected = field.second_quantization(permutation)
    actual = compiler.execute(field.swap_word(0, 5), 6)
    assert np.max(abs(actual-expected)) < 2e-11
    # A particle between the endpoints changes the exchange phase. A bit swap
    # followed by CZ on the endpoints alone would incorrectly return +1 here.
    assert expected[(1 << 3)+(1 << 5), (1 << 0)+(1 << 3)] == -1


def test_full_car_and_two_layer_local_continuity():
    report = field.algebra_and_continuity()
    assert report["all_CAR_and_continuity_identities_exact"]
    assert report["total_corrected_code_dimension"] == 46656
    assert report["CAR_representation_multiplicity"] == 729
    assert report["code_directions_discarded_per_carrier"] == 0
    # All six-code spectator directions support a CAR qubit. Direct matrix
    # verification prevents silently replacing I6 by a rank-two projector.
    c = np.kron(np.eye(3), np.array([[0, 1], [0, 0]]))
    assert np.array_equal(c@c.T+c.T@c, np.eye(6))


def test_every_lattice_momentum_and_both_orientations_share_one_clock():
    result = field.dispersion()
    rows = result["modes"]
    assert len(rows) == 6
    assert {(row["momentum_integer"], row["orientation"]) for row in rows} == {
        (k, orientation) for k in (-1, 0, 1) for orientation in (-1, 1)}
    assert result["charged_retarded_anticommutator_response_residual"] < 3e-11
    assert result["nontrivial_spatial_dispersion"]
    assert not result["physical_chirality_or_3plus1_spacetime_claim"]


@pytest.mark.parametrize("law,alphabet,positive", [("C", 8, 4), ("E", 16, 6)])
def test_complete_native_record_tensor_factors_keep_raw_zero_outputs(law, alphabet, positive):
    result = field.record_factors(law)
    assert len(result["accepted_seams"]) == 30
    assert result["complete_native_record_addresses_per_accepted_cycle"] == 30*alphabet**6
    for row in result["accepted_seams"]:
        branches = row["local_record_factor"]
        assert len(branches) == alphabet
        weights = [F(branch["weight_exact"]) for branch in branches]
        assert sum(weights) == 1
        assert sum(w > 0 for w in weights) == positive
        assert any(branch["abort"] and F(branch["weight_exact"]) == 0 for branch in branches)
    # The complete one-cycle selector/native probability includes both halt
    # addresses; six independent reference-entangled factors cause no loss.
    mass = sum(F(row["selector_weight_exact"])*sum(F(b["weight_exact"])
               for b in row["local_record_factor"])**6 for row in result["accepted_seams"])
    mass += sum(F(row["weight_exact"]) for row in result["terminal_selector_addresses"])
    assert mass == 1


def test_incident_clock_keeps_all_waiting_and_terminal_paths():
    for cap in field.CAPS:
        histories = field.clock_histories(cap)
        assert sum(F(row["weight_exact"]) for row in histories) == 1
        assert {row["flag"] for row in histories} == {"success", "selector_timeout", "deadline_timeout"}
        successes = [row for row in histories if row["flag"] == "success"]
        assert {row["elapsed_cycles"] for row in successes} == {1, 2}
        assert all(row["incident_ticks"] == 1 for row in successes)
        assert any(row["elapsed_cycles"] == 0 and row["flag"] == "selector_timeout" for row in histories)
        for row in histories:
            assert row["elapsed_cycles"] == sum(address < 30 for address in row["selector_addresses"])
            assert F(row["weight_exact"]) == F(1, 32**len(row["selector_addresses"]))
            assert all(address < 30 for address in row["selector_addresses"][:-1])


def test_two_read_feedback_changes_future_joint_law_and_keeps_reference():
    report, arrays = field.detector_experiment()
    outcomes = report["complete_detector_outputs"]
    assert {(r["first_read"], r["second_read"]): F(r["unconditional_probability_exact"])
            for r in outcomes} == {(0, 0): F(225, 512), (0, 1): F(), (1, 0): F(), (1, 1): F(225, 512)}
    assert np.allclose(report["no_feedback_detector_probabilities_before_selector_weights"], [.5, 0, .5, 0], atol=3e-11)
    assert arrays["initial"].shape == (64, 2)
    assert np.allclose(arrays["initial"].conj().T@arrays["initial"], np.eye(2)/2)
    for result in (0, 1):
        conditional = arrays[f"after_first_{result}"]
        # The apparatus learns one local occupation and changes the quantum
        # reference correlations. Native recovery alone would preserve rank2.
        assert np.linalg.matrix_rank(conditional, tol=1e-10) == 1
    probability = sum(F(r["unconditional_probability_exact"]) for r in outcomes)
    probability += F(report["first_selector_timeout_probability_exact"])
    probability += 2*F(report["second_selector_timeout_probability_each_first_read_exact"])
    assert probability == 1
    weighted = [value for name, value in arrays.items()
                if name.startswith("weighted_terminal_") or name.startswith("selector_timeout_")]
    assert abs(sum(float(np.vdot(value, value).real) for value in weighted)-1) < 4e-11
    assert np.allclose(report["identity_walk_control_detector_probabilities_before_selector_weights"], [.5, .5, 0, 0], atol=3e-11)


def test_both_observers_read_one_shared_tape_and_source_advance():
    result = field.joint_clock_histories()
    assert result["one_source_tape_for_both_observers"]
    histories = result["histories"]
    assert len(histories) == 2+30*32
    assert sum(F(row["weight_exact"]) for row in histories) == 1
    seen_disagreement = False
    for row in histories:
        accepted = [address for address in row["selector_addresses"] if address < 30]
        assert len(accepted) == row["elapsed_cycles"]
        expected = []
        for cap in field.CAPS:
            expected.append(next((cycle for cycle, address in enumerate(accepted, 1)
                                  if set(cap).intersection(field.record_code.SEAMS[address])), None))
        assert row["first_incident_cycle_by_observer"] == expected
        seen_disagreement |= expected[0] != expected[1]
    assert seen_disagreement


def test_field_low_bit_wiring_preserves_every_native_spectator_and_syndrome_bit():
    word = field.decoded_register_walk_word()
    assert {target for event in word for target in event["targets"]} == {0, 9, 18, 27, 36, 45}
    assert len(word) == len(field.walk_word())


def test_source_executed_clock_channel_retains_every_terminal_gram_block():
    clock = field.clock_report()
    joint = clock["joint_clock_histories"]
    actual = joint["source_executed_stopped_channel"]
    assert actual["complete_clock_block_count"] == 14
    assert actual["all_work_input_columns_per_clock_block"] == 64
    assert actual["complete_retained_stopped_operator_error"] < 5e-11
    assert actual["complete_stopped_gram_error"] < 1e-10
    assert sum(F(row["weight_exact"]) for row in actual["terminal_source_characteristics"]) == 1
    assert {row["elapsed_cycles"] for row in actual["terminal_source_characteristics"]} == {0, 1, 2}


@pytest.mark.parametrize("phase", [1, -1])
def test_clock_channel_rejects_identity_and_wrong_global_phase(monkeypatch, phase):
    original = field.source_walk()
    false = np.eye(64) if phase == 1 else -original
    monkeypatch.setattr(field, "source_walk", lambda: false)
    with pytest.raises(ArithmeticError, match="stopped source clock"):
        field.clock_report()


@pytest.mark.parametrize("law,archive", [("C", 48), ("E", 60)])
def test_physical_service_counts_are_explicit_without_fabricated_error_values(law, archive):
    supplied_pulses = {"native_pulses": 17, "proper_code_load_events": 17,
                       "proper_code_export_events": 17}
    services = field.declared_services(law, supplied_pulses)
    assert services["maximum_archive_bits"] == archive
    assert services["maximum_archive_and_clock_bits"] == archive+6
    assert services["native_attempts_and_coin_record_pairs"] == 12
    coefficients = services["physical_half_diamond_error_sum_coefficients"]
    assert coefficients["fair_classical_selector_bit"] == 10
    assert coefficients["local_field_number_read_and_commit"] == 2
    assert coefficients["native_coin_or_raw_record_commit"] == 24
    assert coefficients["proper_code_coherent_load"] == 17
    assert not services["physical_epsilon_values_available"]
    assert not services["calibrated_physical_duration_energy_or_total_error_available"]


def test_cached_source_walk_cannot_be_mutated_by_caller():
    with pytest.raises(ValueError):
        field.source_walk().flags.writeable = True


@pytest.mark.parametrize("value", [[0., 1., 2., 3., 4., 5.], [False, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 4]])
def test_permutation_rejects_aliases_and_nonpermutations(value):
    with pytest.raises(ValueError):
        field.second_quantization(value)


@pytest.mark.parametrize("cap", [(False,), (0, 0), (12,), (), "01"])
def test_clock_access_region_validation(cap):
    with pytest.raises(ValueError):
        field.clock_histories(cap)

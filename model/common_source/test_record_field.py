"""Independent all-sector transport, finite records and state-update checks."""
from fractions import Fraction as F
from copy import deepcopy
import math

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
    assert field._field_block_word(word) == field.walk_word()
    targets = list(range(6))
    targets[2], targets[3] = targets[3], targets[2]
    assert np.max(abs(field.source_feedback()-field.second_quantization(targets))) < 2e-11


@pytest.mark.parametrize("action", ["walk", "feedback"])
@pytest.mark.parametrize("fault", ["reordered", "syndrome"])
def test_actual_decoded_register_word_is_checked_before_tensor_factorization(monkeypatch, action, fault):
    name = f"decoded_register_{action}_word"
    execute = getattr(field, f"source_{action}")
    original = getattr(field, name)()
    if fault == "reordered":
        # For feedback the exchange 2<->3 itself is symmetric, so swap its
        # target2 with external target0 to obtain a wrong but unitary operation.
        a, b = (0, 9) if action == "walk" else (0, 18)
        wrong = [{**event, "targets": [b if j == a else a if j == b else j
                                       for j in event["targets"]]} for event in original]
    else:
        wrong = deepcopy(original)
        wrong[0]["targets"][0] += 1
    monkeypatch.setattr(field, name, lambda: wrong)
    execute.cache_clear()
    try:
        with pytest.raises((ArithmeticError, ValueError), match="source |spectator or syndrome"):
            execute()
    finally:
        execute.cache_clear()


def test_field_blocks_lower_the_actual_word_order_not_just_its_target_set():
    original = field.decoded_register_walk_word()
    wrong = [{**event, "targets": [9 if j == 0 else 0 if j == 9 else j
                                   for j in event["targets"]]} for event in original]
    assert {j for event in wrong for j in event["targets"]} == {0, 9, 18, 27, 36, 45}
    actual = compiler.execute(field._field_block_word(wrong), 6)
    assert np.max(abs(actual.conj().T@actual-np.eye(64))) < 4e-11
    assert np.linalg.norm(actual-field.second_quantization(field.permutation()), 2) > 1.9


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


@pytest.mark.parametrize("field_name,bad", [
    ("elapsed_cycles", True),
    ("first_incident_cycle_by_observer", [-1, None]),
    ("first_incident_cycle_by_observer", [True, None]),
    ("first_incident_cycle_by_observer", [3, None]),
    ("first_incident_cycle_by_observer", [None]),
    ("global_terminal_flag", "success"),
    ("path_count", True), ("path_count", 0),
    ("weight_exact", True), ("weight_exact", .0625),
    ("weight_exact", "1/0"), ("weight_exact", "2/32"),
])
def test_stopped_clock_rejects_original_metadata_before_execution(monkeypatch, field_name, bad):
    groups = field._joint_clock_groups(field.joint_clock_histories()["histories"])
    groups[0][field_name] = bad
    def forbidden():
        raise AssertionError("source execution preceded clock input validation")
    monkeypatch.setattr(field, "source_walk", forbidden)
    with pytest.raises(ValueError):
        field.stopped_clock_execution(groups)


@pytest.mark.parametrize("fault", ["future_hit", "duplicate", "false_count", "exchanged_weights"])
def test_normalized_clock_blocks_do_not_authorize_false_tape_metadata(monkeypatch, fault):
    groups = field._joint_clock_groups(field.joint_clock_histories()["histories"])
    if fault == "future_hit":
        row = next(group for group in groups if group["elapsed_cycles"] == 0)
        row["first_incident_cycle_by_observer"] = [1, None]
    elif fault == "duplicate":
        copied = deepcopy(groups[0])
        copied["weight_exact"] = groups[0]["weight_exact"] = str(F(copied["weight_exact"])/2)
        groups.append(copied)
    elif fault == "false_count":
        groups[0]["path_count"] += 1
    else:
        groups[0]["weight_exact"], groups[-1]["weight_exact"] = groups[-1]["weight_exact"], groups[0]["weight_exact"]
    assert sum(F(group["weight_exact"]) for group in groups) == 1
    monkeypatch.setattr(field, "source_walk", lambda: pytest.fail("invalid tape reached execution"))
    with pytest.raises(ValueError):
        field.stopped_clock_execution(groups)


def test_caller_cannot_replace_owned_native_budget_by_a_self_consistent_zero(monkeypatch):
    owned = {"encoder_log_error": 1e-10, "decoder_log_error": 1e-10}
    monkeypatch.setattr(field.record_code, "source_composition_report", lambda: owned)
    supplied = {"encoder_log_error": 0., "decoder_log_error": 0.}
    with pytest.raises(ValueError, match="owned source calculation"):
        field.source_composition_report(supplied)


@pytest.mark.parametrize("law,archive", [("C", 48), ("E", 60)])
def test_physical_service_counts_are_explicit_without_fabricated_error_values(law, archive):
    supplied_pulses = {"native_pulses": 17, "proper_code_load_events": 17,
                       "proper_code_export_events": 17}
    services = field.declared_services(law, supplied_pulses)
    assert services["maximum_archive_bits"] == archive
    assert services["maximum_archive_and_clock_bits"] == archive+6
    assert services["maximum_proper_code_export_history_bits"] == 5
    assert services["maximum_archive_clock_and_export_history_bits"] == archive+11
    assert services["native_attempts_and_coin_record_pairs"] == 12
    coefficients = services["physical_half_diamond_error_sum_coefficients"]
    assert coefficients["fair_classical_selector_bit"] == 10
    assert coefficients["local_field_number_read_and_commit"] == 2
    assert coefficients["native_coin_or_raw_record_commit"] == 24
    assert coefficients["proper_code_coherent_load"] == 17
    assert coefficients["first_export_failure_position_commit"] == 1
    assert coefficients["classical_storage_bit_for_whole_window"] == archive+11
    assert not services["physical_epsilon_values_available"]
    assert not services["calibrated_physical_duration_energy_or_total_error_available"]


def test_export_menu_covers_inverse_macros_all_mixing_sizes_and_native_words():
    report = field._proper_code_export_diagnostic()
    counts = report["pulse_family_counts"]
    assert len(counts) == 18
    assert all(counts[f"mixing_{count}_{direction}"] > 0
               for count in (1, 4, 5, 6) for direction in ("forward", "inverse"))
    assert all(counts[f"native_{law}_coin{coin}"] > 0 for law in ("C", "E") for coin in (0, 1))
    assert counts["routing_macro_forward"] == counts["routing_macro_inverse"] > 0
    assert report["actualpulse_menu_count"] > 0
    assert 0 < report["maximum_complement_operator_norm"] < 2e-11
    assert 0 <= report["maximum_success_norm_excess"] < 2e-11
    assert not report["numeric_complement_is_set_to_zero"]


def test_correct_success_block_does_not_hide_computed_export_failure(monkeypatch):
    original = compiler.native_pulse
    base = field._proper_code_export_diagnostic()
    leak = 1e-12
    def leaky(kind, angle):
        full = original(kind, angle)
        if kind == "I" and angle == math.pi/2:
            retained = full[:2, :2].copy()
            full[2:, :2] = 0
            full[2, 0] = leak
            np.testing.assert_array_equal(full[:2, :2], retained)
        return full
    monkeypatch.setattr(compiler, "native_pulse", leaky)
    changed = field._proper_code_export_diagnostic()
    assert changed["maximum_complement_operator_norm"] == leak
    assert field._first_export_failure_bound(100, changed) > 1000*field._first_export_failure_bound(100, base)


def test_two_pulse_complete_failure_operator_has_orthogonal_position_blocks():
    # Independent literal complete dilation. Both failures must be retained;
    # replacing their stack by one branch loses the sqrt(number) factor.
    theta = .2
    success, failure = math.cos(theta)*np.eye(2), math.sin(theta)*np.eye(2)
    failed = np.concatenate((failure, failure@success))
    actual = float(np.linalg.norm(failed, 2))
    expected = math.sqrt(-math.expm1(4*math.log(math.cos(theta))))
    bound = field._first_export_failure_bound(2, {
        "maximum_complement_operator_norm": math.sin(theta),
        "maximum_success_norm_excess": 0.})
    assert abs(actual-expected) < 1e-15
    assert math.sin(theta) < actual <= bound
    gram = (success@success).T@(success@success)+failed.T@failed
    np.testing.assert_allclose(gram, np.eye(2), atol=1e-15)
    # A nonunitary floating success norm must inflate later failure prefixes.
    nu, leak = 1.01, .03
    failed = np.concatenate((leak*np.eye(2), leak*nu*np.eye(2)))
    inflated = field._first_export_failure_bound(2, {
        "maximum_complement_operator_norm": leak,
        "maximum_success_norm_excess": nu-1})
    assert math.sqrt(2)*leak < np.linalg.norm(failed, 2) <= inflated
    assert np.linalg.norm(np.kron(failed, np.eye(3)), 2) <= inflated


@pytest.mark.parametrize("fault", ["leak", "unitary_leak", "success_norm", "nan"])
def test_export_diagnostic_rejects_bad_full_source_pulses(monkeypatch, fault):
    original = compiler.native_pulse
    def bad(kind, angle):
        full = original(kind, angle)
        if kind == "I" and angle == math.pi/2:
            if fault == "leak":
                full[2, 0] = 1e-8
            elif fault == "unitary_leak":
                rotation = np.eye(6)
                theta = 1e-6
                rotation[0, 0] = rotation[2, 2] = math.cos(theta)
                rotation[0, 2], rotation[2, 0] = -math.sin(theta), math.sin(theta)
                changed = rotation@full
                # Full unitarity and the successful-block tolerance still
                # pass; only checking the complementary columns exposes it.
                assert np.linalg.norm(changed.conj().T@changed-np.eye(6), 2) < 1e-13
                assert np.linalg.norm(changed[:2, :2]-full[:2, :2], 2) < 1e-12
                full = changed
            elif fault == "success_norm":
                full[:2, :2] *= 1+1e-8
            else:
                full[2, 0] = np.nan
        return full
    monkeypatch.setattr(compiler, "native_pulse", bad)
    with pytest.raises(ArithmeticError, match="proper-code export"):
        field._proper_code_export_diagnostic()


def test_cached_source_walk_cannot_be_mutated_by_caller():
    with pytest.raises(ValueError):
        field.source_walk().flags.writeable = True
    with pytest.raises(ValueError):
        field.source_feedback().flags.writeable = True


@pytest.mark.parametrize("value", [[0., 1., 2., 3., 4., 5.], [False, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 4]])
def test_permutation_rejects_aliases_and_nonpermutations(value):
    with pytest.raises(ValueError):
        field.second_quantization(value)


@pytest.mark.parametrize("cap", [(False,), (0, 0), (12,), (), "01"])
def test_clock_access_region_validation(cap):
    with pytest.raises(ValueError):
        field.clock_histories(cap)

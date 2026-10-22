"""Independent endpoint, loop and mutation checks for native current transport."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch

from model.observer import native, native_loop_readout as subject


def direct_event(loads, seam, coin):
    a, b = seam
    low, odd = divmod(loads[a] + loads[b], 2)
    final_a = low + odd * coin
    return loads[a] - final_a


class NativeLoopReadoutTests(unittest.TestCase):
    def test_all_native_records_select_an_exact_gradient(self):
        seams, faces = subject.carrier()
        cycles, chords = subject.cycle_basis(seams)
        self.assertEqual((len(seams), len(faces), len(cycles), len(chords)),
                         (30, 20, 19, 19))
        for loads in native.states():
            measured = subject.read_contrasts(loads, seams)
            for a, b in seams:
                independent = direct_event(loads, (a, b), 0) + direct_event(loads, (a, b), 1)
                self.assertEqual(measured[a, b], independent)
                self.assertEqual(independent, loads[a] - loads[b])
            self.assertTrue(all(subject.integer_holonomy(measured, loop) == 0
                                for loop in cycles))
            self.assertTrue(all(subject.integer_holonomy(measured, face + (face[0],)) == 0
                                for face in faces))

    def test_chord_flux_control_detects_a_loop_but_is_not_native(self):
        seams, _ = subject.carrier()
        cycles, chords = subject.cycle_basis(seams)
        loads = tuple(int(i in (0, 1)) for i in range(12))
        measured = subject.read_contrasts(loads, seams)
        altered = dict(measured)
        altered[chords[0]] += 1
        self.assertEqual([i for i, loop in enumerate(cycles)
                          if subject.integer_holonomy(altered, loop)], [0])
        self.assertNotEqual(altered[chords[0]], loads[chords[0][0]] - loads[chords[0][1]])

    def test_source_branch_mutation_fails_before_physical_promotion(self):
        original = native.native
        loads = tuple(int(i in (0, 1)) for i in range(12))
        seams, _ = subject.carrier()
        _, chords = subject.cycle_basis(seams)

        def changed(source, seam, coin):
            final, record = original(source, seam, coin)
            if source == loads and seam == chords[0] and coin == 0:
                return final, (coin, record[1] + 1)
            return final, record

        with patch.object(native, "native", side_effect=changed):
            with self.assertRaisesRegex(ArithmeticError, "not a port-potential gradient"):
                subject.run()

    def test_binary_transfer_is_exactly_reversible(self):
        seams, _ = subject.carrier()
        counts, (i, j, k) = subject.binary_transfer(seams)
        self.assertEqual(len(counts), 66)
        self.assertEqual(sum(bool(counts[a][b]) for a in range(66)
                             for b in range(a + 1, 66)), 300)
        self.assertTrue(all(counts[a][b] == counts[b][a]
                            for a in range(66) for b in range(66)))
        self.assertEqual(counts[i][j] * counts[j][k] * counts[k][i],
                         counts[i][k] * counts[k][j] * counts[j][i])

    def test_source_executed_path_eraser_separates_record_phase_from_holonomy(self):
        row = subject.same_endpoint_phase_readout()
        self.assertEqual(row["returning_tape_codes"], [0, 21])
        self.assertEqual(set(row["all_four_committed_word_probabilities"].values()), {"1/4"})
        self.assertEqual(row["source_executed_eraser"]["base"]
                         ["exact_returning_plus_minus_probabilities"], ["1/2", "0"])
        self.assertEqual(row["source_executed_eraser"]["first_coin_phase"]
                         ["exact_returning_plus_minus_probabilities"], ["0", "1/2"])
        self.assertTrue(row["first_coin_scalar_phase_leaves_every_committed_CP_branch_unchanged"])
        self.assertEqual(row["complete_C_E_adaptive_branches_checked"], {"C": 312, "E": 312})
        self.assertTrue(row["reversed_seam_pair_has_unit_wilson_product_for_any_U1_link_assignment"])
        self.assertFalse(row["carrier_face_loop_or_nontrivial_link_holonomy"])
        self.assertFalse(row["physical_defect_or_mass_identification"])
        self.assertTrue(all(x["all_outcome_probability_mass_within_2e_10"]
                            and x["source_readout_probability_residual_below_2e_10"]
                            for x in row["source_executed_eraser"].values()))

    def test_broken_eraser_or_missing_record_phase_fails_exact_contrast(self):
        with patch.object(subject.circuit, "basis_swap_word", return_value=[]):
            with self.assertRaisesRegex(ArithmeticError, "phase contrast"):
                subject.same_endpoint_phase_readout()
        original = subject.compiler.pulse

        def missing_phase(kind, angle, targets):
            if kind == "Z" and targets == [0]:
                return original("Z", 0, targets)
            return original(kind, angle, targets)

        with patch.object(subject.compiler, "pulse", side_effect=missing_phase):
            with self.assertRaisesRegex(ArithmeticError, "phase contrast"):
                subject.same_endpoint_phase_readout()

    def test_registered_equal_length_face_routes_keep_every_native_history(self):
        row = subject.registered_face_interferometer()
        self.assertEqual(row["closed_registered_route"], [0, 1, 5, 11, 0])
        self.assertEqual(row["adjacent_registered_faces"], [[0, 5, 1], [0, 11, 5]])
        self.assertEqual(row["two_equal_length_routes"], [[0, 1, 5], [0, 11, 5]])
        self.assertEqual(len(row["all_eight_native_histories"]), 8)
        self.assertEqual({entry["committed_probability"] for entry in row["all_eight_native_histories"]},
                         {"1/8"})
        self.assertEqual(row["successful_history_codes"], [36, 85])
        self.assertEqual(row["success_probability_before_eraser"], "1/4")
        self.assertEqual(row["other_terminal_probability"], "3/4")
        self.assertFalse(row["source_selected_defect_holonomy"])
        self.assertFalse(row["native_repair_source_executes_individual_link_unitaries"])
        for case in row["readout_cases"].values():
            for read in case["quadratures"].values():
                all_outcomes = read["all_endpoint_and_read_code_probabilities"]
                self.assertEqual(sum(subject.Q(value) for bins in all_outcomes.values()
                                     for value in bins.values()), 1)
                self.assertEqual(sum(subject.Q(value) for endpoint, bins in all_outcomes.items()
                                     if endpoint != "2,5" for value in bins.values()),
                                 subject.Q(3, 4))
        receipt = json.loads((Path(__file__).with_name("results") / "native_loop_readout.json")
                             .read_text(encoding="utf-8"))
        original = row["readout_cases"]["flat"]["quadratures"]["X"]
        compact = receipt["registered_face_interferometer"]["readout_cases"]["flat"]["quadratures"]["X"]
        self.assertEqual(compact["complete_outcome_distribution_sha256"], hashlib.sha256(
            subject.canonical(original["all_endpoint_and_read_code_probabilities"])
            .encode("utf-8")).hexdigest())
        self.assertEqual(compact["complete_probability_mass"], "1")
        self.assertNotIn("all_endpoint_and_read_code_probabilities", compact)
        compact_face = receipt["registered_face_interferometer"]
        self.assertEqual(compact_face["native_history_count"], 8)
        self.assertEqual(compact_face["complete_committed_probability_mass"], "1")
        self.assertEqual(compact_face["complete_native_history_distribution_sha256"], hashlib.sha256(
            subject.canonical(row["all_eight_native_histories"]).encode("utf-8")).hexdigest())

    def test_two_quadratures_resolve_supplied_flux_but_native_character_is_flat(self):
        rows = subject.registered_face_interferometer()["readout_cases"]
        for case in ("flat", "native_current_character"):
            self.assertEqual(rows[case]["closed_loop_quarter_turns_mod_4"], 0)
            self.assertEqual(rows[case]["quadratures"]["X"]["exact_success_endpoint_plus_minus"],
                             ["1/4", "0"])
            self.assertEqual(rows[case]["quadratures"]["Y"]["exact_success_endpoint_plus_minus"],
                             ["1/8", "1/8"])
        for case in ("inserted_link_phase", "port_rephased_inserted_link"):
            self.assertEqual(rows[case]["closed_loop_quarter_turns_mod_4"], 1)
            self.assertEqual(rows[case]["quadratures"]["X"]["exact_success_endpoint_plus_minus"],
                             ["1/8", "1/8"])
            self.assertEqual(rows[case]["quadratures"]["Y"]["exact_success_endpoint_plus_minus"],
                             ["1/4", "0"])

    def test_link_gauge_reversal_and_backtrack_are_independent_path_identities(self):
        links = {(5, 11): 1}
        path = (0, 1, 5, 11, 0)
        self.assertEqual(subject.oriented_link_sum(links, path), 1)
        self.assertEqual(subject.oriented_link_sum(links, path[::-1]), -1)
        self.assertEqual(subject.oriented_link_sum(links, (5, 11, 5)), 0)
        shifted = subject.gauge_shift(links, tuple(range(12)))
        self.assertEqual(subject.oriented_link_sum(shifted, path), 1)
        self.assertEqual(subject.oriented_link_sum(shifted, (0, 1, 5))
                         - subject.oriented_link_sum(shifted, (0, 11, 5)), 1)
        with self.assertRaisesRegex(ValueError, "outside the native carrier"):
            subject.oriented_link_sum(links, (0, 2))
        with self.assertRaisesRegex(ValueError, "twelve integer port phases"):
            subject.gauge_shift(links, (0, 1))

    def test_one_face_idle_control_sees_uniform_minus_invisible_to_even_loop(self):
        row = subject.registered_triangle_control()
        self.assertEqual(row["registered_face"], [0, 5, 1])
        self.assertEqual(row["comparison_loop"], [0, 1, 5, 0])
        self.assertEqual(row["native_attempt_counts"], [2, 2])
        self.assertEqual(row["moving_event_counts_on_selected_arms"], [2, 1])
        self.assertEqual(row["second_attempt_on_short_route"], [2, 3])
        self.assertEqual(row["empty_idle_current_on_both_coins"], 0)
        self.assertEqual(row["selected_path_codes"], [36, 5])
        self.assertEqual(row["unselected_idle_coin_code_at_same_terminal"], 21)
        self.assertEqual(len(row["all_eight_native_histories"]), 8)
        self.assertEqual({item["committed_probability"] for item in row["all_eight_native_histories"]},
                         {"1/8"})
        self.assertTrue(row["uniform_sign_gives_two_face_loop_unit_but_one_face_minus_unit"])
        self.assertTrue(row["idle_phase_assumed_zero_in_declared_coherent_lift"])
        self.assertFalse(row["idle_and_move_count_dynamical_phases_independently_calibrated"])
        self.assertFalse(row["phase_difference_identified_as_geometric_holonomy"])
        self.assertFalse(row["source_selected_one_face_holonomy"])
        self.assertEqual(row["readout_cases"]["flat"]["quadratures"]["X"]
                         ["exact_selected_plus_minus"], ["1/4", "0"])
        self.assertEqual(row["readout_cases"]["uniform_minus_moving_link_control"]
                         ["quadratures"]["X"]["exact_selected_plus_minus"], ["0", "1/4"])
        for case in row["readout_cases"].values():
            for read in case["quadratures"].values():
                bins = read["all_endpoint_and_read_code_probabilities"]
                self.assertEqual(sum(subject.Q(value) for endpoint in bins.values()
                                     for value in endpoint.values()), 1)
                self.assertEqual(sum(subject.Q(value) for code, value in bins["5,7"].items()
                                     if code not in ("0", "1")), subject.Q(1, 8))
        receipt = json.loads((Path(__file__).with_name("results") / "native_loop_readout.json")
                             .read_text(encoding="utf-8"))
        original = row["readout_cases"]["uniform_minus_moving_link_control"]["quadratures"]["X"]
        compact = (receipt["registered_triangle_control"]["readout_cases"]
                   ["uniform_minus_moving_link_control"]["quadratures"]["X"])
        self.assertEqual(compact["complete_outcome_distribution_sha256"], hashlib.sha256(
            subject.canonical(original["all_endpoint_and_read_code_probabilities"])
            .encode("utf-8")).hexdigest())
        self.assertEqual(compact["complete_probability_mass"], "1")
        compact_triangle = receipt["registered_triangle_control"]
        self.assertEqual(compact_triangle["native_history_count"], 8)
        self.assertEqual(compact_triangle["complete_committed_probability_mass"], "1")
        self.assertEqual(compact_triangle["complete_native_history_distribution_sha256"], hashlib.sha256(
            subject.canonical(row["all_eight_native_histories"]).encode("utf-8")).hexdigest())

    def test_missing_declared_phase_or_eraser_fails_flux_contrast(self):
        original = subject.compiler.pulse

        def missing_phase(kind, angle, targets):
            if kind == "Z" and targets == [0] and math.isclose(angle, math.pi / 4):
                return original("Z", 0, targets)
            return original(kind, angle, targets)

        with patch.object(subject.compiler, "pulse", side_effect=missing_phase):
            with self.assertRaisesRegex(ArithmeticError, "quadrature"):
                subject.registered_face_interferometer()
        with patch.object(subject.circuit, "basis_swap_word", return_value=[]):
            with self.assertRaisesRegex(ArithmeticError, "quadrature"):
                subject.registered_face_interferometer()

    def test_receipt_replays_and_keeps_scope(self):
        result = subject.run()
        self.assertEqual(result["schema"], "oph.native_loop_readout.v2")
        self.assertEqual(result["ideal_all_seam_readout"]["branch_records_checked"], 4680)
        self.assertEqual(result["status"],
                         "NATIVE_CURRENT_CONTRAST_FLAT__COHERENT_HOLONOMY_OPEN")
        self.assertEqual(result["source"]["registered_carrier_binding"]
                         ["registered_gram_entries_checked_exactly"], 144)
        self.assertFalse(result["source"]["mass_or_particle_target_data_consumed"])
        self.assertEqual(result["uniform_seam_binary_transfer"]["nonzero_unordered_transition_pairs"], 300)
        self.assertEqual(result["fixed_observer"]["same_endpoint_path_control"]
                         ["separate_retained_record_words"],
                         [[[0, 0], [0, 0]], [[1, -1], [0, 1]]])
        self.assertFalse(result["fixed_observer"]["same_endpoint_path_control"]
                         ["relative_phase_readout_in_committed_tape"])
        self.assertEqual(result["fixed_observer"]["coherent_path_record_control"]
                         ["source_pulse_counts"], {"first_coin_phase": 2,
                                                  "record_eraser": 2707})
        self.assertEqual(result["registered_face_interferometer"]["successful_history_codes"], [36, 85])
        self.assertFalse(result["registered_face_interferometer"]["source_selected_defect_holonomy"])
        self.assertEqual(result["registered_triangle_control"]["selected_path_codes"], [36, 5])
        self.assertFalse(result["registered_triangle_control"]["source_selected_one_face_holonomy"])
        receipt = Path(__file__).with_name("results") / "native_loop_readout.json"
        self.assertEqual(receipt.read_text(encoding="utf-8"), subject.canonical(result))


if __name__ == "__main__":
    unittest.main()

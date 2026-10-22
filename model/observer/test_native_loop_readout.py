"""Independent endpoint, loop and mutation checks for native current transport."""
from __future__ import annotations

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

    def test_receipt_replays_and_keeps_scope(self):
        result = subject.run()
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
        receipt = Path(__file__).with_name("results") / "native_loop_readout.json"
        self.assertEqual(receipt.read_text(encoding="utf-8"), subject.canonical(result))


if __name__ == "__main__":
    unittest.main()

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
        receipt = Path(__file__).with_name("results") / "native_loop_readout.json"
        self.assertEqual(receipt.read_text(encoding="utf-8"), subject.canonical(result))


if __name__ == "__main__":
    unittest.main()

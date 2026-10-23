"""Independent graph, gauge, calibration and actual source-readout checks."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from model.observer import holonomy_reconstruction as h


class HolonomyReconstructionTests(unittest.TestCase):
    def test_every_nineteen_face_menu_recovers_gauge_class(self):
        seams, _, _, _, _ = h.incidence()
        connection = {edge: F(i**3 + 2, 113) for i, edge in enumerate(seams)}
        potentials = [F(i * i, 71) for i in range(12)]
        regauged = {(a, b): phase + potentials[a] - potentials[b]
                    for (a, b), phase in connection.items()}
        phases = h.face_phases(connection)
        expected = h.tree_representative(connection)
        self.assertEqual(h.face_phases(regauged), phases)
        self.assertEqual(h.tree_representative(regauged), expected)
        self.assertEqual(sum(phases) % 1, 0)
        for omitted in range(20):
            menu = [value for i, value in enumerate(phases) if i != omitted]
            self.assertEqual(h.reconstruct_faces(menu, omitted), expected)

    def test_inverse_is_integer_identity_not_float_rank(self):
        self.assertEqual(h.integer_inverse(((1, 1), (0, 1))), ((1, -1), (0, 1)))
        with self.assertRaisesRegex(ValueError, "not integral"):
            h.integer_inverse(((2,),))
        with self.assertRaisesRegex(ValueError, "not independent"):
            h.integer_inverse(((1, 1), (1, 1)))
        with self.assertRaisesRegex(ValueError, "nineteen"):
            h.reconstruct_faces([0] * 18)

    def test_missing_face_is_fixed_by_closure_but_missing_two_are_not(self):
        # A changed chord has opposite flux on exactly its two incident faces;
        # omitting both hides that non-gauge perturbation completely.
        seams, _, _, chords, _ = h.incidence()
        flat = {edge: F(0) for edge in seams}
        changed = dict(flat)
        changed[chords[0]] = F(1, 7)
        face_changes = h.face_phases(changed)
        visible = [i for i, value in enumerate(face_changes) if value]
        self.assertEqual(len(visible), 2)
        self.assertNotEqual(h.tree_representative(flat), h.tree_representative(changed))
        self.assertTrue(all(face_changes[i] == 0 for i in range(20) if i not in visible))

    def test_orientation_cancels_common_offset_and_one_anchor_resolves_sign(self):
        seams = h.incidence()[0]
        connection = {edge: F(i * i, 67) for i, edge in enumerate(seams)}
        phases = h.face_phases(connection)
        forward, reverse = h.paired_orientations(phases, F(9, 17))
        self.assertEqual(tuple((a - b) % 1 for a, b in zip(forward, reverse)),
                         tuple((2 * value) % 1 for value in phases))
        candidates = h.recover_common_offset(forward, reverse)
        self.assertEqual(len(candidates), 2)
        self.assertIn((phases, F(9, 17)), candidates)
        for anchor_face in range(20):
            self.assertEqual(h.recover_common_offset(forward, reverse,
                                                     (anchor_face, phases[anchor_face])),
                             ((phases, F(9, 17)),))
        corrupted = list(reverse)
        corrupted[3] += F(1, 19)
        with self.assertRaisesRegex(ValueError, "common move/idle"):
            h.recover_common_offset(forward, corrupted)
        with self.assertRaisesRegex(ValueError, "calibration"):
            h.recover_common_offset(forward, reverse, (0, phases[0] + F(1, 3)))

    def test_all_coherent_paths_still_cannot_split_arbitrary_seam_dynamics(self):
        seams, _, _, chords, _ = h.incidence()
        connection = {edge: F(i, 31) for i, edge in enumerate(seams)}
        dynamics = {edge: F(i * i, 41) for i, edge in enumerate(seams)}
        total = h.directed_phases(connection, dynamics)
        for edge in chords:
            shifted, shifted_dynamics = dict(connection), dict(dynamics)
            shifted[edge] += F(1, 2)
            shifted_dynamics[edge] += F(1, 2)
            self.assertEqual(h.directed_phases(shifted, shifted_dynamics), total)
            self.assertNotEqual(h.tree_representative(shifted), h.tree_representative(connection))

    def test_quadratures_retain_success_mass_and_reject_no_visibility(self):
        self.assertEqual(h.quadrature_character((F(1, 8), F(1, 8)), (F(1, 4), 0)),
                         (0, 1))
        with self.assertRaisesRegex(ValueError, "success mass"):
            h.quadrature_character((F(1, 4), 0), (F(1, 4), F(1, 4)))
        with self.assertRaisesRegex(ValueError, "pure unit"):
            h.quadrature_character((F(1, 8), F(1, 8)), (F(1, 8), F(1, 8)))
        with self.assertRaisesRegex(ValueError, "success mass"):
            h.quadrature_character((0, 0), (0, 0))

    def test_entire_native_face_menu_uses_actual_routes_and_retains_all_branches(self):
        protocols = h.native_face_menu()
        self.assertEqual(len(protocols), 40)
        self.assertEqual({(row["face"], row["orientation"]) for row in protocols},
                         {(face, orientation) for face in range(20) for orientation in (-1, 1)})
        for row in protocols:
            self.assertEqual(len(row["histories"]), 8)
            self.assertEqual(len({(history["route"], *history["coins"])
                                  for history in row["histories"]}), 8)
            self.assertFalse(set(row["idle_seam"]) & (set(row["routes"][0]) | {row["spectator"]}))

    def test_receipt_replays_actual_source_triangle_readout(self):
        result = h.run()
        self.assertEqual(result["source_executed_triangle_characters"],
                         {"flat": ["1", "0"], "uniform_minus_moving_link_control": ["-1", "0"]})
        self.assertEqual(result["squared_holonomy_unresolved_gauge_classes"], 524288)
        self.assertFalse(result["source_selected_coherent_transport"])
        self.assertTrue(result["all_pass"])
        self.assertEqual(Path(h.DEFAULT_RECEIPT).read_text(encoding="utf-8"), h.canonical(result))


if __name__ == "__main__":
    unittest.main()

"""Exact controls for selecting a phase-sensitive native seam transport."""
from __future__ import annotations

from pathlib import Path
import unittest

from model.observer import native, native_link_selection as selection
from model.observer.source_discrimination import certificate


class NativeLinkSelectionTests(unittest.TestCase):
    def test_full_carrier_committed_map_agrees_with_independent_certificate(self):
        states = native.states()
        for law in ("C", "E"):
            for left in states:
                for right in states:
                    columns_l = selection._branches(left, (0, 1), law)
                    columns_r = selection._branches(right, (0, 1), law)
                    expected_l = tuple(certificate._endpoint_image(left, law, coin)
                                       for coin in (0, 1))
                    expected_r = tuple(certificate._endpoint_image(right, law, coin)
                                       for coin in (0, 1))
                    expected = certificate._matrix_unit_image(expected_l, expected_r)
                    observed = selection._matrix_unit(columns_l, columns_r,
                                                       negative_moving_phase=True)
                    self.assertEqual(observed, expected)

    def test_two_symmetry_uniform_phase_completions_disagree_on_every_face(self):
        result = selection.run()
        self.assertEqual(result["complete_committed_matrix_units_checked"], 2 * 78**2)
        self.assertEqual(result["one_token_moving_columns_checked"], 2 * 30 * 2)
        self.assertEqual(result["all_plus_completion_face_products"], "+1 on all 20 faces")
        self.assertEqual(result["all_minus_moving_branch_completion_face_products"],
                         "-1 on all 20 faces")
        self.assertEqual(result["two_face_0_1_5_11_0_loop_product_in_both_completions"], "+1")
        self.assertFalse(result["two_face_even_loop_detects_uniform_pi_per_triangle"])
        self.assertFalse(result["source_selected_or_observed_native_link_holonomy"])
        self.assertTrue(result["all_pass"])
        self.assertEqual(result["schema"], "oph.native_link_selection.v1")
        receipt = Path(__file__).with_name("results") / "native_link_selection.json"
        self.assertEqual(receipt.read_text(encoding="utf-8"), selection.canonical(result))

    def test_moving_branch_is_partial_isometry_not_two_port_unitary(self):
        for law in ("C", "E"):
            forward = selection.moving_branch((0, 1), 0, law)
            reverse = selection.moving_branch((0, 1), 1, law)
            self.assertNotEqual(forward, reverse)
            self.assertEqual(selection._link_sign(0, 1, law, True), -1)
            self.assertEqual(selection._link_sign(1, 0, law, True), -1)
            self.assertEqual(selection._loop_sign((0, 1, 0), law, True), 1)
        with self.assertRaisesRegex(ValueError, "endpoint"):
            selection.moving_branch((0, 1), 5, "C")


if __name__ == "__main__":
    unittest.main()

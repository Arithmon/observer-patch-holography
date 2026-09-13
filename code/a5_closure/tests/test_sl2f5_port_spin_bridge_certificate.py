from __future__ import annotations

import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import sl2f5_port_spin_bridge_certificate as bridge  # noqa: E402


class SL2F5PortSpinBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = bridge.build_certificate()

    def test_explicit_isomorphism_has_120_elements(self) -> None:
        iso = self.certificate["isomorphism"]
        self.assertEqual(iso["source_elements_mapped"], 120)
        self.assertEqual(iso["distinct_spin_images"], 120)
        self.assertEqual(iso["multiplication_checks"], 120 * 120)
        self.assertTrue(iso["bijective"])

    def test_cover_square_commutes_on_every_source_element(self) -> None:
        square = self.certificate["commuting_cover_square"]
        self.assertEqual(square["checks"], 120)
        self.assertEqual(square["source_kernel"], ["+I", "-I"])
        self.assertEqual(square["spin_kernel"], ["+I2", "-I2"])
        self.assertTrue(square["commutes_on_all_source_elements"])

    def test_bridge_uses_constructive_presentation_not_group_profile_inference(self) -> None:
        source = self.certificate["source"]
        spin = self.certificate["port_spin_lift"]
        boundary = self.certificate["claim_boundary"]
        self.assertEqual(source["presentation_solution_count"], 120)
        self.assertEqual(source["chosen_generator_orders"], [4, 6, 10])
        self.assertEqual(len(spin["base_generator_lift_signs"]), 3)
        self.assertEqual(
            boundary["method_not_used"],
            [
                "same-order inference",
                "element-order-profile classification",
                "uniqueness of the non-split central extension",
            ],
        )

    def test_exact_spinor_representation_and_center_map(self) -> None:
        spin = self.certificate["port_spin_lift"]
        iso = self.certificate["isomorphism"]
        self.assertEqual(spin["order"], 120)
        self.assertEqual(spin["center"], ["+I2", "-I2"])
        self.assertEqual(spin["exact_field"], "Q(sqrt(5), i)")
        self.assertTrue(spin["faithful_two_dimensional_complex_representation"])
        self.assertEqual(iso["center_map"], "+I -> +I2; -I -> -I2")


if __name__ == "__main__":
    unittest.main()

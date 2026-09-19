#!/usr/bin/env python3
"""Tests for SOURCE-CURRENT-TOMOGRAPHY-0 stage 0."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import source_current_tomography_stage0 as audit  # noqa: E402


class SourceCurrentTomographyStage0Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = audit.audit_payload()

    def test_fails_closed_on_current_public_source_information(self) -> None:
        self.assertEqual(self.payload["status"], "INSUFFICIENT_SOURCE_DATA")

    def test_source_bound_constraints_are_preserved(self) -> None:
        constraints = self.payload["source_bound_constraints"]
        self.assertTrue(constraints["carrier_bound"])
        self.assertTrue(constraints["target_blind_impulse_readback"])
        self.assertTrue(constraints["response_sector_signs"])
        self.assertTrue(constraints["orientation_convention"])
        self.assertTrue(constraints["refinement_maps"])

    def test_missing_information_is_current_tomography_not_sector_signs(self) -> None:
        requirements = self.payload["tomography_requirements"]
        self.assertFalse(requirements["ordered_two_sided_response_histories"])
        self.assertFalse(requirements["twelve_source_reconstructed_infinitesimal_generators"])
        self.assertFalse(requirements["source_reconstructed_positive_pairing"])
        self.assertFalse(requirements["source_reconstructed_commutator"])
        self.assertFalse(requirements["same_history_overlap_words"])
        self.assertFalse(requirements["same_current_projective_implementers"])

    def test_oriented_face_bracket_is_a_real_negative_control(self) -> None:
        control = self.payload["negative_control"]
        self.assertEqual(control["oriented_face_count"], 20)
        self.assertEqual(control["proper_action_order"], 60)
        self.assertEqual(control["jacobi_nonzero_coordinate_count"], 240)
        self.assertEqual(control["jacobi_squared_norm"], 240)
        self.assertTrue(control["fails_jacobi"])

    def test_b14_repair_is_not_promoted_to_source_selection(self) -> None:
        control = self.payload["negative_control"]
        self.assertFalse(control["nearest_compact_family_is_source_selection"])
        self.assertFalse(control["metric_is_source_derived"])
        self.assertFalse(control["repair_rule_is_source_derived"])

    def test_claim_boundary_is_non_promoting(self) -> None:
        boundary = self.payload["claim_boundary"]
        self.assertTrue(boundary["does_not_reject_conditional_port_current_algebra"])
        self.assertTrue(boundary["does_not_select_charged_double_triplet_fixture"])
        self.assertTrue(boundary["does_not_make_physical_current_claim"])


if __name__ == "__main__":
    unittest.main()

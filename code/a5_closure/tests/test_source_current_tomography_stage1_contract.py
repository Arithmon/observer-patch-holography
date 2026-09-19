#!/usr/bin/env python3
"""Tests for SOURCE-CURRENT-TOMOGRAPHY-0 Stage 1."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import source_current_tomography_stage1_contract as stage1  # noqa: E402


class SourceCurrentTomographyStage1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = stage1.payload()

    def test_registered_algebra_is_exactly_four_dimensional_and_commutative(self) -> None:
        algebra = self.payload["registered_response_word_algebra"]
        self.assertEqual(algebra["dimension"], 4)
        self.assertTrue(algebra["commutative"])
        self.assertEqual(algebra["commutator_nonzero_count"], 0)
        self.assertTrue(algebra["positive_skew_pairing"])

    def test_exact_deficits_are_quantified(self) -> None:
        gap = self.payload["exact_gap"]
        self.assertEqual(gap["independent_tangent_directions_missing"], 8)
        self.assertEqual(gap["derived_commutator_dimensions_missing"], 11)
        self.assertEqual(gap["nonidentity_rechartings_missing_from_response_words"], 59)
        self.assertTrue(gap["same_word_projective_implementers_missing"])

    def test_order_sensitive_packet_records_both_orders(self) -> None:
        ordered = self.payload["source_packet_contract"]["ordered_histories"]
        self.assertTrue(ordered["required"])
        self.assertTrue(ordered["both_composition_orders_per_unordered_port_pair"])
        self.assertTrue(ordered["same_source_packet"])
        self.assertEqual(ordered["unordered_port_pairs"], 66)\n        self.assertEqual(ordered["ordered_compositions"], 132)

    def test_first_order_and_mixed_order_gates_are_nontrivial(self) -> None:
        contract = self.payload["source_packet_contract"]
        self.assertEqual(contract["first_order_gate"]["real_rank_exactly"], 12)
        self.assertEqual(contract["mixed_order_gate"]["commutator_span_real_rank"], 11)
        self.assertEqual(contract["mixed_order_gate"]["center_dimension"], 1)

    def test_holonomy_gate_requires_all_sixty_rechartings(self) -> None:
        holonomy = self.payload["source_packet_contract"]["holonomy_gate"]
        self.assertEqual(holonomy["proper_rechartings_covered"], 60)
        self.assertTrue(holonomy["same_history_projective_implementers"])
        self.assertTrue(holonomy["internal_modulo_pointwise_centralizer"])

    def test_claim_boundary_remains_non_promoting(self) -> None:
        boundary = self.payload["claim_boundary"]
        self.assertTrue(boundary["bounded_obstruction_only"])
        self.assertTrue(boundary["does_not_assert_no_go_for_order_sensitive_histories"])
        self.assertTrue(boundary["does_not_consume_conditional_current_fixture"])
        self.assertTrue(boundary["does_not_promote_physical_current"])
        self.assertTrue(boundary["abstract_forced_lie_type_theorem_preserved"])

    def test_nonidentifiability_is_an_allowed_successful_scientific_exit(self) -> None:
        exits = self.payload["allowed_exits"]
        self.assertIn("inequivalent", exits["nonidentifiable"])
        self.assertIn("explicit witnesses", exits["nonidentifiable"])


if __name__ == "__main__":
    unittest.main()

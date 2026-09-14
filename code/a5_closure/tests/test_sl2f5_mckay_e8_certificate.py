#!/usr/bin/env python3
"""Regression tests for the exact SL(2,F5) McKay/E8 certificate."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import sl2f5_mckay_e8_certificate as cert  # noqa: E402


class SL2F5McKayE8CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_certificate()

    def test_recovers_complete_irreducible_dimension_multiset(self) -> None:
        recovery = self.payload["irreducible_recovery"]
        self.assertEqual(recovery["irreducible_count"], 9)
        self.assertEqual(sorted(recovery["dimensions"].values()), [1, 2, 2, 3, 3, 4, 4, 5, 6])
        self.assertEqual(recovery["sum_squared_dimensions"], 120)
        self.assertFalse(recovery["hard_coded_character_table"])

    def test_certified_spin_doublet_has_affine_e8_mckay_graph(self) -> None:
        mckay = self.payload["mckay"]
        self.assertGreaterEqual(mckay["affine_e8_isomorphism_count"], 1)
        self.assertEqual(len(mckay["edges"]), 8)
        self.assertTrue(mckay["dimension_vector_is_two_eigenvector"])
        self.assertTrue(mckay["connected_tree"])

    def test_center_parity_splits_five_quotient_and_four_spinorial_irreps(self) -> None:
        signs = self.payload["irreducible_recovery"]["center_signs"]
        self.assertEqual(sum(sign == 1 for sign in signs.values()), 5)
        self.assertEqual(sum(sign == -1 for sign in signs.values()), 4)
        self.assertEqual(signs["spin2"], -1)
        self.assertEqual(signs["spin2_galois"], -1)

    def test_galois_conjugate_doublet_is_faithful_and_also_affine_e8(self) -> None:
        control = self.payload["galois_control"]
        self.assertTrue(control["faithful"])
        self.assertEqual(control["kernel_size"], 1)
        self.assertEqual(control["conjugate_doublet"], "spin2_galois")
        self.assertTrue(control["same_affine_e8_graph_type"])
        self.assertGreaterEqual(control["affine_e8_isomorphism_count"], 1)
        self.assertFalse(control["affine_e8_graph_type_selects_galois_embedding"])
        self.assertTrue(control["sqrt5_sensitive_spin_trace_values"])

    def test_dimension_two_alone_does_not_force_affine_e8(self) -> None:
        control = self.payload["negative_control"]
        self.assertEqual(control["dimension"], 2)
        self.assertEqual(control["affine_e8_isomorphism_count"], 0)
        self.assertTrue(control["passes"])

    def test_character_table_has_nine_derived_conjugacy_classes(self) -> None:
        table = self.payload["character_table"]
        self.assertEqual(len(table), 9)
        self.assertEqual(sum(row["size"] for row in table), 120)
        self.assertTrue(all(len(row["characters"]) == 9 for row in table))


if __name__ == "__main__":
    unittest.main()

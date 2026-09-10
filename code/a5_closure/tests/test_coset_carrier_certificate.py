#!/usr/bin/env python3
"""Regression and adversarial tests for the coset carrier certificate."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import coset_carrier_certificate as cert  # noqa: E402


class CosetCarrierCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path, cls.receipt_path, cls.negative_path = cert.default_paths()
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)
        cls.built, _ = cert.construct(cls.manifest)
        cls.quotient = cls.built.quotient

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_negative_control_bundle_is_exactly_recomputable(self) -> None:
        self.assertEqual(cert.load_json(self.negative_path), cert.negative_control_payload(self.manifest))

    def test_placement_classes(self) -> None:
        self.assertEqual(
            self.expected["classification"]["placement_classes"],
            {"CARRIER": 120, "FACE_BOUNDARY_MISMATCH": 180, "FACES_NOT_TRIANGLES": 300, "COLLAPSED_EDGES": 300},
        )

    def test_local_counts_are_the_same_for_every_placement(self) -> None:
        counts = {cert.local_counts(cert.build_incidence(self.quotient, t)) for t in self.built.verdicts}
        self.assertEqual(counts, {((5,), (3,), (5,), (2,), (2,), (3,))})

    def test_carriers_are_exactly_the_placements_with_xyz_equal_1(self) -> None:
        for t, v in self.built.verdicts.items():
            generators = cert.triangle_generators(self.quotient, t)
            self.assertEqual(v == "CARRIER", bool(generators))
            if generators:
                self.assertEqual(len(generators), 1)

    def test_presentation_supplies_exactly_the_carriers(self) -> None:
        solutions, image = cert.presentation_placements(self.quotient)
        self.assertEqual(solutions, 120)
        self.assertEqual(image, set(self.built.carriers))

    def test_every_carrier_relabels_onto_the_committed_packet(self) -> None:
        committed = cert.relabel_faces(self.built.fixture.oriented_faces, range(12))
        for t in self.built.carriers:
            inc = cert.build_incidence(self.quotient, t)
            phi, stats = cert.relabel(self.quotient, t, inc, self.built.fixture)
            self.assertEqual(stats, {"graph_isomorphisms": 120, "orientation_matching_relabellings": 60})
            self.assertEqual(cert.relabel_faces(cert.oriented_faces(inc), phi), committed)

    def test_boundary_mismatch_family_has_the_committed_graph_and_face_sets(self) -> None:
        committed_sets = {frozenset(f) for f in self.built.fixture.oriented_faces}
        family = [t for t, v in self.built.verdicts.items() if v == "FACE_BOUNDARY_MISMATCH"]
        self.assertEqual(len(family), 180)
        for t in family:
            inc = cert.build_incidence(self.quotient, t)
            isos = cert.graph_isomorphisms(cert.adjacency(inc), self.built.fixture.adjacency)
            self.assertTrue(any({frozenset(phi[v] for v in f) for f in inc.face_vertices} == committed_sets for phi in isos))

    def test_every_negative_control_failed_with_its_expected_code(self) -> None:
        rows = cert.load_json(self.negative_path)["finite_controls"]
        self.assertTrue(rows)
        for row in rows:
            self.assertTrue(row["passed"], row["name"])
            self.assertEqual(row["expected_error"], row["actual_error"], row["name"])

    def test_corrupted_edge_face_incidence_is_rejected(self) -> None:
        t = self.built.carriers[0]
        inc = cert.build_incidence(self.quotient, t)
        self.assertEqual(cert.classify(inc), "CARRIER")
        wrong = next(f for f, verts in enumerate(inc.face_vertices) if not inc.edge_vertices[0] <= verts)
        bad = cert.with_edge_faces(inc, 0, frozenset([min(inc.edge_faces[0]), wrong]))
        self.assertEqual(cert.verdict(bad), "FACE_BOUNDARY_MISMATCH")

    def test_fixture_with_every_face_reversed_is_matched(self) -> None:
        f = self.built.fixture
        reversed_fixture = cert.Fixture(f.adjacency, f.port_action, tuple((p, r, q) for p, q, r in f.oriented_faces))
        t = self.built.carriers[0]
        phi, _ = cert.relabel(self.quotient, t, cert.build_incidence(self.quotient, t), reversed_fixture)
        self.assertEqual(sorted(phi), list(range(12)))

    def test_tampered_receipt_is_rejected(self) -> None:
        tampered = copy.deepcopy(self.expected)
        tampered["classification"]["placement_classes"]["CARRIER"] = 300
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_receipt(self.manifest, tampered)
        self.assertEqual(ctx.exception.code, "RECEIPT_MISMATCH")

    def test_source_block_cannot_name_target_structure(self) -> None:
        bad = copy.deepcopy(self.manifest)
        bad["source"]["group"] = "SL2 icosahedral"
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.validate_manifest(bad)
        self.assertEqual(ctx.exception.code, "SOURCE_FIREWALL")

    def test_missing_lean_marker_fails_closed(self) -> None:
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.lean_block("no declaration here", "def neighbors", "\n\n", "Lean/Screen/PortFrameGram.lean")
        self.assertEqual(ctx.exception.code, "LEAN_FIXTURE_PARSE")


if __name__ == "__main__":
    unittest.main()

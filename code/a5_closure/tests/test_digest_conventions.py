#!/usr/bin/env python3
"""Invariants for the three sha256 conventions used by the package pins.

The `sha256:` tag and the field names name none of the three serializations, so
these tests fix the properties a reader needs in order to check a pin without a
label: every pin in the package resolves to exactly one convention, a digest
that was not produced by any of them resolves to none, and the three
conventions disagree on one intact file.
"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import digest_conventions as dc  # noqa: E402


class ConventionArithmeticTests(unittest.TestCase):
    """The three serializations, on one intact artifact."""

    SELF_SIGNING_ARTIFACT = MODULE_DIR / "manifests" / "global_form_semantic_artifact.json"

    def test_three_conventions_give_three_different_digests(self) -> None:
        path = self.SELF_SIGNING_ARTIFACT
        raw = dc.digest_raw_bytes(path)
        canonical = dc.digest_canonical_json(path)
        stripped = dc.digest_canonical_json_self_digest_stripped(path)
        self.assertEqual(len(stripped), 1, "the artifact signs itself with one field")
        (stripped_digest,) = stripped.values()
        self.assertEqual(
            len({raw, canonical, stripped_digest}),
            3,
            "the three conventions must not collapse onto one value",
        )

    def test_canonical_json_is_sorted_compact_and_newline_free(self) -> None:
        payload = {"b": 1, "a": [2, {"d": 3, "c": 4}]}
        self.assertEqual(
            dc.canonical_json_bytes(payload),
            b'{"a":[2,{"c":4,"d":3}],"b":1}',
        )

    def test_self_digest_field_is_found_without_trusting_its_name(self) -> None:
        fields = dc.self_digest_fields(
            json.loads(self.SELF_SIGNING_ARTIFACT.read_text(encoding="utf-8"))
        )
        self.assertEqual(fields, ("artifact_sha256",))

    def test_raw_bytes_convention_reads_the_committed_bytes(self) -> None:
        # The canonical JSON carries no trailing newline, so a committed file
        # that ends in one cannot have the same digest under both conventions.
        path = self.SELF_SIGNING_ARTIFACT
        self.assertTrue(path.read_bytes().endswith(b"\n"))
        self.assertNotEqual(dc.digest_raw_bytes(path), dc.digest_canonical_json(path))

    def test_tag_is_ignored_on_both_sides(self) -> None:
        path = self.SELF_SIGNING_ARTIFACT
        bare = dc.digest_canonical_json(path)
        self.assertEqual(
            dc.resolve_convention(bare, path),
            dc.resolve_convention(f"sha256:{bare}", path),
        )
        self.assertEqual(
            dc.resolve_convention(bare, path), (dc.CANONICAL_JSON,)
        )


class PackagePinSweepTests(unittest.TestCase):
    """Every pin in this package, classified by serialization alone."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.results = dc.classify_pins()
        cls.summary = dc.summarize(cls.results)

    def test_the_sweep_finds_pins(self) -> None:
        self.assertGreater(self.summary["pins"], 50)

    def test_each_located_pin_resolves_to_exactly_one_convention(self) -> None:
        offenders = [
            (row.pin.location, row.pin.declared_path, row.conventions)
            for row in self.results
            if row.located and len(row.conventions) != 1
        ]
        self.assertEqual(offenders, [])

    def test_no_pin_is_unreproduced(self) -> None:
        self.assertEqual(
            [row.pin.location for row in self.summary["unreproduced"]], []
        )

    def test_no_pin_is_ambiguous(self) -> None:
        self.assertEqual([row.pin.location for row in self.summary["ambiguous"]], [])

    def test_all_three_conventions_are_in_use(self) -> None:
        for name in dc.CONVENTIONS:
            self.assertGreater(
                self.summary["split"][name], 0, f"no pin uses {name}"
            )

    def test_the_split_accounts_for_every_located_pin(self) -> None:
        self.assertEqual(sum(self.summary["split"].values()), self.summary["pins"])

    def test_summary_line_reports_a_clean_sweep(self) -> None:
        line = dc.summary_line(self.summary)
        self.assertRegex(line, r"^\d+ path-paired pins, 0 unreproduced, 0 ambiguous$")


class CorruptedPinTests(unittest.TestCase):
    """A pin no convention produced resolves to none of them."""

    ARTIFACT = MODULE_DIR / "manifests" / "global_form_semantic_artifact.json"

    def test_a_corrupted_pin_resolves_to_no_convention(self) -> None:
        good = dc.digest_canonical_json(self.ARTIFACT)
        # Flip one hex digit, so the value keeps the shape of a pin.
        head = "0" if good[0] != "0" else "1"
        corrupted = head + good[1:]
        self.assertTrue(dc.looks_like_digest(corrupted))
        self.assertEqual(dc.resolve_convention(corrupted, self.ARTIFACT), ())
        self.assertEqual(
            dc.resolve_convention(f"sha256:{corrupted}", self.ARTIFACT), ()
        )

    def test_a_digest_of_a_different_file_resolves_to_no_convention(self) -> None:
        other = MODULE_DIR / "manifests" / "axis_center_descent_reference.json"
        self.assertEqual(
            dc.resolve_convention(dc.digest_canonical_json(other), self.ARTIFACT), ()
        )

    def test_the_indented_serialization_is_the_committed_byte_form(self) -> None:
        # The producers write indented, key-sorted JSON with a trailing newline,
        # so a digest over that serialization is the raw-byte convention rather
        # than a fourth one.  A reader who indents the canonical JSON instead of
        # compacting it lands on the wrong digest of the two JSON conventions.
        document = json.loads(self.ARTIFACT.read_text(encoding="utf-8"))
        indented = (
            json.dumps(document, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
        ).encode("utf-8")
        self.assertEqual(indented, self.ARTIFACT.read_bytes())
        self.assertEqual(
            dc.resolve_convention(dc.sha256_hex(indented), self.ARTIFACT),
            (dc.RAW_BYTES,),
        )
        self.assertNotEqual(
            dc.sha256_hex(indented), dc.digest_canonical_json(self.ARTIFACT)
        )


class ResolverBoundaryTests(unittest.TestCase):
    """The resolver refuses to guess which file a pin names."""

    def test_a_path_outside_the_two_bases_is_not_located(self) -> None:
        self.assertIsNone(
            dc.locate("local_domain/stage1_event_complex.py", (dc.REPO_ROOT, MODULE_DIR))
        )

    def test_both_pin_path_bases_resolve(self) -> None:
        roots = (dc.REPO_ROOT, MODULE_DIR)
        self.assertIsNotNone(
            dc.locate("receipts/axis_center_descent_reference.receipt.json", roots)
        )
        self.assertIsNotNone(
            dc.locate(
                "code/a5_closure/receipts/axis_center_descent_reference.receipt.json",
                roots,
            )
        )

    def test_a_bare_basename_would_collide_with_a_local_artifact(self) -> None:
        # classical_realization_receipt.json and clock_unit_verdict.json pin an
        # upstream producer packet by bare basename.  One of those basenames,
        # matter_attachment_receipt.json, also names a file in manifests/, and
        # that file is a different artifact: resolving the pin against the
        # containing directory would report the local artifact as corrupted.
        collision = MODULE_DIR / "manifests" / "matter_attachment_receipt.json"
        self.assertTrue(collision.is_file())
        verdict = json.loads(
            (MODULE_DIR / "manifests" / "clock_unit_verdict.json").read_text(
                encoding="utf-8"
            )
        )
        pin = verdict["upstream_pins"]["matter_attachment_receipt.json"]
        self.assertEqual(dc.resolve_convention(pin, collision), ())
        self.assertIsNone(
            dc.locate("matter_attachment_receipt.json", (dc.REPO_ROOT, MODULE_DIR))
        )

    def test_a_non_json_file_yields_only_the_raw_byte_digest(self) -> None:
        producer = MODULE_DIR / "baryon_dimension_six_census.py"
        digests = dc.convention_digests(producer)
        self.assertEqual(len(digests[dc.RAW_BYTES]), 1)
        self.assertEqual(digests[dc.CANONICAL_JSON], ())
        self.assertEqual(digests[dc.CANONICAL_JSON_SELF_DIGEST_STRIPPED], ())


class KnownPinTests(unittest.TestCase):
    """The named examples that show the tag names no convention."""

    AXIS_RECEIPT = (
        MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json"
    )

    def test_the_same_receipt_is_pinned_raw_and_canonically(self) -> None:
        census = json.loads(
            (MODULE_DIR / "receipts" / "baryon_dimension_six_census.receipt.json")
            .read_text(encoding="utf-8")
        )
        raw_pin = census["upstream_pins"]["axis_center_receipt"]["sha256"]
        common_ew = json.loads(
            (MODULE_DIR / "manifests" / "common_ew_order_unit_carrier_reference.json")
            .read_text(encoding="utf-8")
        )
        canonical_pin = common_ew["upstream_pins"]["global_form_receipt"]["sha256"]
        self.assertEqual(
            dc.resolve_convention(raw_pin, self.AXIS_RECEIPT), (dc.RAW_BYTES,)
        )
        self.assertEqual(
            dc.resolve_convention(canonical_pin, self.AXIS_RECEIPT),
            (dc.CANONICAL_JSON,),
        )
        self.assertNotEqual(
            dc.normalize_digest(raw_pin), dc.normalize_digest(canonical_pin)
        )

    def test_the_tag_appears_on_two_different_conventions(self) -> None:
        census = json.loads(
            (MODULE_DIR / "receipts" / "baryon_dimension_six_census.receipt.json")
            .read_text(encoding="utf-8")
        )
        raw_pin = census["upstream_pins"]["axis_center_receipt"]["sha256"]
        axis_manifest = json.loads(
            (MODULE_DIR / "manifests" / "axis_center_descent_reference.json")
            .read_text(encoding="utf-8")
        )
        stripped_pin = axis_manifest["global_form_artifact_sha256"]
        self.assertTrue(raw_pin.startswith(dc.DIGEST_TAG))
        self.assertTrue(stripped_pin.startswith(dc.DIGEST_TAG))
        self.assertEqual(
            dc.resolve_convention(stripped_pin, ConventionArithmeticTests.SELF_SIGNING_ARTIFACT),
            (dc.CANONICAL_JSON_SELF_DIGEST_STRIPPED,),
        )

    def test_the_self_documenting_manifest_states_its_convention(self) -> None:
        manifest = json.loads(
            (MODULE_DIR / "manifests" / "directed_seam_repair_reference.json")
            .read_text(encoding="utf-8")
        )
        pins = {
            name: pin
            for name, pin in manifest["upstream_pins"].items()
            if "canonical_json_sha256" in pin
        }
        self.assertEqual(
            sorted(pins),
            ["carrier", "integer_record_counting_mechanism", "undirected_scheduler"],
        )
        for name, pin in pins.items():
            with self.subTest(pin=name):
                target = dc.locate(pin["path"], (dc.REPO_ROOT, MODULE_DIR))
                self.assertIsNotNone(target, pin["path"])
                self.assertEqual(
                    dc.resolve_convention(pin["canonical_json_sha256"], target),
                    (dc.CANONICAL_JSON,),
                    "the field name states the convention, so it must hold",
                )


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression and mutation gates for source-current tomography Stage 2."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest


MODULE_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = MODULE_DIR.parents[1]
sys.path.insert(0, str(MODULE_DIR))

import source_current_order_sensitive_inventory as producer  # noqa: E402
import verify_source_current_order_sensitive_inventory as independent  # noqa: E402


class SourceCurrentOrderSensitiveInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.committed = json.loads(producer.INVENTORY_PATH.read_text(encoding="utf-8"))

    def write_mutant(self, value: dict) -> Path:
        value.pop("inventory_sha256", None)
        value["inventory_sha256"] = producer.canonical_sha256(value)
        temporary = tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", encoding="utf-8", delete=False
        )
        with temporary:
            json.dump(value, temporary, indent=2, sort_keys=True)
            temporary.write("\n")
        self.addCleanup(Path(temporary.name).unlink, missing_ok=True)
        return Path(temporary.name)

    def test_committed_inventory_replays_and_verifies_independently(self) -> None:
        rebuilt = producer.build_inventory()
        self.assertEqual(rebuilt, self.committed)
        producer.verify_inventory(self.committed)
        report = independent.verify(producer.INVENTORY_PATH, REPO_ROOT)
        self.assertEqual(report["verdict"], producer.VERDICT)
        self.assertEqual(report["candidate_count"], 21)
        self.assertEqual(report["qualifying_candidate_count"], 0)
        self.assertEqual(report["response_algebra_dimension"], 4)
        self.assertEqual(report["proper_recharting_count"], 60)
        self.assertEqual(report["b14_jacobi_nonzero_count"], 240)

    def test_stop_condition_is_one_packet_and_fail_closed(self) -> None:
        self.assertEqual(
            self.committed["verdict"],
            "SOURCE_CURRENT_ORDER_SENSITIVE_OBJECT_NOT_PRESENT",
        )
        self.assertFalse(self.committed["positive_stages_authorized"])
        self.assertFalse(self.committed["physical_current_source_bridge_attained"])
        self.assertFalse(self.committed["fixture_used_as_reconstruction_oracle"])
        self.assertFalse(self.committed["target_labels_used"])
        self.assertIn(
            "one candidate packet",
            self.committed["audit_scope"]["rule"],
        )
        self.assertIn(
            "different primitives",
            self.committed["summary"]["cross_packet_noncomposition"],
        )
        for row in self.committed["candidates"]:
            self.assertNotEqual(row["classification"], "QUALIFIES")
            self.assertTrue(row["qualification_failures"])

    def test_decisive_near_candidates_remain_separated(self) -> None:
        rows = {row["candidate_id"]: row for row in self.committed["candidates"]}
        response = rows["registered_adjacency_response"]
        self.assertEqual(response["classification"], "COMMUTATIVE_ONLY")
        self.assertTrue(response["properties"]["reversible"])
        self.assertFalse(response["properties"]["raw_histories_serialized"])
        self.assertFalse(response["properties"]["both_composition_orders_recorded"])

        counting = rows["record_counting_repair_628"]
        self.assertEqual(counting["classification"], "IRREVERSIBLE_ONLY")
        self.assertTrue(counting["properties"]["raw_histories_serialized"])
        self.assertTrue(counting["properties"]["refinement_provenance_present"])
        self.assertFalse(counting["properties"]["reversible"])
        self.assertFalse(counting["properties"]["both_composition_orders_recorded"])

        fixture = rows["declared_port_current_fixture"]
        self.assertEqual(fixture["classification"], "DOWNSTREAM_CONTAMINATED")
        self.assertFalse(fixture["properties"]["source_native"])
        self.assertFalse(fixture["properties"]["target_free"])

    def test_minimal_missing_fields_are_explicit(self) -> None:
        self.assertEqual(
            [row["field"] for row in self.committed["minimal_missing_fields"]],
            [
                "perturbation_families",
                "ordered_mixed_histories",
                "raw_history_custody",
                "common_source_binding",
                "reversibility_witness",
                "source_firewall",
            ],
        )
        ordered = self.committed["minimal_missing_fields"][1]["requirement"]
        self.assertIn("66 unordered", ordered)
        self.assertIn("132 ordered", ordered)

    def test_rehashed_false_qualification_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.committed)
        row = next(
            row
            for row in mutant["candidates"]
            if row["candidate_id"] == "registered_adjacency_response"
        )
        row["properties"] = {
            field: True for field in mutant["qualification_fields"]
        }
        row["qualification_failures"] = []
        row["classification"] = "QUALIFIES"
        path = self.write_mutant(mutant)
        with self.assertRaisesRegex(independent.VerificationError, "classification"):
            independent.verify(path, REPO_ROOT)

    def test_rehashed_candidate_omission_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.committed)
        mutant["candidates"] = mutant["candidates"][:-1]
        mutant["summary"]["candidate_count"] -= 1
        path = self.write_mutant(mutant)
        with self.assertRaisesRegex(independent.VerificationError, "candidate coverage"):
            independent.verify(path, REPO_ROOT)

    def test_rehashed_source_pin_doctoring_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.committed)
        mutant["source_pins"][0]["sha256"] = "sha256:" + "0" * 64
        path = self.write_mutant(mutant)
        with self.assertRaisesRegex(independent.VerificationError, "SHA pin drift"):
            independent.verify(path, REPO_ROOT)

    def test_rehashed_positive_stage_promotion_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.committed)
        mutant["positive_stages_authorized"] = True
        mutant["physical_current_source_bridge_attained"] = True
        path = self.write_mutant(mutant)
        with self.assertRaisesRegex(independent.VerificationError, "positive-stage stop"):
            independent.verify(path, REPO_ROOT)

    def test_producer_does_not_import_conditional_current_implementation(self) -> None:
        source = Path(producer.__file__).read_text(encoding="utf-8")
        self.assertNotIn("import port_current_inner_certificate", source)
        self.assertNotIn("from port_current_inner_certificate", source)


if __name__ == "__main__":
    unittest.main()

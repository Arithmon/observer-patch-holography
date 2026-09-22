#!/usr/bin/env python3
"""Regression and mutation gates for source-current tomography Stage 2."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
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

    def copy_required_repository_tree(self, destination: Path) -> None:
        for relative in producer.AUDITED_DIRECTORIES:
            shutil.copytree(REPO_ROOT / relative, destination / relative)
        required_files = {
            row["path"]
            for key in ("source_pins", "implementation_pins")
            for row in self.committed[key]
        }
        for relative in sorted(required_files):
            target = destination / relative
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO_ROOT / relative, target)

    def test_committed_inventory_replays_and_verifies_independently(self) -> None:
        rebuilt = producer.build_inventory()
        self.assertEqual(rebuilt, self.committed)
        producer.verify_inventory(self.committed)
        report = independent.verify(producer.INVENTORY_PATH, REPO_ROOT)
        self.assertEqual(report["verdict"], producer.VERDICT)
        self.assertEqual(report["candidate_count"], 29)
        self.assertEqual(report["qualifying_candidate_count"], 0)
        self.assertEqual(report["response_algebra_dimension"], 4)
        self.assertEqual(report["proper_recharting_count"], 60)
        self.assertEqual(report["b14_jacobi_nonzero_count"], 240)
        self.assertEqual(report["routing_q13_events"], 19113548)
        self.assertEqual(report["routing_q21_events"], 292722053)
        self.assertEqual(
            report["audited_file_count"],
            self.committed["audited_file_snapshot"]["file_count"],
        )

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

        memory = rows["lean_pair_mean_memory_and_reusable_bus"]
        self.assertEqual(memory["classification"], "IRREVERSIBLE_ONLY")
        self.assertTrue(memory["properties"]["source_native"])
        self.assertTrue(memory["properties"]["raw_histories_serialized"])
        self.assertFalse(memory["properties"]["reversible"])
        self.assertFalse(memory["properties"]["same_twelve_port_carrier"])

        selected = rows["lean_selected_pair_mean_histories"]
        self.assertEqual(selected["classification"], "IRREVERSIBLE_ONLY")
        self.assertTrue(selected["properties"]["raw_histories_serialized"])
        self.assertFalse(selected["properties"]["reversible"])
        self.assertFalse(selected["properties"]["both_composition_orders_recorded"])

        routing = rows["source_read_routing_full_family"]
        self.assertEqual(routing["classification"], "IRREVERSIBLE_ONLY")
        self.assertTrue(routing["properties"]["raw_histories_serialized"])
        self.assertTrue(routing["properties"]["refinement_provenance_present"])
        self.assertFalse(routing["properties"]["reversible"])
        self.assertFalse(routing["properties"]["both_composition_orders_recorded"])
        self.assertFalse(routing["properties"]["twelve_reversible_perturbation_families"])

        temporal = rows["native_temporal_tomography_and_checkpoint_selection"]
        self.assertEqual(temporal["classification"], "IRREVERSIBLE_ONLY")
        self.assertFalse(temporal["properties"]["reversible"])
        self.assertFalse(temporal["properties"]["both_composition_orders_recorded"])

        programs = rows["native_stored_programs_and_accumulator"]
        self.assertEqual(programs["classification"], "IRREVERSIBLE_ONLY")
        self.assertFalse(programs["properties"]["reversible"])

        join = rows["finite_source_operator_join"]
        self.assertEqual(join["classification"], "STATIC_ONLY")
        self.assertFalse(join["properties"]["raw_histories_serialized"])

        fermion = rows["fermionic_hypercharge_current_histories"]
        self.assertEqual(fermion["classification"], "DOWNSTREAM_CONTAMINATED")
        self.assertFalse(fermion["properties"]["target_free"])
        self.assertFalse(fermion["properties"]["source_native"])

        maxwell = rows["maxwell_measurement_adapter"]
        self.assertEqual(maxwell["classification"], "DOWNSTREAM_CONTAMINATED")
        self.assertFalse(maxwell["properties"]["raw_histories_serialized"])

        richest = self.committed["summary"]["richest_integrated_near_candidate"]
        self.assertEqual(richest["candidate_id"], "source_read_routing_full_family")

    def test_integrated_tree_review_is_explicit_and_complete(self) -> None:
        review = self.committed["integrated_tree_review"]
        self.assertEqual(review["surface_count"], 111)
        paths = [row["path"] for row in review["surfaces"]]
        self.assertEqual(paths, sorted(set(paths)))
        decisions = {row["decision"] for row in review["surfaces"]}
        self.assertEqual(
            decisions,
            {"NEW_CANDIDATE", "EXTENDS_EXISTING_CANDIDATE"},
        )
        self.assertEqual(
            review["from_upstream_main_sha"],
            "2d9bd11bc47c56d88a2fbbca22e3cc1be171d3f9",
        )
        self.assertEqual(
            review["through_upstream_main_sha"],
            "afed734528edff214c34d4038a64310544df922d",
        )
        candidate_ids = {row["candidate_id"] for row in self.committed["candidates"]}
        for row in review["surfaces"]:
            self.assertIn(row["candidate_id"], candidate_ids)
            self.assertTrue(row["reason"])
        lean_added = {
            path
            for path in paths
            if path.startswith("Lean/Geometry/") and path.endswith(".lean")
        }
        self.assertEqual(len(lean_added), 44)
        self.assertIn("Lean/Geometry/SourceTemporalTomography.lean", lean_added)
        self.assertIn("Lean/Geometry/SourceReadRouting.lean", lean_added)
        prior = review["prior_reviews"]
        self.assertEqual(len(prior), 1)
        self.assertEqual(prior[0]["surface_count"], 40)
        self.assertEqual(
            prior[0]["through_upstream_main_sha"],
            review["from_upstream_main_sha"],
        )
        self.assertIn(
            "Lean/Screen/WhitneyConeMass.lean",
            {row["path"] for row in prior[0]["surfaces"]},
        )

    def test_audit_scope_extension_is_explicit_and_bound(self) -> None:
        extension = self.committed["audit_scope"]["directories_added_at_this_revision"]
        directories = [row["directory"] for row in extension]
        self.assertEqual(directories, list(producer.AUDITED_DIRECTORIES)[13:])
        self.assertIn("code/source_read_routing", directories)
        self.assertIn("evidence/source_net_causal_poset/routed_read_law", directories)
        snapshot_paths = {row["path"] for row in self.committed["audited_file_snapshot"]["files"]}
        for row in extension:
            self.assertTrue(row["reason"])
            self.assertTrue(
                any(path.startswith(row["directory"] + "/") for path in snapshot_paths),
                row["directory"],
            )
        self.assertIn("code/source_read_routing/specification.json", snapshot_paths)
        self.assertIn("code/sm_fermion_current/current_receipt.json", snapshot_paths)

    def test_rehashed_scope_extension_omission_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.committed)
        mutant["audit_scope"]["directories_added_at_this_revision"].pop()
        path = self.write_mutant(mutant)
        with self.assertRaisesRegex(
            independent.VerificationError, "audit scope extension drift"
        ):
            independent.verify(path, REPO_ROOT)

    def test_rehashed_routing_promotion_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.committed)
        row = next(
            row
            for row in mutant["candidates"]
            if row["candidate_id"] == "source_read_routing_full_family"
        )
        row["properties"] = {field: True for field in mutant["qualification_fields"]}
        row["qualification_failures"] = []
        row["classification"] = "QUALIFIES"
        mutant["summary"]["qualifying_candidate_count"] = 1
        path = self.write_mutant(mutant)
        with self.assertRaisesRegex(independent.VerificationError, "classification"):
            independent.verify(path, REPO_ROOT)

    def test_extended_scope_runtime_content_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_required_repository_tree(root)
            target = root / "code/source_read_routing/specification.json"
            original = target.read_bytes()
            changed = original.replace(b"baseline", b"BASELINE", 1)
            self.assertEqual(len(changed), len(original))
            self.assertNotEqual(changed, original)
            target.write_bytes(changed)
            with self.assertRaisesRegex(
                independent.VerificationError, "AUDITED_CONTENT_DRIFT"
            ) as raised:
                independent.verify(producer.INVENTORY_PATH, root)
            self.assertEqual(raised.exception.code, "AUDITED_CONTENT_DRIFT")

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

    def test_audited_file_snapshot_is_canonical_and_current(self) -> None:
        snapshot = self.committed["audited_file_snapshot"]
        self.assertEqual(snapshot["directories"], list(producer.AUDITED_DIRECTORIES))
        paths = [row["path"] for row in snapshot["files"]]
        self.assertEqual(paths, sorted(set(paths)))
        self.assertEqual(snapshot["file_count"], len(snapshot["files"]))
        self.assertEqual(
            snapshot["files_sha256"],
            producer.canonical_sha256(snapshot["files"]),
        )
        self.assertEqual(snapshot["files"], producer.audited_file_records())
        self.assertEqual(
            snapshot["exclusion_policy"]["exact_paths"],
            list(producer.CONTENT_SNAPSHOT_EXACT_EXCLUSIONS),
        )

    def test_new_file_below_audited_directory_makes_inventory_stale(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_required_repository_tree(root)
            added = root / "code/source_feedback_transport/unreviewed_candidate.json"
            added.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(
                independent.VerificationError, "AUDITED_PATH_DRIFT"
            ):
                independent.verify(producer.INVENTORY_PATH, root)

    def test_removed_file_below_audited_directory_makes_inventory_stale(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_required_repository_tree(root)
            (root / "Lean/Dynamics/CenterSpectral.lean").unlink()
            with self.assertRaisesRegex(
                independent.VerificationError, "AUDITED_PATH_DRIFT"
            ):
                independent.verify(producer.INVENTORY_PATH, root)

    def test_existing_audited_file_content_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_required_repository_tree(root)
            unchanged = independent.verify(producer.INVENTORY_PATH, root)
            self.assertEqual(unchanged["verdict"], producer.VERDICT)
            target = root / "Lean/Dynamics/CenterSpectral.lean"
            original = target.read_bytes()
            changed = original.replace(b"theorem", b"Theorem", 1)
            self.assertEqual(len(changed), len(original))
            self.assertNotEqual(changed, original)
            target.write_bytes(changed)
            with self.assertRaisesRegex(
                independent.VerificationError, "AUDITED_CONTENT_DRIFT"
            ) as raised:
                independent.verify(producer.INVENTORY_PATH, root)
            self.assertEqual(raised.exception.code, "AUDITED_CONTENT_DRIFT")

    def test_existing_audited_file_byte_count_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_required_repository_tree(root)
            target = root / "Lean/Dynamics/CenterSpectral.lean"
            target.write_bytes(target.read_bytes() + b"\n")
            with self.assertRaisesRegex(
                independent.VerificationError, "AUDITED_CONTENT_DRIFT"
            ) as raised:
                independent.verify(producer.INVENTORY_PATH, root)
            self.assertEqual(raised.exception.code, "AUDITED_CONTENT_DRIFT")

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

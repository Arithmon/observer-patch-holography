"""Synthetic reference/coverage/freshness tests; no scientific producer runs."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import build_model_theorem_map as builder


class ModelTheoremMapTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "proof.lean").write_text("theorem exact_join : True := trivial\n")
        self.premises = {"rows": [{"id": "PR-01"}, {"id": "PR-85"}]}
        self.claims = {"claims": [
            {"claim_id": "C-SUPPLIED", "assumptions": ["declared_action", "shared"],
             "imported_results": ["spectral_theorem"],
             "premise_dependencies": {"classification": "explicit_non_consumer",
                                      "consumed": [], "open": [], "boundary": []}},
            {"claim_id": "C-BRANCH", "assumptions": ["shared"],
             "imported_results": ["spectral_theorem", "external_scale"],
             "null_decision_assumptions": ["coverage", "shared"],
             "threshold_kinematic_assumptions": ["dispersion"],
             "premise_dependencies": {"classification": "explicit_edges",
                                      "consumed": ["PR-01"], "open": ["PR-01"],
                                      "boundary": ["PR-85"]}},
        ]}
        self.data = {
            "schema": builder.SCHEMA, "as_of": "2026-10-11",
            "coverage_note": "Navigation only; physical identification is excluded.",
            "groups": [{
                "id": "RE-01", "title": "Supplied scalar | branch",
                "links": sorted(builder.MODEL_IDS), "scope": "active",
                "status": "conditional", "premises": ["PR-01", "PR-85"],
                "constructed": "The stated mathematical implication is proved.",
                "remaining": "Supply an implementation on the same domain.",
                "closure": "Bound the implementation error without circular inference.",
                "sources": [{"path": "proof.lean", "label": "Exact implication",
                             "anchor": "theorem exact_join"}],
                "claims": ["C-SUPPLIED"], "assumptions": ["declared_action", "coverage"],
            }],
        }

    def errors(self, data=None):
        return builder.validate(self.data if data is None else data,
                                self.premises, self.claims, self.root)

    def reject(self, mutate, expected):
        data = deepcopy(self.data)
        mutate(data)
        self.assertIn(expected, "\n".join(self.errors(data)))

    def test_complete_fixture_passes_without_changing_claim_edges(self):
        before = deepcopy(self.claims)
        self.assertEqual(self.errors(), [])
        rendered = builder.render(self.data, builder.inventory(self.data, self.premises, self.claims))
        self.assertEqual(self.claims, before)
        self.assertIn("explicit_non_consumer", rendered)
        self.assertIn("no proof-discharge edges", rendered)
        self.assertIn("not a global completeness", rendered.replace("\n", " "))
        self.assertIn("Supplied scalar \\| branch", rendered)

    def test_missing_and_unknown_canonical_ids_fail_separately(self):
        self.reject(lambda d: d["groups"][0]["premises"].pop(), "omitted canonical premises ['PR-85']")
        self.reject(lambda d: d["groups"][0]["premises"].append("PR-404"), "unknown references ['PR-404']")

    def test_missing_and_unknown_model_criteria(self):
        self.reject(lambda d: d["groups"][0]["links"].remove("M12"), "missing model criteria ['M12']")
        self.reject(lambda d: d["groups"][0]["links"].append("M13"), "unknown references ['M13']")

    def test_duplicate_group_id_fails(self):
        self.reject(lambda d: d["groups"].append(deepcopy(d["groups"][0])), "duplicate group id RE-01")

    def test_shared_premises_between_distinct_groups_are_allowed(self):
        other = deepcopy(self.data["groups"][0])
        other["id"] = "RE-02"
        other["scope"] = "downstream"
        self.data["groups"].append(other)
        self.assertEqual(self.errors(), [])
        self.assertEqual(builder.inventory(self.data, self.premises, self.claims)["covered_premises"], 2)

    def test_duplicate_references_within_group_fail(self):
        for field in ("premises", "links", "claims", "assumptions"):
            with self.subTest(field=field):
                self.reject(lambda d: d["groups"][0][field].append(d["groups"][0][field][0]),
                            f"{field}: duplicate entries")

    def test_source_anchor_is_exact_and_case_sensitive(self):
        self.reject(lambda d: d["groups"][0]["sources"][0].update(anchor="Theorem exact_join"),
                    "exact substring absent")
        self.reject(lambda d: d["groups"][0]["sources"][0].update(anchor=" "), "nonempty exact substring")

    def test_source_fragment_must_name_an_explicit_anchor(self):
        (self.root / "page.md").write_text('# Page\n\n<a id="gravity-rungs"></a>\n\n## Ladder\n')
        self.data["groups"][0]["sources"].append({"path": "page.md#gravity-rungs", "label": "Ladder"})
        self.assertEqual(self.errors(), [])
        rendered = builder.render(self.data, builder.inventory(self.data, self.premises, self.claims))
        self.assertIn("- [Ladder](../page.md#gravity-rungs).", rendered)
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="page.md#ladder"),
                    "anchor #ladder absent from page.md")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="page.md#Gravity_Rungs"),
                    "fragment must be lowercase words joined by hyphens")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="page.md#"),
                    "fragment must be lowercase words joined by hyphens")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="absent.md#gravity-rungs"),
                    "missing file")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="../page.md#gravity-rungs"),
                    "inside the repository")

    def test_source_without_anchor_is_allowed_but_missing_file_is_not(self):
        del self.data["groups"][0]["sources"][0]["anchor"]
        self.assertEqual(self.errors(), [])
        self.reject(lambda d: d["groups"][0]["sources"][0].update(path="absent.lean"), "missing file")

    def test_escaping_paths_and_external_symlink_are_rejected(self):
        for path in ("../proof.lean", "/proof.lean", "dir\\proof.lean"):
            with self.subTest(path=path):
                self.reject(lambda d: d["groups"][0]["sources"][0].update(path=path), "inside the repository")
        with tempfile.TemporaryDirectory() as outside:
            external = Path(outside) / "source.lean"
            external.write_text("theorem exact_join")
            (self.root / "link.lean").symlink_to(external)
            self.reject(lambda d: d["groups"][0]["sources"][0].update(path="link.lean"), "outside repository")

    def test_empty_required_text_and_unknown_keys_fail(self):
        for key in ("title", "constructed", "remaining", "closure"):
            with self.subTest(key=key):
                self.reject(lambda d: d["groups"][0].update({key: "\n "}), f"{key}: expected nonempty string")
        self.reject(lambda d: d.update(coverage_note=""), "coverage_note")
        self.reject(lambda d: d["groups"][0].update(proved=True), "unknown keys ['proved']")
        self.reject(lambda d: d.update(discharge_edges=[]), "unknown keys ['discharge_edges']")

    def test_status_scope_and_date_validation(self):
        self.reject(lambda d: d["groups"][0].update(status="closed"), "closed is not supported")
        self.reject(lambda d: d["groups"][0].update(scope="all_proofs"), "scope: expected")
        self.reject(lambda d: d.update(as_of="2026-13-11"), "expected YYYY-MM-DD")
        self.reject(lambda d: d.update(as_of="20261011"), "expected YYYY-MM-DD")

    def test_optional_references_must_resolve(self):
        self.reject(lambda d: d["groups"][0]["claims"].append("invented"), "claims: unknown references")
        self.reject(lambda d: d["groups"][0]["assumptions"].append("invented"), "assumptions: unknown references")

    def test_extra_assumption_fields_are_independently_inventoried(self):
        inv = builder.inventory(self.data, self.premises, self.claims)
        self.assertEqual(inv["ordinary_assumptions"], 2)
        self.assertEqual(inv["imported_names"], 2)
        self.assertEqual(inv["role_counts"], {"consumed": 1, "open": 1, "boundary": 1})
        self.assertEqual(inv["extra_assumption_fields"]["null_decision_assumptions"],
                         {"claims_with_field": 1, "occurrences": 2, "unique_names": 2, "also_ordinary": 1})
        self.assertEqual(inv["extra_assumption_fields"]["threshold_kinematic_assumptions"]["unique_names"], 1)
        text = builder.render(self.data, inv)
        self.assertIn("`threshold_kinematic_assumptions` | 1 | 1 | 1 | 0", text)
        self.assertIn("cannot be added to obtain a certified total", text)

    def test_registry_shape_and_unknown_role_premise_fail_closed(self):
        self.claims["claims"][1]["premise_dependencies"]["open"].append("PR-404")
        self.assertIn("unknown premises ['PR-404']", "\n".join(self.errors()))
        self.claims["claims"][1]["null_decision_assumptions"] = "not a list"
        self.assertIn("null_decision_assumptions: expected a list", "\n".join(self.errors()))

    def test_duplicate_json_keys_are_not_silently_overwritten(self):
        path = self.root / "duplicate.json"
        path.write_text('{"schema": "first", "schema": "second"}')
        with self.assertRaisesRegex(ValueError, "duplicate JSON key: schema"):
            builder.read_json(path)

    def _write_fixture(self, register=None):
        if register is None:
            register = dict(self.premises, **{builder.SECTION: self.data})
        for relative, value in ((builder.REGISTER, register), (builder.CLAIMS, self.claims)):
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(json.dumps(value))

    def test_map_is_the_model_supply_section_of_the_premise_register(self):
        self.assertEqual(builder.REGISTER, Path("tracking/premise_register.json"))
        self.assertEqual(builder.SECTION, "model_supply")
        live = builder.read_json(builder.ROOT / builder.REGISTER)
        self.assertEqual(live[builder.SECTION]["schema"], builder.SCHEMA)
        self.assertEqual(builder.validate(live[builder.SECTION], live,
                                          builder.read_json(builder.ROOT / builder.CLAIMS),
                                          builder.ROOT), [])

    def test_missing_section_fails_closed_without_writing(self):
        self._write_fixture(register=self.premises)
        output = self.root / builder.OUTPUT
        stderr = io.StringIO()
        with patch.object(builder, "ROOT", self.root), redirect_stderr(stderr):
            self.assertEqual(builder.main([]), 1)
        self.assertIn("missing model_supply section", stderr.getvalue())
        self.assertFalse(output.exists())

    def test_cli_generation_and_read_only_stale_check(self):
        self._write_fixture()
        output = self.root / builder.OUTPUT
        with patch.object(builder, "ROOT", self.root), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.assertEqual(builder.main(["--check"]), 1)
            self.assertFalse(output.exists())
            self.assertEqual(builder.main([]), 0)
            before = output.read_bytes()
            mtime = output.stat().st_mtime_ns
            self.assertEqual(builder.main(["--check"]), 0)
            self.assertEqual(output.stat().st_mtime_ns, mtime)
            output.write_text("stale\n")
            self.assertEqual(builder.main(["--check"]), 1)
            self.assertEqual(output.read_text(), "stale\n")
            self.assertEqual(builder.main([]), 0)
            self.assertEqual(output.read_bytes(), before)

    def test_invalid_map_cannot_overwrite_preexisting_output(self):
        self.data["groups"][0]["status"] = "closed"
        self._write_fixture()
        output = self.root / builder.OUTPUT
        output.parent.mkdir(parents=True)
        output.write_text("preserve me\n")
        with patch.object(builder, "ROOT", self.root), redirect_stderr(io.StringIO()):
            self.assertEqual(builder.main([]), 1)
        self.assertEqual(output.read_text(), "preserve me\n")


if __name__ == "__main__":
    unittest.main()

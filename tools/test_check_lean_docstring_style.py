#!/usr/bin/env python3
"""Focused tests for the kernel-axiom coverage scanner in
tools/check_lean_docstring_style.py: name qualification, receipt-line
recognition, suffix matching and the companion registry check."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

import check_lean_docstring_style as checker  # noqa: E402


SOURCE = """
import Mathlib

namespace A.B

theorem x : True := trivial

section Named
variable (n : ℕ)
lemma y : True := trivial
end Named

namespace C
protected theorem z : True := trivial
@[simp] theorem C.w : True := trivial
end C

noncomputable section
theorem v : True := trivial
end

theorem _root_.Top.t : True := trivial
private theorem hidden : True := trivial

-- theorem in_comment : False := sorry
/-- theorem in_docstring : False := sorry -/
def notCounted : ℕ := 0

end A.B

theorem top : True := trivial

#print axioms A.B.x
#print axioms top
"""

AUDIT = """
import Mathlib.Util.AssertNoSorry

open Lean Elab Command in
elab "audit_test_axioms " n:ident : command => do
  logInfo m!"{n}"

/-- error: rejected [sorryAx] -/
#guard_msgs in
audit_test_axioms sorryAx

audit_test_axioms A.B.x
audit_test_axioms B.y
audit_test_axioms A.B.C.z
audit_test_axioms A.B.v
audit_test_axioms Top.t

namespace Other
#print axioms top
end Other
"""


class ScanTests(unittest.TestCase):
    def test_qualified_declarations(self) -> None:
        scan = checker.scan_module(SOURCE)
        self.assertEqual(
            scan.declared,
            ["A.B.x", "A.B.y", "A.B.C.z", "A.B.C.C.w", "A.B.v", "Top.t", "top"],
        )

    def test_receipt_lines_and_forms(self) -> None:
        scan = checker.scan_module(SOURCE)
        self.assertEqual(scan.arguments, {"A.B.x", "top"})
        audit = checker.scan_module(AUDIT)
        self.assertIn("sorryAx", audit.arguments)
        self.assertIn("B.y", audit.arguments)
        forms = dict(audit.entries)["top"]
        self.assertEqual(forms, ("Other.top", "top"))

    def test_final_component_rule_for_in_file_blocks(self) -> None:
        declared, audited = checker.audit_block(SOURCE)
        self.assertEqual(audited, {"x", "top"})
        self.assertEqual(len(declared), 7)
        self.assertEqual(
            checker.unaudited_names(SOURCE), ["y", "z", "w", "v", "t"]
        )

    def test_suffix_matching(self) -> None:
        arguments = {"B.y", "A.B.x"}
        self.assertTrue(checker.names_audited_by(arguments, "A.B.y"))
        self.assertTrue(checker.names_audited_by(arguments, "A.B.x"))
        self.assertFalse(checker.names_audited_by(arguments, "A.B.C.y"))
        self.assertFalse(checker.names_audited_by(arguments, "AB.y"))

    def test_companion_unaudited(self) -> None:
        audit = checker.scan_module(AUDIT)
        source = checker.scan_module(SOURCE)
        self.assertEqual(checker.companion_unaudited(audit, source), ["A.B.C.C.w"])

    def test_companion_registry_flags_omission(self) -> None:
        scans = {
            "Lean/Audit.lean": checker.scan_module(AUDIT),
            "Lean/Source.lean": checker.scan_module(SOURCE),
        }
        registry = {"Lean/Audit.lean": ["Lean/Source.lean", "Lean/Missing.lean"]}
        with mock.patch.object(checker, "AXIOM_AUDIT_COMPANIONS", registry):
            issues = checker.companion_issues(scans)
        self.assertEqual(len(issues), 2)
        self.assertIn("omits 1 public declaration(s) of Lean/Source.lean: A.B.C.C.w", issues[0])
        self.assertIn("Lean/Missing.lean, which is absent", issues[1])

    def test_cross_coverage_prefers_exact_resolution(self) -> None:
        scans = {
            "Lean/Audit.lean": checker.scan_module(AUDIT),
            "Lean/Source.lean": checker.scan_module(SOURCE),
        }
        covered = checker.cross_coverage(scans)
        self.assertEqual(set(covered), {"Lean/Audit.lean"})
        self.assertEqual(
            covered["Lean/Audit.lean"]["Lean/Source.lean"],
            {"A.B.x", "A.B.y", "A.B.C.z", "A.B.v", "Top.t", "top"},
        )


class CorpusTests(unittest.TestCase):
    def test_registered_companion_mappings_are_complete(self) -> None:
        scans = {}
        for path in checker.lean_sources():
            rel = path.relative_to(checker.ROOT).as_posix()
            scans[rel] = checker.scan_module(
                path.read_text(encoding="utf-8", errors="ignore")
            )
        self.assertEqual(checker.companion_issues(scans), [])

    def test_registered_in_file_blocks_are_complete(self) -> None:
        for rel in checker.AXIOM_BLOCK_COMPLETE:
            text = (checker.ROOT / rel).read_text(encoding="utf-8", errors="ignore")
            self.assertEqual(checker.coverage_issues(rel, text), [], rel)


if __name__ == "__main__":
    unittest.main()

"""Semantic mutation checks for finite record construction and completeness."""
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


@pytest.mark.parametrize("module,before,after,test", [
    ("public_record_csp.py",
     "domains[left] = tuple(x for x in domains[left] if read_left[x] == read_right[x])",
     "domains[left] = domains[left]",
     "test_public_record_gluing_regressions.py::test_inconsistent_self_interface_has_no_global_sections"),
    ("public_record_csp.py", "queue = deque(range(len(arcs)))", "return True\n    queue = deque(range(len(arcs)))",
     "test_public_record_sections_independent.py::test_pairwise_feasible_odd_cycle_has_no_global_record"),
    ("public_record_csp.py", "if codomain is not None and", "if False and",
     "test_public_record_gluing_regressions.py::test_readouts_must_land_in_the_declared_interface_atoms"),
    ("public_record_csp.py", 'return value.replace("%", "%25").replace("|", "%7C").replace("=", "%3D")',
     "return value", "test_public_record_gluing_regressions.py::test_distinct_sections_never_share_a_record_identifier"),
    ("public_record_csp.py", 'raise ValueError("section limit exhausted; no complete record set certified")',
     "return sections", "test_public_record_sections_independent.py::test_exhaustion_does_not_return_a_partial_record_set"),
    ("verify_public_record_sections.py", "return actual == _expected(domains, edges, max_candidates)",
     "return actual <= _expected(domains, edges, max_candidates)",
     "test_public_record_sections_independent.py::test_completeness_replay_rejects_forged_or_incomplete_sections"),
    ("verify_public_record_sections.py", "if possible and consistent(assignment):", "if possible:",
     "test_public_record_sections_independent.py::test_spanning_replay_checks_self_interfaces"),
])
def test_broken_record_constructions_are_detected(tmp_path, module, before, after, test):
    source = Path(__file__).parent
    for name in ["public_record_csp.py", "verify_public_record_sections.py",
                 "correctable_public_record_capacity.py", "test_correctable_public_record_capacity.py",
                 "test_public_record_gluing_regressions.py", "test_public_record_sections_independent.py"]:
        shutil.copy2(source/name, tmp_path/name)
    # #1054 adds a numerical helper imported by the capacity evaluator. This
    # copy is optional on main and required when the independent PRs combine.
    if (source/"checkpoint_channels.py").exists():
        shutil.copy2(source/"checkpoint_channels.py", tmp_path)
    target = tmp_path/module
    original = target.read_text(encoding="utf-8")
    assert before in original
    target.write_text(original.replace(before, after), encoding="utf-8")
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", "--tb=short", "-p", "no:cacheprovider", test],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    output = result.stdout+result.stderr
    assert result.returncode == 1, output
    assert "failed" in output and "ERROR" not in output, output
    assert "AssertionError" in output or "DID NOT RAISE" in output, output

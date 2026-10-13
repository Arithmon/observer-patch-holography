"""Fail-closed tests for the #618 gravity premise ladder.

The ladder is the ``gravity_rungs`` section of ``tracking/premise_register.json``
and renders as a section of ``docs/registers/PREMISE_REGISTER_V3.md``.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_gravity_ladder.py"

spec = importlib.util.spec_from_file_location("build_gravity_ladder", BUILDER)
ladder_tool = importlib.util.module_from_spec(spec)
sys.modules["build_gravity_ladder"] = ladder_tool
spec.loader.exec_module(ladder_tool)


def live_ladder() -> dict:
    return ladder_tool.load_ladder()


def test_ladder_is_the_gravity_rungs_section_of_the_premise_register():
    register = json.loads(ladder_tool.REGISTER_PATH.read_text(encoding="utf-8"))
    assert ladder_tool.SECTION == "gravity_rungs"
    assert register[ladder_tool.SECTION] == live_ladder()
    assert ladder_tool.SURFACE_PATH.relative_to(REPO_ROOT).as_posix() == (
        ladder_tool.GENERATED_SURFACE
    )


def test_live_ladder_validates_and_surface_is_current():
    ladder = live_ladder()
    rungs = ladder_tool.validate(ladder)
    rendered = ladder_tool.render(ladder, rungs)
    committed = ladder_tool.SURFACE_PATH.read_text(encoding="utf-8")
    assert committed.count(rendered) == 1
    assert len(rungs) == 11


def test_missing_section_fails(tmp_path, monkeypatch):
    register = json.loads(ladder_tool.REGISTER_PATH.read_text(encoding="utf-8"))
    del register[ladder_tool.SECTION]
    path = tmp_path / "premise_register.json"
    path.write_text(json.dumps(register), encoding="utf-8")
    monkeypatch.setattr(ladder_tool, "REGISTER_PATH", path)
    with pytest.raises(SystemExit, match="no gravity_rungs section"):
        ladder_tool.load_ladder()


def test_stale_generated_surface_fails():
    ladder = live_ladder()
    ladder["generated_surface"] = "docs/registers/GRAVITY_PREMISE_LADDER.md"
    with pytest.raises(SystemExit, match="generated_surface must equal"):
        ladder_tool.validate(ladder)


def test_check_mode_fails_on_a_stale_page(tmp_path, monkeypatch):
    page = ladder_tool.SURFACE_PATH.read_text(encoding="utf-8")
    stale = tmp_path / "PREMISE_REGISTER_V3.md"
    stale.write_text(page.replace("| 11 | Metrological", "| 12 | Metrological", 1), encoding="utf-8")
    monkeypatch.setattr(ladder_tool, "SURFACE_PATH", stale)
    assert ladder_tool.main(["--check"]) == 1
    monkeypatch.undo()
    assert ladder_tool.main(["--check"]) == 0


def test_every_rung_has_countermodel_and_mutation_test():
    rungs = ladder_tool.validate(live_ladder())
    for row in rungs:
        assert row["theorem_or_countermodel"]
        assert row["mutation_tests"]


def test_status_stronger_than_registry_class_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 7)
    row["status"] = "axiom_forced"
    with pytest.raises(SystemExit, match="not among its registry"):
        ladder_tool.validate(ladder)


def test_unknown_interface_fails():
    ladder = live_ladder()
    ladder["rungs"][3]["premise_interfaces"] = ["no_such_interface"]
    with pytest.raises(SystemExit, match="unknown premise interface"):
        ladder_tool.validate(ladder)


def test_missing_lean_declaration_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 9)
    row["theorem_or_countermodel"][0]["declaration"] = "no_such_declaration"
    with pytest.raises(SystemExit, match="absent from"):
        ladder_tool.validate(ladder)


def test_missing_mutation_test_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 5)
    row["mutation_tests"] = []
    with pytest.raises(SystemExit, match="mutation test is required"):
        ladder_tool.validate(ladder)


def test_pending_identification_requires_live_owner():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 10)
    row["owner_issues"] = []
    with pytest.raises(SystemExit, match="requires owner issues"):
        ladder_tool.validate(ladder)


def test_unknown_assumption_token_fails():
    ladder = live_ladder()
    row = next(r for r in ladder["rungs"] if r["rung"] == 4)
    row["assumption_tokens"] = ["no_such_token"]
    with pytest.raises(SystemExit, match="not in the dictionary"):
        ladder_tool.validate(ladder)

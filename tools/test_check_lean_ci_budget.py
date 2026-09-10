from __future__ import annotations

from tools.check_lean_ci_budget import LAKEFILE, WORKFLOW, validate


def _inputs() -> tuple[str, str]:
    return (
        WORKFLOW.read_text(encoding="utf-8"),
        LAKEFILE.read_text(encoding="utf-8"),
    )


def test_current_wiring_passes() -> None:
    workflow, lakefile = _inputs()
    assert validate(workflow, lakefile) == []


def test_default_ophgap_target_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = lakefile.replace(
        "lean_lib «OphGap» where",
        "@[default_target]\nlean_lib «OphGap» where",
        1,
    )
    assert "OphGap must not be a Lake default target" in validate(workflow, mutated)


def test_job_over_thirty_minutes_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace("    timeout-minutes: 30", "    timeout-minutes: 31", 1)
    assert "Lean CI 'build' job exceeds the 30-minute ceiling" in validate(mutated, lakefile)


def test_golden_job_over_thirty_minutes_is_rejected() -> None:
    workflow, lakefile = _inputs()
    marker = (
        "  golden_sector_psl2f5:\n"
        "    name: Golden-sector PSL2F5 representations\n"
        "    runs-on: ubuntu-latest\n"
        "    timeout-minutes: 30"
    )
    mutated = workflow.replace(marker, marker[:-2] + "31", 1)
    assert (
        "Lean CI 'golden_sector_psl2f5' job exceeds the 30-minute ceiling"
        in validate(mutated, lakefile)
    )


def test_full_lake_cache_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace("            Lean/.lake/build", "            Lean/.lake", 1)
    assert (
        "Lean CI must cache build artifacts, not the full Lake dependency tree"
        in validate(mutated, lakefile)
    )


def test_nonresumable_cache_key_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace("-${{ github.run_attempt }}", "")
    errors = validate(mutated, lakefile)
    assert "Lean CI 'build' cache keys must support resumable re-runs" in errors
    assert "Lean CI 'golden_sector_psl2f5' cache keys must support resumable re-runs" in errors
    assert "Lean CI 'ophgap' cache keys must support resumable re-runs" in errors


def test_golden_budget_change_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace(
        "timeout 23m lake build GoldenSectorPSL2F5Representations",
        "timeout 24m lake build GoldenSectorPSL2F5Representations",
        1,
    )
    assert (
        "the typed golden-sector bridge must retain its 23-minute resumable budget"
        in validate(mutated, lakefile)
    )


def test_golden_core_exclusion_is_required() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace(
        "Screen/OPHScreen.lean|Screen/GoldenSectorPSL2F5Representations.lean) continue ;;",
        "Screen/OPHScreen.lean) continue ;;",
        1,
    )
    assert (
        "the core Lean build must leave the import-only OPHScreen umbrella and typed golden bridge to their dedicated/nightly lanes"
        in validate(mutated, lakefile)
    )


def test_duplicate_branch_and_pr_runs_are_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace("    branches: [main]\n", "", 1)
    assert "Lean CI branch pushes must be limited to main" in validate(mutated, lakefile)


def test_missing_superseded_run_cancellation_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace("  cancel-in-progress: true", "  cancel-in-progress: false", 1)
    assert "Lean CI must cancel superseded runs" in validate(mutated, lakefile)


def test_exhaustive_default_build_is_rejected() -> None:
    workflow, lakefile = _inputs()
    mutated = workflow.replace('lake build "$target"', "lake build", 1)
    errors = validate(mutated, lakefile)
    assert "per-change Lean CI must not run an exhaustive default Lake build" in errors

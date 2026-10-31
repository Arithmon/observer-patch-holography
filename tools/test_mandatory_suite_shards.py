"""Exercise the runner's actual argument parsing and subprocess dispatch."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from tools import run_mandatory_suite as runner


def _standard_steps():
    return [step for step in runner.MANDATORY_STEPS if step[0] not in runner.HEAVY_STEP_TITLES]


@pytest.fixture
def dispatch(monkeypatch):
    commands = []

    def capture(command, *, cwd):
        assert cwd == runner.ROOT
        commands.append(list(command))
        return SimpleNamespace(returncode=0)

    # Keep main(), run_steps(), and their ordering/error checks in the path.
    # Only the child scientific programs are replaced by a dispatch recorder.
    monkeypatch.setattr(runner.subprocess, "run", capture)

    def invoke(*arguments):
        commands.clear()
        monkeypatch.setattr(sys, "argv", ["run_mandatory_suite.py", *arguments])
        runner.main()
        return list(commands)

    return invoke, commands


@pytest.mark.parametrize("count", [1, 2, 3, 7, 9, 10, len(_standard_steps())])
def test_all_standard_shards_dispatch_exact_ordered_coverage(dispatch, capsys, count):
    invoke, _ = dispatch
    expected = [command for _, command in _standard_steps()]
    partitions = []
    for index in range(count):
        partitions.append(invoke("--shard-index", str(index), "--shard-count", str(count)))
        output = capsys.readouterr().out
        assert f"standard shard {index + 1}/{count}:" in output
        assert f"mandatory shard OK ({index + 1}/{count}; standard mode)" in output
        assert "mandatory suite OK" not in output
    assert all(partitions)
    # Cost-weighted cut points: an optimal contiguous partition never exceeds
    # the even share plus the largest single step.
    titles = [title for title, _ in _standard_steps()]
    costs, position = [], 0
    for partition in partitions:
        costs.append(sum(runner.step_cost(t) for t in titles[position:position + len(partition)]))
        position += len(partition)
    largest = max(runner.step_cost(t) for t in titles)
    assert max(costs) <= sum(costs) / count + largest
    # Exact concatenation checks both contiguous order and multiplicity.
    assert [command for partition in partitions for command in partition] == expected


def test_duplicate_titles_and_commands_survive_sharding(monkeypatch, dispatch):
    steps = [
        ("duplicate", ["command", "same"]),
        ("duplicate", ["command", "different"]),
        ("heavy", ["excluded", "only-from-standard"]),
        ("duplicate", ["command", "same"]),
        ("last", ["command", "last"]),
    ]
    monkeypatch.setattr(runner, "MANDATORY_STEPS", steps)
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset({"heavy"}))
    invoke, _ = dispatch
    first = invoke("--shard-index", "0", "--shard-count", "2")
    second = invoke("--shard-index", "1", "--shard-count", "2")
    assert first == [["command", "same"], ["command", "different"]]
    assert second == [["command", "same"], ["command", "last"]]


@pytest.mark.parametrize(
    ("arguments", "mode"),
    [
        ([], "standard"),
        (["--full"], "full"),
        (["--certificates"], "standard-certificates"),
        (["--full", "--certificates"], "full-certificates"),
        (["--certificates-only"], "certificates"),
        (["--certificate-smoke-only"], "smoke"),
    ],
)
def test_unsharded_modes_keep_their_complete_command_lists(dispatch, capsys, arguments, mode):
    expected = {
        "standard": _standard_steps(),
        "full": runner.MANDATORY_STEPS,
        "standard-certificates": _standard_steps() + runner.CERTIFICATE_STEPS,
        "full-certificates": runner.MANDATORY_STEPS + runner.CERTIFICATE_STEPS,
        "certificates": runner.CERTIFICATE_STEPS,
        "smoke": runner.CERTIFICATE_SMOKE_STEPS,
    }[mode]
    invoke, _ = dispatch
    assert invoke(*arguments) == [command for _, command in expected]
    assert "shard OK" not in capsys.readouterr().out


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--shard-index", "0"], "must be supplied together"),
        (["--shard-count", "2"], "must be supplied together"),
        (["--shard-index", "-1", "--shard-count", "2"], "0 <= index"),
        (["--shard-index", "2", "--shard-count", "2"], "0 <= index"),
        (["--shard-index", "0", "--shard-count", "0"], "shards must be nonempty"),
        (["--shard-index", "0", "--shard-count", "-1"], "shards must be nonempty"),
        (["--shard-index", "0", "--shard-count", str(len(_standard_steps()) + 1)], "shards must be nonempty"),
        (["--shard-index", "0.5", "--shard-count", "2"], "invalid int value"),
        (["--shard-index", "0", "--shard-count", "two"], "invalid int value"),
        (["--shard-index", "0", "--shard-count", "2", "--full"], "only for the standard suite"),
        (["--shard-index", "0", "--shard-count", "2", "--certificates"], "only for the standard suite"),
        (["--shard-index", "0", "--shard-count", "2", "--certificates-only"], "only for the standard suite"),
        (["--shard-index", "0", "--shard-count", "2", "--certificate-smoke-only"], "only for the standard suite"),
    ],
)
def test_invalid_shard_arguments_fail_before_any_child_runs(dispatch, capsys, arguments, message):
    invoke, commands = dispatch
    with pytest.raises(SystemExit) as error:
        invoke(*arguments)
    assert error.value.code == 2
    assert message in capsys.readouterr().err
    assert commands == []


def test_empty_standard_list_cannot_report_a_successful_shard(monkeypatch, dispatch, capsys):
    monkeypatch.setattr(runner, "MANDATORY_STEPS", [])
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset())
    invoke, commands = dispatch
    with pytest.raises(SystemExit) as error:
        invoke("--shard-index", "0", "--shard-count", "1")
    assert error.value.code == 2
    assert "shards must be nonempty" in capsys.readouterr().err
    assert commands == []


def test_sharding_keeps_unknown_heavy_step_guard(monkeypatch, dispatch):
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset({"missing step"}))
    invoke, commands = dispatch
    with pytest.raises(SystemExit, match="HEAVY_STEP_TITLES entries missing"):
        invoke("--shard-index", "0", "--shard-count", "2")
    assert commands == []


def test_failed_child_stops_its_shard_and_cannot_report_success(monkeypatch, capsys):
    monkeypatch.setattr(runner, "MANDATORY_STEPS", [
        ("first shard", ["unused"]),
        ("first shard again", ["unused-too"]),
        ("failing child", ["fails"]),
        ("must not run", ["later"]),
    ])
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset())
    commands = []

    def failure(command, *, cwd):
        assert cwd == runner.ROOT
        commands.append(command)
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(runner.subprocess, "run", failure)
    monkeypatch.setattr(sys, "argv", [
        "run_mandatory_suite.py", "--shard-index", "1", "--shard-count", "2",
    ])
    with pytest.raises(SystemExit, match=r"FAILED: failing child \(exit 7\)"):
        runner.main()
    assert commands == [["fails"]]
    assert "shard OK" not in capsys.readouterr().out


def test_step_costs_name_existing_steps():
    titles = {title for title, _ in runner.MANDATORY_STEPS}
    assert set(runner.STEP_COST_SECONDS) <= titles
    assert all(cost > 0 for cost in runner.STEP_COST_SECONDS.values())


def test_shard_bounds_minimize_the_most_expensive_contiguous_shard(monkeypatch):
    monkeypatch.setattr(runner, "STEP_COST_SECONDS", {"a": 9, "b": 1, "c": 1, "d": 8, "e": 1})
    steps = [(title, ["command", title]) for title in "abcde"]
    bounds = runner.shard_bounds(steps, 3)
    assert bounds == [0, 1, 3, 5]  # [a] [b c] [d e]: maximum cost 9, the single largest step
    assert runner.shard_bounds(steps, 5) == [0, 1, 2, 3, 4, 5]
    assert runner.shard_bounds(steps, 1) == [0, 5]

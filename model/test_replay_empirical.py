"""All canonical empirical checks must run, with ordinary failure semantics."""
import subprocess
import sys
from types import SimpleNamespace

import pytest

from . import replay_empirical as worker
from . import verify


def test_fixed_worker_invokes_each_canonical_check_and_restores_argv(monkeypatch):
    original = ["parent", "untouched"]
    monkeypatch.setattr(sys, "argv", original)
    imported, checked = [], []

    def load(name):
        imported.append(name)
        def main():
            checked.append(list(sys.argv))
        return SimpleNamespace(main=main)

    monkeypatch.setattr(worker.importlib, "import_module", load)
    worker.main([])
    assert imported == ["model.primitive_observer.__main__",
                        "model.anharmonic_transfer.__main__",
                        "model.anharmonic_observations.__main__"]
    assert checked == [
        ["model.primitive_observer", "--check",
         str(worker.ROOT / "model/primitive_observer/results.json")],
        ["model.anharmonic_transfer", "--check",
         str(worker.ROOT / "model/anharmonic_transfer/results.json")],
        ["model.anharmonic_observations", "--check",
         str(worker.ROOT / "model/anharmonic_observations/results.json")],
    ]
    assert sys.argv is original


@pytest.mark.parametrize("failed_index", [0, 1, 2])
@pytest.mark.parametrize("failure", [ValueError("replay mismatch"), SystemExit(7)])
def test_failure_stops_checks_and_restores_argv(monkeypatch, failed_index, failure):
    original = ["parent"]
    monkeypatch.setattr(sys, "argv", original)
    calls = []

    def load(name):
        index = len(calls)
        calls.append(name)
        def main():
            if index == failed_index:
                raise failure
        return SimpleNamespace(main=main)

    monkeypatch.setattr(worker.importlib, "import_module", load)
    with pytest.raises(type(failure)) as raised:
        worker.main([])
    assert raised.value is failure
    assert len(calls) == failed_index+1
    assert sys.argv is original


@pytest.mark.parametrize("arguments", [
    ["--module", "unreviewed"], ["--check", "other.json"], ["--write", "replacement.json"],
])
def test_worker_rejects_alternate_modules_receipts_and_writes(monkeypatch, arguments):
    monkeypatch.setattr(worker.importlib, "import_module",
                        lambda _: pytest.fail("unsupported arguments reached a producer"))
    with pytest.raises(SystemExit) as raised:
        worker.main(arguments)
    assert raised.value.code == 2


def test_verifier_preserves_every_replay_in_order_when_grouping(monkeypatch):
    # Compare the expanded child commands with the entire original declaration,
    # so grouping cannot silently drop an empirical check or its neighbors.
    monkeypatch.setattr(sys, "argv", ["verify"])
    monkeypatch.setattr(verify, "verify_artifacts", lambda: 1)
    expanded, children = [], []

    def run(argv, **kwargs):
        children.append(argv)
        assert kwargs["check"] is True
        assert kwargs["cwd"] == verify.ROOT
        if argv == [sys.executable, "-m", "model.replay_empirical"]:
            expanded.extend(worker.EMPIRICAL_REPLAYS)
        else:
            check = argv.index("--check")
            command = None if check == 3 else argv[3]
            expanded.append((argv[2], command, argv[check+1]))

    monkeypatch.setattr(verify.subprocess, "run", run)
    verify.main()
    assert tuple(expanded) == verify.REPLAYS
    assert len(children) == len(verify.REPLAYS)-len(worker.EMPIRICAL_REPLAYS)+1
    assert sum(argv[-1] == "model.replay_empirical" for argv in children) == 1


def test_worker_failure_stops_the_outer_verifier(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["verify"])
    monkeypatch.setattr(verify, "verify_artifacts", lambda: 1)
    children = []

    def run(argv, **kwargs):
        children.append(argv)
        if argv[-1] == "model.replay_empirical":
            raise subprocess.CalledProcessError(1, argv)

    monkeypatch.setattr(verify.subprocess, "run", run)
    with pytest.raises(subprocess.CalledProcessError):
        verify.main()
    assert len(children) == verify.REPLAYS.index(worker.EMPIRICAL_REPLAYS[0])+1
    assert children[-1][-1] == "model.replay_empirical"

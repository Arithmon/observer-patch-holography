"""Receipt parsing checks independent of the expensive source execution."""
import json
from pathlib import Path
import sys

import pytest

from tools.strict_json import DuplicateKeyError, load
from . import __main__ as cli


@pytest.fixture
def replayed_receipt(monkeypatch):
    receipt = load(Path(__file__).with_name("results.json"))
    # The source histories are tested separately. Keep this test focused on
    # the actual CLI parser and comparison of a complete retained receipt.
    monkeypatch.setattr(cli.experiment, "run", lambda *args: receipt)
    return receipt


def test_cli_accepts_unambiguous_matching_receipt(replayed_receipt, tmp_path, monkeypatch, capsys):
    path = tmp_path / "receipt.json"
    path.write_text(cli.experiment.canonical(replayed_receipt), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["source_scalar", "--check", str(path)])
    cli.main()
    assert json.loads(capsys.readouterr().out)["all_pass"] is True


@pytest.mark.parametrize("location", ("root", "branch"))
def test_cli_rejects_duplicate_keys_even_when_last_value_matches(
        replayed_receipt, location, tmp_path, monkeypatch):
    text = cli.experiment.canonical(replayed_receipt)
    if location == "root":
        original = '"all_pass": true'
        duplicate = '"all_pass": false, "all_pass": true'
    else:
        probability = replayed_receipt["levels"][0]["step_refinement"][0]["histories"][0]["branches"][0]["probability"]
        original = f'"probability": {probability}'
        duplicate = f'"probability": -1, {original}'
    assert original in text
    text = text.replace(original, duplicate, 1)
    # The permissive decoder would hide the invalid first value, allowing
    # the complete receipt to compare equal to the fresh execution.
    assert json.loads(text) == replayed_receipt
    path = tmp_path / "ambiguous.json"
    path.write_text(text, encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["source_scalar", "--check", str(path)])
    with pytest.raises(DuplicateKeyError, match="duplicate JSON key"):
        cli.main()

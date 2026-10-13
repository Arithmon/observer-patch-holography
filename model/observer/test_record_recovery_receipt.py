"""Bind complete recovery evidence, source bytes and scope; reject trash."""
from copy import deepcopy
from pathlib import Path
import subprocess
import sys

import pytest

from . import record_recovery_receipt as receipt


@pytest.fixture(scope="module")
def actual():
    return receipt.run()


def leaves(value, path=()):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from leaves(child, path + (key,))
    elif isinstance(value, list):
        for key, child in enumerate(value):
            yield from leaves(child, path + (key,))
    else:
        yield path, value


def test_retained_receipt_replays(actual):
    saved = Path(__file__).parent / "results" / "record_recovery.json"
    receipt.compare_receipt(receipt.read_receipt(saved), actual)


def test_every_exported_value_is_binding(actual):
    count = 0
    for path, value in leaves(actual):
        changed = deepcopy(actual)
        parent = changed
        for key in path[:-1]:
            parent = parent[key]
        parent[path[-1]] = not value if type(value) is bool else "corrupt"
        with pytest.raises(ValueError, match="receipt mismatch"):
            receipt.compare_receipt(changed, actual)
        count += 1
    assert count > 100


@pytest.mark.parametrize("trash", (None, True, 1, "pass", [], {}, {"all_pass": True}))
def test_trash_is_rejected(actual, trash):
    with pytest.raises(ValueError):
        receipt.compare_receipt(trash, actual)


def test_missing_sections_extra_keys_and_wrong_scalar_types_fail(actual):
    for key in actual:
        changed = deepcopy(actual)
        del changed[key]
        with pytest.raises(ValueError):
            receipt.compare_receipt(changed, actual)
    for value in (1, 1.0, "true", None):
        changed = deepcopy(actual)
        changed["scope"]["one_event_encoded_qubit_recoverable"] = value
        with pytest.raises(ValueError):
            receipt.compare_receipt(changed, actual)
    changed = deepcopy(actual)
    changed["extra"] = 0
    with pytest.raises(ValueError):
        receipt.compare_receipt(changed, actual)


@pytest.mark.parametrize("text", ('{"x":1,"x":2}', '{"x":{"a":1,"a":2}}',
                                  '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}'))
def test_duplicate_keys_and_nonfinite_values_fail(tmp_path, text):
    path = tmp_path / "bad.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        receipt.read_receipt(path)


def test_changed_source_digest_cannot_replay(actual, monkeypatch):
    monkeypatch.setattr(receipt, "SOURCE_PATHS", ("model/observer/record_recovery_receipt.py",))
    with pytest.raises(ValueError):
        receipt.compare_receipt(actual, receipt.run())


def test_cli_replays_rejects_trash_and_preserves_existing_files(tmp_path, actual):
    path = tmp_path / "receipt.json"
    path.write_text(receipt.canonical(actual), encoding="utf-8")
    command = [sys.executable, "-m", "model.observer.record_recovery_receipt"]
    ok = subprocess.run(command + ["--check", str(path)], cwd=receipt.ROOT,
                        capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    original = path.read_bytes()
    refused = subprocess.run(command + ["--write", str(path)], cwd=receipt.ROOT,
                             capture_output=True, text=True)
    assert refused.returncode != 0
    assert path.read_bytes() == original
    path.write_text('{"all_pass":true}', encoding="utf-8")
    bad = subprocess.run(command + ["--check", str(path)], cwd=receipt.ROOT,
                         capture_output=True, text=True)
    assert bad.returncode != 0
    assert "receipt mismatch" in bad.stderr

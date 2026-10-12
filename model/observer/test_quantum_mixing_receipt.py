"""Exact replay must reject a wrong theorem, missing sector or altered scope."""
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch
import subprocess
import sys

import pytest

from . import quantum_mixing_receipt as receipt


@pytest.fixture(scope="module")
def actual():
    return receipt.run()


def _leaves(value, path=()):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from _leaves(child, path + (key,))
    elif isinstance(value, list):
        for key, child in enumerate(value):
            yield from _leaves(child, path + (key,))
    else:
        yield path, value


def test_retained_exact_receipt(actual):
    saved = Path(__file__).parent / "results" / "quantum_mixing.json"
    receipt.compare_receipt(receipt.read_receipt(saved), actual)


def test_every_exported_scalar_is_binding(actual):
    mutations = 0
    for path, value in _leaves(actual):
        changed = deepcopy(actual)
        parent = changed
        for key in path[:-1]:
            parent = parent[key]
        if type(value) is bool:
            replacement = not value
        elif type(value) is int:
            replacement = value + 1
        else:
            replacement = "corrupted"
        parent[path[-1]] = replacement
        with pytest.raises(ValueError, match="receipt mismatch"):
            receipt.compare_receipt(changed, actual)
        mutations += 1
    assert mutations > 50


@pytest.mark.parametrize("trash", [None, True, 1, "pass", [], {}, {"all_pass": True}])
def test_trash_never_passes(actual, trash):
    with pytest.raises(ValueError, match="receipt mismatch"):
        receipt.compare_receipt(trash, actual)


def test_structure_types_and_empty_sections_are_binding(actual):
    for key in actual:
        changed = deepcopy(actual)
        del changed[key]
        with pytest.raises(ValueError):
            receipt.compare_receipt(changed, actual)
    for value in (1, 1.0, "true", None):
        changed = deepcopy(actual)
        changed["scope"]["binary_eventual_entanglement_breaking"] = value
        with pytest.raises(ValueError):
            receipt.compare_receipt(changed, actual)
    changed = deepcopy(actual)
    changed["sectors"] = {}
    with pytest.raises(ValueError):
        receipt.compare_receipt(changed, actual)


@pytest.mark.parametrize("content", [
    '{"scope": 1, "scope": 2}',
    '{"nested": {"gap": 1, "gap": 2}}',
    '{"gap": NaN}', '{"gap": Infinity}', '{"gap": -Infinity}',
])
def test_invalid_json_values_and_duplicates(tmp_path, content):
    target = tmp_path / "invalid.json"
    target.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError):
        receipt.read_receipt(target)


def test_nonfinite_and_extra_keys_rejected(actual):
    for value in (float("nan"), float("inf"), 0):
        changed = deepcopy(actual)
        changed["extra"] = value
        with pytest.raises(ValueError):
            receipt.compare_receipt(changed, actual)


def test_spectrum_cannot_pass_with_wrong_incidence():
    cycle = tuple((i, (i + 1) % 12) for i in range(12))
    with patch.object(receipt.quantum_mixing, "seams", return_value=cycle):
        with pytest.raises(ArithmeticError, match="wrong exact spectrum"):
            receipt.carrier_spectrum()


def test_cli_checks_and_refuses_overwrite(tmp_path, actual):
    target = tmp_path / "receipt.json"
    target.write_text(receipt.canonical(actual), encoding="utf-8")
    command = [sys.executable, "-m", "model.observer.quantum_mixing_receipt"]
    ok = subprocess.run(command + ["--check", str(target)], capture_output=True,
                        text=True, cwd=receipt.ROOT)
    assert ok.returncode == 0, ok.stderr
    original = target.read_bytes()
    refused = subprocess.run(command + ["--write", str(target)], capture_output=True,
                             text=True, cwd=receipt.ROOT)
    assert refused.returncode != 0
    assert target.read_bytes() == original
    target.write_text('{"all_pass": true}', encoding="utf-8")
    bad = subprocess.run(command + ["--check", str(target)], capture_output=True,
                         text=True, cwd=receipt.ROOT)
    assert bad.returncode != 0
    assert "receipt mismatch" in bad.stderr

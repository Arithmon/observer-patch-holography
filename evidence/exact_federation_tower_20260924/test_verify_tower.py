"""The tower package verifier passes on the package and fails on mutations."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


def run(package: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(package / "verify_tower.py")], capture_output=True, text=True)


def copy_package(tmp_path: Path) -> Path:
    dst = tmp_path / "package"
    shutil.copytree(HERE, dst, ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache"))
    return dst


def rehash(package: Path) -> None:
    import hashlib

    manifest = json.loads((package / "manifest.json").read_text())
    for rel in manifest["files"]:
        manifest["files"][rel] = hashlib.sha256((package / rel).read_bytes()).hexdigest()
    (package / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def test_package_passes() -> None:
    result = run(HERE)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "TOWER_PACKAGE_PASS" in result.stdout


def test_manifest_detects_edit(tmp_path) -> None:
    package = copy_package(tmp_path)
    level = sorted(package.glob("L*"))[0]
    receipt = level / "receipt.json"
    receipt.write_text(receipt.read_text() + "\n")
    result = run(package)
    assert result.returncode == 1 and "manifest digest differs" in result.stdout


@pytest.mark.parametrize("field", ["V_terminal", "quotient_hash", "attempts_to_balanced_class", "ledger"])
def test_mutated_receipt_fails(tmp_path, field) -> None:
    package = copy_package(tmp_path)
    level = sorted(package.glob("L*"))[0]
    receipt = json.loads((level / "receipt.json").read_text())
    entry = receipt["integer_law"]["entries"][0]
    if field == "ledger":
        entry["V_ledger"][1] = entry["V_ledger"][0] + 1
    elif field == "quotient_hash":
        entry["quotient_hash"] = "0" * 64
    else:
        entry[field] += 1
    (level / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    rehash(package)
    result = run(package)
    assert result.returncode == 1 and "FAIL" in result.stdout, result.stdout


def test_kernel_mutation_fails(tmp_path) -> None:
    package = copy_package(tmp_path)
    level = sorted(package.glob("L*"))[0]
    receipt = json.loads((level / "receipt.json").read_text())
    cell = receipt["response_kernels"]["cells"][0]
    n = next(iter(cell["kernels"]))
    cell["kernels"][n][0][1] += 1e-3
    cell["kernels"][n][1][0] += 1e-3
    (level / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    rehash(package)
    result = run(package)
    assert result.returncode == 1 and "slow-band shares" in result.stdout


def test_missing_host_record_fails(tmp_path) -> None:
    package = copy_package(tmp_path)
    level = sorted(package.glob("L*"))[0]
    for p in level.glob("verify_*.json"):
        p.unlink()
    manifest = json.loads((package / "manifest.json").read_text())
    manifest["files"] = {k: v for k, v in manifest["files"].items() if "/verify_" not in k}
    (package / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    result = run(package)
    assert result.returncode == 1 and "host verification" in result.stdout

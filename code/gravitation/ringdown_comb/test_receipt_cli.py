"""Actual CLI replay rejects altered evidence, even with a matching fresh hash."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
RECEIPT = "integer_k_comb_template_receipt.json"


def run_cli(directory: Path, name: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-W", "error", str(directory / name)],
        cwd=directory, capture_output=True, text=True, timeout=120,
    )


@pytest.fixture
def isolated_instrument(tmp_path):
    for name in ("integer_k_comb_template.py", "verify_integer_k_comb_independent.py"):
        shutil.copyfile(HERE / name, tmp_path / name)
    (tmp_path / "runtime").mkdir()
    shutil.copyfile(HERE / "runtime" / RECEIPT, tmp_path / "runtime" / RECEIPT)
    return tmp_path


def test_producer_cli_rebuilds_canonical_receipt_and_verifier_replays(isolated_instrument):
    directory = isolated_instrument
    path = directory / "runtime" / RECEIPT
    expected = path.read_bytes()
    path.write_bytes(b"stale evidence\n")
    produced = run_cli(directory, "integer_k_comb_template.py")
    assert produced.returncode == 0, produced.stdout + produced.stderr
    assert path.read_bytes() == expected
    assert hashlib.sha256(expected).hexdigest() in produced.stdout
    # A successful replay must not depend on importing the producer.
    (directory / "integer_k_comb_template.py").unlink()
    verified = run_cli(directory, "verify_integer_k_comb_independent.py")
    assert verified.returncode == 0, verified.stdout + verified.stderr
    assert "VERIFIED" in verified.stdout


@pytest.mark.parametrize("mutation", [
    "empty", "truncated", "noncanonical", "duplicate_key", "numeric_field",
    "scientific_boundary", "integer_as_bool", "wrong_tooth",
])
def test_verifier_cli_rejects_altered_receipts(isolated_instrument, mutation):
    directory = isolated_instrument
    path = directory / "runtime" / RECEIPT
    original = path.read_bytes()
    receipt = json.loads(original)
    if mutation == "empty":
        payload = b""
    elif mutation == "truncated":
        payload = original[:-7]
    elif mutation == "noncanonical":
        payload = json.dumps(receipt, indent=2).encode("ascii") + b"\n"
    elif mutation == "duplicate_key":
        payload = b'{"schema":"fabricated",' + original[1:]
    else:
        if mutation == "numeric_field":
            receipt["reference_point_synthetic"]["kappa_per_s_sig40"] = "0E+0"
        elif mutation == "scientific_boundary":
            receipt["boundary"] = "Registered physical prediction."
        elif mutation == "integer_as_bool":
            receipt["reference_point_synthetic"]["m_azimuthal"] = True
        else:
            receipt["universal_ladder"][0]["x_sig40"] = "1.000E+0"
        payload = json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode("ascii") + b"\n"
    assert payload != original
    path.write_bytes(payload)
    verified = run_cli(directory, "verify_integer_k_comb_independent.py")
    assert verified.returncode == 1, verified.stdout + verified.stderr
    assert "MISMATCH" in verified.stdout and "VERIFIED" not in verified.stdout
    # A self-consistent digest does not supply the independent numerical replay.
    assert hashlib.sha256(payload).hexdigest() in verified.stdout
    assert path.read_bytes() == payload

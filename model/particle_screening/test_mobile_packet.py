"""Complete receipt custody and hostile-input tests for the mobile observer."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from tools import strict_json
from . import mobile_packet as packet


@pytest.fixture(scope="module")
def result():
    return packet.report()


def test_fresh_complete_report_matches_retained_receipt(result):
    packet.compare(strict_json.load(packet.OUTPUT), result)
    assert result["exact"]["stationary_reset_probability_per_attempt"] == "1/150"
    assert result["exact"]["raw_census"]["C"]["raw"] == 14400
    assert result["exact"]["raw_census"]["E"]["raw"] == 28800
    assert result["physical_particle_identification"] is False
    assert result["integrated_gate_accepted"] is False
    for name, digest in result["source_sha256"].items():
        assert hashlib.sha256((packet.ROOT / name).read_bytes()).hexdigest() == digest


@pytest.mark.parametrize("damage", ["reset", "count", "scope", "boolean", "custody", "missing"])
def test_comparison_rejects_forged_science_costs_scope_and_custody(result, damage):
    value = deepcopy(result)
    if damage == "reset":
        value["exact"]["stationary_reset_probability_per_attempt"] = "0"
    elif damage == "count":
        value["execution"]["resources"]["laws"]["C"]["maximum_recovery_pulses"] = 0
    elif damage == "scope":
        value["physical_particle_identification"] = True
    elif damage == "boolean":
        value["integrated_gate_accepted"] = 0
    elif damage == "custody":
        value["source_sha256"]["model/particle_screening/mobile.py"] = "0" * 64
    else:
        del value["execution"]["source_replay"]
    with pytest.raises(ValueError, match="differs from complete source replay"):
        packet.compare(value, result)


def _diagnostics_from_receipt(value):
    if type(value) is dict:
        if value.get("kind") == "fresh_numerical_admission":
            return 0.0
        return {key: _diagnostics_from_receipt(v) for key, v in value.items()}
    if type(value) is list:
        return [_diagnostics_from_receipt(v) for v in value]
    return value


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1e-30, 0.01, False, 0])
def test_fresh_gate_rejects_invalid_residuals_instead_of_erasing_them(result, bad):
    diagnostic = _diagnostics_from_receipt(result["execution"])
    diagnostic["source_composition"]["maximum_encoder_operator_error"] = bad
    with pytest.raises(ArithmeticError, match="outside fixed gate"):
        packet.admit(diagnostic)


def test_fixed_gates_require_every_fresh_measurement_and_no_unregistered_float(result):
    diagnostic = _diagnostics_from_receipt(result["execution"])
    assert packet.admit(diagnostic) == result["execution"]
    del diagnostic["source_composition"]["maximum_decoder_operator_error"]
    with pytest.raises(ArithmeticError, match="missing fresh"):
        packet.admit(diagnostic)
    diagnostic = _diagnostics_from_receipt(result["execution"])
    diagnostic["unregistered_error"] = 0.0
    with pytest.raises(ArithmeticError, match="unregistered"):
        packet.admit(diagnostic)


@pytest.mark.parametrize("text", ['{"x":0,"x":1}', '{"x":NaN}'])
def test_receipt_consumer_refuses_duplicate_or_nonfinite_data(tmp_path, text):
    path = tmp_path / "bad.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        packet.compare(strict_json.load(path), {})


def test_real_cli_replays_explicit_receipt_path_and_keeps_bytes(tmp_path, result):
    path = tmp_path / "mobile.json"
    original = json.dumps(result, indent=2)
    path.write_text(original, encoding="utf-8")
    run = subprocess.run(
        [sys.executable, "-m", "model.particle_screening.mobile_packet", "--check", str(path)],
        cwd=packet.ROOT, capture_output=True, text=True, timeout=240,
    )
    assert run.returncode == 0, run.stderr
    assert "complete native flags" in run.stdout
    assert path.read_text(encoding="utf-8") == original

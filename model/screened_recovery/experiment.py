"""Source-bound receipt for the screened q=1 active recovery service."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from ..particle_screening import screening
from . import protocol


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "screened_recovery_results.json"
PREDECESSOR_FILES = (
    "model/particle_screening/screening.py",
    "model/particle_screening/screening_results.json",
    "model/common_source/record_code.py",
    "model/common_source/native_routing.py",
    "model/observer/record_recovery.py",
    "model/observer/source_discrimination/circuit.py",
    "model/observer/source_discrimination/experiment.py",
    "model/observer/source_scalar/compiler.py",
    "model/exchange_sector/sector.py",
)
PACKET_FILES = (
    "model/screened_recovery/__init__.py",
    "model/screened_recovery/protocol.py",
    "model/screened_recovery/experiment.py",
    "model/screened_recovery/test_protocol.py",
    "model/screened_recovery/README.md",
    "model/screened_recovery/campaign.json",
)


def _sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def science_result():
    """Machine-stable exact outputs; numerical diagnostics live beside them."""
    return {
        "schema": "oph.screened_recovery.science.v1",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "service": protocol.exact_service_certificate(),
        "compiled_layout": protocol.compiled_layout_certificate(),
        "resources": protocol.resource_report(),
        "cp03_accepted": False,
        "autonomous_particle_identified": False,
        "rest_energy_identified": False,
        "physical_lifetime_identified": False,
    }


def campaign_binding(science):
    campaign = strict_json.load(HERE/"campaign.json")
    expected = {
        "schema": "oph.screened_recovery.campaign.v1",
        "campaign_id": "screened-recovery-20261028-v1",
        "predecessor_git_commit": "7cef80f58f51f06e70ccda3d4373344ca3993dc6",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "whole_program_frozen": False,
        "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
        "integrated_gate_accepted": False,
        "native_family_changed": False,
        "physical_comparison_count_added": 0,
    }
    for key, value in expected.items():
        if type(campaign.get(key)) is not type(value) or campaign[key] != value:
            raise ValueError("screened-recovery campaign scope changed")
    source = campaign.get("predecessor_source_sha256")
    if (type(source) is not dict or set(source) != set(PREDECESSOR_FILES)
            or any(type(digest) is not str or digest != _sha(path)
                   for path, digest in source.items())):
        raise ValueError("screened-recovery predecessor source changed")
    result_hash = hashlib.sha256(_canonical(science)).hexdigest()
    if campaign.get("science_result_sha256") != result_hash:
        raise ValueError("screened-recovery exact science result changed")
    return {"campaign_id": campaign["campaign_id"],
            "campaign_sha256": _sha("model/screened_recovery/campaign.json"),
            "science_result_sha256": result_hash}


def source_bindings():
    bindings = screening.source_bindings()
    for path in PREDECESSOR_FILES+PACKET_FILES:
        bindings[path] = _sha(path)
    bindings["model/particle_screening/screening_results.json"] = _sha(
        "model/particle_screening/screening_results.json")
    return dict(sorted(bindings.items()))


def run():
    science = science_result()
    return {
        **science,
        "source_replay": protocol.source_replay_report(),
        "source_composition": protocol.source_composition_report(),
        **campaign_binding(science),
        "source_sha256": source_bindings(),
    }


def _compare(expected, actual, path="root"):
    if type(expected) is not type(actual):
        raise ValueError(f"receipt type mismatch at {path}")
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise ValueError(f"receipt key mismatch at {path}")
        for key in expected:
            _compare(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f"receipt length mismatch at {path}")
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _compare(left, right, f"{path}[{i}]")
    elif isinstance(expected, float):
        if not math.isclose(expected, actual, rel_tol=1e-6, abs_tol=1e-12):
            raise ValueError(f"receipt numerical diagnostic mismatch at {path}")
    elif expected != actual:
        raise ValueError(f"receipt mismatch at {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    options.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    elif args.check:
        _compare(strict_json.load(args.check), result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

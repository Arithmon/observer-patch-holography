"""Shared source-bound diagnostic contract; no physical gate promotion."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REQUIRED_SOURCES = frozenset((
    "model/completion/campaign.json", "model/completion/observer.py",
    "model/completion/family.py",
    "model/native_modular_reference/history_transport.py",
    "model/native_modular_reference/history_clock.py",
    "model/native_modular_reference/local_preparation.py",
    "model/observer/adaptive_record_recovery.py", "model/observer/native.py",
    "model/observer/quantum_mixing.py", "model/observer/record_recovery.py",
    "model/observer_time_fields/spatial_dynamics.py",
    "model/relative_observer_geometry/intrinsic_geometry.py",
))


def source_bindings(paths):
    result = {}
    for name in paths:
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError("existing repository source required")
        result[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(result.items()))


def load_contract(path=None):
    value = strict_json.load(HERE / "campaign.json" if path is None else path)
    if value["schema"] != "oph.clock_defect_energy.diagnostic_contract.v1":
        raise ValueError("unknown diagnostic contract")
    for name in ("integrated_gate_accepted", "whole_program_frozen", "native_family_changed"):
        if value[name] is not False:
            raise ValueError("diagnostic contract cannot promote a physical or source claim")
    if (type(value["physical_comparison_count_added"]) is not int
            or value["physical_comparison_count_added"] != 0):
        raise ValueError("this diagnostic has no empirical comparison")
    if (value["scope"] != "FINITE_COMPONENT_DIAGNOSTICS"
            or value["boundary_status"] != "FROZEN_DIAGNOSTIC_CONTRACT"
            or value["predecessor"]["merged_into_main_claimed"] is not False):
        raise ValueError("diagnostic scope or pending dependency changed")
    if type(value["source_sha256"]) is not dict or set(value["source_sha256"]) != REQUIRED_SOURCES:
        raise ValueError("complete predecessor source bindings required")
    if source_bindings(value["source_sha256"]) != value["source_sha256"]:
        raise ValueError("predecessor source changed; review a new diagnostic contract")
    return value


def check_receipt(path, expected):
    """Reject duplicate keys and type-changing edits (including false for zero)."""
    actual = strict_json.load(path)
    canonical = lambda value: json.dumps(value, sort_keys=True, allow_nan=False)
    if canonical(actual) != canonical(expected):
        raise ValueError("diagnostic receipt differs from exact replay")


def binding():
    value = load_contract()
    return {"campaign_id": value["campaign_id"],
            "contract_sha256": hashlib.sha256((HERE / "campaign.json").read_bytes()).hexdigest(),
            "predecessor_git_commit": value["predecessor"]["git_commit"]}

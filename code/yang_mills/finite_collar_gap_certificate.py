#!/usr/bin/env python3
"""Replay exact declared-table arithmetic for a finite calibration family.

This interface cannot certify a physical source or continuum construction.
It shares row validation with the theorem-witness interface; neither validates
the physical meaning or completeness of the supplied source descriptions.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_gap_contract import canonical, evaluate, integer, load, name, rational, require_scope

SCHEMA = "oph.yang_mills.collar_gap_certificate.v1"


def expand_types(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand independent self-repeating types, preserving every influence."""
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be an object")
    if ("type_table" in manifest) == ("calibration_type_family" in manifest):
        raise ValueError("supply exactly one type_table or calibration_type_family")
    if "type_table" in manifest:
        if not isinstance(manifest["type_table"], list):
            raise ValueError("type_table must be a list")
        return copy.deepcopy(manifest["type_table"])
    family = manifest["calibration_type_family"]
    if not isinstance(family, dict):
        raise ValueError("calibration_type_family must be an object")
    count = integer(family.get("count"), minimum=1)
    prefix = name(family.get("id_prefix"))
    template = family.get("template")
    if not isinstance(template, dict) or {"id", "refinement_targets"} & template.keys():
        raise ValueError("family template must omit the generated id and refinement_targets")
    influences = template.get("influences")
    if (not isinstance(influences, list)
            or any(not isinstance(item, dict) or "target_type" in item for item in influences)):
        raise ValueError("family influences must be a list omitting the generated target_type")
    types = []
    for index in range(count):
        entry = copy.deepcopy(template)
        entry["id"] = f"{prefix}{index:03d}"
        entry["refinement_targets"] = [entry["id"]]
        for influence in entry["influences"]:
            influence["target_type"] = entry["id"]
        types.append(entry)
    return types


def validate(manifest: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ValueError("unsupported manifest schema")
    require_scope(manifest, "calibration")
    types = expand_types(manifest)
    bounds = evaluate(types)
    return {
        "schema": SCHEMA,
        "scope": manifest["scope"],
        "physical_clay_receipt": False,
        "active_type_count": len(types),
        **bounds,
        "type_table_sha256": hashlib.sha256(canonical(types)).hexdigest(),
        "type_hashes": {e["id"]: hashlib.sha256(canonical(e)).hexdigest() for e in types},
    }


def certify(manifest_path: Path, output_path: Path) -> None:
    manifest = load(manifest_path)
    receipt = validate(manifest)
    receipt["manifest_sha256"] = hashlib.sha256(canonical(manifest)).hexdigest()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")


def verify(manifest_path: Path, receipt_path: Path) -> None:
    manifest = load(manifest_path)
    expected = load(receipt_path)
    actual = validate(manifest)
    actual["manifest_sha256"] = hashlib.sha256(canonical(manifest)).hexdigest()
    # Serialization distinguishes JSON false/0 and true/1.
    if canonical(actual) != canonical(expected):
        raise ValueError("receipt does not exactly recompute from its manifest")


def main() -> int:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("certify", "verify"):
        command = commands.add_parser(name)
        command.add_argument("--manifest", type=Path, required=True)
        command.add_argument("--output" if name == "certify" else "--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "certify":
        certify(args.manifest, args.output)
    else:
        verify(args.manifest, args.receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

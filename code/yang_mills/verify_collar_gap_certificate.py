#!/usr/bin/env python3
"""Check declared finite row arithmetic in the issue-306 theorem witness.

The supplied source descriptions and rates remain declarations. This checker
cannot promote them to verified source, transfer or continuum evidence.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_gap_contract import evaluate, load, rational as fraction, require_scope


def verify(payload: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, dict) or payload.get("schema") != "oph.yang_mills.collar_gap.v1":
        raise ValueError("unsupported certificate schema")
    require_scope(payload, "theorem_contract_witness")
    bounds = evaluate(payload.get("types"))
    expected = payload.get("expected", {})
    allowed = {"c_floor", "eta_upper", "gap_lower"}
    if not isinstance(expected, dict) or expected.keys() - allowed:
        raise ValueError("expected must contain only c_floor, eta_upper or gap_lower")
    for key, value in expected.items():
        if fraction(value) != fraction(bounds[key]):
            raise ValueError(f"expected {key}={value}, computed {bounds[key]}")
    return {"valid": True, "scope": payload["scope"], "physical_clay_receipt": False,
            "type_count": len(payload["types"]), **bounds}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    payload = load(args.certificate)
    print(json.dumps(verify(payload), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

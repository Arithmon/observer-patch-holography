"""Exact accounting of declared finite collar rows, not source/continuum proof.

Both legacy interfaces use this implementation. Independence comes from the
original-input controls, not from presenting two copies of one formula as two
validators. Descriptive source fields do not establish completeness of kernels,
rates, the refinement family or a physical transfer/continuum construction.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.strict_json import load as _load


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def load(path: Path) -> dict:
    value = _load(path)
    canonical(value)  # Reject non-JSON NaN/Infinity even in descriptive fields.
    if not isinstance(value, dict):
        raise ValueError("certificate must be a JSON object")
    return value


def rational(value: Any) -> Fraction:
    """Read original integer/string evidence without numeric coercion."""
    if type(value) not in (int, str):
        raise ValueError("expected integer or rational string; boolean/floating inputs are forbidden")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid rational {value!r}") from exc


def integer(value: Any, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"expected an exact integer >= {minimum}")
    return value


def name(value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("source names and descriptions must be nonempty strings")
    return value


def names(value: Any, *, nonempty: bool = False) -> list[str]:
    if not isinstance(value, list) or (nonempty and not value):
        raise ValueError("expected a list of source names")
    return [name(item) for item in value]


def require_scope(payload: dict, allowed: str) -> None:
    if payload.get("scope") != allowed:
        raise ValueError(
            f"scope must be {allowed}; physical_source_receipt is unsupported: "
            "these checkers do not verify physical source or continuum evidence"
        )


def evaluate(types: Any) -> dict:
    """Validate every supplied row and return exact declared-table bounds."""
    if not isinstance(types, list) or not types or any(not isinstance(e, dict) for e in types):
        raise ValueError("active types must be a nonempty list of source objects")
    ids = [name(e.get("id")) for e in types]
    if len(set(ids)) != len(ids):
        raise ValueError("collar type ids must be unique")
    known = set(ids)
    rates = []
    row_sums = {}
    required = {"id", "rooted_neighborhood", "boundary_sector_flags", "local_alphabets",
                "potential_templates", "complete_repaired_readback", "conditional_kernel",
                "rate_lower", "influences", "refinement_targets"}
    for entry in types:
        if missing := required - entry.keys():
            raise ValueError(f"{entry['id']}: incomplete source signature: {sorted(missing)}")
        for key in ("rooted_neighborhood", "complete_repaired_readback", "conditional_kernel"):
            name(entry[key])
        for key in ("boundary_sector_flags", "potential_templates", "local_alphabets"):
            names(entry[key], nonempty=key == "local_alphabets")
        targets = names(entry["refinement_targets"], nonempty=True)
        if any(target not in known for target in targets):
            raise ValueError(f"{entry['id']}: refinement leaves the finite active type table")
        rate = rational(entry["rate_lower"])
        if rate <= 0:
            raise ValueError(f"{entry['id']}: rate lower bound must be positive")
        rates.append(rate)
        influences = entry["influences"]
        if not isinstance(influences, list):
            raise ValueError(f"{entry['id']}: influences must be a list")
        total = Fraction(0)
        for item in influences:
            if not isinstance(item, dict) or name(item.get("target_type")) not in known:
                raise ValueError(f"{entry['id']}: influence needs a known target type")
            bound = rational(item.get("upper"))
            multiplicity = integer(item.get("multiplicity", 1))
            if bound < 0:
                raise ValueError(f"{entry['id']}: influence bound must be nonnegative")
            rows = item.get("conditional_rows")
            if (not isinstance(rows, list) or len(rows) != 2
                    or any(not isinstance(row, list) or not row for row in rows)
                    or len(rows[0]) != len(rows[1])):
                raise ValueError(f"{entry['id']}: influence needs two equal-length nonempty rows")
            left, right = ([rational(value) for value in row] for row in rows)
            if min(left + right) < 0 or sum(left) != 1 or sum(right) != 1:
                raise ValueError(f"{entry['id']}: conditional rows must be probability laws")
            tv = sum(abs(a-b) for a, b in zip(left, right, strict=True)) / 2
            if tv > bound:
                raise ValueError(f"{entry['id']}: stated influence is below exact TV distance")
            total += multiplicity * bound
        row_sums[entry["id"]] = total
    c_floor, eta = min(rates), max(row_sums.values())
    if eta >= 1:
        raise ValueError(f"Dobrushin upper bound must be < 1, got {eta}")
    return {"c_floor": str(c_floor), "eta_upper": str(eta),
            "approximate_tensorization_upper": str(1 / (1-eta)),
            "gap_lower": str(c_floor * (1-eta)),
            "row_sums": {key: str(value) for key, value in sorted(row_sums.items())}}

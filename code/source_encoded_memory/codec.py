"""Strict serialization and input custody; contains no execution logic."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPPORT = "code/source_routing/support_w12_l3.json"
PINS = (
    SUPPORT,
    "Lean/lean-toolchain", "Lean/lakefile.lean", "Lean/lake-manifest.json",
    "Lean/ObserverPatchHolography/ScalarSeamRepair.lean",
    "Lean/Geometry/SourceEncodedMemory.lean",
    "Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean",
    *(f"code/source_encoded_memory/{p}" for p in
      ("__init__.py", "codec.py", "build.py", "verify.py")),
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate key: {key}")
        result[key] = value
    return result


def reject(value):
    raise ValueError(f"nonintegral JSON number: {value}")


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=unique, parse_float=reject,
                      parse_constant=reject)


def rational(value):
    if not isinstance(value, str):
        raise ValueError("rational must be a string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as e:
        raise ValueError("invalid rational") from e
    if str(result) != value:
        raise ValueError("noncanonical rational")
    return result


def equal(actual, expected, label):
    if canonical(actual) != canonical(expected):
        raise ValueError(f"{label} mismatch")


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS}

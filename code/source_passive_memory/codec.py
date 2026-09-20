"""Serialization only; builder and verifier have separate state machines."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPPORT = "code/source_routing/support_w12_l3.json"
PROOFS = ("SourceNonlinearRecord", "SourcePassiveReset", "SourcePassiveMemoryBudget")
PIN_PATHS = (
    SUPPORT,
    "Lean/lean-toolchain", "Lean/lakefile.lean", "Lean/lake-manifest.json",
    "Lean/ObserverPatchHolography/ScalarSeamRepair.lean",
    "Lean/Geometry/SourceRecordProtection.lean",
    *(f"Lean/Geometry/{name}.lean" for name in PROOFS),
    "Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean",
    "code/source_passive_memory/__init__.py",
    "code/source_passive_memory/codec.py",
    "code/source_passive_memory/build.py",
    "code/source_passive_memory/verify.py",
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def reject_number(value):
    raise ValueError(f"non-integral JSON number: {value}")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate key: {key}")
        result[key] = value
    return result


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=unique_object,
                      parse_float=reject_number, parse_constant=reject_number)


def rational(value):
    if not isinstance(value, str):
        raise ValueError("rational must be a canonical string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("invalid rational") from error
    if str(result) != value:
        raise ValueError("noncanonical rational")
    return result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_pins():
    return {path: digest(ROOT/path) for path in PIN_PATHS}


def equal(actual, expected, label):
    # Type-sensitive: Python's True == 1 is not accepted as a receipt match.
    if canonical(actual) != canonical(expected):
        raise ValueError(f"{label} mismatch")

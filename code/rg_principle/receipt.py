"""Regenerate compact evidence; verification replays every circuit input."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "rg_principle"

from .certificates import algebra_certificate, clock_certificate, resource_certificate
from .checker import check_program
from .compiler import cases

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "code/rg_principle/receipt.json"
SOURCES = [f"code/rg_principle/{name}" for name in
           ("compiler.py", "checker.py", "certificates.py", "receipt.py", "DERIVATION.md",
            "CONTRACT.md", "test_rg_principle.py")]
SOURCES += ["Lean/Geometry/RecordGluingPrinciple.lean", "Lean/Geometry/RecordGluingPrincipleAxiomAudit.lean"]


def compute():
    return {"schema": "oph-rg-principle-v1", "sources": {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCES},
        "circuits": {name: check_program(program, name in ("retained_copy", "clean_conjunction"))
                     for name, program in cases()},
        "clock": clock_certificate(), "algebra": algebra_certificate(),
        "resources": resource_certificate()}


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key}")
        result[key] = value
    return result


def verify(path):
    submitted = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=no_duplicates,
                           parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
    actual = compute()
    # Exact canonical encoding also separates bool from int and float from int.
    if json.dumps(submitted, sort_keys=True, allow_nan=False) != json.dumps(actual, sort_keys=True, allow_nan=False):
        raise ValueError("receipt differs from source-bound independent replay")
    return actual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "verify"))
    parser.add_argument("--path", type=Path, default=EVIDENCE)
    args = parser.parse_args()
    if args.action == "build":
        result = compute()
        header = json.dumps({k: v for k, v in result.items() if k != "circuits"}, indent=2, sort_keys=True)
        rows = ["    " + json.dumps(name) + ": " + json.dumps(row, sort_keys=True)
                for name, row in sorted(result["circuits"].items())]
        rendered = header[:-2] + ',\n  "circuits": {\n' + ",\n".join(rows) + "\n  }\n}\n"
        args.path.write_text(rendered, encoding="utf-8", newline="\n")
    else:
        result = verify(args.path)
    print(f"{args.action}: {len(result['circuits'])} circuits; "
          f"{sum(c['executed_gates'] for c in result['circuits'].values())} replayed gates; "
          f"{sum(c['interventions'] for c in result['circuits'].values())} localized interventions")


if __name__ == "__main__":
    main()

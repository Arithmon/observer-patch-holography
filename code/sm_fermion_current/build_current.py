"""Produce exact finite fermionic hypercharge-current evidence."""
import argparse
import hashlib
from pathlib import Path

import current

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "current_receipt.json"
PINS = (
    "code/sm_fermion_current/current.py",
    "code/sm_fermion_current/build_current.py",
    "code/sm_fermion_current/verify_current.py",
    "code/sm_fermion_current/test_current.py",
    "code/sm_fermion_current/README.md",
    "code/sm_local_action/jet_action.py",
    "code/sm_abelian_reduction/pilot.py",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/A5FamilyBand.lean",
)


def produce():
    result = current.produce()
    result["source_pins"] = {path: {"bytes": len((ROOT/path).read_bytes()),
                                     "sha256": hashlib.sha256((ROOT/path).read_bytes()).hexdigest()}
                             for path in PINS}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = produce()
    data = current.canonical(result)
    if len(data) > 90_000_000:
        raise ValueError("receipt exceeds bounded local artifact budget")
    if args.write:
        args.output.write_bytes(data)
        print(f"wrote {args.output}: {len(data)} bytes")
    else:
        print(f"computed {len(data)} bytes; use --write to retain receipt")
    print("events", [len(run["events"]) for run in result["runs"]])


if __name__ == "__main__":
    main()

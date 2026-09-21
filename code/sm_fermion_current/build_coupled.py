"""Build the bounded exact coupled semiclassical transport receipt."""
import argparse
import hashlib
from pathlib import Path

import coupled
import current

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "coupled_receipt.json"
PINS = (
    "code/sm_fermion_current/coupled.py",
    "code/sm_fermion_current/build_coupled.py",
    "code/sm_fermion_current/verify_coupled.py",
    "code/sm_fermion_current/test_coupled.py",
    "code/sm_fermion_current/coupled_README.md",
    "code/sm_fermion_current/current.py",
    "code/sm_fermion_current/verify_current.py",
    "code/sm_fermion_current/current_receipt.json",
    "code/sm_local_action/jet_action.py",
    "code/sm_abelian_reduction/pilot.py",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/A5FamilyBand.lean",
)


def produce():
    result = coupled.produce()
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
    raw = current.canonical(result)
    if len(raw) > 20_000_000:
        raise ValueError("declared receipt byte budget exceeded")
    if args.write:
        args.output.write_bytes(raw)
        print(f"wrote {args.output}: {len(raw)} bytes")
    else:
        print(f"computed {len(raw)} bytes; use --write to retain")
    print("prefix+continuation events", [len(run["prefix"]["events"])+len(run["events"])
                                         for run in result["runs"]])


if __name__ == "__main__":
    main()

"""Build the exact quantized-link fermion current receipt."""
import argparse
import hashlib
from pathlib import Path

import current
import quantum_link

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "quantum_link_receipt.json"
PINS = (
    "code/sm_fermion_current/quantum_link.py",
    "code/sm_fermion_current/build_quantum_link.py",
    "code/sm_fermion_current/verify_quantum_link.py",
    "code/sm_fermion_current/test_quantum_link.py",
    "code/sm_fermion_current/README.md",
    "code/sm_fermion_current/current.py",
    "code/sm_fermion_current/verify_current.py",
    "code/sm_local_action/jet_action.py",
    "code/sm_abelian_reduction/pilot.py",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/A5FamilyBand.lean",
)


def produce():
    result = quantum_link.produce()
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
    if len(data) > 10_000_000:
        raise ValueError("quantum-link receipt exceeds bounded artifact budget")
    if args.write:
        args.output.write_bytes(data)
        print(f"wrote {args.output}: {len(data)} bytes")
    else:
        print(f"computed {len(data)} bytes; use --write to retain receipt")
    print("branches", result["final_state"]["branch_census"]["branch_count"])
    print("mean electric", result["observables"]["mean_electric_change"])


if __name__ == "__main__":
    main()

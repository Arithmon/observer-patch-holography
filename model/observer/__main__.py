"""Run finite observer calculations or seeded sequential trials."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import capacity, native, simulate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("native", "capacity", "run"):
        subparser = commands.add_parser(command)
        if command == "run":
            subparser.add_argument("--trials", type=int, default=20_000)
            subparser.add_argument("--seed", type=int, default=20261012)
        receipts = subparser.add_mutually_exclusive_group()
        receipts.add_argument("--write", type=Path, help="create a new receipt; never overwrite")
        receipts.add_argument("--check", type=Path, help="require byte-for-byte canonical replay")
    args = parser.parse_args()
    module = {"native": native, "capacity": capacity, "run": simulate}[args.command]
    result = module.run(args.trials, args.seed) if args.command == "run" else module.run()
    encoded = module.canonical(result)
    if args.write:
        with args.write.open("x") as output:
            output.write(encoded)
    if args.check and args.check.read_text() != encoded:
        raise SystemExit("Receipt mismatch")
    if args.command == "native":
        summary = {"all_pass": result["all_pass"], "checks": len(result["checks"]),
                   "histories": len(result["main_histories"])}
    elif args.command == "capacity":
        summary = {"M0_binary": result["families"]["binary_66"]["reachable_terminal_records"],
                   "M0_full": result["families"]["full_78"]["reachable_terminal_records"]}
    else:
        summary = {"trials": result["trials"],
                   "reference_matches": result["execution"]["exact_path_reference_matches"],
                   "brier": {key: value["sample_decimal"] for key, value in result["brier"].items()}}
    summary["result_sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()

"""Replay or write the exact record-law classification receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from tools import strict_json

from . import family


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def main():
    parser = argparse.ArgumentParser(description=family.__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", type=Path, help="replay the retained receipt exactly")
    mode.add_argument("--write", type=Path, help="write a new receipt, refusing overwrite")
    args = parser.parse_args()
    result = family.run()
    if args.check:
        if canonical(strict_json.load(args.check)) != canonical(result):
            raise ValueError("the record-law receipt does not replay")
        print("Admissible record laws, ladder, completion and interference forms replayed.")
        return
    payload = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n"
    if args.write:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()

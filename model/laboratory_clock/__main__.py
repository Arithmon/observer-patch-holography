"""Replay the clock attachment or check a supplied endpoint contract."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from tools import strict_json

from .attachment import analyse
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", type=Path, help="recompute and compare a retained receipt")
    modes.add_argument("--write", type=Path, help="write a new receipt without overwriting")
    modes.add_argument("--data", type=Path, help="test exact endpoint bounds from a supplied JSON file")
    args = parser.parse_args()
    result = analyse(strict_json.load(args.data)) if args.data else run()
    if args.check:
        if strict_json.load(args.check) != result:
            raise ValueError("the laboratory-clock attachment receipt does not replay")
        print("Affine and operation-duration bounds, native histories and host controls replayed.")
        return
    payload = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n"
    if args.write:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()

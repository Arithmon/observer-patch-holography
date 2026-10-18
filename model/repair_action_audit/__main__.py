"""Replay or retain the native repair action audit."""
import argparse
import json
from pathlib import Path

from tools import strict_json
from .construction import canonical, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", type=Path)
    mode.add_argument("--write", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        if canonical(strict_json.load(args.check)) != canonical(result):
            raise ValueError("repair action audit receipt differs")
        print("Native repair action classification and phase controls replayed.")
        return
    payload = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n"
    if args.write:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()

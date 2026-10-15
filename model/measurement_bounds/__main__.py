"""Generate or replay complete native measurement bounds."""
import argparse
import json
from pathlib import Path

from tools import strict_json
from ..occupation_control.__main__ import compare
from .experiment import produce


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = produce()
    if args.check:
        expected = strict_json.load(args.check)
        compare(expected, expected)
        compare(result, expected)
        print('Complete measurement bounds and fresh operator certificates replayed.')
    elif args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
        print(f'Wrote {args.output}')
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()

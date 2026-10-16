"""Replay a finite prepared sequential native-law discrimination protocol."""
import argparse
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical
from .construction import compare, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare(strict_json.load(args.check), result)
        print('Finite sequential preparation and native prediction replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()

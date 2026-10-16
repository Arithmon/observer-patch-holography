"""Replay source unitary probes and the two-order public-data constraint check."""
import argparse
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical
from .experiment import run, compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare(strict_json.load(args.check), result)
        print('Source unitary probes and two-order characteristic constraints replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()

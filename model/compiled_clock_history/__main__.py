"""Reproduce the complete source-compiled observer history."""
import argparse
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical
from .receipt import compare_receipt, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', type=Path)
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(result, strict_json.loads(args.check.read_text(encoding='utf-8')))
        print('Complete source-compiled clock history reproduced.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()

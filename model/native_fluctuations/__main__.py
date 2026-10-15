"""Replay the native fluctuation action and signed-current predictions."""
import argparse
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if result.get('all_pass') is not True:
        raise ArithmeticError('native fluctuation construction failed')
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
        print('Native fluctuation receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()

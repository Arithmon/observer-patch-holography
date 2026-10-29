"""Recompute signed regional CAR maps and complete source-history joins."""
import argparse
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical
from ..interacting_field.replay import compare_field_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if not result['all_pass']:
        raise ArithmeticError('regional field checks failed: '+str(result['checks']))
    # The strict replay format requires native JSON scalar types in memory.
    compare_field_receipt(strict_json.loads(canonical(result)), result)
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_field_receipt(strict_json.load(args.check), result)
    print(json.dumps({'all_pass': result['all_pass'],
                     'checks': result['checks'],
                     'nested_cube_comparison': result['nested_cube_comparison']}, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

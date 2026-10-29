"""Recompute the complete field-taste receipt."""
import argparse
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical,compare_receipt
from .experiment import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path)
    mode.add_argument('--check',type=Path)
    args=parser.parse_args(); result=run()
    if not result['all_pass']:
        raise ArithmeticError(str(result['checks']))
    if args.write:
        with args.write.open('x',encoding='utf8') as out:
            out.write(canonical(result))
    if args.check:
        compare_receipt(strict_json.load(args.check),result)
    print(json.dumps({'all_pass':result['all_pass'],'checks':result['checks']},indent=2))


if __name__=='__main__':
    main()

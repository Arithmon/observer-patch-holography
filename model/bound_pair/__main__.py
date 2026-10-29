"""Produce or independently recompute the bound-pair receipt."""
import argparse
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical,compare_receipt
from .experiment import run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path)
    mode.add_argument('--check',type=Path)
    args=p.parse_args()
    result=run()
    if args.write:
        with args.write.open('x',encoding='utf8') as f:
            f.write(canonical(result))
    if args.check:
        compare_receipt(strict_json.load(args.check),result)
    print(json.dumps({k:result[k] for k in ('scope','mass','integrated_gate_accepted')},indent=2))


if __name__=='__main__':
    main()

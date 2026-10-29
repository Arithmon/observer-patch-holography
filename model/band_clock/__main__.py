"""Recompute the complete-field record-rate receipt."""
import argparse
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical
from .replay import compare_band_receipt
from .experiment import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path)
    mode.add_argument('--check',type=Path)
    args=parser.parse_args(); result=run()
    if args.write:
        args.write.write_text(canonical(result))
    elif args.check:
        compare_band_receipt(strict_json.load(args.check),result)
    else:
        print(canonical(result))


if __name__=='__main__':
    main()

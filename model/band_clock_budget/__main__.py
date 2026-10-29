"""Replay the counted full-field schedule without changing its premises."""
import argparse
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(canonical(result))
    else:
        compare_receipt(strict_json.load(args.check), result)
        print('Complete reduced full-field schedule and source custody verified.')


if __name__ == '__main__':
    main()

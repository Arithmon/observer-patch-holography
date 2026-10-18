"""Execute complete source identification or verify its retained receipt."""
import argparse
import json
from pathlib import Path

from tools import strict_json
from .experiment import canonical, compare_receipt, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', type=Path)
    modes.add_argument('--write', type=Path)
    args = parser.parse_args()
    expected = strict_json.load(args.check) if args.check else None
    result = run()
    if args.check:
        compare_receipt(expected, result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    print(json.dumps({'all_pass': result['all_pass'], 'rank': result['rank']['rank'],
                      'source_probes': 36, 'continuing_histories': len(result['continuing_histories']),
                      'physical_repair_law_selected': False}, indent=2))


if __name__ == '__main__':
    main()

"""Replay the common interacting source and its complete observer records."""
import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path, help='create a receipt, refusing overwrite')
    group.add_argument('--check', type=Path, help='replay exact fields and floating tolerances')
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding='utf-8')), result)
    print(json.dumps({'all_pass': result['all_pass'], 'levels': [
        {key: row[key] for key in ('dimension', 'scalar_feedback_response',
                                   'disconnected_scalar_feedback_response')}
        for row in result['levels']], 'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()

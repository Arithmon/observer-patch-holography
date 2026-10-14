"""Execute or replay the observer-record reconstruction comparison."""
import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding='utf-8')), result)
    print(json.dumps({'all_pass': result['all_pass'],
        'levels': [{'occupation': row['occupation_per_normal_mode'],
                    'mass': row['recovered']['mass'],
                    'coupling': row['recovered']['coupling'],
                    'spatial_dimension': row['recovered']['spatial_dimension'],
                    'withheld_probability_error': row['maximum_withheld_probability_error']}
                   for row in result['levels']], 'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()

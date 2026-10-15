"""Execute or replay the interacting observer-record recovery benchmark."""
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
        'quartic_g': result['recovered']['quartic'],
        'bare_mass_squared': result['recovered']['bare_mass_squared'],
        'levels': [{'tensor_dimension': row['tensor_dimension'],
                    'prediction_tv': max(row['prediction_history_tv_by_phase']),
                    'ideal_prediction_bound': row['ideal_inferred_to_source_history_tv_upper'],
                    'quadratic_deletion_tv': min(row['quadratic_deletion_history_tv_by_phase'])}
                   for row in result['levels']], 'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()

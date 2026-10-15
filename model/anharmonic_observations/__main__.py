"""Execute and score the fixed retrospective anharmonic observation protocol."""
import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical
from .experiment import predict, evaluate
from .replay import compare_receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    parser.add_argument('--predictions-write', type=Path,
                        help='save the complete prediction payload before scoring (exclusive write)')
    args = parser.parse_args()
    prediction = predict(progress=print)
    if args.predictions_write:
        with args.predictions_write.open('x') as handle:
            handle.write(canonical(prediction))
    print('Predictions fixed; now reading positive-duration observations.', flush=True)
    result = evaluate(prediction)
    if args.write:
        with args.write.open('x') as handle:
            handle.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text()), result)
    print(json.dumps({'all_pass': result['all_pass'],
                      'preparation': result['calibration']['parameters'],
                      'refinements': result['numerical_refinements'],
                      'scores': {row['name']: row['comparison']['combined'] for row in result['runs']}}, indent=2))


if __name__ == '__main__':
    main()

"""Execute changed-potential trajectories, then score the complete later data."""
import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical
from .experiment import predict, evaluate
from .replay import compare_receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    output = parser.add_mutually_exclusive_group()
    output.add_argument('--write', type=Path)
    output.add_argument('--check', type=Path)
    parser.add_argument('--predictions-write', type=Path)
    args = parser.parse_args()
    prediction = predict(progress=print)
    if args.predictions_write:
        with args.predictions_write.open('x') as handle:
            handle.write(canonical(prediction))
    print('All predictions fixed; opening positive-time derived measurements.', flush=True)
    result = evaluate(prediction)
    if args.write:
        with args.write.open('x') as handle:
            handle.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text()), result)
    print(json.dumps({'all_pass': result['all_pass'], 'refinements': result['numerical_refinements'],
                      'scores': {k: v['metrics'] for k, v in result['scores'].items()}}, indent=2))


if __name__ == '__main__':
    main()

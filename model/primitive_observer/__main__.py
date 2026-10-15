"""Execute the joined primitive observer and its fixed empirical comparison."""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import digest, run


def compare(expected, actual):
    """Validate each local seal, then compare its primitive numeric values."""
    for result in (expected, actual):
        if (result['all_pass'] is not True or not result['checks']
                or any(value is not True for value in result['checks'].values())):
            raise ValueError('the finite composition contract failed')
        laboratory = result['laboratory']
        if digest(laboratory['prediction']) != laboratory['prediction_sha256_before_scoring']:
            raise ValueError('laboratory prediction seal does not bind its own inputs')
    left, right = deepcopy(expected), deepcopy(actual)
    # Source contractions have bounded portable floating differences; their
    # exact local byte seals are meaningful only with their own payload.
    for result in (left, right):
        result['laboratory'].pop('prediction_sha256_before_scoring')
    compare_receipt(left, right)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run(progress=print)
    compare(result, result)
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare(strict_json.load(args.check), result)
    print(json.dumps({'all_pass': result['all_pass'],
                      'observations': result['laboratory']['scores']['metrics']['combined']['observations'],
                      'position_estimator_rmse': result['laboratory']['scores']['metrics']['combined']['position_estimator_rmse'],
                      'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()

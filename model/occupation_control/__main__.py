"""Reproduce retained occupation-control results without native compilation."""
import argparse
import json
import math
from pathlib import Path

from tools import strict_json

from .experiment import run


def produce():
    return {'schema': 'oph.occupation_control.results.v1',
            'purpose': 'unrestricted-versus-finite apparatus error control; no empirical fit',
            'runs': [run(dimension=d, case='a', steps=50, record_steps=[0, 2, 10, 25, 50])
                     for d in (32, 64)]}


def compare(actual, expected, path='root'):
    if type(actual) is not type(expected):
        raise ValueError(f'{path}: type changed')
    if isinstance(actual, dict):
        if actual.keys() != expected.keys():
            raise ValueError(f'{path}: keys changed')
        for key in actual:
            compare(actual[key], expected[key], f'{path}.{key}')
    elif isinstance(actual, list):
        if len(actual) != len(expected):
            raise ValueError(f'{path}: length changed')
        for index, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, f'{path}[{index}]')
    elif isinstance(actual, float):
        if not math.isfinite(actual) or not math.isfinite(expected):
            raise ValueError(f'{path}: nonfinite numerical result')
        if abs(actual-expected) > 2e-10+2e-8*abs(expected):
            raise ValueError(f'{path}: numerical result changed: {actual} versus {expected}')
    elif path.endswith('acceptance_exact_binary_sum'):
        # Each run independently certifies its binary polynomial and exact
        # readout sum. Different platforms may produce different last bits.
        from fractions import Fraction
        if abs(Fraction(actual)-Fraction(expected)) > Fraction(2e-10):
            raise ValueError(f'{path}: acceptance changed')
    elif actual != expected:
        raise ValueError(f'{path}: value changed')


def read_receipt(path):
    receipt = strict_json.load(path)
    # Reuse the comparison's recursive finite-value checks, including JSON
    # exponent overflow, before this receipt can enter a numerical replay.
    compare(receipt, receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.output and args.check:
        parser.error('choose output or check')
    result = produce()
    if args.check:
        compare(result, read_receipt(args.check))
        print('Occupation-control replay passed; every fresh bound was recomputed.')
    elif args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        print(f'Wrote {args.output}')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

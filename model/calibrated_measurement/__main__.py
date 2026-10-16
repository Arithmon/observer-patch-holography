"""Replay the packet or analyze independently supplied complete count tables."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from tools import strict_json
from ..common_source.output_certificate import compare_receipt
from ..occupation_control.__main__ import compare
from .bounds import analyze_counts, json_safe
from .experiment import fresh_interface, produce


def compare_generated(expected, result):
    """Only execution diagnostics tolerate ordinary platform roundoff."""
    expected, result = dict(expected), dict(result)
    for key in ('source_expectation_execution_diagnostic', 'source_detector_metadata'):
        compare(result.pop(key), expected.pop(key))
    compare_receipt(expected, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', type=Path)
    mode.add_argument('--counts', type=Path, help='JSON counts, service_mean_errors and failure_probability')
    args = parser.parse_args()
    if args.counts:
        supplied = strict_json.load(args.counts)
        if type(supplied) is not dict or set(supplied) != {'counts', 'service_mean_errors', 'failure_probability'}:
            raise ValueError('counts, all service obligations and a failure probability are required')
        _, interface, _, _ = fresh_interface()
        if (type(supplied['service_mean_errors']) is not dict
                or any(value is not None and type(value) is not str for value in supplied['service_mean_errors'].values())
                or type(supplied['failure_probability']) is not str):
            raise ValueError('service errors and failure probability must be exact rational strings, or null for missing services')
        services = {key: None if value is None else F(value) for key, value in supplied['service_mean_errors'].items()}
        result = json_safe(analyze_counts(supplied['counts'], interface, service_mean_errors=services,
                                         failure_probability=F(supplied['failure_probability'])))
    else:
        result = produce()
    if args.check:
        expected = strict_json.load(args.check)
        # Metadata's floating residuals are diagnostic; its hashes, resources
        # and types stay strict. All confidence decisions/radii remain exact.
        compare_generated(expected, result)
        print('Fresh native interface and exact calibrated count certificates replayed.')
    elif args.output:
        with args.output.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
        print(f'Wrote {args.output}')
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()

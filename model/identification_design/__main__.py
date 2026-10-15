"""Execute or replay the inexpensive action-identification design."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from tools import strict_json
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
        compare_receipt(strict_json.load(args.check), result)
    print(json.dumps({'all_pass': result['all_pass'],
                      'added_occupations': result['response_design']['added_occupations'],
                      'native_pulses_per_coherent_trial': result['native_coherent_trials'][0]['resources']['native_pulses'],
                      'coherent_total_variation': result['coherent_total_variation'],
                      'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()

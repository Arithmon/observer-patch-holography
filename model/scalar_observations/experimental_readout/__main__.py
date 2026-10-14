"""Retrospective published-data comparison; no network access during replay."""
import argparse
import json
from pathlib import Path

from ...observer.scalar_response.experiment import canonical, compare_receipt
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
        compare_receipt(json.loads(args.check.read_text()), result)
    summary = result['heldout_comparison']
    print(json.dumps({'all_pass': result['all_pass'],
        'data_kind': result['scope']['data_kind'],
        'calibration_points': len(result['calibration']['axis_means']),
        'withheld_points': summary['withheld_points'],
        'withheld_spin_mean_rmse': summary['spin_mean_rmse'],
        'zero_response_spin_mean_rmse': result['zero_response_baseline']['spin_mean_rmse'],
        'conditional_box_violations': summary['box_violation_count'],
        'coordinate_metadata_discrepancy_resolved': False,
        'seven_mode_dynamics_tested': False}, indent=2))


if __name__ == '__main__':
    main()

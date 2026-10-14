"""Replay spatial interacting observation histories and analytic certificates."""
import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path, help='create a receipt, refusing overwrite')
    group.add_argument('--check', type=Path, help='replay a retained receipt')
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text()), result)
    print(json.dumps({'all_pass': result['all_pass'], 'levels': [
        {key: row[key] for key in ('sites', 'two_particle_dimension',
            'interaction_deletion_history_tv', 'field_feedback_deletion_history_tv')}
        for row in result['levels']], 'scope': result['scope']}, indent=2))


if __name__ == '__main__':
    main()

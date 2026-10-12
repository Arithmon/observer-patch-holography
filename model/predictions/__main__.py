"""Replay predictions or evaluate a complete, independently supplied count file."""
import argparse
import json
from pathlib import Path

from . import protocol as p


def _pairs(items):
    result = {}
    for name, value in items:
        if name in result:
            raise ValueError(f'duplicate JSON field: {name}')
        result[name] = value
    return result


def read_json(path):
    def invalid(value):
        raise ValueError(f'nonfinite JSON value: {value}')
    return json.loads(Path(path).read_text(encoding='utf-8'),
                      object_pairs_hook=_pairs, parse_constant=invalid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_mutually_exclusive_group()
    commands.add_argument('--check', type=Path, help='replay the retained software receipt exactly')
    commands.add_argument('--write', type=Path, help='write a new software receipt, refusing overwrite')
    commands.add_argument('--template', type=Path, help='write an unfilled independent-count template')
    commands.add_argument('--evaluate', type=Path, help='score every supplied count against the fixed protocol')
    commands.add_argument('--demo', choices=('C', 'E'), help='generate labeled synthetic counts and verdict')
    args = parser.parse_args()
    plan = p.checked_protocol()
    if args.template:
        with args.template.open('x', encoding='utf-8') as stream:
            stream.write(p.canonical(p.data_template(plan)))
        result = {'status': 'unfilled_template_written', 'path': str(args.template),
                  'protocol_sha256': p.digest(plan), 'empirical_confirmation': False}
    elif args.evaluate:
        result = p.evaluate(read_json(args.evaluate), plan)
    elif args.demo:
        data = p.synthetic_counts(args.demo, 20261012, plan)
        result = {'data': data, 'verdict': p.evaluate(data, plan)}
    else:
        result = p.run()
        if args.check and result != read_json(args.check):
            raise ValueError('retained prediction receipt does not reproduce exactly')
        if args.write:
            with args.write.open('x', encoding='utf-8') as stream:
                stream.write(p.canonical(result))
        result = {'software_checks_pass': result['software_checks_pass'],
                  'protocol_sha256': result['protocol_sha256'],
                  'no_data_verdict': result['no_data_verdict'],
                  'synthetic_controls': {name: row['verdict']['candidate_verdicts']
                                         for name, row in result['synthetic_controls'].items()},
                  'retrospective_endpoint_audit': result['retrospective_endpoint_audit'],
                  'empirical_confirmation': False}
    print(p.canonical(result), end='')


if __name__ == '__main__':
    main()

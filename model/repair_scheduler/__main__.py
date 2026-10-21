"""Replay the source-derived scheduler and complete continuing observer trial."""
import argparse
import math
from pathlib import Path

from tools import strict_json
from model.observer.scalar_response.experiment import canonical, compare_receipt
from .construction import public_records, run


def _validate_receipt_domains(receipt):
    """Probability/error diagnostics have domains beyond replay closeness.

    This packet has no signed floating diagnostics. Complete numerical laws
    use the same 2e-10 absolute tolerance as replay for their mass sums; the
    structurally impossible public record remains exactly zero. The shared
    comparator still checks every type, label, count and exact bound.
    """
    def walk(value, path='root'):
        if type(value) is float:
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'Nonfinite or negative diagnostic at {path}')
        elif type(value) is dict:
            for key, child in value.items():
                if key in ('probability', 'timeout_probability_diagnostic'):
                    probability(child, f'{path}.{key}')
                walk(child, f'{path}.{key}')
        elif type(value) is list:
            for index, child in enumerate(value):
                walk(child, f'{path}[{index}]')

    def probability(value, path):
        if type(value) is not float or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'Probability outside [0,1] at {path}')

    def grid(value, path):
        if type(value) is not list or not value or any(type(row) is not list or not row for row in value):
            raise ValueError(f'Nonempty probability grid required at {path}')
        values = [entry for row in value for entry in row]
        for entry in values:
            probability(entry, path)
        return values

    def rows(field):
        value = receipt.get(field, [])
        if type(value) is not list or any(type(row) is not dict for row in value):
            raise ValueError(f'Receipt rows required at {field}')
        return value

    if type(receipt) is not dict:
        raise ValueError('A receipt object is required')
    walk(receipt)
    order = [list(record) for record in public_records()]
    impossible = order.index([0, 0, 1])
    for row in rows('continuing_histories'):
        values = grid(row.get('joint_probabilities'), 'joint_probabilities')
        joint = row['joint_probabilities']
        if (row.get('public_record_order') != order or len(joint) != len(order)
                or any(len(line) != len(order) for line in joint)):
            raise ValueError('The complete joint law requires the declared public record order')
        if any(joint[impossible][i] != 0 or joint[i][impossible] != 0 for i in range(len(order))):
            raise ValueError('The impossible public record must have exactly zero joint mass')
        if abs(math.fsum(values)-1) > 2e-10:
            raise ValueError('Complete joint probabilities do not sum to one')
    for row in rows('auxiliary_executions'):
        values = grid(row.get('probabilities_by_round_then_address'), 'auxiliary probabilities')
        timeout = row.get('timeout_probability_diagnostic')
        probability(timeout, 'auxiliary timeout')
        if abs(math.fsum((*values, timeout))-1) > 2e-10:
            raise ValueError('Accepted auxiliary probabilities and timeout do not sum to one')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    _validate_receipt_domains(result)
    if args.check:
        retained = strict_json.load(args.check)
        compare_receipt(retained, result)
        _validate_receipt_domains(retained)
        print('Repair-derived scheduler and complete histories replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()

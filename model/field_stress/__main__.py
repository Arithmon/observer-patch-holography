"""Execute/check the finite same-field energy and entropy receipt."""
import argparse
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = parser.add_mutually_exclusive_group()
    choices.add_argument('--write',type=Path)
    choices.add_argument('--check',type=Path)
    args = parser.parse_args()
    result = run()
    if not result['all_pass']:
        raise ArithmeticError('same-field energy/entropy checks failed')
    if args.write:
        with args.write.open('x',encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()),result)
    print(json.dumps({key:result[key] for key in ('all_pass','energy_continuity','premise_attachment')},indent=2))


if __name__ == '__main__':
    main()

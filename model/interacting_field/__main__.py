"""Execute or independently replay the three-box interacting CAR receipt."""
import argparse
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    choice=parser.add_mutually_exclusive_group()
    choice.add_argument('--write',type=Path,help='write a new receipt; never overwrite evidence')
    choice.add_argument('--check',type=Path,help='recompute and compare every retained field')
    args=parser.parse_args()
    result=run()
    if not result['all_pass']:
        raise ArithmeticError('interacting field checks failed')
    if args.write:
        with args.write.open('x',encoding='utf-8') as stream:stream.write(canonical(result))
    if args.check:
        # Reject duplicate keys/nonfinite constants before numeric comparison.
        expected=strict_json.loads(args.check.read_text(encoding='utf-8'))
        compare_receipt(expected,result)
    print(json.dumps({'all_pass':result['all_pass'],
        'levels':[{k:r[k] for k in ('sites','CAR_modes','two_particle_dimension',
            'interaction_deletion_history_TV','source_reference_history_TV',
            'ideal_source_to_infinite_lattice_complete_history_bound')} for r in result['levels']],
        'connected_four_point':result['non_Gaussian_witness']['connected_four_point'],
        'net_kinetic_node_chirality':result['net_node_chirality']},indent=2,allow_nan=False))


if __name__=='__main__':main()

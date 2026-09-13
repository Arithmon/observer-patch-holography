"""Exact minimal regional collar certificate for the existing q=5 scalar action.

This is an algebraic access certificate, not a quantum-outcome simulation.
The companion verifier uses the independently assembled sqrt(5) parent model
and checks rank factorizations rather than repeating this elimination.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / 'source_scalar_execution'))
from scalar_execution_algebra import Q
from source_scalar_execution import model, TAU

OUTPUT = HERE / 'regional_time_slice_receipt.json'
PARENTS = (
    'code/source_scalar_execution/scalar_execution_algebra.py',
    'code/source_scalar_execution/source_scalar_execution.py',
    'code/source_scalar_execution/verify_source_scalar_execution.py',
    'paper/tex_fragments/SOURCE_SCALAR_EXECUTION.tex',
    'paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex',
)
SCOPE = {
    'same_full_q5_mass_and_force_operator': True,
    'exact_algebraic_rank_and_reconstruction': True,
    'all_real_weyl_parameters_required': True,
    'classical_records_are_not_weyl_access': True,
    'quantum_outcomes_attested': False,
    'source_selected_action_or_regions': False,
    'physical_time_slice_or_clock': False,
    'continuum_regional_net': False,
}


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'),
                       allow_nan=False) + '\n').encode('ascii')


def regions():
    sites = list(product(range(4), repeat=3))
    selectors = {
        'empty': lambda s: False,
        'corner': lambda s: s == (0, 0, 0),
        'interior_singleton': lambda s: s == (1, 1, 1),
        'central_cube': lambda s: all(i in (1, 2) for i in s),
        'left_half': lambda s: s[0] < 2,
        'interior_line': lambda s: s[1:] == (1, 1),
        'checkerboard': lambda s: sum(s) % 2 == 0,
        'whole': lambda s: True,
    }
    return {name: [i for i, s in enumerate(sites) if predicate(s)]
            for name, predicate in selectors.items()}


def rref_witness(matrix, width):
    n = len(matrix)
    a = [row[:] for row in matrix]
    transform = [[Q(int(i == j)) for j in range(n)] for i in range(n)]
    pivots = []
    for j in range(width):
        k = len(pivots)
        pivot = next((i for i in range(k, n) if a[i][j] != Q()), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        transform[k], transform[pivot] = transform[pivot], transform[k]
        divisor = a[k][j]
        a[k] = [x / divisor for x in a[k]]
        transform[k] = [x / divisor for x in transform[k]]
        for i in range(n):
            if i != k and a[i][j] != Q():
                factor = a[i][j]
                a[i] = [x - factor*y for x, y in zip(a[i], a[k])]
                transform[i] = [x - factor*y for x, y in zip(transform[i], transform[k])]
        pivots.append(j)
    return a[:len(pivots)], transform[:len(pivots)], pivots


def encode(matrix):
    return [[x.encode() for x in row] for row in matrix]


def region_certificate(name, indices, masses, operator):
    outside = [i for i in range(64) if i not in indices]
    b = [[operator[i].get(j, Q()) for j in outside] for i in indices]
    c, e, pivots = rref_witness(b, len(outside))
    d = [[row[j] for j in pivots] for row in b]
    active = [j for j in outside if any(operator[i].get(j, Q()) != Q() for i in indices)]
    rank = len(c)
    counterexamples = []
    for column, pivot in enumerate(pivots):
        # C has identity pivot columns. Deleting row column makes this
        # exterior perturbation invisible, but its local velocity is nonzero.
        v = [TAU / 2 * row[pivot] for row in b]
        counterexamples.append({
            'deleted_collar_row': column,
            'exterior_unit_site': outside[pivot],
            'local_velocity_shift_Qphi': [x.encode() for x in v],
        })
    return {
        'name': name, 'region_sites': indices, 'exterior_sites': outside,
        'coordinate_collar_sites': active, 'minimal_linear_collar_dimension': rank,
        'two_layer_real_dimension': 2*len(indices),
        'two_layer_plus_collar_real_dimension': 2*len(indices)+rank,
        'no_collar_recovers_original_regional_algebra': rank == 0,
        'collar_rows_Qphi': encode(c), 'collar_from_force_rows_Qphi': encode(e),
        'force_from_collar_rows_Qphi': encode(d), 'pivot_columns': pivots,
        'q_next_cross_bracket_over_i_hbar_Qphi': [(TAU/masses[i]).encode() for i in indices],
        'missing_row_controls': counterexamples,
    }


def build():
    _, _, masses, rows, _, _, _ = model()
    operator = [dict(row) for row in rows]
    for i in range(64):
        for j in range(64):
            if masses[i]*operator[i].get(j, Q()) != masses[j]*operator[j].get(i, Q()):
                raise ValueError('full mass-weighted symmetry failed')
    return {
        'schema': 'oph.source-scalar-regional-time-slice.v1',
        'parent_sha256': {path: sha256((ROOT/path).read_bytes()).hexdigest() for path in PARENTS},
        'scope': SCOPE, 'sites': 64, 'tau_Qphi': TAU.encode(),
        'mass_Qphi': [m.encode() for m in masses],
        'regions': [region_certificate(name, r, masses, operator) for name, r in regions().items()],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    receipt = build()
    if args.write:
        OUTPUT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps([{k: r[k] for k in ('name', 'minimal_linear_collar_dimension',
                      'two_layer_real_dimension')} for r in receipt['regions']], indent=2))


if __name__ == '__main__':
    main()

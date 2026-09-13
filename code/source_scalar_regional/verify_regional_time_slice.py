"""Independent exact verification of regional collar rank and Weyl brackets.

The producer and its Q(phi) arithmetic are never imported. The preexisting
independent scalar verifier assembles the full stiffness from edge energies
in Q(sqrt(5)). Rank is checked through two factorization identities and an
identity minor, with no elimination algorithm shared with the producer.
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
from verify_source_scalar_execution import R, TAU, model, parse, require, load, equal

OUTPUT = HERE/'regional_time_slice_receipt.json'
PARENTS = (
    'code/source_scalar_execution/scalar_execution_algebra.py',
    'code/source_scalar_execution/source_scalar_execution.py',
    'code/source_scalar_execution/verify_source_scalar_execution.py',
    'paper/tex_fragments/SOURCE_SCALAR_EXECUTION.tex',
    'paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex',
)


def matrix(value, nr, nc, label):
    require(type(value) is list and len(value) == nr, label+' row count')
    require(all(type(row) is list and len(row) == nc for row in value), label+' column count')
    return [[parse(x) for x in row] for row in value]


def mm(a, b, nr, nk, nc):
    return [[sum((a[i][k]*b[k][j] for k in range(nk)), R())
             for j in range(nc)] for i in range(nr)]


def declarations():
    # Direct integer coordinate descriptions, independently expressed.
    grid = list(product(range(4), repeat=3))
    return (
        ('empty', []), ('corner', [0]), ('interior_singleton', [21]),
        ('central_cube', [i for i, x in enumerate(grid) if min(x) >= 1 and max(x) <= 2]),
        ('left_half', list(range(32))),
        ('interior_line', [5, 21, 37, 53]),
        ('checkerboard', [i for i, x in enumerate(grid) if (x[0]+x[1]+x[2]) % 2 == 0]),
        ('whole', list(range(64))),
    )


def verify(receipt, root=ROOT):
    equal(sorted(receipt), sorted(('schema', 'parent_sha256', 'scope', 'sites',
                                   'tau_Qphi', 'mass_Qphi', 'regions')), 'top schema')
    equal(receipt['schema'], 'oph.source-scalar-regional-time-slice.v1', 'schema identifier')
    equal(receipt['sites'], 64, 'site count')
    equal(receipt['tau_Qphi'], TAU.encode(), 'same action tick')
    equal(receipt['parent_sha256'], {p: sha256((root/p).read_bytes()).hexdigest() for p in PARENTS},
          'parent source identity')
    equal(receipt['scope'], {
        'same_full_q5_mass_and_force_operator': True,
        'exact_algebraic_rank_and_reconstruction': True,
        'all_real_weyl_parameters_required': True,
        'classical_records_are_not_weyl_access': True,
        'quantum_outcomes_attested': False,
        'source_selected_action_or_regions': False,
        'physical_time_slice_or_clock': False,
        'continuum_regional_net': False,
    }, 'scope must not be promoted')
    masses, rows, _, _, _ = model(root)
    equal(receipt['mass_Qphi'], [m.encode() for m in masses], 'same mass matrix')
    operator = [dict(row) for row in rows]
    require(all(m.sign() > 0 for m in masses), 'positive masses')
    # Exact CCR of actual q_i at the next split layer, retaining all 64
    # coordinates: [q_i^+,q_j^+]/(i hbar) = 0, not a compressed claim.
    for i in range(64):
        for j in range(64):
            aij = operator[i].get(j, R())
            aji = operator[j].get(i, R())
            require(aij/masses[j] == aji/masses[i], 'full next-field commuting bracket')
    expected_regions = declarations()
    require(type(receipt['regions']) is list and len(receipt['regions']) == len(expected_regions),
            'retain every declared region')
    summary = []
    for entry, (name, inside) in zip(receipt['regions'], expected_regions):
        equal(sorted(entry), sorted(('name', 'region_sites', 'exterior_sites', 'coordinate_collar_sites',
            'minimal_linear_collar_dimension', 'two_layer_real_dimension',
            'two_layer_plus_collar_real_dimension', 'no_collar_recovers_original_regional_algebra',
            'collar_rows_Qphi', 'collar_from_force_rows_Qphi', 'force_from_collar_rows_Qphi',
            'pivot_columns', 'q_next_cross_bracket_over_i_hbar_Qphi', 'missing_row_controls')),
              'region schema')
        equal(entry['name'], name, 'region name/order')
        equal(entry['region_sites'], inside, 'declared regional support')
        outside = [i for i in range(64) if i not in inside]
        equal(entry['exterior_sites'], outside, 'complete complement')
        nr, nc = len(inside), len(outside)
        k = entry['minimal_linear_collar_dimension']
        require(type(k) is int and 0 <= k <= min(nr, nc), 'collar rank range/type')
        b = [[operator[i].get(j, R()) for j in outside] for i in inside]
        c = matrix(entry['collar_rows_Qphi'], k, nc, 'collar')
        e = matrix(entry['collar_from_force_rows_Qphi'], k, nr, 'forward factor')
        d = matrix(entry['force_from_collar_rows_Qphi'], nr, k, 'reverse factor')
        pivots = entry['pivot_columns']
        require(type(pivots) is list and len(pivots) == k and
                all(type(j) is int and 0 <= j < nc for j in pivots) and
                pivots == sorted(set(pivots)), 'independent pivot list')
        # These certify both row-space inclusions. The identity minor
        # certifies rank(C)=k and hence rank(B)=k without numeric tolerance.
        require(mm(d, c, nr, k, nc) == b, 'full exterior force not reconstructed')
        require(mm(e, b, k, nr, nc) == c, 'collar not in exterior force row space')
        require([[row[j] for j in pivots] for row in c] ==
                [[R(int(i == j)) for j in range(k)] for i in range(k)], 'identity rank minor')
        active = [j for j in outside if any(operator[i].get(j, R()) != R() for i in inside)]
        equal(entry['coordinate_collar_sites'], active, 'all incident exterior sites')
        equal(entry['two_layer_real_dimension'], 2*nr, 'two-layer dimension')
        equal(entry['two_layer_plus_collar_real_dimension'], 2*nr+k, 'collar dimension')
        equal(entry['no_collar_recovers_original_regional_algebra'], k == 0, 'no-collar boundary')
        equal(entry['q_next_cross_bracket_over_i_hbar_Qphi'], [(TAU/masses[i]).encode() for i in inside],
              'nonzero local cross-time CCR')
        controls = entry['missing_row_controls']
        require(type(controls) is list and len(controls) == k, 'every minimal collar row has a control')
        for deleted, (control, pivot) in enumerate(zip(controls, pivots)):
            velocity = [TAU/2*row[pivot] for row in b]
            require(any(v != R() for v in velocity), 'control has nonzero local velocity')
            for j in range(k):
                require(c[j][pivot] == R(int(j == deleted)), 'remaining readouts see perturbation')
            for i, v in enumerate(velocity):
                require(TAU*v-TAU**2/2*b[i][pivot] == R(), 'next-field cancellation')
            equal(control, {'deleted_collar_row': deleted, 'exterior_unit_site': outside[pivot],
                            'local_velocity_shift_Qphi': [v.encode() for v in velocity]}, 'control witness')
        summary.append({'region': name, 'sites': nr, 'coordinate_collar_sites': len(active),
                        'minimal_linear_collar': k})
    return {'verified': True, 'region_count': len(summary),
            'missing_row_counterexamples': sum(r['minimal_linear_collar'] for r in summary),
            'regions': summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), indent=2))


if __name__ == '__main__':
    main()

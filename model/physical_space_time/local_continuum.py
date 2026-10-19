"""Compact native spatial Green limits and shrinking action range.

The source topology and its ancestry dictionary are the existing ones.  The
result upgrades their common weak Green subsequence to operator-norm
convergence, with an injective compact limit and a positive inverse generator.
It does not assert a local Dirichlet generator or a physical spatial metric.
The all-level proof is in OCCUPATION_FIELD.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path
import sympy as sp

from . import continuum_field as continuum
from . import source_green as source
from . import uniform_green as uniform
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*continuum.SOURCE_PATHS,
    'model/physical_space_time/local_continuum.py'))))
GREEN_BOUND = uniform.GREEN_BOUND
PATCH_POINCARE_CONSTANT = F(26)
COARSE_ENERGY_TAIL_CONSTANT = F(624)


def _level(value):
    if type(value) is not int or value < 0:
        raise ValueError('an original nonnegative integer is required')
    return value


def compactness_tail_squared(coarse_level):
    """Squared norm bound for (I-P_k)G_l, at every l>=k."""
    return COARSE_ENERGY_TAIL_CONSTANT*GREEN_BOUND/4**_level(coarse_level)


def covariance_eigenvalue_tail(coarse_level):
    """Bound on positive eigenvalue number 20*4^k and every later one."""
    return COARSE_ENERGY_TAIL_CONSTANT/4**_level(coarse_level)


def separated_action_bound(level, cell_distance):
    """Rescaled Schur bilinear bound for unit vectors at cell distance d>=1."""
    level = _level(level)
    if type(cell_distance) is not int or cell_distance < 1:
        raise ValueError('a positive integer cell distance is required')
    return F(4**level, 2**cell_distance)


def macroscopic_action_exponent(level, separation):
    """Base-two exponent of the separated-support action upper bound.

    Separation is intrinsic distance on the auxiliary unit-circumradius flat
    polyhedron.  Adjacent centroids are less than 2^-l apart; each fine cell
    has diameter less than 2*2^-l.  This is a comparison length, not metres.
    An exponent avoids allocating a doubly exponential denominator.
    """
    level = _level(level)
    if type(separation) is not F or separation <= 0:
        raise ValueError('a positive exact Fraction separation is required')
    lower = separation*2**level-4
    distance = max(1, -((-lower.numerator)//lower.denominator))
    if lower <= 0:
        raise ValueError('the cells are not yet certified disjoint at this level')
    return 2*level-distance


def triangle_patch(depth):
    """Exact upright/downward cells and dual edges in one descendant triangle."""
    depth = _level(depth)
    if depth > 5:
        raise ValueError('the finite stencil control is limited to depth five')
    m = 2**depth
    up = tuple((i, j) for i in range(m) for j in range(m-i))
    down = tuple((i, j) for i in range(m-1) for j in range(m-1-i))
    vertices = tuple(('u', *v) for v in up)+tuple(('d', *v) for v in down)
    index = {v: i for i, v in enumerate(vertices)}
    edges = tuple(sorted(tuple(sorted((index[('d', i, j)], index[('u', x, y)])))
        for i, j in down for x, y in ((i, j), (i+1, j), (i, j+1))))
    return m, vertices, edges


def reflected_square_control(depth):
    """The square has exactly twice each upright-triangle edge energy."""
    m, vertices, edges = triangle_patch(depth)
    def reflect(i, j):
        return (i, j) if i+j < m else (m-1-j, m-1-i)
    upright_edges = Counter()
    for i in range(m):
        for j in range(m-i):
            for x, y in ((i+1, j), (i, j+1)):
                if x+y < m:
                    upright_edges[tuple(sorted(((i, j), (x, y))))] = 2
    reflected_edges = Counter()
    for i in range(m):
        for j in range(m):
            for x, y in ((i+1, j), (i, j+1)):
                if x < m and y < m:
                    reflected_edges[tuple(sorted((reflect(i, j), reflect(x, y))))] += 1
    if reflected_edges != upright_edges:
        raise ArithmeticError('square reflection changed the exact edge multiplicities')
    multiplicity = Counter(reflect(i, j) for i in range(m) for j in range(m))
    if (set(multiplicity) != {(i, j) for kind, i, j in vertices if kind == 'u'}
            or not all(n in (1, 2) for n in multiplicity.values())):
        raise ArithmeticError('upright cells are not covered once or twice')
    return {'depth': depth, 'side_subdivisions': m, 'descendant_cells': len(vertices),
        'internal_dual_edges': len(edges), 'square_edges': sum(reflected_edges.values()),
        'square_edge_energy_is_twice_upright_energy': True,
        'all_upright_cells_have_multiplicity_one_or_two': True}


def source_patch_control(depth):
    """Bind every base-face descendant adjacency to the independent stencil."""
    if _level(depth) > 3:
        raise ValueError('native ancestry binding controls use depths zero through three')
    m, vertices, edges = triangle_patch(depth)
    # The first two slots here are independent combinatorial triangle
    # coordinates; they are not the physical x/y slots of the carrier chart.
    simple = ((0, 0, 0, 0, 0, 0), (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0))
    faces = ((0, 1, 2),)
    for _ in range(depth):
        simple, faces = uniform.refine(simple, faces)
    lookup = {v: i for i, v in enumerate(vertices)}
    labels = []
    for face in faces:
        corners = tuple((simple[v][0]//2, simple[v][1]//2) for v in face)
        i, j = min(x for x, _ in corners), min(y for _, y in corners)
        kind = 'u' if (i, j) in corners else 'd'
        labels.append(lookup[(kind, i, j)])
    _, actual_faces = uniform.affine_prefix(depth)
    actual_edges = uniform.dual_edges(actual_faces)
    expected = set(edges)
    for base in range(20):
        lower, upper = base*m*m, (base+1)*m*m
        internal = {tuple(sorted((labels[a-lower], labels[b-lower])))
                    for a, b in actual_edges if lower <= a < upper and lower <= b < upper}
        if internal != expected:
            raise ArithmeticError('native descendants differ from the triangular patch stencil')
    return {'depth': depth, 'base_patches_checked': 20,
        'cell_ancestry_order_and_all_internal_edges_checked': True}


def exact_patch_energy_control(depth):
    """Independent rational LDL certificate for the full patch variance bound."""
    m, vertices, edges = triangle_patch(depth)
    if depth > 3:
        raise ValueError('exact matrix energy controls are limited to depth three')
    n = len(vertices)
    laplacian = sp.zeros(n)
    for a, b in edges:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    constant = sp.ones(n)/n
    certificate = 26*m*m*laplacian-sp.eye(n)+2*constant
    # Constant direction has eigenvalue one; the other directions certify
    # 26*m² L_patch >= I-P_constants without a floating eigenvalue test.
    _, pivots = certificate.LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(n)):
        raise ArithmeticError('the exact triangular patch Poincare bound failed')
    return {'depth': depth, 'cells': n,
        'variance_energy_multiplier': 26*m*m,
        'rational_LDL_positive_pivots': n}


def native_operator_control(level):
    """Native support, compactness and separated-action finite controls."""
    if _level(level) > 2:
        raise ValueError('native operator controls use levels zero through two')
    graph, _, r, integer_f, denominator = source.matrices(level)
    cell_l = r.T@r
    adjacency = -(cell_l-np.diag(np.diag(cell_l)))
    distances = shortest_path(adjacency, directed=False, unweighted=True)
    numerator_t = 2*denominator*np.eye(len(integer_f), dtype=np.int64)-integer_f
    term = r.copy()
    support_rows = []
    row_norm = max(sum(abs(int(value)) for value in row) for row in numerator_t)
    incidence_column_norm = max(sum(abs(int(value)) for value in column) for column in r.T)
    integer_limit = int(np.iinfo(np.int64).max)
    for power in range(4):
        term_norm = max(abs(int(value)) for value in term.flat)
        if incidence_column_norm*term_norm > integer_limit:
            raise ArithmeticError('exact polynomial coefficient could overflow int64')
        coefficient = r.T@term
        if np.any(coefficient[distances > power+1]):
            raise ArithmeticError('Schur polynomial has a forbidden long-range coefficient')
        support_rows.append({'power': power, 'maximum_cell_range': power+1,
            'integer_polynomial_support_checked': True})
        if power < 3:
            if row_norm*term_norm > integer_limit:
                raise ArithmeticError('exact polynomial power could overflow int64')
            term = numerator_t@term
    schur = r.T@np.linalg.solve(integer_f/denominator, r)/12
    green = source._green(schur)/4**level
    compact_rows = []
    for coarse in range(level+1):
        descendants = 4**(level-coarse)
        projected = np.repeat(green.reshape(20*4**coarse, descendants, graph.cells).mean(axis=1),
                              descendants, axis=0)
        squared = float(np.linalg.norm(green-projected, 2)**2)
        if squared > float(compactness_tail_squared(coarse))+1e-9:
            raise ArithmeticError('native Green violates the coarse compactness bound')
        compact_rows.append({'coarse_level': coarse,
            'tail_norm_squared_diagnostic': squared,
            'tail_norm_squared_upper': str(compactness_tail_squared(coarse))})
    distance = int(distances.max())
    a, b = map(int, np.argwhere(distances == distance)[0])
    coupling = float(abs(4**level*schur[a, b]))
    bound = separated_action_bound(level, distance)
    if coupling > float(bound)+1e-12:
        raise ArithmeticError('native distant action exceeded its exact bound')
    return {'level': level, 'cells': graph.cells,
        'exact_polynomial_support': support_rows, 'compactness_controls': compact_rows,
        'distant_native_pair': [a, b], 'cell_distance': distance,
        'absolute_rescaled_coupling_diagnostic': coupling,
        'absolute_rescaled_coupling_upper': str(bound)}


def non_markov_control():
    """Averages are not a vertex Kron reduction: a positive off-diagonal survives."""
    _, _, r, integer_f, denominator = source.matrices(0)
    ri, fi = sp.Matrix(r), sp.Matrix(integer_f)
    entry = (ri[:, 5].T*fi.inv()*ri[:, 9])[0]*denominator/12
    if entry <= 0:
        raise ArithmeticError('the actual source mean-elimination countercontrol changed')
    return {'level': 0, 'cells': [5, 9], 'exact_positive_Schur_off_diagonal': str(entry),
        'finite_cell_mean_generator_is_a_Markov_generator': False,
        'continuum_locality_or_nonlocality_decided_by_this_control': False}


def run():
    return {'schema': 'oph.native_compact_spatial_generator.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'contract': {'source': 'existing spherical shared-midpoint carrier joins with distinct ports',
            'space': 'mean-zero L2 of normalized auxiliary flat-polyhedron area',
            'observation_dictionary': 'normalized ancestral cell means',
            'limit': 'the same common spatial subsequence as the native Green field',
            'physical_spatial_metric_or_radial_coordinate_supplied_by_this_result': False},
        'exact_patch_controls': [reflected_square_control(k) | source_patch_control(k)
                                 | exact_patch_energy_control(k) for k in range(3)],
        'native_operator_controls': [native_operator_control(k) for k in range(3)],
        'mean_elimination_control': non_markov_control(),
        'all_level_bounds': {
            'patch_variance_energy_multiplier': '26*4^(l-k)',
            'coarse_projection_energy_multiplier': '624*4^(-k)',
            'Green_norm_upper': str(GREEN_BOUND),
            'Green_tail_squared': {str(k): str(compactness_tail_squared(k)) for k in (0, 4, 20)},
            'covariance_eigenvalue_number_20_times_4_power_k_upper': {
                str(k): str(covariance_eigenvalue_tail(k)) for k in (0, 4, 20)},
            'separated_cell_support_bilinear_action': '4^l*2^(-d)*||u||*||v||, d>=1',
            'macroscopic_separation_one_quarter_action_base_two_exponents': {
                str(k): macroscopic_action_exponent(k, F(1, 4)) for k in (5, 8, 12)}},
        'theorem': {'common_Green_subsequence_converges_in_operator_norm': True,
            'limit_is_compact_positive_and_injective_on_mean_zero_L2': True,
            'inverse_is_positive_self_adjoint_spatial_generator': True,
            'finite_rescaled_action_couplings_at_separated_continuum_supports_vanish': True,
            'unique_whole_sequence_generator_selected': False,
            'limit_identified_as_local_Dirichlet_or_round_sphere_Laplacian': False,
            'three_dimensional_space_or_Lorentzian_metric_derived': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native compact spatial-generator bounds and source controls verified.')
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
